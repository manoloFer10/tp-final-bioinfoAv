import re
import time
import threading
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm


def _make_session(retries: int = 5, backoff_factor: float = 1.0) -> requests.Session:
    """Crea una sesión con reintentos automáticos y backoff exponencial."""
    session = requests.Session()
    retry = Retry(
        total=retries,
        backoff_factor=backoff_factor,   # espera 1, 2, 4, 8, 16 s entre reintentos
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


_thread_local = threading.local()

def _get_session() -> requests.Session:
    if not hasattr(_thread_local, "session"):
        _thread_local.session = _make_session()
    return _thread_local.session

def _extraer_accession(gene_id_str: str) -> str | None:
    """Extrae el accession de strings como 'Protein_ID: WP_006487096.1',
    o devuelve el valor directamente si ya es un accession (ej. 'WP_006487096.1')."""
    s = str(gene_id_str).strip()
    match = re.search(r'Protein_ID:\s*(\S+)', s)
    if match:
        return match.group(1)
    # Si no tiene el prefijo, asumir que el valor ya es el accession
    return s if s and s.lower() != 'nan' else None


def _mapear_a_uniprot(accessions: list[str]) -> dict[str, str]:
    """
    Envía los accessions a la API asíncrona de UniProt ID Mapping y devuelve
    un diccionario  {accession_original: uniprot_primary_accession}.

    UniProt distingue entre bases de datos de origen:
      - WP_* / NP_* → "RefSeq_Protein"
      - CAR* / CAO* (EMBL protein) → "EMBL-GenBank-DDBJ_CDS"
    """
    RUN_URL    = "https://rest.uniprot.org/idmapping/run"
    RESULT_URL = "https://rest.uniprot.org/idmapping/results/{job_id}"
    BATCH_SIZE = 100         
    MAX_POLLS  = 5
    POLL_SLEEP = 30           
    wp_accs   = [a for a in accessions if a and a.startswith(("WP_", "NP_", "XP_"))]
    embl_accs = [a for a in accessions if a and not a.startswith(("WP_", "NP_", "XP_"))]

    mapping: dict[str, str] = {}

    def _submit_and_collect(accs_batch: list[str], from_db: str):
        if not accs_batch:
            return
        r = requests.post(RUN_URL, data={
            "ids":  ",".join(accs_batch),
            "from": from_db,
            "to":   "UniProtKB",
        })
        r.raise_for_status()
        job_id = r.json()["jobId"]

        for _ in range(MAX_POLLS):
            time.sleep(POLL_SLEEP)
            res = requests.get(RESULT_URL.format(job_id=job_id) + "?format=json")
            if not res.ok:
                continue
            data = res.json()
            if "results" in data:
                for item in data["results"]:
                    # item["to"] puede ser un dict {"primaryAccession": ...} o un string directo
                    to_val = item["to"]
                    mapping[item["from"]] = to_val["primaryAccession"] if isinstance(to_val, dict) else to_val
                next_url = data.get("nextPage")
                while next_url:
                    nr = requests.get(next_url)
                    nr.raise_for_status()
                    nd = nr.json()
                    for item in nd.get("results", []):
                        to_val = item["to"]
                        mapping[item["from"]] = to_val["primaryAccession"] if isinstance(to_val, dict) else to_val
                    next_url = nd.get("nextPage")
                return  # trabajo completado

        print(f"[WARN] Timeout esperando resultados de UniProt para lote {from_db}.")

    # Procesar en lotes para no saturar la API
    for i in tqdm(range(0, len(wp_accs), BATCH_SIZE), desc="Mapeando WP/NP a UniProt", unit="batch"):
        _submit_and_collect(wp_accs[i:i + BATCH_SIZE], "RefSeq_Protein")
        time.sleep(0.5)

    for i in tqdm(range(0, len(embl_accs), BATCH_SIZE), desc="Mapeando EMBL a UniProt", unit="batch"):
        _submit_and_collect(embl_accs[i:i + BATCH_SIZE], "EMBL-GenBank-DDBJ_CDS")
        time.sleep(0.5)

    return mapping


def _obtener_dominios_pfam(uniprot_id: str) -> list[dict]:
    """
    Consulta la API de InterPro para un UniProt ID y devuelve una lista de dicts:
        {
          'accession': 'PF00001',   # PFAM accession
          'name':      'Domain X',  # nombre legible
          'start':     12,          # posición de inicio en la proteína
          'end':       134          # posición de fin
        }

    Si la proteína no tiene dominios PFAM, devuelve [].

    Endpoint:
      GET https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/UniProt/{uniprot_id}/
    """
    url = (
        f"https://www.ebi.ac.uk/interpro/api/entry/pfam/"
        f"protein/UniProt/{uniprot_id}/?format=json"
    )
    MAX_ATTEMPTS = 4
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            r = _get_session().get(url, timeout=30)
            if r.status_code == 204:   # 204 = No Content → sin dominios PFAM
                return []
            r.raise_for_status()
            data = r.json()
            break
        except requests.exceptions.RequestException as e:
            if attempt == MAX_ATTEMPTS:
                print(f"[ERROR] Fallo al consultar InterPro para {uniprot_id} tras {MAX_ATTEMPTS} intentos: {e}")
                return []
            wait = 2 ** attempt
            print(f"[WARN] Reintento {attempt}/{MAX_ATTEMPTS - 1} para {uniprot_id} en {wait}s ({e})")
            time.sleep(wait)

    dominios = []
    for entry in data.get("results", []):
        pfam_acc  = entry["metadata"]["accession"]   # ej. "PF00696"
        pfam_name = entry["metadata"]["name"]

        # Extraer coordenadas de cada fragmento del dominio
        for prot in entry.get("proteins", []):
            for location in prot.get("entry_protein_locations", []):
                for fragment in location.get("fragments", []):
                    dominios.append({
                        "accession": pfam_acc,
                        "name":      pfam_name,
                        "start":     fragment["start"],
                        "end":       fragment["end"],
                    })

    # Manejar paginación de InterPro (poco probable para una sola proteína, pero robusto)
    next_url = data.get("next")
    while next_url:
        r = _get_session().get(next_url, timeout=30)
        r.raise_for_status()
        data = r.json()
        for entry in data.get("results", []):
            pfam_acc  = entry["metadata"]["accession"]
            pfam_name = entry["metadata"]["name"]
            for prot in entry.get("proteins", []):
                for location in prot.get("entry_protein_locations", []):
                    for fragment in location.get("fragments", []):
                        dominios.append({
                            "accession": pfam_acc,
                            "name":      pfam_name,
                            "start":     fragment["start"],
                            "end":       fragment["end"],
                        })
        next_url = data.get("next")

    return dominios


def obtener_dominios_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Agrega los dominios PFAM a cada proteína en el DataFrame.

    Retorna un nuevo DataFrame con la columna adicional:
      - 'Dominios': lista de dicts con claves 'accession', 'name', 'start', 'end'
                    para cada dominio PFAM encontrado en la proteína.
                    Lista vacía si no se encontraron dominios o no hubo mapeo.

    Pasos:
      1. Extrae los accession IDs de 'gene_id'.
      2. Los mapea a UniProt usando la API de ID Mapping.
      3. Consulta InterPro por cada UniProt accession para recuperar dominios PFAM.
    """
    df = df.copy()

    # Paso 1 — Extraer accessions
    df["_accession"] = df["gene_id"].apply(_extraer_accession)

    # Paso 2 — Mapear a UniProt (batch, eficiente)
    all_accs    = df["_accession"].dropna().unique().tolist()
    uniprot_map = _mapear_a_uniprot(all_accs)          # {acc_ncbi: acc_uniprot}
    df["_uniprot"] = df["_accession"].map(uniprot_map) # NaN si no hubo mapeo

    # Paso 3 — Obtener dominios PFAM en paralelo (ThreadPoolExecutor)
    _cache: dict[str, list[dict]] = {}
    unique_uniprots = df["_uniprot"].dropna().unique().tolist()

    MAX_WORKERS = 8   # EBI tolera bien ~8 conexiones simultáneas
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as pool:
        future_to_uid = {pool.submit(_obtener_dominios_pfam, uid): uid for uid in unique_uniprots}
        with tqdm(total=len(future_to_uid), desc="Obteniendo dominios PFAM", unit="proteína") as pbar:
            for future in as_completed(future_to_uid):
                uid = future_to_uid[future]
                _cache[uid] = future.result()
                pbar.update(1)

    def _fetch_con_cache(uniprot_id) -> list[dict]:
        if pd.isna(uniprot_id):
            return []
        return _cache.get(uniprot_id, [])

    df["Dominios"] = df["_uniprot"].apply(_fetch_con_cache)

    # Limpiar columnas temporales
    df.drop(columns=["_accession", "_uniprot"], inplace=True)
    return df


def obtener_comienzo_fin_dominios(fila_proteina: pd.Series) -> list[tuple[int, int]]:
    """
    Dada una fila del DataFrame que ya contiene la columna 'Dominios'
    (lista de dicts con 'accession', 'name', 'start', 'end'),
    devuelve una lista de tuplas (start, end) en el mismo orden que los dominios.

    Si no hay dominios, devuelve [].

    Ejemplo de uso:
        df["posiciones"] = df.apply(obtener_comienzo_fin_dominios, axis=1)
    """
    dominios = fila_proteina.get("Dominios", [])
    if not dominios:
        return []
    return [(d["start"], d["end"]) for d in dominios]

def obtener_lista_homologos(proteina):
    return []

def generar_MSA(fila_proteina):
    """
    Dada una fila del DataFrame con metadata de la proteína, esta función genera un MSA (Multiple Sequence Alignment).
    """
    msa = None
    return msa
    

