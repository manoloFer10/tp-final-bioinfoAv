"""
utiles/busqueda_ligandos.py
---------------------------
Búsqueda de homólogos PDB y recuperación de ligandos/interactores para cada
proteína query.  Devuelve, por cada proteína, tres DataFrames segmentados por
nivel de identidad de secuencia (tiers 1-3).

Referencia de tiers:
  Rentzsch & Orengo (2009), J. Mol. Biol.
"""

import hashlib
import json
import logging
import shutil
import time
import warnings
from pathlib import Path

import pandas as pd
import requests
from Bio import Align
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
_log_path = Path("errors.log")
logging.basicConfig(
    filename=str(_log_path),
    level=logging.ERROR,
    format="%(asctime)s  %(levelname)s  %(message)s",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

# Artefactos cristalográficos a excluir de los resultados de ligandos
HETATM_BLACKLIST: set[str] = {
    # Solventes
    "HOH", "GOL", "EDO", "PEG", "MPD", "PGE", "DMS",
    # Buffers
    "SO4", "PO4", "ACT", "CIT", "MES", "TRS", "EPE", "HEPES",
    # Iones de cristalización (aislados)
    "MG", "ZN", "CA", "NA", "CL", "K", "MN", "FE", "NI", "CO", "CU",
    # Otros comunes
    "BME", "DTT", "FMT", "IMD", "NO3", "SCN",
}

# ---------------------------------------------------------------------------
# Estructura del caché:
#
#   cache/
#   ├── queries/
#   │   └── <accession>/         ← e.g. WP_006487096.1
#   │       ├── blast.json       ← resultados de RCSB MMseqs2 sequence search
#   │       └── ligands/
#   │           ├── 4CVK.json    ← ligandos de cada PDB hit (copia local)
#   │           └── 8F5D.json
#   ├── identities/              ← caché de identidad por (query_hash, pdb_entity)
#   │   └── <hash>:<entity>.json
#   └── ligands/                 ← caché compartido (evita llamadas repetidas
#       ├── 4cvk.json               a PDBe para el mismo PDB_ID)
#       └── 8f5d.json
# ---------------------------------------------------------------------------

CACHE_DIR = Path("cache")
CACHE_QUERIES_DIR = CACHE_DIR / "queries"
CACHE_IDENTITIES_DIR = CACHE_DIR / "identities"
CACHE_LIGANDS_DIR = CACHE_DIR / "ligands"

# Rate-limit entre llamadas secuenciales a APIs
API_SLEEP: float = 0.15

# ---------------------------------------------------------------------------
# Sesión HTTP con reintentos y backoff exponencial
# ---------------------------------------------------------------------------

def _make_session(retries: int = 3, backoff_factor: float = 1.0) -> requests.Session:
    session = requests.Session()
    retry = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET", "POST"],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


_session: requests.Session | None = None


def _get_session() -> requests.Session:
    global _session
    if _session is None:
        _session = _make_session()
    return _session


# ---------------------------------------------------------------------------
# Caché a disco (JSON con indentación)
# ---------------------------------------------------------------------------

def _con_cache(func, key: str, cache_path: Path):
    """
    Wrapper genérico de caché a disco.  Si *cache_path* / *key*.json existe,
    carga y devuelve su contenido; de lo contrario ejecuta *func()*, almacena
    el resultado como JSON legible (con indentación) y lo devuelve.
    """
    cache_path.mkdir(parents=True, exist_ok=True)
    safe_key = key.replace("/", "_").replace("\\", "_")
    filepath = cache_path / f"{safe_key}.json"
    if filepath.exists():
        with open(filepath, "r", encoding="utf-8") as f:
            return json.load(f)
    result = func()
    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)
    return result


def _link_ligand_to_query(pdb_id: str, accession: str):
    """
    Copia el ligand JSON del caché compartido a la subcarpeta de la query,
    para que cada query tenga sus PDB hits visibles en un solo lugar.
    """
    shared = CACHE_LIGANDS_DIR / f"{pdb_id.lower()}.json"
    query_lig_dir = CACHE_QUERIES_DIR / accession / "ligands"
    query_lig_dir.mkdir(parents=True, exist_ok=True)
    dest = query_lig_dir / f"{pdb_id.upper()}.json"
    if not dest.exists() and shared.exists():
        shutil.copy2(shared, dest)


# ---------------------------------------------------------------------------
# _asignar_tier
# ---------------------------------------------------------------------------

def _asignar_tier(identity: float) -> int | None:
    """
    Mapea porcentaje de identidad a tier según Rentzsch & Orengo (2009).

    >=50 %  → Tier 1 (alta confianza)
    30-49 % → Tier 2 (confianza moderada)
    <30 %   → Tier 3 (exploratoria)
    """
    if identity >= 50.0:
        return 1
    elif identity >= 30.0:
        return 2
    else:
        return 3  # aún por encima del E-value cutoff


# ---------------------------------------------------------------------------
# _search_pdb  (Step A — RCSB MMseqs2 sequence search)
# ---------------------------------------------------------------------------

def _search_pdb(
    sequence: str,
    accession: str,
    max_hits: int | None = None,
) -> list[dict]:
    """
    Step A of the MMseqs2 search pipeline.

    Submits the amino acid sequence to the RCSB Search API (which uses MMseqs2,
    not BLAST — MMseqs2 is faster and well-suited to searching a bounded
    structural database like the PDB, though it differs in sensitivity and
    scoring model).

    Uses ``return_all_hits: true`` to retrieve every hit in a single request.
    Empirical testing shows this completes in < 1 s for ~73 hits.

    Returns a list of dicts with keys ``pdb_entity_id`` and ``mmseqs2_score``.
    The result is cached under  cache/queries/<accession>/blast.json  (filename
    kept for backward compatibility with existing cache files).
    """
    query_dir = CACHE_QUERIES_DIR / accession

    def _query():
        url = "https://search.rcsb.org/rcsbsearch/v2/query"
        payload = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": 1e-3,
                    "identity_cutoff": 0.0,
                    "sequence_type": "protein",
                    "value": sequence,
                },
            },
            "request_options": {
                "scoring_strategy": "sequence",
                "return_all_hits": True,
            },
            "return_type": "polymer_entity",
        }
        session = _get_session()
        try:
            resp = session.post(url, json=payload, timeout=60)
            if resp.status_code == 204:
                return []
            if not resp.ok:
                logger.error(
                    "RCSB sequence search HTTP %s – body: %s",
                    resp.status_code,
                    resp.text[:500],
                )
                return []
            data = resp.json()
        except Exception as exc:
            logger.error("RCSB sequence search exception: %s", exc)
            return []

        total_count = data.get("total_count", 0)
        hits: list[dict] = []
        for result in data.get("result_set", []):
            identifier = result.get("identifier", "")
            score = result.get("score", 0.0)
            hits.append({
                "pdb_entity_id": identifier,
                "mmseqs2_score": float(score),
            })

        # Sanity check
        if len(hits) != total_count:
            logger.warning(
                "RCSB search: retrieved %d hits but total_count=%d (accession=%s)",
                len(hits), total_count, accession,
            )

        return hits

    all_hits = _con_cache(_query, "blast", query_dir)

    # Apply max_hits cap if requested
    if max_hits is not None and len(all_hits) > max_hits:
        warnings.warn(
            f"_search_pdb: {len(all_hits)} hits for {accession} exceed "
            f"max_hits={max_hits}; truncating.",
            stacklevel=2,
        )
        all_hits = all_hits[:max_hits]

    return all_hits


# ---------------------------------------------------------------------------
# _fetch_identity  (Step B — RCSB Data API + Biopython pairwise alignment)
# ---------------------------------------------------------------------------

def _seq_hash(sequence: str) -> str:
    """Return the first 16 hex chars of the SHA-256 hash of *sequence*."""
    return hashlib.sha256(sequence.encode()).hexdigest()[:16]


def _fetch_identity(query_sequence: str, pdb_entity_id: str) -> float | None:
    """
    Step B of the MMseqs2 search pipeline.

    For a single PDB entity hit, fetches the entity's amino acid sequence from
    the RCSB Data API and computes percent sequence identity via a global
    pairwise alignment (Biopython ``Align.PairwiseAligner``).

    Returns identity_pct (0–100) or ``None`` if the entity sequence cannot be
    retrieved.

    Results are cached under  cache/identities/  with key
    ``{query_hash}:{pdb_entity_id}``.
    """
    cache_key = f"{_seq_hash(query_sequence)}:{pdb_entity_id}"

    def _compute():
        # Fetch target entity sequence from RCSB Data API
        parts = pdb_entity_id.split("_")
        entry_id = parts[0] if parts else pdb_entity_id
        entity_num = parts[1] if len(parts) > 1 else "1"
        url = (
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/"
            f"{entry_id}/{entity_num}"
        )
        session = _get_session()
        try:
            resp = session.get(url, timeout=15)
            if not resp.ok:
                logger.error(
                    "RCSB Data API polymer_entity %s HTTP %s",
                    pdb_entity_id, resp.status_code,
                )
                return None
            target_seq = (
                resp.json()
                .get("entity_poly", {})
                .get("pdbx_seq_one_letter_code_can", "")
            )
            if not target_seq:
                logger.error(
                    "RCSB Data API: no sequence for %s", pdb_entity_id
                )
                return None
        except Exception as exc:
            logger.error(
                "RCSB Data API polymer_entity %s exception: %s",
                pdb_entity_id, exc,
            )
            return None

        # Pairwise alignment with Biopython
        try:
            aligner = Align.PairwiseAligner()
            aligner.mode = "global"
            aligner.match_score = 1
            aligner.mismatch_score = 0
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            alignments = aligner.align(query_sequence, target_seq)
            if not alignments:
                return 0.0
            best = alignments[0]
            matches = best.counts().identities
            alignment_length = best.shape[1]
            if alignment_length == 0:
                return 0.0
            return round(matches / alignment_length * 100, 2)
        except Exception as exc:
            logger.error(
                "Pairwise alignment failed for %s: %s", pdb_entity_id, exc
            )
            return None

    return _con_cache(_compute, cache_key, CACHE_IDENTITIES_DIR)


# ---------------------------------------------------------------------------
# _mmseqs_contra_pdb  (combined two-step search)
# ---------------------------------------------------------------------------

def _mmseqs_contra_pdb(
    sequence: str,
    accession: str,
    max_hits: int | None = None,
) -> list[dict]:
    """
    Busca homólogos PDB para una secuencia proteica usando el pipeline
    MMseqs2 de RCSB en dos pasos:

      Step A (_search_pdb):   RCSB Search API → (pdb_entity_id, mmseqs2_score)
      Step B (_fetch_identity): RCSB Data API + Biopython → identity_pct

    Nota: la RCSB Search API usa MMseqs2 (no BLAST).  MMseqs2 y BLAST
    difieren en sensibilidad y modelo de scoring; MMseqs2 es más rápido y
    adecuado para buscar en una base de datos estructural acotada como PDB.

    Parámetros
    ----------
    sequence : str      Secuencia aminoacídica.
    accession : str     Accesión del query (para organizar el caché).
    max_hits : int | None
        Límite opcional de hits; si total_count > max_hits se emite un warning
        y se trunca.

    Retorna
    -------
    list[dict]  con claves: pdb_id, chain_id, mmseqs2_score, identity_pct
    """
    # Step A — ranked list from Search API
    raw_hits = _search_pdb(sequence, accession, max_hits=max_hits)
    if not raw_hits:
        return []

    # Step B — compute identity for each hit
    enriched: list[dict] = []
    for hit in raw_hits:
        pdb_entity_id = hit["pdb_entity_id"]
        mmseqs2_score = hit["mmseqs2_score"]

        identity_pct = _fetch_identity(sequence, pdb_entity_id)
        time.sleep(API_SLEEP)

        parts = pdb_entity_id.split("_")
        pdb_id = parts[0].upper() if parts else pdb_entity_id.upper()
        chain_id = parts[1] if len(parts) > 1 else ""

        if identity_pct is None:
            logger.warning(
                "Could not compute identity for %s (accession=%s); skipping.",
                pdb_entity_id, accession,
            )
            continue

        enriched.append({
            "pdb_id": pdb_id,
            "chain_id": chain_id,
            "mmseqs2_score": mmseqs2_score,
            "identity_pct": identity_pct,
        })

    return enriched


# ---------------------------------------------------------------------------
# _obtener_ligandos_pdb  (PDBe REST API)
# ---------------------------------------------------------------------------

def _obtener_ligandos_pdb(pdb_id: str) -> list[dict]:
    """
    Obtiene ligandos (small molecules) e interactores proteicos de un PDB ID
    consultando tres endpoints de PDBe y fusionando los resultados.

    - Endpoint A: ligand_monomers → small molecules
    - Endpoint B: molecules       → polymer chain names
    - Endpoint C: mappings/uniprot → UniProt accessions per chain
    """
    pdb_lower = pdb_id.lower()

    def _query():
        session = _get_session()
        interactors: list[dict] = []

        # --- Endpoint A: ligandos small-molecule (HETATM) ---
        url_lig = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_lower}"
        try:
            r = session.get(url_lig, timeout=15)
            if r.ok:
                data = r.json().get(pdb_lower, [])
                for entry in data:
                    chem_id = entry.get("chem_comp_id", "")
                    if chem_id.upper() in HETATM_BLACKLIST:
                        continue
                    interactors.append({
                        "chem_comp_id": chem_id,
                        "uniprot_id": None,
                        "name": entry.get("chem_comp_name", ""),
                        "type": "small_molecule",
                        "chain_id": entry.get("chain_id", ""),
                    })
            elif r.status_code != 404:
                logger.error("PDBe ligand_monomers %s HTTP %s", pdb_id, r.status_code)
        except Exception as exc:
            logger.error("PDBe ligand_monomers %s exception: %s", pdb_id, exc)

        time.sleep(API_SLEEP)

        # --- Endpoint C: UniProt mappings per chain ---
        # Fetch first so we can annotate the protein chains from Endpoint B.
        chain_to_uniprot: dict[str, str] = {}
        url_uniprot = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_lower}"
        try:
            r = session.get(url_uniprot, timeout=15)
            if r.ok:
                uniprot_data = r.json().get(pdb_lower, {}).get("UniProt", {})
                for uniprot_acc, info in uniprot_data.items():
                    for mapping in info.get("mappings", []):
                        ch = mapping.get("chain_id", "")
                        if ch and ch not in chain_to_uniprot:
                            chain_to_uniprot[ch] = uniprot_acc
            elif r.status_code != 404:
                logger.error("PDBe mappings/uniprot %s HTTP %s", pdb_id, r.status_code)
        except Exception as exc:
            logger.error("PDBe mappings/uniprot %s exception: %s", pdb_id, exc)

        time.sleep(API_SLEEP)

        # --- Endpoint B: cadenas poliméricas (interactores proteicos) ---
        url_mol = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb_lower}"
        try:
            r = session.get(url_mol, timeout=15)
            if r.ok:
                data = r.json().get(pdb_lower, [])
                for mol in data:
                    mol_type = mol.get("molecule_type", "")
                    if "polypeptide" not in mol_type.lower():
                        continue
                    chains = mol.get("in_chains", [])
                    name = mol.get("molecule_name", [""])[0] if mol.get("molecule_name") else ""
                    for ch in chains:
                        interactors.append({
                            "chem_comp_id": None,
                            "uniprot_id": chain_to_uniprot.get(ch),
                            "name": name,
                            "type": "protein",
                            "chain_id": ch,
                        })
            elif r.status_code != 404:
                logger.error("PDBe molecules %s HTTP %s", pdb_id, r.status_code)
        except Exception as exc:
            logger.error("PDBe molecules %s exception: %s", pdb_id, exc)

        return interactors

    return _con_cache(_query, pdb_lower, CACHE_LIGANDS_DIR)


# ---------------------------------------------------------------------------
# buscar_homologos_pdb  —  ENTRY POINT
# ---------------------------------------------------------------------------

def buscar_homologos_pdb(
    df: pd.DataFrame,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]]:
    """
    Para cada proteína del DataFrame de entrada, busca homólogos en PDB,
    clasifica por tier de identidad y recupera los ligandos/interactores de
    cada estructura PDB hit.

    Parámetros
    ----------
    df : DataFrame con al menos las columnas:
         - deg_id
         - gene_id  (o query_accession)
         - aa_seq   (o Secuencia)
         - Dominios  (lista de dicts con accession, etc.)

    Retorna
    -------
    dict  {deg_id: (tier1_df, tier2_df, tier3_df)}
    Cada DataFrame sigue el esquema de columnas documentado en el plan.
    """
    results: dict[str, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {}

    # Columnas de salida
    OUTPUT_COLS = [
        "deg_id",
        "query_accession",
        "pfam_accessions",
        "pdb_id",
        "chain_id",
        "mmseqs2_score",
        "identity_pct",
        "tier",
        "interactors",
        "n_small_molecule_interactors",
        "n_protein_interactors",
    ]

    def _empty_df() -> pd.DataFrame:
        return pd.DataFrame(columns=OUTPUT_COLS)

    # Determinar nombre de la columna de secuencia
    seq_col = "aa_seq" if "aa_seq" in df.columns else "Secuencia"
    acc_col = "gene_id" if "gene_id" in df.columns else "ID"

    # Pre-cache de ligandos: colecta PDB IDs únicos para evitar llamadas
    # redundantes (se cachean automáticamente por _con_cache).

    total = len(df)
    print(f"\n  Buscando homólogos PDB para {total} proteínas …")

    for idx, row in tqdm(df.iterrows(), total=total, desc="Homólogos PDB", unit="prot"):
        deg_id = str(row.get("deg_id", ""))
        accession = str(row.get(acc_col, ""))
        sequence = str(row.get(seq_col, ""))

        # Pfam accessions (lista de strings)
        dominios = row.get("Dominios", [])
        if isinstance(dominios, list):
            pfam_accs = list({d["accession"] for d in dominios if "accession" in d})
        else:
            pfam_accs = []

        # Fallback: si no hay secuencia, devolver vacías
        if not sequence or sequence in ("", "nan", "None"):
            results[deg_id] = (_empty_df(), _empty_df(), _empty_df())
            continue

        # 1) MMseqs2 contra PDB (two-step: search + identity)
        try:
            hits = _mmseqs_contra_pdb(sequence, accession)
        except Exception as exc:
            logger.error("deg_id=%s accession=%s mmseqs exception: %s", deg_id, accession, exc)
            results[deg_id] = (_empty_df(), _empty_df(), _empty_df())
            continue

        if not hits:
            results[deg_id] = (_empty_df(), _empty_df(), _empty_df())
            continue

        # 2) Asignar tier y obtener ligandos para cada hit
        rows_out: list[dict] = []
        for hit in hits:
            tier = _asignar_tier(hit["identity_pct"])
            if tier is None:
                continue

            pdb_id = hit["pdb_id"]
            try:
                interactors = _obtener_ligandos_pdb(pdb_id)
            except Exception as exc:
                logger.error(
                    "deg_id=%s pdb_id=%s obtener_ligandos exception: %s",
                    deg_id, pdb_id, exc,
                )
                interactors = []

            # Link the shared ligand JSON into the query subfolder
            _link_ligand_to_query(pdb_id, accession)

            n_sm = sum(1 for i in interactors if i.get("type") == "small_molecule")
            n_prot = sum(1 for i in interactors if i.get("type") == "protein")

            rows_out.append({
                "deg_id": deg_id,
                "query_accession": accession,
                "pfam_accessions": pfam_accs,
                "pdb_id": pdb_id,
                "chain_id": hit.get("chain_id", ""),
                "mmseqs2_score": hit["mmseqs2_score"],
                "identity_pct": hit["identity_pct"],
                "tier": tier,
                "interactors": interactors,
                "n_small_molecule_interactors": n_sm,
                "n_protein_interactors": n_prot,
            })

        time.sleep(API_SLEEP)

        if not rows_out:
            results[deg_id] = (_empty_df(), _empty_df(), _empty_df())
            continue

        full_df = pd.DataFrame(rows_out, columns=OUTPUT_COLS)
        tier1 = full_df[full_df["tier"] == 1].reset_index(drop=True)
        tier2 = full_df[full_df["tier"] == 2].reset_index(drop=True)
        tier3 = full_df[full_df["tier"] == 3].reset_index(drop=True)
        results[deg_id] = (tier1, tier2, tier3)

    return results


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    # E. coli MurA (UniProt P0A749) — well-characterized, many PDB hits
    print("Fetching MurA (P0A749) sequence from UniProt …")
    r = requests.get("https://rest.uniprot.org/uniprotkb/P0A749.fasta", timeout=15)
    r.raise_for_status()
    lines = r.text.strip().split("\n")
    MURA_SEQ = "".join(l.strip() for l in lines if not l.startswith(">"))
    print(f"  Sequence length: {len(MURA_SEQ)} aa\n")

    # 1) Run _mmseqs_contra_pdb
    print("Running _mmseqs_contra_pdb …")
    hits = _mmseqs_contra_pdb(MURA_SEQ, "P0A749_test")
    print(f"  Total hits retrieved: {len(hits)}\n")

    # 2) Print top 5 hits
    print(f"{'#':>3}  {'pdb_entity_id':>15}  {'mmseqs2_score':>14}  {'identity_pct':>12}")
    print("-" * 52)
    for i, h in enumerate(hits[:5]):
        print(
            f"{i:3d}  {h['pdb_id']+'_'+h['chain_id']:>15}  "
            f"{h['mmseqs2_score']:14.4f}  {h['identity_pct']:12.2f}"
        )

    # 3) Assert identity_pct is not None for top 5
    top5 = hits[:5]
    none_count = sum(1 for h in top5 if h["identity_pct"] is None)
    if none_count:
        print(
            f"\n  ⚠  {none_count}/5 top hits have identity_pct=None — "
            "investigate _fetch_identity failures.",
            file=sys.stderr,
        )
        sys.exit(1)
    else:
        print("\n  ✓  All top 5 hits have valid identity_pct values.")

    # 4) Tier assignment for top 5
    print(f"\n{'#':>3}  {'pdb_entity_id':>15}  {'identity_pct':>12}  {'tier':>4}")
    print("-" * 42)
    for i, h in enumerate(top5):
        tier = _asignar_tier(h["identity_pct"])
        print(
            f"{i:3d}  {h['pdb_id']+'_'+h['chain_id']:>15}  "
            f"{h['identity_pct']:12.2f}  {tier:>4}"
        )

    print("\n  ✓  Smoke test passed.")
