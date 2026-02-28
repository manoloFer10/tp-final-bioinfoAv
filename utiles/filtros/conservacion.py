import pandas as pd
import argparse
import subprocess
import os
import zipfile
import shutil

# Accessions GCA de los 15 proteomas de Burkholderia cenocepacia
GCA_ACCESSIONS = [
    "GCA_001718895.1", "GCA_001999785.1", "GCA_001606135.1",
    "GCA_001999825.1", "GCA_039728355.1", "GCA_021535065.1",
    "GCA_018223805.1", "GCA_025643575.1", "GCA_014357995.1",
    "GCA_018228625.1", "GCA_018223765.1", "GCA_023374855.1",
    "GCA_023374835.1", "GCA_023374335.1", "GCA_001606115.1",
]

PROTEOMAS_DIR = "data/conservacion/15proteomas"
DB_PATH       = "data/conservacion/otros_db"
QUERY_FAA     = "data/protein.faa"


def _descargar_proteomas(proteomas_dir: str = PROTEOMAS_DIR):
    """Descarga los 15 proteomas de NCBI si no existen ya."""
    os.makedirs(proteomas_dir, exist_ok=True)
    if len([f for f in os.listdir(proteomas_dir) if f.endswith('.faa')]) >= len(GCA_ACCESSIONS):
        print("  Proteomas ya descargados. Salteando.")
        return

    print("  Descargando proteomas de NCBI...")
    zip_path = os.path.join(proteomas_dir, "ncbi_dataset.zip")
    tmp_dir  = os.path.join(proteomas_dir, "ncbi_tmp")

    subprocess.run([
        "datasets", "download", "genome", "accession",
        *GCA_ACCESSIONS,
        "--include", "protein",
        "--filename", zip_path
    ], check=True)

    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(tmp_dir)

    for root, dirs, files in os.walk(tmp_dir):
        for file in files:
            if file == "protein.faa":
                cepa = os.path.basename(root)
                shutil.copy(os.path.join(root, file), os.path.join(proteomas_dir, f"{cepa}.faa"))

    os.remove(zip_path)
    shutil.rmtree(tmp_dir)
    print(f"  Descarga completa: {len(os.listdir(proteomas_dir))} proteomas.")


def _armar_db(proteomas_dir: str = PROTEOMAS_DIR, db_path: str = DB_PATH):
    """Tagguea los .faa por cepa y crea la DB BLAST si no existe."""
    if os.path.exists(db_path + ".phr"):
        print("  DB BLAST ya existe. Salteando.")
        return

    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    tagged_dir      = os.path.join(os.path.dirname(db_path), "tagged")
    tagged_combined = os.path.join(os.path.dirname(db_path), "proteomas_tagged.faa")
    os.makedirs(tagged_dir, exist_ok=True)

    faa_files = [f for f in os.listdir(proteomas_dir) if f.endswith(".faa")]
    if not faa_files:
        raise FileNotFoundError(f"No hay archivos .faa en '{proteomas_dir}'.")

    print(f"  Armando DB BLAST con {len(faa_files)} cepas...")
    with open(tagged_combined, "w") as out_faa:
        for faa in faa_files:
            tag         = faa.replace(".faa", "")
            tagged_file = os.path.join(tagged_dir, f"{tag}.tagged.faa")
            with open(os.path.join(proteomas_dir, faa)) as fh, open(tagged_file, "w") as tf:
                for line in fh:
                    tf.write(f">{tag}|{line[1:]}" if line.startswith(">") else line)
            with open(tagged_file) as tf:
                out_faa.write(tf.read())

    subprocess.run([
        "makeblastdb", "-in", tagged_combined,
        "-dbtype", "prot", "-out", db_path
    ], check=True)
    print(f"  DB creada en {db_path}.")


def filtrar_conservacion(proteinas: list[str]) -> list[str]:
    """
    Dado una lista de IDs de proteínas, descarga los 15 proteomas de Burkholderia cenocepacia
    de NCBI (si no existen), arma la DB BLAST (si no existe), corre BLASTp y retorna
    las IDs conservadas en >= 14/15 cepas.

    Criterios de conservación:
    - Identidad >= 50%
    - Cobertura (length/qlen) >= 0.70
    - E-value <= 1e-10
    - Presente en >= 14 de 15 cepas
    """
    output_dir = "data/conservacion/outputs"
    os.makedirs(output_dir, exist_ok=True)

    blast_out      = os.path.join(output_dir, "blast_out.tsv")
    conteo_path    = os.path.join(output_dir, "conteo_por_proteina.tsv")
    presencia_path = os.path.join(output_dir, "presencia_unica.tsv")

    _descargar_proteomas(PROTEOMAS_DIR)
    _armar_db(PROTEOMAS_DIR, DB_PATH)

    print("  Corriendo BLASTp...")
    subprocess.run([
        "blastp", "-query", QUERY_FAA,
        "-db", DB_PATH,
        "-evalue", "1e-10",
        "-num_threads", "10",
        "-outfmt", "6 qseqid sseqid pident length qlen evalue bitscore",
        "-out", blast_out
    ], check=True)

    cols     = ["qseqid", "sseqid", "pident", "length", "qlen", "evalue", "bitscore"]
    df_blast = pd.read_csv(blast_out, sep="\t", names=cols)

    df_blast = df_blast[
        (df_blast["pident"] >= 50.0) &
        (df_blast["length"] / df_blast["qlen"] >= 0.70) &
        (df_blast["evalue"] <= 1e-10)
    ]

    df_blast["cepa"] = df_blast["sseqid"].str.split("|").str[0]
    presencia = df_blast[["qseqid", "cepa"]].drop_duplicates()
    presencia.to_csv(presencia_path, sep="\t", index=False, header=False)

    conteo_df = presencia.groupby("qseqid")["cepa"].nunique().reset_index()
    conteo_df.columns = ["qseqid", "n_cepas"]
    conteo_df = conteo_df.sort_values("n_cepas", ascending=False)
    conteo_df.to_csv(conteo_path, sep="\t", index=False)

    conservadas = conteo_df[conteo_df["n_cepas"] >= 14]["qseqid"].tolist()

    print(f"  Hits tras filtros: {len(df_blast)}")
    print(f"  Proteínas únicas con hits: {len(conteo_df)}")
    print(f"  Conservadas (>=14/15 cepas): {len(conservadas)}")

    return [p for p in conservadas if p in proteinas]


def procesar_conservacion(df):
    """
    Dado un dataframe con columna 'ID',
    devuelve el dataframe filtrado a las proteínas conservadas en otras cepas.
    """
    df = df[df["ID"].isin(filtrar_conservacion(df["ID"].tolist()))]
    return df


if __name__ == "__main__":
    # se llama así: python filtros/conservacion.py proteinas_a_testear.csv
    # (siempre y cuando estén en el root del proyecto)
    parser = argparse.ArgumentParser(description="Filtrar proteínas por conservación en otras cepas.")
    parser.add_argument("archivo_testeo", help="path al CSV con las proteínas a testear.")
    args = parser.parse_args()

    if args.archivo_testeo:
        df = pd.read_csv(args.archivo_testeo)
        print("="*50)
        print(f"Se leyeron {len(df)} proteínas del archivo {args.archivo_testeo}.")
        df = procesar_conservacion(df)
        print("="*50)
        print(f"Se filtraron {len(df)} proteínas:")
        print(df["ID"].tolist())
    else:
        print("No se proporcionó un archivo de testeos. Por favor, use archivo_testeo para especificar el archivo CSV con las proteínas a testear.")
        print(f'Se filtraron {len(df)} proteínas:')
        print(df['ID'].tolist())
    else:
        print("No se proporcionó un archivo de testeos. Por favor, use --archivo_testeo para especificar el archivo CSV con las proteínas a testear.")
