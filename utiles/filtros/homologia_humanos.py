import pandas as pd
import argparse
import subprocess
import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Umbrales para considerar homología
EVALUE_UMBRAL = 1e-5
IDENTIDAD_UMBRAL = 30

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DB_PATH = os.path.join(BASE_DIR, "data", "proteoma_humano", "proteoma_humano_db")

def filtrar_homologia_humanos(proteinas: list[str], hits_blast: set):
    """
    Dado una lista de IDs y un set de IDs con homología,
    devuelve solo los IDs que NO tienen homología con humanos.
    """
    return [p for p in proteinas if p not in hits_blast]

def procesar_homologia_humanos(df):
    """
    El objetivo de esta función es, dado un dataset, filtrarlo por las proteínas que definamos tienen homología. 
    En el futuro puede cambiar, tal vez le pasamos un float con el porcentaje de homología, o algo así.
    Por ahora, lo importante es que esta función tome un dataset y devuelva un dataset filtrado por homología.
    """
    # 1. Crear archivo FASTA temporal
    os.makedirs("temp", exist_ok=True)
    fasta_temporal = "temp/proteinas_consulta.fasta"
    records = []
    for _, row in df.iterrows():
        record = SeqRecord(
            Seq(row['Secuencia']),
            id=row['ID'],
            description=""
        )
        records.append(record)
    SeqIO.write(records, fasta_temporal, "fasta")

    # 2. Correr BLAST
    blast_output = "temp/blast_resultados.txt"
    cmd = [
        "blastp",
        "-query", fasta_temporal,
        "-db", DB_PATH,
        "-out", blast_output,
        "-outfmt", "6 qseqid sseqid pident evalue",
        "-evalue", str(EVALUE_UMBRAL),
        "-num_threads", "4"
    ]
    print("Corriendo BLAST, esto puede tardar varios minutos...")
    subprocess.run(cmd, check=True)

    # 3. Leer resultados
    blast_df = pd.read_csv(
        blast_output,
        sep="\t",
        names=["qseqid", "sseqid", "pident", "evalue"]
    )

    # 4. Aplicar umbrales y obtener IDs con homología
    hits_filtrados = blast_df[
        (blast_df['evalue'] < EVALUE_UMBRAL) &
        (blast_df['pident'] > IDENTIDAD_UMBRAL)
    ]
    ids_con_homologia = set(hits_filtrados['qseqid'].unique())

    # 5. Filtrar usando filtrar_homologia_humanos
    ids_sin_homologia = filtrar_homologia_humanos(df['ID'].tolist(), ids_con_homologia)
    df = df[df['ID'].isin(ids_sin_homologia)]

    # 6. Imprimir resumen
    print("="*50)
    print("FILTRO DE HOMOLOGÍA CON HUMANOS")
    print("="*50)
    print(f"Criterios utilizados:")
    print(f"  - E-value < {EVALUE_UMBRAL}")
    print(f"  - Identidad > {IDENTIDAD_UMBRAL}%")
    print(f"Proteínas descartadas por homología: {len(ids_con_homologia)}")
    print(f"Proteínas conservadas: {len(df)}")
    print("="*50)

    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filtrar proteínas por homología con humanos.")
    parser.add_argument("archivo_testeo", help="path al CSV con las proteínas a testear.")
    args = parser.parse_args()
    if args.archivo_testeo:
        df = pd.read_csv(args.archivo_testeo)
        print('='*50)
        print(f'Se leyeron {len(df)} proteínas del archivo {args.archivo_testeo}.')
        df = procesar_homologia_humanos(df)
        print('='*50)
        print(f'Se filtraron {len(df)} proteínas:')
        print(df['ID'].tolist())
    else:
        print("No se proporcionó un archivo de testeos. Por favor, use --archivo_testeo para especificar el archivo CSV con las proteínas a testear.")