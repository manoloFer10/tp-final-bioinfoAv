import os
import gzip
import shutil
import subprocess
import urllib.request

# Rutas
DATA_DIR = "data/proteoma_humano"
FASTA_GZ = os.path.join(DATA_DIR, "UP000005640_9606.fasta.gz")
FASTA = os.path.join(DATA_DIR, "UP000005640_9606.fasta")
DB_PATH = os.path.join(DATA_DIR, "proteoma_humano_db")

# URL del proteoma humano en UniProt
URL = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/UP000005640/UP000005640_9606.fasta.gz"

def descargar_proteoma():
    print("Descargando proteoma humano desde UniProt...")
    os.makedirs(DATA_DIR, exist_ok=True)
    urllib.request.urlretrieve(URL, FASTA_GZ)
    print(f"Descargado en {FASTA_GZ}")

def descomprimir_proteoma():
    print("Descomprimiendo proteoma humano...")
    with gzip.open(FASTA_GZ, 'rb') as f_in:
        with open(FASTA, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    print(f"Descomprimido en {FASTA}")

def crear_base_datos():
    print("Creando base de datos BLAST...")
    cmd = [
        "makeblastdb",
        "-in", FASTA,
        "-dbtype", "prot",
        "-out", DB_PATH
    ]
    subprocess.run(cmd, check=True)
    print(f"Base de datos creada en {DB_PATH}")

if __name__ == "__main__":
    print("="*50)
    print("SETUP - TP BIOINFORMÁTICA")
    print("="*50)

    if os.path.exists(DB_PATH + ".pin"):
        print("La base de datos ya existe, saltando setup.")
    else:
        descargar_proteoma()
        descomprimir_proteoma()
        crear_base_datos()
        print("="*50)
        print("Setup completado exitosamente.")
        print("Ya podés correr: python main.py")
        print("="*50)