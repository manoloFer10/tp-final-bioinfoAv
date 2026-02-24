import argparse
import pandas as pd
from Bio import SeqIO 
from utiles.filtros.esencialidad import procesar_esencialidad
from utiles.filtros.conservacion import procesar_conservacion
from utiles.filtros.homologia_humanos import procesar_homologia_humanos
from utiles.dominios import obtener_dominios_df


def main(args):
    archivo_proteinas = args.archivo_proteinas
    seqs = []
    ids = []
    descripciones = []
    for record in SeqIO.parse(archivo_proteinas, "fasta"):
        seqs.append(str(record.seq))
        ids.append(record.id)
        descripciones.append(record.description)

    df = pd.DataFrame({
        'ID': ids,
        'Descripcion': descripciones,
        'Secuencia': seqs
    })

    print(f'Se leyeron {len(df)} proteínas del archivo {archivo_proteinas}.')
    
    print('='*20, ' Filtrando esencialidad ', '='*20)
    df = procesar_esencialidad(df)
    print(f'Hay {len(df)} proteínas luego del filtro de esencialidad.')
    print('='*20, ' Filtrando conservación ', '='*20)
    df = procesar_conservacion(df)
    print(f'Hay {len(df)} proteínas luego del filtro de conservación.')
    print('='*20, ' Filtrando homología con humanos ', '='*20)
    df = procesar_homologia_humanos(df)
    print(f'Hay {len(df)} proteínas luego del filtro de homología con humanos.')
    print('='*50)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline de filtrado de proteínas.")
    parser.add_argument("archivo_proteinas", help="path al archivo FASTA con las proteínas a filtrar.")
    args = parser.parse_args()
    main(args)