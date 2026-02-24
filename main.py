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
    
    print('='*20, ' FILTRADO ', '='*20)
    df = procesar_esencialidad(df)
    print(' Esencialidad: ',f'Hay {len(df)} proteínas luego del filtro de esencialidad.')
    #df = procesar_conservacion(df)
    print(' Conservación: ',f'Hay {len(df)} proteínas luego del filtro de conservación.')
    #df = procesar_homologia_humanos(df)
    print(' Homología Humanos: ',f'Hay {len(df)} proteínas luego del filtro de homología con humanos.')

    print("\n",'='*20, ' DOMINIOS Y MSA ', '='*20)
    df = obtener_dominios_df(df)
    #mostramos algunas estadísticas básicas sobre los dominios encontrados:
    total_dominios = df['Dominios'].apply(len).sum()
    print(f'Total de dominios PFAM encontrados: {total_dominios}')
    #mostramos los dominios encontrados para 2 proteínas aleatorias:
    print("\nEjemplo de dominios encontrados para 2 proteínas:")
    for i in range(min(2, len(df))):
        print(f"Proteína: {df.iloc[i]['gene_id']}")
        print(f"Descripción: {df.iloc[i]['function']}")
        print("Dominios PFAM encontrados:")
        for dominio in df.iloc[i]['Dominios']:
            print(f"  - Accession: {dominio['accession']}, Nombre: {dominio['name']}, Coordenadas: {dominio['start']}-{dominio['end']}")
        print("-"*40)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline de filtrado de proteínas.")
    parser.add_argument("--archivo_proteinas", required=False, default="data\genoma_burkholderia_cenocepacia\protein.faa", help="path al archivo FASTA con las proteínas a filtrar.")
    args = parser.parse_args()
    main(args)