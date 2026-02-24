import argparse
import pickle
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from utiles.filtros.esencialidad import procesar_esencialidad
from utiles.filtros.conservacion import procesar_conservacion
from utiles.filtros.homologia_humanos import procesar_homologia_humanos
from utiles.dominios import obtener_dominios_df
from utiles.busqueda_ligandos import buscar_homologos_pdb


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
    
    # ----- FILTRADO DE PROTEÍNAS DE INTERÉS -----
    print('='*20, ' FILTRADO ', '='*20)
    df = procesar_esencialidad(df)
    print(' Esencialidad: ',f'Hay {len(df)} proteínas luego del filtro de esencialidad.')
    #df = procesar_conservacion(df)
    print(' Conservación: ',f'Hay {len(df)} proteínas luego del filtro de conservación.')
    #df = procesar_homologia_humanos(df)
    print(' Homología Humanos: ',f'Hay {len(df)} proteínas luego del filtro de homología con humanos.')

    # ----- BÚSQUEDA DE DOMINIOS CONSERVADOS Y ARMADO MSA -----
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

    # ----- BÚSQUEDA DE HOMÓLOGOS PDB Y LIGANDOS -----
    print("\n", '='*20, ' HOMÓLOGOS PDB Y LIGANDOS ', '='*20)
    homologos = buscar_homologos_pdb(df)

    # Persistir resultados (pickle por la estructura anidada de DataFrames)
    output_path = Path("output")
    output_path.mkdir(exist_ok=True)
    with open(output_path / "homologos_pdb.pkl", "wb") as f:
        pickle.dump(homologos, f)
    print(f"  Resultados guardados en {output_path / 'homologos_pdb.pkl'}")

    # Resumen por proteína
    print("\n  Resumen de homólogos PDB por proteína:")
    print(f"  {'deg_id':<20} {'T1':>4} {'T2':>4} {'T3':>4}  {'T1+lig':>6} {'T2+lig':>6} {'T3+lig':>6}")
    print("  " + "-"*60)
    for deg_id, (t1, t2, t3) in homologos.items():
        def _con_lig(tier_df):
            if tier_df.empty:
                return 0
            return int((tier_df["n_small_molecule_interactors"] + tier_df["n_protein_interactors"] > 0).sum())
        n1, n2, n3 = len(t1), len(t2), len(t3)
        l1, l2, l3 = _con_lig(t1), _con_lig(t2), _con_lig(t3)
        if n1 + n2 + n3 > 0:
            print(f"  {deg_id:<20} {n1:>4} {n2:>4} {n3:>4}  {l1:>6} {l2:>6} {l3:>6}")

    total_con_hits = sum(1 for (t1, t2, t3) in homologos.values() if len(t1) + len(t2) + len(t3) > 0)
    print(f"\n  Proteínas con al menos un homólogo PDB: {total_con_hits}/{len(homologos)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline de filtrado de proteínas.")
    parser.add_argument("--archivo_proteinas", required=False, default="data\genoma_burkholderia_cenocepacia\protein.faa", help="path al archivo FASTA con las proteínas a filtrar.")
    args = parser.parse_args()
    main(args)