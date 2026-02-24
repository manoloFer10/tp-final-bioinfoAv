# existe https://academic.oup.com/nar/article/49/D1/D677/5937083, 
# una base de datos de genes esenciales que parece que es bastante usada en la literatura.
# la extraccion de los datos de cenocepacia está en data/tubic_esenciales_bacterias.

import pandas as pd 
import argparse

def filtrar_esencialidad(proteinas: list[str]):
    # esta función debería filtrar las proteínas por si son esenciales o no...
    # debería tomar los ids de entrada y devolver los que están en la base de datos de esenciales. 
    filtradas = []
    return filtradas

def procesar_esencialidad(df):
    df = df[df['ID'].isin(filtrar_esencialidad(df['ID'].tolist()))]
    return df

if __name__ == "__main__":
    # esto se usa para ejecutar el script desde la línea de comandos, pasando un archivo CSV con las proteínas a testear.
    # el CSV debería tener una columna llamada 'Proteinas' con los nombres de las proteínas a testear.
    # se llama así: python filtros/esencialidad.py --archivo_testeo proteinas_a_testear.csv
    # (siempre y cuando estén en el root del proyecto)
    parser = argparse.ArgumentParser(description="Filtrar proteínas por esencialidad.")
    parser.add_argument("archivo_testeo", help="path al CSV con las proteínas a testear.")
    args = parser.parse_args()
    if args.archivo_testeo:
        df = pd.read_csv(args.archivo_testeo)
        print('='*50)
        print(f'Se leyeron {len(df)} proteínas del archivo {args.archivo_testeo}.')
        df = procesar_esencialidad(df)
        print('='*50)
        print(f'Se filtraron {len(df)} proteínas:')
        print(df['ID'].tolist())
    else:
        print("No se proporcionó un archivo de testeos. Usar --archivo_testeo para especificar el archivo CSV con las proteínas a testear.")