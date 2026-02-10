import pandas as pd
import argparse

def filtrar_conservacion(proteinas: list[str]):
    # esta función debería filtrar las proteínas por su nivel de conservación en otras cepas...
    filtradas = []
    return filtradas

if __name__ == "__main__":
    # esto se usa para ejecutar el script desde la línea de comandos, pasando un archivo CSV con las proteínas a testear.
    # el CSV debería tener una columna llamada 'Proteinas' con los nombres de las proteínas a testear.
    # se llama así: python filtros/conservacion.py --archivo_testeo proteinas_a_testear.csv
    # (siempre y cuando estén en el root del proyecto)
    parser = argparse.ArgumentParser(description="Filtrar proteínas por conservación en otras cepas.")
    parser.add_argument("archivo_testeo", help="path al CSV con las proteínas a testear.")
    args = parser.parse_args()
    if args.archivo_testeo:
        df = pd.read_csv(args.archivo_testeo)
        proteinas = list(df['Proteinas'])
        print('='*50)
        print(f'Se leyeron {len(proteinas)} proteínas del archivo {args.archivo_testeo}.')
        filtradas = filtrar_conservacion(proteinas)
        print('='*50)
        print(f'Se filtraron {len(filtradas)} proteínas:')
        print(filtradas)
    else:
        print("No se proporcionó un archivo de testeos. Por favor, use --archivo_testeo para especificar el archivo CSV con las proteínas a testear.")