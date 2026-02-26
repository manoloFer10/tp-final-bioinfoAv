# Filtros de Proteínas - Bioinformática Avanzada

Herramienta de bioinformática para filtrar y analizar proteínas basándose en diferentes criterios: conservación en otras cepas, esencialidad y homología con proteínas humanas.

## Requisitos

- Python 3.11+
- Conda (Miniconda o Anaconda)
- BLAST+ instalado en el sistema. Descargarlo desde:
  https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  Para Windows elegir el archivo `ncbi-blast-...-win64.exe`

## Instalación

### 1. Clonar o descargar el repositorio

```bash
git clone <url-del-repositorio>
cd tp-final-bioinfoAv
```

### 2. Crear el ambiente de Conda

```bash
conda env create -f environment.yml
```

Este comando crea un ambiente llamado `bioinfo` con todas las dependencias necesarias.

### 3. Activar el ambiente

```bash
conda activate bioinfo
```

Debes activar este ambiente cada vez que trabajes con los scripts.

### 4. Setup inicial para que funcione el filtro de homología con humanos
Antes de correr el proyecto por primera vez, ejecutá el script de setup
que descarga el proteoma humano y crea la base de datos BLAST:
```bash
python setup.py
```
Este script:
- Descarga el proteoma humano desde UniProt (~500MB, puede tardar varios minutos)
- Descomprime el archivo
- Crea la base de datos BLAST

## Uso de los Filtros

Todos los scripts de filtrado se encuentran en la carpeta `filtros/` y utilizan archivos CSV como entrada. El archivo CSV debe contener una columna llamada `Proteinas` con los nombres de las proteínas a filtrar.

### Formato del archivo de entrada

El archivo CSV debe tener la siguiente estructura:

```
Proteinas
proteína_1
proteína_2
proteína_3
...
```

### Ejecutar los filtros

#### 1. Filtro de Conservación

Filtra las proteínas según su nivel de conservación en otras cepas bacterianas.

```bash
python filtros/conservacion.py archivo_con_proteinas.csv
```

**Salida:** Lista de proteínas filtradas por su grado de conservación.

#### 2. Filtro de Esencialidad

Filtra las proteínas según su grado de esencialidad para la viabilidad de la bacteria.

```bash
python filtros/esencialidad.py archivo_con_proteinas.csv
```

**Salida:** Lista de proteínas filtradas por esencialidad.

#### 3. Filtro de Homología con Humanos

Filtra las proteínas que tienen homología con proteínas humanas.

```bash
python filtros/homologia_humanos.py archivo_con_proteinas.csv
```

**Salida:** Lista de proteínas con homología identificada en proteomas humanos.

## Estructura del Proyecto

```
tp-final-bioinfoAv/
├── README.md                    # Este archivo
├── environment.yml              # Dependencias del proyecto
├── main.py                      # Script principal (en desarrollo)
└── filtros/
    ├── conservacion.py          # Filtro de conservación
    ├── esencialidad.py          # Filtro de esencialidad
    └── homologia_humanos.py     # Filtro de homología con humanos
```


