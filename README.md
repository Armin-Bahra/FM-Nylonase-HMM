# TFM-Nylonase-HMM
FM-Nylonase-HMM

Repositorio del Trabajo Fin de Máster «Diseño de un modelo oculto de Markov para la identificación de nuevas nylonasas».

Autor: Armin Bahraminia De La Fuente
Tutora: Tania Díaz
Máster Universitario en Bioinformática, Universidad Europea de Madrid (octubre de 2026)
Resumen

Se construye un perfil HMM (HMMER 3.4) a partir de cinco secuencias NylC curadas, se valida con controles positivo y negativo y se usa para buscar posibles candidatas a NylC en UniProtKB mediante el servidor HMMER-EBI. Las candidatas se seleccionan por similitud de secuencia con el perfil (umbral de 88,1 bits) y no constituyen una demostración de actividad nylonasa.

Estructura del repositorio
FM-Nylonase-HMM/
├── README.md
├── data/        Datos de entrada (secuencias, alineamientos, matrices, controles, tablas de HMMER-EBI)
├── results/     Resultados (perfil HMM, tablas, figuras)
└── scripts/     Scripts de análisis (01 a 05)

Los scripts usan rutas relativas (data/ y results/), por lo que deben ejecutarse desde la raíz del repositorio.

Scripts (orden de uso)
Script	Apartado del TFM	Qué hace
scripts/02_heatmaps_identidad.R	4.3 y 4.4	Dibuja los heatmaps de identidad de secuencia a partir de las matrices .pim de Clustal Omega (14 secuencias con 0 y 5 iteraciones, y 8 NylC).
scripts/05_arbol_agrupamiento.py	4.3 y 4.4	Dibuja los árboles de agrupamiento de Clustal Omega (.phylotree). No es un análisis filogenético.
scripts/03_triada.py	4.4	Comprueba la conservación del sitio de autocorte (Asn266/Thr267) y de la tríada catalítica (Thr267, Asp306, Asp308) en las NylC, tomando Q79F77 como referencia de numeración.
scripts/01_hmm_construccion_y_validacion.R	4.5	Construye el perfil HMM (hmmbuild), lo evalúa contra los controles positivo y negativo (hmmsearch, en local mediante WSL) y calcula el margen para fijar el umbral.
scripts/04_analisis_candidatos.py	4.5	Procesa la tabla de HMMER-EBI, filtra por el umbral (88,1 bits), clasifica la anotación de UniProt, cuenta especies y géneros y genera tablas y figuras.
Ejemplos de uso
bash
# Heatmaps (R)
Rscript scripts/02_heatmaps_identidad.R

# Árbol de agrupamiento (Python)
python scripts/05_arbol_agrupamiento.py data/arbol_14seq_5iter.phylotree "Título" results/arbol_14seq_5iter.png

# Tríada catalítica (Python)
python scripts/03_triada.py data/alineamiento_8nylc.aln-fasta results/triada_NylC.tsv

# Perfil HMM y controles (R; requiere HMMER 3.4 en WSL)
Rscript scripts/01_hmm_construccion_y_validacion.R

# Análisis de candidatas (Python)
python scripts/04_analisis_candidatos.py data/<tabla_hmmer_ebi>.tsv results/candidatas 88.1
Flujo de trabajo
Obtención de secuencias (PlasticDB, UniProt y revisión bibliográfica): 14 proteínas (NylA, NylB y NylC).
Alineamiento múltiple con Clustal Omega y matrices de identidad (heatmaps y árboles).
Alineamiento de las 8 NylC y comprobación de la tríada catalítica.
Construcción del perfil HMM con 5 NylC y validación con 3 NylC de control positivo y 8.364 secuencias de control negativo.
Búsqueda con HMMER-EBI en UniProtKB y análisis descriptivo de las candidatas.

El diagrama completo está en el Anexo 3 de la memoria.

Versiones del software
Software	Versión
HMMER	3.4
Clustal Omega	1.2.4
Python	3.12.10
Biopython	1.88
pandas	3.1.0
matplotlib	3.11.2
R	4.5.2
tidyverse	2.0.0
pheatmap	1.0.13
PyMOL	3.1.7.2
Bases de datos
PlasticDB (https://plasticdb.org/proteins): consulta del 19 de septiembre de 2026.
UniProt versión 2026_03: selección de secuencias y controles negativos (Swiss-Prot, secuencias revisadas).
UniProtKB versión 2025_01: búsqueda final con HMMER-EBI (8 de octubre de 2026).
Limitaciones
El perfil se construyó con pocas secuencias muy similares entre sí (eff_nseq = 0,49).
El control positivo es pequeño (3 secuencias).
Las anotaciones de las candidatas son automáticas y no prueban función.
No hay validación experimental.
Cómo citar

Bahraminia De La Fuente, A. (2026). Diseño de un modelo oculto de Markov para la identificación de nuevas nylonasas [Trabajo Fin de Máster, Universidad Europea de Madrid]. https://github.com/Armin-Bahra/FM-Nylonase-HMM
