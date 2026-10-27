"""
04_analisis_candidatos.py
------------------------------------------------------------------------
QUÉ HACE (apartado 4.5 del TFM, "Análisis de las candidatas")
    Procesa la tabla de resultados de hmmsearch descargada de HMMER-EBI
    (TSV) y hace un análisis DESCRIPTIVO de las candidatas:
      1. Deja una sola fila por proteína (el TSV trae una fila por dominio).
      2. Comprueba que todas superan el umbral de puntuación (88,1 bits).
      3. Clasifica cada proteína según su ANOTACIÓN en UniProtKB (campo
         "Description") mediante expresiones regulares.
      4. Agrupa las proteínas por género, a partir del campo "Species".
      5. Guarda dos tablas CSV y tres figuras.
    NO hace ningún análisis estadístico ni comprueba la tríada catalítica
    en las candidatas.

USO (desde la terminal)
    python 04_analisis_candidatos.py tabla.tsv carpeta_salida [umbral] [etiqueta]

    tabla.tsv         TSV descargado de HMMER-EBI
    carpeta_salida    carpeta donde se guardan tablas y figuras (se crea)
    umbral            (opcional) umbral en bits; por defecto 88.1
    etiqueta          (opcional) texto que se añade al nombre de la tabla de
                      candidatas; por defecto "88.1bits". Para la búsqueda por
                      E-value conviene usar, por ejemplo, "evalue", y una
                      carpeta de salida distinta.

ENTRADA: columnas del TSV de HMMER-EBI que se usan
    Target Name, Target Accession, Target Length, Score (bit-score de la
    secuencia completa), E-value, Domain Index, Description, Species, Tax-Id

SALIDAS (en carpeta_salida)
    candidatos_<etiqueta>.csv         una fila por proteína candidata
    top50_candidatos.csv              las 50 con mayor puntuación
    fig_histograma_puntuaciones.png   distribución de puntuaciones
    fig_categorias_anotacion.png      tipo de proteína según su anotación
    fig_generos_frecuentes.png        20 géneros/grupos más frecuentes

REQUISITOS
    Python 3, pandas y matplotlib (en el TFM: Python 3.13.16, pandas 3.0.5,
    matplotlib 3.11.2).   Instalación:  pip install pandas matplotlib
------------------------------------------------------------------------
"""
import os
import re
import sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")            # modo sin ventana: guarda las figuras directamente en fichero
import matplotlib.pyplot as plt

# --- 1. Argumentos -----------------------------------------------------
if len(sys.argv) < 3:
    sys.exit("Uso: python 04_analisis_candidatos.py tabla.tsv carpeta_salida [umbral] [etiqueta]")
TSV = sys.argv[1]                                   # tabla de HMMER-EBI
OUT = sys.argv[2]                                   # carpeta de salida
UMBRAL = float(sys.argv[3]) if len(sys.argv) > 3 else 88.1   # umbral en bits
ETIQUETA = sys.argv[4] if len(sys.argv) > 4 else "88.1bits"   # sufijo del nombre de la tabla
os.makedirs(OUT, exist_ok=True)                     # crea la carpeta si no existe

# --- 2. Lectura y una fila por proteína --------------------------------
# quotechar='"' porque HMMER-EBI entrecomilla algunos campos de texto.
d = pd.read_csv(TSV, sep="\t", quotechar='"')
n_filas = len(d)

# El TSV tiene una fila por dominio encontrado; una proteína con varios
# dominios aparece en varias filas. "Score" es la puntuación de la secuencia
# completa (igual en todas sus filas). Se ordena por "Domain Index" y se
# conserva la primera fila de cada proteína ("Target Name").
s = d.sort_values("Domain Index").drop_duplicates("Target Name").copy()

# Control: el umbral ya se aplicó en el formulario de HMMER-EBI, por lo que
# este filtro no debería eliminar nada; se deja como comprobación.
n_antes = len(s)
s = s[s["Score"] >= UMBRAL].copy()
print(f"Filas del TSV: {n_filas} | proteínas distintas: {n_antes} | "
      f"con puntuación >= {UMBRAL}: {len(s)}")

# Limpia los campos de texto (quita comillas sobrantes y valores vacíos)
s["Description"] = s["Description"].fillna("").str.strip('"')
s["Species"] = s["Species"].fillna("").str.strip('"')


# --- 3. Clasificación según la anotación de UniProtKB ------------------
# IMPORTANTE: es una clasificación AUTOMÁTICA basada solo en el texto de la
# anotación. No se ha verificado la función real de las proteínas.
def categoria(x):
    """Categoría general (5 categorías) según palabras de la anotación."""
    t = x.lower()
    if re.search(r"6-aminohexanoat|nylon|nylc", t):
        return "Nylonasa / 6-aminohexanoato hidrolasa"
    if re.search(r"peptidase|dmpa|d-esterase|s58|aminopeptidase", t):
        return "Peptidasa (S58, P1, T4, L-aminopeptidasa)"
    if re.search(r"uncharacterized|hypothetical|putative|unannotated|predicted|^$", t):
        return "Sin caracterizar / putativa"
    if "hydrolase" in t:
        return "Hidrolasa (otras)"
    return "Otras"


def tipo(x):
    """Tipo más detallado (familias de peptidasas, etc.)."""
    t = x.lower()
    if re.search(r"6-aminohexanoat|nylon|nylc", t):
        return "6-aminohexanoato hidrolasa"
    if re.search(r"d-esterase|l-aminopeptidase", t):
        return "L-aminopeptidasa / D-esterasa"
    if re.search(r"\bp1\b", t):
        return "Peptidasa familia P1"
    if re.search(r"s58|dmpa", t):
        return "Peptidasa familia S58"
    if re.search(r"\bt4\b", t):
        return "Peptidasa familia T4"
    if re.search(r"peptidase|aminopeptidase", t):
        return "Otras peptidasas"
    if re.search(r"uncharacterized|hypothetical|putative|unannotated|predicted|^$", t):
        return "Sin caracterizar / putativa"
    if "hydrolase" in t:
        return "Otras hidrolasas"
    return "Otras"


# --- 4. Género o grupo a partir del nombre de especie -------------------
def grupo(sp):
    """Primera palabra del nombre de especie (el género). Casos especiales:
    'Candidatus Género' se mantiene junto, y 'X bacterium' (bacterias sin
    género asignado) se agrupa como 'X bacterium'."""
    w = sp.split()
    if not w:
        return "NA"
    if w[0] == "Candidatus" and len(w) > 1:
        return "Candidatus " + w[1]
    if len(w) > 1 and w[1] == "bacterium":
        return w[0] + " bacterium"
    return w[0]


s["Categoria"] = s["Description"].map(categoria)
s["Tipo"] = s["Description"].map(tipo)
s["Genero_o_grupo"] = s["Species"].map(grupo)

# --- 5. Tablas ---------------------------------------------------------
cols = ["Target Name", "Target Accession", "Score", "E-value", "Target Length",
        "Description", "Categoria", "Tipo", "Species", "Genero_o_grupo", "Tax-Id"]
tabla = s[cols].sort_values("Score", ascending=False)       # de mayor a menor puntuación
tabla.to_csv(f"{OUT}/candidatos_{ETIQUETA}.csv", index=False)  # todas las candidatas
tabla.head(50).to_csv(f"{OUT}/top50_candidatos.csv", index=False)  # las 50 mejores

# Resumen numérico en pantalla (son las cifras que se citan en Resultados)
print("Candidatos:", len(tabla))
print("Especies distintas:", tabla["Species"].nunique(),
      "| géneros/grupos:", tabla["Genero_o_grupo"].nunique())
print(tabla["Categoria"].value_counts())
print(tabla["Genero_o_grupo"].value_counts().head(20))

# --- 6. Estilo común de las figuras ------------------------------------
AZUL, GRIS, TXT = "#2a6f97", "#8a96a3", "#222222"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                     "axes.edgecolor": GRIS})


def limpiar(ax):
    """Quita los bordes superior y derecho del gráfico (estética)."""
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    ax.tick_params(colors=TXT)


# --- 7. Figura 1: histograma de puntuaciones ---------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.2))
bins = list(range(75, 775, 25))                     # intervalos de 25 bits
ax.hist(tabla["Score"], bins=bins, color=AZUL, edgecolor="white", linewidth=1)
ax.axvline(UMBRAL, color=TXT, lw=1, ls="--")        # línea vertical en el umbral
ax.annotate(f"Umbral = {UMBRAL} bits".replace(".", ","), (UMBRAL, ax.get_ylim()[1] * 0.92),
            xytext=(8, 0), textcoords="offset points", color=TXT, fontsize=10)
ax.set_yscale("log")                                # escala logarítmica en el eje Y
ax.set_xlabel("Puntuación (bits)")
ax.set_ylabel("Número de proteínas (escala log)")
limpiar(ax)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_histograma_puntuaciones.png", dpi=300)
plt.close(fig)

# --- 8. Figura 2: tipo de proteína según la anotación (3 grupos) --------
NARANJA, GRIS_C = "#e07a2f", "#b4bcc5"


def grupo3(t):
    """Agrupa los 'tipos' en solo 3 grupos para la figura."""
    if t.startswith("6-amino"):
        return "6-aminohexanoato hidrolasa"
    if "eptidasa" in t or "esterasa" in t:
        return "Peptidasas"
    return "Otras hidrolasas,\nsin caracterizar y otras"


g3 = tabla["Tipo"].map(grupo3).value_counts()
g3 = g3.reindex(["Otras hidrolasas,\nsin caracterizar y otras",
                 "6-aminohexanoato hidrolasa", "Peptidasas"])   # orden de las barras
colores = {"Peptidasas": AZUL, "6-aminohexanoato hidrolasa": NARANJA}
fig, ax = plt.subplots(figsize=(7.5, 3.2))
b = ax.barh(g3.index, g3.values, color=[colores.get(k, GRIS_C) for k in g3.index],
            height=0.6)
for r, v in zip(b, g3.values):
    # Porcentaje con coma decimal y SIN espacio antes del % (p. ej. 41,0%)
    pc = f"{v / len(tabla) * 100:.1f}".replace(".", ",")
    ax.text(v + g3.max() * 0.012, r.get_y() + r.get_height() / 2,
            f"{v:,}".replace(",", ".") + f" ({pc}%)",
            va="center", color=TXT, fontsize=10.5)
ax.set_xlabel("Número de proteínas")
ax.set_xlim(0, g3.max() * 1.38)
ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{int(v):,}".replace(",", ".")))
limpiar(ax)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_categorias_anotacion.png", dpi=300)
plt.close(fig)

# --- 9. Figura 3: géneros/grupos más frecuentes -------------------------
g = tabla["Genero_o_grupo"].value_counts().head(20).sort_values()   # 20 más frecuentes
fig, ax = plt.subplots(figsize=(7.5, 6.0))
b = ax.barh(g.index, g.values, color=AZUL, height=0.65)
for r, v in zip(b, g.values):
    ax.text(v + g.max() * 0.01, r.get_y() + r.get_height() / 2, str(v),
            va="center", color=TXT, fontsize=9)
ax.set_xlabel("Número de proteínas")
ax.set_xlim(0, g.max() * 1.1)
limpiar(ax)
fig.tight_layout()
fig.savefig(f"{OUT}/fig_generos_frecuentes.png", dpi=300)
plt.close(fig)
