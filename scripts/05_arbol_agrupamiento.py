"""
05_arbol_agrupamiento.py
------------------------------------------------------------------------
QUÉ HACE (apartado 4.3 del TFM)
    Dibuja el "árbol de agrupamiento" que da Clustal Omega (fichero
    .phylotree, formato Newick) con el nombre de cada proteína coloreado
    según su familia: NylA (naranja), NylB (azul) y NylC (verde).

    El árbol NO es un análisis filogenético en sentido estricto (no hay
    modelo de sustitución ni soporte estadístico); solo representa las
    relaciones de similitud entre las secuencias.

USO (desde la terminal)
    python 05_arbol_agrupamiento.py arbol.phylotree "Título" salida.png

    arbol.phylotree   fichero .phylotree de Clustal Omega
    "Título"          título que se escribe sobre la figura
    salida.png        imagen de salida

EJEMPLO
    python scripts/05_arbol_agrupamiento.py data/analisis_identidad_seq/14seq_5iter/arbol_14seq_5iter.phylotree \
        "Árbol de agrupamiento de las 14 proteínas (5 iteraciones)" results/arboles/arbol_14seq_5iter.png

REQUISITOS
    Python 3, Biopython y matplotlib (en el TFM: Python 3.13.16,
    Biopython 1.88, matplotlib 3.11.2).
    Instalación:  pip install biopython matplotlib

FORMATO DE LAS ETIQUETAS DEL ÁRBOL
    Las etiquetas del .phylotree tienen la forma  ACCESIÓN|NOMBRE|ORGANISMO
    (p. ej. "P07062|NYLB2_PAEUR|Paenarthrobacter_ureafaciens"). En la figura
    se muestra  NOMBRE (ACCESIÓN)  para que sea legible.
------------------------------------------------------------------------
"""
import sys
import os
from Bio import Phylo                      # lectura y dibujo de árboles
import matplotlib
matplotlib.use("Agg")                      # guarda la figura directamente en fichero
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# --- 1. Argumentos -----------------------------------------------------
if len(sys.argv) < 4:
    sys.exit('Uso: python 05_arbol_agrupamiento.py arbol.phylotree "Título" salida.png')
archivo_newick, titulo, archivo_salida = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(os.path.dirname(archivo_salida) or ".", exist_ok=True)   # crea la carpeta de salida si no existe

# --- 2. Colores por familia ---------------------------------------
# Para cada familia: (color, texto de la leyenda).
COLORES = {
    "NylA": ("#e08214", "NylA"),
    "NylB": ("#3366cc", "NylB"),
    "NylC": ("#1a7a1a", "NylC"),
}


def familia(nombre):
    """Devuelve la familia a partir del NOMBRE de la proteína
    (sin distinguir mayúsculas de minúsculas).
      - empieza por NYLA -> NylA ; empieza por NYLB -> NylB
      - NYLC..., NYL50, TVGC -> NylC (familia NylC y homólogos)
    Si no encaja en ninguna, devuelve None (se dibuja en negro y se avisa)."""
    n = nombre.upper()
    if n.startswith("NYLA"):
        return "NylA"
    if n.startswith("NYLB"):
        return "NylB"
    if n.startswith("NYLC") or n.startswith("NYL50") or n.startswith("TVGC"):
        return "NylC"
    return None


# --- 3. Leer el árbol y preparar las etiquetas -------------------------
tree = Phylo.read(archivo_newick, "newick")

color_de_etiqueta = {}                     # etiqueta mostrada -> color
for clado in tree.get_terminals():         # solo las hojas (las proteínas)
    partes = clado.name.split("|")         # ACCESIÓN | NOMBRE | ORGANISMO
    accesion = partes[0]
    nombre = partes[1] if len(partes) > 1 else partes[0]
    fam = familia(nombre)
    if fam is None:
        print("AVISO: no se reconoce la familia de", nombre)
    clado.name = f"{nombre} ({accesion})"  # etiqueta que se verá en la figura
    color_de_etiqueta[clado.name] = COLORES[fam][0] if fam else "black"

# --- 4. Dibujar --------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 7))
Phylo.draw(tree, axes=ax, do_show=False,
           label_colors=lambda etiqueta: color_de_etiqueta.get(etiqueta, "black"),
           branch_labels=None)             # sin números sobre las ramas

# Margen a la derecha para que las etiquetas largas no se salgan del marco
x0, x1 = ax.get_xlim()
ax.set_xlim(x0, x1 + 0.22 * (x1 - x0))

ax.set_ylabel("")                          # sin título ni marcas en el eje vertical
ax.set_yticks([])
ax.set_xlabel("Longitud de rama (distancia estimada por Clustal Omega)", fontsize=10)
ax.set_title(titulo, fontsize=13, fontweight="bold", pad=15)

# Leyenda de colores
elementos = [mlines.Line2D([0], [0], color=c, lw=3, label=txt)
             for c, txt in COLORES.values()]
ax.legend(handles=elementos, loc="upper left", fontsize=9,
          frameon=True, title="Familia")

plt.tight_layout()
plt.savefig(archivo_salida, dpi=300, bbox_inches="tight")
plt.close()
print("Generado:", archivo_salida)
