# ==========================================
# ÁRBOL GUÍA Y ÁRBOL FILOGENÉTICO COLOREADOS POR SUPERFAMILIA
# TFM - Diseño de un modelo oculto de Markov para la identificación de nuevas nylonasas
# ==========================================

from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ------------------------------------------
# Configuración común: colores por superfamilia
# ------------------------------------------
color_map = {
    "NYLC2": ("#66c266", "NylC (putativa)"),
    "NYLC":  ("#1a7a1a", "NylC"),
    "NYLA":  ("#e08214", "NylA"),
    "NYLB":  ("#3366cc", "NylB"),
}

def get_color(name):
    if name is None:
        return "black"
    for key, (color, label) in color_map.items():
        if key in name:
            return color
    return "black"

def dibujar_arbol(archivo_newick, titulo, archivo_salida, separador="|"):
    tree = Phylo.read(archivo_newick, "newick")

    # Reformatear etiquetas a "Nombre (Accession)"
    for clade in tree.find_clades():
        if clade.name:
            partes = clade.name.split(separador, 1)
            if len(partes) == 2:
                accession, nombre = partes
                clade.name = f"{nombre} ({accession})"

    fig, ax = plt.subplots(figsize=(11, 7))
    Phylo.draw(tree, axes=ax, do_show=False,
               label_colors=lambda name: get_color(name),
               branch_labels=None)

    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_xlabel("Distancia (sustituciones por posición)", fontsize=10)
    ax.set_title(titulo, fontsize=13, fontweight="bold", pad=15)

    legend_elements = [mlines.Line2D([0], [0], color=c, lw=3, label=lbl)
                        for c, lbl in dict.fromkeys(color_map.values())]
    ax.legend(handles=legend_elements, loc="upper left", fontsize=9,
              frameon=True, title="Superfamilia")

    plt.tight_layout()
    plt.savefig(archivo_salida, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Generado: {archivo_salida}")

# ------------------------------------------
# 1. Árbol guía (guide tree)
# ------------------------------------------
dibujar_arbol(
    archivo_newick="clustalo_8secuencias_guide.tree",
    titulo="Árbol guía para las 8 proteínas inicialmente seleccionadas",
    archivo_salida="arbol_guia_8secuencias.png",
    separador="|"
)

# ------------------------------------------
# 2. Árbol filogenético
# ------------------------------------------
dibujar_arbol(
    archivo_newick="clustalo_8secuencias_filogenetico.tree",
    titulo="Árbol filogenético para las 8 proteínas inicialmente seleccionadas",
    archivo_salida="arbol_filogenetico_8secuencias.png",
    separador="_"
)
