"""
03_triada.py
------------------------------------------------------------------------
QUÉ HACE (apartado 4.4 del TFM)
    Comprueba si las proteínas de un alineamiento múltiple conservan el
    sitio de autocorte (Asn-266 / Thr-267) y los residuos Asp-306 y Asp-308
    de la tríada catalítica de NylC, aunque estén en posiciones distintas
    de cada secuencia.

IDEA
    Las posiciones 266, 267, 306 y 308 están descritas para la proteína de
    referencia Q79F77 (NylC de Paenarthrobacter ureafaciens). En las demás
    secuencias esos mismos residuos pueden tener otro número (porque cada
    proteína puede tener un extremo N-terminal más largo o más corto, o
    inserciones/deleciones). El alineamiento dice qué residuo de cada
    secuencia ocupa la misma COLUMNA que cada posición de la referencia.

USO (desde la terminal)
    python scripts/03_triada.py data/analisis_identidad_seq/8nylc_0iter/alineamiento_8nylc_0iter.aln-fasta results/triada/triada_NylC.tsv [id_referencia]

    alineamiento.aln-fasta  alineamiento múltiple (FASTA alineado o Clustal)
    salida.tsv              tabla de resultados (se crea)
    id_referencia           (opcional) texto que identifica la referencia en
                            las cabeceras; por defecto busca "Q79F77" o
                            "NYLC_PAEUR"

SALIDA (TSV, una línea por secuencia y sitio)
    secuencia  sitio  residuo  posicion  conservado
      - sitio       Asn266, Thr267, Asp306 o Asp308 (numeración de Q79F77)
      - residuo     aminoácido que hay en esa columna (o "-" si es hueco)
      - posicion    número de ese residuo en SU PROPIA secuencia ("-" si hueco)
      - conservado  SI si el residuo es el esperado (N, T, D, D); si no, NO

REQUISITOS
    Python 3 y Biopython (en el TFM: Python 3.13.16, Biopython 1.88).
    Instalación:  pip install biopython
------------------------------------------------------------------------
"""
import sys
import os
from Bio import AlignIO   # lectura de alineamientos múltiples

# --- 1. Argumentos de la línea de comandos -----------------------------
# sys.argv[0] es el nombre del script; [1] y [2] son los dos argumentos
# obligatorios; [3] es opcional.
if len(sys.argv) < 3:
    sys.exit("Uso: python 03_triada.py alineamiento.aln-fasta salida.tsv [id_referencia]")
ruta_aln, ruta_out = sys.argv[1], sys.argv[2]
ref_id = sys.argv[3].lower() if len(sys.argv) > 3 else None
os.makedirs(os.path.dirname(ruta_out) or ".", exist_ok=True)   # crea la carpeta de salida si no existe

# --- 2. Sitios que se comprueban ---------------------------------------
# nombre del sitio: (aminoácido esperado, posición en la proteína de referencia)
# Asn266/Thr267: sitio de autocorte; Thr267, Asp306 y Asp308: tríada catalítica.
SITIOS = {"Asn266": ("N", 266), "Thr267": ("T", 267),
          "Asp306": ("D", 306), "Asp308": ("D", 308)}

# --- 3. Leer el alineamiento -------------------------------------------
# Se prueba primero como FASTA alineado (la salida .aln-fasta de Clustal
# Omega) y, si falla, como formato Clustal clásico.
aln = None
for fmt in ("fasta", "clustal"):
    try:
        aln = AlignIO.read(ruta_aln, fmt)
        break
    except Exception:
        pass
if aln is None:
    sys.exit("No se pudo leer el alineamiento (FASTA o Clustal).")

# --- 4. Localizar la secuencia de referencia ---------------------------
# Se busca la primera secuencia cuyo identificador contenga alguno de los
# textos de 'claves' (sin distinguir mayúsculas de minúsculas).
claves = [ref_id] if ref_id else ["q79f77", "nylc_paeur"]
ref = next((r for r in aln if any(k in r.id.lower() for k in claves)), None)
if ref is None:
    sys.exit("No encuentro la referencia. IDs: " + ", ".join(r.id for r in aln))
print("Referencia:", ref.id)

# --- 5. Tabla: posición en la referencia -> columna del alineamiento ----
# Se recorre la referencia columna a columna. 'n' cuenta los residuos reales
# (los huecos "-" no cuentan). col_de[266] será la columna donde está el
# residuo 266 de la referencia.
col_de, n = {}, 0
for col, aa in enumerate(str(ref.seq)):
    if aa != "-":
        n += 1
        col_de[n] = col

# --- 6. Comprobar cada secuencia en cada sitio -------------------------
total, conservados = 0, 0
with open(ruta_out, "w", encoding="utf-8") as f:
    f.write("secuencia\tsitio\tresiduo\tposicion\tconservado\n")
    for r in aln:                      # cada secuencia del alineamiento
        s = str(r.seq)                 # su secuencia alineada (con huecos)
        for nombre, (aa_esp, pos_ref) in SITIOS.items():
            c = col_de[pos_ref]        # columna del sitio
            aa = s[c]                  # residuo de esta secuencia en esa columna
            # Posición del residuo en la propia secuencia = nº de residuos
            # (sin huecos) hasta esa columna incluida. Si hay hueco, "-".
            pos = len(s[:c + 1].replace("-", "")) if aa != "-" else "-"
            ok = (aa == aa_esp)        # ¿es el aminoácido esperado?
            total += 1
            conservados += ok
            f.write(f"{r.id}\t{nombre}\t{aa}\t{pos}\t{'SI' if ok else 'NO'}\n")

# Resumen en pantalla
print(f"Secuencias: {len(aln)} | sitios evaluados: {total} | conservados: {conservados}")
print("Guardado en", ruta_out)
