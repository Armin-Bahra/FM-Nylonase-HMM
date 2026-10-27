# =====================================================================
# 02_heatmaps_identidad.R
# ---------------------------------------------------------------------
# QUÉ HACE ESTE SCRIPT (apartados 4.3 y 4.4 del TFM)
#   Dibuja mapas de calor (heatmaps) de la identidad de secuencia entre
#   proteínas, a partir de las matrices de identidad (.pim) que genera
#   Clustal Omega. Cada celda muestra el % de identidad entre dos proteínas.
#
#   IMPORTANTE: NO se hace agrupamiento jerárquico. Las proteínas aparecen
#   en el mismo orden que en el fichero .pim de Clustal Omega (cluster_rows
#   y cluster_cols = FALSE).
#
# REQUISITOS
#   - R con el paquete pheatmap (en el TFM: R 4.5.2, pheatmap 1.0.13).
#   - Carpetas "data/" (entradas) y "results/" (salidas) en el directorio
#     de trabajo de R.
#
# ENTRADAS (data/analisis_identidad_seq/<grupo>/): ficheros .pim de Clustal Omega
#   matriz_8nylc.pim          8 NylC (0 iteraciones)
#   matriz_14seq_5iter.pim    14 secuencias, 5 iteraciones
#   matriz_14seq_0iter.pim    14 secuencias, 0 iteraciones
#
# SALIDAS (carpeta results/heatmaps/): tres imágenes PNG
#   heatmap_NylC.png, heatmap_14seq_0iter.png, heatmap_14seq_5iter.png
#
# FORMATO DE UN .pim (para entender el código de lectura)
#   Cabecera con líneas que empiezan por "#", una línea en blanco, y luego
#   una línea por proteína:
#        1: P07062|NYLB2_PAEUR|Paenarthrobacter_ureafaciens  100.00  93.11 ...
#   es decir: número con ":", nombre completo, y n valores de identidad
#   (n = número de proteínas).
# =====================================================================

library(pheatmap)

# ---------------------------------------------------------------------
# 1. PALETA DE COLORES
# ---------------------------------------------------------------------
# Se construye una paleta de 100 colores (uno por cada punto porcentual) con
# más resolución donde están las identidades altas, que es donde hay más
# diferencias entre proteínas parecidas.
colores <- c(
  colorRampPalette(c("#fff5eb", "#fdae6b"))(60),   # 0-60 %: de crema a naranja
  colorRampPalette(c("#fdae6b", "#a1d99b"))(25),   # 60-85 %: de naranja a verde claro
  colorRampPalette(c("#a1d99b", "#00682b"))(15)    # 85-100 %: de verde claro a verde oscuro
)

# ---------------------------------------------------------------------
# 2. FUNCIÓN QUE LEE UN .pim Y DIBUJA EL HEATMAP
# ---------------------------------------------------------------------
# Argumentos:
#   ruta_pim : ruta del fichero .pim de Clustal Omega
#   titulo   : título que se escribe sobre el heatmap
#   ruta_png : ruta de la imagen PNG de salida
hacer_heatmap <- function(ruta_pim, titulo, ruta_png) {

  # Lee todas las líneas y se queda solo con las de datos: quita las líneas
  # vacías y las que empiezan por "#" (cabecera de Clustal).
  lineas <- readLines(ruta_pim)
  lineas <- lineas[lineas != "" & !grepl("^#", lineas)]

  n <- length(lineas)               # número de proteínas
  nombres <- character(n)           # aquí se guardará el nombre de cada proteína
  matriz <- matrix(NA, n, n)        # matriz n x n que se irá rellenando

  for (i in seq_along(lineas)) {
    # Separa la línea en "palabras" por los espacios
    partes <- strsplit(trimws(lineas[i]), "\\s+")[[1]]
    # Los n ÚLTIMOS elementos son los porcentajes de identidad
    matriz[i, ] <- as.numeric(tail(partes, n))
    # El nombre está entre el "1:" inicial y los n valores finales
    nombres[i] <- paste(partes[2:(length(partes) - n)], collapse = " ")
  }

  # Etiqueta corta: el 2.º campo del nombre separado por "|"
  # (p. ej. "P07062|NYLB2_PAEUR|Paenarthrobacter..." -> "NYLB2_PAEUR")
  etiquetas <- sapply(nombres, function(x) strsplit(x, "\\|")[[1]][2])

  # Seguridad: si dos proteínas tuvieran la misma etiqueta, el heatmap
  # mezclaría filas. Se avisa y se detiene.
  stopifnot(!anyDuplicated(etiquetas))

  rownames(matriz) <- etiquetas
  colnames(matriz) <- etiquetas

  # Números que se escriben dentro de las celdas, con una cifra decimal y
  # COMA decimal (98,6), según la tipografía española. Al pasar una matriz
  # de texto a display_numbers, pheatmap la usa tal cual.
  numeros <- matrix(sub("\\.", ",", sprintf("%.1f", matriz)), nrow = n)

  pheatmap(matriz,
           display_numbers = numeros,   # texto de cada celda
           number_color = "black",      # color del texto
           color = colores,             # paleta definida arriba
           breaks = seq(0, 100, length.out = 101),  # escala fija 0-100 % (comparable entre heatmaps)
           main = titulo,
           fontsize_number = 7,         # tamaño de los números dentro de las celdas
           fontsize = 9,                # tamaño del resto del texto
           cluster_rows = FALSE,        # SIN agrupamiento jerárquico: se mantiene el orden del .pim
           cluster_cols = FALSE,
           filename = ruta_png,         # guarda directamente en PNG
           width = 10, height = 8)      # tamaño en pulgadas
}

# Crea la carpeta de resultados si no existe
dir.create("results/heatmaps", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 3. HEATMAPS QUE SE GENERAN
# ---------------------------------------------------------------------
# 3.1 Las 8 NylC (apartado 4.4), alineadas con 0 iteraciones
hacer_heatmap("data/analisis_identidad_seq/8nylc_0iter/matriz_8nylc_0iter.pim",
              "Identidad de secuencia de las NylC (0 iteraciones combinadas)",
              "results/heatmaps/heatmap_NylC.png")

# 3.2 Las 14 secuencias NylA, NylB y NylC (apartado 4.3)
# Matrices de identidad de las 14 secuencias (5 y 0 iteraciones combinadas)
# Identificadores originales de Clustal Omega (22/09/2026), por trazabilidad:
#   5 iteraciones: clustalo-I20260922-121545-0458-96962559-p2m  (resultado final)
#   0 iteraciones: clustalo-I20260922-120454-0681-17317983-p1m
RUTA_5ITER <- "data/analisis_identidad_seq/14seq_5iter/matriz_14seq_5iter.pim"
RUTA_0ITER <- "data/analisis_identidad_seq/14seq_0iter/matriz_14seq_0iter.pim"

hacer_heatmap(RUTA_5ITER,
              "Identidad de secuencia de NylA, NylB y NylC (5 iteraciones)",
              "results/heatmaps/heatmap_14seq_5iter.png")

hacer_heatmap(RUTA_0ITER,
              "Identidad de secuencia de NylA, NylB y NylC (0 iteraciones)",
              "results/heatmaps/heatmap_14seq_0iter.png")
