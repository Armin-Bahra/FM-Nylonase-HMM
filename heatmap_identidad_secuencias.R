# ==========================================
# HEATMAP DE IDENTIDAD DE SECUENCIA CON AGRUPAMIENTO JERÁRQUICO
# TFM - Diseño de un modelo oculto de Markov para la identificación de nuevas nylonasas
# ==========================================

library(tidyverse)
library(pheatmap)

# ------------------------------------------
# 1. Leer el archivo .pim descargado de Clustal Omega
# ------------------------------------------
lineas <- readLines("clustalo-I20260903-012720-0791-96518903-p1m.pim")
lineas <- lineas[lineas != "" & !grepl("^#", lineas)]
n <- length(lineas)
nombres <- character(n)
matriz <- matrix(NA, n, n)

for (i in seq_along(lineas)) {
  partes <- strsplit(trimws(lineas[i]), "\\s+")[[1]]
  numeros <- as.numeric(tail(partes, n))
  nombre_completo <- paste(partes[2:(length(partes) - n)], collapse = " ")
  nombres[i] <- nombre_completo
  matriz[i, ] <- numeros
}

# ------------------------------------------
# 2. Crear etiquetas con nombre y código de acceso
# ------------------------------------------
etiquetas <- unname(sapply(nombres, function(x) {
  campos <- strsplit(x, "\\|")[[1]]
  paste0(campos[2], "\n(", campos[1], ")")
}))

rownames(matriz) <- etiquetas
colnames(matriz) <- etiquetas

# ------------------------------------------
# 3. Generar el heatmap con dendrograma
# ------------------------------------------
pheatmap(matriz,
         display_numbers = TRUE,
         number_format = "%.1f",
         number_color = "white",
         color = colorRampPalette(c("red", "yellow", "darkgreen"))(100),
         breaks = seq(0, 100, length.out = 101),
         main = "Heatmap con dendrograma de las 8 proteínas seleccionadas\n\n",
         fontsize_number = 8,
         fontsize = 9,
         treeheight_row = 100,
         treeheight_col = 100,
         cutree_rows = 4,
         cutree_cols = 4,
         filename = "heatmap_identidad_secuencias.png",
         width = 11, height = 10)
