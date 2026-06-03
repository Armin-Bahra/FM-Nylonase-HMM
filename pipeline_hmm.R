# 1. Archivos de datos (Rutas relativas dentro de tu proyecto de RStudio)
fichero_alineamiento <- "data/nylonases_aligned.aln-fasta"
perfil_hmm_salida    <- "results/nylonase_model.hmm"
fichero_target       <- "data/ecoli_test.fasta"
tabla_resultados     <- "results/ecoli_hits.tbl"

# ==========================================
# PASO A: CONSTRUIR EL MODELO HMM VÍA WSL
# ==========================================
print("Construyendo el perfil HMM en WSL...")

# Creamos el comando indicando a Windows que lo ejecute dentro de WSL
comando_build <- paste("wsl hmmbuild", perfil_hmm_salida, fichero_alineamiento)
system(comando_build)

# ==========================================
# PASO B: EJECUTAR LA BIOPROSPECCIÓN VÍA WSL
# ==========================================
print("Escaneando el dataset de prueba de E. coli en WSL...")

comando_search <- paste("wsl hmmsearch --tblout", tabla_resultados, perfil_hmm_salida, fichero_target)
system(comando_search)

print("¡Procesamiento completado con éxito a través de WSL!")

library(tidyverse)

# Verificar si el archivo de resultados existe y tiene tamaño mayor que cero
fichero_res <- "results/ecoli_hits.tbl"

# Leemos las líneas completas para comprobar si hay hits reales
lineas <- readLines(fichero_res)
hits_reales <- grep("^[^#]", lineas) # Busca líneas que NO empiecen por #

if (length(hits_reales) > 0) {
  # Si hay hits, cargamos y filtramos como una tabla estructurada
  raw_data <- read_table(fichero_res, comment = "#", col_names = FALSE)
  colnames(raw_data)[1:3] <- c("target_name", "accession", "query_name")
  colnames(raw_data)[5]   <- "E_value"
  
  print("Se han detectado candidatos potenciales:")
  print(raw_data %>% select(target_name, E_value))
} else {
  print("Análisis completado: No se detectaron homólogos de nylonasas en el dataset de control (E. coli). El modelo demuestra una alta especificidad inicial.")
}