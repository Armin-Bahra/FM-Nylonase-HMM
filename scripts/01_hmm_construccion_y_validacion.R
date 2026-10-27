# =====================================================================
# 01_hmm_construccion_y_validacion.R
# ---------------------------------------------------------------------
# QUÉ HACE ESTE SCRIPT (apartado 4.5 del TFM)
#   1. Construye el perfil HMM (hmmbuild) a partir del alineamiento de las
#      5 NylC de construcción (Q79F77, Q1EPR4, Q1EPR5, NylC1 y NylC2).
#   2. Busca el perfil contra el CONTROL POSITIVO (C+: NylC5, Nyl50 y TvgC)
#      con hmmsearch y guarda la tabla de resultados (--tblout).
#   3. Prepara el CONTROL NEGATIVO (C-): une los FASTA descargados de
#      UniProt (Swiss-Prot, revisadas), quita las secuencias repetidas y
#      las de NylC, y CUENTA cuántas había en cada fichero, cuántas estaban
#      repetidas y cuántas se descartaron.
#   4. Busca el perfil contra el C- con hmmsearch (filtros desactivados y
#      umbral de notificación muy laxo, para ver la puntuación de TODAS las
#      secuencias, incluso las que casi no se parecen al perfil).
#   5. Calcula el margen entre la puntuación MÁXIMA del C- y la MÍNIMA del
#      C+. Entre ambas se fija el umbral de confianza (bit-score).
#
# REQUISITOS
#   - R con el paquete tidyverse.
#   - HMMER 3.4 instalado dentro de WSL (Ubuntu). Cada comando de HMMER se
#     lanza desde R con el prefijo "wsl ".
#   - Carpetas "data/" (entradas) y "results/" (salidas) en el directorio
#     de trabajo de R (la raíz del repositorio).
#
# ENTRADAS (dentro de data/)
#   analisis_identidad_seq/5nylc_construccion/alineamiento_5nylc_construccion.aln-fasta
#                                    alineamiento Clustal Omega de las 5 NylC
#   controles/positivo/control_positivo.fasta   NylC5, Nyl50 y TvgC
#   controles/negativo/uniprotkb_lipase_..._2026_10_08.fasta
#   controles/negativo/uniprotkb_aminotransferase_..._2026_10_08.fasta
#   controles/negativo/uniprotkb_dmpA_OR_bapA_..._2026_10_08.fasta (S58 y afines)
#
# SALIDAS (en results/hmm/)
#   nylonase_model.hmm                      perfil HMM construido
#   control_positivo_hits.tbl / control_negativo_hits.tbl   tablas de hmmsearch
#   control_negativo_conteos.csv            recuentos del control negativo
#   control_negativo_hits_ordenados.csv     puntuaciones del C-, de mayor a menor
#   versiones.txt                           versiones del software usado
#   (y data/controles/negativo/control_negativo_completo.fasta, el FASTA unificado)
#
# Ejecutar SIEMPRE desde la carpeta raíz del repositorio.
# =====================================================================

library(tidyverse)   # readr, dplyr, tibble... (en el TFM: tidyverse 2.0.0)

# ---------------------------------------------------------------------
# 0. AJUSTES GENERALES
# ---------------------------------------------------------------------
# Prefijo para ejecutar HMMER. En tu ordenador (Windows) HMMER está dentro
# de WSL, por eso es "wsl ". Si algún día lo ejecutas directamente en
# Linux, cámbialo por "" (cadena vacía).
prefijo <- "wsl "

# Rutas de entrada y salida
fichero_alineamiento <- "data/analisis_identidad_seq/5nylc_construccion/alineamiento_5nylc_construccion.aln-fasta"
perfil_hmm_salida    <- "results/hmm/nylonase_model.hmm"
fichero_target       <- "data/controles/positivo/control_positivo.fasta"
tabla_resultados     <- "results/hmm/control_positivo_hits.tbl"

# Crea la carpeta de resultados si no existe (no da error si ya existe)
dir.create("results/hmm", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------
# 0.1 REGISTRO DE VERSIONES (para el apartado de Metodología y el README)
# ---------------------------------------------------------------------
# Guarda en results/versiones.txt la versión de R, de tidyverse y de HMMER
# tal como se han usado de verdad en esta ejecución.
version_hmmer <- system(paste0(prefijo, "hmmsearch -h"), intern = TRUE)[2]
versiones <- c(
  R.version.string,
  paste("tidyverse", as.character(packageVersion("tidyverse"))),
  version_hmmer
)
writeLines(versiones, "results/hmm/versiones.txt")
print(versiones)

# ---------------------------------------------------------------------
# PASO 1. CONSTRUIR EL PERFIL HMM (hmmbuild)
# ---------------------------------------------------------------------
# hmmbuild <salida.hmm> <alineamiento>
# Usa los parámetros por defecto. El modelo se nombra automáticamente a
# partir del nombre del archivo del alineamiento (sin extensión); por eso ese
# nombre aparece como "Query Name" en las tablas de HMMER-EBI. (El perfil
# usado en el TFM se construyó con el archivo "perfilhmm_aligned"; si se
# reconstruye con el nombre nuevo, solo cambia esa etiqueta, no el modelo.)
print("Construyendo el perfil HMM en WSL...")

comando_build <- paste0(prefijo, "hmmbuild ", perfil_hmm_salida, " ", fichero_alineamiento)
system(comando_build)

# ---------------------------------------------------------------------
# PASO 2. COMPROBAR QUE EL HMM SE HA CREADO
# ---------------------------------------------------------------------
# Si hmmbuild falla (p. ej. ruta mal escrita) el fichero no existe y el
# script se detiene aquí con un mensaje claro, en vez de seguir y dar
# errores confusos más abajo.
if (!file.exists(perfil_hmm_salida)) {
  stop("ERROR: No se ha podido crear el archivo del perfil HMM.")
} else {
  print("Perfil HMM creado correctamente.")
}

# ---------------------------------------------------------------------
# PASO 3. CONTROL POSITIVO (C+): hmmsearch del perfil contra NylC5, Nyl50 y TvgC
# ---------------------------------------------------------------------
# hmmsearch --tblout <tabla> <perfil.hmm> <secuencias.fasta>
#   --tblout : guarda una tabla con una línea por secuencia encontrada
# El resto de parámetros son los POR DEFECTO de HMMER (umbral de notificación
# E-value <= 10).
print("Analizando el control positivo (C+)...")

comando_search <- paste0(prefijo, "hmmsearch --tblout ", tabla_resultados, " ",
                         perfil_hmm_salida, " ", fichero_target)
system(comando_search)

# ---------------------------------------------------------------------
# PASO 4. COMPROBAR QUE HAY TABLA DE RESULTADOS
# ---------------------------------------------------------------------
if (!file.exists(tabla_resultados)) {
  stop("ERROR: No se ha generado el archivo de resultados.")
}

# ---------------------------------------------------------------------
# PASO 5. LEER Y MOSTRAR LOS HITS DEL C+
# ---------------------------------------------------------------------
# En el formato --tblout las líneas que empiezan por "#" son comentarios.
# Las demás son resultados ("hits"). Contamos cuántas hay.
lineas <- readLines(tabla_resultados)
hits_reales <- grep("^[^#]", lineas)   # posiciones de las líneas que NO empiezan por #

if (length(hits_reales) > 0) {

  # Lee la tabla ignorando los comentarios. Las columnas de --tblout son:
  #   1 target name | 2 accession | 3 query name | 4 accession |
  #   5 E-value (secuencia completa) | 6 score (secuencia completa) | ...
  raw_data <- read_table(tabla_resultados, comment = "#", col_names = FALSE)
  colnames(raw_data)[1] <- "target_name"
  colnames(raw_data)[2] <- "accession"
  colnames(raw_data)[3] <- "query_name"
  colnames(raw_data)[5] <- "E_value"
  colnames(raw_data)[6] <- "score"

  print("======================================")
  print("HITS DETECTADOS EN EL CONTROL POSITIVO")
  print("======================================")
  print(raw_data %>% select(target_name, accession, query_name, E_value, score))

} else {
  print("No se detectaron hits en el control positivo.")
}

# ---------------------------------------------------------------------
# PASO 6. CONTROL NEGATIVO (C-)
# ---------------------------------------------------------------------
# Tres FASTA descargados de UniProtKB (solo Swiss-Prot, "reviewed") el
# 08/10/2026: lipasas, aminotransferasas y peptidasas S58 y afines
# (consulta "dmpA OR bapA"). Ajusta los nombres si los tuyos son distintos.
ficheros_neg <- c(
  "data/controles/negativo/uniprotkb_lipase_AND_reviewed_true_2026_10_08.fasta",
  "data/controles/negativo/uniprotkb_aminotransferase_AND_reviewed_2026_10_08.fasta",
  "data/controles/negativo/uniprotkb_dmpA_OR_bapA_AND_reviewed_tru_2026_10_08.fasta"
)
fichero_target_neg   <- "data/controles/negativo/control_negativo_completo.fasta"   # FASTA unificado (se crea aquí)
tabla_resultados_neg <- "results/hmm/control_negativo_hits.tbl"

# --- 6.1 Juntar los FASTA, quitar repetidas y contar ------------------
# Una misma proteína puede salir en dos búsquedas (p. ej. en "lipase" y en
# "dmpA OR bapA"). Se identifica por el primer "palabra" de la cabecera
# (p. ej. "sp|P13398|NYLA_PAEUR") y solo se conserva la primera aparición.
# Además se eliminan las entradas de NylC (identificador con "NYLC"), para
# que el control negativo no contenga la familia que queremos detectar.

# Lee cada fichero por separado (para poder contar cuántas secuencias tiene)
lineas_por_fichero <- lapply(ficheros_neg, readLines)
n_por_fichero <- sapply(lineas_por_fichero, function(l) sum(grepl("^>", l)))  # nº de cabeceras ">" = nº de secuencias
names(n_por_fichero) <- basename(ficheros_neg)

# Junta todas las líneas en un único vector
lineas_fasta <- unlist(lineas_por_fichero)

# Posición de cada cabecera y de la última línea de cada secuencia
cab    <- grep("^>", lineas_fasta)               # líneas que son cabeceras
inicio <- cab                                     # cada secuencia empieza en su cabecera...
fin    <- c(cab[-1] - 1, length(lineas_fasta))    # ...y termina justo antes de la siguiente

# Identificador de cada secuencia (primer campo de la cabecera, sin ">")
ids <- sub("^>(\\S+).*", "\\1", lineas_fasta[cab])

# A qué fichero pertenece cada secuencia (en el mismo orden que 'cab')
fichero_de_cada <- rep(basename(ficheros_neg), times = n_por_fichero)

# Marcas lógicas (TRUE/FALSE) para cada secuencia
es_duplicada <- duplicated(ids)        # TRUE si el id ya había aparecido antes
es_nylc      <- grepl("NYLC", ids)     # TRUE si es una entrada de NylC
mantener     <- !es_duplicada & !es_nylc

# Escribe el FASTA unificado solo con las secuencias que se conservan
bloques <- unlist(Map(function(i, f) lineas_fasta[i:f],
                      inicio[mantener], fin[mantener]))
writeLines(bloques, fichero_target_neg)

# --- Tabla de recuentos (para Resultados y Anexos) ---
conteos <- tibble(
  concepto = c(
    paste("Secuencias descargadas en", names(n_por_fichero)),
    "Total antes de quitar repetidas",
    paste("Repetidas eliminadas (ya vistas en un fichero anterior), localizadas en", names(n_por_fichero)),
    "Entradas de NylC eliminadas (no repetidas)",
    "Total del control negativo (secuencias únicas)",
    paste("Conservadas procedentes de", names(n_por_fichero))
  ),
  n = c(
    as.integer(n_por_fichero),
    length(ids),
    as.integer(sapply(names(n_por_fichero), function(f) sum(es_duplicada & fichero_de_cada == f))),
    sum(!es_duplicada & es_nylc),
    sum(mantener),
    as.integer(sapply(names(n_por_fichero), function(f) sum(mantener & fichero_de_cada == f)))
  )
)
print(conteos, n = Inf)
write_csv(conteos, "results/hmm/control_negativo_conteos.csv")
cat("Secuencias en el control negativo:", sum(mantener), "\n")

# --- 6.2 hmmsearch contra el C- ---------------------------------------
# Parámetros NO por defecto (se deben mencionar en la Metodología):
#   --max    desactiva TODOS los filtros heurísticos de HMMER (más sensible,
#            más lento): cada secuencia se puntúa de forma completa.
#   -E 1000  umbral de notificación E-value <= 1000 (el de por defecto es 10),
#            para que aparezcan también secuencias con similitud muy débil.
# Objetivo: conocer la puntuación MÁS ALTA que alcanza cualquier secuencia
# del control negativo.
print("Analizando el control negativo (C-)...")

comando_search_neg <- paste0(prefijo, "hmmsearch --max -E 1000 --tblout ",
                             tabla_resultados_neg, " ", perfil_hmm_salida,
                             " ", fichero_target_neg)
system(comando_search_neg)

if (!file.exists(tabla_resultados_neg)) {
  stop("ERROR: No se ha generado el archivo de resultados del control negativo.")
}

# --- 6.3 Leer las tablas ----------------------------------------------
# La última columna de --tblout (descripción) puede tener espacios, y eso
# confunde a read_table. Por eso se leen solo los campos 1, 5 y 6
# (nombre, E-value y puntuación), que nunca tienen espacios.
leer_tblout <- function(f) {
  l <- readLines(f)
  l <- l[!grepl("^#", l)]                 # descarta comentarios
  if (length(l) == 0) {                   # sin hits: devuelve tabla vacía
    return(tibble(target_name = character(),
                  E_value = numeric(), score = numeric()))
  }
  campos <- strsplit(trimws(l), "\\s+")   # separa cada línea por espacios
  tibble(
    target_name = sapply(campos, `[`, 1),
    E_value     = as.numeric(sapply(campos, `[`, 5)),
    score       = as.numeric(sapply(campos, `[`, 6))
  ) %>% arrange(desc(score))              # de mayor a menor puntuación
}

neg <- leer_tblout(tabla_resultados_neg)
pos <- leer_tblout(tabla_resultados)

# Comprobación de seguridad: TODAS las secuencias del control positivo deben
# aparecer en la tabla. Si alguna no se detecta, no tiene puntuación y el
# "mínimo del C+" calculado más abajo sería engañosamente alto.
n_pos_esperadas <- sum(grepl("^>", readLines(fichero_target)))
cat("Secuencias del C+ detectadas:", nrow(pos), "de", n_pos_esperadas, "\n")
if (nrow(pos) < n_pos_esperadas) {
  warning("OJO: alguna secuencia del C+ no se ha detectado; el umbral calculado no es fiable.")
}

print("10 puntuaciones más altas del control negativo:")
print(head(neg, 10))
write_csv(neg, "results/hmm/control_negativo_hits_ordenados.csv")

# Cuántas secuencias del C- ha puntuado HMMER (con -E 1000 pueden ser
# menos que las del FASTA: las que no se parecen nada no se notifican)
cat("Secuencias del C- con puntuación en la tabla:", nrow(neg),
    "de", sum(mantener), "\n")

# --- 6.4 Umbral de confianza ------------------------------------------
# Umbral = valor entre la puntuación máxima del C- y la mínima del C+.
# Si hay separación (min C+ > max C-), cualquier valor en ese margen
# distingue ambos controles. En el TFM se eligió la puntuación mínima del
# C+ (la de NylC5).
max_neg <- if (nrow(neg) > 0) max(neg$score) else NA
min_pos <- min(pos$score)

cat("\nPuntuación máxima C-: ", max_neg, "bits\n")
cat("Puntuación mínima C+: ", min_pos, "bits\n")
cat("Margen de separación: ", max_neg, "-", min_pos, "bits\n")
