#!/usr/bin/env Rscript
# ==========================================
# 03_Marcadores.R
# Calcula FindAllMarkers para cada columna de clustering en meta.data
# y guarda SOLO la tabla completa (sin Top10) por resolución/modelo.
# Pensado para usarse después de 01_Resultados_clustering.R y
# Clustering_joinmerged.R, sobre cualquier RDS *_joinmerged_clustering.rds.
# ==========================================

library(Seurat)
library(dplyr)

# --- CONFIGURACIÓN ---
# Cambiar solo esta ruta al RDS que quieras revisar
merged_rds  <- path.expand("~/Documents/Master/TFM/Resultados_Algoritmos_v1/thin_75_joinmerged_clustering.rds")
base_salida <- path.expand("~/Documents/Master/TFM/Resultados_Algoritmos_v1")

# --- CARGA DEL OBJETO Y DETECCIÓN DE COLUMNAS DE CLUSTERING ---
merged <- readRDS(merged_rds)
nms    <- colnames(merged@meta.data)

leiden_cols <- nms[grepl("^Leiden_.*res[0-9]", nms)]
leiden_cols <- leiden_cols[order(as.numeric(gsub(".*res([0-9.]+)$", "\\1", leiden_cols)))]
ist_cols    <- nms[grepl("^IST_.*nclusts[0-9]", nms) & !grepl("_conf$", nms)]
ist_cols    <- ist_cols[order(as.numeric(gsub(".*nclusts([0-9]+).*", "\\1", ist_cols)))]

columnas_clustering <- c(leiden_cols, ist_cols)
if (length(columnas_clustering) == 0) {
  stop("No se encontraron columnas de clustering en meta.data (Leiden_/IST_).")
}

# Inferir la versión (original / thin_75) a partir de las columnas
version_tag <- NA_character_
if (length(leiden_cols) > 0) {
  base_leiden <- sub("^Leiden_", "", leiden_cols[1])
  version_tag <- strsplit(base_leiden, "_res", fixed = TRUE)[[1]][1]
} else if (length(ist_cols) > 0) {
  base_ist <- sub("^IST_", "", ist_cols[1])
  version_tag <- strsplit(base_ist, "_nclusts", fixed = TRUE)[[1]][1]
} else {
  version_tag <- "desconocida"
}

dir_version <- file.path(base_salida, paste0("Resultados_", version_tag))
if (!dir.exists(dir_version)) dir.create(dir_version, recursive = TRUE)

cat("RDS:", basename(merged_rds),
    "-> versión detectada:", version_tag,
    "| Carpeta versión:", dir_version, "\n")
cat("Columnas de clustering:", paste(columnas_clustering, collapse = ", "), "\n\n")

# --- FindAllMarkers POR COLUMNA (solo tabla completa) ---
for (res_actual in columnas_clustering) {
  nombre_archivo <- paste0("Marcadores_", res_actual, ".txt")
  ruta_salida    <- file.path(dir_version, nombre_archivo)

  if (file.exists(ruta_salida)) {
    cat("[FindAllMarkers]", res_actual, "ya existe en", dir_version, "- omitiendo.\n")
    next
  }

  cat("========================================\n")
  cat("[FindAllMarkers] Versión:", version_tag,
      "| Columna:", res_actual, "\n")
  cat("Salida:", ruta_salida, "\n")
  cat("Hora de inicio:", format(Sys.time(), "%H:%M:%S"), "\n")

  Idents(merged) <- res_actual
  marcadores <- FindAllMarkers(
    merged,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25
  )

  write.table(
    marcadores,
    file      = ruta_salida,
    sep       = "\t",
    quote     = FALSE,
    row.names = FALSE
  )

  cat("-> Completado: tabla completa guardada para", res_actual, "\n\n")
}

cat("Todos los marcadores han sido calculados y exportados.\n")
