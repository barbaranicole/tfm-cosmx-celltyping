#!/usr/bin/env Rscript
# ==========================================
# INICIALIZACIÓN Y LIBRERÍAS
# ==========================================
options(warn = 1) 

suppressPackageStartupMessages({
  library(Seurat)
})

# ==========================================
# 1. LECTURA DE ARGUMENTOS
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  stop("Error. Uso: Rscript 01_Seurat_Clustering.R <input_preprocesado.rds> <dims> <resolution> <directorio_salida>")
}

in_rds     <- args[1]
dims_val   <- as.numeric(args[2])
res_val    <- as.numeric(args[3])
out_dir    <- args[4]

# Nos aseguramos de que la carpeta de salida exista
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==========================================
# 2. GENERACIÓN DINÁMICA DE NOMBRES
# ==========================================
# Extraemos el nombre base del archivo de entrada (ej. "data_subset" de "data_subset.rds")
nombre_base <- tools::file_path_sans_ext(basename(in_rds))

# Creamos los nombres de salida con los parámetros incrustados
nombre_salida <- sprintf("%s_dims%d_res%s", nombre_base, dims_val, res_val)

tiempo_actual <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_txt       <- file.path(out_dir, paste0(nombre_salida, "_", tiempo_actual, ".txt"))
log_file      <- file.path(out_dir, paste0("Log_", nombre_salida, "_", tiempo_actual, ".txt"))

# ==========================================
# 3. CONFIGURACIÓN DEL LOG (SINK)
# ==========================================
con <- file(log_file, open = "wt")
sink(con)
sink(con, type = "message")

cat("========================================\n")
cat(" REPORTE DE EJECUCIÓN: Seurat Clustering \n")
cat(" Fecha:", as.character(Sys.time()), "\n")
cat("========================================\n")
cat("[PARÁMETROS UTILIZADOS]\n")
cat(" -> Input RDS  :", in_rds, "\n")
cat(" -> Dims (PCA) : 1 a", dims_val, "\n")
cat(" -> Resolution :", res_val, "\n")
cat(" -> Output TXT :", out_txt, "\n")
cat(" -> Log File   :", log_file, "\n")
cat("========================================\n")

# ==========================================
# 4. CARGA DEL OBJETO
# ==========================================
cat("\n[Seurat] 1. Cargando objeto Seurat...\n")
data <- readRDS(in_rds)

cat("\n========================================\n")
cat("[PASO] Objeto Seurat inicial\n")
cat("Hora:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  -> Tipo de objeto :", class(data)[1], "\n")
cat("  -> Dimensiones    :", dim(data)[1], "filas x", dim(data)[2], "columnas\n")
cat("  -> Peso en RAM    :", format(object.size(data), units = "MB"), "\n")
cat("========================================\n")

# ==========================================
# 5. CLUSTERING (FindNeighbors & FindClusters)
# ==========================================
cat("\n[Seurat] 2. Calculando grafo de vecinos (FindNeighbors)...\n")
# Nota: Asumimos que RunPCA ya se ejecutó en el objeto de entrada
set.seed(357)
data <- FindNeighbors(data, dims = 1:dims_val)

cat("\n[Seurat] 3. Identificando clusters (FindClusters)...\n")
data <- FindClusters(data, resolution = res_val, algorithm = 4)

# ==========================================
# 6. RESUMEN DE CLUSTERS DETECTADOS
# ==========================================
cat("\n========================================\n")
cat("[RESUMEN DE CLUSTERING]\n")

clusters_detectados <- length(unique(data$seurat_clusters))
cat(" -> Total de clusters detectados:", clusters_detectados, "\n\n")
cat(" -> Distribución de células por cluster:\n")

tabla_clusters <- table(data$seurat_clusters)
for (i in seq_along(tabla_clusters)) {
  cat(sprintf("    - Cluster %-2s : %5d células\n", names(tabla_clusters)[i], tabla_clusters[i]))
}
cat("========================================\n")


# ==========================================
# 7. GUARDADO DEL OBJETO
# ==========================================
cat("\n[Seurat] 4. Guardando objeto final...\n")

resultados <- data.frame(
   cell_ID = colnames(data),
   label_seurat = data$seurat_clusters,
   stringsAsFactors = FALSE
)

write.table(resultados, file = out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[Seurat] ¡Proceso completado con éxito en:", out_txt, "!\n")

# ==========================================
# 8. CIERRE DEL LOG
# ==========================================
sink(type = "message")
sink()
close(con)