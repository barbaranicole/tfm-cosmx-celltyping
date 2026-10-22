#!/usr/bin/env Rscript
# ==========================================
# INICIALIZACIÓN Y LIBRERÍAS
# ==========================================
options(warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(SC3)
  library(SingleCellExperiment)
  library(scater)
  library(matrixStats)
})

# ==========================================
# 1. LECTURA DE ARGUMENTOS
# ==========================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Error. Uso: Rscript 03_SC3.R <input_preprocesado.rds> <range_num_clusters> <n_cores> <directorio_salida>")
}

in_rds    <- args[1]
k_input   <- args[2]
n_cores   <- as.numeric(args[3])
out_dir   <- args[4]

k_valores <- eval(parse(text = k_input))

# Nos aseguramos de que la carpeta de salida exista
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==========================================
# 2. GENERACIÓN DINÁMICA DE NOMBRES
# ==========================================
nombre_base    <- tools::file_path_sans_ext(basename(in_rds))
k_label        <- paste(k_valores, collapse = "-")
nombre_salida  <- sprintf("%s_k%s_ncores%d", nombre_base, k_label, n_cores)
tiempo_actual  <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_txt        <- file.path(out_dir, paste0(nombre_salida, "_", tiempo_actual, ".txt"))
log_file       <- file.path(out_dir, paste0("Log_", nombre_salida, "_", tiempo_actual, ".txt"))

# ==========================================
# 3. CONFIGURACIÓN DEL LOG (SINK)
# ==========================================
con <- file(log_file, open = "wt")
sink(con)
sink(con, type = "message")

cat("========================================\n")
cat(" REPORTE DE EJECUCIÓN: SC3 Clustering \n")
cat(" Fecha:", as.character(Sys.time()), "\n")
cat("========================================\n")
cat("[PARÁMETROS UTILIZADOS]\n")
cat(" -> Input RDS   :", in_rds, "\n")
cat(" -> k (clusters):", paste(k_valores, collapse = ", "), "\n")
cat(" -> n_cores     :", n_cores, "\n")
cat(" -> gene_filter : FALSE\n")
cat(" -> svm_num_cells : 5000\n")
cat(" -> Output TXT  :", out_txt, "\n")
cat(" -> Log File    :", log_file, "\n")
cat("========================================\n")

# ==========================================
# 4. CARGA DEL OBJETO Y CONVERSIÓN A SCE
# ==========================================
cat("\n[SC3] 1. Cargando objeto Seurat y convirtiendo a SingleCellExperiment...\n")
input_seurat <- readRDS(in_rds)
assay_name   <- DefaultAssay(input_seurat)
sce_subset   <- as.SingleCellExperiment(input_seurat, assay = assay_name)
counts(sce_subset)    <- as.matrix(counts(sce_subset))
logcounts(sce_subset) <- as.matrix(logcounts(sce_subset))
rowData(sce_subset)$feature_symbol <- rownames(sce_subset)

cat("\n========================================\n")
cat("[PASO] Objeto SingleCellExperiment inicial\n")
cat("Hora:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  -> Tipo de objeto :", class(sce_subset)[1], "\n")
cat("  -> Dimensiones    :", nrow(sce_subset), "genes x", ncol(sce_subset), "células\n")
cat("  -> Peso en RAM   :", format(object.size(sce_subset), units = "MB"), "\n")
cat("========================================\n")

# ==========================================
# 5. CLUSTERING (SC3)
# ==========================================
cat("\n[SC3] 2. Preparando y ejecutando SC3 (sc3_prepare, dists, transfs, kmeans, consens, svm)...\n")
set.seed(357)
sce_subset <- sc3_prepare(sce_subset,
                          n_cores = n_cores,
                          gene_filter = FALSE,
                          svm_num_cells = 5000)
sce_subset <- sc3_calc_dists(sce_subset)
sce_subset <- sc3_calc_transfs(sce_subset)
sce_subset <- sc3_kmeans(sce_subset, ks = k_valores)
sce_subset <- sc3_calc_consens(sce_subset)
sce_subset <- sc3_run_svm(sce_subset, ks = k_valores)

# ==========================================
# 6. RESUMEN Y ETIQUETAS PARA TODOS LOS k
# ==========================================
# Extraer etiquetas para cada k del rango (ej. sc3_8_clusters, sc3_9_clusters, ...)
resultados <- data.frame(cell_ID = colnames(sce_subset),
                         stringsAsFactors = FALSE,
                         row.names = NULL)
for (k in k_valores) {
  col_name   <- paste0("sc3_", k, "_clusters")
  resultados[[paste0("label_sc3_k", k)]] <- colData(sce_subset)[[col_name]]
}

cat("\n========================================\n")
cat("[RESUMEN DE CLUSTERING]\n")
for (k in k_valores) {
  etiquetas_k <- resultados[[paste0("label_sc3_k", k)]]
  n_clust     <- length(unique(etiquetas_k))
  cat(" -> k =", k, ": ", n_clust, " clusters detectados\n")
  tabla_k <- table(etiquetas_k)
  for (i in seq_along(tabla_k)) {
    cat(sprintf("      Cluster %-2s : %5d células\n", names(tabla_k)[i], tabla_k[i]))
  }
  cat("\n")
}
cat("========================================\n")

# ==========================================
# 7. GUARDADO DEL RESULTADO
# ==========================================
cat("\n[SC3] 3. Guardando etiquetas (todas las k)...\n")
write.table(resultados, file = out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[SC3] ¡Proceso completado con éxito en:", out_txt, "!\n")

# ==========================================
# 8. CIERRE DEL LOG
# ==========================================
sink(type = "message")
sink()
close(con)
