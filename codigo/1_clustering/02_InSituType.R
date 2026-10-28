#!/usr/bin/env Rscript
# ==============================================================================
# 02_InSituType.R
# ------------------------------------------------------------------------------
# Descripción: Clustering semisupervisado con InSituType. Combina perfiles de
#              referencia fijos con un número de perfiles libres (n_clusts) para
#              asignar cada célula a una clase, con su probabilidad.
#
# Ejecución:   Servidor del laboratorio (por requerimientos de memoria).
# Uso:         Rscript 02_InSituType.R <input.rds> <neg_vector.rds> <reference.csv>
#                               <n_clusts> <dir_salida>
#              En el TFM: n_clusts = 6, 8 y 10; referencia IO_1k.
# Entradas:    input.rds       Mismo objeto Seurat que Leiden, con las capas de
#                              conteos unidas (JoinLayers), formato requerido
#                              por InSituType (*_joinmerged.rds).
#              neg_vector.rds  Media de sondas negativas por célula (generado
#                              en el control de calidad del laboratorio).
#              reference.csv   Perfiles de referencia (IO_1k.profiles.csv).
# Salidas:     <dir_salida>/<input>_nclusts<n>_<fecha>.txt  cell_ID, clase y
#                                                           probabilidad
#              <dir_salida>/Log_<input>_nclusts<n>_<fecha>.txt  registro
# Siguiente:   2_anotacion/03_Unir_resultados.R
# ==============================================================================

# ==========================================
# INICIALIZACIÓN Y LIBRERÍAS
# ==========================================
options(warn = 1)

suppressPackageStartupMessages({
  library(Seurat)
  library(InSituType)
})

# ==========================================
# 1. LECTURA DE ARGUMENTOS
# ==========================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Error. Uso: Rscript 02_InSituType.R <input.rds> <neg_vector.rds> <reference.csv> <n_clusts> <dir_salida>")
}

in_rds    <- args[1]
in_neg    <- args[2]
ref_path  <- args[3]
n_clusts  <- as.numeric(args[4])
out_dir   <- args[5]

# Nos aseguramos de que la carpeta de salida exista
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==========================================
# 2. GENERACIÓN DINÁMICA DE NOMBRES
# ==========================================
nombre_base   <- tools::file_path_sans_ext(basename(in_rds))
nombre_salida <- sprintf("%s_nclusts%d", nombre_base, n_clusts)
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
cat(" REPORTE DE EJECUCIÓN: InSituType \n")
cat(" Fecha:", as.character(Sys.time()), "\n")
cat("========================================\n")
cat("[PARÁMETROS UTILIZADOS]\n")
cat(" -> Input RDS (Seurat) :", in_rds, "\n")
cat(" -> Input neg_vector   :", in_neg, "\n")
cat(" -> Reference (CSV)   :", ref_path, "\n")
cat(" -> n_clusts          :", n_clusts, "\n")
cat(" -> Output TXT        :", out_txt, "\n")
cat(" -> Log File          :", log_file, "\n")
cat("========================================\n")

# ==========================================
# 4. CARGA DE LOS TRES OBJETOS DE ENTRADA
# ==========================================
cat("\n[IST] 1. Cargando objeto Seurat...\n")
input_seurat <- readRDS(in_rds)

cat("\n========================================\n")
cat("[PASO] Objeto Seurat (input)\n")
cat("Hora:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  -> Tipo de objeto :", class(input_seurat)[1], "\n")
cat("  -> Dimensiones    :", dim(input_seurat)[1], "genes x", dim(input_seurat)[2], "células\n")
cat("  -> Peso en RAM    :", format(object.size(input_seurat), units = "MB"), "\n")
cat("========================================\n")

cat("\n[IST] 2. Cargando vector de negativos...\n")
neg_vector_final <- readRDS(in_neg)

cat("\n========================================\n")
cat("[PASO] neg_vector_final\n")
cat("Hora:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  -> Tipo          :", class(neg_vector_final)[1], "\n")
cat("  -> Longitud      :", length(neg_vector_final), "(una entrada por célula)\n")
cat("  -> Resumen       : min =", min(neg_vector_final), ", max =", max(neg_vector_final), ", media =", round(mean(neg_vector_final), 4), "\n")
cat("========================================\n")

cat("\n[IST] 3. Cargando referencia (perfiles)...\n")
reference <- read.csv(ref_path, row.names = 1)
reference <- as.matrix(reference)

cat("\n========================================\n")
cat("[PASO] Objeto reference (perfiles de referencia)\n")
cat("Hora:", format(Sys.time(), "%H:%M:%S"), "\n")
cat("  -> Dimensiones   :", nrow(reference), "genes x", ncol(reference), "tipos/perfiles\n")
cat("  -> Tipos (columnas):", paste(colnames(reference), collapse = ", "), "\n")
cat("  -> Peso en RAM   :", format(object.size(reference), units = "KB"), "\n")
cat("========================================\n")

# ==========================================
# 5. PREPARACIÓN: MATRIZ CÉLULAS x GENES Y GENES COMUNES
# ==========================================
cat("\n[IST] 4. Preparando matriz de conteos (células x genes)...\n")
counts <- GetAssayData(input_seurat, assay = DefaultAssay(input_seurat), layer = "counts")
counts <- as(counts, "dgCMatrix")
counts_matrix_t <- t(counts)
nombres_celulas <- rownames(counts_matrix_t)

common_genes <- intersect(colnames(counts_matrix_t), rownames(reference))
reference    <- reference[common_genes, ]
counts_matrix_t <- counts_matrix_t[, common_genes]

cat("  -> Genes en común (datos y referencia):", length(common_genes), "\n")
cat("  -> Matriz final para InSituType:", nrow(counts_matrix_t), "células x", ncol(counts_matrix_t), "genes\n")

# ==========================================
# 6. CLUSTERING InSituType
# ==========================================
cat("\n[IST] 5. Ejecutando InSituType...\n")
set.seed(357)
ist_semi <- insitutype(
  x = counts_matrix_t,
  neg = neg_vector_final,
  reference_profiles = reference,
  n_clusts = n_clusts
)

ist_labels    <- ist_semi$clust
ist_confianza <- ist_semi$prob
nombres_celulas <- names(ist_semi$clust)

# ==========================================
# 7. RESUMEN DE CLASES DETECTADAS
# ==========================================
cat("\n========================================\n")
cat("[RESUMEN DE CLASIFICACIÓN]\n")
clases_detectadas <- length(unique(ist_labels))
cat(" -> Total de clases detectadas:", clases_detectadas, "\n\n")
cat(" -> Distribución de células por clase:\n")
tabla_clases <- table(ist_labels)
for (i in seq_along(tabla_clases)) {
  cat(sprintf("    - Clase %-12s : %5d células\n", names(tabla_clases)[i], tabla_clases[i]))
}
cat("========================================\n")

# ==========================================
# 8. GUARDADO DEL RESULTADO
# ==========================================
cat("\n[IST] 6. Guardando etiquetas...\n")
resultados <- data.frame(
  cell_ID      = nombres_celulas,
  ist_class    = ist_labels,
  ist_confianza = ist_confianza,
  stringsAsFactors = FALSE
)
write.table(resultados, file = out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[IST] ¡Proceso completado con éxito en:", out_txt, "!\n")

# ==========================================
# 9. CIERRE DEL LOG
# ==========================================
sink(type = "message")
sink()
close(con)
