suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  stop("Error. Uso: Rscript 04_CellAssign.R <input.rds> <reference.csv> <out_labels.txt>")
}

in_rds      <- args[1]
ref_path <- args[2]
out_txt     <- args[4]

## --------------------------------------------- ##
cat("[CellAssign] Configurando entorno de Python (Conda)...\n")
Sys.setenv(RETICULATE_PYTHON = "/opt/homebrew/anaconda3/envs/r-reticulate/bin/python")
Sys.setenv(RETICULATE_CONDA = "/opt/homebrew/anaconda3/bin/conda")
library(reticulate)
use_condaenv("r-reticulate", required = TRUE)

tryCatch({
  tf <- import("tensorflow")
  cat("             TensorFlow cargado correctamente:", tf$version$VERSION, "\n")
}, error = function(e) stop("Error crítico: No se pudo cargar TensorFlow."))

library(cellassign)

## --------------------------------------------- ##
generate_marker_matrix_logfc <- function(ref_matrix, logfc_threshold = 1, min_expr = 1) {
  
  cell_types <- colnames(ref_matrix)
  genes <- rownames(ref_matrix)
  
  # Matriz para guardar los resultados
  fc_matrix <- matrix(0, nrow = length(genes), ncol = length(cell_types))
  rownames(fc_matrix) <- genes
  colnames(fc_matrix) <- cell_types
  
  # Matriz binaria de salida (la Rho que necesita CellAssign)
  binary_matrix <- matrix(0, nrow = length(genes), ncol = length(cell_types))
  rownames(binary_matrix) <- genes
  colnames(binary_matrix) <- cell_types
  
  for (ctype in cell_types) {
    # Definir el grupo "Target" y el grupo "Resto"
    target_expr <- ref_matrix[, ctype]
    
    # El resto son todas las columnas MENOS la actual
    other_types <- ref_matrix[, colnames(ref_matrix) != ctype, drop = FALSE]
    rest_expr <- rowMeans(other_types) # Promedio del resto
    
    # Log2 Fold Change añadir +1 (pseudocount) para evitar log(0)
    log_fc <- log2(target_expr + 1) - log2(rest_expr + 1)
    
    # Guardamos el FC por si quieres inspeccionarlo
    fc_matrix[, ctype] <- log_fc
    
    # Binarización
    # Fold Change positivo alto y tener una expresión mínima absoluta (Evitar marcadores de ruido)
    es_marcador <- (log_fc > logfc_threshold) & (target_expr > min_expr)
    binary_matrix[es_marcador, ctype] <- 1
  }
  
  # Eliminamos genes que sean 0 para todos (no informativos)
  keep_genes <- rowSums(binary_matrix) > 0
  
  final_matrix <- binary_matrix[keep_genes, ]
  
  message("Hecho. Se seleccionaron ", nrow(final_matrix), " genes marcadores.")
  return(final_matrix)
}

## --------------------------------------------- ##
cat("[CellAssign] Cargando objetos...\n")

input_seurat <- readRDS(input_rds)

## --------------------------------------------- ##
cat("[CellAssign] Preparando matriz de conteos y size factors...\n")
# Extraemos dispersa y transponemos (Células x Genes)
counts <- GetAssayData(input_seurat, assay = DefaultAssay(input_seurat), layer = "counts")
counts_matrix_t <- t(counts)
nombres_celulas_totales <- rownames(counts_matrix_t)

# Calculamos factores de tamaño. Para 500k células usar rowSums en dispersa es rapidísimo.
size_factors <- rowSums(counts_matrix_t) / mean(rowSums(counts_matrix_t))

## --------------------------------------------- ##
cat("[CellAssign] Procesando Perfil de Referencia...\n")
reference <- as.matrix(read.csv(ref_path, row.names = 1))
common_genes <- intersect(colnames(counts_matrix_t), rownames(reference))
reference <- reference[common_genes, ]

rho_matrix <- generate_marker_matrix_logfc(reference,
                                           logfc_threshold = 1.0,
                                           min_expr = 5)

# Filtramos la matriz de expresión solo a los marcadores
genes_ca  <- rownames(rho_matrix)
counts_ca <- counts_matrix_t[, genes_cak]

## --------------------------------------------- ##
cat("Ejecutando algoritmo CellAssign (6k)...\n")

set.seed(357)
fit <- cellassign(exprs_obj = counts_ca, 
                  marker_gene_info = rho_matrix, 
                  s = size_factors, 
                  learning_rate = 1e-2, 
                  verbose = FALSE)

# Extracción segura: CellAssign devuelve un vector con los labels en el mismo orden de las filas
ca_labels <- fit$cell_type
names(ca_labels) <- rownames(counts_ca)

## --------------------------------------------- ##
cat("[CellAssign] Extrayendo etiquetas...\n")

# Extraemos las etiquetas usando la lista maestra de células para evitar cualquier desfase
resultados_ca <- data.frame(
  cell_ID = nombres_celulas_totales,
  cellassign_type_6k = as.character(ca_labels_6k[nombres_celulas_totales]),
  cellassign_type_1k = as.character(ca_labels_1k[nombres_celulas_totales]),
  stringsAsFactors = FALSE
)

cat("[CellAssign] Guardando TXT...\n")
write.table(resultados_ca, file = out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[CellAssign] ¡Listo!")
