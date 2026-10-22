suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(ggplot2)
})
## --------------------------------------------- ##
#Argumentos para entrada y salida 
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) stop("Uso: Rscript 01_QC.R <input_raw.rds> <out_subset.rds> <out_neg_vector.rds>")

input_rds <- args[1]
out_rds <- args[2]
out_neg <- args[3]

## --------------------------------------------- ##
cat("[QC] Cargando datos...\n")
#Lectura de datos
datos <- readRDS(input_rds)

#actualizacion del objeto seurat de v4 a v5
datos <- UpdateSeuratObject(datos)

#assay_name <- DefaultAssay(datos) #"Nanostring"
#print(dim(datos))

## --------------------------------------------- ##
cat("[QC] Identificando sondas negativas...\n")

# identificar sondas negativas (Background)
neg_features <- grep("Neg|neg|Ctrl|System", rownames(datos), value = TRUE)

# Extraemos sus conteos y calculamos el promedio por célula
raw_neg_counts <- GetAssayData(datos, layer = "counts")[neg_features, ]
neg_vector_final <- colMeans(as.matrix(raw_neg_counts)) 

#Extraccion de sondas negativas de la matriz 
real_genes <- setdiff(rownames(datos), neg_features)
datos <- subset(datos, features = real_genes)

## --------------------------------------------- ##
cat("[QC] Aplicando filtros de calidad...\n")

#Filtros QC 
min_counts <- 20
min_features <- 20
min_area <- 100  # Para eliminar detritos muy pequeños

# Aplicar el filtrado y guardar en un nuevo objeto (o sobrescribir el actual)
datos_subset <- subset(
  x = datos, 
  subset = nCount_Nanostring >= min_counts & 
    nFeature_Nanostring >= min_features & 
    Area >= min_area)

dim(datos_subset)

## --------------------------------------------- ##
cat("[QC] Normalizando...\n")
datos_subset <- NormalizeData(datos_subset, assay = "Nanostring", normalization.method = "LogNormalize")

## --------------------------------------------- ##
cat("[QC] Guardando resultados...\n")
saveRDS(datos_subset, file = out_rds)
saveRDS(neg_vector_final, file = out_neg)
cat("[QC] ¡Listo!\n")


