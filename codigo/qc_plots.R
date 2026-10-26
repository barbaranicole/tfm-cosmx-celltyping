library(ggplot2)
library(dplyr)
library(patchwork)

# Extraccion de metadatos
merged_rds   <- here::here("Resultados_Algoritmos_v1/original_joinmerged_clustering.rds")
merged <- readRDS(merged_rds)
meta_original <- merged@meta.data
meta_original$Version <- "Original"

merged_rds   <- here::here("Resultados_Algoritmos_v1/thin_75_joinmerged_clustering.rds")
merged_thin75 <- readRDS(merged_rds)
meta_thin <- merged_thin75@meta.data
meta_thin$Version <- "Thin75"

meta_original
dim(merged)
dim(merged_thin75)

#revision de numero de genes con valores 0 
typeof(merged@assays$RNA)
counts_matrix <- merged_thin75@assays$RNA@layers$counts
sum(rownames(counts_matrix)[rowSums(counts_matrix) > 0])

columnas_qc <- c("nCount_RNA", "nFeature_RNA", "Version")

meta_combinada <- rbind(
  meta_original[, columnas_qc],
  meta_thin[, columnas_qc]
)

meta_combinada$Version <- factor(meta_combinada$Version, levels = c("Original", "Thin75"))

#https://hbctraining.github.io/scRNA-seq/lessons/04_SC_quality_control.html

plot_hist_ncount <- ggplot(meta_combinada, aes(x = nCount_RNA, fill = Version, color = Version)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  labs(title = "Distribución de Transcritos",
       x = "nCount_RNA", y = "Count") +
  theme_classic()

plot_hist_nft <- ggplot(meta_combinada, aes(x = nFeature_RNA, fill = Version, color = Version)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  labs(title = "Distribución de Genes",
       x = "nFeature_RNA", y = "Count") +
  theme_classic()

plots <- (plot_hist_ncount + plot_hist_nft) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

dir_plots <- here::here("Plots") 
nombre_qc <- paste0("QC", ".png")
ruta_qc <- file.path(dir_plots, nombre_qc)
png(ruta_qc, width = 6000, height = 3000, res = 300) 
print(plots)
dev.off()
cat("-> QC guardada:", ruta_qc, "\n")

noqc_rds   <- here::here("Resultados_Algoritmos_v1/datos_histograma.rds")
merged_noqc <- readRDS(noqc_rds)

plot_hist_ncount <- ggplot(merged_noqc, aes(x = nCount_RNA, fill = dataset, color = dataset)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  guides(fill=guide_legend(title="Version"), color=guide_legend(title="Version")) +
  labs(title = "Distribución de Transcritos",
       x = "nCount_RNA", y = "Count") +
  theme_classic()

plot_hist_nft <- ggplot(merged_noqc, aes(x = nFeature_RNA, fill = dataset, color = dataset)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  guides(fill=guide_legend(title="Version"), color=guide_legend(title="Version")) +
  labs(title = "Distribución de Genes",
       x = "nFeature_RNA", y = "Count") +
  theme_classic()

plots <- (plot_hist_ncount + plot_hist_nft) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

dir_plots <- here::here("Plots") 
nombre_qc <- paste0("noQC", ".png")
ruta_qc <- file.path(dir_plots, nombre_qc)
png(ruta_qc, width = 6000, height = 3000, res = 300) 
print(plots)
dev.off()
cat("-> QC guardada:", ruta_qc, "\n")

### Clustree 

table(merged@meta.data$IST_original_nclusts6, useNA = "always")
table(merged_thin75@meta.data$IST_thin_75_nclusts6,  useNA = "always")

celulas_comunes <- intersect(rownames(merged@meta.data), rownames(merged_thin75@meta.data))
df_comparacion <- data.frame(
  cell_ID = celulas_comunes,
  comp_1 = merged@meta.data[celulas_comunes, "IST_original_nclusts6"],
  comp_2 = merged_thin75@meta.data[celulas_comunes, "IST_thin_75_nclusts6"] # O como se llame en el thin75
)

diccionario_cortos <- c(
  "B-cell"         = "B",
  "endothelial"    = "Endo",
  "fibroblast"     = "Fib",
  "macrophage"     = "Mac",
  "mast"           = "Mast",
  "mDC"            = "mDC",
  "monocyte"       = "Mono",
  "neutrophil"     = "Neu",
  "NK"             = "NK",
  "pDC"            = "pDC",
  "plasmablast"    = "Plas",
  "T CD4 memory"   = "T4m",
  "T CD4 naive"    = "T4n",
  "T CD8 memory"   = "T8m",
  "T CD8 naive"    = "T8n",
  "Treg"           = "Treg"
)

traducir_etiquetas <- function(vector_original) {
  traduccion <- unname(diccionario_cortos[vector_original])
  vector_final <- ifelse(is.na(traduccion), vector_original, traduccion)
  vector_final[vector_final == "undefined" | is.na(vector_final)] <- "Unk"
  
  return(vector_final)
}

df_comparacion$comp_1 <- traducir_etiquetas(df_comparacion$comp_1)
df_comparacion$comp_2 <- traducir_etiquetas(df_comparacion$comp_2)

plot_clustree_comparativo <- clustree(df_comparacion, prefix = "comp_") +
  labs(title = "Estabilidad del Clustering: Original vs Thin75",
       subtitle = "InSituType (k=6)")

print(plot_clustree_comparativo)


#Plots de densidad
plot_densidad_ncount <- ggplot(meta_combinada, aes(x = nCount_RNA, fill = Version, color = Version)) +
  geom_density(aes(y= after_stat(count)), alpha = 0.5) +
  labs(title = "Distribución de Transcritos",
       x = "nCount_RNA", y = "Densidad") +
  theme(plot.title = element_text(face = "bold")) +
  theme_classic()

plot_densidad_nft <- ggplot(meta_combinada, aes(x = nFeature_RNA, fill = Version, color = Version)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribución de Genes",
       x = "nFeature_RNA", y = "Densidad") +
  theme_classic()

plot_densidad_ncount + plot_densidad_nft


#https://sydneybiox.github.io/SingleCellPlus/qc.html
#Plot conjunto
# media de nCount_RNA y nFeature para Original y Thin75
meta_medias <- meta_combinada |>
  group_by(Version) |>
  summarise(mean_nCount = mean(nCount_RNA, na.rm = TRUE), mean_nFtr = mean(nFeature_RNA, na.rm = TRUE))

#nCount
plot_nCount <- ggplot(meta_combinada, aes(x = nCount_RNA, colour = Version, fill = Version)) +
  # Histograma y densidad
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white", bins = 30) +
  geom_density(alpha = 0.4) +
  geom_vline(data = meta_medias, #media
             aes(xintercept = mean_nCount),
             colour = "blue", linetype = "dashed", linewidth = 1.2) +
  
  labs(title = "Distribución de Transcritos: Original vs Thin75",
       x = "Total de Transcritos (nCount_RNA)",
       y = "Densidad") +
  
  # paneles
  facet_wrap(~ Version, ncol = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "lightgray"), # Titulo de cada plot
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

print(plot_nCount)

#nFeature
plot_nFeature <- ggplot(meta_combinada, aes(x = nFeature_RNA, colour = Version, fill = Version)) +
  # Histograma y densidad
  geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white", bins = 30) +
  geom_density(alpha = 0.4) +
  geom_vline(data = meta_medias, #media
             aes(xintercept = mean_nFtr),
             colour = "blue", linetype = "dashed", linewidth = 1.2) +
  
  labs(title = "Distribución de Transcritos: Original vs Thin75",
       x = "Total de Transcritos (nCount_RNA)",
       y = "Densidad") +
  
  # paneles
  facet_wrap(~ Version, ncol = 2) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.background = element_rect(fill = "lightgray"), # Titulo de cada plot
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "none"
  )

print(plot_nFeature)

meta_original$log10GenesPerUMI <- log10(meta_original$nFeature_RNA) / log10(meta_original$nCount_RNA)
meta_thin$log10GenesPerUMI <- log10(meta_thin$nFeature_RNA) / log10(meta_thin$nCount_RNA)
meta_combinada$log10GenesPerUMI <- c(meta_original$log10GenesPerUMI, meta_thin$log10GenesPerUMI)

plot_densidad <- ggplot(meta_combinada, aes(x = log10GenesPerUMI, fill = Version, color = Version)) +
  geom_density(alpha = 0.5) +
  theme_classic() +
  labs(title = "Comparación de Complejidad",
       x = "log10GenesPerUMI", y = "Densidad") +
  theme(plot.title = element_text(face = "bold"))
