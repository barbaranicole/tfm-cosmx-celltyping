# ==============================================================================
# 05_Figuras_QC.R
# ------------------------------------------------------------------------------
# Descripción: Figuras de control de calidad: distribución del número de
#              transcritos (nCount_RNA) y de genes (nFeature_RNA) por célula,
#              comparando el escenario original con el degradado (thin_75),
#              antes y después del control de calidad.
#
# Entradas:    Resultados_Algoritmos_v1/<escenario>_joinmerged_clustering.rds
#              Resultados_Algoritmos_v1/datos_histograma.rds (datos previos al QC)
# Salidas:     Plots/QC.png    Histogramas tras el QC (original vs. thin_75)
#              Plots/noQC.png  Histogramas antes del QC (original vs. degradado)
# Requiere:    04_Anadir_clustering.R
# Siguiente:   06_Marcadores.R
# ==============================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

# ==========================================
# 1. CARGA DE DATOS
# ==========================================
# Extraccion de metadatos
merged_rds   <- here::here("Resultados_Algoritmos_v1/original_joinmerged_clustering.rds")
merged <- readRDS(merged_rds)
meta_original <- merged@meta.data
meta_original$Version <- "Original"

merged_rds   <- here::here("Resultados_Algoritmos_v1/thin_75_joinmerged_clustering.rds")
merged_thin75 <- readRDS(merged_rds)
meta_thin <- merged_thin75@meta.data
meta_thin$Version <- "Thin75"
# Revisiones rápidas (solo se muestran en consola) 
meta_original
dim(merged)
dim(merged_thin75)

#revision de numero de genes con valores 0 
typeof(merged@assays$RNA)
counts_matrix <- merged_thin75@assays$RNA@layers$counts
sum(rownames(counts_matrix)[rowSums(counts_matrix) > 0])

# ==========================================
# 2. FIGURA QC: DATOS TRAS EL CONTROL DE CALIDAD (QC.png)
# ==========================================
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

# ==========================================
# 3. FIGURA QC: DATOS ANTES DEL CONTROL DE CALIDAD (noQC.png)
# ==========================================

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

# ==========================================
# 4. EXPLORATORIO: DENSIDADES Y COMPLEJIDAD (log10GenesPerUMI)
# ==========================================

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
