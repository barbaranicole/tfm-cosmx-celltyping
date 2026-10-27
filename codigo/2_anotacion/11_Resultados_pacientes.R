library(Seurat)
library(ggplot2)
library(scales)
library(tidyverse)

#### datos pre control de calidad ####

base <- here::here("Resultados_Algoritmos_v1")
base_salida <- file.path(base, "Resultados_etiquetado")

noqc_rds   <- here::here("Resultados_Algoritmos_v1/datos_histograma.rds")
merged_noqc <- readRDS(noqc_rds)
head(merged_noqc)
merged_noqc %>% group_by(dataset) %>% summarise(
  media = mean(nCount_RNA, na.rm = TRUE),
  mediana = median(nCount_RNA, na.rm = TRUE),
  cuartil_1 = quantile(nCount_RNA, probs = 0.25, na.rm = TRUE), # Q1 (25%)
  cuartil_3 = quantile(nCount_RNA, probs = 0.75, na.rm = TRUE)  # Q3 (75%)
)

merged_noqc %>% group_by(dataset) %>% summarise(
  media = mean(nFeature_RNA, na.rm = TRUE),
  mediana = median(nFeature_RNA, na.rm = TRUE),
  cuartil_1 = quantile(nFeature_RNA, probs = 0.25, na.rm = TRUE), # Q1 (25%)
  cuartil_3 = quantile(nFeature_RNA, probs = 0.75, na.rm = TRUE)  # Q3 (75%)
)

merged_noqc$paciente <- str_extract(rownames(merged_noqc), "Lung\\d+")
merged_noqc$Id <- str_extract(rownames(merged_noqc), "Lung\\d+_\\d+_\\d+")
merged_noqc <- merged_noqc %>%
  mutate(dataset = ifelse(dataset == "Originnal", "Original", dataset))
merged_noqc$dataset <- factor(merged_noqc$dataset, levels = c("Original", "Degradado"))

plot_hist_ncount <- ggplot(merged_noqc, aes(x = log(nCount_RNA+1), fill = paciente, color = paciente)) +
  geom_histogram(alpha = 0.5, bins = 50, position = "identity") +
  guides(fill=guide_legend(title="Version"), color=guide_legend(title="Version")) +
  labs(title = "Distribución de Transcritos",
       x = "nCount_RNA", y = "Count") +
  theme_classic() +
  facet_grid(rows = vars(dataset),
             cols = vars(paciente))

plot_hist_ncount

plot_hist_nft <- ggplot(merged_noqc, aes(x = nFeature_RNA, fill = paciente, color = paciente)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  guides(fill=guide_legend(title="Version"), color=guide_legend(title="Version")) +
  labs(title = "Distribución de Genes",
       x = "nFeature_RNA", y = "Count") +
  theme_classic() +
  facet_grid(rows = vars(dataset),
             cols = vars(paciente),
             scales = 'free')
plot_hist_nft

#### datos post control de calidad ####

merged_rds   <- here::here("Resultados_Algoritmos_v1/original_joinmerged_clustering.rds")
merged <- readRDS(merged_rds)
meta_original <- merged@meta.data
meta_original$Version <- "Original"

merged_rds   <- here::here("Resultados_Algoritmos_v1/thin_75_joinmerged_clustering.rds")
merged_thin75 <- readRDS(merged_rds)
meta_thin <- merged_thin75@meta.data
meta_thin$Version <- "Thin75"

columnas_qc <- c("nCount_RNA", "nFeature_RNA", "Version", "orig.ident")

meta_combinada <- bind_rows(
  Original = meta_original[, columnas_qc] %>% rownames_to_column("cell_ID"),
  Thin75 = meta_thin[, columnas_qc] %>% rownames_to_column("cell_ID"),
  .id = "Version" 
)

meta_combinada$Version <- factor(meta_combinada$Version, levels = c("Original", "Thin75"))

plot_hist_ncount <- ggplot(meta_combinada, aes(x = log(nCount_RNA+1), fill = orig.ident, color = orig.ident)) +
  geom_histogram(alpha = 0.5, bins = 50, position = "identity") +
  guides(fill=guide_legend(title="Version"), color=guide_legend(title="Version")) +
  labs(title = "Distribución de Transcritos",
       x = "nCount_RNA", y = "Count") +
  theme_classic() +
  facet_grid(rows = vars(Version),
             cols = vars(orig.ident))
plot_hist_ncount

plot_hist_nft <- ggplot(meta_combinada, aes(x = nFeature_RNA, fill = orig.ident, color = orig.ident)) +
  geom_histogram(alpha = 0.5, bins = 50) +
  guides(fill=guide_legend(title="Version"), color=guide_legend(title="Version")) +
  labs(title = "Distribución de Genes",
       x = "nFeature_RNA", y = "Count") +
  theme_classic() +
  facet_grid(rows = vars(Version),
             cols = vars(orig.ident))
plot_hist_nft

### Anadir clustering ###
resultados <- read.table(here::here("Resultados_Algoritmos_v1/clustering_etiquetado.txt"), 
                         header = TRUE, 
                         sep = "\t",
                         row.names = 1)

#extraccion de cada escenario
original <- colnames(resultados)[!str_detect(colnames(resultados), ".thin_75.")]
thin_75 <- colnames(resultados)[str_detect(colnames(resultados), ".thin_75.")]
meta_original <- resultados[, c(original)]
meta_thin_75 <- resultados[, c(thin_75)]

# pivot para unir con los datos combinados (ncount y nfeature)
meta_original <- meta_original %>%
  rownames_to_column("cell_ID") %>%
  mutate(across(-cell_ID, as.character)) %>% 
  pivot_longer(
    cols = -cell_ID,                                
    names_to = c("Metodo", "Resolucion"),             
    names_pattern = "(.*)_original_(.*)",           
    values_to = "Cluster"                          
  )
head(meta_original)

meta_thin_75 <- meta_thin_75 %>%
  rownames_to_column("cell_ID") %>%
  mutate(across(-cell_ID, as.character)) %>% 
  pivot_longer(
    cols = -cell_ID,                                
    names_to = c("Metodo", "Resolucion"),             
    names_pattern = "(.*)_thin_75_(.*)",           
    values_to = "Cluster"                          
  )

#Join 
meta_original_part <- meta_combinada %>%
  filter(Version == "Original") %>%
  left_join(meta_original, by = "cell_ID")

meta_thin_part <- meta_combinada %>%
  filter(Version == "Thin75") %>% 
  left_join(meta_thin_75, by = "cell_ID")

head(meta_thin_part)
meta_final <- bind_rows(meta_original_part, meta_thin_part)

head(meta_final)
meta_final$orig.ident <- factor(meta_final$orig.ident, levels = c("Lung5_Rep3", "Lung6", "Lung9_Rep1", "Lung12", "Lung13"))
meta_final$Metodo <- factor(meta_final$Metodo, levels = c("Leiden", "IST"))
meta_final$Resolucion <- factor(meta_final$Resolucion, levels = c("res0.8", "res1", "res1.2", "nclusts6", "nclusts8", "nclusts10"))

#### Plot ####
paleta_maestra <- c(
  # --- LINAJE LINFOIDE (Azules, pizarras y púrpuras maduros) ---
  "cel. B"  = "#4DBBD5",  
  "Plas"    = "#7A3B69", 
  "Linf. T" = "#3C5488",  
  "T4"      = "#5C6BC0",  
  "T8"      = "#0072B5", 
  "Treg"    = "#F39B7F",  
  "NK"      = "#8A4198",  
  
  # --- LINAJE MIELOIDE (Verdes, tierras y cálidos) ---
  "Mieloide"= "#91D1C2", 
  "Mac"     = "#3B8132",  
  "Mono"    = "#B2DF9A",  
  "cDc"     = "#F26052",  
  "pDC"     = "#C49A00",
  "Mast"    = "#8B5A2B",  
  "Neu"     = "#B09C85",  
  "Mac/Fib" = "#4D6C7B",  
  
  # --- ESTROMA / MICROAMBIENTE ---
  "Endo"    = "#D15396",  
  "Fib"     = "#E68AC8",  
  
  # --- TUS GRUPOS DEL RECODE (Gama principal de alto impacto) ---
  "Tumor"   = "#C22D3C",  
  "Estroma" = "#00A087",  
  "Epit"    = "#005083" 
)


# Original - comparacion entre pacientes y entre metodos 
meta_final_original <- meta_final %>%
                        filter(Version == "Original")

# Pacientes
bar_paciente_original <- ggplot(meta_final_original, aes(y = Resolucion, 
                                fill = Cluster
                                )) +
                          geom_bar() +
                          scale_fill_manual(values=paleta_maestra[unique(meta_final_original$Cluster)]) +
                          scale_x_continuous(breaks = scales::breaks_width(20000)) +
                          theme_classic() +
                          facet_grid(cols = vars(orig.ident), scales="free") +
                          labs(title = "Proporción de clústers por paciente",
                               subtitle = paste0("Original"))

filename <- file.path(base_salida, "clusters_paciente_original.png")
ggsave(filename = filename, 
       plot = bar_paciente_original, 
       bg = "white",
       width = 12,      
       height = 7,      
       dpi = 300)

#Metodos
bar_metodos_original <- ggplot(meta_final_original, aes(y = Resolucion, 
                          fill = Cluster)) +
                        geom_bar() +
                        scale_fill_manual(values=paleta_maestra[unique(meta_final_original$Cluster)]) +
                        scale_x_continuous(breaks = scales::breaks_width(50000)) +
                        theme_classic() +
                        labs(title = "Proporción de clústers por metodo",
                             subtitle = paste0("Original"))

filename <- file.path(base_salida, "clusters_metodo_original.png")
ggsave(filename = filename, 
       plot = bar_metodos_original,
       bg = "white",
       width = 10,      
       height = 7,      
       dpi = 300)

### -------------------------------------------------------------------------------- ##

# Degradado 
meta_final_thin <- meta_final %>%
  filter(Version == "Thin75")

bar_paciente_thin <- ggplot(meta_final_thin, aes(y = Resolucion, 
                                fill = Cluster)) +
                      geom_bar() +
                      scale_fill_manual(values=paleta_maestra[unique(meta_final_thin$Cluster)]) +
                      scale_x_continuous(breaks = scales::breaks_width(20000)) +
                      theme_classic() +
                      facet_grid(cols = vars(orig.ident), 
                                 scales="free") +
                      labs(title = "Proporción de clústers por paciente",
                           subtitle = paste0("Degradado (Thin75)"))

filename <- file.path(base_salida, "clusters_paciente_degradado.png")
ggsave(filename = filename, 
       plot = bar_paciente_thin, 
       bg = "white",
       width = 12,      
       height = 7,      
       dpi = 300)

bar_metodos_thin <- ggplot(meta_final_thin, aes(y = Resolucion, 
                                fill = Cluster)) +
                    geom_bar() +
                    scale_fill_manual(values=paleta_maestra[unique(meta_final_thin$Cluster)]) +
                    scale_x_continuous(breaks = scales::breaks_width(50000)) +
                    theme_classic() +
                    labs(title = "Proporción de clústers por metodo",
                         subtitle = paste0("Degradado"))

filename <- file.path(base_salida, "clusters_metodo_degradado.png")
ggsave(filename = filename, 
       plot = bar_metodos_thin,
       bg = "white",
       width = 10,      
       height = 7,      
       dpi = 300)

### -------------------------------------------------------------------------------- ##

#Escenarios 
bar_escenarios <- ggplot(meta_final, aes(y = Resolucion, 
                                fill = Cluster)) +
                  geom_bar() +
                  scale_fill_manual(values=paleta_maestra[unique(meta_final$Cluster)]) +
                  scale_x_continuous(breaks = scales::breaks_width(50000)) +
                  theme_classic() + 
                  facet_grid(cols = vars(Version), scales="free") +
                  labs(title = "Proporción de clústers por método",
                       subtitle = paste0("Escenarios"))

filename <- file.path(base_salida, "clusters_escenarios.png")
ggsave(filename = filename, 
       plot = bar_escenarios,
       bg = "white",
       width = 12,      
       height = 7,      
       dpi = 300)
### -------------------------------------------------------------------------------- ##


library(dplyr)
library(tidyr)
library(clustree)

# 1. Filtramos y preparamos el dataframe
df_comparacion <- meta_final %>%
  filter(
    Version == "Original",
    (Metodo == "Leiden" & Resolucion == "res1.2") | 
      (Metodo == "IST" & Resolucion == "nclusts10")
  ) %>%
  # CORRECCIÓN: Usar estrictamente números después del prefijo "k_"
  mutate(Nivel_Clustree = case_when(
    Metodo == "Leiden" ~ "k_1",
    Metodo == "IST" ~ "k_2"
  )) %>%
  select(cell_ID, nCount_RNA, nFeature_RNA, orig.ident, Nivel_Clustree, Cluster) %>%
  pivot_wider(
    names_from = Nivel_Clustree,
    values_from = Cluster
  )

# 2. Generamos el gráfico
# Le decimos a clustree que busque las columnas que empiezan por "k_"
clustree_comp <- clustree(df_comparacion, prefix = "k_", layout = "sugiyama",) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    #legend.key.height = unit(0.2, "cm"),
    legend.position = "right" ) +
  ggplot2::guides(
    color = "none")

ruta_clustree_comp <- file.path(base_salida, paste0("Leiden_IST_clustree_comparativo.png"))
ggsave(filename = ruta_clustree_comp, 
       plot = clustree_comp, 
       width = 10, 
       height = 6, 
       dpi = 300, 
       bg = "white")
