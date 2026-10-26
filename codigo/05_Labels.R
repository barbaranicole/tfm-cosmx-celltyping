library(dplyr)

# etiquetas QMDs
leiden <- readRDS(here::here("Resultados_Algoritmos_v1/temp_leiden_labels.rds"))
ist    <- readRDS(here::here("Resultados_Algoritmos_v1/temp_ist_labels.rds"))

etiquetas_original <- full_join(leiden$orig, ist$orig, by = "cell_ID")
etiquetas_thin     <- full_join(leiden$thin, ist$thin, by = "cell_ID")
resultados_finales <- full_join(etiquetas_original, etiquetas_thin, by = "cell_ID")

names(resultados_finales) <- gsub("_grupo", "", names(resultados_finales)) #quitar el sufijo "_grupo"

write.table(resultados_finales, 
            here::here("Resultados_Algoritmos_v1/clustering_etiquetado.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
