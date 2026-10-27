# ==============================================================================
# 05_Labels.R
# ------------------------------------------------------------------------------
# Descripción: Une las etiquetas de tipo celular asignadas en los reportes de
#              Leiden e InSituType (ambos escenarios) en una sola tabla. Elimina
#              el sufijo "_grupo" de los nombres de columna.
#
# Entradas:    Resultados_Algoritmos_v1/temp_leiden_labels.rds
#              Resultados_Algoritmos_v1/temp_ist_labels.rds
# Salidas:     Resultados_Algoritmos_v1/clustering_etiquetado.txt
#              Columnas: cell_ID y una columna por método, escenario y parámetro
#              con el tipo celular asignado.
# Requiere:    Haber renderizado 04_Resultados_Leiden.qmd y 04_Resultados_IST.qmd
# Siguiente:   06_Confusion.R, resultados_pacientes.R
# ==============================================================================

library(dplyr)

# ==========================================
# 1. CARGA DE LAS ETIQUETAS DE CADA MÉTODO
# ==========================================

# etiquetas manuales extraidas de los QMDs
leiden <- readRDS(here::here("Resultados_Algoritmos_v1/temp_leiden_labels.rds"))
ist    <- readRDS(here::here("Resultados_Algoritmos_v1/temp_ist_labels.rds"))

# ==========================================
# 2. UNIÓN POR cell_ID Y GUARDADO
# ==========================================

etiquetas_original <- full_join(leiden$orig, ist$orig, by = "cell_ID")
etiquetas_thin     <- full_join(leiden$thin, ist$thin, by = "cell_ID")
resultados_finales <- full_join(etiquetas_original, etiquetas_thin, by = "cell_ID")

names(resultados_finales) <- gsub("_grupo", "", names(resultados_finales)) #quitar el sufijo "_grupo"

write.table(resultados_finales, 
            here::here("Resultados_Algoritmos_v1/clustering_etiquetado.txt"), 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
