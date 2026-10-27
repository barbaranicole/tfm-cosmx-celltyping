# ==============================================================================
# 06_Confusion.R
# ------------------------------------------------------------------------------
# Descripción: Compara las etiquetas de tipo celular entre escenarios, parámetros
#              y métodos mediante matrices de confusión, y calcula la
#              concordancia global entre todas las combinaciones con el índice
#              de Rand ajustado (ARI).
#
# Entradas:    Resultados_Algoritmos_v1/clustering_etiquetado.txt
# Salidas:     Resultados_Algoritmos_v1/Matriz_confusion/<col1>_vs_<col2>_combinado.png
#              Resultados_Algoritmos_v1/Matriz_confusion/ARI.png
# Requiere:    05_Labels.R
# Funciones:   guardar_matriz() y mat_confusion() de funciones.R
# ==============================================================================

library(dplyr)
library(pheatmap)
library(stringr)
library(ggplot2)
library(mclust)

# ==========================================
# 0. CONFIGURACIÓN Y LECTURA DE ETIQUETAS
# ==========================================

box::use(./funciones[mat_confusion])
box::use(./funciones[guardar_matriz])

resultados <- read.table(here::here("Resultados_Algoritmos_v1/clustering_etiquetado.txt"), 
                      header = TRUE, 
                      sep = "\t",
                      row.names = 1)

base <- here::here("Resultados_Algoritmos_v1")
base_salida <- file.path(base, "Matriz_confusion")


# ==========================================
# 1. ORIGINAL VS. DEGRADADO (mismo método y parámetro)
# ==========================================

# InSituType
guardar_matriz(resultados, "IST_original_nclusts6", "IST_thin_75_nclusts6", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_original_nclusts8", "IST_thin_75_nclusts8", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_original_nclusts10", "IST_thin_75_nclusts10", base_salida, mostrar_conteo = FALSE)

# Leiden
guardar_matriz(resultados, "Leiden_original_res0.8", "Leiden_thin_75_res0.8", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "Leiden_original_res1", "Leiden_thin_75_res1", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "Leiden_original_res1.2", "Leiden_thin_75_res1.2", base_salida, mostrar_conteo = FALSE)

# ==========================================
# 2. ENTRE PARÁMETROS DEL MISMO MÉTODO: ESCENARIO ORIGINAL
# ==========================================

guardar_matriz(resultados, "Leiden_original_res0.8", "Leiden_original_res1", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "Leiden_original_res0.8", "Leiden_original_res1.2", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "Leiden_original_res1", "Leiden_original_res1.2", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_original_nclusts6", "IST_original_nclusts8", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_original_nclusts6", "IST_original_nclusts10", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_original_nclusts8", "IST_original_nclusts10", base_salida, mostrar_conteo = FALSE)

guardar_matriz(resultados, "IST_original_nclusts10", "IST_original_nclusts8", base_salida, mostrar_conteo = FALSE)

# ==========================================
# 3. ENTRE PARÁMETROS DEL MISMO MÉTODO: ESCENARIO DEGRADADO
# ==========================================

guardar_matriz(resultados, "Leiden_thin_75_res0.8", "Leiden_thin_75_res1", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "Leiden_thin_75_res0.8", "Leiden_thin_75_res1.2", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "Leiden_thin_75_res1", "Leiden_thin_75_res1.2", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_thin_75_nclusts6", "IST_thin_75_nclusts8", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_thin_75_nclusts6", "IST_thin_75_nclusts10", base_salida, mostrar_conteo = FALSE)
guardar_matriz(resultados, "IST_thin_75_nclusts8", "IST_thin_75_nclusts10", base_salida, mostrar_conteo = FALSE)

# ==========================================
# 4. ENTRE MÉTODOS (Leiden res 1.2 vs. InSituType 10 perfiles)
# ==========================================

guardar_matriz(resultados, "Leiden_original_res1.2", "IST_original_nclusts10", base_salida, mostrar_conteo = FALSE)

# ==========================================
# 5. CONCORDANCIA GLOBAL (ARI)
# ==========================================
# 1. Crear la matriz inicial (llena de 1s para la diagonal)
columnas <- colnames(resultados)
matriz_ari <- matrix(1, nrow = length(columnas), ncol = length(columnas),
                     dimnames = list(columnas, columnas))

# 2. Calcular los valores ARI para cada par de columnas
for (i in 1:(length(columnas) - 1)) {
  for (j in (i + 1):length(columnas)) {
    # print(paste(columnas[i], columnas[j])) # Opcional: para ver el progreso
    valor_ari <- adjustedRandIndex(resultados[[columnas[i]]], resultados[[columnas[j]]])
    
    matriz_ari[i, j] <- valor_ari
    matriz_ari[j, i] <- valor_ari
  }
}

# 3. Preparar las anotaciones estéticas (Método y Escenario)
metodo <- ifelse(grepl("IST", columnas), "InSituType", "Leiden")
escenario <- ifelse(grepl("original", columnas), "Original", "Degradado")

anotaciones <- data.frame(
  Método = metodo,
  Escenario = escenario
)
# Es vital que el nombre de las filas coincida con la matriz_ari original
rownames(anotaciones) <- columnas 

# 4. Crear las etiquetas cortas (solo el hiperparámetro)
etiquetas_cortas <- gsub(".*nclusts|.*res", "", columnas)
etiquetas_cortas <- ifelse(metodo == "InSituType", 
                           paste("Perfiles:", etiquetas_cortas), 
                           paste("Res:", etiquetas_cortas))

# 5. Definir colores para las barras de anotación
colores_anotacion <- list(
  Método = c(InSituType = "#005083", Leiden = "#F26052"),
  Escenario = c(Original = "#4DBBD5", Degradado = "#7A3B69")
)

# 6. Generar el heatmap
color <- colorRampPalette(c("white", "cadetblue2", "turquoise4"))(5)
breaks <- c(0.6, 0.7, 0.8, 0.9, 1)

file_png <- paste0("ARI.png")
filename <- file.path(base_salida, file_png)

mat_plot <- pheatmap(matriz_ari, 
         display_numbers = TRUE,
         number_format = "%.3f",
         color = color,
         breaks = breaks,
         annotation_row = anotaciones,     # Agrega la barra lateral
         annotation_col = anotaciones,     # Agrega la barra superior
         annotation_colors = colores_anotacion, # Aplica los colores
         labels_row = etiquetas_cortas,    # Cambia el texto de la fila
         labels_col = etiquetas_cortas,    # Cambia el texto de la columna
         main = "Concordancia global (ARI) entre métodos y escenarios",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         # --- AJUSTES DE LEYENDA Y MÁRGENES ---
         fontsize = 10,                 # Tamaño base de fuente (ajusta si es necesario)
         fontsize_number = 12,           # Tamaño de los números dentro de las celdas
         annotation_legend = TRUE,      # Asegura que se muestre la leyenda de anotaciones
         legend = TRUE,                 # Muestra la leyenda de colores del ARI
         # fontsize_col y fontsize_row pueden ayudar a dar más espacio a la derecha
         fontsize_row = 10,
         fontsize_col = 10,
         annotation_names_row = FALSE,  # Quita los textos de abajo a la izquierda
         annotation_names_col = FALSE,  # Quita los textos de arriba a la derecha
         filename = filename,  # Le pasas la ruta aquí
         width = 12,           # Ancho en pulgadas
         height = 8           # Alto en pulgadas
         )

file_png <- paste0("ARI.png")
filename <- file.path(base_salida, file_png)

# Nota: pheatmap ya guarda ARI.png arriba (8 pulgadas de alto); este ggsave
# lo sobrescribe con 10 pulgadas de alto. La versión final es esta.

ggsave(filename = filename, 
       plot = mat_plot, 
       bg = "white",
       width = 12,      
       height = 10,      
       dpi = 300)

# ==========================================
# 6. ANÁLISIS EXPLORATORIO: ARI SIN CÉLULAS TUMORALES
# ==========================================
# Recalcula el ARI excluyendo las células etiquetadas como "Tumor" en alguna
# columna. Solo se muestra en pantalla (no se guarda) y no forma parte de los
# resultados finales.

resultados_sin_tumor <- resultados %>% filter_all(all_vars(. != 'Tumor'))
dim(resultados)
dim(resultados_sin_tumor)
columnas <- colnames(resultados)
matriz_ari <- matrix(1, nrow = length(columnas), ncol = length(columnas),
                     dimnames = list(columnas, columnas))

for (i in 1:(length(columnas) - 1)) {
  for (j in (i + 1):length(columnas)) {
    print(paste(columnas[i], columnas[j]))
    valor_ari <- adjustedRandIndex(resultados_sin_tumor[[columnas[i]]], resultados_sin_tumor[[columnas[j]]])
    
    matriz_ari[i, j] <- valor_ari
    matriz_ari[j, i] <- valor_ari
  }
}

color = colorRampPalette(c("white", "steelblue"))(6)
breaks = c(0.5, 0.6, 0.7, 0.8, 0.9, 1)

pheatmap(matriz_ari, 
         display_numbers = TRUE, # Muestra el valor numérico en cada celda
         number_format = "%.4f",
         color = color,
         breaks = breaks,
         main = "Matriz de Similitud (ARI) entre métodos de Clustering - Sin tumor")
