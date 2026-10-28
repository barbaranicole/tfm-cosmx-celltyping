# ==============================================================================
# funciones.R
# ------------------------------------------------------------------------------
# Descripción: Funciones auxiliares para el etiquetado manual de clusters y la
#              comparación entre escenarios (original vs. degradado).
#
# Uso:         - 07_Resultados_Leiden.qmd y 08_Resultados_InSituType.qmd la
#                cargan con source("funciones.R").
#              - 10_Matrices_confusion.R importa solo mat_confusion() y
#                guardar_matriz() con box::use() (marcadas con #' @export).
# Requiere:    Las funciones usadas desde los qmd asumen que ya están cargados
#              Seurat, ggplot2, patchwork, dplyr, tidyr, knitr, kableExtra y DT.
# ==============================================================================

box::use(
  stringr[str_match, str_detect, str_remove], # funciones de texto
  pheatmap[pheatmap],                        
  grDevices[colorRampPalette],
  ggplot2[ggsave],
  gridExtra[arrangeGrob],
  grid[textGrob, gpar]
)

# ==========================================
# 1. DICCIONARIO Y PALETA DE COLORES, encima de diccionario_cortos
# ==========================================
# Abreviaturas de los tipos celulares del perfil IO_1k (usadas en tablas y figuras).
diccionario_cortos <- c(
  "B.cell"                        = "B-cell",
  "Conventional.dendritic.cell"   = "cDc",
  "Endothelial"                   = "Endo",
  "Fibroblast"                    = "Fib",
  "Macrophage"                    = "Mac",
  "Mast.cell"                     = "Mast",
  "Monocyte"                      = "Mono",
  "Neutrophil"                    = "Neu",
  "NK.cell"                       = "NK",
  "Plasmacytoid.dendritic.cell"   = "pDC",
  "Myeloid.dendritic.cell"        = "mDC",
  "Plasma"                        = "Plas",
  "Plasmablast"                   = "Plasb",
  "T.cell.CD4"                    = "T4",
  "T.cell.CD8"                    = "T8",
  "T.cell.regulatory"             = "Treg",
  "undefined"                     = "*"
)

# Abrevia etiquetas usando diccionario_cortos. Las que no están en el
# diccionario se mantienen; NA y "undefined" pasan a "*".
#    - vector_original: vector de etiquetas
# Devuelve: vector de caracteres con las etiquetas abreviadas.
traducir_etiquetas <- function(vector_original) {
  vector_original[is.na(vector_original)] <- "*"
  traduccion <- unname(diccionario_cortos[vector_original])
  vector_final <- ifelse(is.na(traduccion), vector_original, traduccion)
  vector_final[vector_final == "undefined" | is.na(vector_final)] <- "*"
  return(vector_final)
}

# Colores fijos por tipo celular / grupo, para que coincidan entre figuras.
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
  
  # --- GRUPOS DEL RECODE (Gama principal de alto impacto) ---
  "Tumor"   = "#C22D3C",  
  "Estroma" = "#00A087",  
  "Epit"    = "#005083" 
)

# ==========================================
# 2. FIGURAS Y TABLAS PARA LOS INFORMES (07 y 08)
# ==========================================
# UMAP coloreado por una columna de clusters usando paleta_maestra. Las etiquetas
# sin color asignado reciben colores extra de scales::hue_pal().
#    - obj: Seurat
#    - col_cluster: columna a usar como identidad
#    - mi_base_size: tamaño base del texto
#    - pt_leyenda: tamaño de los puntos en la leyenda
# Devuelve: objeto ggplot.
plot_umap_ist <- function(obj, col_cluster, mi_base_size = 12, pt_leyenda = 3) {
  Idents(obj) <- col_cluster
  etiquetas_presentes <- levels(as.factor(obj@meta.data[[col_cluster]]))
  colores_usar <- paleta_maestra
  faltantes <- setdiff(etiquetas_presentes, names(colores_usar))
  if(length(faltantes) > 0) {
    colores_extra <- scales::hue_pal()(length(faltantes)) # Genera colores extra
    names(colores_extra) <- faltantes
    colores_usar <- c(colores_usar, colores_extra)
  }
  
  DimPlot(obj, cols = colores_usar, reduction = "umap", group.by = col_cluster) +
    theme_classic(base_size = mi_base_size) +
    theme(legend.margin = margin(t = 0, r = 10, b = 0, l = -10)) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = pt_leyenda)))
}

# Tabla muestra (orig.ident) x cluster con "porcentaje por fila (conteo)".
#    - obj: Seurat
#    - col_cluster: columna de clusters
# Devuelve: tabla kable con estilo.
tabla_conteos_ist <- function(obj, col_cluster) {
  tabla_base <- table(obj$orig.ident, obj[[col_cluster]][,1])
  tabla_prop <- prop.table(tabla_base, margin = 1) * 100
  matriz_texto <- matrix(
    sprintf("%.2f%% (%d)", tabla_prop, tabla_base),
    nrow = nrow(tabla_base),
    ncol = ncol(tabla_base),
    dimnames = dimnames(tabla_base) # Recuperamos los nombres de filas y columnas
  )
  wide <- as.data.frame.matrix(matriz_texto)
  kable(wide, caption = paste("Distribución Celular: Porcentaje (Conteo) -", col_cluster)) |>
    kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))
}

# Lee la tabla de marcadores (06_Marcadores.R), filtra avg_log2FC > 1 y
# p_val_adj < 0.05, y ordena por avg_log2FC dentro de cada cluster.
#    - ruta_markers: ruta del archivo de marcadores
# Devuelve: tabla interactiva (DT) con un cluster por columna.
resumen_marcadores_ist <- function(ruta_markers) {
  markers <- read.delim(ruta_markers)
  markers <- markers[(markers$avg_log2FC > 1) & (markers$p_val_adj < 0.05), ]
  markers$cluster <- traducir_etiquetas(as.character(markers$cluster))
  markers_sort <- markers |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC), .by_group = TRUE) |>
    mutate(rank = row_number()) |>
    ungroup() |>
    select(rank, cluster, gene) |>
    pivot_wider(names_from = cluster, values_from = gene) |>
    arrange(rank)

  datatable(markers_sort,
            extensions = 'Buttons',
            rownames = FALSE,
            class = 'cell-border stripe compact', #Emula "striped" y "condensed"
            options = list(
              dom = 'Bfrtip',
              buttons = c('copy', 'excel', 'csv'),
              pageLength = 20,
              lengthMenu = c(5, 10, 15, 30), 
              autoWidth = TRUE,
              scrollX = TRUE # Vital por si tienes muchos clústeres y la tabla se hace muy ancha
            ),
            caption = "Genes Marcadores por Clúster")
}

# Heatmap de la expresión media por cluster (AverageExpression), escalada por
# fila, con todos los genes ordenados por varianza.
#    - obj: Seurat
#    - col_cluster: columna a usar como identidad
#    - titulo: título del gráfico (por defecto "Perfil promedio - <col_cluster>")
#    - ruta_archivo: si se indica, pheatmap guarda la figura en esa ruta
perfil_heatmap_ist <- function(obj, col_cluster, titulo = NULL, ruta_archivo = NA) {
  Idents(obj) <- col_cluster
  perfiles_seurat <- AverageExpression(obj, assays = "RNA", return.seurat = FALSE)
  mat <- as.matrix(perfiles_seurat$RNA)
  var_genes <- apply(mat, 1, var)
  genes_top <- names(sort(var_genes, decreasing = TRUE))
  mat_filt <- mat[genes_top, ]
  
  pheatmap(mat_filt,
           scale = "row",
           cluster_cols = TRUE,
           cluster_rows = TRUE,
           color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
           main = titulo %||% paste("Perfil promedio -", col_cluster),
           show_rownames = FALSE,
           fontsize_row = 8,
           fontsize_col = 12,
           filename = ruta_archivo) # La clave para que se guarde correctamente
}

# Boxplots de Max.PanCK y Mean.PanCK por cluster, ordenados por mediana.
# PanCK es un marcador epitelial que ayuda a identificar clusters tumorales.
#    - obj: Seurat
#    - col_cluster: columna de clusters
# Devuelve: dos gráficos combinados con patchwork.
boxplots_panck_ist <- function(obj, col_cluster) {
  meta <- obj@meta.data
  p1 <- ggplot(meta, aes(x = reorder(.data[[col_cluster]], Max.PanCK, median),
                         y = Max.PanCK)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    theme(plot.title=element_text(size=10),
          axis.title=element_text(size=8),
          axis.text.x = element_text(size=5, angle = 90, hjust = 1)) +
    labs(title = "Max PanCK", x = "Cluster")
  
  p2 <- ggplot(meta, aes(x = reorder(.data[[col_cluster]], Mean.PanCK, median),
                         y = Mean.PanCK)) +
    geom_boxplot(outlier.shape = NA) +
    theme_minimal() +
    theme(plot.title=element_text(size=10),
          axis.title=element_text(size=8),
          axis.text.x = element_text(size=5, angle = 90, hjust = 1)) +
    labs(title = "Mean PanCK", x = "Cluster")
  p1 + p2
}

# Clustree entre las etiquetas del escenario original (comp_1) y del degradado
# (comp_2), emparejadas por cell_ID.
#    - obj_orig, obj_thin: objetos Seurat original y degradado
#    - col_orig, col_thin: columnas de etiquetas a comparar
#    - k_val: valor que se muestra en el subtítulo
#    - prop_filter: oculta las conexiones por debajo de esa proporción
# Nota: el subtítulo está fijo como "InSituType (k=...)", por lo que las
#       figuras de 07_Resultados_Leiden.qmd lo muestran aunque k_val sea
#       la resolución de Leiden. Se deja así para no alterar las figuras
#       publicadas.
plot_clustree_comparativo_ist <- function(obj_orig, obj_thin, col_orig, col_thin, k_val, prop_filter = 0.1) {
  df_original <- data.frame(
    cell_ID = rownames(obj_orig@meta.data),
    comp_1  = as.character(obj_orig@meta.data[[col_orig]])
  )

  df_thin <- data.frame(
    cell_ID = rownames(obj_thin@meta.data),
    comp_2  = as.character(obj_thin@meta.data[[col_thin]])
  )

  df_comparacion <- dplyr::full_join(df_original, df_thin, by = "cell_ID")
  df_comparacion$comp_1 <- traducir_etiquetas(df_comparacion$comp_1)
  df_comparacion$comp_2 <- traducir_etiquetas(df_comparacion$comp_2)
  p <- clustree::clustree(df_comparacion, prefix = "comp_", layout = "sugiyama", prop_filter = prop_filter) +
    ggplot2::labs(title = "Estabilidad del Clustering: Original vs Degradado",
                  subtitle = paste0("InSituType (k=", k_val, ")")) +
    ggplot2::scale_size_continuous(range = c(2, 8), breaks = c(50000, 150000)) +
    ggplot2::theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      legend.position = "right",
      legend.box.spacing = unit(0.2, "cm")) + 
    ggplot2::guides(
      color = "none"
      #edge_colour = ggplot2::guide_colorbar(barheight = ggplot2::unit(3, "cm"), 
      #                                      barwidth = ggplot2::unit(0.4, "cm"))
    )

  p
}

# ==========================================
# 3. RECODIFICACIÓN Y GUARDADO DE FIGURAS
# ==========================================
# Agrupa clusters en categorías (p. ej. Tumor, Estroma) en una nueva columna;
# los clusters que no están en `mapa` conservan su etiqueta.
#    - col_cluster: nombre de la columna original (ej. "IST_original_nclusts6")
#    - nueva_col: nombre de la nueva columna (ej. "IST_original_nclusts6_grupo")
#    - mapa: named vector, p.ej. c("a"="Tumor", "b"="Inmune*")

recode_cluster_grupos <- function(obj, col_cluster, nueva_col, mapa) {
  labels_orig <- as.character(obj@meta.data[[col_cluster]])
  labels_new  <- labels_orig
  es_mapeable <- labels_orig %in% names(mapa)
  labels_new[es_mapeable] <- mapa[labels_new[es_mapeable]]
  obj@meta.data[[nueva_col]] <- factor(labels_new)
  obj
}

# Guardar un UMAP a disco usando plot_umap_ist ya definido.
#    - obj: Seurat
#    - col_cluster: columna a usar como identidad
#    - ruta_png: ruta completa del archivo PNG
guardar_umap_ist <- function(obj, col_cluster, ruta_png) {
  png(ruta_png, width = 1400, height = 1200, res = 400)
  p <- plot_umap_ist(obj, col_cluster, mi_base_size = 4, pt_leyenda = 1.5)
  print(p)
  dev.off()
}

# Guardar boxplots de PanCK a disco (usa boxplots_panck_ist).
#    - obj: Seurat
#    - col_cluster: columna a usar como identidad
#    - ruta_png: ruta completa del archivo PNG
guardar_panck_ist <- function(obj, col_cluster, ruta_png) {
  png(ruta_png, width = 1600, height = 1000, res = 300)
  print(boxplots_panck_ist(obj, col_cluster))
  dev.off()
}

# ==========================================
# 4. MATRICES DE CONFUSIÓN 
# ==========================================
# Matrices de confusión entre dos columnas de etiquetas: porcentaje por fila
# en ambos sentidos y, opcionalmente, conteos absolutos. Título y ejes se
# deducen del nombre de las columnas (<Leiden|IST>_<original|thin_75>_<param>).
#    - obj: data.frame / metadata con ambas columnas
#    - col1, col2: columnas a comparar
#    - titulo: título manual (si es NULL se genera automáticamente)
#    - mostrar_conteo: añade la matriz de conteos
# Devuelve: grob combinado (gridExtra::arrangeGrob).
#' @export
mat_confusion <- function(obj, col1, col2, titulo = NULL, mostrar_conteo = TRUE) {
  
  # --- 1. LÓGICA DE EXTRACCIÓN Y TÍTULOS ---
  parse_col <- function(x) {
    partes <- str_match(x, "^(Leiden|IST)_(original|thin_75)_(.*)$") 
    metodo <- ifelse(partes[2] == "IST", "InSituType", "Leiden")
    version <- ifelse(partes[3] == "original", "Original", "Degradado")
    res_raw <- partes[4]
    
    if (str_detect(res_raw, "res")) {
      res_fmt <- paste0("res = ", str_remove(res_raw, "res"))
    } else if (str_detect(res_raw, "nclusts")) {
      res_fmt <- paste0(str_remove(res_raw, "nclusts"), " perfiles libres")
    } else {
      res_fmt <- res_raw
    }
    return(list(metodo = metodo, version = version, res = res_fmt))
  }
  
  p1 <- parse_col(col1)
  p2 <- parse_col(col2)
  
  # Determinamos el título principal si no se pasa uno manual
  if (is.null(titulo)) {
    if (p1$metodo == p2$metodo && p1$res == p2$res && p1$version != p2$version) {
      titulo_auto <- sprintf("Comparación Escenarios: %s (%s)", p1$metodo, p1$res)
    } else if (p1$metodo == p2$metodo && p1$res != p2$res && p1$version == p2$version) {
      titulo_auto <- sprintf("Comparación Hiperparámetros: %s (%s)", p1$metodo, p1$version)
    } else {
      titulo_auto <- sprintf("Comparación Cruzada: %s vs %s", p1$metodo, p2$metodo)
    }
  } else {
    titulo_auto <- titulo
  }
  
  # --- 2. LÓGICA DE ETIQUETAS DINÁMICAS PARA LOS EJES ---
  # Caso A: Mismo método e hiperparámetro -> Comparamos escenarios (Original vs Degradado)
  if (p1$metodo == p2$metodo && p1$res == p2$res) {
    label1 <- p1$version
    label2 <- p2$version
  } 
  # Caso B: Mismo método y escenario -> Comparamos hiperparámetros (res 0.8 vs 1.2)
  else if (p1$metodo == p2$metodo && p1$version == p2$version) {
    label1 <- p1$res
    label2 <- p2$res
  } 
  # Caso C: Diferentes métodos (Caso borde) -> Método + Hiperparámetro
  else {
    label1 <- paste(p1$metodo, p1$res)
    label2 <- paste(p2$metodo, p2$res)
  }
  
  # --- 3. CÁLCULO DE LAS MATRICES ---
  mat_counts <- table(Filas = obj[[col1]], Columnas = obj[[col2]])
  
  mat_1_vs_2 <- prop.table(mat_counts, margin = 1) * 100
  mat_1_vs_2 <- round(as.matrix(mat_1_vs_2), 1)
  
  mat_counts_inv <- table(Filas = obj[[col2]], Columnas = obj[[col1]])
  mat_2_vs_1 <- prop.table(mat_counts_inv, margin = 1) * 100
  mat_2_vs_1 <- round(as.matrix(mat_2_vs_1), 1)
  
  # --- 4. GENERACIÓN DE HEATMAPS (Sin subtítulos internos) ---
  col_pct <- colorRampPalette(c("white", "cadetblue2", "turquoise4"))(11)
  breaks_pct <- seq(0, 100, 10)
  
  # main = "" elimina el subtitulo de normalización interno de pheatmap
  hm1 <- pheatmap(mat_1_vs_2, display_numbers = TRUE, number_format = "%.1f",
                  color = col_pct, breaks = breaks_pct,
                  cluster_rows = FALSE, cluster_cols = FALSE, 
                  main = "", silent = TRUE) 
  
  hm2 <- pheatmap(mat_2_vs_1, display_numbers = TRUE, number_format = "%.1f",
                  color = col_pct, breaks = breaks_pct,
                  cluster_rows = FALSE, cluster_cols = FALSE, 
                  main = "", silent = TRUE)
  
  # --- 5. COMBINACIÓN CON ETIQUETAS PERSONALIZADAS ---
  formato_ejes <- gpar(fontsize = 12, fontface = "bold")
  
  plot_hm1 <- arrangeGrob(hm1$gtable, 
                          left = textGrob(label1, rot = 90, gp = formato_ejes),
                          bottom = textGrob(label2, gp = formato_ejes))
  
  plot_hm2 <- arrangeGrob(hm2$gtable, 
                          left = textGrob(label2, rot = 90, gp = formato_ejes),
                          bottom = textGrob(label1, gp = formato_ejes))
  
  if (mostrar_conteo) {
    col_cnt <- colorRampPalette(c("#F7FBFF", "lightskyblue2", "steelblue4"))(50)
    hm3 <- pheatmap(mat_counts, display_numbers = TRUE, number_format = "%.0f",
                    color = col_cnt, cluster_rows = FALSE, cluster_cols = FALSE, 
                    main = "", silent = TRUE)
    
    plot_hm3 <- arrangeGrob(hm3$gtable, 
                            left = textGrob(label1, rot = 90, gp = formato_ejes),
                            bottom = textGrob(label2, gp = formato_ejes))
    
    plot_combinado <- arrangeGrob(plot_hm1, plot_hm2, plot_hm3, ncol = 3, 
                                  top = textGrob(titulo_auto, gp = gpar(fontsize = 16, fontface = "bold")))
  } else {
    plot_combinado <- arrangeGrob(plot_hm1, plot_hm2, ncol = 2, 
                                  top = textGrob(titulo_auto, gp = gpar(fontsize = 16, fontface = "bold")))
  }
  
  return(plot_combinado)
}

# Guarda la figura de mat_confusion() como
# <carpeta_salida>/<col1>_vs_<col2>_combinado.png (20 o 14 pulgadas de ancho).
#' @export
guardar_matriz <- function(obj, col1, col2, carpeta_salida, mostrar_conteo = TRUE) {
  
  # Pasamos la preferencia a la función principal
  mat_plot <- mat_confusion(obj, col1, col2, mostrar_conteo = mostrar_conteo)
  
  file_png <- paste0(col1, "_vs_", col2, "_combinado.png")
  filename <- file.path(carpeta_salida, file_png)
  
  # Ajuste dinámico del ancho: 20 si son tres gráficos, 14 si son dos
  ancho_plot <- ifelse(mostrar_conteo, 20, 14)
  
  ggsave(filename = filename, 
         plot = mat_plot, 
         bg = "white",
         width = ancho_plot,      
         height = 7,      
         dpi = 300)
}

