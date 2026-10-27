# ==============================================================================
# 01_Resultados_clustering.R
# ------------------------------------------------------------------------------
# Descripción: Une en una sola tabla los resultados de clustering de Leiden e
#              InSituType de ambos escenarios (original y thin_75). Cada columna
#              corresponde a un método, escenario y parámetro.
#
# Entradas:    Resultados_Algoritmos_v1/Leiden/<escenario>/*_res<r>_<fecha>.txt
#              Resultados_Algoritmos_v1/Insitutype/<escenario>/*_nclusts<n>_<fecha>.txt
#              (lista fija en `lista_archivos`; si se vuelven a ejecutar los
#              métodos, hay que actualizar los nombres con la nueva fecha)
# Salidas:     Resultados_Algoritmos_v1/clustering_resultados.txt
#              Columnas: cell_ID, Leiden_<escenario>_res<r>,
#              IST_<escenario>_nclusts<n> e IST_<escenario>_nclusts<n>_conf
# Requiere:    Resultados de Metodos/02_Seurat.R y Metodos/03_IST.R
# Siguiente:   02_Clustering_joinmerged.R
# ==============================================================================

library(dplyr)

# ==========================================
# 0. CONFIGURACIÓN
# ==========================================
base_dir   <- here::here("Resultados_Algoritmos_v1/")
out_txt    <- here::here("Resultados_Algoritmos_v1/clustering_resultados.txt")
versiones  <- c("original", "thin_75")

# Lista de archivos: path, metodo (Leiden / IST), version (original / thin_75)
# Leiden: se usan las resoluciones 0.8, 1 y 1.2 (la 0.6 se descartó).
# Solo archivos de resultados (no Log_*)
lista_archivos <- list(
  # Leiden - original
  list(path = file.path(base_dir, "Leiden/original/original_merged_dims30_res0.8_20260224_121152.txt"),  metodo = "Leiden", version = "original"),
  list(path = file.path(base_dir, "Leiden/original/original_merged_dims30_res1_20260224_121713.txt"),   metodo = "Leiden", version = "original"),
  list(path = file.path(base_dir, "Leiden/original/original_merged_dims30_res1.2_20260406_142020.txt"),   metodo = "Leiden", version = "original"),
  # Leiden - thin_75
  list(path = file.path(base_dir, "Leiden/thin_75/thin_75_merged_dims30_res0.8_20260303_121326.txt"),  metodo = "Leiden", version = "thin_75"),
  list(path = file.path(base_dir, "Leiden/thin_75/thin_75_merged_dims30_res1_20260303_121912.txt"),   metodo = "Leiden", version = "thin_75"),
  list(path = file.path(base_dir, "Leiden/thin_75/thin_75_merged_dims30_res1.2_20260406_141316.txt"),   metodo = "Leiden", version = "thin_75"),
  
  # InSituType - original
  list(path = file.path(base_dir, "Insitutype/original/original_joinmerged_nclusts6_20260309_131933.txt"),  metodo = "IST", version = "original"),
  list(path = file.path(base_dir, "Insitutype/original/original_joinmerged_nclusts8_20260309_134002.txt"),  metodo = "IST", version = "original"),
  list(path = file.path(base_dir, "Insitutype/original/original_joinmerged_nclusts10_20260309_140111.txt"), metodo = "IST", version = "original"),
  # InSituType - thin_75
  list(path = file.path(base_dir, "Insitutype/thin_75/thin_75_joinmerged_nclusts6_20260309_142304.txt"),  metodo = "IST", version = "thin_75"),
  list(path = file.path(base_dir, "Insitutype/thin_75/thin_75_joinmerged_nclusts8_20260309_143857.txt"),  metodo = "IST", version = "thin_75"),
  list(path = file.path(base_dir, "Insitutype/thin_75/thin_75_joinmerged_nclusts10_20260309_145723.txt"),   metodo = "IST", version = "thin_75")
)

# ==========================================
# 1. LEER CADA TXT Y CONSTRUIR TABLA CON COLUMNAS NOMBRADAS
# ==========================================
tablas <- list()

for (i in seq_along(lista_archivos)) {
  entrada <- lista_archivos[[i]]
  path    <- entrada$path
  metodo  <- entrada$metodo
  version <- entrada$version

  if (!file.exists(path)) {
    cat("[aviso] No existe:", path, "\n")
    next
  }

  d <- read.delim(path, stringsAsFactors = FALSE)
  if (any(duplicated(d$cell_ID))) d <- d[!duplicated(d$cell_ID), ]

  if (metodo == "Leiden") {
    sufijo <- regmatches(path, regexpr("res[0-9\\.]+", path))
    nom_col <- paste0("Leiden_", version, "_", sufijo)
    tablas[[i]] <- d %>% select(cell_ID) %>% mutate(!!nom_col := d$label_seurat)
  } else {
    sufijo <- regmatches(path, regexpr("nclusts[0-9]+", path))
    nom_class <- paste0("IST_", version, "_", sufijo)
    nom_conf  <- paste0("IST_", version, "_", sufijo, "_conf")
    tablas[[i]] <- d %>%
      select(cell_ID) %>%
      mutate(!!nom_class := d$ist_class, !!nom_conf := as.numeric(d$ist_confianza))
  }
  cat("  Leído:", basename(path), "->", paste(names(tablas[[i]])[-1], collapse = ", "), "\n")
}

tablas <- tablas[!sapply(tablas, is.null)]
if (length(tablas) == 0) stop("No se leyó ningún archivo.")

# ==========================================
# 2. UNIR TODAS LAS TABLAS POR cell_ID (full join)
# ==========================================
merged_df <- tablas[[1]]
for (j in seq_along(tablas)[-1]) {
  merged_df <- full_join(merged_df, tablas[[j]], by = "cell_ID")
}
merged_df <- merged_df %>% arrange(cell_ID)

# ==========================================
# 3. REORDENAR COLUMNAS: cell_ID, Leiden por versión, IST por versión (original, thin_75)
# ==========================================
nms <- names(merged_df)
id_col <- "cell_ID"
orden <- id_col
for (v in versiones) {
  # Leiden: columnas Leiden_<version>_res*
  leiden_v <- sort(nms[grepl(paste0("^Leiden_", v, "_"), nms)])
  orden <- c(orden, leiden_v)
}
for (v in versiones) {
  # IST: class y conf por nclusts (orden numérico)
  pref <- paste0("^IST_", v, "_")
  ist_v      <- nms[grepl(pref, nms) & !grepl("_conf$", nms)]
  ist_v      <- ist_v[order(as.numeric(gsub(paste0("IST_", v, "_nclusts"), "", ist_v)))]
  ist_v_conf <- nms[grepl(paste0(pref, ".*_conf$"), nms)]
  ist_v_conf <- ist_v_conf[order(as.numeric(gsub(paste0("IST_", v, "_nclusts|_conf"), "", ist_v_conf)))]
  orden <- c(orden, c(rbind(ist_v, ist_v_conf)))
}
orden <- orden[orden %in% nms]
merged_df <- merged_df[, orden, drop = FALSE]

# ==========================================
# 4. GUARDAR TXT
# ==========================================
dir_out <- dirname(out_txt)
if (nchar(dir_out) > 0 && !dir.exists(dir_out)) dir.create(dir_out, recursive = TRUE)
write.table(merged_df, file = out_txt, sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nResultado guardado:", out_txt, "\n")
cat("  Filas:", nrow(merged_df), "| Columnas:", ncol(merged_df), "\n")
