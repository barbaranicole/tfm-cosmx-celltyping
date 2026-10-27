# ==============================================================================
# 02_Clustering_joinmerged.R
# ------------------------------------------------------------------------------
# Descripción: Añade los clusters de Leiden e InSituType a la metadata del
#              objeto Seurat de cada escenario (original y thin_75). A partir
#              de aquí, ambos métodos siguen el mismo procesamiento y análisis.
#
# Entradas:    Resultados_Algoritmos_v1/clustering_resultados.txt
#              Resultados_Algoritmos_v1/<escenario>_joinmerged.rds (no se modifica)
# Salidas:     Resultados_Algoritmos_v1/<escenario>_joinmerged_clustering.rds
# Requiere:    01_Resultados_clustering.R
# Siguiente:   03_Marcadores.R
# ==============================================================================

library(Seurat)

# ==========================================
# 0. CONFIGURACIÓN
# ==========================================

base_dir    <- here::here("Resultados_Algoritmos_v1/")
merged_txt  <- here::here("Resultados_Algoritmos_v1/clustering_resultados.txt")
versiones   <- c("original", "thin_75")
# Nombre del RDS de salida:
sufijo_salida <- "joinmerged_clustering"

if (!file.exists(merged_txt)) stop("Ejecuta antes 01_Resultados_clustering.R para generar ", merged_txt)
# ==========================================
# 1. LECTURA DE LA TABLA DE CLUSTERING
# ==========================================

annot_full <- read.delim(merged_txt, stringsAsFactors = FALSE)
if (any(duplicated(annot_full$cell_ID))) annot_full <- annot_full[!duplicated(annot_full$cell_ID), ]

# ==========================================
# 2. AÑADIR CLUSTERS A CADA OBJETO SEURAT
# ==========================================

for (version in versiones) {
  rds_path <- file.path(base_dir, paste0(version, "_joinmerged.rds"))
  if (!file.exists(rds_path)) {
    cat("[aviso] No existe:", rds_path, "\n");next
  }

  # Columnas de esta versión (Leiden_<v>_* e IST_<v>_*)
  pref_leiden <- paste0("Leiden_", version, "_")
  pref_ist    <- paste0("IST_", version, "_")
  nms         <- names(annot_full)
  cols        <- c("cell_ID",
                  nms[startsWith(nms, pref_leiden)],
                  nms[startsWith(nms, pref_ist)])

  if (length(cols) < 2) {
    cat(basename(rds_path), ": sin columnas para", version, "\n"); next 
  }

  obj   <- readRDS(rds_path)
  cells <- colnames(obj)
  if (length(cells) == 0) cells <- Cells(obj)

  cols_meta <- setdiff(cols, "cell_ID")
  nuevos    <- setdiff(cols_meta, colnames(obj@meta.data))
  if (length(nuevos) == 0) { cat(basename(rds_path), ": ya actualizado.\n"); next }

  idx     <- match(cells, annot_full$cell_ID)
  meta_df <- annot_full[idx, nuevos, drop = FALSE]
  rownames(meta_df) <- cells

  obj <- AddMetaData(obj, metadata = meta_df)
  out_path <- file.path(base_dir, paste0(version, "_", sufijo_salida, ".rds"))
  saveRDS(obj, out_path)
  cat("Guardado:", basename(out_path), "(", length(nuevos), "columnas:", paste(nuevos, collapse = ", "), ")\n")
}

cat("Listo.\n")
