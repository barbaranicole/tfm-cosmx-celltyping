# ==========================================
# Agregar resultados de clustering a cada joinmerged.rds
# Requiere haber ejecutado antes Resultados_clustering.R (TXT unificado).
# Lee <version>_joinmerged.rds (no se modifica) y guarda el resultado en
# <version>_joinmerged_con_clustering.rds.
# ==========================================

library(Seurat)

# --- CONFIGURACIÓN ---
base_dir    <- path.expand("~/Documents/Master/TFM/Resultados_Algoritmos_v1/")
merged_txt  <- path.expand("~/Documents/Master/TFM/Resultados_Algoritmos_v1/clustering_resultados.txt")
versiones   <- c("original", "thin_75")
# Nombre del RDS de salida:
sufijo_salida <- "joinmerged_clustering"

if (!file.exists(merged_txt)) stop("Ejecuta antes Resultados_clustering.R para generar ", merged_txt)

# Leer TXT unificado una sola vez
annot_full <- read.delim(merged_txt, stringsAsFactors = FALSE)
if (any(duplicated(annot_full$cell_ID))) annot_full <- annot_full[!duplicated(annot_full$cell_ID), ]

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
