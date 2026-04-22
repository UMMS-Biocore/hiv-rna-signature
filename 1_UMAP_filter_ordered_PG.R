# ==== Wei et al. (Immunity 2023) — scRNA-only pipeline to UMAP (Harmony on RNA PCs) ====
# INPUT : Please put your input location
# OUTPUT: Please put your output location

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(readr); library(stringr)
  library(dplyr);  library(purrr);  library(ggplot2)
})
set.seed(1234)
options(stringsAsFactors = FALSE)

# ----- Paths -----
# PROJECT_ROOT must be set in the environment. No hardcoded paths.
project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  stop("PROJECT_ROOT environment variable is not set. ",
       "Set it to the absolute path of your project directory, e.g.:\n",
       "  PROJECT_ROOT=/path/to/project Rscript ", "this_script.R\n",
       "Or use scripts/run_pipeline.sh, which handles this for you.")
}
input_dir    <- file.path(project_root, "GSE239909_RAW")
out_dir      <- file.path(project_root, "Outputs_Code1_20251130")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- I/O helpers -----
read_mtx_geo <- function(matrix_path, features_path, barcodes_path) {
  m <- Matrix::readMM(matrix_path)
  feats <- suppressWarnings(readr::read_tsv(features_path, col_names = FALSE, show_col_types = FALSE))
  bcs   <- suppressWarnings(readr::read_tsv(barcodes_path,  col_names = FALSE, show_col_types = FALSE))[[1]]
  gene_col <- if (ncol(feats) >= 2) 2 else 1  # prefer symbols if present
  gnames <- as.character(feats[[gene_col]])
  bad <- is.na(gnames) | !nzchar(gnames)
  if (any(bad)) gnames[bad] <- as.character(feats[[1]])[bad]
  gnames <- make.unique(gnames)
  stopifnot(nrow(m) == length(gnames), ncol(m) == length(bcs))
  rownames(m) <- gnames; colnames(m) <- bcs
  m
}
ggsave_tiff <- function(path, plot, w=7, h=6, dpi=300) {
  # NOTE: using no compression per your preference
  ggsave(path, plot = plot, device = "tiff", width = w, height = h, units = "in", dpi = dpi, compression = "none")
}

# ----- File discovery (top-level GEO MTX triples) -----
mtx_files <- unique(c(
  Sys.glob(file.path(input_dir, "*_matrix.mtx.gz")),
  Sys.glob(file.path(input_dir, "*_matrix.mtx"))
))
stopifnot("No *_matrix.mtx[.gz] files found in input_dir" = length(mtx_files) > 0)

samples <- purrr::map_dfr(mtx_files, function(m) {
  base1 <- sub("_matrix\\.mtx(\\.gz)?$", "", basename(m))
  d     <- dirname(m)
  feat <- c(file.path(d, paste0(base1,"_features.tsv.gz")),
            file.path(d, paste0(base1,"_features.tsv")))
  bc   <- c(file.path(d, paste0(base1,"_barcodes.tsv.gz")),
            file.path(d, paste0(base1,"_barcodes.tsv")))
  feat <- feat[file.exists(feat)][1]; bc <- bc[file.exists(bc)][1]
  if (length(feat)==0 || length(bc)==0) return(tibble())
  tibble(
    mtx = m, features = feat, barcodes = bc,
    gsm = sub("^(GSM\\d+)_.*$", "\\1", basename(m)),
    sid = sub("^GSM\\d+_([A-Za-z0-9-]+)_matrix.*$", "\\1", basename(m))
  )
}) |> dplyr::distinct()
stopifnot("Could not pair features/barcodes for any matrix file" = nrow(samples) > 0)
write.csv(samples, file.path(out_dir, "found_mtx_triples.csv"), row.names = FALSE)

# ----- Build per-sample objects + paper-like RNA QC -----
objs <- vector("list", nrow(samples))
for (i in seq_len(nrow(samples))) {
  row <- samples[i,]
  message("Reading ", row$sid, " (", row$gsm, ")")
  mat <- read_mtx_geo(row$mtx, row$features, row$barcodes)
  colnames(mat) <- paste(row$sid, colnames(mat), sep = "_")  # globally unique
  
  so <- CreateSeuratObject(counts = mat, project = paste0("Wei_", row$sid),
                           min.cells = 3, min.features = 200)
  so$sample_id <- row$sid; so$GSM <- row$gsm
  
  # percent.mt (safe if MT- genes exist; else set 0)
  if (any(grepl("^MT-", rownames(so)))) {
    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")
  } else {
    so$percent.mt <- 0
  }
  
  # ---- QC thresholds as in Wei et al. (RNA-only) ----
  # keep cells with: mt% < 25; nFeature_RNA > 200; 500 < nCount_RNA < 10000
  so <- subset(so, subset = percent.mt < 25 & nFeature_RNA > 200 &
                 nCount_RNA > 500 & nCount_RNA < 10000)
  
  objs[[i]] <- so
}
qc_counts <- data.frame(sample = vapply(objs, \(x) unique(x$sample_id), character(1)),
                        cells_after_QC = sapply(objs, ncol))
write.csv(qc_counts, file.path(out_dir, "cells_per_sample_after_QC.csv"), row.names = FALSE)

# ----- Merge - RNA normalization - HVGs - Scale - PCA -----
combined <- Reduce(function(a,b) merge(a, b), objs)
DefaultAssay(combined) <- "RNA"

combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                           scale.factor = 1e4, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst",
                                 nfeatures = 3000, verbose = FALSE)

# do NOT regress nCount_RNA or percent.mt here
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)

# ----- Batch correction + graph/UMAP/clusters + save/plots -----
use_harmony <- TRUE  # paper uses Harmony for RNA batch correction
DefaultAssay(combined) <- "RNA"
if (!"pca" %in% Reductions(combined)) {
  combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
}

reduction_used <- "pca"
if (use_harmony) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("The 'harmony' package is required when use_harmony=TRUE. Install with install.packages('harmony').")
  }
  combined <- harmony::RunHarmony(combined, group.by.vars = "sample_id")
  reduction_used <- "harmony"
}

graph_name <- paste0(reduction_used, "_snn")
umap_name  <- paste0(reduction_used, "_umap")

combined <- FindNeighbors(combined, reduction = reduction_used, dims = 1:30,
                          graph.name = graph_name, verbose = FALSE)
combined <- RunUMAP(combined, reduction = reduction_used, dims = 1:30,
                    reduction.name = umap_name, reduction.key = "UMAP_", verbose = FALSE)

# Seurat SLM clustering at resolution = 0.6 (paper picked via Clustree)
combined <- FindClusters(combined, graph.name = graph_name, resolution = 0.6, verbose = FALSE)

# ----- Plot & save -----
p1 <- DimPlot(combined, reduction = umap_name, group.by = "seurat_clusters",
              label = TRUE, repel = TRUE) +
  ggtitle(sprintf("RNA UMAP • Clusters (%s)", if (use_harmony) "Harmony" else "PCA"))
ggsave_tiff(file.path(out_dir, sprintf("UMAP_clusters_%s.tiff", reduction_used)), p1)

# Dump metadata & object
write.csv(combined@meta.data,
          file.path(out_dir, sprintf("combined_metadata_%s.csv", reduction_used)),
          row.names = FALSE)

saveRDS(
  combined,
  file.path(out_dir, sprintf("Wei_scRNA_RNA_%s_UMAP.rds",
                             if (use_harmony) "Harmony" else "PCA"))
)

# ----- Session info for reproducibility -----
sink(file.path(out_dir, sprintf("sessionInfo_%s.txt", reduction_used)))
print(sessionInfo())
sink()

message("✅ Done. Reduction=", reduction_used, " • Outputs saved to: ", out_dir)
