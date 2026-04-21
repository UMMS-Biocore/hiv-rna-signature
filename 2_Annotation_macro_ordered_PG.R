#####################################################################################
# ==== Code 2: Macro programs (RNA-only) via direct per-cell scoring — forced argmax ====
# INPUT : RDS from Code 1 (UMAP/Harmony step)
# OUTPUT: Annotated Seurat object + macro tables/plot

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(tibble); library(ggplot2); library(Matrix)
})
set.seed(1234)
options(stringsAsFactors = FALSE)

# ----- Paths (EDIT THESE) -----
project_root <- "{Your_Folder}"
in_dir  <- file.path(project_root, "Outputs_Code1_20251130")  # where Code 1 wrote the RDS
out_dir <- file.path(project_root, "Outputs_Code2_20251130")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- Load input from Code 1 -----
# If you used PCA instead of Harmony in Code 1, change the file name accordingly.
combined <- readRDS(file.path(in_dir, "Wei_scRNA_RNA_Harmony_UMAP.rds"))

stopifnot(exists("combined"), inherits(combined, "Seurat"))
DefaultAssay(combined) <- "RNA"
Idents(combined) <- "seurat_clusters"

# Join layers for Seurat v5 (single "data" matrix)
if (inherits(combined[["RNA"]], "Assay5")) {
  if (!"data" %in% Layers(combined[["RNA"]])) {
    combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)
  }
  combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
}

# --- Robust marker sets (4 macros) ---
macro_panels <- list(
  Proliferation = c("MKI67","TOP2A","BIRC5","CENPF","TPX2","ASPM","RRM2","CCNB1","CCNB2","CDK1",
                    "HJURP","BUB1","KIF11","KIF15","KIF23","KIF2C","PCLAF","TYMS","UHRF1","PRR11"),
  Cytotoxic     = c("GZMB","PRF1","NKG7","CTSW","GNLY","KLRD1","HCST","ZEB2","FGFBP2"),
  AP1_TNF       = c("FOS","FOSB","JUN","JUNB","JUND","DUSP1","DUSP2","TNFAIP3","NFKBIA","ZFP36",
                    "ZFP36L1","ZFP36L2","IER2","IER5","GADD45B","KLF2","KLF6","RELB","NFKBIZ","EGR1"),
  IRF_IFN       = c("ISG15","IFIT1","IFIT3","IFITM1","IFITM2","IFITM3","MX1","MX2","OAS1","OAS2",
                    "OAS3","OASL","RSAD2","USP18","IRF1","IRF7","IRF9","STAT1","STAT2","BST2","XAF1","HERC6","CXCL10")
)

# Keep only genes that exist
present <- rownames(combined)
macro_panels <- lapply(macro_panels, function(g) intersect(unique(g), present))
stopifnot("No macro marker genes found in this object." =
            sum(vapply(macro_panels, length, integer(1))) > 0)

# --- Pull normalized log-expression matrix (brace style to avoid else issues) ---
if (inherits(combined[["RNA"]], "Assay5")) {
  M <- LayerData(combined[["RNA"]], layer = "data")
} else {
  M <- GetAssayData(combined, assay = "RNA", slot = "data")
}

# --- Per-cell average scores for each program ---
score_df <- lapply(names(macro_panels), function(nm){
  genes <- macro_panels[[nm]]
  if (length(genes) == 0) return(rep(NA_real_, ncol(M)))
  Matrix::colMeans(as.matrix(M[genes, , drop = FALSE]))
}) |> setNames(names(macro_panels)) |> as.data.frame()

rownames(score_df) <- colnames(M)
colnames(score_df) <- paste0("MAC_", names(score_df))  # MAC_Proliferation, ...

# Attach scores to metadata
combined <- AddMetaData(combined, metadata = score_df)

# --- Per-cluster means and Z across clusters (per program) ---
md <- combined@meta.data %>% mutate(cluster = as.character(Idents(combined)))
mean_df <- md %>%
  group_by(cluster) %>%
  summarise(across(starts_with("MAC_"), \(x) mean(x, na.rm = TRUE)), .groups = "drop")

z_df_macro <- mean_df
for (cn in grep("^MAC_", names(mean_df), value = TRUE)) {
  mu <- mean(mean_df[[cn]], na.rm = TRUE); sdv <- sd(mean_df[[cn]], na.rm = TRUE)
  z_df_macro[[paste0("Z_", cn)]] <- if (!is.na(sdv) && sdv > 0) (mean_df[[cn]] - mu)/sdv else 0
}

# --- Build program Z table (no NSE on 'cluster'); forced argmax (no 'Other') ---
getZ <- function(tag) {
  z <- paste0("Z_MAC_", tag)
  if (z %in% names(z_df_macro)) z_df_macro[[z]] else rep(-Inf, nrow(z_df_macro))
}

prog_macro <- tibble::tibble(
  cluster       = as.character(z_df_macro$cluster),
  Proliferation = getZ("Proliferation"),
  Cytotoxic     = getZ("Cytotoxic"),
  `AP-1/TNF`    = getZ("AP1_TNF"),
  `IRF/IFN`     = getZ("IRF_IFN")
)

assign_macro <- function(irow){
  v <- as.numeric(irow[-1])
  v[!is.finite(v)] <- -Inf
  names(v) <- colnames(prog_macro)[-1]
  names(v)[which.max(v)]
}

macro_tbl <- data.frame(
  cluster = prog_macro$cluster,
  Macro   = apply(prog_macro, 1, assign_macro),
  row.names = NULL
)

# --- Attach macro labels to cells (exactly 4 labels) ---
lab_map <- setNames(macro_tbl$Macro, macro_tbl$cluster)
cell_macro <- unname(lab_map[as.character(Idents(combined))])
combined <- AddMetaData(
  combined,
  metadata = setNames(data.frame(annotation_macro = cell_macro, row.names = colnames(combined)),
                      "annotation_macro")
)
combined$annotation_macro <- factor(
  combined$annotation_macro,
  levels = c("Proliferation","Cytotoxic","AP-1/TNF","IRF/IFN")
)

# --- Pick UMAP for plotting (brace style) ---
if ("harmony_umap" %in% names(combined@reductions)) {
  umap_name <- "harmony_umap"
} else if ("pca_umap" %in% names(combined@reductions)) {
  umap_name <- "pca_umap"
} else {
  umap_name <- Reductions(combined)[[1]]
}

# --- Plot & save (no LZW compression) ---
pal_macro <- c("Proliferation"="#F8766D","Cytotoxic"="#7CAE00","AP-1/TNF"="#00BFC4","IRF/IFN"="#C77CFF")
p_macro <- DimPlot(combined, reduction = umap_name, group.by = "annotation_macro",
                   label = FALSE, repel = TRUE, cols = pal_macro) +
  ggtitle("Macro programs (RNA-only): Proliferation / Cytotoxic / AP-1/TNF / IRF/IFN") +
  theme(
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 20))
  

write.csv(macro_tbl, file.path(out_dir, "macro_annotation_by_cluster.csv"), row.names = FALSE)
write.csv(prog_macro, file.path(out_dir, "macro_program_cluster_Z.csv"), row.names = FALSE)
ggsave(file.path(out_dir, "UMAP_macro_programs.tiff"), p_macro,
       width = 6, height = 5, dpi = 300, compression = "none")

# ----- Save output for next step -----
saveRDS(
  combined,
  file.path(out_dir, "Wei_scRNA_RNA_Annotated_UMAP.rds")
)

# Console visibility
print(macro_tbl)
print(table(combined$annotation_macro))
#####################################################################################
