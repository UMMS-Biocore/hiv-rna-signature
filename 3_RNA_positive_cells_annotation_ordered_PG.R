#####################################################################################
# ==== Code 3: Append HIV_RNA annotations (continue from annotated UMAP object) ====
# INPUT : RDS from Code 2 → "Wei_scRNA_RNA_Annotated_UMAP.rds"
# OUTPUT: Updated object + HIV_RNA summary + UMAP highlighting RNA+ cells

suppressPackageStartupMessages({
  library(Seurat); library(readr); library(dplyr); library(stringr); library(ggplot2)
})
set.seed(1234)
options(stringsAsFactors = FALSE)

# ----- Paths -----
project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(project_root)) {
  stop("PROJECT_ROOT environment variable is not set. ",
       "Set it to the absolute path of your project directory, or use scripts/run_pipeline.sh.")
}
input_dir    <- file.path(project_root, "GSE239909_RAW")
in_dir       <- file.path(project_root, "Outputs_Code2_20251130")
out_dir      <- file.path(project_root, "Outputs_Code3_20251130")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- Load input from Code 2 -----
combined <- readRDS(file.path(in_dir, "Wei_scRNA_RNA_Annotated_UMAP.rds"))

DefaultAssay(combined) <- "RNA"

# ---- Read HIV RNA calls (TXT with columns: cell | hiv_rna | copies) ----
txt_path <- if (file.exists(file.path(input_dir, "bc_RNA_HIV.txt"))) {
  file.path(input_dir, "bc_RNA_HIV.txt")
} else {
  file.path(project_root, "bc_RNA_HIV.txt")
}
stopifnot("bc_RNA_HIV.txt not found in input_dir or project_root" = file.exists(txt_path))

hiv_tbl <- suppressMessages(readr::read_tsv(txt_path, col_types = cols(.default = "c")))
names(hiv_tbl) <- tolower(names(hiv_tbl))
stopifnot("TXT must contain columns: cell, hiv_rna" = all(c("cell","hiv_rna") %in% names(hiv_tbl)))
if (!"copies" %in% names(hiv_tbl)) hiv_tbl$copies <- NA_character_

# Fix common OCR typo: YW1O_ -> YW10_
hiv_tbl$cell <- gsub("^YW1O_", "YW10_", hiv_tbl$cell)

# Normalize & enforce ≥2 reads rule for RNA+
hiv_tbl <- hiv_tbl %>%
  mutate(
    cell    = trimws(cell),
    hiv_rna = tolower(trimws(hiv_rna)),
    copies  = suppressWarnings(as.integer(trimws(copies))),
    rna_pos = ifelse(!is.na(copies), copies >= 2,
                     hiv_rna %in% c("positive","pos","rna+","yes","true"))
  ) %>%
  filter(rna_pos) %>%
  group_by(cell) %>% slice_max(order_by = copies, n = 1, with_ties = FALSE) %>% ungroup()

# ---- Exact-match ONLY to avoid cross-sample collisions ----
obj_cols  <- colnames(combined)
n_txt     <- nrow(hiv_tbl)
exact_idx <- match(hiv_tbl$cell, obj_cols, nomatch = 0L)
n_found   <- sum(exact_idx != 0L)
n_missing <- n_txt - n_found

# Prepare/append meta columns
if (!"HIV_RNA"         %in% colnames(combined@meta.data)) combined$HIV_RNA         <- NA_character_
if (!"HIV_Copies"      %in% colnames(combined@meta.data)) combined$HIV_Copies      <- NA_integer_
if (!"HIV_RNA_source"  %in% colnames(combined@meta.data)) combined$HIV_RNA_source  <- NA_character_

if (n_found > 0) {
  matched_obj <- obj_cols[exact_idx[exact_idx != 0L]]
  combined@meta.data[matched_obj, "HIV_RNA"]        <- "positive"
  combined@meta.data[matched_obj, "HIV_Copies"]     <- hiv_tbl$copies[exact_idx != 0L]
  combined@meta.data[matched_obj, "HIV_RNA_source"] <- "bc_RNA_HIV.txt (exact, copies>=2)"
}

# ---- Report counts + per-sample breakdown (if sample_id exists) ----
summary_lines <- c(
  paste("TXT path:", txt_path),
  paste("Total RNA+ rows in TXT after Copies>=2 & dedup:", n_txt),
  paste("RNA+ cells found in Seurat (exact match):", n_found),
  paste("Unmatched in Seurat:", n_missing)
)
writeLines(summary_lines, file.path(out_dir, "HIV_RNA_import_summary.txt"))
if ("sample_id" %in% colnames(combined@meta.data) && n_found > 0) {
  by_sample <- table(combined@meta.data$sample_id[obj_cols %in% hiv_tbl$cell])
  capture.output(by_sample, file = file.path(out_dir, "HIV_RNA_import_summary.txt"), append = TRUE)
}
cat(paste(summary_lines, collapse = "\n"), "\n")

# ---- Plot HIV_RNA on Harmony UMAP (keep original plotting logic) ----
red <- "harmony_umap"
emb <- Embeddings(combined, reduction = red)[, 1:2]
colnames(emb) <- c("x", "y")

df <- as.data.frame(emb)
df$HIV_RNA <- combined$HIV_RNA

# counts summary
n_total <- nrow(df)
n_pos   <- sum(tolower(trimws(as.character(df$HIV_RNA))) == "positive", na.rm = TRUE)
cat(sprintf("HIV_RNA == 'positive': %s of %s (%.1f%%)\n",
            format(n_pos, big.mark=","), format(n_total, big.mark=","),
            100 * n_pos / n_total))

# plot: grey background + big red positives
p <- ggplot(df, aes(x, y)) +
  geom_point(data = df, color = "grey80", size = 0.3, alpha = 0.35) +
  geom_point(data = df %>% filter(tolower(trimws(as.character(HIV_RNA))) == "positive"),
             color = "red", size = 0.5, alpha = 0.6) +
  labs(title = "HIV_RNA", x = "UMAP_1", y = "UMAP_2") +
  theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 20),
    legend.position = "none"
  )

# save high-res TIFF (no LZW compression)
ggsave(
  filename = file.path(out_dir, "HIV_RNA_harmonyUMAP.tiff"),
  plot = p, device = "tiff",
  width = 5, height = 5, units = "in",
  dpi = 300, compression = "none"
)

# also print overall count on the full object
n_total <- ncol(combined)
n_pos   <- sum(tolower(trimws(as.character(combined$HIV_RNA))) == "positive", na.rm = TRUE)
cat(sprintf("HIV_RNA == 'positive': %s of %s (%.1f%%)\n",
            format(n_pos, big.mark = ","), format(n_total, big.mark = ","),
            100 * n_pos / n_total))

# ---- Save updated object for next step ----
saveRDS(combined, file.path(out_dir, "Wei_scRNA_with_HIV_RNA_exact_annotated.rds"))

# sanity check
print(table(HIV_RNA = combined$HIV_RNA, useNA = "ifany"))
#####################################################################################
