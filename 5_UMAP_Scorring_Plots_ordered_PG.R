#####################################################################################
# ======================= Code 5: Prevalence + Triptych UMAPs =======================
# INPUT : Code 4 output ONLY → "Wei_scRNA_with_profiles_scored.rds"
# OUTPUT: Prevalence CSVs + UMAP triptychs (z ∈ {1.0, 1.5, 2.0})

suppressPackageStartupMessages({
  library(Seurat)    # for Embeddings / DimPlot
  library(dplyr); library(tidyr); library(ggplot2)
  library(scales); library(viridis); library(patchwork)
})

set.seed(1234)
options(stringsAsFactors = FALSE)

# ----- Paths (EDIT THESE) ----------------------------------------------------------
project_root <- "{Your_Folder}"
in_dir_code4   <- file.path(project_root, "Outputs_Code4_20251130")
out_dir        <- file.path(project_root, "Outputs_Code5_20251130")
if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }

# ----- Load input (only from Code 4) ----------------------------------------------
if (!exists("combined")) {
  rds_path <- file.path(in_dir_code4, "Wei_scRNA_with_profiles_scored.rds")
  stopifnot("Expected input from Code 4 not found." = file.exists(rds_path))
  combined <- readRDS(rds_path)
}
stopifnot(inherits(combined, "Seurat"))

# ----- Basic setup ----------------------------------------------------------------
macro_levels <- c("Proliferation","Cytotoxic","AP-1/TNF","IRF/IFN")
combined$annotation_macro <- factor(as.character(combined$annotation_macro), levels = macro_levels)

# Legend labels
dataset_labels <- c("our" = "Bulk RNA-seq profile",
                    "article" = "Wei et al. RNA+ profile")

# Colors for prevalence bars
cols_manual <- c("Wei et al. RNA+ profile" = "#33638D",
                 "Bulk RNA-seq profile"    = "#35B779")

# ----- Helper: prevalence plot ----------------------------------------------------
plot_prevalence_auto <- function(prev_tbl, thr) {
  ytop <- min(1, max(prev_tbl$pct, na.rm = TRUE) * 1.25)
  ggplot(prev_tbl, aes(x = annotation_macro, y = pct, fill = dataset)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.6) +
    geom_text(
      aes(label = sprintf("%s\n(n=%s)", percent(pct, 0.1), comma(pos))),
      position = position_dodge(width = 0.75),
      vjust = -0.25, size = 3.4
    ) +
    scale_y_continuous(labels = percent_format(accuracy = 1),
                       limits = c(0, ytop),
                       expand = expansion(mult = c(0, 0.05))) +
    scale_fill_manual(values = cols_manual, name = NULL) +
    labs(
      x = NULL, y = "Profile+ within macro (%)",
      title = sprintf("Profile+ prevalence by macro (z ≥ %.1f)", thr)
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 25, hjust = 1),
      panel.grid.minor = element_blank()
    )
}

# -------- main loop: ONLY prevalence & tables (composition removed) ---------------
# sanity checks for required columns (will be checked per-threshold too)
req_cols_soft <- c("annotation_macro","our_profile_z","article_profile_z")
stopifnot(all(req_cols_soft %in% colnames(combined@meta.data)))

for (thr in c(1.0, 1.5, 2.0)) {
  our_col <- paste0("our_data", sprintf("%.1f", thr))
  art_col <- paste0("article_data", sprintf("%.1f", thr))
  stopifnot(all(c(our_col, art_col) %in% names(combined@meta.data)))
  
  md <- combined@meta.data %>%
    transmute(
      annotation_macro = droplevels(annotation_macro),
      our_call = .data[[our_col]],
      art_call = .data[[art_col]]
    ) %>%
    filter(!is.na(annotation_macro))
  
  long <- md %>%
    mutate(our = our_call == "Profile+",
           article = art_call == "Profile+") %>%
    select(annotation_macro, our, article) %>%
    pivot_longer(c(our, article), names_to = "dataset", values_to = "is_pos") %>%
    mutate(dataset = recode(dataset, !!!dataset_labels))
  
  # Prevalence within macro
  prev_tbl <- long %>%
    group_by(annotation_macro, dataset) %>%
    summarise(pos = sum(is_pos, na.rm = TRUE),
              total = n(),
              pct = ifelse(total > 0, pos / total, 0),
              .groups = "drop")
  
  p_prev <- plot_prevalence_auto(prev_tbl, thr)
  ggsave(
    file.path(out_dir, sprintf("PREVALENCE_macro_Bulk_vs_Wei_z%s.tiff",
                               gsub("\\.", "_", sprintf("%.1f", thr)))),
    p_prev, width = 9.5, height = 6.2, dpi = 300, compression = "none"
  )
  
  # Save tidy table
  write.csv(prev_tbl, file.path(out_dir, sprintf("TABLE_prevalence_Bulk_vs_Wei_z%s.csv",
                                                 gsub("\\.", "_", sprintf("%.1f", thr)))),
            row.names = FALSE)
}

################################    UMAPs   ###########################################################
## ========================= Triptych UMAPs (macro + two profiles) ===================================
stopifnot(all(c("our_profile_z","article_profile_z","annotation_macro") %in% colnames(combined@meta.data)))

umap_name <- if ("harmony_umap" %in% names(combined@reductions)) {
  "harmony_umap"
} else if ("pca_umap" %in% names(combined@reductions)) {
  "pca_umap"
} else {
  Reductions(combined)[[1]]
}

pal_macro <- c(
  "Proliferation" = "#F8766D",
  "Cytotoxic"     = "#7CAE00",
  "AP-1/TNF"      = "#00BFC4",
  "IRF/IFN"       = "#C77CFF"
)

pt_macro <- 0.5
pt_bg    <- 0.06
alpha_bg <- 0.06
pt_hi    <- 0.14
alpha_hi <- 0.65

col_bg       <- "#C0C0C0"
col_our      <- "#2CA25F"
col_article  <- "#2B8CBE"

emb <- as.data.frame(Embeddings(combined, reduction = umap_name))
colnames(emb)[1:2] <- c("UMAP_1","UMAP_2")
df_all <- cbind(emb, combined@meta.data[, c("annotation_macro","our_profile_z","article_profile_z")])

centroids <- df_all %>%
  dplyr::filter(!is.na(annotation_macro)) %>%
  dplyr::group_by(annotation_macro) %>%
  dplyr::summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

make_triptych_soft <- function(thr, out_dir = getwd()) {
  thr_tag <- format(thr, nsmall = 1)
  df <- df_all %>%
    dplyr::mutate(
      our_call = ifelse(our_profile_z >= thr, "Profile+", "Other"),
      art_call = ifelse(article_profile_z >= thr, "Profile+", "Other")
    )
  n_our <- sum(df$our_call == "Profile+", na.rm = TRUE)
  n_art <- sum(df$art_call == "Profile+", na.rm = TRUE)
  
  old_id <- Idents(combined)
  on.exit({ try({Idents(combined) <<- old_id}, silent = TRUE) }, add = TRUE)
  Idents(combined) <- "annotation_macro"
  
  p_macro <- DimPlot(
    combined, reduction = umap_name,
    cols = pal_macro, pt.size = pt_macro, label = FALSE
  ) +
    geom_text(
      data = centroids,
      aes(UMAP_1, UMAP_2, label = annotation_macro),
      size = 3, color = "black"
    ) +
    ggtitle("Annotated clusters UMAP") +
    theme(
      plot.title = element_text(size = 16),
      legend.title = element_blank(),
      legend.position = "right"
    )
  
  our_pos_df <- df %>% dplyr::filter(our_call == "Profile+") %>% as.data.frame()
  df_df      <- as.data.frame(df)
  p_our <- ggplot() +
    geom_point(data = df_df,      aes(UMAP_1, UMAP_2), color = col_bg,
               size = pt_bg, alpha = alpha_bg, shape = 16) +
    geom_point(data = our_pos_df, aes(UMAP_1, UMAP_2), color = col_our,
               size = pt_hi, alpha = alpha_hi, shape = 16) +
    labs(x = "UMAP_1", y = "UMAP_2") +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      panel.grid = element_blank(),
      axis.title = element_text(size = 20),
      axis.text  = element_text(size = 20),
      axis.line  = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    ggtitle(sprintf("Bulk RNA-seq profile (z ≥ %s)\nProfile+ n = %s",
                    thr_tag, scales::comma(n_our)))
  
  art_pos_df <- df %>% dplyr::filter(art_call == "Profile+") %>% as.data.frame()
  p_art <- ggplot() +
    geom_point(data = df_df,      aes(UMAP_1, UMAP_2), color = col_bg,
               size = pt_bg, alpha = alpha_bg, shape = 16) +
    geom_point(data = art_pos_df, aes(UMAP_1, UMAP_2), color = col_article,
               size = pt_hi, alpha = alpha_hi, shape = 16) +
    labs(x = "UMAP_1", y = "UMAP_2") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      panel.grid = element_blank(),
      axis.title = element_text(size = 20),
      axis.text  = element_text(size = 20),
      axis.line  = element_line(colour = "black", linewidth = 0.5),
      axis.ticks = element_line(colour = "black", linewidth = 0.5)
    ) +
    ggtitle(sprintf("Wei et al. RNA+ profile (z ≥ %s)\nProfile+ n = %s",
                    thr_tag, scales::comma(n_art)))
  
  g <- (p_macro | p_our | p_art) + plot_layout(widths = c(1, 1, 1))
  if (!dir.exists(out_dir)) { dir.create(out_dir, recursive = TRUE) }
  fn <- file.path(out_dir, sprintf("UMAP_triptych_soft_z%s.tiff", gsub("\\.", "", thr_tag)))
  ggsave(fn, g, width = 16, height = 5.5, dpi = 300, compression = "none")
  message("Saved: ", fn)
  invisible(g)
}

for (thr in c(1.0, 1.5, 2.0)) {
  make_triptych_soft(thr, out_dir)
}
#####################################################################################
