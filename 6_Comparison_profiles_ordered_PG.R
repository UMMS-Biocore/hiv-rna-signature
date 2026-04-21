#####################################################################################
# ========== Code 6: Agreement metrics & viridis-styled plots ======================
# INPUT : Code 4 output ONLY → "Wei_scRNA_with_profiles_scored.rds"
# OUTPUT: Metrics CSVs + agreement/venn/bar/UMAP plots

suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(viridis)
  library(ggvenn); library(readr); library(tidyr)
})
set.seed(1234); options(stringsAsFactors = FALSE)

# ----- Paths (EDIT THESE) ----------------------------------------------------------
project_root <- "/Users/Paula/Desktop/Wei_2023_Analysis"
in_dir_code4 <- file.path(project_root, "Outputs_Code4_20251130")
# If out_dir not set earlier (e.g., by Code 5), define it here:
if (!exists("out_dir")) out_dir <- file.path(project_root, "Outputs_Code6_20251130")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- Load input (only if not already in memory) ---------------------------------
if (!exists("combined")) {
  rds_path <- file.path(in_dir_code4, "Wei_scRNA_with_profiles_scored.rds")
  stopifnot("Expected input from Code 4 not found." = file.exists(rds_path))
  combined <- readRDS(rds_path)
}
stopifnot(inherits(combined, "Seurat"))

# Pick an available UMAP reduction name already derived in your script (fallback safe)
umap_name <- if ("harmony_umap" %in% names(combined@reductions)) "harmony_umap" else {
  if ("pca_umap" %in% names(combined@reductions)) "pca_umap" else Reductions(combined)[[1]]
}

# ---------------- Helpers ----------------------------------------------------------
.jaccard <- function(a, b) {
  inter <- sum(a & b, na.rm = TRUE); uni <- sum(a | b, na.rm = TRUE)
  if (uni == 0) return(NA_real_); inter / uni
}
.mcc_from_tbl <- function(tbl) {
  rn <- rownames(tbl); cn <- colnames(tbl)
  stopifnot(identical(rn, c("Profile+","Other")), identical(cn, c("Profile+","Other")))
  tp <- as.numeric(tbl["Profile+","Profile+"]); tn <- as.numeric(tbl["Other","Other"])
  fp <- as.numeric(tbl["Profile+","Other"]);    fn <- as.numeric(tbl["Other","Profile+"])
  denom <- sqrt( (tp+fp)*(tp+fn)*(tn+fp)*(tn+fn) )
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  (tp*tn - fp*fn) / denom
}
.kappa_from_tbl <- function(tbl) {
  n <- sum(tbl); if (n == 0) return(NA_real_)
  p0 <- (tbl["Profile+","Profile+"] + tbl["Other","Other"]) / n
  pr <- rowSums(tbl) / n; pc <- colSums(tbl) / n; pe <- sum(pr * pc)
  if (1 - pe == 0) return(NA_real_); (p0 - pe) / (1 - pe)
}
.acc_from_tbl <- function(tbl) {
  n <- sum(tbl); if (n == 0) return(NA_real_)
  (tbl["Profile+","Profile+"] + tbl["Other","Other"]) / n
}
.f1_from_tbl <- function(tbl) {
  tp <- as.numeric(tbl["Profile+","Profile+"]); fp <- as.numeric(tbl["Profile+","Other"])
  fn <- as.numeric(tbl["Other","Profile+"])
  prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  rec  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  if (is.na(prec) || is.na(rec) || (prec + rec) == 0) return(NA_real_)
  2 * prec * rec / (prec + rec)
}

# Ensure factor levels
lv <- c("Profile+","Other")
combined$our_data1.0     <- factor(combined$our_data1.0,     levels = lv)
combined$our_data1.5     <- factor(combined$our_data1.5,     levels = lv)
combined$our_data2.0     <- factor(combined$our_data2.0,     levels = lv)
combined$article_data1.0 <- factor(combined$article_data1.0, levels = lv)
combined$article_data1.5 <- factor(combined$article_data1.5, levels = lv)
combined$article_data2.0 <- factor(combined$article_data2.0, levels = lv)

# ---------------- Metrics for each threshold --------------------------------------
thrs <- c("1.0","1.5","2.0")
metrics <- list()

for (thr in thrs) {
  our_col <- paste0("our_data", thr)
  art_col <- paste0("article_data", thr)
  stopifnot(our_col %in% colnames(combined@meta.data),
            art_col %in% colnames(combined@meta.data))
  
  ok    <- complete.cases(combined@meta.data[, c(our_col, art_col)])
  our_v <- droplevels(combined@meta.data[[our_col]][ok])
  art_v <- droplevels(combined@meta.data[[art_col]][ok])
  
  tbl <- table(our_v, art_v)
  full_tbl <- matrix(0, nrow = 2, ncol = 2, dimnames = list(lv, lv))
  full_tbl[rownames(tbl), colnames(tbl)] <- tbl
  tbl <- full_tbl
  
  our_log <- our_v == "Profile+"; art_log <- art_v == "Profile+"
  
  metrics[[thr]] <- data.frame(
    threshold    = thr,
    n_compared   = sum(tbl),
    TP           = as.numeric(tbl["Profile+","Profile+"]),
    TN           = as.numeric(tbl["Other","Other"]),
    FP           = as.numeric(tbl["Profile+","Other"]),
    FN           = as.numeric(tbl["Other","Profile+"]),
    Accuracy     = .acc_from_tbl(tbl),
    Kappa        = .kappa_from_tbl(tbl),
    MCC          = .mcc_from_tbl(tbl),
    Jaccard      = .jaccard(our_log, art_log),
    F1           = .f1_from_tbl(tbl),
    stringsAsFactors = FALSE
  )
}

metrics_df <- bind_rows(metrics)
readr::write_csv(metrics_df, file.path(out_dir, "ProfileSimilarity_metrics_by_threshold.csv"))
message("✅ Metrics saved → ", file.path(out_dir, "ProfileSimilarity_metrics_by_threshold.csv"))

# ---------------- Continuous agreement (z-score correlations) ---------------------
z_df <- data.frame(
  our_profile_z      = combined$our_profile_z,
  article_profile_z  = combined$article_profile_z
) %>% mutate(ok = complete.cases(.))
pear_r  <- suppressWarnings(cor(z_df$our_profile_z[z_df$ok], z_df$article_profile_z[z_df$ok], method = "pearson"))
spear_r <- suppressWarnings(cor(z_df$our_profile_z[z_df$ok], z_df$article_profile_z[z_df$ok], method = "spearman"))

# Scatter
p_scatter <- ggplot(z_df[z_df$ok, ], aes(x = article_profile_z, y = our_profile_z)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
  theme_classic(base_size = 12) +
  labs(title = sprintf("Profile z-scores (r=%.2f, ρ=%.2f)", pear_r, spear_r),
       x = "Article profile z", y = "Our profile z") +
  scale_color_viridis_d()

ggsave(file.path(out_dir, "ProfileSimilarity_scatter_zscores.png"),
       p_scatter, width = 6.2, height = 5.8, dpi = 300)

# Bland–Altman
ba_df <- z_df[z_df$ok, ]
ba_df$mean <- rowMeans(ba_df[, c("our_profile_z","article_profile_z")])
ba_df$diff <- ba_df$our_profile_z - ba_df$article_profile_z
mu  <- mean(ba_df$diff); sdv <- sd(ba_df$diff)
p_ba <- ggplot(ba_df, aes(x = mean, y = diff)) +
  geom_point(alpha = 0.35, size = 0.7) +
  geom_hline(yintercept = mu, linewidth = 0.8) +
  geom_hline(yintercept = mu + 1.96*sdv, linetype = "dashed") +
  geom_hline(yintercept = mu - 1.96*sdv, linetype = "dashed") +
  theme_classic(base_size = 12) +
  labs(title = "Bland–Altman: Our − Article (z-scores)",
       x = "Mean z-score", y = "Difference (Our − Article)") +
  scale_color_viridis_d()

ggsave(file.path(out_dir, "ProfileSimilarity_bland_altman.png"),
       p_ba, width = 6.2, height = 5.8, dpi = 300)

# ---------------- Per-threshold visuals ------------------------------------------
make_per_threshold_plots <- function(thr) {
  our_col <- paste0("our_data", thr)
  art_col <- paste0("article_data", thr)
  nm_tag  <- gsub("\\.", "", thr)
  
  ok    <- complete.cases(combined@meta.data[, c(our_col, art_col)])
  our_v <- droplevels(combined@meta.data[[our_col]][ok])
  art_v <- droplevels(combined@meta.data[[art_col]][ok])
  
  tbl <- table(our_v, art_v)
  full_tbl <- matrix(0, nrow = 2, ncol = 2, dimnames = list(lv, lv))
  full_tbl[rownames(tbl), colnames(tbl)] <- tbl
  tbl <- as.table(full_tbl)
  
  agreement <- rep(NA_character_, nrow(combined@meta.data))
  agreement[ok & our_v == "Profile+" & art_v == "Profile+"] <- "Both Profile+"
  agreement[ok & our_v == "Profile+" & art_v == "Other"]    <- "Our+ only"
  agreement[ok & our_v == "Other"    & art_v == "Profile+"] <- "Article+ only"
  agreement[ok & our_v == "Other"    & art_v == "Other"]    <- "Both Other"
  combined$agreement_tmp__ <- factor(agreement, levels = c("Both Profile+","Our+ only","Article+ only","Both Other"))
  
  p_umap <- DimPlot(combined, reduction = umap_name, group.by = "agreement_tmp__", pt.size = 0.3) +
    scale_color_viridis_d(na.value = "grey80") +
    labs(title = paste0("Agreement on UMAP (z ≥ ", thr, ")"), color = "Agreement") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(file.path(out_dir, sprintf("ProfileSimilarity_umap_agreement_thr%s.png", nm_tag)),
         p_umap, width = 7.0, height = 6.2, dpi = 300)
  
  cells_all <- rownames(combined@meta.data)[ok]
  our_pos   <- cells_all[our_v == "Profile+"]
  art_pos   <- cells_all[art_v == "Profile+"]
  p_venn <- ggvenn(list(Our = our_pos, Article = art_pos),
                   fill_color = viridis(2), show_percentage = FALSE) +
    ggtitle(paste0("Overlap of Profile+ (z ≥ ", thr, ")"))
  
  ggsave(file.path(out_dir, sprintf("ProfileSimilarity_venn_thr%s.png", nm_tag)),
         p_venn, width = 6.0, height = 6.0, dpi = 300)
  
  df_tbl <- as.data.frame(tbl); colnames(df_tbl) <- c("Our","Article","Count")
  p_bar <- ggplot(df_tbl, aes(x = Our, y = Count, fill = Article)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_viridis_d() +
    theme_classic(base_size = 12) +
    labs(title = paste0("Agreement 2×2 (z ≥ ", thr, ")"),
         x = "Our call", y = "Cell count", fill = "Article call")
  
  ggsave(file.path(out_dir, sprintf("ProfileSimilarity_bar_thr%s.png", nm_tag)),
         p_bar, width = 6.2, height = 5.6, dpi = 300)
  
  combined@meta.data$agreement_tmp__ <- NULL
}
invisible(lapply(thrs, make_per_threshold_plots))

message("✅ Plots saved to: ", out_dir)

# ---------------- Optional: README ------------------------------------------------
readme_txt <- c(
  "# Profile Similarity Outputs",
  sprintf("- Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "- Metrics CSV: ProfileSimilarity_metrics_by_threshold.csv",
  "- Continuous plots: ProfileSimilarity_scatter_zscores.png, ProfileSimilarity_bland_altman.png",
  "- Per-threshold plots:",
  "  * UMAP agreement overlays: ProfileSimilarity_umap_agreement_thr10/15/20.png",
  "  * Venn overlaps: ProfileSimilarity_venn_thr10/15/20.png",
  "  * 2x2 bar plots: ProfileSimilarity_bar_thr10/15/20.png"
)
writeLines(readme_txt, file.path(out_dir, "ProfileSimilarity_README.txt"))

###################################################################################
# ============== Macro-cluster correlation panels (Profile+ only) =================
suppressPackageStartupMessages({ library(ggplot2); library(viridis); library(dplyr); library(readr) })
meta <- combined@meta.data

# Choose the macro-cluster column (auto-detect if 'macro_col' not supplied)
if (!exists("macro_col") || !(macro_col %in% colnames(meta))) {
  candidates <- grep("macro|major|broad|annot", colnames(meta), ignore.case = TRUE, value = TRUE)
  candidates <- unique(c(candidates,
                         intersect(c("macro_cluster","macro_clusters","macro_annot",
                                     "macro_annotation","annotation_macro","broad_class",
                                     "major_class","lineage","macro","annotation"), colnames(meta))))
  pick <- NA_character_
  for (nm in candidates) {
    if (length(unique(na.omit(meta[[nm]]))) == 4) { pick <- nm; break }
  }
  if (is.na(pick)) {
    for (nm in colnames(meta)) {
      if (length(unique(na.omit(meta[[nm]]))) == 4) { pick <- nm; break }
    }
  }
  if (is.na(pick)) stop("Could not auto-detect the macro cluster column. Set macro_col <- '<colname>' and re-run.")
  macro_col <- pick
}
message("Using macro cluster column: ", macro_col)

make_macro_panel <- function(thr) {
  our_col <- paste0("our_data", thr)
  art_col <- paste0("article_data", thr)
  stopifnot(all(c(our_col, art_col) %in% names(meta)))
  
  lv <- c("Profile+","Other")
  meta[[our_col]] <- factor(meta[[our_col]], levels = lv)
  meta[[art_col]] <- factor(meta[[art_col]], levels = lv)
  
  keep <- (meta[[our_col]] == "Profile+" | meta[[art_col]] == "Profile+") &
    !is.na(combined$our_profile_z) & !is.na(combined$article_profile_z)
  
  df <- data.frame(
    macro = factor(meta[[macro_col]], levels = levels(as.factor(meta[[macro_col]]))),
    our_z = combined$our_profile_z,
    art_z = combined$article_profile_z,
    keep = keep
  )
  df <- df[df$keep & !is.na(df$macro), , drop = FALSE]
  if (nrow(df) == 0) stop("No Profile+ cells at threshold ", thr, " for any macro cluster.")
  
  xlim <- as.numeric(quantile(df$art_z, c(0.01, 0.99), na.rm = TRUE))
  ylim <- as.numeric(quantile(df$our_z, c(0.01, 0.99), na.rm = TRUE))
  
  labels <- df %>%
    group_by(macro) %>%
    summarise(
      n   = n(),
      r   = suppressWarnings(cor(our_z, art_z, method = "pearson")),
      rho = suppressWarnings(cor(our_z, art_z, method = "spearman")),
      .groups = "drop"
    ) %>%
    mutate(
      lab = paste0("n=", format(n, big.mark=","), "\n r=", sprintf("%.2f", r),
                   ",  \u03C1=", sprintf("%.2f", rho)),
      x = xlim[2] - 0.05 * diff(xlim),
      y = ylim[2] - 0.06 * diff(ylim)
    )
  
  p <- ggplot(df, aes(x = art_z, y = our_z)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", contour = TRUE, alpha = 0.65) +
    geom_point(alpha = 0.18, size = 0.55, color = "grey10") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9) +
    scale_fill_viridis_c(option = "plasma", name = "Density") +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
    facet_wrap(~ macro, ncol = 2) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0("Profile+ only (z \u2265 ", thr, ") — annotated macro clusters"),
      subtitle = "Cells called Profile+ by our OR article signature, per macro cluster",
      x = "Article profile z", y = "Our profile z"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    geom_text(data = labels, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3.6)
  
  nm_tag <- gsub("\\.", "", thr)
  outfile <- file.path(out_dir, sprintf("ProfileSimilarity_MacroClusters_thr%s_ProfilePlus_scatter.png", nm_tag))
  ggsave(outfile, p, width = 12, height = 9, dpi = 300)
  message("✅ Saved: ", outfile)
  
  readr::write_csv(labels |> select(macro, n, r, rho),
                   file.path(out_dir, sprintf("ProfileSimilarity_MacroClusters_thr%s_metrics.csv", nm_tag)))
}

invisible(lapply(c("1.0","1.5","2.0"), make_macro_panel))
#####################################################################################
