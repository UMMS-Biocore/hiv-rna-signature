#############################################################################################
# ========================== Code 7: HIV_RNA+ composition & correlations ===================
# INPUT : Code 4 output ONLY → "Wei_scRNA_with_profiles_scored.rds"
# OUTPUT: HIV_RNA+ composition barplots + correlation panels (overall & per macro)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr); library(tidyr); library(ggplot2)
  library(cowplot); library(viridis); library(readr)
})

set.seed(1234)
options(stringsAsFactors = FALSE)

# ----- Paths (EDIT THESE) ----------------------------------------------------------
project_root <- "/Users/Paula/Desktop/Wei_2023_Analysis"
in_dir_code4 <- file.path(project_root, "Outputs_Code4_20251130")
out_dir      <- file.path(project_root, "Outputs_Code7_20251130")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- Load input (only if not already in memory) ---------------------------------
if (!exists("combined")) {
  rds_path <- file.path(in_dir_code4, "Wei_scRNA_with_profiles_scored.rds")
  stopifnot("Expected input from Code 4 not found." = file.exists(rds_path))
  combined <- readRDS(rds_path)
}
stopifnot(inherits(combined, "Seurat"))

# ---------- helpers ----------
norm_str <- function(x) tolower(trimws(as.character(x)))

# ---------- inputs ----------
meta <- combined@meta.data
colnames(meta)  # (for quick console glance)

is_pos   <- norm_str(meta$HIV_RNA) %in% c("positive","pos","rna+","positive_rna")
our_cols <- c("our_data1.0","our_data1.5","our_data2.0")
art_cols <- c("article_data1.0","article_data1.5","article_data2.0")
stopifnot(all(our_cols %in% colnames(meta)),
          all(art_cols %in% colnames(meta)),
          "annotation_macro" %in% colnames(meta))

# macro clusters
macros <- levels(factor(meta$annotation_macro))
group_levels <- c("Bulk RNA-seq Profile +",
                  "Wei et al. RNA+ Profile +",
                  "Both",
                  "Neither")

# core counter: returns tibble for one mask (subset) and one cutoff
count_for <- function(mask, our_col, art_col, cut_label) {
  total  <- sum(mask & is_pos, na.rm = TRUE)
  our_hit <- norm_str(meta[[our_col]]) == "profile+"
  art_hit <- norm_str(meta[[art_col]]) == "profile+"
  
  our_hit <- our_hit & mask & is_pos
  art_hit <- art_hit & mask & is_pos
  
  n_our_only <- sum(our_hit & !art_hit, na.rm = TRUE)
  n_art_only <- sum(!our_hit & art_hit, na.rm = TRUE)
  n_both     <- sum(our_hit &  art_hit, na.rm = TRUE)
  n_neither  <- total - (n_our_only + n_art_only + n_both)
  
  ns  <- c(n_our_only, n_art_only, n_both, n_neither)
  pct <- 100 * ns / ifelse(total == 0, NA, total)
  
  tibble(
    cutoff = cut_label,
    group  = factor(group_levels, levels = group_levels),
    n   = ns,
    pct = pct
  )
}

# plotting: bars + viridis (no error bars)
make_bar <- function(df, title_txt) {
  ggplot(df, aes(x = group, y = pct, fill = group)) +
    geom_col(width = 0.72, color = "black") +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%s)", pct, scales::comma(n))),
              vjust = -0.15, lineheight = 0.95, size = 3.3) +
    scale_fill_viridis_d(option = "D") +
    scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = title_txt,
      subtitle = "Percentage distribution of Profile+ calls in HIV_RNA+ cells",
      x = NULL, y = "Percentage (%)"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(angle = 15, hjust = 1))
}

# ---------- BUILD & SAVE ----------
cut_labels <- c("z ≥ 1.0","z ≥ 1.5","z ≥ 2.0")
for (i in seq_along(cut_labels)) {
  cutlab  <- cut_labels[i]
  our_col <- our_cols[i]
  art_col <- art_cols[i]
  
  # overall (all HIV_RNA+ cells)
  overall_tbl <- count_for(mask = rep(TRUE, nrow(meta)), our_col, art_col, cutlab)
  p_overall <- make_bar(overall_tbl, paste0("HIV_RNA+ — ", cutlab, " (All macro clusters)"))
  ggsave(file.path(out_dir, paste0("HIVRNApos_pct_overall_", gsub("[ ≥]", "", cutlab), ".tiff")),
         p_overall, device = "tiff", width = 6, height = 5, units = "in", dpi = 300, compression = "none")
  
  # per-macro plots (2x2 grid per cutoff)
  macro_plots <- vector("list", length(macros)); names(macro_plots) <- macros
  for (m in macros) {
    tbl_m <- count_for(mask = meta$annotation_macro == m, our_col, art_col, cutlab)
    macro_plots[[m]] <- make_bar(tbl_m, paste0("HIV_RNA+ — ", cutlab, " — ", m))
  }
  
  grid_fig <- plot_grid(plotlist = macro_plots, ncol = 2)
  ggsave(file.path(out_dir, paste0("HIVRNApos_pct_macros_", gsub("[ ≥]", "", cutlab), ".tiff")),
         grid_fig, device = "tiff", width = 12, height = 10, units = "in", dpi = 300, compression = "none")
}

#################################################################################################
## =================== Correlations in HIV_RNA+ cells (Profile+ union) ===================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(viridis); library(readr) })

stopifnot(exists("combined"), inherits(combined, "Seurat"))
stopifnot(exists("out_dir") && dir.exists(out_dir))

meta <- combined@meta.data
norm_str <- function(x) tolower(trimws(as.character(x)))

# ---- columns ----
our_z <- "our_profile_z"; art_z <- "article_profile_z"
stopifnot(all(c(our_z, art_z, "HIV_RNA") %in% colnames(meta)))

# ---- HIV_RNA+ mask ----
is_rna_pos <- norm_str(meta$HIV_RNA) %in% c("positive","pos","rna+","positive_rna")

# ---- macro cluster column (auto-detect if not supplied) ----
if (!exists("macro_col") || !(macro_col %in% colnames(meta))) {
  candidates <- grep("macro|major|broad|annot", colnames(meta), ignore.case = TRUE, value = TRUE)
  candidates <- unique(c(candidates,
                         intersect(c("macro_cluster","macro_clusters","macro_annot",
                                     "macro_annotation","annotation_macro","broad_class",
                                     "major_class","lineage","macro","annotation"), colnames(meta))))
  pick <- NA_character_
  for (nm in candidates) if (length(unique(na.omit(meta[[nm]]))) == 4) { pick <- nm; break }
  if (is.na(pick)) for (nm in colnames(meta)) if (length(unique(na.omit(meta[[nm]]))) == 4) { pick <- nm; break }
  if (is.na(pick)) stop("Could not auto-detect the macro cluster column. Set macro_col <- '<colname>' and re-run.")
  macro_col <- pick
}
message("Using macro cluster column: ", macro_col)

# ---- aesthetics helpers (same style) ----
make_corr_plot <- function(df, title_txt) {
  pear  <- suppressWarnings(cor(df$our_z, df$art_z, method = "pearson"))
  spear <- suppressWarnings(cor(df$our_z, df$art_z, method = "spearman"))
  ggplot(df, aes(x = art_z, y = our_z)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                    contour = TRUE, alpha = 0.65) +
    geom_point(alpha = 0.22, size = 0.65, color = "grey10") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "deepskyblue3") +
    scale_fill_viridis_c(option = "plasma", name = "Density") +
    theme_minimal(base_size = 14) +
    labs(
      title = sprintf("%s — n=%s", title_txt, format(nrow(df), big.mark=",")),
      subtitle = sprintf("Pearson r = %.2f   (Spearman \u03C1 = %.2f)", pear, spear),
      x = "Article profile z-score",
      y = "Our profile z-score"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

# ---- builder per threshold (overall + macro panel), HIV_RNA+ + Profile+ union ----
make_corr_hivpos <- function(thr) {
  our_call <- paste0("our_data", thr)
  art_call <- paste0("article_data", thr)
  stopifnot(all(c(our_call, art_call) %in% names(meta)))
  
  keep <- (meta[[our_call]] == "Profile+" | meta[[art_call]] == "Profile+") & is_rna_pos &
    !is.na(meta[[our_z]]) & !is.na(meta[[art_z]])
  
  df <- data.frame(
    macro = factor(meta[[macro_col]]),
    our_z = meta[[our_z]],
    art_z = meta[[art_z]],
    keep  = keep
  ) %>% filter(keep, !is.na(macro))
  
  if (nrow(df) == 0) { warning("No HIV_RNA+ Profile+ cells at threshold ", thr, ". Skipping."); return(invisible(NULL)) }
  
  xlim <- as.numeric(quantile(df$art_z, c(0.01, 0.99), na.rm = TRUE))
  ylim <- as.numeric(quantile(df$our_z, c(0.01, 0.99), na.rm = TRUE))
  
  # overall plot
  p_overall <- make_corr_plot(df, paste0("HIV_RNA+ Profile+ (z \u2265 ", thr, ")"))
  ggsave(file.path(out_dir, sprintf("ProfileSimilarity_corr_overall_thr%s_HIVRNApos_ProfilePlus.png",
                                    gsub("\\.", "", thr))),
         p_overall, width = 7, height = 6, dpi = 300)
  
  # macro faceted panel + per-facet labels
  lab_df <- df %>%
    group_by(macro) %>%
    summarise(
      n = n(),
      r = suppressWarnings(cor(our_z, art_z, method = "pearson")),
      rho = suppressWarnings(cor(our_z, art_z, method = "spearman")),
      .groups = "drop"
    ) %>%
    mutate(lab = paste0("n=", format(n, big.mark=","), "\n r=", sprintf("%.2f", r),
                        ",  \u03C1=", sprintf("%.2f", rho)),
           x = xlim[2] - 0.05 * diff(xlim),
           y = ylim[2] - 0.06 * diff(ylim))
  
  p_macro <- ggplot(df, aes(x = art_z, y = our_z)) +
    stat_density_2d(aes(fill = after_stat(level)), geom = "polygon",
                    contour = TRUE, alpha = 0.65) +
    geom_point(alpha = 0.18, size = 0.55, color = "grey10") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, color = "deepskyblue3") +
    scale_fill_viridis_c(option = "plasma", name = "Density") +
    coord_cartesian(xlim = xlim, ylim = ylim, expand = TRUE) +
    facet_wrap(~ macro, ncol = 2) +
    theme_minimal(base_size = 13) +
    labs(
      title = paste0("HIV_RNA+ Profile+ (z \u2265 ", thr, ") — macro clusters"),
      subtitle = "Cells called Profile+ by our OR article signature",
      x = "Article profile z-score", y = "Our profile z-score"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    ) +
    geom_text(data = lab_df, aes(x = x, y = y, label = lab),
              inherit.aes = FALSE, hjust = 1, vjust = 1, size = 3.6)
  
  ggsave(file.path(out_dir, sprintf("ProfileSimilarity_corr_macros_thr%s_HIVRNApos_ProfilePlus.png",
                                    gsub("\\.", "", thr))),
         p_macro, width = 12, height = 9, dpi = 300)
}

invisible(lapply(c("1.0","1.5","2.0"), make_corr_hivpos))
message("✅ HIV_RNA+ plots saved to: ", out_dir)
#############################################################################################
