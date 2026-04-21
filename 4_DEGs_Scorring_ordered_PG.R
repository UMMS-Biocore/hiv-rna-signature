#####################################################################################
# ======================= Code 4: Scoring (our + article) ===========================
# INPUT : Code 3 output ONLY → "Wei_scRNA_with_HIV_RNA_exact_annotated.rds"
# OUTPUT: Seurat RDS with our/article scores & z thresholds

suppressPackageStartupMessages({
  library(Seurat); library(readxl); library(readr); library(dplyr); library(Matrix)
})
set.seed(1234)
options(stringsAsFactors = FALSE)

# ----- Paths (EDIT THESE) ----------------------------------------------------------
project_root <- "{Your_Folder}"
# Excel inputs:
deg_path      <- file.path(project_root, "Clayton_DEG_List.xlsx")    # EDIT if stored elsewhere
article_path  <- file.path(project_root, "Article_DEG_list.xlsx")    # EDIT if stored elsewhere
# RDS inputs/outputs:
in_dir_code3  <- file.path(project_root, "Outputs_Code3_20251130")
out_dir       <- file.path(project_root, "Outputs_Code4_20251130")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----- Load input RDS (CODE 3 ONLY) ----------------------------------------------
rds3 <- file.path(in_dir_code3, "Wei_scRNA_with_HIV_RNA_exact_annotated.rds")
stopifnot("Expected input from Code 3 not found." = file.exists(rds3))
combined <- readRDS(rds3)

stopifnot(exists("combined"), inherits(combined, "Seurat"))
DefaultAssay(combined) <- "RNA"

# Seurat v5 hygiene: ensure a single 'data' layer
if (inherits(combined[["RNA"]], "Assay5")) {
  if (!"data" %in% Layers(combined[["RNA"]])) {
    combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)
  }
  combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
}

# ======================= Scorring =======================
# (1) Read OUR DEG table -----------------------------------------------------------
deg <- read_excel(deg_path, .name_repair = "minimal") %>%
  rename(GeneID = `Gene ID`, log2FoldChange = `log2FoldChange`, padj = `padj`)
print(names(deg))
print(head(deg, 10))

stopifnot(exists("deg"), all(c("GeneID","log2FoldChange") %in% names(deg)))

# Drop stale columns (including any old OUR_UP*/OUR_DOWN*)
drop_cols <- grep("^(OUR_UP|OUR_DOWN)", colnames(combined@meta.data), value = TRUE)
drop_cols <- unique(c(drop_cols,
                      "our_up_score","our_down_score","our_score",
                      "our_profile_z","our_profile_call","our_data1.0","our_data1.5","our_data2.0"))
if (length(drop_cols) > 0) {
  combined@meta.data[, drop_cols] <- list(NULL)
}

# 1) Build UP/DOWN lists from your DEG tibble
deg_tbl <- deg %>%
  transmute(
    gene = trimws(as.character(.data[["GeneID"]])),
    lfc  = as.numeric(.data[["log2FoldChange"]])
  ) %>%
  filter(!is.na(gene), nzchar(gene), !is.na(lfc))

up_genes   <- unique(deg_tbl$gene[deg_tbl$lfc >  0])
down_genes <- unique(deg_tbl$gene[deg_tbl$lfc <  0])

# 2) Map to object & keep genes detected in at least 1 cell
M <- if (inherits(combined[["RNA"]], "Assay5")) {
  LayerData(combined[["RNA"]], layer = "data")
} else {
  GetAssayData(combined, assay = "RNA", slot = "data")
}
present <- rownames(M)
keep_expr <- function(g) { g <- intersect(g, present); g[Matrix::rowSums(M[g, , drop = FALSE] > 0) > 0] }
up_use   <- keep_expr(up_genes)
down_use <- keep_expr(down_genes)

cat(sprintf("After guards: UP=%d, DOWN=%d\n", length(up_use), length(down_use)))

# 3) Run AddModuleScore with robust column capture
pick_added <- function(obj, before, prefix){
  added <- setdiff(colnames(obj@meta.data), before)
  added <- added[grepl(paste0("^", prefix), added)]
  if (length(added) > 0) added[1] else NA_character_
}

# UP score
before <- colnames(combined@meta.data)
if (length(up_use) > 0) {
  combined <- AddModuleScore(combined, features = list(up_use), name = "OUR_UP",
                             assay = "RNA", nbin = 24, ctrl = 100, seed = 1234)
}
up_col <- pick_added(combined, before, "OUR_UP")
combined$our_up_score <- if (!is.na(up_col)) combined@meta.data[[up_col]] else NA_real_

# DOWN score
before <- colnames(combined@meta.data)
if (length(down_use) > 0) {
  combined <- AddModuleScore(combined, features = list(down_use), name = "OUR_DOWN",
                             assay = "RNA", nbin = 24, ctrl = 100, seed = 1234)
}
down_col <- pick_added(combined, before, "OUR_DOWN")
combined$our_down_score <- if (!is.na(down_col)) combined@meta.data[[down_col]] else NA_real_

cat(sprintf("Captured columns → UP: %s | DOWN: %s\n",
            ifelse(is.na(up_col),"NA",up_col), ifelse(is.na(down_col),"NA",down_col)))

# 4) One signed score (UP − DOWN) and z-standardize across cells (fallback safe)
fallback_signed_z <- function() {
  g_all <- unique(c(up_use, down_use))
  if (length(g_all) == 0) {
    return(rep(0, ncol(combined)))
  } else {
    combined <<- ScaleData(combined, features = g_all, do.center = TRUE, do.scale = TRUE, verbose = FALSE)
    Z <- if (inherits(combined[["RNA"]], "Assay5")) {
      LayerData(combined[["RNA"]], layer = "scale.data")
    } else {
      GetAssayData(combined, assay = "RNA", slot = "scale.data")
    }
    Zu <- if (length(up_use))   Z[up_use,   , drop = FALSE] else NULL
    Zd <- if (length(down_use)) Z[down_use, , drop = FALSE] else NULL
    signed <- do.call(rbind, c(list(Zu), if (!is.null(Zd)) list(-Zd)))
    prof <- if (!is.null(signed)) Matrix::colMeans(as.matrix(signed), na.rm = TRUE) else rep(0, ncol(combined))
    sdv <- stats::sd(prof, na.rm = TRUE)
    if (is.finite(sdv) && sdv > 0) as.numeric(scale(prof)) else rep(0, length(prof))
  }
}

if (all(is.na(combined$our_up_score)) && all(is.na(combined$our_down_score))) {
  warning("Both OUR_UP and OUR_DOWN scores are NA; using per-gene signed z-average fallback.")
  combined$our_score     <- NA_real_
  combined$our_profile_z <- fallback_signed_z()
} else {
  if (all(is.na(combined$our_down_score))) {
    message("DOWN score missing → using UP only.")
    combined$our_score <- combined$our_up_score
  } else if (all(is.na(combined$our_up_score))) {
    message("UP score missing → using -DOWN only.")
    combined$our_score <- -combined$our_down_score
  } else {
    combined$our_score <- combined$our_up_score - combined$our_down_score
  }
  sdv <- stats::sd(combined$our_score, na.rm = TRUE)
  if (is.finite(sdv) && sdv > 0) {
    combined$our_profile_z <- as.numeric(scale(combined$our_score))
  } else {
    combined$our_profile_z <- fallback_signed_z()
  }
}

# 5) Calls at z ≥ {1.0, 1.5, 2.0}
bin_call <- function(z, thr) factor(ifelse(z >= thr, "Profile+", "Other"),
                                    levels = c("Profile+","Other"))
combined$our_data1.0 <- bin_call(combined$our_profile_z, 1.0)
combined$our_data1.5 <- bin_call(combined$our_profile_z, 1.5)
combined$our_data2.0 <- bin_call(combined$our_profile_z, 2.0)
combined$our_profile_call <- combined$our_profile_z > 0

cat("\nour_profile_z summary:\n"); print(summary(combined$our_profile_z))
cat("\nCounts @ z ≥ 1.0:\n"); print(table(combined$our_data1.0))
cat("\nCounts @ z ≥ 1.5:\n"); print(table(combined$our_data1.5))
cat("\nCounts @ z ≥ 2.0:\n"); print(table(combined$our_data2.0))
#####################################################################################################

# ===================== Article signature: UP/DOWN → signed score → z → thresholds
DefaultAssay(combined) <- "RNA"

# --- Seurat v5 hygiene: ensure a single 'data' layer to read ---
if (inherits(combined[["RNA"]], "Assay5")) {
  if (!"data" %in% Layers(combined[["RNA"]])) {
    combined <- NormalizeData(combined, normalization.method = "LogNormalize",
                              scale.factor = 1e4, verbose = FALSE)
  }
  combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
}

# --------- 0) Read Article_DEG_list.xlsx with explicit columns ----------
# robust numeric parser (handles decimal commas & Unicode minus)
to_num <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  s <- as.character(x)
  s <- gsub("\\s+", "", s)
  s <- chartr(",", ".", s)                                  # decimal commas → dot
  s <- gsub("[\u2212\u2010-\u2014]", "-", s)               # Unicode minus/dashes → "-"
  suppressWarnings(as.numeric(s))
}
clean_sym <- function(x) {
  y <- trimws(as.character(x))
  gsub("\\s+", "", y)
}

art_raw <- read_excel(article_path, .name_repair = "minimal")
nm <- names(art_raw)
stopifnot(all(c("Gene","Average log2 fold change") %in% nm))

art_tbl <- art_raw %>%
  transmute(
    gene = clean_sym(.data[["Gene"]]),
    lfc  = to_num(.data[["Average log2 fold change"]])
  ) %>%
  filter(!is.na(gene), nzchar(gene), !is.na(lfc))

# --------- 1) Build UP/DOWN lists ----------
article_up_genes   <- unique(art_tbl$gene[art_tbl$lfc >  0])
article_down_genes <- unique(art_tbl$gene[art_tbl$lfc <  0])

# --------- 2) Map to object & keep genes detected in ≥1 cell ----------
M <- if (inherits(combined[["RNA"]], "Assay5")) LayerData(combined[["RNA"]], layer="data") else
  GetAssayData(combined, assay="RNA", slot="data")
present <- rownames(M)

keep_expr <- function(g) {
  g <- intersect(g, present)
  if (!length(g)) return(character(0))
  g[Matrix::rowSums(M[g, , drop = FALSE] > 0) > 0]
}
art_up_use   <- keep_expr(article_up_genes)
art_down_use <- keep_expr(article_down_genes)

message(sprintf("Article mapping after guards: UP=%d, DOWN=%d", length(art_up_use), length(art_down_use)))
stopifnot("Article UP or DOWN set is empty after mapping; cannot proceed with exact pipeline." =
            length(art_up_use) > 0 && length(art_down_use) > 0)

# --------- Clean previous article_* columns we will recreate ----------
drop_cols <- grep("^(ART_UP|ART_DOWN)", colnames(combined@meta.data), value = TRUE)
drop_cols <- unique(c(drop_cols,
                      "article_up_score","article_down_score","article_score",
                      "article_profile_z","article_profile_call","article_data1.0","article_data1.5","article_data2.0"))
if (length(drop_cols) > 0) {
  combined@meta.data[, drop_cols] <- list(NULL)
}

# helper to capture the column AddModuleScore creates (e.g., ART_UP1)
pick_added <- function(obj, before, prefix){
  added <- setdiff(colnames(obj@meta.data), before)
  added <- added[grepl(paste0("^", prefix), added)]
  if (length(added) > 0) added[1] else NA_character_
}

# --------- 3) Run AddModuleScore twice (UP, DOWN) ----------
# UP
before <- colnames(combined@meta.data)
combined <- AddModuleScore(combined, features = list(art_up_use), name = "ART_UP",
                           assay = "RNA", nbin = 24, ctrl = 100, seed = 1234)
art_up_col <- pick_added(combined, before, "ART_UP")
stopifnot(!is.na(art_up_col))
combined$article_up_score <- combined@meta.data[[art_up_col]]

# DOWN
before <- colnames(combined@meta.data)
combined <- AddModuleScore(combined, features = list(art_down_use), name = "ART_DOWN",
                           assay = "RNA", nbin = 24, ctrl = 100, seed = 1234)
art_down_col <- pick_added(combined, before, "ART_DOWN")
stopifnot(!is.na(art_down_col))
combined$article_down_score <- combined@meta.data[[art_down_col]]

# --------- 4) Single signed score (UP − DOWN) & z-score across cells ----------
combined$article_score <- combined$article_up_score - combined$article_down_score
sdv <- stats::sd(combined$article_score, na.rm = TRUE)
stopifnot(is.finite(sdv), sdv > 0)  # exact pipeline; no fallback here
combined$article_profile_z <- as.numeric(scale(combined$article_score))

# --------- 5) Binary calls at z ≥ 1.0, 1.5, 2.0 ----------
bin_call <- function(z, thr) factor(ifelse(z >= thr, "Profile+", "Other"),
                                    levels = c("Profile+","Other"))
combined$article_data1.0 <- bin_call(combined$article_profile_z, 1.0)
combined$article_data1.5 <- bin_call(combined$article_profile_z, 1.5)
combined$article_data2.0 <- bin_call(combined$article_profile_z, 2.0)
combined$article_profile_call <- combined$article_profile_z > 0

# --------- Quick summaries ----------
cat("\narticle_profile_z summary:\n"); print(summary(combined$article_profile_z))
cat("\nArticle counts @ z ≥ 1.0:\n"); print(table(combined$article_data1.0))
cat("\nArticle counts @ z ≥ 1.5:\n"); print(table(combined$article_data1.5))
cat("\nArticle counts @ z ≥ 2.0:\n"); print(table(combined$article_data2.0))

# ===== Save updated Seurat object (RDS) =====
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
stopifnot(exists("combined"), inherits(combined, "Seurat"))

rds_name <- sprintf("Wei_scRNA_with_profiles_%s.rds",
                    format(Sys.time(), "%Y%m%d_%H%M%S"))
rds_path <- file.path(out_dir, rds_name)
saveRDS(combined, rds_path)
# Stable copy for chaining Code 5
saveRDS(combined, file.path(out_dir, "Wei_scRNA_with_profiles_scored.rds"))

message("✅ Saved timestamped: ", rds_path)
message("✅ Saved stable copy: ", file.path(out_dir, "WWei_scRNA_with_profiles_scored.rds"))
#####################################################################################
