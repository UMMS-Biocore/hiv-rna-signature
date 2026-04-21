# HIV RNA+ Profile Scoring and Comparison in Wei et al. scRNA-seq

This repository contains an ordered 7-step R workflow for processing the Wei et al. (Immunity 2023) single-cell RNA-seq dataset, annotating macro programs, importing HIV RNA+ labels, and comparing two gene-signature profiles:
- A bulk RNA-seq derived profile (`our` profile)
- The Wei et al. article-derived RNA+ profile (`article` profile)

## Biological and analytical goal
The workflow classifies cells as `Profile+` or `Other` using signature scoring and strict z-score thresholds, then compares agreement between profiles globally, by macro program, and within HIV RNA+ cells.

Core scoring logic implemented in the code:
1. Split each signature into `UP` and `DOWN` gene sets.
2. Compute separate module scores for each set with `AddModuleScore`.
3. Build one signed score per cell: `UP - DOWN`.
4. Standardize signed scores across all cells: `z = scale(signed_score)`.
5. Call profile-positive at increasing strictness: `z >= 1.0`, `z >= 1.5`, `z >= 2.0`.

This makes direct comparison possible between the bulk RNA-seq profile and the Wei et al. profile under the same scoring framework.

## Repository structure
- `1_UMAP_filter_ordered_PG.R`
- `2_Annotation_macro_ordered_PG.R`
- `3_RNA_positive_cells_annotation_ordered_PG.R`
- `4_DEGs_Scorring_ordered_PG.R`
- `5_UMAP_Scorring_Plots_ordered_PG.R`
- `6_Comparison_profiles_ordered_PG.R`
- `7_RNAplus_comparison_ordered_PG.R`

## End-to-end pipeline

### 1) Build Seurat object and UMAP (`1_UMAP_filter_ordered_PG.R`)
- Reads GEO matrix triples (`*_matrix.mtx[.gz]`, features, barcodes).
- Creates per-sample Seurat objects.
- Applies RNA QC filters:
  - `percent.mt < 25`
  - `nFeature_RNA > 200`
  - `500 < nCount_RNA < 10000`
- Merges samples and runs:
  - LogNormalize
  - HVG selection (`nfeatures = 3000`)
  - PCA (`npcs = 50`)
  - Harmony batch correction (`group.by.vars = sample_id`)
  - neighbor graph + clustering (`resolution = 0.6`)
  - UMAP
- Saves combined object and metadata.

### 2) Macro-program annotation (`2_Annotation_macro_ordered_PG.R`)
- Scores four macro programs using curated marker sets:
  - `Proliferation`
  - `Cytotoxic`
  - `AP-1/TNF`
  - `IRF/IFN`
- Computes per-cell average expression score for each program.
- Computes per-cluster z-scores per program.
- Assigns one macro label per cluster using forced argmax.
- Writes cluster-level program tables and annotated UMAP.

### 3) Import HIV RNA+ calls (`3_RNA_positive_cells_annotation_ordered_PG.R`)
- Loads `bc_RNA_HIV.txt` (expects `cell`, `hiv_rna`, optional `copies`).
- Enforces RNA+ rule using `copies >= 2` when available.
- Matches barcodes to Seurat cells by exact cell ID only.
- Adds metadata:
  - `HIV_RNA`
  - `HIV_Copies`
  - `HIV_RNA_source`
- Writes import summaries and an HIV RNA UMAP highlight plot.

### 4) Signature scoring: bulk profile + article profile (`4_DEGs_Scorring_ordered_PG.R`)
Inputs:
- `Clayton_DEG_List.xlsx` for bulk profile (`Gene ID`, `log2FoldChange`)
- `Article_DEG_list.xlsx` for Wei profile (`Gene`, `Average log2 fold change`)

For each profile:
- Build `UP` genes (`log2FC > 0`) and `DOWN` genes (`log2FC < 0`).
- Keep genes present in the object and detected in at least one cell.
- Compute module scores (`OUR_UP/OUR_DOWN`, `ART_UP/ART_DOWN`).
- Compute signed score:
  - `our_score = our_up_score - our_down_score`
  - `article_score = article_up_score - article_down_score`
- Standardize to z-score across cells:
  - `our_profile_z`
  - `article_profile_z`
- Binary calls at three thresholds:
  - `our_data1.0`, `our_data1.5`, `our_data2.0`
  - `article_data1.0`, `article_data1.5`, `article_data2.0`

Note: for the `our` profile, the script includes a fallback signed-z approach if both module scores are unavailable.

### 5) Prevalence summaries and UMAP triptychs (`5_UMAP_Scorring_Plots_ordered_PG.R`)
For each threshold (`1.0`, `1.5`, `2.0`):
- Computes Profile+ prevalence by macro cluster for both datasets.
- Saves prevalence tables and bar plots.
- Creates triptych UMAPs:
  - macro annotations
  - bulk profile highlights
  - article profile highlights

### 6) Agreement metrics and concordance plots (`6_Comparison_profiles_ordered_PG.R`)
Per threshold:
- Confusion matrix components (`TP`, `TN`, `FP`, `FN`)
- Similarity metrics:
  - Accuracy
  - Cohen's kappa
  - MCC
  - Jaccard
  - F1

Also generates:
- continuous z-score scatter and Bland-Altman plots
- UMAP agreement overlays
- Venn overlap plots
- bar plots of 2x2 agreement
- macro-specific correlation panels

### 7) HIV RNA+ specific comparisons (`7_RNAplus_comparison_ordered_PG.R`)
Within HIV RNA+ cells:
- Computes composition across:
  - `Bulk RNA-seq Profile +`
  - `Wei et al. RNA+ Profile +`
  - `Both`
  - `Neither`
- Generates overall and per-macro composition bar plots.
- Produces correlation plots (overall and faceted by macro) for Profile+ union cells across all z thresholds.

## Key metadata fields produced
- Macro annotation: `annotation_macro`
- HIV RNA labels: `HIV_RNA`, `HIV_Copies`, `HIV_RNA_source`
- Bulk profile:
  - `our_up_score`, `our_down_score`, `our_score`, `our_profile_z`
  - `our_data1.0`, `our_data1.5`, `our_data2.0`
- Article profile:
  - `article_up_score`, `article_down_score`, `article_score`, `article_profile_z`
  - `article_data1.0`, `article_data1.5`, `article_data2.0`

## Software dependencies
Main R packages used across scripts:
- `Seurat`
- `harmony`
- `Matrix`
- `dplyr`, `tibble`, `tidyr`, `purrr`, `stringr`
- `readr`, `readxl`
- `ggplot2`, `patchwork`, `cowplot`, `viridis`, `scales`, `ggvenn`

## Required external inputs
- GEO matrix files from `GSE239909_RAW`.
- HIV RNA call table: `bc_RNA_HIV.txt`.
- DEG Excel files:
  - `Clayton_DEG_List.xlsx`
  - `Article_DEG_list.xlsx`

## Reproducibility notes
- All scripts set `set.seed(1234)`.
- Path variables are currently hard-coded and should be updated before running.
- Scripts are designed to run in numeric order (`1` to `7`).

## Known caveats
- Filenames use `Scorring` spelling in multiple script names.
- Script 4 prints one output filename with an apparent typo (`WWei...`) in a message string, while writing the correct stable file (`Wei_scRNA_with_profiles_scored.rds`).

