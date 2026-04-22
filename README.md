# HIV-infected CD4+ T cells that survive NK cell attack — source code and processed single-cell data

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19699921.svg)](https://doi.org/10.5281/zenodo.19699921)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Data License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

This repository contains the source code, inputs, and derived outputs for the secondary single-cell RNA-seq analysis in Grasberger et al., *"HIV-infected CD4+ T cells that survive NK cell attack evade detection by gain of inhibitory MHC-I and loss of activating TRAIL-R2 expression."* The pipeline re-processes the [Wei et al. (Immunity 2023) CITE-seq dataset](https://doi.org/10.1016/j.immuni.2023.10.002) (GEO: [GSE239909](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239909)) of HIV RNA+ CD4+ T cells from people with HIV, and projects onto it a bulk RNA-seq-derived gene signature of *in vitro* HIV-infected CD4+ T cells that survive NK cell attack (the *our* profile, derived from the associated study) alongside the Wei et al. HIV RNA+ signature (the *article* profile).

Both profiles are scored on every cell using a common framework: split each signature into `UP` and `DOWN` gene sets, module-score each with `AddModuleScore`, combine into a signed `UP − DOWN` score, z-standardize across cells, and call cells *Profile+* at increasing thresholds (`z ≥ 1.0`, `1.5`, `2.0`). Profile concordance is then quantified globally, per macro-program, and within the HIV RNA+ subset.

## Repository structure

Source (tracked here):
- [1_UMAP_filter_ordered_PG.R](1_UMAP_filter_ordered_PG.R)
- [2_Annotation_macro_ordered_PG.R](2_Annotation_macro_ordered_PG.R)
- [3_RNA_positive_cells_annotation_ordered_PG.R](3_RNA_positive_cells_annotation_ordered_PG.R)
- [4_DEGs_Scorring_ordered_PG.R](4_DEGs_Scorring_ordered_PG.R)
- [5_UMAP_Scorring_Plots_ordered_PG.R](5_UMAP_Scorring_Plots_ordered_PG.R)
- [6_Comparison_profiles_ordered_PG.R](6_Comparison_profiles_ordered_PG.R)
- [7_RNAplus_comparison_ordered_PG.R](7_RNAplus_comparison_ordered_PG.R)
- [Dockerfile](Dockerfile), [scripts/run_pipeline.sh](scripts/run_pipeline.sh), [scripts/download_data.sh](scripts/download_data.sh), [scripts/zenodo_upload.sh](scripts/zenodo_upload.sh)

Data (archived on Zenodo):
- `GSE239909_RAW/` — raw 10x matrix triples mirrored from GEO
- `Outputs_Code1..7_20251130/` — per-step Seurat RDS objects, metadata, tables, and figure panels
- `bc_RNA_HIV.txt`, `Clayton_DEG_List.xlsx`, `Article_DEG_list.xlsx` — non-GEO inputs

## Pipeline

### 1. Build Seurat object and UMAP ([1_UMAP_filter_ordered_PG.R](1_UMAP_filter_ordered_PG.R))
Reads GEO matrix triples (`*_matrix.mtx[.gz]`, features, barcodes) and builds per-sample Seurat objects. RNA QC filters: `percent.mt < 25`, `nFeature_RNA > 200`, `500 < nCount_RNA < 10000`. Samples are merged, log-normalized, HVG-selected (`nfeatures = 3000`), PCA'd (`npcs = 50`), integrated with Harmony (`group.by.vars = sample_id`), clustered at resolution 0.6, and embedded with UMAP.

### 2. Macro-program annotation ([2_Annotation_macro_ordered_PG.R](2_Annotation_macro_ordered_PG.R))
Four canonical transcriptional programs — Proliferation, Cytotoxic, AP-1/TNF, IRF/IFN — are scored with curated marker sets. Per-cluster z-scores are computed and a single macro label is assigned per cluster via forced argmax.

### 3. Import HIV RNA+ calls ([3_RNA_positive_cells_annotation_ordered_PG.R](3_RNA_positive_cells_annotation_ordered_PG.R))
`bc_RNA_HIV.txt` is loaded (`cell`, `hiv_rna`, optional `copies`). The RNA+ call is enforced with `copies ≥ 2` when copies are present; barcodes are matched to Seurat cells by exact cell ID. Metadata fields `HIV_RNA`, `HIV_Copies`, and `HIV_RNA_source` are written.

### 4. Signature scoring — *our* profile + *article* profile ([4_DEGs_Scorring_ordered_PG.R](4_DEGs_Scorring_ordered_PG.R))
Inputs: `Clayton_DEG_List.xlsx` (bulk RNA-seq DEGs underlying the *our* profile) and `Article_DEG_list.xlsx` (Wei et al.-derived DEGs for the *article* profile). For each signature, `UP` and `DOWN` gene sets are built from the sign of the log2 fold change, restricted to genes present and detected in the object, module-scored, combined into a signed score, and z-standardized:
- `our_score = our_up_score − our_down_score` → `our_profile_z`
- `article_score = article_up_score − article_down_score` → `article_profile_z`

Binary Profile+ calls are made at three z thresholds: `our_data1.0`, `our_data1.5`, `our_data2.0`, and the matching `article_data*`.

### 5. Prevalence summaries and UMAP triptychs ([5_UMAP_Scorring_Plots_ordered_PG.R](5_UMAP_Scorring_Plots_ordered_PG.R))
At each threshold, Profile+ prevalence is computed per macro-cluster for both datasets, with prevalence tables and bar plots. Triptych UMAPs display macro annotations, *our* profile highlights, and *article* profile highlights side by side.

### 6. Concordance metrics ([6_Comparison_profiles_ordered_PG.R](6_Comparison_profiles_ordered_PG.R))
Per threshold, the two profiles are compared via confusion-matrix components (`TP`, `TN`, `FP`, `FN`) and similarity metrics (accuracy, Cohen's κ, MCC, Jaccard, F1). Supporting plots: continuous z-score scatter, Bland-Altman, UMAP agreement overlays, Venn overlaps, 2×2 bar plots, and macro-specific correlation panels.

### 7. HIV RNA+ composition and correlations ([7_RNAplus_comparison_ordered_PG.R](7_RNAplus_comparison_ordered_PG.R))
Within HIV RNA+ cells, composition is computed across Bulk+, Wei+, Both, and Neither, with overall and per-macro bar plots. Correlation panels are produced for Profile+ union cells across all z thresholds.

## Metadata fields produced
- `annotation_macro` ∈ {Proliferation, Cytotoxic, AP-1/TNF, IRF/IFN}
- `HIV_RNA`, `HIV_Copies`, `HIV_RNA_source`
- `our_up_score`, `our_down_score`, `our_score`, `our_profile_z`, `our_data1.0`, `our_data1.5`, `our_data2.0`
- `article_up_score`, `article_down_score`, `article_score`, `article_profile_z`, `article_data1.0`, `article_data1.5`, `article_data2.0`

## Software dependencies
`Seurat`, `harmony`, `Matrix`, `dplyr`, `tibble`, `tidyr`, `purrr`, `stringr`, `readr`, `readxl`, `ggplot2`, `patchwork`, `cowplot`, `viridis`, `scales`, `ggvenn`. Exact package versions used to produce the archived outputs are captured per step in the corresponding `Outputs_CodeN_*/sessionInfo_*.txt`.

## Configuration

Scripts read the project path from the `PROJECT_ROOT` environment variable:

```bash
export PROJECT_ROOT=/absolute/path/to/your/project
```

Inside the Docker image, `PROJECT_ROOT` defaults to `/mnt/project` (the mount point).

## How to run

### Docker (recommended)
```bash
docker build -t hiv-rna-signature:latest .

docker run --rm \
  -v /absolute/path/to/project:/mnt/project \
  hiv-rna-signature:latest all

# Or a single step:
docker run --rm \
  -v /absolute/path/to/project:/mnt/project \
  hiv-rna-signature:latest 4
```

Pass `all`, a digit `1`..`7`, or an exact script filename as the final argument.

### Native R
```bash
export PROJECT_ROOT=/absolute/path/to/your/project
bash scripts/run_pipeline.sh all       # steps 1-7
bash scripts/run_pipeline.sh 4         # step 4 only
Rscript 1_UMAP_filter_ordered_PG.R     # also works directly
```

## Data availability

All bulk data is archived on Zenodo with a permanent DOI. The record is self-contained: raw 10x matrices mirrored from GEO GSE239909, every per-step Seurat RDS and figure panel, and the three non-GEO inputs (`bc_RNA_HIV.txt`, `Clayton_DEG_List.xlsx`, `Article_DEG_list.xlsx`).

| Resource | Link |
|---|---|
| Zenodo concept DOI (always latest version) | https://doi.org/10.5281/zenodo.19699921 |
| Zenodo version DOI (v1.0.0) | https://doi.org/10.5281/zenodo.19699922 |
| GEO raw data (canonical upstream) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE239909 |
| Wei et al. 2023, *Immunity* | https://doi.org/10.1016/j.immuni.2023.10.002 |
| GitHub source | https://github.com/UMMS-Biocore/hiv-rna-signature |

### Reconstitute the working tree
```bash
cp .env.example .env     # then edit to set ZENODO_RECORD_ID
export PROJECT_ROOT=/absolute/path/to/project
cd "$PROJECT_ROOT"
bash scripts/download_data.sh all       # GEO + Zenodo
# bash scripts/download_data.sh geo     # GEO only
# bash scripts/download_data.sh zenodo  # Zenodo only
```

`download_data.sh` verifies every file against the MD5 returned by the Zenodo API. The repo-level `MANIFEST.sha256` covers every shipped file for post-download verification:

```bash
shasum -a 256 -c MANIFEST.sha256
```

## Citation

Use the GitHub "Cite this repository" button (powered by [CITATION.cff](CITATION.cff)), or:

> Grasberger PE, Sondrini AR, Glidden N, Modica A, Pushlar N, Bedir S, Bromfield T, Gentling S, Cheema K, Kucukural A, Ozdemir M, Zapp M, Bosque A, Leyre L, Shulkin A, Piechocka-Trocha A, Jones RB, Clayton KL. *HIV-infected CD4+ T cells that survive NK cell attack evade detection by gain of inhibitory MHC-I and loss of activating TRAIL-R2 expression — source code and processed single-cell data.* Zenodo. https://doi.org/10.5281/zenodo.19699921

Please also cite the Wei et al. (2023) publication and GEO accession GSE239909 when reusing the raw 10x matrices.

## License

- Source code: [MIT](LICENSE)
- Data artifacts on Zenodo: [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/)

