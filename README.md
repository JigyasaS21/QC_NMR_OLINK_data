# MGBB multi-omics QC: NMR metabolomics and Olink proteomics

Quality-control and imputation pipelines for the Mass General Brigham Biobank (MGBB) **Nightingale NMR metabolomics** and **Olink proteomics** data. Each pipeline takes the delivered data and produces a clean, analysis-ready table (one row per sample, one column per metabolite or protein, no missing values), together with logs of every value or sample that was flagged or removed.

Downstream analyses (association tests, signatures, prediction models) are not part of this repository.

---

| Pipeline | Script | Guide |
|---|---|---|
| NMR metabolomics (Nightingale) | [`nmr/impute_no_panel_info_NMR_MGBB.R`](nmr/impute_no_panel_info_NMR_MGBB.R) | [`nmr/README_NMR_MGBB.md`](nmr/README_NMR_MGBB.md) |
| Olink proteomics (NPX) | [`olink/impute_no_panel_info_Olink_MGBB.R`](olink/impute_no_panel_info_Olink_MGBB.R) | [`olink/README_Olink_MGBB.md`](olink/README_Olink_MGBB.md) |

---

## What each pipeline does

**NMR metabolomics** (input already standardized: mean 0, SD 1 per metabolite)
1. Checks that the input is standardized
2. Detects placeholder (filler) values from failed measurements and sets them to missing
3. Removes metabolites and samples with too much missing data
4. Fills remaining gaps by k-nearest-neighbour imputation (k = 10)
5. Writes the cleaned table, a QC status for every sample, and QC logs

**Olink proteomics** (NPX, log2 scale)
1. Keeps participant samples only (removes Olink control wells)
2. Collapses duplicate runs of the same specimen to one sample
3. Applies Olink's own QC flags (`SampleQC`, `AssayQC`) and logs values below the limit of detection
4. Removes undetectable proteins and proteins/samples with too much missing data
5. Imputes any remaining gaps (k = 10) and standardizes each protein (z-score; rank-inverse-normal version also written)
6. Writes the cleaned tables and QC logs

| | NMR | Olink |
|---|---|---|
| Input | 47,991 samples × 249 metabolites | 1,892 samples × 1,034 proteins |
| **Output** | **47,868 samples × 249 metabolites** | **1,878 samples × 1,024 proteins** |
| Main output file | `MGBB_nmr_imputed_final.csv` | `MGBB_olink_imputed_final.csv` |

See each pipeline's guide for details, thresholds and output files.

---

## Quick start

**Requirements:** R ≥ 4.1 with `data.table` (CRAN) and `impute` (Bioconductor).

```r
install.packages("data.table")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("impute")
```

**Run a pipeline:**
1. Open the script and edit the input and output paths in the `CONFIG` section at the top.
2. Run it in a fresh R session:
   ```bash
   Rscript nmr/impute_no_panel_info_NMR_MGBB.R
   Rscript olink/impute_no_panel_info_Olink_MGBB.R
   ```
3. Check the console summary, then use the `*_imputed_final.csv` file for analysis.

---

## Data

**No data are stored in this repository**, only code and documentation. The input and output files contain participant-level data and are kept on MGB-approved storage. File locations are listed in each pipeline's guide.
