# MGBB Proteomics (Olink): Data Cleaning README

## Overview

This document describes the cleaning and imputation of the MGBB Olink proteomics data: the data, the location of the input and output files, what the pipeline does, and the thresholds used. For the code, see `impute_no_panel_info_Olink_MGBB.R`.

> **Input scale:** Olink NPX values are already on a log2 scale, so no further log transform is applied. Olink's own QC labels (`SampleQC`, `AssayQC`, `LOD`) are used directly.

---

## 1. Requirements

- **R** 4.1 or newer
- **R packages:** `data.table` (CRAN) and `impute` (Bioconductor)

```r
install.packages("data.table")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("impute")
```

**Memory:** the whole export is loaded at once. A file of about 2 million rows (≈1,900 samples × ≈1,000 proteins) needs roughly 4–8 GB of RAM.

## Data summary

| | |
|---|---|
| Platform | Olink Reveal panel (NPX, log2 scale) |
| Delivery | `T2D_Reveal_NPX_2026-07-06` |
| Format | Long: one row per sample × protein (2,161,060 rows) |
| Proteins (assays) | 1,034 |
| Participant samples | 1,892 (plus Olink control wells, removed) |
| Duplicate samples | 13 pairs: the same specimen run twice (IDs ending `_1` / `_2`) |
| Measurements failing Olink QC | 2,742 |
| Measurements below limit of detection (LOD) | 274,732 (14%) |

**After cleaning:** 1,878 samples × 1,024 proteins, with no missing values.
- 13 duplicate pairs collapsed to one sample each.
- 10 proteins removed (9 below LOD in ≥ 99% of samples, 1 missing in > 40%) and 1 sample removed (missing > 60%).
- 1,877 samples map to a `Biobank_Subject_Id` through the sample manifest.

---

## File locations

**Input (raw data):**
`/data/dgag/mgbb_proteomics/raw_data/T2D_Reveal_NPX_2026-07-06.csv`

**Sample manifest** (Olink `SampleID` ↔ `Biobank_Subject_Id`), same folder:
`NGHXXXXX-SampleManifest_template_05.19.26MGBB.xlsx` (columns `Specimen Barcode` ↔ `Cust. Subj.`)

**Script:** `impute_no_panel_info_Olink_MGBB.R`

**Outputs:** `/data/dgag/mgbb_proteomics/OLINK Results/`

| File | Contents |
|---|---|
| `MGBB_olink_imputed_final.csv` | **Final cleaned dataset, z-scored. Use this file for analysis.** 1,878 samples: `SampleID`, 1,024 proteins, and 3 quality columns (`n_values_imputed`, `pct_values_imputed`, `duplicate_pair_correlation`). |
| `MGBB_olink_imputed_final_ranknorm.csv` | Same data, rank-inverse-normal transformed instead of z-scored (for comparison or sensitivity analyses). |
| `MGBB_olink_feature_columns.txt` | List of the 1,024 protein column names. |
| `MGBB_olink_flagged_qc_values.csv` | Every measurement that failed Olink QC or was below LOD (sample, protein, original value, reason). |
| `MGBB_olink_qc_dropped_items.csv` | Every protein or sample removed, with the reason. |
| `MGBB_olink_duplicate_pair_correlations.csv` | One row per duplicate pair: correlation, QC-issue counts, and which run was kept and why. |
| `MGBB_olink_duplicate_samples_flagged.csv` | Both original runs of every duplicate pair, marked kept or dropped. |
| `MGBB_olink_qc_summary.csv` | Counts at each step. |

---

## Pipeline steps

1. **Read the data** and keep only participant samples (`SampleType == "SAMPLE"`). Olink control wells are removed.
2. **Collapse duplicate samples.** Each pair of runs of the same specimen (`ID_1`, `ID_2`, same plate) becomes one sample. The run with fewer QC problems is kept, unless the other run is clearly more typical of the cohort and more deeply sequenced. The `_1`/`_2` suffix is removed so the ID matches the manifest.
3. **Apply Olink QC.** Measurements failing `SampleQC` or `AssayQC` are set to missing. Values below LOD are kept and logged.
4. **Reshape** to one row per sample, one column per protein.
5. **Remove undetectable proteins:** those below LOD in ≥ 99% of samples.
6. **Missingness QC.** Proteins missing in more than 40% of samples are dropped, then samples missing more than 60% of proteins.
7. **KNN imputation** of any remaining missing values (k = 10). None remained in this delivery.
8. **Standardize** each protein (z-score). A rank-inverse-normal version is also written.
9. **Add quality columns** and write the outputs.

---

## Thresholds used

| Parameter | Value | Meaning |
|---|---|---|
| QC pass | `PASS` only | Measurements with `SampleQC` or `AssayQC` other than `PASS` (including `WARN`) are set to missing |
| Below LOD | kept | Values below LOD are kept (and logged), not set to missing |
| Undetectable protein | ≥ 99% below LOD | Protein removed |
| Protein missingness cutoff | 40% | Proteins missing in more than 40% of samples are dropped |
| Sample missingness cutoff | 60% | Samples missing more than 60% of proteins are dropped |
| KNN neighbours (k) | 10 | Samples used to impute a missing value (fixed seed: 12345) |
| Duplicate override | BH q < 0.05, ≥ 60% of proteins agree, deeper sequencing | All three needed to keep the run with more QC problems |
| Low duplicate agreement | r < 0.80 | Flagged only; does not change which run is kept |

---

## Key decisions and known limitations

- **Duplicate pairs are true duplicates** (the same specimen run twice), so each pair is collapsed to one sample rather than dropped. 2 of the 13 pairs agree poorly between runs (r < 0.80); the correlation is kept in `duplicate_pair_correlation` so they can be excluded from sensitive analyses.
- **`WARN` counts as a QC failure.** In this delivery all such values were in the removed protein or sample, so this had no effect on the final data.
- **Below-LOD values are kept.** 149 proteins are below LOD in ≥ 50% of samples. Interpret results for these proteins with care.
- **z-score is the main output.** Rank-normalization spreads near-noise below-LOD values over a large part of the scale, so it is provided for comparison only.
- **The Olink cohort is enriched for T2D** (about 87% with T2D at blood draw), which matters for downstream analyses.