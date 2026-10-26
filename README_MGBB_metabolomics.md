# MGBB Metabolomics (NMR) — Data Cleaning README

## Overview

This document describes the cleaning and imputation pipeline applied to the
MGBB Nightingale NMR metabolomics data, the thresholds used, and the location
of the raw and cleaned files. For the actual code, see `impute_no_panel_info_NMR_MGBB.R`.

> ⚠️ **Important: this script assumes your data are already normalized.**
> Each metabolite must already be standardized (mean 0, SD 1). Negative values are therefore expected and are **not** treated as errors, and **no** log transform or new standardization is applied. The script checks this when it loads your file and warns you if the data don't look standardized.
> **Do not use it on raw concentrations** (mmol/L, g/L, %).
 
---

## Data summary
 
| | |
|---|---|
| Platform | Nightingale Health NMR |
| Samples | 47,991 (one row per sample) |
| Metabolites | 249: lipoprotein subclasses (size, particle concentration, lipid composition), cholesterol, triglycerides, phospholipids, apolipoproteins, fatty acids, amino acids, glycolysis-related metabolites, ketone bodies, creatinine, albumin, GlycA |
| ID columns | `Sample id`, `Biobank Subject ID` |
| Scale | Already standardized: all 249 columns have mean 0 and SD 1 |
| Missing values in the input | 69,286 (0.6% of all values) |
 
**After cleaning:** 47,868 samples × 249 metabolites.
 
| QC outcome | Samples |
|---|---|
| Passed QC (`kept`) | 47,847 |
| One group of measures was placeholder, imputed and kept (`kept_partial_placeholder`) | 21 |
| Whole sample was placeholder, removed (`dropped_placeholder`) | 123 |

---
 
## File locations
 
**Input (raw data):**
`/data/dgag/mgbb_metabolomics/raw_data/1_NGT_QC_final_240801.csv`
Wide format, 47,991 rows × 251 columns (2 ID columns + 249 metabolites).
 
**Script:** `impute_no_panel_info_NMR_MGBB.R`
 
**Outputs:** `/data/dgag/mgbb_metabolomics/NMR_Results`
 
| File | Contents |
|---|---|
| `MGBB_nmr_imputed_final.csv` | **Final cleaned, imputed dataset. Use this file for analysis.** 47,868 samples: 2 ID columns, 249 metabolites, and 3 quality columns (`n_values_imputed`, `pct_values_imputed`, `placeholder_sample`). |
| `MGBB_nmr_sample_qc_status.csv` | Every one of the 47,991 input samples with its QC outcome (`kept`, `kept_partial_placeholder`, `dropped_placeholder`, `dropped_missingness`), placeholder score and % imputed. |
| `MGBB_nmr_feature_columns.txt` | List of the 249 metabolite column names. |
| `MGBB_nmr_flagged_values.csv` | Every placeholder value that was set to missing (sample, metabolite, original value, reason). |
| `MGBB_nmr_qc_dropped_items.csv` | Every sample or metabolite removed, with the reason. |
| `MGBB_nmr_qc_summary.csv` | Counts at each step. |
 
---
 
## Pipeline steps
 
1. **Read the data** and check that each metabolite is standardized (mean ≈ 0, SD ≈ 1). The ID columns are kept as-is.
2. **Flag placeholder samples.** Some samples contain filler values from a failed measurement: many samples carry exactly the same value for a metabolite. A sample is flagged when it holds such repeated values in at least 12 metabolites, the size of a whole group of measures failing together. Only those values are set to missing. A few repeated values per sample are normal ("too low to measure" floor values) and are kept.
3. **Missingness QC.** Metabolites missing in more than 40% of samples are dropped, then samples missing more than 60% of metabolites. Samples that failed almost completely are removed here.
4. **Zero-variance check.** Metabolites with no variation are dropped. None were.
5. **KNN imputation.** Remaining missing values are filled using the 10 most similar samples (k = 10).
6. **Add quality columns** (`n_values_imputed`, `pct_values_imputed`, `placeholder_sample`) and write the outputs.
No log transform and no re-standardization are applied.
 
---
 
## Thresholds used
 
| Parameter | Value | Meaning |
|---|---|---|
| Repeated value | ≥ 20 samples | A value counts as "repeated" if at least 20 samples share it exactly in a metabolite |
| Placeholder cutoff | ≥ 12 metabolites | Chosen automatically at the gap in the data: no sample had exactly 11 |
| Metabolite missingness cutoff | 40% | Metabolites missing in more than 40% of samples are dropped |
| Sample missingness cutoff | 60% | Samples missing more than 60% of metabolites are dropped |
| KNN neighbours (k) | 10 | Samples used to impute a missing value (fixed seed: 12345) |
 
---
 
## Key decisions and known limitations
 
- **The input was already standardized.** All 249 columns have mean 0 and SD 1, and the same values appear in `df_clin_met_mbgg.csv`. It is not known whether a log transform was applied before standardization.
- **Placeholder values:** there is no codebook confirming what the repeated values represent. Flagging samples with repeated values in ≥ 12 metabolites is a data-driven rule, not a certainty.
- **28 samples have more than 30% of their values imputed, 4 have more than 40%, and the maximum is 49.0%.** These samples are kept in the dataset rather than excluded. Their imputed-value percentage is visible in the final file via `pct_values_imputed`, so this is transparent to anyone using the data downstream (e.g. exclude > 30% for sensitive analyses).
- **123 samples were removed** because nearly the whole measurement failed: they held placeholder values in about two-thirds of all metabolites. They are listed in `MGBB_nmr_sample_qc_status.csv` (`dropped_placeholder`).
- **No plate or batch identifier** exists for this NMR data. `DRAW_DATE` is the best available proxy if a batch check is needed.