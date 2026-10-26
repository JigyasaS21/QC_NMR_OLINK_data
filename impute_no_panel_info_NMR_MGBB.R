# ---------------------------------------------------------------------------
# impute_no_panel_info_NMR_MGBB.R
#
# Quality control and imputation for Nightingale NMR metabolomics data
# (wide format: one row per sample, one column per metabolite).
#
# INPUT ASSUMPTION -- DATA ARE ALREADY NORMALIZED
#   This script assumes the input values are already standardized (z-scores:
#   each metabolite has mean 0 and SD 1). Negative values are therefore
#   expected and are NOT treated as errors, and NO log transform or new
#   standardization is applied. The script checks this assumption on loading
#   and warns if the columns do not look standardized. Do not use this script
#   on raw concentrations (mmol/L, g/L, %).
#
# STEPS
#   1. Read data; check that the input looks standardized
#   2. Detect placeholder samples: samples holding values that exactly match a
#      value shared by many samples, in at least PLACEHOLDER_MIN_COLS
#      metabolites. This catches whole failed samples and samples where one
#      group of measures failed (e.g. all fatty-acid or all HDL-cholesterol
#      measures). Only those repeated cells are set to NA; the sample's other
#      values are kept. A few repeated values per sample (e.g. a detection floor
#      in some metabolites) are normal and kept.
#      The cutoff is chosen automatically at the first gap (a count with no
#      samples) above the normal range; set PLACEHOLDER_MIN_COLS to a number
#      to override. Always check the printed distribution.
#   3. Missingness QC: drop metabolites, then samples, with too many missing
#      values; drop constant metabolites
#   4. Impute remaining missing values by k-nearest neighbours
#   5. Write the imputed table plus per-sample quality columns and QC logs
#
# OUTPUT (in OUT_DIR, prefixed with OUT_PREFIX; .csv by default, .tsv if OUT_FORMAT = "tsv")
#   <prefix>_imputed_final.csv      ID columns + metabolites + quality columns:
#                                     n_values_imputed, pct_values_imputed,
#                                     placeholder_sample (1 = flagged in step 2)
#   <prefix>_feature_columns.txt    metabolite columns in the final table
#   <prefix>_flagged_values.csv     every cell set to NA in step 2
#   <prefix>_qc_dropped_items.csv   metabolites / samples removed, with reason
#   <prefix>_sample_qc_status.csv   EVERY input sample with its QC outcome:
#                                     kept | kept_partial_placeholder |
#                                     dropped_placeholder | dropped_missingness
#   <prefix>_qc_summary.csv         counts at each step
#
# Suggested downstream use: exclude samples with high pct_values_imputed
# (e.g. > 30%) or placeholder_sample == 1 from sensitive analyses.
#
# REQUIRES: R >= 4.1, data.table (CRAN), impute (Bioconductor)
#   install.packages("data.table")
#   BiocManager::install("impute")
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(impute)
})

# =============================================================================
# CONFIG -- edit these
# =============================================================================
INPUT_FILE <- "/Path/to/input/file"     # wide file, comma- or tab-separated
ID_COLS    <- c("Sample id", "Biobank Subject ID")   # identifier columns, passed through unchanged;
# the FIRST is the sample-level ID (must be unique)
OUT_DIR    <- "/Path/to/OutputDir"
OUT_PREFIX <- "MGBB_nmr"
OUT_FORMAT <- "csv"                       

# placeholder-sample detection
PLACEHOLDER_MIN_REPEAT <- 20   # a value is "repeated" if >= this many samples share it exactly in a metabolite
PLACEHOLDER_MIN_COLS   <- "auto"   # a sample is a placeholder if it holds repeated values in >= this many metabolites.
# "auto" = first gap above the normal range (see below); or give a number, e.g. 12
PLACEHOLDER_AUTO_START_QUANTILE <- 0.995   # gap search starts above the count reached by 99.5% of samples (the normal range)
PLACEHOLDER_AUTO_MIN_START      <- 5       # ... and never below this count
PLACEHOLDER_AUTO_FALLBACK       <- 12      # used if no gap is found (smallest Nightingale measure group that fails together)

# missingness QC
FEATURE_MISS_MAX <- 0.40   # drop metabolites missing in > 40% of samples
SAMPLE_MISS_MAX  <- 0.60   # then drop samples missing > 60% of remaining metabolites

# imputation
KNN_K    <- 10             # neighbours for impute.knn
KNN_SEED <- 12345          # fixed seed for reproducibility

# check of the input assumption (warn if columns are not close to mean 0, SD 1)
STANDARDIZED_TOL <- 0.05

# =============================================================================
# Helpers
# =============================================================================
check_standardized <- function(dat, feature_cols, tol) {
  mat <- as.matrix(dat[, ..feature_cols])
  col_means <- colMeans(mat, na.rm = TRUE)
  col_sds   <- apply(mat, 2, sd, na.rm = TRUE)
  message(sprintf("  input check: median column mean %.4f | median column SD %.4f | %.1f%% of values negative",
                  median(col_means), median(col_sds), 100 * mean(mat < 0, na.rm = TRUE)))
  off <- abs(col_means) > tol | abs(col_sds - 1) > tol
  if (any(off)) {
    warning(sum(off), " of ", length(feature_cols), " metabolite columns are not close to mean 0 / SD 1. ",
            "This script assumes standardized (z-scored) input -- check the input file before continuing.")
  }
  invisible(off)
}

# first count with no samples, searching upward from just above the normal range
choose_placeholder_cutoff <- function(n_rep,
                                      start_q  = PLACEHOLDER_AUTO_START_QUANTILE,
                                      min_start = PLACEHOLDER_AUTO_MIN_START,
                                      fallback = PLACEHOLDER_AUTO_FALLBACK) {
  counts <- tabulate(n_rep + 1L, nbins = max(n_rep) + 1L)       # counts[k + 1] = samples with k repeated values
  start  <- max(min_start, as.integer(ceiling(quantile(n_rep, start_q))))
  if (start >= max(n_rep)) {
    message("  auto cutoff: no samples above the normal range -> no placeholders; using ", fallback)
    return(fallback)
  }
  empty <- which(counts[(start + 1L):(max(n_rep) + 1L)] == 0L)
  if (!length(empty)) {
    warning("No clear gap found in the repeated-value distribution; using PLACEHOLDER_AUTO_FALLBACK = ", fallback,
            ". Check the printed distribution and set PLACEHOLDER_MIN_COLS by hand if needed.")
    return(fallback)
  }
  gap <- start + empty[1] - 1L                                  # first count with no samples
  message(sprintf("  auto cutoff: normal range up to ~%d repeated values; first empty count = %d -> cutoff %d",
                  start, gap, gap + 1L))
  gap + 1L
}

flag_placeholder_samples <- function(dat, feature_cols, id_col, min_repeat, min_cols) {
  rep_mat <- matrix(FALSE, nrow(dat), length(feature_cols))
  for (j in seq_along(feature_cols)) {
    x <- dat[[feature_cols[j]]]
    common <- data.table(v = x)[!is.na(v), .N, by = v][N >= min_repeat, v]
    if (length(common)) rep_mat[, j] <- !is.na(x) & x %in% common
  }
  n_rep <- rowSums(rep_mat)
  if (identical(min_cols, "auto")) min_cols <- choose_placeholder_cutoff(n_rep)
  message("  number of metabolites holding a repeated value, per sample:")
  print(table(cut(n_rep, unique(c(-1, 0, 2, 4, 6, 8, 10, min_cols - 1, min_cols, 15, 20, 40, 100, length(feature_cols))))))
  rows <- which(n_rep >= min_cols)
  message(sprintf("  cutoff: >= %d repeated metabolites -> %d placeholder sample(s)", min_cols, length(rows)))
  flagged <- rbindlist(lapply(rows, function(i) {
    cols <- feature_cols[rep_mat[i, ]]
    data.table(sample_id = dat[[id_col]][i], metabolite = cols,
               original_value = unlist(dat[i, ..cols], use.names = FALSE),
               reason = sprintf("placeholder sample: repeated values in %d metabolites", n_rep[i]))
  }))
  for (i in rows) for (cn in feature_cols[rep_mat[i, ]]) set(dat, i = i, j = cn, value = NA_real_)
  if (!nrow(flagged)) flagged <- data.table(sample_id = character(0), metabolite = character(0),
                                            original_value = numeric(0), reason = character(0))
  list(dat = dat, flagged = flagged, placeholder_ids = dat[[id_col]][rows], n_rep = n_rep)
}

knn_impute_matrix <- function(mat, k, seed) {
  if (!anyNA(mat)) { message("  no missing values; imputation skipped"); return(mat) }
  fit <- impute.knn(t(mat), k = k, rowmax = 0.9, colmax = 0.9, rng.seed = seed)   # impute.knn expects features x samples
  out <- t(fit$data)
  dimnames(out) <- dimnames(mat)
  out
}

# =============================================================================
# 1. Read data and check the input assumption
# =============================================================================
if (!file.exists(INPUT_FILE)) stop("Input file not found: ", INPUT_FILE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_DIR)) stop("Output folder could not be created: ", OUT_DIR, " -- set OUT_DIR to a folder you can write to")
out_path <- function(name) file.path(OUT_DIR, paste0(OUT_PREFIX, "_", name))
OUT_SEP <- if (OUT_FORMAT == "tsv") "\t" else ","
write_out <- function(x, name) fwrite(x, out_path(paste0(name, ".", OUT_FORMAT)), sep = OUT_SEP, na = "NA")

message("Reading ", INPUT_FILE)
dat <- fread(INPUT_FILE)
missing_ids <- setdiff(ID_COLS, names(dat))
if (length(missing_ids)) stop("ID column(s) not found: ", paste(missing_ids, collapse = ", "))
if (anyDuplicated(dat[[ID_COLS[1]]])) stop("Sample ID column '", ID_COLS[1], "' has duplicated values")

feature_cols <- setdiff(names(dat), ID_COLS)
non_numeric  <- feature_cols[!vapply(dat[, ..feature_cols], is.numeric, logical(1))]
if (length(non_numeric)) stop("Non-numeric metabolite column(s): ", paste(head(non_numeric, 10), collapse = ", "),
                              " -- add them to ID_COLS or remove them")
message("  ", nrow(dat), " samples x ", length(feature_cols), " metabolites")
n_samples_start  <- nrow(dat)
n_features_start <- length(feature_cols)
n_missing_input  <- sum(is.na(dat[, ..feature_cols]))
status <- copy(dat[, ..ID_COLS])                     # one row per input sample, for the status file
check_standardized(dat, feature_cols, STANDARDIZED_TOL)

# =============================================================================
# 2. Placeholder samples
# =============================================================================
message("Detecting placeholder samples ...")
ph <- flag_placeholder_samples(dat, feature_cols, ID_COLS[1], PLACEHOLDER_MIN_REPEAT, PLACEHOLDER_MIN_COLS)
dat <- ph$dat
status[, n_repeated_metabolites := ph$n_rep]
status[, placeholder_sample := as.integer(get(ID_COLS[1]) %in% ph$placeholder_ids)]
write_out(ph$flagged, "flagged_values")
message("  placeholder samples: ", length(ph$placeholder_ids), " | cells set to NA: ", nrow(ph$flagged))

# =============================================================================
# 3. Missingness QC
# =============================================================================
qc_log <- list()

feat_miss <- vapply(dat[, ..feature_cols], function(x) mean(is.na(x)), numeric(1))
drop_feat <- names(feat_miss)[feat_miss > FEATURE_MISS_MAX]
if (length(drop_feat)) {
  qc_log$features_missing <- data.table(item = drop_feat, item_type = "metabolite",
                                        reason = sprintf("missing in > %.0f%% of samples", 100 * FEATURE_MISS_MAX))
  dat <- dat[, setdiff(names(dat), drop_feat), with = FALSE]
  feature_cols <- setdiff(feature_cols, drop_feat)
}
message("Metabolites dropped for missingness: ", length(drop_feat))

row_miss <- rowSums(is.na(dat[, ..feature_cols])) / length(feature_cols)
drop_rows <- which(row_miss > SAMPLE_MISS_MAX)
if (length(drop_rows)) {
  drop_ids <- dat[[ID_COLS[1]]][drop_rows]
  qc_log$samples_missing <- data.table(item = as.character(drop_ids), item_type = "sample",
                                       reason = ifelse(drop_ids %in% ph$placeholder_ids,
                                                       sprintf("placeholder sample: after removing placeholder values, missing > %.0f%% of metabolites", 100 * SAMPLE_MISS_MAX),
                                                       sprintf("missing > %.0f%% of metabolites", 100 * SAMPLE_MISS_MAX)))
  dat <- dat[-drop_rows]
}
message("Samples dropped for missingness: ", length(drop_rows))

feat_var <- vapply(dat[, ..feature_cols], function(x) var(x, na.rm = TRUE), numeric(1))
drop_const <- names(feat_var)[is.na(feat_var) | feat_var == 0]
if (length(drop_const)) {
  qc_log$features_constant <- data.table(item = drop_const, item_type = "metabolite", reason = "zero variance")
  dat <- dat[, setdiff(names(dat), drop_const), with = FALSE]
  feature_cols <- setdiff(feature_cols, drop_const)
}
message("Constant metabolites dropped: ", length(drop_const))
message("After QC: ", nrow(dat), " samples x ", length(feature_cols), " metabolites")

write_out(if (length(qc_log)) rbindlist(qc_log) else data.table(item = character(0), item_type = character(0), reason = character(0)),
          "qc_dropped_items")

# =============================================================================
# 4. KNN imputation (values stay on the input, standardized scale)
# =============================================================================
block <- as.matrix(dat[, ..feature_cols])
rownames(block) <- dat[[ID_COLS[1]]]
n_imputed   <- rowSums(is.na(block))
pct_imputed <- 100 * n_imputed / length(feature_cols)
message(sprintf("Imputed share per sample: mean %.2f%% | > 30%%: %d samples | > 40%%: %d samples",
                mean(pct_imputed), sum(pct_imputed > 30), sum(pct_imputed > 40)))

message("KNN imputation (k = ", KNN_K, ") ...")
block_imp <- knn_impute_matrix(block, KNN_K, KNN_SEED)
stopifnot(!anyNA(block_imp))

# =============================================================================
# 5. Write outputs
# =============================================================================
final <- cbind(dat[, ..ID_COLS], as.data.table(block_imp))
final[, n_values_imputed   := n_imputed]
final[, pct_values_imputed := pct_imputed]
final[, placeholder_sample := as.integer(get(ID_COLS[1]) %in% ph$placeholder_ids)]
setcolorder(final, c(ID_COLS, feature_cols, "n_values_imputed", "pct_values_imputed", "placeholder_sample"))

write_out(final, "imputed_final")
writeLines(feature_cols, out_path("feature_columns.txt"))

kept_ids <- final[[ID_COLS[1]]]
status[final, pct_values_imputed := i.pct_values_imputed, on = ID_COLS[1]]
status[, qc_status := fifelse(!get(ID_COLS[1]) %in% kept_ids,
                              fifelse(placeholder_sample == 1L, "dropped_placeholder", "dropped_missingness"),
                              fifelse(placeholder_sample == 1L, "kept_partial_placeholder", "kept"))]
write_out(status, "sample_qc_status")
cat("\n=== Sample QC status ===\n"); print(status[, .N, by = qc_status])

summary_dt <- data.table(
  metric = c("samples_input", "metabolites_input", "missing_cells_input",
             "placeholder_samples", "cells_set_to_NA_placeholder",
             "metabolites_dropped_missingness", "samples_dropped_missingness", "metabolites_dropped_constant",
             "samples_final", "metabolites_final", "cells_imputed",
             "samples_over_30pct_imputed", "samples_over_40pct_imputed"),
  value  = c(n_samples_start, n_features_start, n_missing_input,
             length(ph$placeholder_ids), nrow(ph$flagged),
             length(drop_feat), length(drop_rows), length(drop_const),
             nrow(final), length(feature_cols), sum(n_imputed),
             sum(pct_imputed > 30), sum(pct_imputed > 40)))
write_out(summary_dt, "qc_summary")

cat("\n=== QC summary ===\n"); print(summary_dt)
message("Outputs written to ", OUT_DIR)