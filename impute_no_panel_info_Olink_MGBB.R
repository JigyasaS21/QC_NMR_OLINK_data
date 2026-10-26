# ---------------------------------------------------------------------------
# impute_no_panel_info_Olink_MGBB.R
#
# MGBB Olink proteomics (long-format NPX export): QC, duplicate-pair
# resolution, KNN impute, z-score (+ rank-inverse-normal for comparison).
#
# Pipeline:
#   1. Read long NPX export; drop control samples
#   1.2 Detect + resolve "_N"-suffixed duplicate pairs (keep one, flag both)
#   1.5 Assay-name collision check
#   1.6 Flag Olink QC failures (-> NA) and below-LOD (kept, flagged)
#   1.7 Pivot to wide; 1.8 drop proteins >= 99% below LOD
#   2. Missingness QC; 3. KNN impute -> z-score (+ RINT); 4. write
# Requires: data.table; impute (Bioconductor)
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(impute)
})

# =============================================================================
# CONFIG
# =============================================================================
DATA_TYPE <- "olink"

# --- Input -------------------------------------------------------------------
INPUT_TSV    <- "/Path/to/Input/file"
INPUT_FORMAT <- "long"
USE_PRELOADED_OBJECT <- FALSE   
COVARIATE_COLS <- c("SampleID")

LONG_ID_COL      <- "SampleID"
LONG_FEATURE_COL <- "Assay"      # checked for cross-panel collisions (Section 1.5)
LONG_VALUE_COL   <- "NPX"
LONG_VISIT_COL   <- NULL
LONG_VISIT_VALUE <- NULL

# --- Output ------------------------------------------------------------------
OLINK_DIR        <- "/Path/to/OutputDir"
OUT_TSV          <- file.path(OLINK_DIR, "MGBB_olink_imputed_final.csv")
OUT_TSV_RANKNORM <- file.path(OLINK_DIR, "MGBB_olink_imputed_final_ranknorm.csv")
OUT_FEATURES     <- file.path(OLINK_DIR, "MGBB_olink_feature_columns.txt")
QC_SUMMARY_TXT   <- file.path(OLINK_DIR, "MGBB_olink_qc_summary.csv")
QC_DROPPED_TSV   <- file.path(OLINK_DIR, "MGBB_olink_qc_dropped_items.csv")
QC_FLAGGED_TSV   <- file.path(OLINK_DIR, "MGBB_olink_flagged_qc_values.csv")
DUP_CORR_TSV     <- file.path(OLINK_DIR, "MGBB_olink_duplicate_pair_correlations.csv")
DUP_FLAG_TSV     <- file.path(OLINK_DIR, "MGBB_olink_duplicate_samples_flagged.csv")
OUT_SEP          <- ","   # outputs are comma-separated (.csv)

# --- QC / imputation parameters ----------------------------------------------
PROTEIN_MISS_MAX      <- 0.40
SAMPLE_MISS_MAX       <- 0.60
KNN_K                 <- 10
QC_PASS_VALUES        <- c("PASS")   # SampleQC/AssayQC labels counted as passing; WARN/FAIL/NA -> NA.
# Add "WARN" here to keep WARN measurements instead.
NA_BELOW_LOD          <- FALSE       # below-LOD values are flagged but kept
LOD_EXCLUDE_THRESHOLD <- 0.99        # drop proteins below LOD in >= 99% of samples
LOG_METHOD <- "none"                 # NPX is already log2
LOG_OFFSET <- NULL; SKIP_LOG_FEATURES <- character(0); SKIP_LOG_PATTERN <- NULL

# --- Duplicate-pair handling (fully automatic; no hard-coded IDs) -------------
# Magda confirmed "_N"-suffixed pairs are the same specimen run twice. Each clean
# pair is collapsed to ONE row: default keep = fewer QC issues; overridden only if
# the kept-member population-outlier test (BH q < alpha), effect size and
# sequencing depth ALL point the other way. Correlation is a quality flag only.
BRIDGING_ID_SUFFIX_PATTERN      <- "_[0-9]+$"
BRIDGING_PLATE_COL              <- "PlateID"
DUP_LOW_CORR_THRESHOLD          <- 0.80   # flag-only; never decides which member is kept
HANDLE_DUPLICATE_PAIRS          <- TRUE
AUTO_RESOLVE_WITH_TRUST_SIGNALS <- TRUE
OUTLIER_TEST_ALPHA              <- 0.05
OUTLIER_EFFECT_SIZE_MIN         <- 0.60
MIN_SHARED_PROTEINS_FOR_TEST    <- 30

dir.create(OLINK_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Helper functions
# =============================================================================
pivot_long_to_wide <- function(long, id_col, feature_col, value_col, visit_col = NULL, visit_value = NULL) {
  long <- as.data.table(long)
  if (!is.null(visit_col)) { stopifnot(visit_col %in% names(long)); long <- long[get(visit_col) == visit_value] }
  dcast(long, as.formula(paste(id_col, "~", feature_col)), value.var = value_col)
}

flag_olink_qc <- function(long, id_col, feature_col, value_col, na_below_lod = NA_BELOW_LOD,
                          pass_values = QC_PASS_VALUES) {
  long <- as.data.table(long)
  stopifnot(all(c("SampleQC", "AssayQC", "LOD") %in% names(long)))
  fail_idx <- which(is.na(long$SampleQC) | !long$SampleQC %in% pass_values |
                      is.na(long$AssayQC)  | !long$AssayQC  %in% pass_values)
  below_lod_idx <- which(!is.na(long[[value_col]]) & !is.na(long$LOD) & long[[value_col]] < long$LOD)
  flagged <- data.table(sample_id = character(0), column = character(0), original_value = numeric(0), reason = character(0))
  if (length(fail_idx)) {
    flagged <- rbind(flagged, data.table(
      sample_id = long[[id_col]][fail_idx], column = long[[feature_col]][fail_idx],
      original_value = long[[value_col]][fail_idx],
      reason = paste0("failed Olink QC (SampleQC=", long$SampleQC[fail_idx], ", AssayQC=", long$AssayQC[fail_idx], ")")))
    set(long, i = fail_idx, j = value_col, value = NA_real_)
  }
  if (length(below_lod_idx)) {
    flagged <- rbind(flagged, data.table(
      sample_id = long[[id_col]][below_lod_idx], column = long[[feature_col]][below_lod_idx],
      original_value = long[[value_col]][below_lod_idx], reason = "below limit of detection (LOD)"))
    if (na_below_lod) set(long, i = below_lod_idx, j = value_col, value = NA_real_)
  }
  message("  flagged ", length(fail_idx), " row(s) failing Olink QC (labels other than ",
          paste(pass_values, collapse = "/"), " -> NA), ", length(below_lod_idx), " row(s) below LOD",
          if (na_below_lod) " (-> NA)" else " (kept, flagged)")
  list(long = long, flagged = flagged)
}

detect_bridging_samples <- function(long, id_col, plate_col = NULL, suffix_pattern = "_[0-9]+$") {
  long <- as.data.table(long)
  ids <- unique(long[[id_col]])
  suffixed_ids <- ids[grepl(suffix_pattern, ids)]
  grp <- data.table(id = suffixed_ids, base_id = sub(suffix_pattern, "", suffixed_ids))
  multi_base <- grp[, .N, by = base_id][N > 1, base_id]
  candidates <- grp[base_id %in% multi_base]
  if (!nrow(candidates)) return(data.table(id = character(0), base_id = character(0), bare_exists = logical(0),
                                           n_plates = integer(0), auto_exclude = logical(0), needs_review = logical(0)))
  candidates[, bare_exists := base_id %in% ids]
  if (!is.null(plate_col) && plate_col %in% names(long)) {
    plate_lookup <- long[, .(plate = get(plate_col)[1]), by = c(id_col)]
    setnames(plate_lookup, id_col, "id")
    candidates <- merge(candidates, plate_lookup, by = "id", all.x = TRUE)
    candidates <- merge(candidates, candidates[, .(n_plates = uniqueN(plate)), by = base_id], by = "base_id", all.x = TRUE)
  } else candidates[, n_plates := NA_integer_]
  candidates[, auto_exclude := !bare_exists & (is.na(n_plates) | n_plates == 1)]
  candidates[, needs_review := !bare_exists & !is.na(n_plates) & n_plates > 1]
  candidates[]
}

knn_impute_matrix <- function(mat, k = 10) {
  if (!anyNA(mat)) { message("  no missing values, skipping impute.knn"); return(mat) }
  fit <- impute.knn(t(mat), k = k, rowmax = 0.9, colmax = 0.9, rng.seed = 12345)
  out <- t(fit$data); dimnames(out) <- dimnames(mat); out
}
zscore_matrix <- function(mat) scale(mat, center = TRUE, scale = TRUE)
rank_inverse_normal <- function(x, c = 3/8) {          # Blom; as in UKB-PPP
  r <- rank(x, na.last = "keep", ties.method = "average"); n <- sum(!is.na(x))
  qnorm((r - c) / (n - 2 * c + 1))
}
rank_normalize_matrix <- function(mat) { out <- apply(mat, 2, rank_inverse_normal); dimnames(out) <- dimnames(mat); out }
sort_by_suffix_number <- function(ids) ids[order(as.integer(gsub(".*_", "", ids)))]

# =============================================================================
# 1. Read data
# =============================================================================
if (USE_PRELOADED_OBJECT && exists("olink_data")) {
  message("Using already-loaded 'olink_data' object (", nrow(olink_data), " rows)")
  raw <- as.data.table(olink_data)
} else {
  if (!file.exists(INPUT_TSV)) stop("Input file not found: ", INPUT_TSV, " -- update INPUT_TSV in CONFIG")
  message("Reading ", INPUT_TSV, " ...")
  raw <- fread(INPUT_TSV)
}
message("  ", nrow(raw), " rows | sample types: ", paste(unique(raw$SampleType), collapse = ", "))
n_before_type_filter <- nrow(raw)
raw <- raw[SampleType == "SAMPLE"]
message("  dropped ", n_before_type_filter - nrow(raw), " control row(s); ", uniqueN(raw$SampleID), " SAMPLE IDs")

# ---- 1.2 duplicate pairs ------------------------------------------------------
pop_stats <- raw[!is.na(NPX), .(pop_mean = mean(NPX), pop_sd = sd(NPX)), by = Assay]
outlier_z <- function(id) {
  d <- merge(raw[SampleID == id, .(Assay, NPX)], pop_stats, by = "Assay")[pop_sd > 0]
  setNames(abs((d$NPX - d$pop_mean) / d$pop_sd), d$Assay)
}
issue_count <- function(qc, aqc, npx, lod) {
  sum(is.na(qc) | !qc %in% QC_PASS_VALUES | is.na(aqc) | !aqc %in% QC_PASS_VALUES |
        (!is.na(npx) & !is.na(lod) & npx < lod))
}

ids_to_keep <- character(0); ids_to_drop <- character(0); needs_review_ids <- character(0)
kept_corr_lookup <- data.table(SampleID = character(0), duplicate_pair_correlation = numeric(0))

if (!HANDLE_DUPLICATE_PAIRS) {
  message("HANDLE_DUPLICATE_PAIRS = FALSE -- every sample treated independently.")
} else {
  bridging <- detect_bridging_samples(raw, LONG_ID_COL, BRIDGING_PLATE_COL, BRIDGING_ID_SUFFIX_PATTERN)
  dup_candidates   <- bridging[auto_exclude == TRUE]
  needs_review_ids <- bridging[needs_review == TRUE, id]
  n_same  <- uniqueN(dup_candidates$base_id)
  n_multi <- uniqueN(sub(BRIDGING_ID_SUFFIX_PATTERN, "", needs_review_ids))
  if (!n_same && !n_multi) {
    message("No duplicate-pair-like SampleIDs detected -- nothing to collapse.")
  } else {
    message(n_same + n_multi, " candidate duplicate base ID(s): ", n_same, " single-plate (auto-resolved), ",
            n_multi, " multi-plate (left un-collapsed, needs_review)")
  }
  if (length(needs_review_ids)) warning(length(needs_review_ids), " suffixed SampleID(s) span >1 plate -- NOT collapsed, review: ",
                                        paste(unique(sub(BRIDGING_ID_SUFFIX_PATTERN, "", needs_review_ids)), collapse = ", "))
  
  if (nrow(dup_candidates)) {
    group_sizes <- dup_candidates[, .N, by = base_id]
    clean_pairs <- group_sizes[N == 2, base_id]
    irregular   <- group_sizes[N != 2, base_id]
    if (length(irregular)) warning(length(irregular), " duplicate group(s) with != 2 members -- NOT collapsed: ",
                                   paste(irregular, collapse = ", "))
    if (length(clean_pairs)) {
      dup_corr_results <- rbindlist(lapply(clean_pairs, function(b) {
        sibs <- sort_by_suffix_number(dup_candidates[base_id == b, id])
        m <- merge(raw[SampleID == sibs[1], .(Assay, NPX, SampleQC, AssayQC, LOD)],
                   raw[SampleID == sibs[2], .(Assay, NPX, SampleQC, AssayQC, LOD)], by = "Assay", suffixes = c("_a", "_b"))
        data.table(base_id = b, id_a = sibs[1], id_b = sibs[2], n_proteins = nrow(m),
                   correlation = round(suppressWarnings(cor(m$NPX_a, m$NPX_b, use = "complete.obs")), 4),
                   n_issues_a = issue_count(m$SampleQC_a, m$AssayQC_a, m$NPX_a, m$LOD_a),
                   n_issues_b = issue_count(m$SampleQC_b, m$AssayQC_b, m$NPX_b, m$LOD_b))
      }))
      message("Pairs: ", nrow(dup_corr_results), " | median correlation ", round(median(dup_corr_results$correlation), 4))
      low <- dup_corr_results[correlation < DUP_LOW_CORR_THRESHOLD, base_id]
      if (length(low)) message(length(low), " pair(s) correlate below ", DUP_LOW_CORR_THRESHOLD,
                               " (flag only, still collapsed): ", paste(low, collapse = ", "))
      
      dup_corr_results[, `:=`(default_keep = fifelse(n_issues_a <= n_issues_b, "a", "b"),
                              outlier_test_p = NA_real_, outlier_effect = NA_real_,
                              depth_a = NA_real_, depth_b = NA_real_, depth_agrees = NA)]
      dup_corr_results[, signal_keep := default_keep]
      if (AUTO_RESOLVE_WITH_TRUST_SIGNALS) {
        for (i in seq_len(nrow(dup_corr_results))) {
          row <- dup_corr_results[i]
          za <- outlier_z(row$id_a); zb <- outlier_z(row$id_b); shared <- intersect(names(za), names(zb))
          if (length(shared) >= MIN_SHARED_PROTEINS_FOR_TEST) {
            wt <- suppressWarnings(wilcox.test(za[shared], zb[shared], paired = TRUE))
            diffs <- za[shared] - zb[shared]
            keep_sig <- if (mean(diffs) > 0) "b" else "a"
            eff <- if (keep_sig == "b") mean(diffs > 0) else 1 - mean(diffs > 0)
            d_a <- median(raw[SampleID == row$id_a, Count], na.rm = TRUE)
            d_b <- median(raw[SampleID == row$id_b, Count], na.rm = TRUE)
            set(dup_corr_results, i, "outlier_test_p", wt$p.value)
            set(dup_corr_results, i, "outlier_effect", eff)
            set(dup_corr_results, i, "signal_keep", keep_sig)
            set(dup_corr_results, i, "depth_a", d_a); set(dup_corr_results, i, "depth_b", d_b)
            set(dup_corr_results, i, "depth_agrees", (keep_sig == "a" && d_a >= d_b) || (keep_sig == "b" && d_b >= d_a))
          }
        }
      }
      dup_corr_results[, outlier_test_q := p.adjust(outlier_test_p, method = "BH")]
      dup_corr_results[, use_outlier_override := !is.na(outlier_test_q) & outlier_test_q < OUTLIER_TEST_ALPHA &
                         !is.na(outlier_effect) & outlier_effect >= OUTLIER_EFFECT_SIZE_MIN &
                         !is.na(depth_agrees) & depth_agrees & signal_keep != default_keep]
      fwrite(dup_corr_results, DUP_CORR_TSV, sep = OUT_SEP)
      if (any(dup_corr_results$use_outlier_override))
        message(sum(dup_corr_results$use_outlier_override), " automatic override(s): ",
                paste(dup_corr_results[use_outlier_override == TRUE, base_id], collapse = ", "))
      
      for (i in seq_len(nrow(dup_corr_results))) {
        row <- dup_corr_results[i]
        k <- if (row$use_outlier_override) row$signal_keep else row$default_keep
        keep_id <- if (k == "a") row$id_a else row$id_b
        drop_id <- if (k == "a") row$id_b else row$id_a
        ids_to_keep <- c(ids_to_keep, keep_id); ids_to_drop <- c(ids_to_drop, drop_id)
        kept_corr_lookup <- rbind(kept_corr_lookup, data.table(
          SampleID = sub(BRIDGING_ID_SUFFIX_PATTERN, "", keep_id), duplicate_pair_correlation = row$correlation))
      }
      fwrite(merge(dup_candidates[, .(id, base_id, kept = id %in% ids_to_keep)],
                   dup_corr_results[, .(base_id, correlation)], by = "base_id", all.x = TRUE), DUP_FLAG_TSV, sep = OUT_SEP)
    }
  }
}
n_before_dup <- nrow(raw)
raw <- raw[!SampleID %in% ids_to_drop]
raw[SampleID %in% ids_to_keep, SampleID := sub(BRIDGING_ID_SUFFIX_PATTERN, "", SampleID)]
message("  collapsed ", length(ids_to_keep), " pair(s); dropped ", n_before_dup - nrow(raw), " sibling row(s)")

# ---- 1.5 assay-name collisions ----------------------------------------------
dup_assays <- unique(raw[, .(Assay, Panel)])[, .N, by = Assay][N > 1, Assay]
if (length(dup_assays)) stop("Assay name(s) on >1 Panel: ", paste(head(dup_assays, 10), collapse = ", "),
                             ". Switch LONG_FEATURE_COL to 'OlinkID'.")

# ---- 1.6 Olink QC flags --------------------------------------------------------
message("Flagging rows failing Olink QC / below LOD ...")
flag_result <- flag_olink_qc(raw, LONG_ID_COL, LONG_FEATURE_COL, LONG_VALUE_COL)
raw <- flag_result$long; flagged_values <- flag_result$flagged
fwrite(flagged_values, QC_FLAGGED_TSV, sep = OUT_SEP)

# ---- 1.7 pivot (guard: one row per SampleID x Assay) ---------------------------
n_dup_cells <- raw[, .N, by = .(SampleID, Assay)][N > 1, .N]
if (n_dup_cells) stop(n_dup_cells, " SampleID x Assay combination(s) appear more than once -- pivot would aggregate.",
                      " Check for un-collapsed duplicates (needs_review / irregular groups).")
dat <- pivot_long_to_wide(raw, LONG_ID_COL, LONG_FEATURE_COL, LONG_VALUE_COL, LONG_VISIT_COL, LONG_VISIT_VALUE)
rm(raw); invisible(gc())
feature_cols <- setdiff(names(dat), COVARIATE_COLS)
message("  ", nrow(dat), " samples x ", length(feature_cols), " proteins")
n_features_start <- length(feature_cols); n_samples_start <- nrow(dat)
qc_log <- list()

# ---- 1.8 essentially undetectable proteins -----------------------------------
lod_n <- flagged_values[reason == "below limit of detection (LOD)", .N, by = column]
lod_n[, pct_below_lod := N / nrow(dat)]
drop_lod <- lod_n[pct_below_lod >= LOD_EXCLUDE_THRESHOLD, column]
message("Proteins >= ", LOD_EXCLUDE_THRESHOLD * 100, "% below LOD (excluded): ", length(drop_lod),
        if (length(drop_lod)) paste0(" -- ", paste(drop_lod, collapse = ", ")) else "")
if (length(drop_lod)) {
  qc_log[["below_lod_exclusion"]] <- data.table(item = drop_lod, item_type = "feature",
                                                reason = sprintf("%.0f%%+ of samples below LOD", LOD_EXCLUDE_THRESHOLD * 100))
  dat <- dat[, setdiff(names(dat), drop_lod), with = FALSE]; feature_cols <- setdiff(feature_cols, drop_lod)
}

# =============================================================================
# 2. Missingness QC
# =============================================================================
feat_miss <- vapply(dat[, ..feature_cols], function(x) mean(is.na(x)), numeric(1))
drop_feat <- names(feat_miss)[feat_miss > PROTEIN_MISS_MAX]
message("Features with >", PROTEIN_MISS_MAX * 100, "% missingness: ", length(drop_feat),
        if (length(drop_feat)) paste0(" -- ", paste(drop_feat, collapse = ", ")) else "")
if (length(drop_feat)) {
  qc_log[["feat_miss"]] <- data.table(item = drop_feat, item_type = "feature",
                                      reason = sprintf("missingness > %.0f%%", PROTEIN_MISS_MAX * 100))
  dat <- dat[, setdiff(names(dat), drop_feat), with = FALSE]; feature_cols <- setdiff(feature_cols, drop_feat)
}
row_miss <- rowSums(is.na(dat[, ..feature_cols])) / length(feature_cols)
drop_row_idx <- which(row_miss > SAMPLE_MISS_MAX)
message("Samples with >", SAMPLE_MISS_MAX * 100, "% missingness: ", length(drop_row_idx))
if (length(drop_row_idx)) {
  qc_log[["sample_miss"]] <- data.table(item = as.character(dat[[COVARIATE_COLS[1]]][drop_row_idx]),
                                        item_type = "sample", reason = sprintf("missingness > %.0f%%", SAMPLE_MISS_MAX * 100))
  dat <- dat[-drop_row_idx]
}
feat_var <- vapply(dat[, ..feature_cols], function(x) var(x, na.rm = TRUE), numeric(1))
drop_zerovar <- names(feat_var)[is.na(feat_var) | feat_var == 0]
if (length(drop_zerovar)) {
  qc_log[["zero_var"]] <- data.table(item = drop_zerovar, item_type = "feature", reason = "zero variance")
  dat <- dat[, setdiff(names(dat), drop_zerovar), with = FALSE]; feature_cols <- setdiff(feature_cols, drop_zerovar)
}
message("After QC: ", nrow(dat), " samples x ", length(feature_cols), " features")

fwrite(data.table(
  metric = c("features_start", "features_dropped_below_lod", "features_dropped_missingness", "features_dropped_zero_variance",
             "features_final", "samples_start", "samples_dropped_missingness", "samples_final",
             "rows_flagged_failed_qc_or_below_lod", "duplicate_pairs_collapsed_to_one_row",
             "duplicate_rows_dropped_as_sibling", "duplicate_ids_flagged_for_manual_review"),
  value  = c(n_features_start, length(drop_lod), length(drop_feat), length(drop_zerovar), length(feature_cols),
             n_samples_start, length(drop_row_idx), nrow(dat), nrow(flagged_values),
             length(ids_to_keep), length(ids_to_drop), length(needs_review_ids))), QC_SUMMARY_TXT, sep = OUT_SEP)
if (length(qc_log)) fwrite(rbindlist(qc_log), QC_DROPPED_TSV, sep = OUT_SEP)

# =============================================================================
# 3. KNN impute -> z-score (+ RINT)
# =============================================================================
block <- as.matrix(dat[, ..feature_cols]); rownames(block) <- dat[[COVARIATE_COLS[1]]]
n_missing_before_impute   <- rowSums(is.na(block))
pct_missing_before_impute <- 100 * n_missing_before_impute / length(feature_cols)
block_imp  <- knn_impute_matrix(block, k = KNN_K)   # NPX already log2 -> no log step
stopifnot(!anyNA(block_imp))
block_z    <- zscore_matrix(block_imp)
block_rank <- rank_normalize_matrix(block_imp)

# =============================================================================
# 4. Write outputs
# =============================================================================
write_final <- function(mat, path) {
  out <- cbind(dat[, ..COVARIATE_COLS], as.data.table(mat))
  out[, n_values_imputed := n_missing_before_impute]
  out[, pct_values_imputed := pct_missing_before_impute]
  out[, duplicate_pair_correlation := NA_real_]
  if (nrow(kept_corr_lookup)) out[kept_corr_lookup, duplicate_pair_correlation := i.duplicate_pair_correlation, on = "SampleID"]
  fwrite(out, path, sep = OUT_SEP, na = "NA")
  message("wrote ", path, " (", nrow(out), " x ", ncol(out), ")")
}
write_final(block_z, OUT_TSV)
write_final(block_rank, OUT_TSV_RANKNORM)
writeLines(colnames(block_z), OUT_FEATURES)

cat("\n=== Final sample count reconciliation ===\n")
cat("SAMPLE IDs after control filter:     ", n_samples_start + length(ids_to_drop), "\n")
cat("  minus duplicate siblings dropped:  -", length(ids_to_drop), "\n", sep = "")
cat("  minus missingness-QC samples:      -", length(drop_row_idx), "\n", sep = "")
cat("  = final:                            ", nrow(dat), "\n")