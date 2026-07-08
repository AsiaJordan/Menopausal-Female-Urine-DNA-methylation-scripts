#!/usr/bin/env Rscript
#
# preprocess_menopause_urine.R
#
# SeSAMe QCDPB preprocessing for the menopause urine EPIC v2 pilot (n=8).
# Saves:
#   betas_full.csv      — all EPIC v2 probes × samples (for EpiDish / limma)
#   qc_report.csv       — per-sample QC metrics
#   sentrix_info.csv    — chip/row/col from filenames (batch structure check)
#
# Usage:
#   Rscript preprocess_menopause_urine.R \
#       --idat_dir  /path/to/idats \
#       --output_dir /path/to/output
#
# Requires: sesame, sesameData, BiocParallel

cat("=== Menopause Urine EPIC v2: SeSAMe Preprocessing ===\n\n")

# ── Args ──────────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
idat_dir   <- NULL
output_dir <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--idat_dir"   && i < length(args)) { idat_dir   <- args[i+1]; i <- i+2
  } else if (args[i] == "--output_dir" && i < length(args)) { output_dir <- args[i+1]; i <- i+2
  } else { i <- i+1 }
}

if (is.null(idat_dir) || !dir.exists(idat_dir)) {
  cat("Usage: Rscript preprocess_menopause_urine.R --idat_dir /path/to/idats\n")
  quit(status = 1)
}
if (is.null(output_dir)) output_dir <- file.path(idat_dir, "sesame_out")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ── Packages ───────────────────────────────────────────────────────────────
for (pkg in c("sesame", "sesameData", "BiocParallel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
suppressPackageStartupMessages({
  library(sesame)
  library(BiocParallel)
})

sesameDataCache()   # downloads manifests on first run; cached thereafter

# ── Discover IDATs (non-recursive: only the root of idat_dir) ─────────────
grn_files <- list.files(idat_dir, pattern = "_Grn\\.idat(\\.gz)?$",
                        recursive = FALSE, full.names = TRUE)

if (length(grn_files) == 0) {
  cat("ERROR: No *_Grn.idat files found directly in", idat_dir, "\n")
  quit(status = 1)
}

idat_prefixes <- sub("_Grn\\.idat(\\.gz)?$", "", grn_files)
sample_names  <- basename(idat_prefixes)
n_samples     <- length(idat_prefixes)
cat(sprintf("Found %d samples:\n", n_samples))
for (s in sample_names) cat("  ", s, "\n")
cat("\n")

# ── Parse Sentrix IDs for batch structure check ────────────────────────────
# Format: {SentrixID}_{RowCol}  e.g. 207758450138_R03C01
sentrix_df <- data.frame(
  sample     = sample_names,
  sentrix_id = sub("_R\\d+C\\d+$", "", sample_names),
  row_col    = regmatches(sample_names, regexpr("R\\d+C\\d+$", sample_names)),
  stringsAsFactors = FALSE
)
cat("Sentrix chip structure (batch check):\n")
print(sentrix_df)
cat("\n")

# ── Run SeSAMe QCDPB prep ─────────────────────────────────────────────────
cat("Running openSesame (QCDPB prep)...\n")
cat("  Prep: pOOBAH → NOOB → dye-bias correction → BMIQ\n\n")

n_cpus <- min(n_samples, max(1L, parallel::detectCores() - 1L))
BPPARAM <- MulticoreParam(n_cpus)

# Process each sample, get SigDF (not just betas) for QC
sdfs <- bplapply(idat_prefixes, function(pref) {
  tryCatch(
    openSesame(pref, prep = "QCDPB", func = NULL),
    error = function(e) { message("FAILED: ", basename(pref), " — ", e$message); NULL }
  )
}, BPPARAM = BPPARAM)
names(sdfs) <- sample_names

# ── QC metrics ─────────────────────────────────────────────────────────────
cat("\nCalculating QC metrics...\n")
qc_rows <- lapply(seq_along(sdfs), function(si) {
  sname <- sample_names[si]
  sdf   <- sdfs[[si]]

  if (is.null(sdf)) {
    return(data.frame(sample=sname, qc_flag="FAIL_PREPROCESSING",
                      frac_dt=NA, mean_intensity=NA, mean_beta=NA,
                      frac_na=NA, RGratio=NA, RGdistort=NA, stringsAsFactors=FALSE))
  }

  qc_det  <- tryCatch(as.data.frame(sesameQC_calcStats(sdf, "detection")),  error=function(e) NULL)
  qc_int  <- tryCatch(as.data.frame(sesameQC_calcStats(sdf, "intensity")),  error=function(e) NULL)
  qc_bet  <- tryCatch(as.data.frame(sesameQC_calcStats(sdf, "betas")),      error=function(e) NULL)
  qc_dye  <- tryCatch(as.data.frame(sesameQC_calcStats(sdf, "dyeBias")),    error=function(e) NULL)

  frac_dt   <- if (!is.null(qc_det)) qc_det$frac_dt          else NA
  mean_int  <- if (!is.null(qc_int)) qc_int$mean_intensity    else NA
  frac_na   <- if (!is.null(qc_bet)) qc_bet$frac_na_cg        else NA

  flag <- "PASS"
  if (!is.na(frac_dt)  && frac_dt  < 0.80) flag <- "WARN_LOW_DETECTION"
  if (!is.na(mean_int) && mean_int < 500)   flag <- "WARN_LOW_INTENSITY"
  if (!is.na(frac_na)  && frac_na  > 0.10)  flag <- "WARN_HIGH_MISSINGNESS"

  cat(sprintf("  %s  frac_dt=%.3f  mean_int=%.0f  frac_na=%.4f  [%s]\n",
              sname,
              if (is.na(frac_dt)) -1 else frac_dt,
              if (is.na(mean_int)) -1 else mean_int,
              if (is.na(frac_na)) -1 else frac_na,
              flag))

  data.frame(
    sample        = sname,
    sentrix_id    = sentrix_df$sentrix_id[si],
    row_col       = sentrix_df$row_col[si],
    qc_flag       = flag,
    frac_dt       = round(if (is.na(frac_dt))  NA else frac_dt,  4),
    mean_intensity= round(if (is.na(mean_int)) NA else mean_int, 1),
    mean_beta     = round(if (!is.null(qc_bet)) qc_bet$mean_beta_cg  else NA, 4),
    median_beta   = round(if (!is.null(qc_bet)) qc_bet$median_beta_cg else NA, 4),
    frac_na       = round(if (is.na(frac_na))  NA else frac_na,  4),
    RGratio       = round(if (!is.null(qc_dye)) qc_dye$RGratio    else NA, 4),
    RGdistort     = round(if (!is.null(qc_dye)) qc_dye$RGdistort  else NA, 4),
    mean_oob_grn  = round(if (!is.null(qc_int)) qc_int$mean_oob_grn else NA, 1),
    mean_oob_red  = round(if (!is.null(qc_int)) qc_int$mean_oob_red else NA, 1),
    stringsAsFactors = FALSE
  )
})
qc_df <- do.call(rbind, qc_rows)

n_pass <- sum(qc_df$qc_flag == "PASS")
n_warn <- sum(grepl("^WARN", qc_df$qc_flag))
n_fail <- sum(grepl("^FAIL", qc_df$qc_flag))
cat(sprintf("\nQC summary: %d PASS, %d WARN, %d FAIL\n\n", n_pass, n_warn, n_fail))

# ── Extract full beta matrix ───────────────────────────────────────────────
cat("Extracting beta values (all EPIC v2 probes)...\n")
beta_list <- lapply(sdfs, function(sdf) {
  if (is.null(sdf)) return(NULL)
  getBetas(sdf)
})

# Find common probes across all passing samples
passing <- !sapply(beta_list, is.null)
all_probes <- Reduce(intersect, lapply(beta_list[passing], names))
cat(sprintf("  Common probes across passing samples: %d\n", length(all_probes)))

# Build matrix: probes × samples
beta_mat <- matrix(NA_real_, nrow = length(all_probes), ncol = n_samples,
                   dimnames = list(all_probes, sample_names))
for (si in seq_along(beta_list)) {
  if (!is.null(beta_list[[si]])) {
    b <- beta_list[[si]]
    common <- intersect(names(b), all_probes)
    beta_mat[common, si] <- b[common]
  }
}

cat(sprintf("  Beta matrix: %d probes × %d samples\n", nrow(beta_mat), ncol(beta_mat)))
cat(sprintf("  Overall missingness: %.2f%%\n\n",
            100 * mean(is.na(beta_mat))))

# ── Save outputs ───────────────────────────────────────────────────────────
# betas_full.csv stays at the project root since every downstream script
# reads it as the master input. QC tables go to 01_qc/.
qc_dir <- file.path(output_dir, "01_qc")
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

cat("Saving outputs to", output_dir, "...\n")

beta_df <- as.data.frame(beta_mat)
beta_df <- cbind(probe = rownames(beta_df), beta_df)
write.csv(beta_df, file.path(output_dir, "betas_full.csv"), row.names = FALSE)
cat("  betas_full.csv\n")

write.csv(qc_df, file.path(qc_dir, "qc_report.csv"), row.names = FALSE)
cat("  01_qc/qc_report.csv\n")

write.csv(sentrix_df, file.path(qc_dir, "sentrix_info.csv"), row.names = FALSE)
cat("  01_qc/sentrix_info.csv\n")

cat("\n=== Done ===\n")
cat(sprintf("Beta matrix saved: %d probes × %d samples\n", nrow(beta_mat), ncol(beta_mat)))
cat("Next steps:\n")
cat("  1. Check qc_report.csv — flag any WARN/FAIL samples\n")
cat("  2. Check sentrix_info.csv — samples split across 2 chips (batch effect risk)\n")
cat("  3. Run EpiDish deconvolution on betas_full.csv\n")
cat("  4. Run limma differential methylation\n")
cat("  5. Run exhaustive permutation test (C(8,3)=56 permutations)\n")
