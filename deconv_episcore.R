#!/usr/bin/env Rscript
#
# deconv_episcore.R
#
# EpiSCORE deconvolution of urine EPIC v2 data (n=8) using mrefBladder.m
# (most tissue-appropriate reference for urine samples), plus IC-specific
# CpG overlap analysis using the centEpiFibIC EpiDish reference.
#
# Outputs:
#   episcore_bladder_fractions.csv
#   episcore_summary.txt
#   ic_cpg_overlap.txt
#   episcore_vs_epidish.png
#
# Usage:
#   Rscript deconv_episcore.R <betas_full.csv> <output_dir>

suppressPackageStartupMessages(library(ggplot2))

# Shared project theme + palettes.
.script_dir <- (function() {
  args <- commandArgs(trailingOnly = FALSE)
  fn <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(fn)) return(dirname(normalizePath(fn)))
  of <- try(sys.frame(1)$ofile, silent = TRUE)
  if (!inherits(of, "try-error") && !is.null(of))
    return(dirname(normalizePath(of)))
  "."
})()
source(file.path(.script_dir, "_style.R"))

args       <- commandArgs(trailingOnly = TRUE)
betas_path <- if (length(args) >= 1) args[1] else
  "menopause_urine_out/betas_full.csv"
out_dir    <- if (length(args) >= 2) args[2] else
  "menopause_urine_out"
dmp_path   <- file.path(path.expand(out_dir), "dmp_list.csv")

betas_path <- path.expand(betas_path)
out_dir    <- path.expand(out_dir)

# ── Section subfolders ───────────────────────────────────────────────────────
qc_dir      <- file.path(out_dir, "01_qc")
deconv_dir  <- file.path(out_dir, "02_deconvolution")
dm_dir      <- file.path(out_dir, "03_diff_meth")
fig_dir     <- file.path(out_dir, "04_figures")
pub_fig_dir <- file.path(out_dir, "05_pub_figures")
for (d in c(qc_dir, deconv_dir, dm_dir, fig_dir, pub_fig_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("=== EpiSCORE + IC CpG Overlap: Menopause Urine Pilot ===\n\n")

# ── Phenotype ──────────────────────────────────────────────────────────────
pheno <- data.frame(
  sample_id = c("207758450138_R03C01","207758450138_R04C01",
                "207758450138_R07C01","207758450138_R08C01",
                "208356980025_R02C01","208356980025_R03C01",
                "208356980025_R05C01","208356980025_R06C01"),
  subject   = c("FMLS_007","FMLS_009","FMLS_005","FMLS_006",
                "FMLS_011","FMLS_015","FMLS_002","FMLS_010"),
  status    = c("Post","Pre","Post","Post","Pre","Pre","Pre","Pre"),
  age       = c(61, 31, 62, 66, 45, 23, 27, 41),
  flag      = c("","OC_pill","","ovary_removed","PREGNANT","OC_pill","",""),
  stringsAsFactors = FALSE
)

# ── Load betas ─────────────────────────────────────────────────────────────
cat("Loading beta matrix...\n")
beta_df <- read.csv(betas_path, row.names = 1, check.names = FALSE)
beta_df <- beta_df[, pheno$sample_id, drop = FALSE]
colnames(beta_df) <- pheno$subject
beta_mat <- as.matrix(beta_df)
cat(sprintf("  %d probes x %d samples loaded\n\n", nrow(beta_mat), ncol(beta_mat)))

# EPIC v2 suffix stripping: all probes are cg00XXXXX_TCXX etc.
array_probes   <- rownames(beta_mat)
array_base_ids <- sub("_[^_]+$", "", array_probes)
base_to_full   <- setNames(array_probes, array_base_ids)
cat(sprintf("  Probe ID examples (raw):    %s\n", paste(head(array_probes, 3), collapse=", ")))
cat(sprintf("  Probe ID examples (base):   %s\n\n", paste(head(array_base_ids, 3), collapse=", ")))

# ══════════════════════════════════════════════════════════════════════════
# 1. EpiSCORE — Bladder reference (mrefBladder.m)
# ══════════════════════════════════════════════════════════════════════════
cat("=== 1. EpiSCORE Deconvolution (Bladder reference) ===\n\n")

if (!requireNamespace("EpiSCORE", quietly = TRUE)) {
  cat("  Installing EpiSCORE from GitHub...\n")
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes", repos = "https://cloud.r-project.org")
  remotes::install_github("aet21/EpiSCORE")
}
if (!requireNamespace("EpiDISH", quietly = TRUE)) {
  BiocManager::install("EpiDISH", ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(EpiSCORE)
  library(EpiDISH)
})
cat("  EpiSCORE and EpiDISH loaded\n\n")

# ── Step 1a: Map EPIC v2 probes → Entrez gene IDs via probeInfo450k ────────
cat("  Step 1a: Mapping EPIC v2 probes to genes (TSS-proximal CpGs)...\n")
pi450   <- probeInfo450k.lv      # built-in EpiSCORE data object
# probeID, EID (Entrez), GeneGroup: 1=TSS200, 2=TSS1500, 3=body, etc.
idx     <- match(array_base_ids, pi450$probeID)
matched <- !is.na(idx)

our_eid <- rep(NA_integer_, length(array_probes))
our_gg  <- rep(NA_integer_, length(array_probes))
our_eid[matched] <- pi450$EID[idx[matched]]
our_gg[matched]  <- pi450$GeneGroup[idx[matched]]

n_matched  <- sum(matched)
n_promoter <- sum(!is.na(our_eid) & our_gg %in% c(1, 2))
n_genes    <- length(unique(our_eid[!is.na(our_eid) & our_gg %in% c(1, 2)]))

cat(sprintf("  Probes matched to 450K probeInfo:  %d / %d (%.1f%%)\n",
            n_matched, length(array_probes), 100*n_matched/length(array_probes)))
cat(sprintf("  Promoter CpGs (TSS200 + TSS1500):  %d\n", n_promoter))
cat(sprintf("  Unique Entrez genes covered:        %d\n\n", n_genes))

# ── Step 1b: Aggregate betas per gene (mean over promoter CpGs) ────────────
cat("  Step 1b: Aggregating betas per gene...\n")
pmask          <- !is.na(our_eid) & our_gg %in% c(1, 2)
promoter_idx   <- which(pmask)
promoter_eids  <- our_eid[pmask]
unique_genes   <- sort(unique(promoter_eids))

gene_beta <- matrix(NA_real_, nrow = length(unique_genes), ncol = ncol(beta_mat),
                    dimnames = list(as.character(unique_genes), colnames(beta_mat)))

for (g in unique_genes) {
  gp <- promoter_idx[promoter_eids == g]
  gene_beta[as.character(g), ] <- if (length(gp) == 1) beta_mat[gp, ] else
    colMeans(beta_mat[gp, , drop = FALSE], na.rm = TRUE)
}
gene_beta[is.nan(gene_beta)] <- NA

# Impute remaining NAs with row median
na_genes <- rowSums(is.na(gene_beta))
if (any(na_genes > 0)) {
  cat(sprintf("  Imputing NAs in %d genes (row median)\n", sum(na_genes > 0)))
  gene_meds <- apply(gene_beta, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(gene_beta))) {
    nas <- is.na(gene_beta[i, ]); if (any(nas)) gene_beta[i, nas] <- gene_meds[i]
  }
}
cat(sprintf("  Gene-level beta matrix: %d genes x %d samples\n\n",
            nrow(gene_beta), ncol(gene_beta)))

# ── Step 1c: Load bladder reference and run wRPC ───────────────────────────
cat("  Step 1c: Running wRPC with mrefBladder.m...\n")
data(mrefBladder.m)
ref_genes   <- intersect(rownames(mrefBladder.m), rownames(gene_beta))
ref_missing <- setdiff(rownames(mrefBladder.m), rownames(gene_beta))

cat(sprintf("  Bladder reference: %d genes total\n", nrow(mrefBladder.m)))
cat(sprintf("  Reference genes in our data: %d / %d (%.1f%%)\n",
            length(ref_genes), nrow(mrefBladder.m),
            100*length(ref_genes)/nrow(mrefBladder.m)))
cat(sprintf("  Reference genes missing from our data: %d\n\n", length(ref_missing)))
cat(sprintf("  Cell types in reference: %s\n\n",
            paste(colnames(mrefBladder.m), collapse = ", ")))

episcore_res <- tryCatch(
  wRPC(data   = gene_beta[ref_genes, ],
       ref.m  = mrefBladder.m[ref_genes, ],
       useW   = TRUE, wth = 0.4, maxit = 200),
  error = function(e) {
    cat("  ERROR in wRPC:", e$message, "\n")
    NULL
  }
)

if (is.null(episcore_res) || is.null(episcore_res$estF)) {
  cat("  EpiSCORE returned no results.\n")
  frac_es <- NULL
} else {
  frac_es <- as.data.frame(episcore_res$estF)
  frac_es$subject <- rownames(frac_es)
  frac_es <- merge(frac_es, pheno[, c("subject","status","age","flag")], by = "subject")

  cell_cols <- setdiff(colnames(frac_es), c("subject","status","age","flag"))

  # Sort by status then age for display
  frac_ord <- frac_es[order(frac_es$status, frac_es$age), ]

  cat("\n--- EpiSCORE Bladder fractions (per sample) ---\n")
  print(format(frac_ord[, c("subject","status","age","flag", cell_cols)],
               digits = 3, nsmall = 3), row.names = FALSE)

  cat("\n--- Mean fractions by menopausal status ---\n")
  for (grp in c("Pre","Post")) {
    sub  <- frac_es[frac_es$status == grp, cell_cols, drop = FALSE]
    vals <- colMeans(sub, na.rm = TRUE)
    sds  <- apply(sub, 2, sd, na.rm = TRUE)
    cat(sprintf("\n  %s-menopausal (n=%d):\n", grp, nrow(sub)))
    for (ct in cell_cols) {
      cat(sprintf("    %-20s  mean=%.4f  SD=%.4f\n", ct, vals[ct], sds[ct]))
    }
  }

  # Statistical tests on each cell type
  cat("\n--- Wilcoxon tests (Post > Pre) per cell type ---\n")
  for (ct in cell_cols) {
    pre  <- frac_es[frac_es$status == "Pre",  ct]
    post <- frac_es[frac_es$status == "Post", ct]
    wt   <- tryCatch(wilcox.test(post, pre, alternative = "greater", exact = FALSE),
                     error = function(e) NULL)
    if (!is.null(wt)) {
      cat(sprintf("  %-20s  W=%.0f  p=%.4f\n", ct, wt$statistic, wt$p.value))
    }
  }

  # Check for IC-like column
  ic_col <- grep("(?i)(immun|leuko|IC|stromal)", cell_cols, value = TRUE, perl = TRUE)
  if (length(ic_col) == 0) ic_col <- cell_cols[1]  # fallback: first column
  cat(sprintf("\n  Note: treating '%s' as primary immune/stromal proxy for downstream use\n\n",
              paste(ic_col, collapse = "+")))

  write.csv(frac_es, file.path(deconv_dir, "episcore_bladder_fractions.csv"), row.names = FALSE)
  cat("  Saved: episcore_bladder_fractions.csv\n\n")

  # ── Write detailed summary ─────────────────────────────────────────────
  sink_file <- file.path(deconv_dir, "episcore_summary.txt")
  sink(sink_file)
  cat("=== EpiSCORE Deconvolution: Menopause Urine Pilot ===\n\n")
  cat(sprintf("Reference: mrefBladder.m (%d cell types)\n", length(cell_cols)))
  cat(sprintf("Reference cell types: %s\n\n", paste(colnames(mrefBladder.m), collapse = ", ")))
  cat(sprintf("Reference genes used: %d / %d (%.1f%%)\n\n",
              length(ref_genes), nrow(mrefBladder.m),
              100*length(ref_genes)/nrow(mrefBladder.m)))

  cat("Per-sample fractions:\n")
  print(frac_ord[, c("subject","status","age","flag", cell_cols)], row.names = FALSE, digits = 4)

  cat("\nMean by status:\n")
  for (grp in c("Pre","Post")) {
    sub  <- frac_es[frac_es$status == grp, cell_cols, drop = FALSE]
    vals <- colMeans(sub, na.rm = TRUE)
    sds  <- apply(sub, 2, sd, na.rm = TRUE)
    cat(sprintf("\n%s-menopausal (n=%d):\n", grp, nrow(sub)))
    for (ct in cell_cols) {
      cat(sprintf("  %-20s  mean=%.4f  SD=%.4f\n", ct, vals[ct], sds[ct]))
    }
  }

  cat("\nWilcoxon tests (Post > Pre):\n")
  for (ct in cell_cols) {
    pre  <- frac_es[frac_es$status == "Pre",  ct]
    post <- frac_es[frac_es$status == "Post", ct]
    wt   <- tryCatch(wilcox.test(post, pre, alternative = "greater", exact = FALSE),
                     error = function(e) NULL)
    if (!is.null(wt))
      cat(sprintf("  %-20s  W=%.0f  p=%.4f\n", ct, wt$statistic, wt$p.value))
  }
  sink()
  cat(sprintf("  Saved: episcore_summary.txt\n\n"))
}

# ══════════════════════════════════════════════════════════════════════════
# 2. IC-specific CpG overlap with DMPs
# ══════════════════════════════════════════════════════════════════════════
cat("=== 2. IC-specific CpG Overlap with Pre vs Post DMPs ===\n\n")

if (!file.exists(dmp_path)) {
  cat("  dmp_list.csv not found — skipping. Run age_cpg_enrichment.R first.\n\n")
} else {
  tt <- read.csv(dmp_path, row.names = 1)
  # Sort by p-value ascending (limma output should already be sorted but ensure)
  tt <- tt[order(tt$P.Value), ]
  cat(sprintf("  Loaded %d probes from limma (Pre vs Post)\n", nrow(tt)))
  cat(sprintf("  Top DMP: %s  logFC=%.3f  p=%.2e\n\n",
              rownames(tt)[1], tt$logFC[1], tt$P.Value[1]))

  # ── Load centEpiFibIC reference ──────────────────────────────────────────
  suppressPackageStartupMessages(data(centEpiFibIC.m))
  cat(sprintf("  centEpiFibIC reference: %d probes x %d cell types (%s)\n\n",
              nrow(centEpiFibIC.m), ncol(centEpiFibIC.m),
              paste(colnames(centEpiFibIC.m), collapse = ", ")))

  # IC specificity score: IC beta minus best non-IC beta
  ic_score  <- centEpiFibIC.m[, "IC"] - pmax(centEpiFibIC.m[, "Epi"],
                                               centEpiFibIC.m[, "Fib"])
  epi_score <- centEpiFibIC.m[, "Epi"] - pmax(centEpiFibIC.m[, "IC"],
                                                centEpiFibIC.m[, "Fib"])

  ic_specific   <- names(ic_score[ic_score   >  0.3])
  ic_moderate   <- names(ic_score[ic_score   >  0.1 & ic_score <= 0.3])
  epi_specific  <- names(epi_score[epi_score >  0.3])   # negative control
  epi_moderate  <- names(epi_score[epi_score >  0.1 & epi_score <= 0.3])

  cat(sprintf("  IC-specific  probes (IC score >  0.3):  %d\n", length(ic_specific)))
  cat(sprintf("  IC-moderate  probes (IC score 0.1-0.3): %d\n", length(ic_moderate)))
  cat(sprintf("  Epi-specific probes (Epi score > 0.3):  %d\n", length(epi_specific)))
  cat(sprintf("  Epi-moderate probes (Epi score 0.1-0.3):%d\n\n", length(epi_moderate)))

  # Map reference probe IDs to our EPIC v2 DMP universe (suffix-aware)
  dmp_base     <- sub("_[^_]+$", "", rownames(tt))
  all_dmp_ids  <- rownames(tt)

  find_dmp_probes <- function(ref_ids) {
    all_dmp_ids[dmp_base %in% ref_ids]
  }

  ic_in_array   <- find_dmp_probes(ic_specific)
  icm_in_array  <- find_dmp_probes(ic_moderate)
  epi_in_array  <- find_dmp_probes(epi_specific)
  epim_in_array <- find_dmp_probes(epi_moderate)

  cat("  Probe overlap with our EPIC v2 array:\n")
  cat(sprintf("    IC-specific  → %d / %d (%.1f%%) in array\n",
              length(ic_in_array),   length(ic_specific),   100*length(ic_in_array)/max(1,length(ic_specific))))
  cat(sprintf("    IC-moderate  → %d / %d (%.1f%%) in array\n",
              length(icm_in_array),  length(ic_moderate),   100*length(icm_in_array)/max(1,length(ic_moderate))))
  cat(sprintf("    Epi-specific → %d / %d (%.1f%%) in array\n",
              length(epi_in_array),  length(epi_specific),  100*length(epi_in_array)/max(1,length(epi_specific))))
  cat(sprintf("    Epi-moderate → %d / %d (%.1f%%) in array\n\n",
              length(epim_in_array), length(epi_moderate),  100*length(epim_in_array)/max(1,length(epi_moderate))))

  # ── Fisher enrichment tests ──────────────────────────────────────────────
  run_fisher <- function(sig_probes, target_probes, universe, label) {
    n11 <- sum(sig_probes %in% target_probes)
    n10 <- length(sig_probes) - n11
    n01 <- length(target_probes) - n11
    n00 <- length(universe) - n11 - n10 - n01
    if (n00 < 0) n00 <- 0
    mat <- matrix(c(n11, n10, n01, n00), 2)
    ft  <- fisher.test(mat, alternative = "greater")
    list(label = label, n_sig = length(sig_probes), n_target = length(target_probes),
         n_overlap = n11, obs_pct = 100*n11/max(1,length(sig_probes)),
         bg_pct = 100*length(target_probes)/length(universe),
         OR = as.numeric(ft$estimate), p = ft$p.value)
  }

  print_fisher_row <- function(r) {
    cat(sprintf("  %-38s n_sig=%6d  overlap=%4d  obs=%.2f%%  bg=%.2f%%  OR=%6.2f  p=%.4g\n",
                r$label, r$n_sig, r$n_overlap, r$obs_pct, r$bg_pct, r$OR, r$p))
  }

  p_thresholds <- c(0.05, 0.01, 0.001)
  top_ns       <- c(500, 1000, 5000, 10000)

  cat("--- Fisher enrichment: IC-specific probes in DMPs ---\n")
  cat(sprintf("  %-38s %-14s %-13s %-13s %-13s %-9s %s\n",
              "Threshold", "n_DMPs", "overlap", "obs%", "bg%", "OR", "p"))
  for (thr in p_thresholds) {
    sig <- rownames(tt)[tt$P.Value < thr]
    print_fisher_row(run_fisher(sig, ic_in_array, all_dmp_ids,
                                sprintf("IC-specific (p<%.3f)", thr)))
  }
  for (n in top_ns) {
    sig <- rownames(tt)[seq_len(min(n, nrow(tt)))]
    print_fisher_row(run_fisher(sig, ic_in_array, all_dmp_ids,
                                sprintf("IC-specific (top %d)", n)))
  }

  cat("\n--- Fisher enrichment: IC-moderate probes in DMPs ---\n")
  for (thr in p_thresholds) {
    sig <- rownames(tt)[tt$P.Value < thr]
    print_fisher_row(run_fisher(sig, icm_in_array, all_dmp_ids,
                                sprintf("IC-moderate (p<%.3f)", thr)))
  }
  for (n in top_ns) {
    sig <- rownames(tt)[seq_len(min(n, nrow(tt)))]
    print_fisher_row(run_fisher(sig, icm_in_array, all_dmp_ids,
                                sprintf("IC-moderate (top %d)", n)))
  }

  cat("\n--- Fisher enrichment: Epi-specific probes (negative control) ---\n")
  for (thr in p_thresholds) {
    sig <- rownames(tt)[tt$P.Value < thr]
    print_fisher_row(run_fisher(sig, epi_in_array, all_dmp_ids,
                                sprintf("Epi-specific (p<%.3f)", thr)))
  }
  for (n in top_ns) {
    sig <- rownames(tt)[seq_len(min(n, nrow(tt)))]
    print_fisher_row(run_fisher(sig, epi_in_array, all_dmp_ids,
                                sprintf("Epi-specific (top %d)", n)))
  }

  # ── Direction of IC-specific probes among top DMPs ───────────────────────
  cat("\n--- Direction check: IC-specific probes among top DMPs ---\n")
  cat("  (positive logFC = higher methylation in Post vs Pre)\n\n")
  for (n in c(500, 1000, 5000)) {
    topn     <- rownames(tt)[seq_len(min(n, nrow(tt)))]
    topn_ic  <- topn[topn %in% ic_in_array]
    if (length(topn_ic) > 0) {
      lfc      <- tt[topn_ic, "logFC"]
      n_pos    <- sum(lfc > 0)
      n_neg    <- sum(lfc < 0)
      med_lfc  <- median(lfc)
      bt       <- tryCatch(binom.test(n_pos, length(topn_ic), 0.5, alternative="greater"),
                           error = function(e) NULL)
      bt_p     <- if (!is.null(bt)) bt$p.value else NA
      cat(sprintf("  Top %5d DMPs → %d IC-specific:  logFC>0=%d  logFC<0=%d  median=%.3f  binomial_p(>0)=%.4g\n",
                  n, length(topn_ic), n_pos, n_neg, med_lfc, bt_p))
    } else {
      cat(sprintf("  Top %5d DMPs → 0 IC-specific probes\n", n))
    }
  }

  # ── Detailed look at top IC-specific DMPs ───────────────────────────────
  cat("\n--- Top IC-specific probes by DMP rank ---\n")
  ic_dmp_rows <- which(all_dmp_ids %in% ic_in_array)
  if (length(ic_dmp_rows) > 0) {
    top_ic <- tt[ic_dmp_rows[seq_len(min(20, length(ic_dmp_rows)))],
                 c("logFC","AveExpr","t","P.Value","adj.P.Val"), drop=FALSE]
    top_ic$IC_score <- ic_score[dmp_base[ic_dmp_rows[seq_len(min(20, length(ic_dmp_rows)))]]]
    top_ic$rank     <- ic_dmp_rows[seq_len(min(20, length(ic_dmp_rows)))]
    cat(sprintf("  %-20s %7s %8s %7s %10s %10s %8s %6s\n",
                "probe", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "IC_score", "rank"))
    for (i in seq_len(nrow(top_ic))) {
      cat(sprintf("  %-20s %7.3f %8.3f %7.3f %10.2e %10.2e %8.3f %6d\n",
                  rownames(top_ic)[i], top_ic$logFC[i], top_ic$AveExpr[i],
                  top_ic$t[i], top_ic$P.Value[i], top_ic$adj.P.Val[i],
                  top_ic$IC_score[i], top_ic$rank[i]))
    }
  } else {
    cat("  No IC-specific probes found in DMP list.\n")
  }

  # ── Save full overlap summary ────────────────────────────────────────────
  all_results <- list()
  for (thr in p_thresholds) {
    sig <- rownames(tt)[tt$P.Value < thr]
    for (probe_set in list(list(ic_in_array, "IC_specific"), list(epi_in_array, "Epi_specific"))) {
      r <- run_fisher(sig, probe_set[[1]], all_dmp_ids,
                      sprintf("%s_p%.3f", probe_set[[2]], thr))
      all_results[[length(all_results)+1]] <- data.frame(
        probe_set = probe_set[[2]], threshold = sprintf("p<%s",thr),
        n_sig = r$n_sig, n_target = r$n_target, n_overlap = r$n_overlap,
        obs_pct = r$obs_pct, bg_pct = r$bg_pct, OR = r$OR, p = r$p,
        stringsAsFactors = FALSE)
    }
  }
  overlap_df <- do.call(rbind, all_results)
  write.csv(overlap_df, file.path(dm_dir, "ic_cpg_overlap.csv"), row.names = FALSE)
  cat(sprintf("\n  Saved: ic_cpg_overlap.csv\n\n"))
}

# ══════════════════════════════════════════════════════════════════════════
# 3. Comparison plot: EpiSCORE vs EpiDish cell fractions
# ══════════════════════════════════════════════════════════════════════════
cat("=== 3. EpiSCORE vs EpiDish comparison plot ===\n")

epidish_path <- file.path(deconv_dir, "epidish_fractions.csv")
if (!is.null(frac_es) && file.exists(epidish_path)) {
  frac_ed <- read.csv(epidish_path)

  # EpiDish IC column
  ed_plot <- data.frame(
    subject = frac_ed$subject,
    status  = frac_ed$status,
    value   = frac_ed$IC,
    method  = "EpiDish\n(centEpiFibIC)",
    cell    = "IC (Immune)",
    stringsAsFactors = FALSE
  )

  # EpiSCORE: identify all cell columns
  es_cell_cols <- setdiff(colnames(frac_es), c("subject","status","age","flag"))
  es_long <- do.call(rbind, lapply(es_cell_cols, function(ct) {
    data.frame(subject = frac_es$subject,
               status  = frac_es$status,
               value   = frac_es[[ct]],
               method  = "EpiSCORE\n(mrefBladder)",
               cell    = ct,
               stringsAsFactors = FALSE)
  }))

  # Combined plot: EpiDish IC vs each EpiSCORE cell type
  p1 <- ggplot(ed_plot, aes(x = status, y = value, colour = status, shape = status)) +
    geom_jitter(width = 0.1, size = 2.6) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.35,
                 colour = "grey20", linewidth = 0.4) +
    scale_colour_manual(values = STATUS_COLORS) +
    short_labs(title = "EpiDish IC",
               x = NULL, y = "Estimated fraction") +
    project_theme(base_size = 10) + theme(legend.position = "none")

  p2 <- ggplot(es_long, aes(x = status, y = value, colour = status, shape = status)) +
    geom_jitter(width = 0.1, size = 2.6) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.35,
                 colour = "grey20", linewidth = 0.4) +
    scale_colour_manual(values = STATUS_COLORS) +
    facet_wrap(~ cell, scales = "free_y") +
    short_labs(title = "EpiSCORE bladder",
               x = NULL, y = "Estimated fraction") +
    project_theme(base_size = 10) + theme(legend.position = "none")

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    g <- gridExtra::arrangeGrob(p1, p2, ncol = 2, widths = c(1, 2))
    ggsave(file.path(deconv_dir, "episcore_vs_epidish.png"), g,
           width = 7.2, height = 3.8, dpi = 300, bg = "white")
  } else {
    ggsave(file.path(deconv_dir, "episcore_bladder_fractions.png"), p2,
           width = 6.5, height = 3.8, dpi = 300, bg = "white")
  }
  cat("  Saved: episcore_vs_epidish.png\n")
} else if (!is.null(frac_es)) {
  es_cell_cols <- setdiff(colnames(frac_es), c("subject","status","age","flag"))
  es_long <- do.call(rbind, lapply(es_cell_cols, function(ct) {
    data.frame(subject = frac_es$subject, status = frac_es$status,
               value = frac_es[[ct]], cell = ct, stringsAsFactors = FALSE)
  }))
  p <- ggplot(es_long, aes(x = status, y = value, colour = status)) +
    geom_jitter(width = 0.1, size = 2.6) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.35,
                 colour = "grey20", linewidth = 0.4) +
    scale_colour_manual(values = STATUS_COLORS) +
    facet_wrap(~ cell, scales = "free_y") +
    short_labs(x = NULL, y = "Fraction") +
    project_theme(base_size = 10) + theme(legend.position = "none")
  ggsave(file.path(deconv_dir, "episcore_bladder_fractions.png"), p,
         width = 6.5, height = 3.8, dpi = 300, bg = "white")
  cat("  Saved: episcore_bladder_fractions.png\n")
} else {
  cat("  No EpiSCORE results to plot.\n")
}

cat("\n=== Done ===\n")
