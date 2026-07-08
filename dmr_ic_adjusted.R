#!/usr/bin/env Rscript
#
# dmr_ic_adjusted.R
#
# DMR analysis of urine EPIC v2 data (n=8, Pre vs Post) using DMRcate,
# with and without IC (immune cell) fraction adjustment.
#
# Models compared:
#   Unadjusted:  ~ status
#   IC-adjusted: ~ IC + status
#
# Outputs:
#   dmr_unadjusted.csv
#   dmr_ic_adjusted.csv
#   dmr_summary.txt
#   dmr_top_regions.png
#
# Usage:
#   Rscript dmr_ic_adjusted.R <betas_full.csv> <output_dir>

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

cat("=== DMR Analysis (unadjusted vs IC-adjusted) ===\n\n")

# ── Install / load packages ───────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
for (pkg in c("limma", "DMRcate")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
suppressPackageStartupMessages({
  library(limma)
  library(DMRcate)
})
cat(sprintf("  DMRcate version: %s\n\n", as.character(packageVersion("DMRcate"))))

# ── Phenotype ─────────────────────────────────────────────────────────────
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
pheno$status <- factor(pheno$status, levels = c("Pre","Post"))

# ── Load IC fractions ─────────────────────────────────────────────────────
ic_path <- file.path(deconv_dir, "epidish_fractions.csv")
if (!file.exists(ic_path))
  stop("epidish_fractions.csv not found — run epidish_menopause_urine.R first")
ic_df <- read.csv(ic_path)
pheno <- merge(pheno, ic_df[, c("subject","IC")], by = "subject")
# Restore original sample order
pheno <- pheno[match(c("FMLS_007","FMLS_009","FMLS_005","FMLS_006",
                        "FMLS_011","FMLS_015","FMLS_002","FMLS_010"),
                     pheno$subject), ]
cat("Phenotype + IC fractions:\n")
print(pheno[, c("subject","status","age","flag","IC")], row.names = FALSE)
cat("\n")

# ── Load betas ────────────────────────────────────────────────────────────
cat("Loading beta matrix...\n")
beta_df  <- read.csv(betas_path, row.names = 1, check.names = FALSE)
beta_df  <- beta_df[, pheno$sample_id, drop = FALSE]
colnames(beta_df) <- pheno$subject
beta_mat <- as.matrix(beta_df)

# Remove probes with any NA (DMRcate dislikes them)
na_rows  <- rowSums(is.na(beta_mat)) > 0
beta_mat <- beta_mat[!na_rows, ]
cat(sprintf("  %d probes x %d samples (%d probes with NAs removed)\n\n",
            nrow(beta_mat), ncol(beta_mat), sum(na_rows)))

# ── Suffix-stripped beta matrix ───────────────────────────────────────────
# EPIC v2 probes: cg00XXXXX_TC21 etc. DMRcate EPIC/450K needs plain cg IDs.
# For each base ID keep the probe with lowest variance (most stable signal).
cat("Building suffix-stripped beta matrix (deduplication by lowest row variance)...\n")
base_ids <- sub("_[^_]+$", "", rownames(beta_mat))
dup_base <- duplicated(base_ids)
cat(sprintf("  %d raw probes → %d unique base IDs (%d duplicates removed)\n\n",
            nrow(beta_mat), sum(!dup_base), sum(dup_base)))

# Among duplicated base IDs, keep the one with lowest row variance (more stable)
row_vars  <- apply(beta_mat, 1, var, na.rm = TRUE)
ord       <- order(row_vars)          # ascending: lowest variance first
beta_sorted  <- beta_mat[ord, ]
base_sorted  <- base_ids[ord]
keep         <- !duplicated(base_sorted)
beta_stripped <- beta_sorted[keep, ]
rownames(beta_stripped) <- base_sorted[keep]
cat(sprintf("  Stripped matrix: %d probes x %d samples\n\n",
            nrow(beta_stripped), ncol(beta_stripped)))

# ── M-value matrices ─────────────────────────────────────────────────────
to_M <- function(b) {
  M <- log2(pmax(b, 1e-6) / pmax(1 - b, 1e-6))
  M[!is.finite(M)] <- NA
  M
}
M_stripped <- to_M(beta_stripped)

# ── Limma fit helper ──────────────────────────────────────────────────────
# Returns list(fit2 = MArrayLM, tt = topTable df).
run_limma <- function(M_mat, design, coef_name, label = "") {
  cat(sprintf("  [%s] Running limma on %d probes...\n", label, nrow(M_mat)))
  fit  <- lmFit(M_mat, design)
  fit2 <- eBayes(fit, robust = TRUE)
  tt   <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "none")
  cat(sprintf("  [%s] p<0.05: %d (%.1f%%)  p<0.01: %d  p<0.001: %d  FDR<0.05: %d  FDR<0.20: %d\n",
              label,
              sum(tt$P.Value < 0.05),   100*mean(tt$P.Value < 0.05),
              sum(tt$P.Value < 0.01),   sum(tt$P.Value < 0.001),
              sum(tt$adj.P.Val < 0.05), sum(tt$adj.P.Val < 0.20)))
  list(fit2 = fit2, tt = tt)
}

# ── cpg.annotate wrapper ──────────────────────────────────────────────────
# DMRcate ≥ 3.x takes the M-value matrix as `object` and fits the model
# internally from design + coef.
annotate_fit <- function(M_mat, design, coef_name, label = "") {
  for (atype in c("EPICv2", "EPICv1", "EPIC", "450K")) {
    res <- tryCatch({
      cat(sprintf("  [%s] cpg.annotate arraytype='%s'...\n", label, atype))
      cpg.annotate(
        datatype      = "array",
        object        = M_mat,
        what          = "M",
        arraytype     = atype,
        analysis.type = "differential",
        design        = design,
        coef          = coef_name,
        fdr           = 0.05
      )
    }, error = function(e) {
      cat(sprintf("  [%s] %s failed: %s\n", label, atype, conditionMessage(e)))
      NULL
    })
    if (!is.null(res)) {
      cat(sprintf("  [%s] Annotation succeeded with arraytype='%s'\n", label, atype))
      return(res)
    }
  }
  cat(sprintf("  [%s] All arraytypes failed.\n", label))
  NULL
}

# ══════════════════════════════════════════════════════════════════════════
# Designs
# ══════════════════════════════════════════════════════════════════════════
design_unadj <- model.matrix(~ status,    data = pheno)
design_adj   <- model.matrix(~ IC + status, data = pheno)
coef <- "statusPost"

cat("=== Model 1: Unadjusted (~ status) ===\n")
res_u  <- run_limma(M_stripped, design_unadj, coef, "Unadjusted")
anno_u <- annotate_fit(M_stripped, design_unadj, coef, "Unadjusted")

cat("\n=== Model 2: IC-adjusted (~ IC + status) ===\n")
res_a  <- run_limma(M_stripped, design_adj, coef, "IC-adjusted")
anno_a <- annotate_fit(M_stripped, design_adj, coef, "IC-adjusted")

# ══════════════════════════════════════════════════════════════════════════
# DMRcate
# ══════════════════════════════════════════════════════════════════════════
run_dmrcate <- function(anno, label) {
  if (is.null(anno)) {
    cat(sprintf("  [%s] Skipping — annotation unavailable.\n", label))
    return(NULL)
  }
  cat(sprintf("  [%s] Running dmrcate() (lambda=1000, C=2, min.cpgs=3)...\n", label))
  tryCatch({
    dmrcate(anno, lambda = 1000, C = 2, min.cpgs = 3, pcutoff = 0.05)
  }, error = function(e) {
    cat(sprintf("  [%s] dmrcate() failed: %s\n", label, conditionMessage(e)))
    NULL
  })
}

cat("\n=== DMRcate region calling ===\n")
dmr_u <- run_dmrcate(anno_u, "Unadjusted")
dmr_a <- run_dmrcate(anno_a, "IC-adjusted")

# ── Extract ranges ────────────────────────────────────────────────────────
# DMRcate ≥ 3.x renamed the effect-size columns from meanbetafc/maxbetafc
# to meandiff/maxdiff. Rename back so the rest of this script (and any
# downstream consumers) keeps working with either version.
normalise_dmr_cols <- function(df) {
  if ("meandiff" %in% colnames(df) && !"meanbetafc" %in% colnames(df))
    colnames(df)[colnames(df) == "meandiff"] <- "meanbetafc"
  if ("maxdiff"  %in% colnames(df) && !"maxbetafc"  %in% colnames(df))
    colnames(df)[colnames(df) == "maxdiff"]  <- "maxbetafc"
  df
}

extract_dmrs <- function(dmr_obj, label) {
  if (is.null(dmr_obj)) return(NULL)
  for (genome in c("hg38", NULL)) {
    df <- tryCatch(
      as.data.frame(if (is.null(genome)) extractRanges(dmr_obj)
                    else extractRanges(dmr_obj, genome = genome)),
      error = function(e) NULL
    )
    if (!is.null(df)) {
      df <- normalise_dmr_cols(df)
      df <- df[order(abs(df$meanbetafc), decreasing = TRUE), ]
      df$model <- label

      # Report total candidates, HMFDR-significant, and manuscript-filter
      # (|Δβ|≥0.10 AND Stouffer<0.05) counts side by side.
      n_total <- nrow(df)
      n_hmfdr <- if ("HMFDR" %in% colnames(df))
                   sum(df$HMFDR < 0.05, na.rm = TRUE) else NA_integer_
      n_ms    <- if (all(c("Stouffer","meanbetafc") %in% colnames(df)))
                   sum(df$Stouffer < 0.05 & abs(df$meanbetafc) >= 0.10,
                       na.rm = TRUE) else NA_integer_
      cat(sprintf("  [%s] %d candidates, %d HMFDR<0.05, %d at manuscript filter (|Δβ|≥0.10 AND Stouffer<0.05)\n",
                  label, n_total, n_hmfdr, n_ms))
      return(df)
    }
  }
  cat(sprintf("  [%s] extractRanges failed.\n", label))
  NULL
}

df_u <- extract_dmrs(dmr_u, "unadjusted")
df_a <- extract_dmrs(dmr_a, "ic_adjusted")

# ══════════════════════════════════════════════════════════════════════════
# Print results
# ══════════════════════════════════════════════════════════════════════════
print_dmrs <- function(df, label, n = 40) {
  if (is.null(df) || nrow(df) == 0) {
    cat(sprintf("\n=== %s DMRs (0 at manuscript filter) ===\n  None.\n", label))
    return(invisible())
  }
  n_total <- nrow(df)
  n_ms    <- if (all(c("Stouffer","meanbetafc") %in% colnames(df)))
               sum(df$Stouffer < 0.05 & abs(df$meanbetafc) >= 0.10,
                   na.rm = TRUE) else n_total
  cat(sprintf("\n=== %s DMRs (%d at manuscript filter of %d candidates) ===\n",
              label, n_ms, n_total))

  sig <- if (all(c("Stouffer","meanbetafc") %in% colnames(df)))
           df[!is.na(df$Stouffer) & df$Stouffer < 0.05 &
              abs(df$meanbetafc) >= 0.10, ] else df
  if (nrow(sig) == 0) {
    cat("  No DMRs meet |Δβ|≥0.10 AND Stouffer<0.05. Top candidates by |Δβ|:\n")
    sig <- df
  }

  cols <- intersect(c("seqnames","start","end","width","no.cpgs",
                       "meanbetafc","maxbetafc","min_smoothed_fdr",
                       "Stouffer","HMFDR","Fisher","overlapping.genes"),
                    colnames(sig))
  top <- head(sig[, cols, drop = FALSE], n)
  for (col in colnames(top))
    if (is.numeric(top[[col]])) top[[col]] <- signif(top[[col]], 3)
  print(top, row.names = FALSE)

  cat("\n  Summary (manuscript-filter regions only):\n")
  cat(sprintf("    Regions: %d (of %d candidates)\n", nrow(sig), n_total))
  if ("no.cpgs" %in% colnames(sig))
    cat(sprintf("    CpGs/region: median=%d  range=%d–%d\n",
                as.integer(median(sig$no.cpgs)), min(sig$no.cpgs), max(sig$no.cpgs)))
  if ("meanbetafc" %in% colnames(sig)) {
    cat(sprintf("    Post > Pre (meanbetafc > 0): %d\n", sum(sig$meanbetafc > 0, na.rm=TRUE)))
    cat(sprintf("    Pre > Post (meanbetafc < 0): %d\n", sum(sig$meanbetafc < 0, na.rm=TRUE)))
    cat(sprintf("    |mean beta FC|: %.3f – %.3f\n",
                min(abs(sig$meanbetafc),na.rm=TRUE), max(abs(sig$meanbetafc),na.rm=TRUE)))
  }
  if ("overlapping.genes" %in% colnames(sig)) {
    genes <- unique(trimws(unlist(strsplit(sig$overlapping.genes, ";|,"))))
    genes <- genes[nchar(genes) > 0]
    cat(sprintf("    Unique overlapping genes: %d\n", length(genes)))
    if (length(genes) <= 40) cat(sprintf("    %s\n", paste(sort(genes), collapse=", ")))
  }
}

print_dmrs(df_u, "Unadjusted")
print_dmrs(df_a, "IC-adjusted")

# Cross-model comparison using the manuscript's DMR filter so counts are
# directly comparable to paragraph 88 of the manuscript ("6,072 DMRs").
ms_sig <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df[0, ])
  if (!all(c("Stouffer","meanbetafc") %in% colnames(df))) return(df)
  df[!is.na(df$Stouffer) & df$Stouffer < 0.05 &
     abs(df$meanbetafc) >= 0.10, ]
}
sig_u <- ms_sig(df_u)
sig_a <- ms_sig(df_a)

if (!is.null(df_u) && !is.null(df_a)) {
  cat("\n=== Cross-model comparison (|Δβ|≥0.10 AND Stouffer<0.05) ===\n")
  cat(sprintf("  Unadjusted:  %d DMRs (of %d candidates)\n",
              nrow(sig_u), nrow(df_u)))
  cat(sprintf("  IC-adjusted: %d DMRs (of %d candidates)\n",
              nrow(sig_a), nrow(df_a)))
  if ("meanbetafc" %in% colnames(sig_u))
    cat(sprintf("  Unadjusted:  %d Hyper / %d Hypo (%.1f%% Hyper)\n",
                sum(sig_u$meanbetafc > 0), sum(sig_u$meanbetafc < 0),
                100 * mean(sig_u$meanbetafc > 0)))
  if ("meanbetafc" %in% colnames(sig_a))
    cat(sprintf("  IC-adjusted: %d Hyper / %d Hypo (%.1f%% Hyper)\n",
                sum(sig_a$meanbetafc > 0), sum(sig_a$meanbetafc < 0),
                100 * mean(sig_a$meanbetafc > 0)))

  parse_genes <- function(df) {
    if (!"overlapping.genes" %in% colnames(df)) return(character(0))
    unique(trimws(unlist(strsplit(df$overlapping.genes, ";|,"))))
  }
  gu <- parse_genes(sig_u)[nchar(parse_genes(sig_u)) > 0]
  ga <- parse_genes(sig_a)[nchar(parse_genes(sig_a)) > 0]
  cat(sprintf("\n  Genes shared between models:        %d\n",
              length(intersect(gu, ga))))
  cat(sprintf("  Genes unadjusted only (lost on IC):  %d\n",
              length(setdiff(gu, ga))))
  cat(sprintf("  Genes IC-adjusted only (newly seen): %d\n",
              length(setdiff(ga, gu))))
}

# ── Save outputs ──────────────────────────────────────────────────────────
if (!is.null(df_u)) {
  write.csv(df_u, file.path(dm_dir, "dmr_unadjusted.csv"), row.names = FALSE)
  cat("\nSaved: dmr_unadjusted.csv\n")
}
if (!is.null(df_a)) {
  write.csv(df_a, file.path(dm_dir, "dmr_ic_adjusted.csv"), row.names = FALSE)
  cat("Saved: dmr_ic_adjusted.csv\n")
}

# ── Save summary to file ──────────────────────────────────────────────────
sink(file.path(dm_dir, "dmr_summary.txt"))
cat("=== DMR Analysis: Menopause Urine Pilot (EPIC v2, n=8) ===\n\n")
cat("Models:\n")
cat("  1. Unadjusted:  ~ status\n")
cat("  2. IC-adjusted: ~ IC + status\n\n")
cat("Probe-level limma (on suffix-stripped, deduplicated matrix):\n")
for (nm in c("Unadjusted","IC-adjusted")) {
  r <- if (nm == "Unadjusted") res_u$tt else res_a$tt
  cat(sprintf("  %-12s p<0.05: %d (%.1f%%)  p<0.01: %d  p<0.001: %d  FDR<0.05: %d  FDR<0.20: %d\n",
              nm,
              sum(r$P.Value < 0.05), 100*mean(r$P.Value < 0.05),
              sum(r$P.Value < 0.01), sum(r$P.Value < 0.001),
              sum(r$adj.P.Val < 0.05), sum(r$adj.P.Val < 0.20)))
}
cat("\n")
print_dmrs(df_u, "Unadjusted")
print_dmrs(df_a, "IC-adjusted")
sink()
cat("Saved: dmr_summary.txt\n")

# ── Lollipop plot: top IC-adjusted DMRs ──────────────────────────────────
# Restrict to manuscript-filter DMRs (|Δβ|≥0.10 AND Stouffer<0.05) so the
# plot shows regions comparable to the manuscript's published top hits.
# Fall back to all candidates if nothing survives.
plot_src <- if (!is.null(df_a) &&
                all(c("Stouffer","meanbetafc") %in% colnames(df_a)))
              df_a[!is.na(df_a$Stouffer) & df_a$Stouffer < 0.05 &
                   abs(df_a$meanbetafc) >= 0.10, ] else df_a
if (is.null(plot_src) || nrow(plot_src) == 0) plot_src <- df_a

if (!is.null(plot_src) && nrow(plot_src) > 0 &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  top_n <- min(25, nrow(plot_src))
  top   <- head(plot_src[order(abs(plot_src$meanbetafc), decreasing = TRUE), ],
                top_n)
  coord <- paste0(top$seqnames, ":",
                  formatC(top$start, format = "d", big.mark = ","))
  if ("overlapping.genes" %in% colnames(top)) {
    top$label <- ifelse(nchar(trimws(top$overlapping.genes)) > 0,
                        sub(";.*|,.*", "", trimws(top$overlapping.genes)),
                        coord)
  } else {
    top$label <- coord
  }
  # Disambiguate repeated gene names by appending coordinates — multiple
  # DMRs can hit the same gene's first listed name.
  dupes <- top$label[duplicated(top$label) | duplicated(top$label, fromLast = TRUE)]
  if (length(dupes)) {
    needs_coord <- top$label %in% unique(dupes)
    top$label[needs_coord] <- paste0(top$label[needs_coord], " (",
                                       coord[needs_coord], ")")
  }
  top$label <- factor(top$label, levels = top$label[order(top$meanbetafc)])

  p <- ggplot(top, aes(x = meanbetafc, y = label, colour = meanbetafc > 0)) +
    geom_segment(aes(x = 0, xend = meanbetafc, yend = label), linewidth = 0.6) +
    geom_point(size = 2.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    scale_colour_manual(values = c("TRUE"  = unname(STATUS_COLORS["Post"]),
                                   "FALSE" = unname(STATUS_COLORS["Pre"])),
                        labels = c("TRUE"="Post > Pre","FALSE"="Pre > Post"),
                        name   = NULL) +
    short_labs(x = "Mean beta FC (Post - Pre)", y = NULL) +
    project_theme(base_size = 10) +
    theme(legend.position = "bottom")

  h <- max(3.5, 0.22 * top_n + 1.0)
  ggsave(file.path(dm_dir, "dmr_top_regions.png"), p,
         width = 5, height = h, dpi = 300, bg = "white")
  cat("Saved: dmr_top_regions.png\n")
}

# ── DMRplot for top region ────────────────────────────────────────────────
if (!is.null(dmr_a) && !is.null(df_a) && nrow(df_a) > 0) {
  # DMRcate ≥ 3.x renamed DMRplot → DMR.plot. Pick whichever exists.
  dmr_plot_fn <- if (exists("DMR.plot", where = "package:DMRcate",
                            inherits = FALSE)) DMRcate::DMR.plot else
                 if (exists("DMRplot",  where = "package:DMRcate",
                            inherits = FALSE)) DMRcate::DMRplot  else NULL

  if (is.null(dmr_plot_fn)) {
    cat("DMRplot skipped: neither DMR.plot nor DMRplot exported by DMRcate.\n")
  } else tryCatch({
    top_gene <- if ("overlapping.genes" %in% colnames(df_a))
      sub(";.*|,.*", "", trimws(df_a$overlapping.genes[1])) else "top_DMR"
    top_gene <- gsub("[^A-Za-z0-9_]", "_", top_gene)
    fname    <- sprintf("dmrplot_%s.png", top_gene)
    png(file.path(dm_dir, fname), width = 1000, height = 550, res = 120)
    # Use positional args for the DMResults object — DMRcate ≥ 3.x renamed
    # the first argument from `cpgr` to `ranges`.
    dmr_plot_fn(
      dmr_a, dmr = 1,
      phen.col = ifelse(pheno$status == "Post",
                        unname(STATUS_COLORS["Post"]),
                        unname(STATUS_COLORS["Pre"])),
      genome   = "hg38"
    )
    dev.off()
    cat(sprintf("Saved: %s\n", fname))
  }, error = function(e) cat(sprintf("DMRplot skipped: %s\n", e$message)))
}

cat("\n=== Done ===\n")
