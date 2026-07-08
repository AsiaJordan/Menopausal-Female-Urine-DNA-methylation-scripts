#!/usr/bin/env Rscript
#
# neutrophil_dmp_overlap.R
#
# Two complementary analyses:
#
# 1. REFERENCE-BASED: For CpGs in centDHSbloodDMC.m (333 probes), does the
#    neutrophil reference beta value predict the limma logFC (Post vs Pre)?
#    If DMPs are neutrophil-driven, probes with high Neutro beta should be
#    hypermethylated in Post (positive logFC), while CD4T/NK/B-high probes
#    should trend the other way.
#
# 2. DATA-DRIVEN: Correlate each probe's beta values with the IC fraction
#    across 8 samples. IC-correlated probes = neutrophil-influenced probes.
#    Show these dominate the DMP list.
#
# Outputs:
#   neutrophil_reference_scatter.png   (Neutro ref beta vs logFC)
#   neutrophil_ic_correlation.png      (IC correlation vs logFC)
#   neutrophil_dmp_overlap.txt
#
# Usage:
#   Rscript neutrophil_dmp_overlap.R <betas_full.csv> <output_dir>

args       <- commandArgs(trailingOnly = TRUE)
betas_path <- path.expand(if (length(args) >= 1) args[1] else
  "menopause_urine_out/betas_full.csv")
out_dir    <- path.expand(if (length(args) >= 2) args[2] else
  "menopause_urine_out")

# ── Section subfolders ───────────────────────────────────────────────────────
qc_dir      <- file.path(out_dir, "01_qc")
deconv_dir  <- file.path(out_dir, "02_deconvolution")
dm_dir      <- file.path(out_dir, "03_diff_meth")
fig_dir     <- file.path(out_dir, "04_figures")
pub_fig_dir <- file.path(out_dir, "05_pub_figures")
for (d in c(qc_dir, deconv_dir, dm_dir, fig_dir, pub_fig_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("=== Neutrophil-specific CpG overlap with menopausal DMPs ===\n\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(EpiDISH)
})
if (!requireNamespace("patchwork", quietly=TRUE))
  install.packages("patchwork", repos="https://cloud.r-project.org")
suppressPackageStartupMessages(library(patchwork))

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

# ── Phenotype + IC fractions ──────────────────────────────────────────────
pheno <- data.frame(
  sample_id = c("207758450138_R03C01","207758450138_R04C01",
                "207758450138_R07C01","207758450138_R08C01",
                "208356980025_R02C01","208356980025_R03C01",
                "208356980025_R05C01","208356980025_R06C01"),
  subject   = c("FMLS_007","FMLS_009","FMLS_005","FMLS_006",
                "FMLS_011","FMLS_015","FMLS_002","FMLS_010"),
  status    = c("Post","Pre","Post","Post","Pre","Pre","Pre","Pre"),
  stringsAsFactors = FALSE
)
ic_df <- read.csv(file.path(deconv_dir, "epidish_fractions.csv"))
pheno <- merge(pheno, ic_df[, c("subject","IC")], by="subject")
pheno <- pheno[match(c("FMLS_007","FMLS_009","FMLS_005","FMLS_006",
                        "FMLS_011","FMLS_015","FMLS_002","FMLS_010"),
                     pheno$subject), ]

# ── Load DMP lists ────────────────────────────────────────────────────────
tt_unadj <- read.csv(file.path(dm_dir, "dmp_list.csv"),          row.names=1)
tt_adj   <- read.csv(file.path(dm_dir, "limma_IC_adjusted.csv"), row.names=1)
tt_unadj <- tt_unadj[order(tt_unadj$P.Value), ]
tt_adj   <- tt_adj[order(tt_adj$P.Value), ]

all_probes <- rownames(tt_unadj)
base_ids   <- sub("_[^_]+$", "", all_probes)
cat(sprintf("Unadjusted limma: %d probes\nIC-adjusted limma: %d probes\n\n",
            nrow(tt_unadj), nrow(tt_adj)))

# ── Load blood reference ──────────────────────────────────────────────────
cat("Loading centDHSbloodDMC.m...\n")
data(centDHSbloodDMC.m)
cat(sprintf("  %d CpGs x %d cell types: %s\n",
            nrow(centDHSbloodDMC.m), ncol(centDHSbloodDMC.m),
            paste(colnames(centDHSbloodDMC.m), collapse=", ")))

# Map reference CpGs to our EPIC v2 array
ref_base    <- rownames(centDHSbloodDMC.m)
matched_idx <- match(ref_base, base_ids)
in_array    <- !is.na(matched_idx)
cat(sprintf("  %d / %d reference CpGs in our array (%.1f%%)\n\n",
            sum(in_array), length(ref_base),
            100*sum(in_array)/length(ref_base)))

ref_probes_full <- all_probes[matched_idx[in_array]]   # full probe IDs
ref_sub         <- centDHSbloodDMC.m[in_array, ]       # matched reference rows

# ══════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: Reference-based — Neutro beta vs limma logFC
# ══════════════════════════════════════════════════════════════════════════
cat("=== Analysis 1: Neutro reference beta vs limma logFC ===\n\n")

# Build data frame of reference beta values + limma stats
ref_df <- as.data.frame(ref_sub)
ref_df$probe    <- ref_probes_full
ref_df$logFC    <- tt_unadj[ref_probes_full, "logFC"]
ref_df$logFC_adj <- tt_adj[ref_probes_full, "logFC"]
ref_df$pval     <- tt_unadj[ref_probes_full, "P.Value"]

# Which cell type has the highest beta? (defines cell-type assignment)
ct_cols <- colnames(centDHSbloodDMC.m)
ref_df$dominant_ct <- ct_cols[apply(ref_sub, 1, which.max)]
ref_df$dominant_ct <- factor(ref_df$dominant_ct,
                              levels = c("Neutro","CD4T","CD8T","Mono","NK","B","Eosino"))

cat("Cell-type assignments in reference (dominant beta):\n")
print(table(ref_df$dominant_ct))
cat("\n")

# Spearman correlation: Neutro beta vs logFC
cor_neutro <- cor(ref_df$Neutro, ref_df$logFC,    method="spearman")
cor_cd4    <- cor(ref_df$CD4T,   ref_df$logFC,    method="spearman")
cor_mono   <- cor(ref_df$Mono,   ref_df$logFC,    method="spearman")
cor_nk     <- cor(ref_df$NK,     ref_df$logFC,    method="spearman")

cat("Spearman r(reference_beta, unadjusted_logFC):\n")
cat(sprintf("  Neutro: r = %+.3f  (positive = neutro-high CpGs are Post>Pre)\n", cor_neutro))
cat(sprintf("  CD4T:   r = %+.3f\n", cor_cd4))
cat(sprintf("  Mono:   r = %+.3f\n", cor_mono))
cat(sprintf("  NK:     r = %+.3f\n\n", cor_nk))

cor_neutro_adj <- cor(ref_df$Neutro, ref_df$logFC_adj, method="spearman")
cat(sprintf("Spearman r(Neutro, IC-adjusted logFC): r = %+.3f\n\n", cor_neutro_adj))

# ── Plot 1: Neutro reference beta vs logFC, coloured by dominant cell type ─
ct_colours <- CELL_COLORS[c("Neutro","CD4T","CD8T","Mono","NK","B","Eosino")]

p_scatter <- ggplot(ref_df, aes(x = Neutro, y = logFC, colour = dominant_ct)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.6, alpha = 0.85) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
              colour = "black", linewidth = 0.5, alpha = 0.15) +
  scale_colour_manual(values = ct_colours, name = NULL) +
  annotate("text", x = Inf, y = Inf,
           label = sprintf("italic(r) == %+.2f", cor_neutro),
           hjust = 1.1, vjust = 1.5, size = 2.8, parse = TRUE) +
  short_labs(
    x = "Neutrophil reference beta",
    y = "logFC (Post - Pre)"
  ) +
  project_theme(base_size = 10)

# ── Plot 2: All 7 cell types in one faceted scatter ───────────────────────
ref_long <- do.call(rbind, lapply(ct_cols, function(ct) {
  data.frame(
    probe   = ref_df$probe,
    ct      = ct,
    ref_beta = ref_df[[ct]],
    logFC   = ref_df$logFC,
    stringsAsFactors = FALSE
  )
}))
ref_long$ct <- factor(ref_long$ct, levels = names(ct_colours))

cor_labels <- vapply(ct_cols, function(ct) {
  r <- cor(ref_df[[ct]], ref_df$logFC, method="spearman")
  sprintf("r = %+.2f", r)
}, character(1))

p_facet <- ggplot(ref_long, aes(x = ref_beta, y = logFC, colour = ct)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(size = 1.4, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "black",
              linewidth = 0.5, linetype = "solid") +
  geom_text(data = data.frame(ct = factor(ct_cols, levels=names(ct_colours)),
                               lab = cor_labels,
                               ref_beta = Inf, logFC = Inf),
            aes(x = ref_beta, y = logFC, label = lab),
            hjust = 1.1, vjust = 1.5, size = 2.6,
            colour = "black", inherit.aes = FALSE) +
  scale_colour_manual(values = ct_colours, guide = "none") +
  facet_wrap(~ct, ncol=4) +
  short_labs(
    x = "Reference beta",
    y = "logFC (Post - Pre)"
  ) +
  project_theme(base_size = 10)

# ══════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: Data-driven — IC correlation across samples
# ══════════════════════════════════════════════════════════════════════════
cat("=== Analysis 2: Data-driven IC-probe correlations ===\n\n")

cat("Loading beta matrix...\n")
beta_df  <- read.csv(betas_path, row.names=1, check.names=FALSE)
beta_df  <- beta_df[, pheno$sample_id, drop=FALSE]
colnames(beta_df) <- pheno$subject
beta_mat <- as.matrix(beta_df)
beta_mat <- beta_mat[rowSums(is.na(beta_mat)) == 0, ]
cat(sprintf("  %d x %d (NAs removed)\n\n", nrow(beta_mat), ncol(beta_mat)))

# Spearman correlation of each probe with IC fraction
# (IC ≈ neutrophil fraction given Post = ~98% neutrophils)
cat("Computing probe-IC Spearman correlations (may take ~30s)...\n")
ic_vec <- pheno$IC
ic_rank <- rank(ic_vec)
ic_r <- apply(beta_mat, 1, function(x) {
  if (var(x, na.rm=TRUE) < 1e-10) return(0)
  cor(rank(x), ic_rank, method="pearson")  # Spearman via ranks
})
cat(sprintf("  Done. IC correlation range: %.3f to %.3f\n\n", min(ic_r), max(ic_r)))

# Link IC correlation to limma logFC and DMP rank
ic_df2 <- data.frame(
  probe    = names(ic_r),
  ic_r     = ic_r,
  logFC    = tt_unadj[names(ic_r), "logFC"],
  pval     = tt_unadj[names(ic_r), "P.Value"],
  stringsAsFactors = FALSE
)
ic_df2 <- ic_df2[!is.na(ic_df2$logFC), ]
ic_df2$dmp_rank <- rank(ic_df2$pval)

# Correlation between IC-r and logFC
cat(sprintf("Spearman r(IC_correlation, logFC): %+.4f\n", 
            cor(ic_df2$ic_r, ic_df2$logFC, method="spearman")))

# What fraction of top DMPs have |IC_r| > 0.7?
thresholds_ic <- c(0.5, 0.7, 0.9)
top_ns <- c(500, 1000, 5000, 10000)
bg_frac <- sapply(thresholds_ic, function(t) mean(abs(ic_df2$ic_r) > t))

cat("\nFraction of probes with |IC_r| above threshold:\n")
cat(sprintf("  Background (all probes):  |r|>0.5: %.1f%%  |r|>0.7: %.1f%%  |r|>0.9: %.1f%%\n\n",
            100*bg_frac[1], 100*bg_frac[2], 100*bg_frac[3]))

cat(sprintf("  %-12s  %6s  %6s  %6s\n", "DMP set", "|r|>0.5", "|r|>0.7", "|r|>0.9"))
for (n in top_ns) {
  topn <- ic_df2$probe[order(ic_df2$pval)][seq_len(n)]
  fracs <- sapply(thresholds_ic, function(t)
    100 * mean(abs(ic_df2[topn, "ic_r"]) > t, na.rm=TRUE))
  cat(sprintf("  top %-8d  %5.1f%%  %5.1f%%  %5.1f%%\n", n, fracs[1], fracs[2], fracs[3]))
}
for (pt in c(0.05, 0.01, 0.001)) {
  sig  <- ic_df2$probe[ic_df2$pval < pt]
  fracs <- sapply(thresholds_ic, function(t)
    100 * mean(abs(ic_df2[sig, "ic_r"]) > t, na.rm=TRUE))
  cat(sprintf("  p<%-9s  %5.1f%%  %5.1f%%  %5.1f%%\n", pt, fracs[1], fracs[2], fracs[3]))
}

# ── Plot 3: IC correlation vs logFC (hexbin density) ─────────────────────
# Subsample for plotting speed
set.seed(42)
samp_idx <- sample(nrow(ic_df2), min(100000, nrow(ic_df2)))
ic_plot  <- ic_df2[samp_idx, ]

p_ic <- ggplot(ic_plot, aes(x = ic_r, y = logFC)) +
  geom_hex(bins = 80, aes(fill = after_stat(log10(count)))) +
  scale_fill_gradient(low = "#deebf7", high = "#08306b",
                      name = expression(log[10]*" n")) +
  geom_hline(yintercept = 0,  linetype="dashed", colour="grey50") +
  geom_vline(xintercept = 0,  linetype="dashed", colour="grey50") +
  geom_vline(xintercept =  0.7, linetype="dotted",
             colour = unname(STATUS_COLORS["Post"]), linewidth=0.5) +
  geom_vline(xintercept = -0.7, linetype="dotted",
             colour = unname(STATUS_COLORS["Pre"]),  linewidth=0.5) +
  annotate("text", x= 0.85, y= max(ic_plot$logFC, na.rm=TRUE)*0.9,
           label="IC-correlated", size=2.6,
           colour = unname(STATUS_COLORS["Post"])) +
  annotate("text", x=-0.85, y= max(ic_plot$logFC, na.rm=TRUE)*0.9,
           label="IC-anticorrelated", size=2.6,
           colour = unname(STATUS_COLORS["Pre"])) +
  short_labs(
    x = "r(probe beta, IC fraction)",
    y = "logFC (Post - Pre)"
  ) +
  project_theme(base_size = 10)

# ── Plot 4: IC|r| distribution in top DMPs vs background ─────────────────
top5k_mask <- ic_df2$probe %in% ic_df2$probe[order(ic_df2$pval)][1:5000]
ic_df2$group <- ifelse(top5k_mask, "Top 5000 DMPs", "Background")
ic_df2$group <- factor(ic_df2$group, levels=c("Background","Top 5000 DMPs"))

p_dist <- ggplot(ic_df2, aes(x = abs(ic_r), fill = group, colour = group)) +
  geom_density(alpha = 0.45, linewidth = 0.5) +
  scale_fill_manual(values = c("Background"="#aec7e8",
                                "Top 5000 DMPs" = unname(STATUS_COLORS["Post"]))) +
  scale_colour_manual(values = c("Background" = unname(STATUS_COLORS["Pre"]),
                                  "Top 5000 DMPs" = "#9e0000")) +
  short_labs(
    x = "|r(probe beta, IC fraction)|",
    y = "Density",
    fill = NULL, colour = NULL
  ) +
  project_theme(base_size = 10) +
  theme(legend.position = "top")

# ── Combine all 4 plots ───────────────────────────────────────────────────
combined <- (p_scatter | p_facet) / (p_ic | p_dist)

ggsave(file.path(fig_dir, "neutrophil_dmp_overlap.png"),
       combined, width=7.5, height=6, dpi=300, bg = "white")
cat("\nSaved: neutrophil_dmp_overlap.png\n")

# Also save individual plots
ggsave(file.path(fig_dir, "neutro_ref_scatter.png"),   p_scatter, width=4.5, height=3.3, dpi=300, bg = "white")
ggsave(file.path(fig_dir, "neutro_ref_facet.png"),     p_facet,   width=7.2, height=4,   dpi=300, bg = "white")
ggsave(file.path(fig_dir, "ic_correlation_logfc.png"), p_ic,      width=4.5, height=3.3, dpi=300, bg = "white")
ggsave(file.path(fig_dir, "ic_correlation_dist.png"),  p_dist,    width=4.5, height=3.3, dpi=300, bg = "white")
cat("Saved individual plots.\n\n")

cat("=== Done ===\n")
