#!/usr/bin/env Rscript
#
# age_confound_fmls002.R
#
# 1. Checks FMLS_002 beta distribution vs all other samples (outlier diagnosis)
# 2. Tests whether age explains PC1 as well as menopausal status does
#    Key: if age drives the signal, within-Pre group PC1 should increase
#    monotonically with age. If it doesn't, biology (menopause) is the driver.

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

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

cat("=== FMLS_002 QC + Age Confound Assessment ===\n\n")

# ── Phenotype ──────────────────────────────────────────────────────────────
pheno <- data.frame(
  sample_id = c("207758450138_R03C01","207758450138_R04C01",
                "207758450138_R07C01","207758450138_R08C01",
                "208356980025_R02C01","208356980025_R03C01",
                "208356980025_R05C01","208356980025_R06C01"),
  subject   = c("FMLS_007","FMLS_009","FMLS_005","FMLS_006",
                "FMLS_011","FMLS_015","FMLS_002","FMLS_010"),
  status    = c("Post","Pre","Post","Post","Pre","Pre","Pre","Pre"),
  chip      = c("207758450138","207758450138","207758450138","207758450138",
                "208356980025","208356980025","208356980025","208356980025"),
  age       = c(61, 31, 62, 66, 45, 23, 27, 41),
  flag      = c("","OC_pill","","ovary_removed","PREGNANT","OC_pill","",""),
  stringsAsFactors = FALSE
)

# ── Load betas ─────────────────────────────────────────────────────────────
cat("Loading beta matrix...\n")
beta_df <- read.csv(betas_path, row.names = 1, check.names = FALSE)
beta_df <- beta_df[, pheno$sample_id, drop = FALSE]
colnames(beta_df) <- pheno$subject
cat(sprintf("  %d probes x %d samples\n\n", nrow(beta_df), ncol(beta_df)))

# ── 1. FMLS_002 beta distribution ─────────────────────────────────────────
cat("=== 1. FMLS_002 beta distribution ===\n")

# Summary stats for each sample
for (s in colnames(beta_df)) {
  b <- beta_df[[s]]
  b <- b[!is.na(b)]
  cat(sprintf("  %-10s  n=%d  mean=%.4f  median=%.4f  sd=%.4f  frac_na=%.4f\n",
              s, length(b), mean(b), median(b), sd(b),
              mean(is.na(beta_df[[s]]))))
}

# Beta density plot - all samples overlaid
cat("\n")
melt_betas <- function(beta_df, pheno, n_probes = 50000) {
  set.seed(42)
  idx <- sample(nrow(beta_df), min(n_probes, nrow(beta_df)))
  do.call(rbind, lapply(colnames(beta_df), function(s) {
    b <- beta_df[idx, s]
    df <- data.frame(beta = b[!is.na(b)], subject = s,
                     stringsAsFactors = FALSE)
    merge(df, pheno[, c("subject","status","age","flag")], by = "subject")
  }))
}

cat("Building density data (50K probe sample)...\n")
long_df <- melt_betas(beta_df, pheno)

# Highlight FMLS_002
long_df$highlight <- ifelse(long_df$subject == "FMLS_002", "FMLS_002", "Other")
long_df$label     <- paste0(long_df$subject,
                            ifelse(long_df$flag != "", paste0(" [",long_df$flag,"]"), ""),
                            " (", long_df$status, ", age ", long_df$age, ")")

p_density <- ggplot(long_df, aes(x = beta, colour = label,
                                  linewidth = highlight)) +
  geom_density(adjust = 0.5, alpha = 0.9) +
  scale_linewidth_manual(values = c("FMLS_002" = 1.8, "Other" = 0.6),
                         guide = "none") +
  scale_colour_manual(
    values = c(
      "FMLS_002 (Pre, age 27)"           = "#e6550d",
      "FMLS_005 (Post, age 62)"          = "#d62728",
      "FMLS_006 [ovary_removed] (Post, age 66)" = "#9467bd",
      "FMLS_007 (Post, age 61)"          = "#8c564b",
      "FMLS_009 [OC_pill] (Pre, age 31)" = "#17becf",
      "FMLS_010 (Pre, age 41)"           = "#1f77b4",
      "FMLS_011 [PREGNANT] (Pre, age 45)"= "#bcbd22",
      "FMLS_015 [OC_pill] (Pre, age 23)" = "#2ca02c"
    ), name = NULL) +
  short_labs(x = "Beta value (0 = unmethylated, 1 = methylated)",
             y = "Density") +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
  project_theme(base_size = 10) +
  theme(legend.text      = element_text(size = 7),
        legend.key.size  = unit(0.35, "cm"),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.margin    = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0))

ggsave(file.path(qc_dir, "fmls002_beta_density.png"), p_density,
       width = 7.2, height = 5.2, dpi = 300, bg = "white")
cat("  Saved: fmls002_beta_density.png\n\n")

# ── 2. Age vs menopausal status as driver of PC1 ──────────────────────────
cat("=== 2. Age vs menopause: which drives PC1? ===\n\n")

# Recompute PCA (top 10K variable probes)
row_vars   <- apply(beta_df, 1, var, na.rm = TRUE)
top_probes <- order(row_vars, decreasing = TRUE)[1:10000]
betas_hv   <- beta_df[top_probes, ]

betas_imp <- apply(betas_hv, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE); x
})
betas_imp[is.na(betas_imp)] <- 0.5

pca_res <- prcomp(betas_imp, center = TRUE, scale. = FALSE)
pca_df  <- as.data.frame(pca_res$x[, 1:3])
pca_df$subject <- rownames(pca_df)
pca_df  <- merge(pca_df, pheno, by = "subject")
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

cat("PC1 scores by sample (sorted by age):\n")
pca_sorted <- pca_df[order(pca_df$age), c("subject","status","age","flag","PC1")]
print(pca_sorted, row.names = FALSE)
cat("\n")

# Model 1: PC1 ~ age (continuous)
m_age    <- lm(PC1 ~ age, data = pca_df)
r2_age   <- summary(m_age)$r.squared

# Model 2: PC1 ~ status (binary)
m_status <- lm(PC1 ~ status, data = pca_df)
r2_status <- summary(m_status)$r.squared

cat(sprintf("R² of PC1 ~ age:    %.4f (%.1f%%)\n", r2_age,    r2_age*100))
cat(sprintf("R² of PC1 ~ status: %.4f (%.1f%%)\n", r2_status, r2_status*100))
cat("\n")

# Within-Pre group: does PC1 correlate with age?
pre_df <- pca_df[pca_df$status == "Pre", ]
cor_pre <- cor(pre_df$PC1, pre_df$age, method = "spearman")
cat(sprintf("Within Pre group — Spearman correlation of PC1 with age: r=%.4f\n", cor_pre))
cat("  (If age drives signal: expect r strongly negative — older Pre → lower PC1,\n")
cat("   closer to Post. If r is near 0 or positive, age is NOT the driver.)\n\n")

cat("Within Pre group (age vs PC1):\n")
print(pre_df[order(pre_df$age), c("subject","age","flag","PC1")], row.names = FALSE)
cat("\n")

# Within-Post group: does PC1 correlate with age?
post_df <- pca_df[pca_df$status == "Post", ]
cor_post <- cor(post_df$PC1, post_df$age, method = "spearman")
cat(sprintf("Within Post group — Spearman correlation of PC1 with age: r=%.4f\n", cor_post))
cat("\n")

# ── 3. PC1 vs age scatter plot ─────────────────────────────────────────────
p_age <- ggplot(pca_df, aes(x = age, y = PC1, colour = status, label = subject)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              linetype = "dashed", linewidth = 0.5) +
  # Within-group trend lines
  geom_smooth(aes(group = status, colour = status),
              method = "lm", se = FALSE, linewidth = 0.6, linetype = "dotted") +
  geom_point(size = 2.6) +
  geom_text(vjust = -0.9, size = 2.6) +
  scale_colour_manual(values = STATUS_COLORS, name = "Status") +
  # Expand axes so the labels at the corners aren't clipped.
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  # Two-line annotation anchored to the top-left interior so it never
  # collides with the right panel edge (previous version clipped to
  # "R²(age) = 0.68  R²(s..." because the string ran off the frame).
  annotate("text", x = -Inf, y = Inf,
           label = sprintf(
             "italic(R)^2*'(age)' == %.2f ~~ italic(R)^2*'(status)' == %.2f",
             r2_age, r2_status),
           size = 2.6, hjust = -0.05, vjust = 1.6, colour = "grey30",
           parse = TRUE) +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf(
             "italic(r)[Pre] == %.2f ~~ italic(r)[Post] == %.2f",
             cor_pre, cor_post),
           size = 2.6, hjust = -0.05, vjust = 3.2, colour = "grey30",
           parse = TRUE) +
  short_labs(
    x = "Age (years)",
    y = sprintf("PC1 (%.1f%% variance)", var_exp[1])
  ) +
  project_theme(base_size = 10)

ggsave(file.path(qc_dir, "pc1_vs_age.png"), p_age,
       width = 5.0, height = 3.8, dpi = 300, bg = "white")
cat("  Saved: pc1_vs_age.png\n")

# ── Summary ────────────────────────────────────────────────────────────────
cat("\n=== SUMMARY ===\n")
cat(sprintf("  Age alone explains %.1f%% of PC1 variance\n", r2_age * 100))
cat(sprintf("  Status alone explains %.1f%% of PC1 variance\n", r2_status * 100))
cat(sprintf("  Within-Pre age-PC1 correlation: r=%.3f\n", cor_pre))

if (abs(cor_pre) < 0.3) {
  cat("  -> Within Pre, PC1 does NOT track age → age is not the primary driver\n")
} else if (cor_pre < -0.3) {
  cat("  -> Within Pre, older Pre samples shift toward Post → age may confound\n")
} else {
  cat("  -> Within Pre, older Pre samples shift AWAY from Post → age not confounding\n")
}
