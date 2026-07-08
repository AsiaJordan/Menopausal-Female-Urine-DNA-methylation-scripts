#!/usr/bin/env Rscript
#
# epidish_menopause_urine.R
#
# Cell-type deconvolution of urine EPIC v2 methylation (n=8 pilot).
#
# Reference panels:
#   centEpiFibIC.m  — 3-component (Epithelial / Fibroblast / Immune)
#                     Epithelial ≈ urothelial/transitional cells
#                     IC         ≈ leukocytes (central to our hypothesis)
#   centDHSbloodDMC.m — 7 blood cell types (B, NK, CD4T, CD8T, Mono, Neu, Eos)
#                       gives finer breakdown of the immune component
#
# Outputs:
#   epidish_fractions.csv       — all estimated fractions, both references
#   epidish_barplot.png         — stacked bar chart by sample
#   epidish_prepost.png         — Pre vs Post IC fraction comparison
#   epidish_pc1_correlation.png — IC fraction vs PC1 scatter
#   limma_IC_adjusted.csv       — DMPs after adjusting for IC fraction
#   epidish_summary.txt         — text summary for paper/response
#
# Usage:
#   Rscript epidish_menopause_urine.R \
#       <betas_full.csv> <output_dir>

suppressPackageStartupMessages({
  library(ggplot2)
  library(limma)
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
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── Section subfolders ───────────────────────────────────────────────────────
qc_dir      <- file.path(out_dir, "01_qc")
deconv_dir  <- file.path(out_dir, "02_deconvolution")
dm_dir      <- file.path(out_dir, "03_diff_meth")
fig_dir     <- file.path(out_dir, "04_figures")
pub_fig_dir <- file.path(out_dir, "05_pub_figures")
for (d in c(qc_dir, deconv_dir, dm_dir, fig_dir, pub_fig_dir))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

cat("=== EpiDish Cell-Type Deconvolution: Menopause Urine ===\n\n")

# ── Install EpiDISH if needed ──────────────────────────────────────────────
if (!requireNamespace("EpiDISH", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install("EpiDISH", ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages(library(EpiDISH))

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

# ── Load full beta matrix ──────────────────────────────────────────────────
cat("Loading beta matrix...\n")
beta_df <- read.csv(betas_path, row.names = 1, check.names = FALSE)
beta_df <- beta_df[, pheno$sample_id, drop = FALSE]
colnames(beta_df) <- pheno$subject
cat(sprintf("  %d probes x %d samples\n\n", nrow(beta_df), ncol(beta_df)))

beta_mat <- as.matrix(beta_df)   # probes x samples

# ── EPIC v2 probe suffix handling ──────────────────────────────────────────
# EPIC v2 SeSAMe output has all probes suffixed (e.g. cg00574958_TC21).
# EpiDish references use plain 450K/EPIC cg IDs. Build a lookup table so we
# can subset the beta matrix to the reference probes by base cg ID.
array_probes   <- rownames(beta_mat)
array_base_ids <- sub("_[^_]+$", "", array_probes)  # strip last _XXXX
# named vector: base_id → suffixed probe id
base_to_full   <- setNames(array_probes, array_base_ids)

subset_to_ref <- function(ref_matrix) {
  ref_ids <- rownames(ref_matrix)
  # which ref probes have a match in our array (by base ID)?
  matched_base <- ref_ids[ref_ids %in% array_base_ids]
  if (length(matched_base) == 0) stop("No probe overlap after suffix stripping")

  full_ids <- base_to_full[matched_base]   # suffixed IDs for row subsetting
  beta_sub <- beta_mat[full_ids, , drop = FALSE]
  rownames(beta_sub) <- matched_base       # restore plain cg IDs for EpiDish

  ref_sub  <- ref_matrix[matched_base, , drop = FALSE]

  cat(sprintf("  Reference probes: %d | In array (after suffix strip): %d (%.1f%%)\n",
              nrow(ref_matrix), length(matched_base),
              100 * length(matched_base) / nrow(ref_matrix)))

  # Impute NAs with probe median
  na_count <- sum(is.na(beta_sub))
  if (na_count > 0) {
    probe_meds <- apply(beta_sub, 1, median, na.rm = TRUE)
    for (i in seq_len(nrow(beta_sub))) {
      nas <- is.na(beta_sub[i, ])
      if (any(nas)) beta_sub[i, nas] <- probe_meds[i]
    }
    cat(sprintf("  Imputed %d NAs\n", na_count))
  }

  list(beta = beta_sub, ref = ref_sub)
}

# ── 1. centEpiFibIC.m — 3-compartment (primary analysis) ──────────────────
cat("=== 1. EpiDish: centEpiFibIC.m (Epi / Fib / Immune) ===\n")
data(centEpiFibIC.m)
d3 <- subset_to_ref(centEpiFibIC.m)

res3 <- epidish(beta.m = d3$beta, ref.m = d3$ref, method = "RPC")
frac3 <- as.data.frame(res3$estF)
colnames(frac3) <- c("Epi", "Fib", "IC")
frac3$subject <- rownames(frac3)
frac3 <- merge(frac3, pheno[, c("subject","status","age","flag","chip")],
               by = "subject")

cat("\nEstimated cell-type fractions (centEpiFibIC):\n")
print(frac3[order(frac3$status, frac3$age),
            c("subject","status","age","flag","Epi","Fib","IC")],
      row.names = FALSE, digits = 3)
cat("\n")

# Group means
cat("Mean fractions by menopausal status:\n")
for (g in c("Pre","Post")) {
  sub <- frac3[frac3$status == g, ]
  cat(sprintf("  %s (n=%d): Epi=%.3f  Fib=%.3f  IC=%.3f\n",
              g, nrow(sub), mean(sub$Epi), mean(sub$Fib), mean(sub$IC)))
}
cat("\n")

# ── 2. centDHSbloodDMC.m — 7 blood cell types (immune breakdown) ──────────
cat("=== 2. EpiDish: centDHSbloodDMC.m (7 blood cell types) ===\n")
data(centDHSbloodDMC.m)
d7 <- subset_to_ref(centDHSbloodDMC.m)

res7 <- epidish(beta.m = d7$beta, ref.m = d7$ref, method = "RPC")
frac7 <- as.data.frame(res7$estF)
frac7$subject <- rownames(frac7)
frac7 <- merge(frac7, pheno[, c("subject","status","age","flag")],
               by = "subject")

cat("\nEstimated blood cell fractions:\n")
print(frac7[order(frac7$status, frac7$age), ], row.names = FALSE, digits = 3)
cat("\n")

cat("Mean fractions by menopausal status:\n")
blood_cols <- setdiff(colnames(frac7), c("subject","status","age","flag"))
for (g in c("Pre","Post")) {
  sub <- frac7[frac7$status == g, blood_cols, drop = FALSE]
  vals <- colMeans(sub)
  cat(sprintf("  %s: %s\n", g,
              paste(sprintf("%s=%.3f", names(vals), vals), collapse="  ")))
}
cat("\n")

# ── 3. Statistical tests: Pre vs Post ─────────────────────────────────────
cat("=== 3. Pre vs Post comparison ===\n\n")

pre_IC  <- frac3$IC[frac3$status == "Pre"]
post_IC <- frac3$IC[frac3$status == "Post"]

cat(sprintf("IC fraction — Pre:  mean=%.4f  median=%.4f  range=[%.4f, %.4f]\n",
            mean(pre_IC), median(pre_IC), min(pre_IC), max(pre_IC)))
cat(sprintf("IC fraction — Post: mean=%.4f  median=%.4f  range=[%.4f, %.4f]\n",
            mean(post_IC), median(post_IC), min(post_IC), max(post_IC)))

wt <- wilcox.test(post_IC, pre_IC, alternative = "greater", exact = FALSE)
tt_ic <- t.test(post_IC, pre_IC, alternative = "greater")
cat(sprintf("\nWilcoxon one-sided (Post > Pre): W=%.0f, p=%.4f\n",
            wt$statistic, wt$p.value))
cat(sprintf("t-test one-sided   (Post > Pre): t=%.3f, p=%.4f\n\n",
            tt_ic$statistic, tt_ic$p.value))

# ── 4. Correlation of IC with PC1 ─────────────────────────────────────────
cat("=== 4. IC fraction vs PC1 ===\n")

# Recompute PCA (top 10K variable probes)
row_vars   <- apply(beta_mat, 1, var, na.rm = TRUE)
top_probes <- order(row_vars, decreasing = TRUE)[1:10000]
betas_hv   <- beta_mat[top_probes, ]
betas_imp  <- apply(betas_hv, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE); x
})
betas_imp[is.na(betas_imp)] <- 0.5

pca_res <- prcomp(betas_imp, center = TRUE, scale. = FALSE)
pca_df  <- as.data.frame(pca_res$x[, 1:3])
pca_df$subject <- rownames(pca_df)
pca_df  <- merge(pca_df, frac3[, c("subject","IC","Epi","Fib",
                                    "status","age","flag")], by = "subject")
var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)

cor_IC_PC1    <- cor(pca_df$IC, pca_df$PC1, method = "spearman")
cor_IC_PC1_p  <- cor(pca_df$IC, pca_df$PC1, method = "pearson")
ct_pearson    <- cor.test(pca_df$IC, pca_df$PC1, method = "pearson")
cat(sprintf("Spearman r(IC, PC1) = %.4f\n", cor_IC_PC1))
cat(sprintf("Pearson  r(IC, PC1) = %.4f  (p = %.2e, r^2 = %.4f)\n",
            cor_IC_PC1_p, ct_pearson$p.value, cor_IC_PC1_p^2))
cat(sprintf("PC1 explains %.1f%% of total variance\n\n", var_exp[1]))

# ── 5. Limma adjusted for IC fraction ─────────────────────────────────────
cat("=== 5. Limma adjusted for IC fraction ===\n")

# Impute full matrix for limma
betas_all_t <- apply(beta_mat, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE); x
})
betas_all_t[is.na(betas_all_t)] <- 0.5
beta_for_limma <- t(betas_all_t)   # probes x samples

# Merge IC fractions into pheno (in correct sample order)
pheno2 <- merge(pheno, frac3[, c("subject","IC")], by = "subject")
pheno2 <- pheno2[match(colnames(beta_for_limma), pheno2$subject), ]

pheno2$status_f <- factor(pheno2$status, levels = c("Pre","Post"))
design_adj <- model.matrix(~ IC + status_f, data = pheno2)
cat("Design matrix (IC-adjusted):\n")
print(design_adj)
cat("\n")

fit_adj  <- lmFit(beta_for_limma, design_adj)
fit_adj  <- eBayes(fit_adj, trend = TRUE)
tt_adj   <- topTable(fit_adj, coef = "status_fPost",
                     number = Inf, sort.by = "P", confint = TRUE)

cat(sprintf("IC-adjusted limma — probes nominal p<0.05: %d (%.1f%%)\n",
            sum(tt_adj$P.Value < 0.05),
            100 * mean(tt_adj$P.Value < 0.05)))
cat(sprintf("IC-adjusted limma — probes FDR<0.05:      %d\n",
            sum(tt_adj$adj.P.Val < 0.05)))
cat(sprintf("IC-adjusted limma — probes FDR<0.20:      %d\n\n",
            sum(tt_adj$adj.P.Val < 0.20)))

cat("Top 10 IC-adjusted DMPs:\n")
print(head(tt_adj[, c("logFC","t","P.Value","adj.P.Val")], 10))
cat("\n")

write.csv(tt_adj, file.path(dm_dir, "limma_IC_adjusted.csv"))
cat("  Saved: limma_IC_adjusted.csv\n\n")

# ── 6. Plots ───────────────────────────────────────────────────────────────
cat("Generating plots...\n")

# 6a. Stacked bar chart (centEpiFibIC)
long3 <- do.call(rbind, lapply(c("Epi","Fib","IC"), function(ct) {
  data.frame(
    subject   = frac3$subject,
    status    = frac3$status,
    age       = frac3$age,
    flag      = frac3$flag,
    cell_type = ct,
    fraction  = frac3[[ct]],
    stringsAsFactors = FALSE
  )
}))
long3$cell_type <- factor(long3$cell_type, levels = c("Epi","Fib","IC"))
long3$label <- paste0(long3$subject,
                      ifelse(long3$flag != "", paste0("\n[", long3$flag, "]"), ""),
                      "\n(", long3$status, ", age ", long3$age, ")")
long3$label <- factor(long3$label,
                       levels = unique(long3$label[order(long3$status,
                                                         long3$age)]))

p_bar <- ggplot(long3, aes(x = label, y = fraction, fill = cell_type)) +
  geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = CELL_COLORS[c("Epi","Fib","IC")],
                    name   = NULL) +
  facet_grid(~ status, scales = "free_x", space = "free_x") +
  short_labs(x = NULL, y = "Estimated cell fraction") +
  project_theme(base_size = 10) +
  theme(axis.text.x = element_text(size = 7, angle = 15, hjust = 1))

ggsave(file.path(deconv_dir, "epidish_barplot.png"), p_bar,
       width = 7.5, height = 4, dpi = 300, bg = "white")
cat("  Saved: epidish_barplot.png\n")

# 6b. IC fraction Pre vs Post
p_IC <- ggplot(frac3, aes(x = status, y = IC, colour = status,
                            label = subject)) +
  geom_jitter(width = 0.05, size = 2.6) +
  geom_text(vjust = -0.9, size = 2.6, show.legend = FALSE) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.35,
               colour = "grey20", linewidth = 0.4) +
  scale_colour_manual(values = STATUS_COLORS, guide = "none") +
  annotate("text", x = 1.5, y = max(frac3$IC) * 1.05,
           label = sprintf("italic(p) == %.3f", wt$p.value),
           size = 3, colour = "grey20", parse = TRUE) +
  short_labs(x = NULL, y = "IC fraction (EpiDish)") +
  project_theme(base_size = 10)

ggsave(file.path(deconv_dir, "epidish_prepost.png"), p_IC,
       width = 3.3, height = 3.5, dpi = 300, bg = "white")
cat("  Saved: epidish_prepost.png\n")

# 6c. IC vs PC1 scatter
p_IC_PC1 <- ggplot(pca_df, aes(x = IC, y = PC1, colour = status,
                                 label = subject)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40",
              fill = "grey85",
              linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2.6) +
  geom_text(vjust = -1.0, size = 2.4, show.legend = FALSE) +
  scale_colour_manual(values = STATUS_COLORS, name = "Status") +
  # Expand axes so labels at the extremes aren't clipped by the panel edge.
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.10))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.10))) +
  annotate("text",
           x = min(pca_df$IC), y = min(pca_df$PC1),
           label = sprintf(
             "italic(r) == %.3f ~~ italic(r)^2 == %.3f",
             cor_IC_PC1_p, cor_IC_PC1_p^2),
           size = 3, hjust = 0, vjust = 0, colour = "grey20",
           parse = TRUE) +
  short_labs(
    x = "IC fraction (EpiDish)",
    y = sprintf("PC1 (%.1f%% variance)", var_exp[1])
  ) +
  project_theme(base_size = 10)

ggsave(file.path(deconv_dir, "epidish_pc1_correlation.png"), p_IC_PC1,
       width = 4.8, height = 3.5, dpi = 300, bg = "white")
cat("  Saved: epidish_pc1_correlation.png\n\n")

# ── 7. Save combined fractions & text summary ──────────────────────────────
frac_out <- merge(frac3[, c("subject","status","age","flag","Epi","Fib","IC")],
                  frac7, by = c("subject","status","age","flag"))
write.csv(frac_out, file.path(deconv_dir, "epidish_fractions.csv"), row.names = FALSE)
cat("  Saved: epidish_fractions.csv\n")

sink(file.path(deconv_dir, "epidish_summary.txt"))
cat("=== EpiDish Deconvolution: Menopause Urine Pilot ===\n\n")
cat("Reference: centEpiFibIC.m (Epithelial / Fibroblast / Immune cells)\n\n")
cat("Estimated IC (immune/leukocyte) fractions:\n")
print(frac3[order(frac3$status, frac3$age),
            c("subject","status","age","flag","Epi","Fib","IC")],
      row.names = FALSE, digits = 4)
cat(sprintf("\nPre-menopausal  IC: mean=%.4f  SD=%.4f\n",
            mean(pre_IC), sd(pre_IC)))
cat(sprintf("Post-menopausal IC: mean=%.4f  SD=%.4f\n",
            mean(post_IC), sd(post_IC)))
cat(sprintf("\nWilcoxon test (Post > Pre): W=%.0f, p=%.4f\n",
            wt$statistic, wt$p.value))
cat(sprintf("t-test         (Post > Pre): t=%.3f, p=%.4f\n",
            tt_ic$statistic, tt_ic$p.value))
cat(sprintf("\nSpearman r(IC, PC1) = %.4f\n", cor_IC_PC1))
cat(sprintf("(PC1 = %.1f%% of total methylation variance)\n\n", var_exp[1]))
cat("IC-adjusted limma (~ IC + status):\n")
cat(sprintf("  Probes p<0.05:  %d (%.1f%%)\n",
            sum(tt_adj$P.Value < 0.05),
            100 * mean(tt_adj$P.Value < 0.05)))
cat(sprintf("  Probes FDR<0.05:%d\n", sum(tt_adj$adj.P.Val < 0.05)))
cat(sprintf("  Probes FDR<0.20:%d\n\n", sum(tt_adj$adj.P.Val < 0.20)))
cat("INTERPRETATION:\n")
if (mean(post_IC) > mean(pre_IC)) {
  cat(sprintf("  Post-menopausal IC fraction is HIGHER (%.3f vs %.3f).\n",
              mean(post_IC), mean(pre_IC)))
  cat("  This is CONSISTENT with the paper's hypothesis:\n")
  cat("  increased urinary leukocytes in postmenopausal women.\n")
} else {
  cat(sprintf("  Pre-menopausal IC fraction is HIGHER (%.3f vs %.3f).\n",
              mean(pre_IC), mean(post_IC)))
  cat("  IC composition does not explain the menopausal methylation signal.\n")
}
if (abs(cor_IC_PC1) > 0.7) {
  cat(sprintf("\n  IC fraction strongly correlates with PC1 (r=%.3f).\n",
              cor_IC_PC1))
  cat("  Cell composition partially co-varies with menopausal status on PC1.\n")
} else {
  cat(sprintf("\n  IC fraction weakly correlates with PC1 (r=%.3f).\n",
              cor_IC_PC1))
  cat("  PC1 is not primarily driven by immune cell composition.\n")
}
sink()
cat("  Saved: epidish_summary.txt\n")

cat("\n=== Done ===\n")
cat("Key outputs:\n")
cat("  epidish_barplot.png         — stacked cell fractions per sample\n")
cat("  epidish_prepost.png         — IC fraction Pre vs Post\n")
cat("  epidish_pc1_correlation.png — IC vs PC1\n")
cat("  limma_IC_adjusted.csv       — DMPs with IC as covariate\n")
cat("  epidish_summary.txt         — text summary\n")
