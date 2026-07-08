#!/usr/bin/env Rscript
#
# batch_check_menopause.R
#
# Assesses whether chip (batch) or menopausal status (biology) drives
# the primary variance structure in the urine EPIC v2 data.
#
# Key diagnostic: FMLS_009 is the only premenopausal sample on chip 1
# (with 3 postmenopausal). If chip dominates, it clusters with post.
# If biology dominates, it clusters with the 4 pre samples on chip 2.
#
# Outputs:
#   pca_batch_biology.png   — PCA coloured by chip AND status
#   hclust_heatmap.png      — hierarchical clustering of samples
#   variance_explained.csv  — PC variances
#   batch_diagnosis.txt     — plain-text verdict

suppressPackageStartupMessages({
  library(ggplot2)
  # ComplexHeatmap gives us per-cell text colouring (white on dark,
  # black on light) which pheatmap cannot do natively.
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
    BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
  if (!requireNamespace("circlize", quietly = TRUE))
    install.packages("circlize", repos = "https://cloud.r-project.org")
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

# Shared project theme + palettes. Resolve script dir whether run via
# Rscript or source(), falling back to CWD.
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

# ── Paths ──────────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
betas_path  <- if (length(args) >= 1) args[1] else
  "menopause_urine_out/betas_full.csv"
out_dir     <- if (length(args) >= 2) args[2] else
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

cat("=== Batch vs Biology Diagnostic ===\n\n")

# ── Phenotype table (from sample sheet) ───────────────────────────────────
pheno <- data.frame(
  sample_id  = c("207758450138_R03C01","207758450138_R04C01",
                 "207758450138_R07C01","207758450138_R08C01",
                 "208356980025_R02C01","208356980025_R03C01",
                 "208356980025_R05C01","208356980025_R06C01"),
  subject    = c("FMLS_007","FMLS_009","FMLS_005","FMLS_006",
                 "FMLS_011","FMLS_015","FMLS_002","FMLS_010"),
  status     = c("Post","Pre","Post","Post",
                 "Pre","Pre","Pre","Pre"),
  chip       = c("207758450138","207758450138","207758450138","207758450138",
                 "208356980025","208356980025","208356980025","208356980025"),
  age        = c(61, 31, 62, 66, 45, 23, 27, 41),
  flag       = c("","OC_pill","","ovary_removed",
                 "PREGNANT","OC_pill","",""),
  stringsAsFactors = FALSE
)

cat("Sample phenotypes:\n")
print(pheno[, c("subject","status","chip","age","flag")])
cat("\n")

# ── Load betas ─────────────────────────────────────────────────────────────
cat("Loading beta matrix...\n")
beta_df <- read.csv(betas_path, row.names = 1, check.names = FALSE)
cat(sprintf("  %d probes × %d samples\n\n", nrow(beta_df), ncol(beta_df)))

# Align columns to pheno order
beta_df <- beta_df[, pheno$sample_id, drop = FALSE]
colnames(beta_df) <- pheno$subject   # rename to FMLS IDs for readability

# ── Select high-variance probes ────────────────────────────────────────────
cat("Selecting top variable probes...\n")
row_vars  <- apply(beta_df, 1, var, na.rm = TRUE)
top_n     <- 10000
top_probes <- order(row_vars, decreasing = TRUE)[1:min(top_n, sum(!is.na(row_vars)))]
betas_hv   <- beta_df[top_probes, ]
cat(sprintf("  Using top %d variable probes\n\n", nrow(betas_hv)))

# ── PCA ────────────────────────────────────────────────────────────────────
cat("Running PCA...\n")
# Impute NAs with row median before PCA
betas_imp <- apply(betas_hv, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
})  # returns samples × probes
betas_imp[is.na(betas_imp)] <- 0.5  # fallback

pca_res <- prcomp(betas_imp, center = TRUE, scale. = FALSE)
pca_df  <- as.data.frame(pca_res$x[, 1:3])
pca_df$subject <- rownames(pca_df)
pca_df  <- merge(pca_df, pheno, by = "subject")

var_exp <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 2)
cat(sprintf("  PC1: %.1f%%  PC2: %.1f%%  PC3: %.1f%%\n\n",
            var_exp[1], var_exp[2], var_exp[3]))

write.csv(data.frame(PC = paste0("PC", seq_along(var_exp)),
                     variance_pct = var_exp),
          file.path(qc_dir, "variance_explained.csv"), row.names = FALSE)

# PCA plot: shape = chip, colour = status, label = subject.
# All diagnostic narrative lives in the figure caption / response letter —
# the plot itself only carries variance-explained in the axis labels.
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2,
                             colour = status, shape = chip,
                             label = subject)) +
  geom_point(size = 3.2, stroke = 0.9) +
  geom_text(vjust = -1.0, size = 2.8, show.legend = FALSE) +
  scale_colour_manual(values = STATUS_COLORS,
                      name = "Status") +
  scale_shape_manual(values = CHIP_SHAPES,
                     name = "Chip",
                     labels = c("207758450138" = "...450138",
                                "208356980025" = "...980025")) +
  # Expand axes so FMLS_005 / FMLS_009 labels at the corners don't clip.
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.12))) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  short_labs(
    x = sprintf("PC1 (%.1f%% variance)", var_exp[1]),
    y = sprintf("PC2 (%.1f%% variance)", var_exp[2])
  ) +
  project_theme(base_size = 10)

ggsave(file.path(qc_dir, "pca_batch_biology.png"), p_pca,
       width = 5.8, height = 4.2, dpi = 300, bg = "white")
cat("  Saved: pca_batch_biology.png\n")

# ── Hierarchical clustering heatmap ───────────────────────────────────────
cat("Running hierarchical clustering...\n")

# Sample distance matrix on top variable probes
dist_mat <- dist(t(as.matrix(betas_imp)), method = "euclidean")
hc       <- hclust(dist_mat, method = "complete")

# Annotation for heatmap: status, chip, plus the three hormonal
# confounders (OC, pregnancy, oophorectomy) and age as a gradient.
ann_col <- data.frame(
  Status       = pheno$status,
  Chip         = pheno$chip,
  Age          = pheno$age,
  OC_pill      = ifelse(pheno$flag == "OC_pill",      "Yes", "No"),
  Pregnant     = ifelse(pheno$flag == "PREGNANT",     "Yes", "No"),
  Oophorectomy = ifelse(pheno$flag == "ovary_removed","Yes", "No"),
  row.names    = pheno$subject,
  check.names  = FALSE
)
ann_colours <- list(
  Status       = c(Post = "#d62728", Pre = "#1f77b4"),
  Chip         = c("207758450138" = "#ff7f0e",
                   "208356980025" = "#2ca02c"),
  # ComplexHeatmap continuous annotation expects a colorRamp2 function
  Age          = colorRamp2(c(20, 70), c("#f2f0f7", "#54278f")),
  OC_pill      = c(No = "grey92", Yes = "#6a4b9b"),
  Pregnant     = c(No = "grey92", Yes = "#e67e22"),
  Oophorectomy = c(No = "grey92", Yes = "#2c3e50")
)

# Correlation matrix of samples
cor_mat <- cor(t(as.matrix(betas_imp)), use = "pairwise.complete.obs")
rownames(cor_mat) <- colnames(cor_mat) <- pheno$subject

# Standard divergent palette over the full correlation range.
hm_col <- colorRamp2(c(-1, 0, 1), c("#4575b4", "white", "#d73027"))

# Switch text colour by cell darkness — black on pale cells, white on
# strongly coloured ones.
cell_text_fun <- function(j, i, x, y, w, h, fill) {
  val <- cor_mat[i, j]
  txt_col <- if (abs(val) > 0.4) "white" else "#111111"
  grid.text(sprintf("%.3f", val), x, y,
            gp = gpar(fontsize = 8, col = txt_col))
}

# Top annotation
top_ann <- HeatmapAnnotation(
  df = ann_col,
  col = ann_colours,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontsize = 9),
  gap = unit(0.6, "mm"),
  simple_anno_size = unit(4, "mm"),
  show_annotation_name = TRUE
)

hm <- Heatmap(
  cor_mat,
  name = "correlation",
  col  = hm_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows = "complete",
  clustering_method_columns = "complete",
  top_annotation = top_ann,
  cell_fun = cell_text_fun,
  column_title = "Sample-to-sample methylation correlation (top 10K variable probes)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 9),
                               labels_gp = gpar(fontsize = 8))
)

png(file.path(qc_dir, "hclust_heatmap.png"),
    width = 1100, height = 950, res = 120)
draw(hm, heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)
dev.off()
cat("  Saved: hclust_heatmap.png\n\n")

# ── Diagnostic text ────────────────────────────────────────────────────────
# Where does FMLS_009 fall on PC1?
fmls009_pc1 <- pca_df$PC1[pca_df$subject == "FMLS_009"]
post_pc1    <- pca_df$PC1[pca_df$status  == "Post"]
pre_pc1     <- pca_df$PC1[pca_df$status  == "Pre" & pca_df$subject != "FMLS_009"]

dist_to_post <- abs(fmls009_pc1 - mean(post_pc1))
dist_to_pre  <- abs(fmls009_pc1 - mean(pre_pc1))

verdict <- if (dist_to_post < dist_to_pre) {
  "FMLS_009 is CLOSER to the Post group on PC1 → chip effect likely dominates"
} else {
  "FMLS_009 is CLOSER to the Pre group on PC1 → biology likely dominates"
}

diag_text <- c(
  "=== Batch vs Biology Diagnostic ===",
  "",
  "FMLS_009: Pre-menopausal, age 31, OC pill — on chip 207758450138 with 3 Post samples",
  "",
  sprintf("PC1 values:"),
  sprintf("  FMLS_009 (Pre, chip1):  %.4f", fmls009_pc1),
  sprintf("  Post mean (chip1):      %.4f", mean(post_pc1)),
  sprintf("  Pre mean  (chip2):      %.4f", mean(pre_pc1)),
  "",
  sprintf("Distance FMLS_009 → Post centroid: %.4f", dist_to_post),
  sprintf("Distance FMLS_009 → Pre centroid:  %.4f", dist_to_pre),
  "",
  paste("VERDICT:", verdict),
  "",
  "Additional confounds present:",
  "  FMLS_011: PREGNANT (chip 2, Pre group)",
  "  FMLS_009: OC pill  (chip 1, the diagnostic Pre sample)",
  "  FMLS_015: OC pill  (chip 2, Pre group)",
  "  FMLS_006: One ovary removed/clipped (chip 1, Post group)"
)

writeLines(diag_text, file.path(qc_dir, "batch_diagnosis.txt"))
cat(paste(diag_text, collapse = "\n"), "\n")
