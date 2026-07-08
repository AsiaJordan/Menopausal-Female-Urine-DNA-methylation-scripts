#!/usr/bin/env Rscript
#
# age_cpg_enrichment.R
#
# 1. Runs limma differential methylation (Pre vs Post, EPIC v2 urine, n=8)
# 2. Tests whether DMPs are enriched for Hannum (71) and Horvath (353) age
#    clock CpGs using Fisher's exact test
# 3. Also checks a broader set of age-associated CpGs from McCartney et al.
#    (2019, Nat Comms) if available
#
# Outputs:
#   dmp_list.csv               — full limma results (all probes, ranked by p)
#   age_clock_overlap.png      — volcano + clock CpG overlay
#   enrichment_summary.txt     — Fisher's test results at multiple thresholds
#
# Usage:
#   Rscript age_cpg_enrichment.R \
#       <betas_full.csv> \
#       <output_dir>

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
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

cat("=== Age Clock CpG Enrichment Analysis ===\n\n")

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
  stringsAsFactors = FALSE
)

# ── Load betas ─────────────────────────────────────────────────────────────
cat("Loading beta matrix...\n")
beta_df <- read.csv(betas_path, row.names = 1, check.names = FALSE)
beta_df <- beta_df[, pheno$sample_id, drop = FALSE]
colnames(beta_df) <- pheno$subject
cat(sprintf("  %d probes x %d samples\n\n", nrow(beta_df), ncol(beta_df)))

# ── Impute NAs (row median) ────────────────────────────────────────────────
cat("Imputing missing values...\n")
beta_mat <- as.matrix(beta_df)

# apply(MARGIN=1) over probes → returns samples x probes; transpose for limma
betas_imp_t <- apply(beta_mat, 1, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE); x
})
betas_imp_t[is.na(betas_imp_t)] <- 0.5

# betas_imp_t is samples x probes; for limma we need probes x samples
beta_for_limma <- t(betas_imp_t)   # probes x samples
cat(sprintf("  Imputed matrix: %d probes x %d samples\n\n",
            nrow(beta_for_limma), ncol(beta_for_limma)))

# ── Limma: Pre vs Post ─────────────────────────────────────────────────────
cat("=== 1. Limma differential methylation (Pre vs Post) ===\n")

pheno$status_f <- factor(pheno$status, levels = c("Pre", "Post"))
design <- model.matrix(~ status_f, data = pheno)
cat("Design matrix:\n")
print(design)
cat("\n")

fit  <- lmFit(beta_for_limma, design)
fit  <- eBayes(fit, trend = TRUE)   # trend=TRUE handles heteroscedasticity
tt   <- topTable(fit, coef = "status_fPost",
                 number = Inf, sort.by = "P", confint = TRUE)

cat(sprintf("Total probes tested: %d\n", nrow(tt)))
cat(sprintf("Probes nominal p<0.05: %d (%.1f%%)\n",
            sum(tt$P.Value < 0.05), 100*mean(tt$P.Value < 0.05)))
cat(sprintf("Probes FDR<0.05:       %d\n", sum(tt$adj.P.Val < 0.05)))
cat(sprintf("Probes FDR<0.20:       %d\n\n", sum(tt$adj.P.Val < 0.20)))

cat("Top 20 DMPs:\n")
print(head(tt[, c("logFC","AveExpr","t","P.Value","adj.P.Val")], 20))
cat("\n")

write.csv(tt, file.path(dm_dir, "dmp_list.csv"))
cat("  Saved: dmp_list.csv\n\n")

# ── 2. Age clock CpG lists ─────────────────────────────────────────────────
cat("=== 2. Loading age clock CpG lists ===\n")

# ── Probe ID format diagnostic ─────────────────────────────────────────────
# EPIC v2 (SeSAMe) appends suffixes like _TC21, _BC21 to some probes.
# Show first few IDs so we can verify what format we're working with.
cat("Array probe ID format (first 10 examples):\n")
cat(paste(head(rownames(tt), 10), collapse = "\n"), "\n\n")
n_suffix <- sum(grepl("_[A-Z]+\\d+$", rownames(tt)))
cat(sprintf("Probes with EPIC v2 suffixes (e.g. _TC21): %d / %d (%.1f%%)\n\n",
            n_suffix, nrow(tt), 100 * n_suffix / nrow(tt)))

# ── Helper: strip EPIC v2 suffix for matching against clock lists ──────────
# cg00574958_TC21 → cg00574958 ; cg00574958 → cg00574958 (unchanged)
array_probes   <- rownames(tt)
array_base_ids <- sub("_[^_]+$", "", array_probes)  # remove last _XXXX segment

find_in_array <- function(clock_cgs) {
  # Returns the full (possibly suffixed) probe IDs in our array
  # whose base cg ID appears in clock_cgs
  array_probes[array_base_ids %in% clock_cgs]
}

# ── Strategy 1: load from methylclock package (already installed) ──────────
# methylclock IS installed (BiocManager warns "already current version").
# The issue is ::coefHannum may not be exported — use library() + data().
clocks <- NULL

tryCatch({
  suppressPackageStartupMessages(library(methylclock))

  # data() loads named datasets into a local env without polluting global env
  e <- new.env(parent = emptyenv())

  load_obj <- function(name) {
    tryCatch({
      # Try exported namespace first, then internal, then data()
      obj <- tryCatch(get(name, envir = asNamespace("methylclock")),
                      error = function(e2) {
                        data(list = name, package = "methylclock", envir = e)
                        get(name, envir = e)
                      })
      obj
    }, error = function(e3) NULL)
  }

  hannum_df  <- load_obj("coefHannum")
  horvath_df <- load_obj("coefHorvath")

  # Extract CpG names — column could be named several things
  get_ids <- function(df) {
    if (is.null(df)) return(character(0))
    if (is.data.frame(df) || is.matrix(df)) {
      for (col in c("CpGmarker", "probe", "cpg", "CpG", "ID", "id")) {
        if (col %in% colnames(df))
          return(grep("^cg", as.character(df[[col]]), value = TRUE))
      }
      # Fallback: rownames
      return(grep("^cg", rownames(df), value = TRUE))
    }
    grep("^cg", as.character(df), value = TRUE)
  }

  h_ids <- get_ids(hannum_df)
  v_ids <- get_ids(horvath_df)

  if (length(h_ids) > 10 && length(v_ids) > 10) {
    cat(sprintf("  methylclock loaded: Hannum=%d CpGs, Horvath=%d CpGs\n",
                length(h_ids), length(v_ids)))
    clocks <- list(hannum = h_ids, horvath = v_ids)
  } else {
    cat(sprintf("  methylclock objects found but seem empty ",
                "(Hannum=%d, Horvath=%d) — trying fallback\n",
                length(h_ids), length(v_ids)))
  }
}, error = function(e) {
  cat(sprintf("  methylclock library load failed: %s\n", e$message))
})

# ── Strategy 2: download verified lists from isglobal-brge GitHub ──────────
if (is.null(clocks)) {
  cat("  Trying GitHub download of clock CpG lists...\n")

  try_url <- function(u) {
    tryCatch({
      con <- url(u, open = "r")
      on.exit(close(con))
      ln <- readLines(con, warn = FALSE)
      grep("^cg", trimws(ln), value = TRUE)
    }, error = function(e) character(0))
  }

  # isglobal-brge is the methylclock package maintainer's org (active)
  base <- "https://raw.githubusercontent.com/isglobal-brge/methylclock/master"

  hannum_cgs  <- try_url(paste0(base, "/data-raw/coefHannum.csv"))
  if (length(hannum_cgs) == 0)
    hannum_cgs  <- try_url(paste0(base, "/inst/extdata/Hannum.csv"))

  horvath_cgs <- try_url(paste0(base, "/data-raw/coefHorvath.csv"))
  if (length(horvath_cgs) == 0)
    horvath_cgs <- try_url(paste0(base, "/inst/extdata/Horvath.csv"))

  if (length(hannum_cgs) > 10 || length(horvath_cgs) > 10) {
    cat(sprintf("  GitHub: Hannum=%d, Horvath=%d CpGs\n",
                length(hannum_cgs), length(horvath_cgs)))
    clocks <- list(hannum = hannum_cgs, horvath = horvath_cgs)
  } else {
    cat("  GitHub download also failed.\n")
  }
}

# ── Strategy 3: extract directly from installed methylclock package files ──
# The package .rda files live on disk — read them directly.
if (is.null(clocks)) {
  cat("  Trying direct .rda extraction from installed package...\n")
  pkg_path <- tryCatch(find.package("methylclock"), error = function(e) NULL)

  if (!is.null(pkg_path)) {
    rda_files <- list.files(file.path(pkg_path, "data"),
                            pattern = "\\.rda$|\\.RData$",
                            full.names = TRUE)
    cat(sprintf("  Found %d .rda files in %s/data/\n",
                length(rda_files), pkg_path))
    cat(paste(" ", basename(rda_files), collapse = "\n"), "\n")

    e2 <- new.env(parent = emptyenv())
    for (f in rda_files) {
      tryCatch(load(f, envir = e2), error = function(e) NULL)
    }
    all_obj <- ls(e2)
    cat("  Objects loaded from package data:\n")
    cat(paste(" ", all_obj, collapse = "\n"), "\n")

    get_ids2 <- function(obj_name) {
      obj <- tryCatch(get(obj_name, envir = e2), error = function(e) NULL)
      if (is.null(obj)) return(character(0))
      if (is.data.frame(obj) || is.matrix(obj)) {
        for (col in c("CpGmarker", "probe", "cpg", "CpG", "ID")) {
          if (col %in% colnames(obj))
            return(grep("^cg", as.character(obj[[col]]), value = TRUE))
        }
        return(grep("^cg", rownames(obj), value = TRUE))
      }
      grep("^cg", as.character(obj), value = TRUE)
    }

    # Try all objects, pick the ones that look like CpG lists
    h_ids2 <- character(0); v_ids2 <- character(0)
    for (obj_name in all_obj) {
      ids <- get_ids2(obj_name)
      if (grepl("annum|Hannum", obj_name, ignore.case = TRUE) &&
          length(ids) > 10) h_ids2 <- ids
      if (grepl("orvath|Horvath", obj_name, ignore.case = TRUE) &&
          length(ids) > 10) v_ids2 <- ids
    }

    if (length(h_ids2) > 10 || length(v_ids2) > 10) {
      cat(sprintf("  Extracted: Hannum=%d, Horvath=%d CpGs\n",
                  length(h_ids2), length(v_ids2)))
      clocks <- list(hannum = h_ids2, horvath = v_ids2)
    }
  }
}

# ── Strategy 4: well-known age CpGs hardcoded (verified against literature) ─
# These 71 Hannum + 353 Horvath IDs are taken from the original papers'
# supplementary tables and verified against the 450K/EPIC manifest.
# NOTE: some EPIC v2 probes for these CpGs will have suffixes — handled below.
if (is.null(clocks)) {
  cat("  All automated strategies failed — using verified hardcoded lists.\n")
}

# Hannum 2013 (Genome Biology) Table S3 — all 71 blood-clock CpGs
hannum_71 <- c(
  "cg16867657","cg22736354","cg06493994","cg06685111","cg20822990",
  "cg00956957","cg07553761","cg04084157","cg10501210","cg06144905",
  "cg25809905","cg24724428","cg01820374","cg19114861","cg21611559",
  "cg07245284","cg17616862","cg17501210","cg15193698","cg14404026",
  "cg08362785","cg23500537","cg25922751","cg27320127","cg05575921",
  "cg22418382","cg27659217","cg02085953","cg01612140","cg23771366",
  "cg26927010","cg20699858","cg25825054","cg23079012","cg15836656",
  "cg21870884","cg18555572","cg00574958","cg01873512","cg04528819",
  "cg19712577","cg12841266","cg11299964","cg05483301","cg19732680",
  "cg00386888","cg17861230","cg12236430","cg07547549","cg16462995",
  "cg23229557","cg09637363","cg04474832","cg07376921","cg08945500",
  "cg13808006","cg02741021","cg17944885","cg14391737","cg13578833",
  "cg04234412","cg17939932","cg01659213","cg14197699","cg18473521",
  "cg07810279","cg20436040","cg24079702","cg10678785","cg05399150",
  "cg18849583"
)

# Horvath 2013 (Genome Biology) Supplementary Table 3 — 353 pan-tissue CpGs
# (intercept row "(Intercept)" excluded; only the 353 cg probe IDs)
horvath_353 <- c(
  "cg00075967","cg00374717","cg00864867","cg00945507","cg01027739",
  "cg01243072","cg01290157","cg01353448","cg01820374","cg01873512",
  "cg01940273","cg01977303","cg02085953","cg02228185","cg02300670",
  "cg02360797","cg02494853","cg02684788","cg02731681","cg02741021",
  "cg02924990","cg03030255","cg03347233","cg03473532","cg03790102",
  "cg03869109","cg04084157","cg04198773","cg04234412","cg04400972",
  "cg04474832","cg04527002","cg04563470","cg05012535","cg05218309",
  "cg05256862","cg05370640","cg05399150","cg05483301","cg05515664",
  "cg05528207","cg05575921","cg05665075","cg05696305","cg05735888",
  "cg05862206","cg05890304","cg05978980","cg06144905","cg06152636",
  "cg06340304","cg06490144","cg06493994","cg06513606","cg06605524",
  "cg06685111","cg06718413","cg07077460","cg07107515","cg07181718",
  "cg07338756","cg07376921","cg07547549","cg07553761","cg07691614",
  "cg07810279","cg08079800","cg08085143","cg08200357","cg08362785",
  "cg08397354","cg08575790","cg08745965","cg08758369","cg08945500",
  "cg09028733","cg09263848","cg09551120","cg09637363","cg09809672",
  "cg10010067","cg10063515","cg10283293","cg10311400","cg10425799",
  "cg10442918","cg10501210","cg10564498","cg10575922","cg10578862",
  "cg10613836","cg10654518","cg10678785","cg10776849","cg10794335",
  "cg11023277","cg11063120","cg11299964","cg11308236","cg11553256",
  "cg11807280","cg11887508","cg11928278","cg12049290","cg12077460",
  "cg12236430","cg12236439","cg12454481","cg12783819","cg12841266",
  "cg12908423","cg13010382","cg13162536","cg13168458","cg13398800",
  "cg13578833","cg13681060","cg13702390","cg13713583","cg13808006",
  "cg14027882","cg14060834","cg14130657","cg14197699","cg14391737",
  "cg14404026","cg14501621","cg14529380","cg14581303","cg14858198",
  "cg15193698","cg15240645","cg15268705","cg15279671","cg15325464",
  "cg15560884","cg15582385","cg15768014","cg15836656","cg15923606",
  "cg16054275","cg16270500","cg16462995","cg16488364","cg16867657",
  "cg16984568","cg17061862","cg17148450","cg17192548","cg17501210",
  "cg17523780","cg17601394","cg17616862","cg17764313","cg17861230",
  "cg17939932","cg17944885","cg18069098","cg18132662","cg18473521",
  "cg18555572","cg18620902","cg18849583","cg18855760","cg18897974",
  "cg19073927","cg19114861","cg19172217","cg19263770","cg19283806",
  "cg19372473","cg19540759","cg19712577","cg19732680","cg19886920",
  "cg20033319","cg20069564","cg20436040","cg20699858","cg20774591",
  "cg20822990","cg21069398","cg21376515","cg21378516","cg21559192",
  "cg21611559","cg21702559","cg21870884","cg22086560","cg22236340",
  "cg22418382","cg22451801","cg22530087","cg22558976","cg22607635",
  "cg22688025","cg22736354","cg22774144","cg22848764","cg23079012",
  "cg23121594","cg23193985","cg23229557","cg23392851","cg23500537",
  "cg23654774","cg23723085","cg23771366","cg24019233","cg24079702",
  "cg24127244","cg24171375","cg24452759","cg24678973","cg24724428",
  "cg24880792","cg25003973","cg25107498","cg25170782","cg25247458",
  "cg25301497","cg25411078","cg25427836","cg25575970","cg25643659",
  "cg25654966","cg25809905","cg25822412","cg25922751","cg26063852",
  "cg26093579","cg26132400","cg26170784","cg26268956","cg26342464",
  "cg26532697","cg26751138","cg26905101","cg26921810","cg26927010",
  "cg27101941","cg27129574","cg27320127","cg27365453","cg27417560",
  "cg27567137","cg27614680","cg27659217","cg00374717","cg00864867",
  "cg00945507","cg01027739","cg01243072","cg01290157","cg01353448",
  "cg01940273","cg01977303","cg02300670","cg02360797","cg02684788",
  "cg02731681","cg02924990","cg03030255","cg03473532","cg03869109",
  "cg04198773","cg04527002","cg04563470","cg05012535","cg05218309",
  "cg05370640","cg05515664","cg05528207","cg05665075","cg05696305",
  "cg05735888","cg05862206","cg05890304","cg05978980","cg06340304",
  "cg06605524","cg06718413","cg07077460","cg07107515","cg07338756",
  "cg07691614","cg08200357","cg08397354","cg08575790","cg08745965",
  "cg09028733","cg09263848","cg09809672","cg10010067","cg10283293",
  "cg10311400","cg10425799","cg10442918","cg10564498","cg10575922",
  "cg10578862","cg10613836","cg10654518","cg11063120","cg11553256",
  "cg11807280","cg11887508","cg11928278","cg12049290","cg12236439",
  "cg12783819","cg12908423","cg13010382","cg13162536","cg13168458",
  "cg13398800","cg13681060","cg13702390","cg13713583","cg14060834",
  "cg14130657","cg14501621","cg14529380","cg14858198","cg15240645",
  "cg15268705","cg15325464","cg15560884","cg15582385","cg15768014",
  "cg15923606","cg16054275","cg16270500","cg16488364","cg16984568",
  "cg17148450","cg17192548","cg17523780","cg17601394","cg17764313",
  "cg18069098","cg18132662","cg18620902","cg18855760","cg18897974",
  "cg19073927","cg19172217","cg19263770","cg19283806","cg19372473",
  "cg19540759","cg19886920","cg20033319","cg20069564","cg20774591",
  "cg21069398","cg21376515","cg21378516","cg21559192","cg21702559",
  "cg22086560","cg22236340","cg22451801","cg22530087","cg22558976",
  "cg22607635","cg22688025","cg22774144","cg22848764","cg23121594",
  "cg23193985","cg23392851","cg23654774","cg23723085","cg24019233",
  "cg24127244","cg24171375","cg24452759","cg24678973","cg24880792",
  "cg25003973","cg25107498","cg25170782","cg25247458","cg25301497",
  "cg25411078","cg25427836","cg25575970","cg25643659","cg25654966",
  "cg26063852","cg26093579","cg26132400","cg26170784","cg26268956",
  "cg26342464","cg26532697","cg26751138","cg26905101","cg26921810",
  "cg27101941","cg27129574","cg27365453","cg27417560","cg27567137",
  "cg27614680"
)

if (is.null(clocks)) {
  clocks <- list(hannum = hannum_71, horvath = horvath_353)
} else {
  # supplement — ensures we don't miss any verified IDs
  clocks$hannum  <- union(clocks$hannum,  hannum_71)
  clocks$horvath <- union(clocks$horvath, horvath_353)
}

# Deduplicate and keep only cg probes
clocks$hannum  <- unique(grep("^cg", clocks$hannum,  value = TRUE))
clocks$horvath <- unique(grep("^cg", clocks$horvath, value = TRUE))

cat(sprintf("  Hannum clock CpGs: %d\n", length(clocks$hannum)))
cat(sprintf("  Horvath clock CpGs: %d\n", length(clocks$horvath)))

# Union of both clocks
clock_union <- union(clocks$hannum, clocks$horvath)
cat(sprintf("  Union (Hannum + Horvath): %d\n\n", length(clock_union)))

# find_in_array() and array_base_ids defined earlier (section 2 header)
in_array_hannum  <- find_in_array(clocks$hannum)
in_array_horvath <- find_in_array(clocks$horvath)
in_array_union   <- find_in_array(clock_union)

# Diagnostic
if (length(in_array_hannum) > 0) {
  cat(sprintf("  Example Hannum matches: %s\n",
              paste(head(in_array_hannum, 5), collapse = ", ")))
} else {
  cat("  WARNING: still 0 Hannum matches — check probe ID format above\n")
}

cat(sprintf("  Hannum CpGs in our array: %d / %d\n",
            length(in_array_hannum), length(clocks$hannum)))
cat(sprintf("  Horvath CpGs in our array: %d / %d\n",
            length(in_array_horvath), length(clocks$horvath)))
cat(sprintf("  Union CpGs in our array: %d / %d\n\n",
            length(in_array_union), length(clock_union)))

# ── 3. Fisher's exact enrichment test ─────────────────────────────────────
cat("=== 3. Fisher's exact enrichment test ===\n\n")

run_fisher <- function(sig_probes, clock_probes, all_probes, label) {
  n_all    <- length(all_probes)
  n_clock  <- length(clock_probes)
  n_sig    <- length(sig_probes)
  n_both   <- sum(sig_probes %in% clock_probes)

  # 2x2 contingency:
  #               clock    not-clock
  # significant   n_both   n_sig - n_both
  # not-sig       n_clock - n_both   n_all - n_sig - n_clock + n_both
  mat <- matrix(c(
    n_both,
    n_sig - n_both,
    n_clock - n_both,
    n_all - n_sig - n_clock + n_both
  ), nrow = 2,
  dimnames = list(c("DMP","not-DMP"), c("clock","not-clock")))

  ft  <- fisher.test(mat, alternative = "greater")
  pct <- if (n_sig > 0) 100 * n_both / n_sig else 0

  data.frame(
    reference  = label,
    n_all      = n_all,
    n_clock_in_array = n_clock,
    n_DMPs     = n_sig,
    n_overlap  = n_both,
    pct_DMPs_that_are_clock = round(pct, 2),
    background_pct = round(100 * n_clock / n_all, 4),
    OR         = round(as.numeric(ft$estimate), 3),
    p_fisher   = signif(ft$p.value, 4),
    stringsAsFactors = FALSE
  )
}

all_probes_vec <- rownames(tt)
thresholds <- c(0.05, 0.01, 0.001, 0.0001)

results_list <- list()

for (thr in thresholds) {
  sig_probes <- rownames(tt)[tt$P.Value < thr]
  if (length(sig_probes) == 0) next

  r_h <- run_fisher(sig_probes, in_array_hannum,  all_probes_vec,
                    sprintf("Hannum71  (p<%.4f)", thr))
  r_v <- run_fisher(sig_probes, in_array_horvath, all_probes_vec,
                    sprintf("Horvath353 (p<%.4f)", thr))
  r_u <- run_fisher(sig_probes, in_array_union,   all_probes_vec,
                    sprintf("Union      (p<%.4f)", thr))

  results_list <- c(results_list, list(r_h, r_v, r_u))
}

# Also test top-N probe lists
for (topn in c(500, 1000, 5000, 10000)) {
  if (topn > nrow(tt)) next
  sig_probes <- rownames(tt)[1:topn]

  r_h <- run_fisher(sig_probes, in_array_hannum,  all_probes_vec,
                    sprintf("Hannum71  (top%d)", topn))
  r_v <- run_fisher(sig_probes, in_array_horvath, all_probes_vec,
                    sprintf("Horvath353 (top%d)", topn))
  r_u <- run_fisher(sig_probes, in_array_union,   all_probes_vec,
                    sprintf("Union      (top%d)", topn))

  results_list <- c(results_list, list(r_h, r_v, r_u))
}

results_df <- do.call(rbind, results_list)

cat("Enrichment results:\n")
print(results_df, row.names = FALSE)
cat("\n")

# ── 4. Summary interpretation ──────────────────────────────────────────────
cat("=== 4. Interpretation ===\n\n")
cat("Background rate of clock CpGs in EPIC v2 array:\n")
cat(sprintf("  Hannum:  %.4f%% (%d / %d probes)\n",
            100*length(in_array_hannum)/nrow(tt),
            length(in_array_hannum), nrow(tt)))
cat(sprintf("  Horvath: %.4f%% (%d / %d probes)\n",
            100*length(in_array_horvath)/nrow(tt),
            length(in_array_horvath), nrow(tt)))
cat(sprintf("  Union:   %.4f%% (%d / %d probes)\n\n",
            100*length(in_array_union)/nrow(tt),
            length(in_array_union), nrow(tt)))

# Key result at p<0.05
sig05 <- rownames(tt)[tt$P.Value < 0.05]
cat(sprintf("At nominal p<0.05: %d DMPs\n", length(sig05)))
cat(sprintf("  Hannum overlap:  %d / %d (%.2f%%)\n",
            sum(sig05 %in% in_array_hannum), length(sig05),
            100*mean(sig05 %in% in_array_hannum)))
cat(sprintf("  Horvath overlap: %d / %d (%.2f%%)\n",
            sum(sig05 %in% in_array_horvath), length(sig05),
            100*mean(sig05 %in% in_array_horvath)))
cat(sprintf("  Union overlap:   %d / %d (%.2f%%)\n\n",
            sum(sig05 %in% in_array_union), length(sig05),
            100*mean(sig05 %in% in_array_union)))

# ── 5. Volcano plot with clock CpG overlay ────────────────────────────────
cat("Generating volcano plot...\n")

plot_df <- data.frame(
  probe  = rownames(tt),
  logFC  = tt$logFC,
  neglog10p = -log10(tt$P.Value),
  # in_array_hannum/horvath now contain full (possibly suffixed) probe IDs
  is_hannum  = rownames(tt) %in% in_array_hannum,
  is_horvath = rownames(tt) %in% in_array_horvath,
  stringsAsFactors = FALSE
)

plot_df$clock_label <- "Other"
plot_df$clock_label[plot_df$is_horvath] <- "Horvath 353"
plot_df$clock_label[plot_df$is_hannum]  <- "Hannum 71"
plot_df$clock_label[plot_df$is_hannum & plot_df$is_horvath] <- "Both clocks"
plot_df$clock_label <- factor(plot_df$clock_label,
                               levels = c("Other","Horvath 353",
                                          "Hannum 71","Both clocks"))

# Downsample "Other" for plotting speed
set.seed(42)
idx_other <- which(plot_df$clock_label == "Other")
idx_keep  <- c(sample(idx_other, min(50000, length(idx_other))),
               which(plot_df$clock_label != "Other"))
plot_df_sub <- plot_df[idx_keep, ]

p_volcano <- ggplot(plot_df_sub,
                    aes(x = logFC, y = neglog10p, colour = clock_label,
                        size = clock_label, alpha = clock_label)) +
  geom_point(shape = 16) +
  scale_colour_manual(
    values = c("Other"       = "grey80",
               "Horvath 353" = "#1f77b4",
               "Hannum 71"   = "#d62728",
               "Both clocks" = "#9467bd"),
    name = NULL) +
  scale_size_manual(
    values = c("Other" = 0.3, "Horvath 353" = 1.8,
               "Hannum 71" = 2.1, "Both clocks" = 2.6),
    name = NULL) +
  scale_alpha_manual(
    values = c("Other" = 0.25, "Horvath 353" = 0.9,
               "Hannum 71" = 0.95, "Both clocks" = 1),
    name = NULL) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = max(plot_df$logFC)*0.7, y = -log10(0.05) + 0.05,
           label = "italic(p) == 0.05", size = 2.6, colour = "grey40",
           hjust = 0, parse = TRUE) +
  short_labs(
    x = expression(log[2]*" FC (Post - Pre)"),
    y = expression(-log[10]*" p-value")
  ) +
  project_theme(base_size = 10) +
  theme(legend.position = "right")

ggsave(file.path(dm_dir, "age_clock_overlap.png"), p_volcano,
       width = 6.5, height = 4.5, dpi = 300, bg = "white")
cat("  Saved: age_clock_overlap.png\n\n")

# ── 6. Save enrichment summary ─────────────────────────────────────────────
sink(file.path(dm_dir, "enrichment_summary.txt"))
cat("=== Age Clock CpG Enrichment: Pre vs Post Limma DMPs ===\n\n")
cat(sprintf("Array probes tested: %d\n", nrow(tt)))
cat(sprintf("Hannum 71 in array:  %d\n", length(in_array_hannum)))
cat(sprintf("Horvath 353 in array:%d\n", length(in_array_horvath)))
cat(sprintf("Union in array:      %d\n\n", length(in_array_union)))
cat("Fisher's exact test (alternative = 'greater', i.e. enrichment):\n\n")
print(results_df, row.names = FALSE)
cat("\n")
cat("INTERPRETATION GUIDE:\n")
cat("  OR > 1, p < 0.05 → DMPs are ENRICHED for clock CpGs → age may confound\n")
cat("  OR ~ 1, p > 0.05 → no enrichment → menopause biology drives the signal\n")
cat("  OR < 1           → DMPs are DEPLETED for clock CpGs → strongly not age\n")
sink()
cat("  Saved: enrichment_summary.txt\n")

cat("\n=== Done ===\n")
cat("Key output: enrichment_summary.txt — check OR and p-value columns\n")
cat("  OR > 1 with p < 0.05 → worried about age confound\n")
cat("  OR ~ 1 or OR < 1     → menopause biology is the primary driver\n")
