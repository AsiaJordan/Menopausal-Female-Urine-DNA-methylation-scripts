#!/usr/bin/env Rscript
#
# map_clock_cpgs_to_epicv2.R
#
# Finds EPIC v2 equivalents for Hannum/Horvath clock CpGs that are absent
# from the EPIC v2 array by direct probe ID matching.
#
# Strategy:
#   1. Identify which clock CpGs are "missing" from EPIC v2 (no probe ID match)
#   2. Get their genomic coordinates from the 450K annotation (hg19)
#   3. Lift over hg19 → hg38 using UCSC chain file
#   4. Match against EPIC v2 probe positions (hg38) within 2 bp window
#      (CpG sites have a C and G; probe can target either strand)
#   5. Report EPIC v2 probe IDs covering the same CpG sites
#   6. Re-run Fisher enrichment with extended probe set
#
# Outputs:
#   clock_cpg_v2_mapping.csv      — mapping table (old ID → v2 equivalents)
#   enrichment_extended.txt       — Fisher results with extended probe set
#
# Usage:
#   Rscript map_clock_cpgs_to_epicv2.R \
#       <betas_full.csv> <output_dir>

suppressPackageStartupMessages(library(limma))

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

cat("=== Mapping clock CpGs to EPIC v2 equivalents ===\n\n")

# ── Install packages ───────────────────────────────────────────────────────
bioc_pkgs <- c(
  "IlluminaHumanMethylation450kanno.ilmn12.hg19",  # 450K positions (hg19)
  "IlluminaHumanMethylationEPICv2anno.20a1.hg38",  # EPIC v2 positions (hg38)
  "rtracklayer",                                     # liftOver
  "GenomicRanges"
)
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    cat(sprintf("Installing %s...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}
suppressPackageStartupMessages({
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  library(rtracklayer)
  library(GenomicRanges)
})

# ── Clock CpG lists (same as age_cpg_enrichment.R) ────────────────────────
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
  "cg27567137","cg27614680","cg27659217"
)

clock_union <- unique(c(hannum_71, horvath_353))
cat(sprintf("Clock CpGs: Hannum=%d, Horvath=%d, Union=%d\n\n",
            length(hannum_71), length(horvath_353), length(clock_union)))

# ── Load EPIC v2 array probes (from betas_full.csv header only) ────────────
cat("Reading EPIC v2 probe IDs from beta matrix...\n")
con <- file(betas_path, "r")
header <- readLines(con, n = 1)
close(con)
# probe IDs are row names (first column after header), read them separately
array_probes_raw <- read.csv(betas_path, nrows = 1, row.names = 1,
                             check.names = FALSE)
# Actually just read probe names column
probe_col <- read.csv(betas_path, colClasses = c(NA, rep("NULL", 8)),
                      check.names = FALSE)
array_probes   <- probe_col[, 1]
array_base_ids <- sub("_[^_]+$", "", array_probes)
cat(sprintf("  %d probes in array\n", length(array_probes)))
cat(sprintf("  %d unique base IDs\n\n", length(unique(array_base_ids))))

# Which clock CpGs are already found by direct base ID match?
found_direct   <- clock_union[clock_union %in% array_base_ids]
missing_clocks <- clock_union[!clock_union %in% array_base_ids]

cat(sprintf("Clock CpGs found by direct ID match: %d / %d\n",
            length(found_direct), length(clock_union)))
cat(sprintf("Clock CpGs MISSING from EPIC v2:     %d / %d\n\n",
            length(missing_clocks), length(clock_union)))

if (length(missing_clocks) == 0) {
  cat("All clock CpGs found — no mapping needed.\n")
  quit(save = "no")
}

# ── Get 450K positions (hg19) for missing clock CpGs ──────────────────────
cat("=== Step 1: 450K positions for missing clock CpGs (hg19) ===\n")

anno450 <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450_sub <- anno450[rownames(anno450) %in% missing_clocks,
                        c("chr","pos","strand")]
found_in_450k <- rownames(anno450_sub)
not_in_450k   <- missing_clocks[!missing_clocks %in% found_in_450k]

cat(sprintf("  Missing CpGs found in 450K annotation: %d\n",
            nrow(anno450_sub)))
if (length(not_in_450k) > 0) {
  cat(sprintf("  Not in 450K at all (truly absent): %d\n", length(not_in_450k)))
  cat(paste("   ", not_in_450k, collapse = "\n"), "\n")
}
cat("\n")

# ── Build GRanges for 450K positions ──────────────────────────────────────
gr_450k <- GRanges(
  seqnames = anno450_sub$chr,
  ranges   = IRanges(start = anno450_sub$pos, width = 1),
  probe_id = rownames(anno450_sub)
)

# ── Liftover hg19 → hg38 ──────────────────────────────────────────────────
cat("=== Step 2: LiftOver hg19 → hg38 ===\n")

chain_path <- file.path(out_dir, "hg19ToHg38.over.chain.gz")

if (!file.exists(chain_path)) {
  cat("  Downloading hg19ToHg38 chain file...\n")
  chain_url <- "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
  tryCatch(
    download.file(chain_url, chain_path, quiet = TRUE, method = "curl"),
    error = function(e) {
      download.file(chain_url, chain_path, quiet = TRUE, method = "wget")
    }
  )
}

chain  <- import.chain(chain_path)
gr_38  <- liftOver(gr_450k, chain)

# liftOver returns a GRangesList — keep only 1:1 mappings
mapped_idx  <- which(lengths(gr_38) == 1)
unmapped    <- gr_450k$probe_id[lengths(gr_38) != 1]
gr_38_flat  <- unlist(gr_38[mapped_idx])

cat(sprintf("  Mapped: %d / %d probes\n", length(gr_38_flat), nrow(anno450_sub)))
if (length(unmapped) > 0)
  cat(sprintf("  Unmapped (multi-map or dropout): %d — %s\n",
              length(unmapped), paste(unmapped, collapse = ", ")))
cat("\n")

# ── Get EPIC v2 probe positions (hg38) ────────────────────────────────────
cat("=== Step 3: EPIC v2 probe positions (hg38) ===\n")

anno_v2 <- getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
cat(sprintf("  EPIC v2 annotation: %d probes\n", nrow(anno_v2)))

gr_v2 <- GRanges(
  seqnames = anno_v2$chr,
  ranges   = IRanges(start = anno_v2$pos, width = 1),
  probe_id = rownames(anno_v2)
)

# ── Coordinate match: 450K (hg38) ↔ EPIC v2 (hg38) within 2 bp ───────────
cat("=== Step 4: Coordinate matching (window ±2 bp) ===\n\n")

# Expand 450K positions by ±2 bp to catch strand-shifted positions
gr_38_expanded <- gr_38_flat + 2

hits <- findOverlaps(gr_38_expanded, gr_v2)

# Build mapping table
map_df <- data.frame(
  clock_cpg      = gr_38_flat$probe_id[queryHits(hits)],
  v2_probe_id    = gr_v2$probe_id[subjectHits(hits)],
  v2_base_id     = sub("_[^_]+$", "", gr_v2$probe_id[subjectHits(hits)]),
  chr_hg38       = as.character(seqnames(gr_v2)[subjectHits(hits)]),
  pos_hg38       = start(gr_v2)[subjectHits(hits)],
  stringsAsFactors = FALSE
)

# Flag whether the v2 probe is the same cg ID (suffix) or a genuinely new ID
map_df$same_base_id <- map_df$v2_base_id == map_df$clock_cpg

# Exclude v2 probes that already matched by direct ID (suffix match already handled)
# Only report probes where the v2 base ID is DIFFERENT from the clock probe ID
map_df_new <- map_df[!map_df$same_base_id, ]
map_df_same <- map_df[map_df$same_base_id, ]

cat(sprintf("Coordinate matches found: %d total hits\n", nrow(map_df)))
cat(sprintf("  Same base cg ID (suffix variants already in array): %d\n",
            nrow(map_df_same)))
cat(sprintf("  Different cg ID (genuinely new v2 probe for same site): %d\n\n",
            nrow(map_df_new)))

if (nrow(map_df_new) > 0) {
  cat("New EPIC v2 equivalents for clock CpGs:\n")
  print(map_df_new[, c("clock_cpg","v2_probe_id","chr_hg38","pos_hg38")],
        row.names = FALSE)
  cat("\n")
} else {
  cat("No genuinely new EPIC v2 probe IDs found for missing clock CpGs.\n")
  cat("The missing CpGs are truly absent from EPIC v2 (probe retired, no replacement).\n\n")
}

# Summary by clock
missing_hannum  <- hannum_71[!hannum_71 %in% array_base_ids]
missing_horvath <- horvath_353[!horvath_353 %in% array_base_ids]

recovered_hannum  <- missing_hannum[missing_hannum %in% map_df_new$clock_cpg]
recovered_horvath <- missing_horvath[missing_horvath %in% map_df_new$clock_cpg]
truly_gone_hannum  <- missing_hannum[!missing_hannum %in% map_df_new$clock_cpg]
truly_gone_horvath <- missing_horvath[!missing_horvath %in% map_df_new$clock_cpg]

cat("=== Summary ===\n\n")
cat(sprintf("Hannum 71:\n"))
cat(sprintf("  Found by direct ID:   %d\n",
            sum(hannum_71 %in% array_base_ids)))
cat(sprintf("  Recovered via coords: %d\n", length(recovered_hannum)))
cat(sprintf("  Truly gone from v2:   %d\n\n", length(truly_gone_hannum)))

cat(sprintf("Horvath 353:\n"))
cat(sprintf("  Found by direct ID:   %d\n",
            sum(horvath_353 %in% array_base_ids)))
cat(sprintf("  Recovered via coords: %d\n", length(recovered_horvath)))
cat(sprintf("  Truly gone from v2:   %d\n\n", length(truly_gone_horvath)))

if (length(truly_gone_hannum) > 0) {
  cat("Hannum CpGs with no EPIC v2 equivalent:\n")
  cat(paste(" ", truly_gone_hannum, collapse = "\n"), "\n\n")
}

# ── Save mapping table ─────────────────────────────────────────────────────
write.csv(map_df, file.path(dm_dir, "clock_cpg_v2_mapping.csv"),
          row.names = FALSE)
cat(sprintf("Saved: clock_cpg_v2_mapping.csv\n\n"))

# ── Re-run Fisher enrichment with recovered probes if any found ────────────
if (nrow(map_df_new) > 0) {
  cat("=== Re-running Fisher with extended probe set ===\n\n")

  # Load DMP results from age_cpg_enrichment.R
  dmp_path <- file.path(dm_dir, "dmp_list.csv")
  if (!file.exists(dmp_path)) {
    cat("dmp_list.csv not found — skipping extended Fisher test.\n")
    cat("Run age_cpg_enrichment.R first, then re-run this script.\n")
  } else {
    tt <- read.csv(dmp_path, row.names = 1)

    # Extended probe sets (direct + coordinate-recovered)
    new_v2_probes <- map_df_new$v2_probe_id

    # Find these new probes in our array (they may also have suffixes)
    array_probes_raw2 <- rownames(tt)
    array_base2 <- sub("_[^_]+$", "", array_probes_raw2)

    # For recovered Hannum
    hannum_extended_base <- c(hannum_71,
                               map_df_new$clock_cpg[map_df_new$clock_cpg %in%
                                                       missing_hannum])
    # Map new v2 IDs to what's in our array
    extra_v2_in_array <- array_probes_raw2[
      sub("_[^_]+$", "", array_probes_raw2) %in%
        sub("_[^_]+$", "", new_v2_probes)]

    in_array_hannum_ext <- c(
      array_probes_raw2[array_base2 %in% hannum_71],
      extra_v2_in_array[sub("_[^_]+$", "", extra_v2_in_array) %in%
                          sub("_[^_]+$", "", map_df_new$v2_probe_id[
                            map_df_new$clock_cpg %in% missing_hannum])]
    )
    in_array_hannum_ext <- unique(in_array_hannum_ext)

    cat(sprintf("Extended Hannum set in array: %d (was %d)\n",
                length(in_array_hannum_ext),
                sum(array_base2 %in% hannum_71)))

    # Fisher at top-N thresholds
    run_fisher_simple <- function(sig, clock, all_p) {
      n_both <- sum(sig %in% clock)
      mat <- matrix(c(n_both, length(sig)-n_both,
                      length(clock)-n_both,
                      length(all_p)-length(sig)-length(clock)+n_both), 2)
      ft <- fisher.test(mat, alternative = "greater")
      c(n_sig=length(sig), n_clock=length(clock), overlap=n_both,
        OR=round(as.numeric(ft$estimate),3), p=signif(ft$p.value,4))
    }

    all_p <- rownames(tt)
    for (topn in c(500, 1000, 5000, 10000)) {
      r <- run_fisher_simple(all_p[1:topn], in_array_hannum_ext, all_p)
      cat(sprintf("  top%d: overlap=%d, OR=%.2f, p=%.4f\n",
                  topn, r["overlap"], r["OR"], r["p"]))
    }
  }
}

cat("\n=== Done ===\n")
