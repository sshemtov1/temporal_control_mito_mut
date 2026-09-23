# ============================================================================
#  mtDNA tissue duplex-seq: per-tissue mutation-frequency processing + log2FC,
#  then (after a manual VEP run) merge log2FC onto the VEP annotation and analyze.
# ----------------------------------------------------------------------------
#  Loops over tissues (heart, liver, spleen, muscle, intestine) and compares  
#  aged conditions (16_mo_rec or 16_mo_unrec) to the 2_mo baseline. 
#
#  STAGE 1  process_tissue():
#     - baseline 2_mo: read one .dcs.txt per mouse, match on CHROM/POS/REF/ALT/
#       VARTYPE, average COV & MUT.FREQ, reproducibility-filter.
#     - for EACH aged condition present (16_mo_rec and/or 16_mo_unrec): same, then
#       merge with the 2_mo baseline and write <tissue>_<cond>_all.csv (includes 
#       Log2FC_16mo_vs_2mo). Conditions whose folder is absent are skipped, so a 
#       tissue with only 16_mo_rec just produces the rec output.
#
#  --- MANUAL: for each <tissue>_<cond>_all.csv, format as VCF, run VEP on the
#      Ensembl website, and save the result as <tissue>_<cond>_VEP_anno.txt in
#      output_dir. ---
#
#  STAGE 2  merge_vep():
#     left-join Log2FC_16mo_vs_2mo onto each <tissue>_<cond>_VEP_anno.txt ->
#     <tissue>_<cond>_anno_with_log2FC.txt. Pairs whose anno file isn't present
#     yet are skipped.
#
#  STAGE 3  run_selection_pipeline():
#     mtDNA mutation analysis pipeline. Fed the <tissue>_<cond>_anno_with_log2FC.txt 
#     files from Stage 2, it writes, into output/selection/:
#       <prefix>_impact_medians.csv               (per tissue x condition)
#       selection_coefficient_by_tissue_all.csv   (one combined summary)
#     The tissue x condition table is built automatically from the same config
#     below, and only files that exist are run.
# ============================================================================

library(openxlsx)
library(dplyr)
library(tidyr)
library(readr)

# ###########################################################################
#                               ##  CONFIG    ##
# ###########################################################################
# data_root holds one subfolder per tissue, each with 2_mo/, 16_mo_rec/ and
# (optionally) 16_mo_unrec/ folders of .dcs.txt files. 

data_root  <- "path/to/data"
output_dir <- file.path(data_root, "output")

# baseline timepoint, and the aged conditions compared against it. Each aged
# condition gets its own output; add/remove entries to change the set.
baseline <- list(dir = "2_mo", tag = "2mo")
aged <- list(dir = "16_mo_rec",   tag = "16mo_rec")

# optionally can include unrecombined group too
#aged <- list(
#  rec   = list(dir = "16_mo_rec",   tag = "16mo_rec"),
#  unrec = list(dir = "16_mo_unrec", tag = "16mo_unrec")
#)

# Reproducibility filter: keep a variant seen in ALL mice at a timepoint, then
# trim technical outliers with a frequency-scaled range tolerance.
range_floor <- 0.005   # absolute floor: protects low-freq variants near noise
k_rel       <- 1.0     # relative tolerance: among-mouse spread up to 100% of mean

write_master_xlsx <- TRUE   # also save the per-timepoint master workbook (1 sheet/mouse)

# One entry per tissue: its subfolder under data_root, and the filename prefix
# of its .dcs.txt files (the only thing that differed between the originals).
tissues <- list(
  Heart     = list(subdir = "Heart",     file_prefix = "^Heart"),
  Liver     = list(subdir = "Liver",     file_prefix = "^Liver"),
  Spleen    = list(subdir = "Spleen",    file_prefix = "^Duplex_Spleen"),
  Muscle    = list(subdir = "Muscle",    file_prefix = "^Hom"),
  Intestine = list(subdir = "Intestine", file_prefix = "^Hom")
)
# ###########################################################################

key_cols <- c("CHROM", "POS", "REF", "ALT", "VARTYPE")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helpers --------------------------------------------------------------

# read every <prefix>*.dcs.txt in a timepoint folder; one data frame per mouse,
# with its MUT.FREQ and COV columns renamed uniquely (MUT.FREQ_<mouse>, COV_<mouse>)
read_mouse_files <- function(tp_dir, file_prefix) {
  files <- list.files(tp_dir, pattern = paste0(file_prefix, ".*\\.txt$"))
  if (length(files) == 0) return(NULL)
  labels <- sub("\\.dcs\\.txt$", "", files)
  setNames(lapply(seq_along(files), function(i) {
    df <- read.table(file.path(tp_dir, files[i]), header = TRUE, sep = "\t",
                     stringsAsFactors = FALSE)
    names(df)[names(df) == "MUT.FREQ"] <- paste0("MUT.FREQ_", labels[i])
    names(df)[names(df) == "COV"]      <- paste0("COV_",      labels[i])
    df
  }), labels)
}

# match mice on key columns 
match_timepoint <- function(mouse_list) {
  trimmed <- lapply(mouse_list, function(df)
    df[, c(key_cols, grep("^(MUT\\.FREQ_|COV_)", names(df), value = TRUE)), drop = FALSE])
  matched  <- Reduce(function(a, b) merge(a, b, by = key_cols, all = TRUE), trimmed)
  cov_cols <- grep("^COV_",        names(matched), value = TRUE)
  mf_cols  <- grep("^MUT\\.FREQ_", names(matched), value = TRUE)
  matched$AVG.COV      <- if (length(cov_cols)) rowMeans(matched[, cov_cols, drop = FALSE], na.rm = TRUE) else NA_real_
  matched$AVG.MUT.FREQ <- rowMeans(matched[, mf_cols, drop = FALSE], na.rm = TRUE)
  matched
}

# keep variants detected in ALL mice, within the frequency-scaled range tolerance
reproducibility_filter <- function(df) {
  rep_cols <- grep("^MUT\\.FREQ_", names(df), value = TRUE)
  m     <- df[, rep_cols, drop = FALSE]
  n_det <- rowSums(!is.na(m))
  rng   <- suppressWarnings(apply(m, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
  mn    <- rowMeans(m, na.rm = TRUE)
  tol   <- pmax(range_floor, k_rel * mn)
  df[n_det == length(rep_cols) & rng <= tol, , drop = FALSE]
}

# suffix every non-key column, so two timepoints can be merged without clashes
rename_with_suffix <- function(df, suffix) {
  other <- setdiff(names(df), key_cols)
  names(df)[names(df) %in% other] <- paste0(other, "_", suffix)
  df
}

# save the per-timepoint master workbook (one sheet per mouse)
save_master <- function(mice, tag, tp_tag) {
  if (!write_master_xlsx) return(invisible(NULL))
  wb <- createWorkbook()
  for (nm in names(mice)) { addWorksheet(wb, nm); writeData(wb, nm, mice[[nm]]) }
  saveWorkbook(wb, file.path(output_dir, sprintf("%s_%s_master.xlsx", tag, tp_tag)), overwrite = TRUE)
}

# read + match + filter one timepoint folder; returns the suffix-renamed filtered
# table, or NULL if the folder has no files
prepare_timepoint <- function(tp_dir, file_prefix, tag, tp_tag) {
  mice <- read_mouse_files(tp_dir, file_prefix)
  if (is.null(mice)) return(NULL)
  message(sprintf("  %-11s: %d mice", basename(tp_dir), length(mice)))
  save_master(mice, tag, tp_tag)
  matched <- match_timepoint(mice)
  write.csv(matched, file.path(output_dir, sprintf("%s_%s_all_rep.csv", tag, tp_tag)), row.names = FALSE)
  rename_with_suffix(reproducibility_filter(matched), tp_tag)
}

# ---- STAGE 1: per-tissue processing --------------------------------------
process_tissue <- function(tissue, cfg) {
  message("==== ", tissue, " ====")
  tdir <- file.path(data_root, cfg$subdir)
  tag  <- tolower(tissue)

  # baseline (2mo): processed once, reused for every aged comparison
  base_filtered <- prepare_timepoint(file.path(tdir, baseline$dir), cfg$file_prefix, tag, baseline$tag)
  if (is.null(base_filtered)) { message("  [MISSING] baseline ", baseline$dir, " - skipping ", tissue); return(invisible(NULL)) }
  a2 <- paste0("AVG.MUT.FREQ_", baseline$tag)

  for (cond in names(aged)) {
    ac <- aged[[cond]]
    aged_filtered <- prepare_timepoint(file.path(tdir, ac$dir), cfg$file_prefix, tag, ac$tag)
    if (is.null(aged_filtered)) { message("  [skip] no ", ac$dir, " folder"); next }

    merged <- merge(base_filtered, aged_filtered, by = key_cols, all = TRUE)
    a16 <- paste0("AVG.MUT.FREQ_", ac$tag)
    merged$Ratio_16mo_vs_2mo  <- merged[[a16]] / merged[[a2]]
    merged$Log2FC_16mo_vs_2mo <- log2(merged$Ratio_16mo_vs_2mo)
    merged <- merged %>% drop_na()          # keep variants present across all replicates

    out <- file.path(output_dir, sprintf("%s_%s_all.csv", tag, cond))
    write.csv(merged, out, row.names = FALSE)
    message(sprintf("    -> %s  (%d variants)", basename(out), nrow(merged)))
  }
  invisible(NULL)
}

message("==== STAGE 1: per-tissue log2FC ====")
for (t in names(tissues)) process_tissue(t, tissues[[t]])

# ---- STAGE 2: merge log2FC onto VEP annotation (after the manual VEP run) --
#   Expects, in output_dir, per tissue x condition:
#     <tissue>_<cond>_all.csv        (Stage 1 output)
#     <tissue>_<cond>_VEP_anno.txt   (VEP web output, tab-delimited)
#   Writes <tissue>_<cond>_anno_with_log2FC.txt. Join is on POS/REF/ALT.
merge_vep <- function(tissue, cond) {
  tag       <- tolower(tissue)
  all_file  <- file.path(output_dir, sprintf("%s_%s_all.csv", tag, cond))
  anno_file <- file.path(output_dir, sprintf("%s_%s_VEP_anno.txt", tag, cond))
  out_file  <- file.path(output_dir, sprintf("%s_%s_anno_with_log2FC.txt", tag, cond))

  if (!file.exists(all_file) || !file.exists(anno_file)) {
    message("  [SKIP] ", tag, " ", cond, " - need ", basename(all_file), " and ", basename(anno_file))
    return(invisible(NULL))
  }
  all_tbl <- read_csv(all_file,  show_col_types = FALSE)
  anno    <- read_tsv(anno_file, show_col_types = FALSE)
  out <- anno %>%
    left_join(all_tbl %>% dplyr::select(POS, REF, ALT, Log2FC_16mo_vs_2mo),
              by = c("POS", "REF", "ALT"))
  write_tsv(out, out_file)
  message("  -> ", basename(out_file))
}

message("\n==== STAGE 2: merge log2FC onto VEP annotation ====")
for (t in names(tissues)) for (cond in names(aged)) merge_vep(t, cond)

# ---- STAGE 3: tissue mtDNA selection pipeline --------------

run_selection_pipeline <- function(input_folder, output_folder, fc_col, samples) {

  NBOOT       <- 2000
  DLOOP_START <- 15423
  set.seed(0)
  NONSYN_RE   <- "missense|stop_gained|frameshift|inframe_insertion|inframe_deletion|stop_lost|start_lost"
  dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

  # Stage 2 always writes a tab-delimited .txt, so just read that.
  read_anno <- function(path)
    read.delim(path, sep = "\t", header = TRUE, quote = "\"",
               stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))

  boot_median_vec <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) return(rep(NA_real_, NBOOT))
    n <- length(x)
    vapply(seq_len(NBOOT), function(i) median(x[sample.int(n, n, replace = TRUE)]), numeric(1))
  }
  ci95 <- function(v) as.numeric(quantile(v, c(0.025, 0.975), na.rm = TRUE, names = FALSE))

  group_summary <- function(df, key, fc, sm, mad) {
    ks <- unique(df[[key]])
    out <- lapply(ks, function(k) {
      x <- df[[fc]][df[[key]] == k]; m <- median(x); ci <- ci95(boot_median_vec(x))
      data.frame(key = k, n_mutations = length(x), median_log2FC = m,
                 median_CI95_low = ci[1], median_CI95_high = ci[2],
                 prop_constrained = mean(x < 0),
                 median_vs_syn = m - sm, selection_score = (m - sm) / mad,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, out)
    names(out)[1] <- key
    out[order(out$median_log2FC), ]
  }

  process_tissue <- function(path, prefix, tissue) {
    a <- read_anno(path)
    a[[fc_col]] <- suppressWarnings(as.numeric(a[[fc_col]]))
    a <- a[is.finite(a[[fc_col]]), , drop = FALSE]
    if (!"DOMAINS"  %in% names(a)) a$DOMAINS  <- NA
    if (!"BLOSUM62" %in% names(a)) a$BLOSUM62 <- NA

    # D-Loop normalization
    posn <- suppressWarnings(as.numeric(a$POS))
    is_dl <- a$SYMBOL %in% c("D-Loop","D-LOOP") |
             a$Consequence %in% c("intergenic","intergenic_variant","D-Loop","D-LOOP","D_Loop") |
             toupper(as.character(a$IMPACT)) %in% c("D-LOOP","INTERGENIC") |
             (posn >= DLOOP_START)
    is_dl[is.na(is_dl)] <- FALSE
    a$Consequence[is_dl] <- "intergenic"; a$SYMBOL[is_dl] <- "D-Loop"; a$IMPACT[is_dl] <- "MODIFIER"

    # synonymous neutral reference
    syn <- a[[fc_col]][a$Consequence == "synonymous_variant"]
    sm  <- median(syn); mad <- median(abs(syn - sm))

    of <- function(name) file.path(output_folder, name)
    fc <- fc_col

    # impact-level median log2FC, bootstrap CI, selection score
    write.csv(group_summary(a, "IMPACT", fc, sm, mad),
              of(sprintf("%s_impact_medians.csv", prefix)), row.names = FALSE)

    # selection coefficient: median(syn) - median(nonsyn), with bootstrap CI
    nsyn <- a[[fc]][grepl(NONSYN_RE, as.character(a$Consequence))]
    s_coef <- median(syn) - median(nsyn)
    bs <- vapply(seq_len(NBOOT), function(i)
            median(syn[sample.int(length(syn), length(syn), TRUE)]) -
            median(nsyn[sample.int(length(nsyn), length(nsyn), TRUE)]), numeric(1))
    r3 <- data.frame(tissue = tissue, n_syn = length(syn), n_nonsyn = length(nsyn),
                     median_syn_log2FC = median(syn), median_nonsyn_log2FC = median(nsyn),
                     s = s_coef, ci_lower = quantile(bs, 0.025, names = FALSE),
                     ci_upper = quantile(bs, 0.975, names = FALSE), stringsAsFactors = FALSE)

    cat(sprintf("  %-12s | syn_median=%.3f syn_MAD=%.3f | s=%.3f\n", tissue, sm, mad, s_coef))
    r3
  }

  cat("Generating impact medians + selection-coefficient summary per tissue...\n")
  A3 <- list()
  for (i in seq_len(nrow(samples))) {
    path <- file.path(input_folder, samples$file[i])
    if (!file.exists(path)) { cat("  MISSING:", samples$file[i], "-- skipped\n"); next }
    A3[[i]] <- process_tissue(path, samples$prefix[i], samples$tissue[i])
  }
  write.csv(do.call(rbind, A3), file.path(output_folder, "selection_coefficient_by_tissue_all.csv"), row.names = FALSE)
  cat("Done. Selection outputs in: ", output_folder, "\n", sep = "")
}

message("\n==== STAGE 3: tissue mtDNA selection pipeline ====")
# Build the tissue x condition table automatically from the config, pointing at
# the Stage 2 outputs, and keep only the files that actually exist.
sel_samples <- do.call(rbind, lapply(names(tissues), function(t) {
  do.call(rbind, lapply(names(aged), function(cond) {
    data.frame(file   = sprintf("%s_%s_anno_with_log2FC.txt", tolower(t), cond),
               prefix = sprintf("%s_%s", toupper(t), cond),
               tissue = sprintf("%s_%s", t, cond),
               stringsAsFactors = FALSE)
  }))
}))
sel_samples <- sel_samples[file.exists(file.path(output_dir, sel_samples$file)), , drop = FALSE]

if (nrow(sel_samples) == 0) {
  message("  [SKIP] no *_anno_with_log2FC.txt files yet - run Stage 2 (post-VEP) first")
} else {
  run_selection_pipeline(
    input_folder  = output_dir,
    output_folder = file.path(output_dir, "selection"),
    fc_col        = "Log2FC_16mo_vs_2mo",
    samples       = sel_samples
  )
}
