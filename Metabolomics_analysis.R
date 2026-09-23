#############################################################################
## Metabolomics differential-abundance analysis (heart + liver)
##   1. Read the data.
##   2. Documented exclusions: samples removed for reasons independent of the
##      metabolomics data (LTreated2: tissue weight abnormality).
##   3. Outlier screen on the remaining samples of each organ - three checks,
##      run once, same rule for heart and liver:
##        a. within-group correlation: leave-one-out mean Pearson correlation
##           (log2 values) to the other samples of the same group;
##           flagged if < median - 3 x MAD
##        b. silhouette (Euclidean, per-metabolite z-scored data, groups as
##           clusters); flagged if < 0
##        c. robust Mahalanobis distance (MCD) on the first 3 principal
##           components of the z-scored data; flagged if > sqrt(chi2_0.975, 3 df)
##      A sample is removed ONLY if it is flagged by all three checks.
##   4. Differential abundance on the final sample set:
##        log2 -> one-way ANOVA -> Tukey HSD -> BH across metabolites
##      Also run without the screen removals (documented exclusions only).
#############################################################################

library(cluster)
library(robustbase)
options(width = 200)

#############################################################################
                                  ## CONFIG  ##                         
#############################################################################
project_dir <- "path/to/metabolomics_files"
infile      <- "metabolomics_filtered.txt"        #suppl file 3, sheet 2
## ###########################################################################

## ---- 1. Input -------------------------------------------------------------
setwd(project_dir)

raw <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)
metabolites <- as.character(raw[[1]])     # metabolite names (duplicates allowed)
raw[[1]] <- NULL
raw <- raw[, colnames(raw) != "" & !grepl("^Unnamed", colnames(raw)), drop = FALSE]
raw <- as.data.frame(lapply(raw, function(x) suppressWarnings(as.numeric(x))),
                     check.names = FALSE)

## ---- 2. Documented exclusions ----------------------------------------------
## Samples excluded for reasons independent of the metabolomics data.
documented_exclusions <- list(
  H = character(0),
  L = c(LTreated2 = "tissue weight abnormality")
)

## ---- Experimental design helper ------------------------------------------
make_design <- function(prefix, exclude = character(0)) {
  cols <- grep(paste0("^", prefix), colnames(raw), value = TRUE)
  cols <- setdiff(cols, exclude)
  grp  <- ifelse(grepl("Untreated", cols), "Untreated",
          ifelse(grepl("Cre", cols),       "Cre", "Treated"))
  data.frame(sample = cols,
             group  = factor(grp, levels = c("Cre", "Untreated", "Treated")),
             stringsAsFactors = FALSE)
}

## ---- 3. Outlier screen -----------------------------------------------------
outlier_screen <- function(design, mcd_seed = 0) {
  cols <- design$sample
  grp  <- as.character(design$group)
  n    <- length(cols)

  ## log2; non-positive -> NA; per-metabolite median imputation;
  ## drop metabolites that are all missing or have zero variance
  X <- as.matrix(raw[, cols, drop = FALSE]); X[X <= 0] <- NA
  L <- log2(X)
  L <- t(apply(L, 1, function(r) { r[is.na(r)] <- median(r, na.rm = TRUE); r }))
  L <- L[rowSums(is.na(L)) == 0, , drop = FALSE]
  L <- L[apply(L, 1, sd) > 0, , drop = FALSE]
  Xs <- t(L)                                                  # samples x metabolites
  psd <- apply(Xs, 2, function(v) sqrt(mean((v - mean(v))^2)))
  Xz <- sweep(sweep(Xs, 2, colMeans(Xs)), 2, psd, "/")        # z-scored

  ## a. within-group correlation (log2 values, not z-scored)
  C   <- cor(t(Xs))
  loo <- sapply(seq_len(n), function(i) mean(C[i, grp == grp[i] & seq_len(n) != i]))
  corr_cut <- median(loo) - 3 * mad(loo)                      # mad() is already x 1.4826

  ## b. silhouette
  sil <- silhouette(as.integer(factor(grp)), dist(Xz))[, "sil_width"]

  ## c. MCD robust distance on PC1-3
  scores <- prcomp(Xz, center = TRUE, scale. = FALSE)$x[, 1:3]
  set.seed(mcd_seed)
  mcd <- covMcd(scores, use.correction = FALSE)
  rd  <- sqrt(mahalanobis(scores, mcd$center, mcd$cov))
  mcd_cut <- sqrt(qchisq(0.975, df = 3))

  out <- data.frame(sample = cols, group = grp,
                    within_grp_corr = round(loo, 3), corr_flag = loo < corr_cut,
                    silhouette = round(sil, 3),      sil_flag  = sil < 0,
                    MCD_dist = round(rd, 2),         mcd_flag  = rd > mcd_cut,
                    stringsAsFactors = FALSE)
  out$n_flags <- out$corr_flag + out$sil_flag + out$mcd_flag
  out$removed <- out$n_flags == 3
  attr(out, "cutoffs") <- c(corr = corr_cut, silhouette = 0, MCD = mcd_cut)
  attr(out, "n_metabolites") <- nrow(L)
  out[order(-out$n_flags, -out$MCD_dist), ]
}

screen_exclusions <- list()
for (organ in c("H", "L")) {
  scr <- outlier_screen(make_design(organ, names(documented_exclusions[[organ]])))
  cut <- attr(scr, "cutoffs")
  cat(sprintf("\n==== Outlier screen: %s (%d samples, %d metabolites) ====\n",
              ifelse(organ == "H", "Heart", "Liver"), nrow(scr), attr(scr, "n_metabolites")))
  if (length(documented_exclusions[[organ]]))
    cat(sprintf("Documented exclusions: %s\n",
                paste(sprintf("%s (%s)", names(documented_exclusions[[organ]]),
                              documented_exclusions[[organ]]), collapse = "; ")))
  cat(sprintf("Cutoffs: correlation < %.3f | silhouette < 0 | MCD > %.3f\n",
              cut["corr"], cut["MCD"]))
  print(scr, row.names = FALSE)
  screen_exclusions[[organ]] <- scr$sample[scr$removed]
  cat(sprintf("Removed (flagged by all three): %s\n",
              if (length(screen_exclusions[[organ]])) paste(screen_exclusions[[organ]], collapse = ", ") else "none"))
  write.csv(scr, sprintf("%s_outlier_screen.csv", ifelse(organ == "H", "heart", "liver")),
            row.names = FALSE)
}

## final exclusion lists = documented exclusions + screen removals
exclude_final <- list(H = c(names(documented_exclusions$H), screen_exclusions$H),
                      L = c(names(documented_exclusions$L), screen_exclusions$L))
## documented exclusions only (for the with/without comparison)
exclude_doc   <- list(H = names(documented_exclusions$H),
                      L = names(documented_exclusions$L))

## ---- 4a. Differential abundance: ANOVA + Tukey ------------------------------
run_de <- function(raw, metabolites, design, organ, min_n = 4) {
  X <- as.matrix(raw[, design$sample, drop = FALSE])
  X[X <= 0] <- NA
  L <- log2(X)
  grp <- design$group

  rows <- lapply(seq_len(nrow(L)), function(i) {
    d <- data.frame(y = as.numeric(L[i, ]), group = grp)
    d <- d[!is.na(d$y), ]
    mu <- tapply(d$y, d$group, mean)
    n  <- tapply(d$y, d$group, length)
    g  <- function(x) unname(ifelse(is.na(mu[x]), NA, mu[x]))
    r <- data.frame(
      Organ = organ, Metabolite = metabolites[i],
      n_Untreated = unname(n["Untreated"]), n_Cre = unname(n["Cre"]), n_Treated = unname(n["Treated"]),
      mean_Untreated_log2 = g("Untreated"), mean_Cre_log2 = g("Cre"), mean_Treated_log2 = g("Treated"),
      F = NA_real_, p_anova = NA_real_,
      log2FC_Untreated_vs_Cre     = g("Untreated") - g("Cre"),
      p_tukey_Untreated_vs_Cre    = NA_real_,
      log2FC_Treated_vs_Cre       = g("Treated") - g("Cre"),
      p_tukey_Treated_vs_Cre      = NA_real_,
      log2FC_Untreated_vs_Treated = g("Untreated") - g("Treated"),
      p_tukey_Untreated_vs_Treated= NA_real_,
      stringsAsFactors = FALSE, check.names = FALSE)

    ok <- all(!is.na(n)) && min(n) >= min_n && nlevels(droplevels(d$group)) == 3
    if (ok) {
      fit <- aov(y ~ group, data = d)
      a   <- summary(fit)[[1]]
      r$F <- a[["F value"]][1]; r$p_anova <- a[["Pr(>F)"]][1]
      tk  <- TukeyHSD(fit)$group
      getp <- function(x, y) {
        k1 <- paste0(x, "-", y); k2 <- paste0(y, "-", x)
        if (k1 %in% rownames(tk)) tk[k1, "p adj"] else tk[k2, "p adj"]
      }
      r$p_tukey_Untreated_vs_Cre     <- getp("Untreated", "Cre")
      r$p_tukey_Treated_vs_Cre       <- getp("Treated", "Cre")
      r$p_tukey_Untreated_vs_Treated <- getp("Untreated", "Treated")
    }
    r
  })
  res <- do.call(rbind, rows)

  res$q_anova_BH                      <- p.adjust(res$p_anova, method = "BH")
  res$q_tukey_Untreated_vs_Cre_BH     <- p.adjust(res$p_tukey_Untreated_vs_Cre, method = "BH")
  res$q_tukey_Treated_vs_Cre_BH       <- p.adjust(res$p_tukey_Treated_vs_Cre, method = "BH")
  res$q_tukey_Untreated_vs_Treated_BH <- p.adjust(res$p_tukey_Untreated_vs_Treated, method = "BH")
  res[["significant_q<0.05"]] <- ifelse(!is.na(res$q_anova_BH) & res$q_anova_BH < 0.05, "yes", "no")

  res[, c("Organ","Metabolite","n_Untreated","n_Cre","n_Treated",
          "mean_Untreated_log2","mean_Cre_log2","mean_Treated_log2",
          "F","p_anova","q_anova_BH","significant_q<0.05",
          "log2FC_Untreated_vs_Cre","p_tukey_Untreated_vs_Cre","q_tukey_Untreated_vs_Cre_BH",
          "log2FC_Treated_vs_Cre","p_tukey_Treated_vs_Cre","q_tukey_Treated_vs_Cre_BH",
          "log2FC_Untreated_vs_Treated","p_tukey_Untreated_vs_Treated","q_tukey_Untreated_vs_Treated_BH")]
}

## ---- 4b. Run & write --------------------------------------------------------
n_sig <- function(res) sum(res[["significant_q<0.05"]] == "yes")
summary_rows <- list()
for (organ in c("H", "L")) {
  nm <- ifelse(organ == "H", "Heart", "Liver")
  fname <- ifelse(organ == "H", "heart", "liver")
  for (set in c("final", "doc_only")) {
    ex  <- if (set == "final") exclude_final[[organ]] else exclude_doc[[organ]]
    des <- make_design(organ, ex)
    ## the "doc_only" run is only needed when the screen removed something
    if (set == "doc_only" && setequal(exclude_final[[organ]], exclude_doc[[organ]])) next
    tag <- if (set == "final") "" else "_without_screen_removals"
    a <- run_de(raw, metabolites, des, nm)
    write.csv(a, sprintf("%s_metabolomics%s.csv", fname, tag), row.names = FALSE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      organ = nm,
      samples = if (set == "final") "final (after screen)" else "documented exclusions only",
      n_samples = nrow(des),
      excluded = if (length(ex)) paste(ex, collapse = ", ") else "none",
      significant = n_sig(a))
  }
}
cat("\n==== Significant metabolites (q < 0.05) ====\n")
print(do.call(rbind, summary_rows), row.names = FALSE)
