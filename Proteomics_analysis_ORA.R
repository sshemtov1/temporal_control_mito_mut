# =====================================================
# Proteomics pipeline: DE (limma) + Reactome ORA
# -----------------------------------------------------
# DE analysis parameters:
#   detection filter (>=4 of 6 in >=1 group), log2, quantile normalization, ~0 + group design, same 3 contrasts, eBayes.
# Enrichment: Reactome ORA on significant DEPs, with the background universe set to all quantified/tested proteins
#   (not the whole genome). Run for all / up / down per contrast.
# =====================================================

# ---- Libraries ------------------------
library(readxl)
library(limma)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)

# =====================================================
# CONFIG
# =====================================================

# #############################################################################
project_dir <- "path/to/proteomics_file"
setwd(project_dir)
# #############################################################################

tissues <- list(
  heart = list(
    file              = "USC-MS-Core_Data_VermulstLab-Sarah_Proteomics-Heart_02-19-2026.xlsx",   ## sheet 2 in Suppl file 2
    abundance_pattern = "^Abundance..F"
  ),
  liver = list(
    file              = "USC-MS-Core_Data_VermulstLab-Sarah_Proteomics-Liver_02-19-2026.xlsx",   ## sheet 1 in Suppl file 2
    abundance_pattern = "^Abundance..F"
  )
)

# Shared experimental design.
group_levels <- c("UBC.CRE", "WT.Treated", "WT.Untreated")
group_design <- factor(
  c(rep("UBC.CRE", 6), rep("WT.Treated", 6), rep("WT.Untreated", 6)),
  levels = group_levels
)

sig_cutoff <- 0.05   # adj.P.Val threshold defining significant DEPs

# =====================================================
# FUNCTIONS
# =====================================================

# ---- Expression matrix ----
prepare_expression <- function(df, abundance_pattern, group, min_detect = 4) {

  sample_cols <- grep(abundance_pattern, colnames(df), value = TRUE)

  if (length(sample_cols) != length(group)) {
    warning(sprintf(
      "%d abundance columns found but group vector has %d entries -- check column order/pattern!",
      length(sample_cols), length(group)))
  }

  expr <- as.data.frame(lapply(df[, sample_cols], as.numeric))
  rownames(expr) <- df$Accession

  keep <- apply(expr, 1, function(x) {
    any(sapply(levels(group), function(g) sum(!is.na(x[group == g])) >= min_detect))
  })
  expr <- expr[keep, ]

  expr <- log2(as.matrix(expr))
  expr[is.infinite(expr)] <- NA
  expr <- expr[complete.cases(expr), ]

  normalizeBetweenArrays(expr, method = "quantile")
}

# ---- PCA plot ----
plot_pca <- function(expr_norm, group, tissue) {
  mat <- t(expr_norm)

  nzv <- apply(mat, 2, function(z) var(z) > 0)
  if (any(!nzv)) {
    message(sprintf("[%s] PCA: dropped %d constant protein(s).", tissue, sum(!nzv)))
    mat <- mat[, nzv, drop = FALSE]
  }

  pca <- prcomp(mat, scale. = TRUE)
  imp <- summary(pca)$importance[2, 1:2] * 100
  pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Group = group)

  p <- ggplot(pca_df, aes(PC1, PC2, color = Group)) +
    geom_point(size = 4) +
    theme_minimal() +
    labs(title = paste(tissue, "PCA"),
         x = sprintf("PC1 (%.1f%%)", imp[1]),
         y = sprintf("PC2 (%.1f%%)", imp[2]))

  print(p)
  ggsave(paste0(tissue, "_PCA.pdf"), p, width = 6, height = 5)
  invisible(p)
}

# ---- limma DE ----
run_limma <- function(expr_norm, group) {
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)

  fit <- lmFit(expr_norm, design)

  cm <- makeContrasts(
    Treated_vs_CRE       = WT.Treated   - UBC.CRE,
    Untreated_vs_CRE     = WT.Untreated - UBC.CRE,
    Untreated_vs_Treated = WT.Untreated - WT.Treated,
    levels = design
  )

  fit2 <- eBayes(contrasts.fit(fit, cm))

  results <- lapply(colnames(cm), function(cn) topTable(fit2, coef = cn, number = Inf))
  names(results) <- colnames(cm)
  results
}

# ---- Attach gene symbols ----
attach_symbols <- function(results_table, mapping) {
  results_table$Accession <- rownames(results_table)
  merge(results_table, mapping, by = "Accession", all.x = TRUE)
}

# ---- Volcano ----
make_volcano <- function(results_table, title,
                         fc_cutoff = 0.58, fdr_cutoff = 0.05,
                         top_n = 60, output_file = NULL) {

  df <- results_table
  df$Gene <- df$Gene.Symbol

  df$Significance <- "NS"
  df$Significance[df$adj.P.Val < fdr_cutoff & df$logFC >  fc_cutoff] <- "Up"
  df$Significance[df$adj.P.Val < fdr_cutoff & df$logFC < -fc_cutoff] <- "Down"

  up_idx   <- which(df$Significance == "Up")
  down_idx <- which(df$Significance == "Down")
  top_up   <- head(up_idx[order(df$adj.P.Val[up_idx])], top_n)
  top_down <- head(down_idx[order(df$adj.P.Val[down_idx])], top_n)
  top_genes <- df[c(top_up, top_down), ]

  p <- ggplot(df, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(color = Significance), alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "gray70")) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 50) +
    theme_minimal() +
    theme(legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")) +
    labs(title = title, x = "log2 Fold Change", y = "-log10 FDR")

  print(p)
  if (!is.null(output_file)) ggsave(output_file, p, width = 7, height = 6)
  invisible(p)
}

# ---- helper: parse clusterProfiler "x/y" ratio to numeric ----
parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"),
         function(p) as.numeric(p[1]) / as.numeric(p[2]))
}

# ---- Reactome ORA with detected-protein universe ----
# sig_symbols     : gene symbols of the significant test set (all/up/down)
# universe_symbols: gene symbols of ALL quantified/tested proteins (background)
run_reactome_ORA <- function(sig_symbols, universe_symbols,
                             label, filename_prefix, showCategory = 20) {

  sig_symbols      <- unique(sig_symbols[!is.na(sig_symbols)])
  universe_symbols <- unique(universe_symbols[!is.na(universe_symbols)])

  if (length(sig_symbols) == 0) {
    cat("No significant genes for", label, "\n"); return(NULL)
  }

  # Map both sets to Entrez
  sig_entrez <- suppressWarnings(
    bitr(sig_symbols, "SYMBOL", "ENTREZID", org.Mm.eg.db))$ENTREZID
  uni_entrez <- suppressWarnings(
    bitr(universe_symbols, "SYMBOL", "ENTREZID", org.Mm.eg.db))$ENTREZID

  if (length(sig_entrez) == 0) {
    cat("No genes mapped for", label, "\n"); return(NULL)
  }

  reactome <- enrichPathway(
    gene         = sig_entrez,
    universe     = uni_entrez,          # <-- detected-protein background
    organism     = "mouse",
    pvalueCutoff = 0.05,
    readable     = TRUE
  )

  if (is.null(reactome) || nrow(reactome@result) == 0) {
    cat("No enrichment for", label, "\n"); return(reactome)
  }

  write.csv(reactome@result, paste0(filename_prefix, ".csv"), row.names = FALSE)

  # ---- plot: fold enrichment, top pathways selected by p.adjust ----
  rr <- reactome@result[reactome@result$p.adjust < 0.05, , drop = FALSE]
  if (nrow(rr) == 0) {
    cat("No pathways pass FDR for", label, "\n"); return(reactome)
  }

  if (!"FoldEnrichment" %in% colnames(rr)) {
    rr$FoldEnrichment <- parse_ratio(rr$GeneRatio) / parse_ratio(rr$BgRatio)
  }
  rr$Score       <- -log10(rr$p.adjust)
  rr$Description <- str_wrap(rr$Description, width = 45)

  rr <- rr[order(rr$p.adjust), ]
  rr <- head(rr, showCategory)

  p <- ggplot(rr,
              aes(x = FoldEnrichment,
                  y = reorder(Description, FoldEnrichment),
                  size = Count,
                  color = Score)) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "#6BAED6", high = "#08306B") +
    theme_minimal(base_size = 12) +
    labs(title = label,
         x = "Fold Enrichment", y = "",
         color = "-log10(p.adj)", size = "Gene count") +
    theme(plot.title  = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 9),
          legend.position = "right")

  print(p)
  ggsave(paste0(filename_prefix, ".pdf"), p, width = 7, height = 6)

  reactome
}

# =====================================================
# MAIN: loop over tissues
# =====================================================

all_results <- list()   # DE tables per tissue/contrast
all_ora     <- list()   # ORA objects per tissue/contrast/direction

for (tissue in names(tissues)) {

  cfg <- tissues[[tissue]]
  message("==== Processing: ", tissue, " ====")

  df <- as.data.frame(read_excel(cfg$file))
  colnames(df) <- make.names(colnames(df))

  req  <- c("Accession", "Gene.Symbol")
  miss <- setdiff(req, colnames(df))
  if (length(miss) > 0) {
    stop(sprintf("[%s] Missing columns: %s\nAvailable: %s",
                 tissue, paste(miss, collapse = ", "),
                 paste(colnames(df), collapse = ", ")))
  }

  expr_norm <- prepare_expression(df, cfg$abundance_pattern, group_design, min_detect = 4)
  plot_pca(expr_norm, group_design, tissue)
  results <- run_limma(expr_norm, group_design)

  mapping  <- unique(df[, c("Accession", "Gene.Symbol")])
  ora_list <- list()

  for (cn in names(results)) {
    res <- attach_symbols(results[[cn]], mapping)
    write.csv(res, paste0(tissue, "_DEP_", cn, ".csv"), row.names = FALSE)

    make_volcano(res, paste(tissue, "-", cn),
                 output_file = paste0(tissue, "_Volcano_", cn, ".pdf"))

    # Background universe = ALL quantified/tested proteins for this contrast
    universe_symbols <- res$Gene.Symbol

    # Significant DEPs
    sig <- res[res$adj.P.Val < sig_cutoff, ]

    ora_list[[cn]] <- list(
      all = run_reactome_ORA(
        sig$Gene.Symbol, universe_symbols,
        paste(tissue, "ORA -", cn),
        paste0(tissue, "_ORA_", cn, "_ALL")),

      up = run_reactome_ORA(
        sig$Gene.Symbol[sig$logFC > 0], universe_symbols,
        paste(tissue, "ORA UP -", cn),
        paste0(tissue, "_ORA_", cn, "_UP")),

      down = run_reactome_ORA(
        sig$Gene.Symbol[sig$logFC < 0], universe_symbols,
        paste(tissue, "ORA DOWN -", cn),
        paste0(tissue, "_ORA_", cn, "_DOWN"))
    )

    results[[cn]] <- res
  }

  all_results[[tissue]] <- results
  all_ora[[tissue]]     <- ora_list
}
