# mtDNA Duplex-Seq & omics analysis pipelines

R analysis code used in manuscript "Temporal control of mitochondrial mutagenesis 
reveals the fate of mtDNA mutations with age" for mitochondrial DNA duplex sequencing
(cell lines and tissues), tissue proteomics and metabolomics analyses. Each pipeline
is a single R script that turns per-sample data into the tables and figures used in the paper. 
Each pipeline has its own README with the details about it.

---

## The four pipelines

| Pipeline | Script | What it does |
|----------|--------|--------------|
| **Duplex — tissue** | `Tissues_Duplex_Analysis.R` | Per-tissue mtDNA mutation frequencies (2 mo vs 16 mo, rec/unrec) → log2FC → VEP annotation → impact medians + selection coefficients |
| **Proteomics** | `Proteomics_analysis_ORA.R` | Differential expression (limma) + Reactome over-representation analysis, heart & liver |
| **Metabolomics** | `Metabolomics_analysis.R` | Outlier screen + differential abundance (ANOVA + Tukey + BH), heart & liver |
| **Duplex — cells** | `Cell_Duplex_Analysis.R` | Per-replicate mtDNA mutation frequencies (pre vs 1 month) → log2FC → VEP annotation → selection analysis |

All four are mouse (*Mus musculus*, mm10 / Ensembl mouse annotation).

---

## What each one is

### Duplex Seq — Tissues (`Tissues_Duplex_Analysis.R`)
Same shape, generalized across tissues (heart, liver, spleen, muscle,
intestine) in one loop. **(1)** Per tissue, matches per-mouse `.dcs.txt` files,
applies a reproducibility filter, and computes 16mo-vs-2mo log2FC for each aged
condition (`rec`, `unrec`) against the shared 2mo baseline. After a
manual **Ensembl VEP** run, **(2)** merges log2FC onto the annotation. **(3)** Writes
per-tissue impact-median tables and a combined selection-coefficient summary
(median syn − median nonsyn, with bootstrap CIs).

### Proteomics (`Proteomics_analysis_ORA.R`)
For heart and liver TMT abundance matrices (3 groups × 6 replicates): detection
filter → log2 → quantile normalization → `limma` differential expression across
three contrasts → volcano plots + PCA, then **Reactome ORA** on the significant
proteins using all quantified proteins as the background universe.

### Metabolomics (`Metabolomics_analysis.R`)
For heart and liver metabolite intensities (3 groups): documented exclusion →
a three-part outlier screen (within-group correlation, silhouette,
robust Mahalanobis distance; a sample is dropped only if all three flag it) →
differential abundance (log2 → one-way ANOVA → Tukey HSD → Benjamini–Hochberg),
run with and without the screen removals.

### Duplex Seq — Cell lines (`Cell_Duplex_Analysis.R`)
Three stages. **(1)** Reads the cleaned `.dcs.txt` files per replicate,
computes each variant's log2FC in mutation frequency vs its own pre timepoint, 
and writes long/wide/summary tables and the shared pre to 1mo sets. After a
manual **Ensembl VEP** run, **(2)** merges the VEP annotation back onto the variants. 
**(3)** Runs an mtDNA selection analysis (impact-level and frequency-adjusted 
functional class tables) used for Fig 6D and 6E. 

---

## Common conventions

- **One config block at the top of each script** — set the data path(s) there;
  nothing else needs editing for a standard layout.
- **The two duplex pipelines share a manual VEP step.** After the first stage
  writes per-sample variant CSVs, you format them as VCFs and run them through
  [Ensembl VEP](https://jun2026.archive.ensembl.org/Mus_musculus/Tools/VEP?db=core) 
  on the web, or however you'd prefer save the results next to the inputs, and run 
  the remainder of the script to finish. However, both are written so the whole
  script can be run safely before *and* after VEP (stages whose inputs aren't
  ready yet are skipped).
- **Outputs are written to an `output/` folder** beside the inputs.

---

## Requirements

- **R ≥ 4.0**
- Per pipeline:

| Pipeline | Packages |
|----------|----------|
| Duplex — tissue | `openxlsx`, `dplyr`, `tidyr`, `readr` |
| Proteomics | CRAN: `readxl`, `ggplot2`, `ggrepel`, `dplyr`, `stringr` · Bioconductor: `limma`, `clusterProfiler`, `ReactomePA`, `org.Mm.eg.db` |
| Metabolomics | `cluster`, `robustbase` |
| Duplex — cells | `dplyr`, `purrr`, `stringr`, `tidyr`, `data.table` |

```r
# CRAN
install.packages(c("dplyr","purrr","stringr","tidyr","readr","openxlsx",
                   "readxl","ggplot2","ggrepel","cluster","robustbase","data.table"))
# Bioconductor (proteomics only)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","clusterProfiler","ReactomePA","org.Mm.eg.db"))
```

The two duplex pipelines also need a manual **Ensembl VEP** run (no install).

---

