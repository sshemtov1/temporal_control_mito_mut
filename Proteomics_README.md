# Proteomics pipeline — DE (limma) + Reactome ORA

Differential-expression and pathway-enrichment analysis of mouse **heart**
and **liver** proteomics (TMT abundances), for three groups (6 replicates
each): `UBC.CRE`, `WT.Treated`, `WT.Untreated`.

For each tissue it filters and normalizes the abundance matrix, runs `limma`
differential expression across three contrasts, and performs **Reactome
over-representation analysis (ORA)** on the significant proteins, using all
quantified proteins in that contrast as the background universe (not the
whole genome).

---

## Requirements

- R (≥ 4.0)
- CRAN packages: `readxl`, `ggplot2`, `ggrepel`, `dplyr`, `stringr`
- Bioconductor packages: `limma`, `clusterProfiler`, `ReactomePA`,
  `org.Mm.eg.db`

```r
install.packages(c("readxl","ggplot2","ggrepel","dplyr","stringr"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","clusterProfiler","ReactomePA","org.Mm.eg.db"))
```

---

## Setup

Edit:

```r
project_dir <- "path/to/proteomics_file"   # folder with the .xlsx file
```

The script `setwd()`s there, so **all inputs are read from and all outputs
are written to that folder.** Point `project_dir` at the folder holding the
two input spreadsheets.

---

## Inputs

Two Excel files (named in the `tissues` list near the top):

| Tissue | File                                                                 |
|--------|----------------------------------------------------------------------|
| heart  | `USC-MS-Core_Data_VermulstLab-Sarah_Proteomics-Heart_02-19-2026.xlsx` |
| liver  | `USC-MS-Core_Data_VermulstLab-Sarah_Proteomics-Liver_02-19-2026.xlsx` |

Each must contain an `Accession` column, a `Gene.Symbol` column, and 18
abundance columns matching the pattern `^Abundance..F` (6 per group, in the
order UBC.CRE → WT.Treated → WT.Untreated). If your column count doesn't
match the 18-entry group vector, the script warns — check the pattern and
column order.

From the supplementary files, it is supplementary file 2, sheets 1 (liver) and 2 (heart)

## Design & thresholds

- **Detection filter:** keep a protein if ≥ 4 of 6 samples are quantified in
  at least one group.
- **Normalization:** log2 → quantile (`normalizeBetweenArrays`).
- **Model:** `~ 0 + group`, `eBayes`.
- **Contrasts:** Treated vs CRE, Untreated vs CRE, Untreated vs Treated.
- **Significance:** `adj.P.Val < 0.05` (`sig_cutoff`); volcano fold-change
  guide line at `logFC = 0.58` (~1.5×).

## Outputs (written to `project_dir`)

Per tissue × contrast:
- `<tissue>_DEP_<contrast>.csv` — full DE table with gene symbols
- `<tissue>_Volcano_<contrast>.pdf` — volcano plot
- `<tissue>_ORA_<contrast>_{ALL,UP,DOWN}.csv` / `.pdf` — Reactome ORA table
  and fold-enrichment dot plot for all / up- / down-regulated DEPs

Per tissue:
- `<tissue>_PCA.pdf` — PCA of normalized abundances

