# Tissue mtDNA log2FC + analysis pipeline

One script (`tissue_log2FC_pipeline.R`) that takes per-mouse duplex-sequencing
mutation calls across tissues (dcs.vcf files, which are generated for every sample 
by the duplex seq raw data processing pipeline), computes a 16mo-vs-2mo log2 fold change 
per variant, and after a manual VEP annotation step, analyzes the annotation. 

Each aged condition (`16_mo_rec`, `16_mo_unrec`) is compared against the shared
`2_mo` baseline **separately**, so every tissue produces one result per
condition.

---

## Requirements

- R (≥ 4.0)
- Packages: `openxlsx`, `dplyr`, `tidyr`, `readr`

```r
install.packages(c("openxlsx","dplyr","tidyr","readr"))
```

- An **Ensembl VEP** run, done manually on the
  [VEP web interface](https://jun2026.archive.ensembl.org/Mus_musculus/Tools/VEP?db=core) between Stages 1 and 2.

---

## Setup

Edit the CONFIG block at the top:

```r
data_root  <- "path/to/data"
output_dir <- file.path(data_root, "output")   
```

Expected input layout — one folder per tissue, each with a folder per
timepoint of `.dcs.txt` files (one file per mouse):

```
data_root/
├── Heart/       { 2_mo/  16_mo_rec/  16_mo_unrec/ }
├── Liver/       { 2_mo/  16_mo_rec/  16_mo_unrec/ }
├── Spleen/      { 2_mo/  16_mo_rec/  16_mo_unrec/ }
├── Muscle/      { 2_mo/  16_mo_rec/  16_mo_unrec/ }
└── Intestine/   { 2_mo/  16_mo_rec/  16_mo_unrec/ }
```

The `tissues` list sets each tissue's subfolder and the filename prefix of its
`.dcs.txt` files:

```r
tissues <- list(
  Heart     = list(subdir = "Heart",     file_prefix = "^Heart"),
  Liver     = list(subdir = "Liver",     file_prefix = "^Liver"),
  Spleen    = list(subdir = "Spleen",    file_prefix = "^Duplex_Spleen"),
  Muscle    = list(subdir = "Muscle",    file_prefix = "^Hom"),
  Intestine = list(subdir = "Intestine", file_prefix = "^Hom")
)
```

Timepoints are set by `baseline` and `aged`. Any tissue missing a condition's folder is
skipped for that condition.

---

## How to run

Run stage 1, then the manual VEP step, and then stages 2 and 3.

### Stage 1 — per-tissue log2FC

For each tissue, for the baseline and each aged condition present: read one
`.dcs.txt` per mouse, match on `CHROM/POS/REF/ALT/VARTYPE`, average COV and
MUT.FREQ, and apply the reproducibility filter. Then merge each aged condition
against the 2mo baseline and compute `Log2FC_16mo_vs_2mo`.

**Reproducibility filter:** keep a variant detected in **all** mice at a
timepoint whose among-mouse range ≤ `pmax(range_floor, k_rel * mean)`, 
a frequency-scaled tolerance (`range_floor = 0.005`, `k_rel = 1.0`). Set
`k_rel = 0` to recover a flat absolute `range ≤ range_floor` cutoff.

Outputs (in `output_dir`):
- `<tissue>_<tp>_master.xlsx` : per-timepoint workbook, one sheet per mouse
  (set `write_master_xlsx <- FALSE` to skip)
- `<tissue>_<tp>_all_rep.csv` : matched per-mouse table with AVG.COV and AVG.MUT.FREQ
- `<tissue>_<cond>_all.csv` : merged 2mo-vs-aged table with `Log2FC_16mo_vs_2mo`
  (`<cond>` = `rec` or `unrec`)

### Manual step — VEP

For each `<tissue>_<cond>_all.csv`, format as VCF, run it through VEP on the
Ensembl website, and save the result as `<tissue>_<cond>_VEP_anno.txt` in
`output_dir`.

### Stage 2 — merge log2FC onto the annotation

Left-joins `Log2FC_16mo_vs_2mo` onto each `<tissue>_<cond>_VEP_anno.txt` by
`POS/REF/ALT` and generates `<tissue>_<cond>_anno_with_log2FC.txt`. 
Annotation files that aren't present are skipped.

### Stage 3 — mtDNA selection analysis

Runs on the `<tissue>_<cond>_anno_with_log2FC.txt` files. The tissue × condition
table is built automatically from the config and filtered to files that exist.
Base R; median log2FC with bootstrap 95% CIs (`NBOOT = 2000`, `set.seed(0)`). 
D-Loop positions (POS ≥ 15423 or intergenic) are normalized to a single class.

Outputs (in `output_dir/selection/`):
- `<PREFIX>_impact_medians.csv` — per tissue × condition:
  impact-level (HIGH/MODERATE/LOW/MODIFIER) median log2FC, bootstrap CI,
  `prop_constrained`, `median_vs_syn`, `selection_score`
- `selection_coefficient_by_tissue_all.csv` — one combined table: per tissue ×
  condition, `median_syn`, `median_nonsyn`, `s = median(syn) − median(nonsyn)`,
  and its bootstrap CI

Input columns Stage 3 needs: `CHROM, POS, REF, ALT, Consequence, IMPACT,
SYMBOL, BIOTYPE, Log2FC_16mo_vs_2mo`.

---

## Key parameters

| Parameter           | Default                    | Meaning                                            |
|---------------------|----------------------------|----------------------------------------------------|
| `range_floor`       | 0.005                      | Stage 1: absolute floor of the range tolerance     |
| `k_rel`             | 1.0                        | Stage 1: range tolerance as a fraction of the mean |
| `write_master_xlsx` | TRUE                       | Stage 1: also write the per-timepoint workbooks    |
| `baseline` / `aged` | 2_mo / 16_mo rec, unrec    | timepoints compared                                |
| `NBOOT`             | 2000                       | Stage 3: bootstrap resamples (inside the wrapper)  |

---


