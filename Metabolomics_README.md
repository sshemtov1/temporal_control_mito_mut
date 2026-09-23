# Metabolomics differential-abundance analysis

Differential-abundance analysis of mouse **heart** and **liver**
metabolomics across three groups: `Cre`, `Untreated`, `Treated`. Includes a
reproducible three-part outlier screen and a with/without-screen comparison.

---

## Requirements

- R (≥ 4.0)
- CRAN packages: `cluster`, `robustbase`

```r
install.packages(c("cluster","robustbase"))
```

---

## Setup

Edit the two lines in the CONFIG block at the top:

```r
project_dir <- "path/to/metabolomics_files"   # folder with the input file
infile      <- "metabolomics_filtered.txt"    # tab-delimited input
```

The script `setwd()`s to `project_dir`, so the input is read from there and
all outputs are written there.

---

## Input

`metabolomics_filtered.txt` — tab-delimited, first column = metabolite names
(duplicates allowed), remaining columns = per-sample intensities. Sample
column names encode organ and group by prefix/keyword: they start with `H`
(heart) or `L` (liver), and contain `Untreated`, `Cre`, or otherwise are
treated as `Treated` (e.g. `HUntreated1`, `LCre3`, `HTreated2`).

## What it does

1. **Read** the intensity matrix (non-numeric coerced to `NA`).
2. **Documented exclusions** — samples removed for reasons independent of the
   data (`LTreated2`: tissue weight abnormality). 
3. **Outlier screen** (per organ, same rule for both), three checks:
   - within-group leave-one-out Pearson correlation, flagged if
     `< median − 3 × MAD`
   - silhouette width (groups as clusters), flagged if `< 0`
   - robust Mahalanobis distance (MCD) on PC1–3, flagged if
     `> sqrt(qchisq(0.975, df = 3))`
   A sample is removed **only if flagged by all three**. The MCD step is
   seeded (`mcd_seed = 0`) so results are reproducible.
4. **Differential abundance** on the final sample set:
   log2 → one-way ANOVA → Tukey HSD → Benjamini–Hochberg across metabolites.
   Also re-run with documented exclusions only (no screen removals), for the
   with/without comparison — but only when the screen actually removed
   something.

Minimum group size for a metabolite to be tested: `min_n = 4`.

## Outputs (written to `project_dir`)

- `heart_outlier_screen.csv`, `liver_outlier_screen.csv` — full screen table
  (per-sample flags, distances, cutoffs)
- `heart_metabolomics.csv`, `liver_metabolomics.csv` — DE results on the
  final sample set (means, log2FC, ANOVA F/p, BH q, Tukey p per contrast)
- `<organ>_metabolomics_without_screen_removals.csv` — the doc-only run,
  written only if the screen removed a sample for that organ
- A console summary of significant metabolites (q < 0.05) per run

---

## Notes

- Contrasts reported: Untreated vs Cre, Treated vs Cre, Untreated vs Treated.
- Significance is called on the **BH-adjusted ANOVA q-value** (`q < 0.05`);
  Tukey p-values (also BH-adjusted) give the per-contrast direction.

