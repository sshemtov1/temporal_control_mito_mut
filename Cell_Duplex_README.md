# Cell line pipeline

R pipeline for a duplex-sequencing experiment comparing mitochondrial DNA
mutation frequencies **pre** vs **1 month** in two cell lines, per replicate:

| Line       | Sample # | Description        |
|------------|----------|--------------------|
| `hom`      | 3823     | Hom control        |
| `mfn1`     | 4748     | Hom;Mfn1 knockout  |

For each variant it computes a per-replicate `log2FC` of mutation frequency
relative to that replicate's own pre timepoint (negative = the mutation was
lost over time), annotates the variants with Ensembl VEP, and tests whether
higher-impact classes are depleted relative to synonymous variants.
Used to generate **Fig 6D and 6E**.

---

## Requirements

- R (≥ 4.0)
- R packages: `dplyr`, `purrr`, `stringr`, `tidyr`, `data.table`
  (`data.table` is auto-installed by the script if missing; install the
  others once with
  `install.packages(c("dplyr","purrr","stringr","tidyr"))`)
- An **Ensembl VEP** run — done manually on the
  [VEP web interface](https://jun2026.archive.ensembl.org/Mus_musculus/Tools/VEP?db=core) between stages 1 and 2.

---

## Setup

Everything is driven by one variable at the top of `duplex_cell_analysis.R`:

```r
project_dir <- "path/to/files"   # <-- the only line you must edit
```

All other paths are built from it, so you normally don't touch them:

```
project_dir/               input .dcs.txt files (in Hom/ and Hom_Mfn1/)
├── output/                stage 1 outputs
├── VEP/                   stage 2 inputs + outputs
└── output2/               stage 3 inputs
    └── output3/           stage 3 outputs
```

---

## How to run

The script runs top to bottom in three stages, but with **two manual steps
in between** (variant → VCF formatting, and the VEP web run). Run stage by
stage rather than all at once.

### Stage 1 — Initial duplex processing

**In:** `Cell_<sample>_<rep>.dcs.txt` under `project_dir/`
(the duplex `.dcs.vcf` files, cleaned up in Excel and saved as tab-delimited
`.txt`). Filename format: `Cell_<sample>_<rep>.dcs.txt`, e.g.
`Cell_3823_1_pre.dcs.txt`. Replicates expected: 1, 2, 3.

**Out (to `output/`):**
- `exp_<line>_within_line_long.csv` — per-replicate long table with
  `MUT_FREQ_pre`, and `MUT_FREQ`/`COV`/`log2FC` for the 1mo timepoint
- `exp_<line>_wide.csv` — one `log2FC` column per replicate (Prism-ready)
- `exp_within_line_summary.csv` — per line × replicate summary
  (n, median log2FC, % negative, Wilcoxon p vs 0)
- `exp_<line>_<rep>_shared_pre_1mo.csv` — variants called in both pre and 1mo

Filtering: a variant is kept only if its **pre** frequency ≥ `MIN_MUT_FREQ`
(1e-6) and pre coverage ≥ `MIN_COV` (100).

### Manual step A — build VCFs for VEP

For each line, take the `within_line_long.csv`, split it into one file per
replicate, save each as `<sample>.csv` (naming them
`hom_rep1`, `hom_rep2`, `hom_rep3`, `hom_mfn1_rep1`, …), format as VCF, and
run each through **VEP on the Ensembl website**. Save each result as
`<sample>_all_anno.txt` and place both the `<sample>.csv` and
`<sample>_all_anno.txt` in `project_dir/VEP/`.

### Stage 2 — Merge VEP annotation back onto the variants

Left-joins each `<sample>.csv` with its `<sample>_all_anno.txt` on
`CHROM, POS, REF, ALT`, keeping every variant and appending the VEP columns
(GO dropped, columns in a fixed order).

**In / Out (both in `VEP/`):** writes `<sample>_all_anno_merge.csv`.

Edit the `samples` vector if your sample names differ:

```r
samples <- c("hom_mfn1_rep1", "hom_mfn1_rep2", "hom_mfn1_rep3",
             "hom_rep1",      "hom_rep2",      "hom_rep3")
```

Missing files are skipped with a `SKIP` message rather than erroring.

### Manual step B — stage the merged files

Move (and rename as needed) the `*_all_anno_merge.csv` files you want to
analyze into `project_dir/output2/`.

### Stage 3 — mtDNA selection analysis

Auto-detects every `.csv` in `output2/`, infers each file's group from its
name (contains `mfn1` → `Mfn1-KO`, otherwise `Control`) and its timepoints
from its `log2FC_*` columns, then writes:

**Out (to `output2/output3/`):**
- `<sample>_<tp>_UNIFIED_impact.csv` — per-replicate impact-level table
  (loss rate, unified survival score, selection score vs synonymous)
- `functional_unified_freqadj.csv` — frequency-adjusted functional-class
  effects, summarized across replicates

Each input file must contain: `CHROM, POS, REF, ALT, MUT_FREQ_pre, IMPACT,
Consequence, SYMBOL, BIOTYPE`, plus `log2FC_<tp>` / `COV_<tp>` 

---

## Key parameters

Set near the top of each stage:

| Parameter       | Default | Meaning                                         |
|-----------------|---------|-------------------------------------------------|
| `MIN_MUT_FREQ`  | 1e-6    | Stage 1: noise floor at pre                     |
| `MIN_COV`       | 100     | Stage 1: coverage floor at pre                  |
| `REPLICATES`    | 1,2,3   | Stage 1: replicates to look for                 |
| `baseline_thr`  | 4e-5    | Stage 3: ignore variants below this pre freq    |
| `dloop_start`   | 15423   | Stage 3: positions ≥ this are the D-loop        |
| `group_keyword` | `mfn1`  | Stage 3: filename marker for the KO group       |

---

## Notes
- Stages 1–2 use `dplyr`; stage 3 uses `data.table`.
