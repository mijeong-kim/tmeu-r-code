# EJS Code Release

This folder contains a GitHub-ready subset of the code used for the paper

> Efficient targeted maximum likelihood estimation of average treatment effects under structured outcome models with unknown error distributions

The release keeps only the scripts and input data needed to reproduce the main simulation results, the imbalanced-assignment benchmarks, the supplementary nonlinear sensitivity check, the NSW empirical analysis, and the figures derived from those outputs.

If you upload this package to GitHub, use the contents of this folder as the repository root.

## Folder layout

- `code/`
  - Simulation, empirical-analysis, and figure-generation scripts.
- `data/`
  - Input data required to run the empirical analysis.
- `results/`
  - Generated `.csv` summaries, replicate files, and LaTeX tables. This folder is kept empty in version control.
- `figures/`
  - Generated `.pdf` figures. This folder is kept empty in version control.

## Core scripts

- `code/run_proposed_repeated_cf_main.R`
  - Balanced-design simulation for the proposed estimator.
- `code/run_simulation_unbalanced_aipw.R`
  - Imbalanced-assignment simulation comparing the proposed estimator with Gaussian OLS, AIPW, BART, and optional forest-based TMLE.
- `code/run_grf_unbalanced_focus.R`
  - Forest-based TMLE benchmark for the heavy-tailed and skewed-mixture imbalanced designs.
- `code/run_nsw_analysis.R`
  - NSW empirical comparison.
- `code/run_saturation_spline_compare.R`
  - Supplementary nonlinear sensitivity analysis using a saturation model approximated by a common spline basis.
- `code/make_main_figures.R`
  - Main simulation and empirical comparison figures.
- `code/make_nsw_supp_figure.R`
  - Supplementary NSW figures.

## R package requirements

The scripts use base R plus the following contributed packages.

- `MASS`
- `grf`
- `bartCause`
- `bartMachine` as a fallback when `bartCause` is unavailable
- `causaldata` for the NSW dataset, with a CSV fallback already included in `data/`
- `splines`

## Typical workflow

Run commands from the root of this folder.

### 1. Main balanced simulation

```bash
Rscript code/run_proposed_repeated_cf_main.R out_prefix=results/main_n500_rcf
```

### 2. Main imbalanced simulation

```bash
Rscript code/run_simulation_unbalanced_aipw.R n=500 reps=1000 split_reps=20 out_prefix=results/unbalanced_aipw_n500_mc1000
Rscript code/run_grf_unbalanced_focus.R n=500 reps=1000 out_prefix=results/grf_unbalanced_focus_n500_mc1000
```

### 3. Supplementary sample-size runs

```bash
Rscript code/run_simulation_unbalanced_aipw.R n=300 reps=1000 split_reps=20 out_prefix=results/unbalanced_aipw_n300_mc1000_bart
Rscript code/run_grf_unbalanced_focus.R n=300 reps=1000 out_prefix=results/grf_unbalanced_focus_n300_mc1000

Rscript code/run_simulation_unbalanced_aipw.R n=1000 reps=1000 split_reps=20 out_prefix=results/unbalanced_aipw_n1000_mc1000_bart
Rscript code/run_grf_unbalanced_focus.R n=1000 reps=1000 out_prefix=results/grf_unbalanced_focus_n1000_mc1000
```

### 4. Supplementary nonlinear sensitivity analysis

```bash
Rscript code/run_saturation_spline_compare.R out_prefix=results/saturation_spline_n500_mc1000
```

### 5. NSW empirical analysis

```bash
Rscript code/run_nsw_analysis.R
```

### 6. Figures

```bash
Rscript code/make_main_figures.R
Rscript code/make_nsw_supp_figure.R
```

## Notes

- Generated outputs are written to `results/` and `figures/` and are intentionally excluded from version control.
- The scripts are direct, cleaned copies of the working analysis scripts used for the manuscript, with project-specific manuscript paths removed.
- If `bartCause` is not installed, `run_nsw_analysis.R` can reuse an existing `results/nsw_table2.csv` row for the BART benchmark.
