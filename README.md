# Forecasting SARS-CoV-2 Variant Expansion and Persistence

This repository contains the analysis workflow for forecasting two complementary dimensions of SARS-CoV-2 variant behavior from early national frequency trajectories:

- **Subsequent growth:** whether a lineage gains more than 5 percentage points after the initial observation window.
- **Long duration:** whether a lineage remains in circulation for more than 90 days.

The workflow harmonizes Pango lineage labels across countries, identifies sustained national circulation, extracts 28 trajectory features from the first 14, 21, or 28 days, and evaluates binary classifiers using a lineage-disjoint temporal test set. The models use routinely collected variant-frequency data and are intended as a reproducible forecasting framework rather than a prospectively validated operational early-warning system.

## Repository layout

```text
.
├── code/
│   ├── 01_select_countries.R       # Rank countries by lineage diversity and sequencing volume
│   ├── 02_select_variants.R        # Harmonize lineages and retain sustained trajectories
│   ├── 03_extract_features.R       # Build 14/21/28-day features, outcomes, and temporal splits
│   ├── feature_functions.R         # Trajectory-feature utilities
│   ├── modeling_comparison.py      # Shared model training, tuning, and evaluation code
│   ├── modeling_duration.py        # Duration-model entry point
│   ├── modeling_growth.py          # Subsequent-growth model entry point
│   ├── make_modeling_outputs.py    # Manuscript tables, performance plots, and SHAP outputs
│   └── draw_variants.R             # Stacked national variant-frequency figures
├── data/                            # Publicly distributable summary datasets
├── figs/                            # Final manuscript figures and editable Figure 1 source
├── reference/ref.bib               # Manuscript bibliography
├── .gitignore                      # Local data and generated-output exclusions
└── README.md
```

## Data availability

The analysis starts from the daily GISAID lineage summary used in the study. GISAID data are governed by the GISAID Database Access Agreement and are not redistributed in this repository. After obtaining authorized access, place the input file at:

```text
reference/UKHSA-UConn-variant-modelling/variant_modelling/data/summary_GISAID_20240918.csv
```

The file must contain `country`, `date`, `lineage`, `denominator`, and `numerator` columns. The `reference/` directory is ignored except for `reference/ref.bib`.

## Software requirements

The workflow was developed with R, Python 3, Pangolin 4.4, and pangolin-data 1.39.

Core R packages are `dplyr`, `data.table`, `jsonlite`, `tidyr`, and `Rcatch22` (or `catch22`). Figure generation additionally uses `ggplot2`, `patchwork`, and `scales`.

Core Python packages are `interpret-core`, `joblib`, `matplotlib`, `numpy`, `openpyxl`, `pandas`, `pygam`, `rdata`, `scikit-learn`, `scipy`, and `xgboost`. One reproducible setup pattern is:

```bash
python -m venv .venv-modeling
source .venv-modeling/bin/activate
python -m pip install --upgrade pip
python -m pip install interpret-core joblib matplotlib numpy openpyxl pandas pygam rdata scikit-learn scipy xgboost
```

Install the core R packages from an R session if they are not already available:

```r
install.packages(c("dplyr", "data.table", "jsonlite", "tidyr", "Rcatch22", "ggplot2", "patchwork", "scales"))
```

`02_select_variants.R` verifies the Pangolin and pangolin-data versions before processing. If Pangolin is not available on `PATH`, provide its executable explicitly:

```bash
export PANGOLIN_BIN=/absolute/path/to/pangolin
```

## Reproducing the analysis

Run the R preprocessing scripts from `code/` in numerical order:

```bash
cd code
Rscript 01_select_countries.R
Rscript 02_select_variants.R
Rscript 03_extract_features.R
cd ..
```

These steps generate the selected-country list, harmonized eligible trajectories, and lineage-disjoint training and test sets. The temporal cutoff is May 5, 2023; all country-level observations of the same lineage remain in one partition.

Fit and evaluate the duration and subsequent-growth classifiers from the repository root:

```bash
python code/modeling_duration.py
python code/modeling_growth.py
```

Full model runs use lineage-grouped five-fold cross-validation and 1,000 lineage-level bootstrap replicates. For a quick environment check, add `--smoke`. Existing completed tasks are reused unless `--overwrite` is supplied.

Generate the manuscript tables, input-window comparisons, country summaries, and SHAP figures after both model runs finish:

```bash
python code/make_modeling_outputs.py --outcome duration
python code/make_modeling_outputs.py --outcome growth
```

To regenerate the grouped national variant-frequency figures:

```bash
cd code
Rscript draw_variants.R
```

## Generated files

Intermediate `.rds` files, model artifacts, bootstrap results, diagnostics, manuscript files, and temporary review outputs are intentionally excluded from Git. Principal generated locations include:

```text
code/train.rds
code/test.rds
code/variants.rds
code/modeling_results/
result/02_group_variants/
```

Final manuscript figures selected for version control are stored in `figs/`. Model runs also write a manifest containing arguments, random seed, software versions, and package versions to their output directory.

## Reproducibility notes

- Missing dates within an otherwise eligible feature window are represented as zero variant share.
- Lineages with a maximum share above 50% during an input window are excluded from that window.
- Model tuning, threshold selection, and Super Learner weight estimation use training data only; the held-out temporal test set is reserved for final evaluation.
- Raw GISAID data, locally generated model objects, and manuscript drafts must not be committed.

## Contact

Yifan (Franky) Zhang  
University of Connecticut, Department of Statistics  
zfranky6@uconn.edu
