# SARS-CoV-2 Variant Early-Growth Prediction Framework

This repository contains code and summary data for predicting SARS-CoV-2 variant outcomes from early genomic growth patterns. The workflow builds variant time series, extracts early-stage features, trains models for peak prevalence and dominance duration, and generates summary figures.

## Repository Structure

```text
.
├── code/                 # Analysis scripts, helper functions, and SHAP notebook
├── data/                 # Input summary datasets used by the pipeline
├── result/figs/          # Tracked output figures
├── .gitignore            # Local/intermediate file ignore rules
├── Project2_Variants_Modeling.Rproj
└── README.md
```

## Code Workflow

The main scripts are in `code/`:

```text
01_load_raw_data.R
02_clean_grouping.R
03_feature_extraction.R
04_duration_modeling.R
04_peak_modeling.R
05_joint_comparison_table.R
SHAP_Analysis.ipynb
```

Run the R scripts in order from the repository root:

```bash
Rscript code/01_load_raw_data.R
Rscript code/02_clean_grouping.R
Rscript code/03_feature_extraction.R
Rscript code/04_duration_modeling.R
Rscript code/04_peak_modeling.R
Rscript code/05_joint_comparison_table.R
```

`SHAP_Analysis.ipynb` generates SHAP feature-importance figures and saves them to `result/figs/`.

## Outputs

Tracked figures are stored in:

```text
result/figs/
```

Local intermediate outputs such as `.rds` files, `code/py_data/`, `result/plots/`, `result/plots2/`, manuscripts, references, and PNAS submission files are ignored by Git.

## Data

The `data/` folder contains summary datasets used by the analysis. Raw genomic data originate from GISAID and are subject to its data-use agreement.

## Contact

Yifan (Franky) Zhang  
University of Connecticut, Department of Statistics
