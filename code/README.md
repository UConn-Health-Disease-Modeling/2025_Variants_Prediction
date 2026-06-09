# Code Overview

This folder contains the analysis scripts for the COVID-19 variant prediction project.

- `01_load_raw_data.R`: load raw variant, hospitalization, infection, and case data.
- `02_clean_grouping.R`: clean the raw data and define variant groupings.
- `03_feature_extraction.R`: construct analytical periods and extract time-series features.
- `04_duration_modeling.R`: train and evaluate models for variant duration categories.
- `04_peak_modeling.R`: train and evaluate models for peak-share categories.
- `05_joint_comparison_table.R`: combine model outputs into comparison tables.
- `SHAP_Analysis.ipynb`: generate SHAP-based feature importance plots for selected models.
- `classify_lineage.R`: helper code for assigning or grouping lineages.
- `functions_ml_utils.R`: shared machine-learning utility functions.
- `funtions_data_io.R`: shared data input/output helper functions.

Local intermediate files such as `.rds` objects and `py_data/` are ignored by Git and are not required to browse the source code.
