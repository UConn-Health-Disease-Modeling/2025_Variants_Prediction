# SARS-CoV-2 Variant Early-Growth Prediction Framework

This repository implements a sequence-only forecasting framework that predicts the peak prevalence and dominance duration of emerging SARS-CoV-2 variants using early genomic growth signals. The pipeline is modular, reproducible, and based on multi-country GISAID data (2020–2024).  
See the accompanying manuscript for methodological details.

## Repository Structure

The top-level directory contains all modeling code, utilities, and processed data objects used throughout the workflow.

```
Project2_Variants_Modeling/
│
├── 01_load_raw_data.R              # Load lineage counts, raw share trajectories
├── 02_clean_grouping.R             # Cleaning, analytic period selection, variant grouping
├── 03_feature_extraction.R         # Handcrafted + catch22 feature generation
├── 04_duration_modeling.R          # Duration (>10%) prediction models
├── 04_peak_modeling.R              # Peak share prediction models
├── 05_joint_comparison_table.R     # Combine predictions and compare with observed outcomes
│
├── classify_lineage.R              # Helper functions for lineage → variant grouping
├── functions_ml_utils.R            # ML helper functions (SuperLearner setup, metrics, SHAP utils)
├── funtions_data_io.R              # Data I/O utilities for loading/saving intermediate objects
├── model_summary.R                 # Summaries, tables, final evaluation outputs
│
├── all_data.rds                    # Full processed dataset (all lineages, all countries)
├── data3_and_dropped.rds           # Cleaned dataset after exclusions (analytic periods)
├── feat_list.rds                   # Saved feature list (handcrafted + catch22)
│
├── peak_results_1117.rds           # Saved peak-share prediction results
├── duration_results_1117.rds       # Saved duration prediction results
```

## Pipeline Overview

1. **Data Loading**  
   Run `01_load_raw_data.R` to import raw lineage-level share data and construct initial time-series matrices.

2. **Cleaning and Variant Grouping**  
   `02_clean_grouping.R` applies:
   - analytic-period filtering  
   - minimum sequencing requirements  
   - lineage → variant grouping rules  
   Output is stored in `data3_and_dropped.rds`.

3. **Feature Engineering**  
   `03_feature_extraction.R` computes:
   - handcrafted epidemiological features (slopes, runs, early maxima)  
   - 22 catch22 time-series descriptors  
   Output is saved in `feat_list.rds`.

4. **Peak and Duration Modeling**  
   - `04_peak_modeling.R`
   - `04_duration_modeling.R`  
   Each script trains GLM, GAM, Elastic Net, SVM, CART, neural networks, and a SuperLearner ensemble.  
   Predictions are stored in:
   - `peak_results_1117.rds`
   - `duration_results_1117.rds`

5. **Joint Comparison and Tables**  
   `05_joint_comparison_table.R` merges predictions with observed metrics, producing the final evaluation tables.

## Running the Pipeline

Execute the full workflow in order:

```
Rscript 01_load_raw_data.R
Rscript 02_clean_grouping.R
Rscript 03_feature_extraction.R
Rscript 04_peak_modeling.R
Rscript 04_duration_modeling.R
Rscript 05_joint_comparison_table.R
```

Intermediate RDS files are generated automatically in the working directory.

## Data Availability

Raw genomic data originate from GISAID and follow its data-use agreement.  
This repository includes only processed, non-sensitive summary datasets.

## Contact

Yifan (Franky) Zhang  
University of Connecticut, Department of Statistics
