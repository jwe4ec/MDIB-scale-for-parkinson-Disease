# MDIB-scale-for-parkinson-Disease

Analysis code for manuscript "Psychometric Evaluation of the Movement Disorders 
Interpretation Bias Scale in Parkinson's Disease"

# Overview

Dandan Tang adapted this code from Jeremy W. Eberle's code (https://github.com/jwe4ec/mdib-hd) 
used in Gibson et al. (2025; https://doi.org/pw62)

# File Organization

Scripts are numbered in the order to be run

Data and results are stored on this repo's corresponding [OSF project](https://osf.io/zt4sd)

```plaintext
.
├── .gitignore   # Used to prevent committing data and results to GitHub
├── README.md
|
├── code/  # Code
|   ├── 1a_define_functions.R
|   ├── 1b_set_flextable_defaults.R
|   ├── 1c_set_officer_properties.R
|   ├── 2_clean_data_compute_item_missingness.R
|   ├── 3_clean_demog_data_create_table.R
|   ├── 4a_define_efa_functions.R
|   ├── 4b_run_efa.R
|   ├── 5_compute_internal_consistency.R
|   ├── 6_compute_test_retest_external_validity_create_tables.R
|   └── 7_run_exploratory_analyses.R
|
├── data/   # Data (stored on OSF)
|   ├── bot_cleaned/
|   |   └── final PD Aim 1 data_deid_2022-12-08_OSF.csv   # Source data
|   ├── further_clean/
|   └── helper/
|
└── results/   # Results (stored on OSF)
    ├── demographics/
    ├── efa_pd/
    └── missingness/
```