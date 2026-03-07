# Replication data and scripts for "Interdependent cooperation distinguishes chimpanzees and humans from bonobos" 

by Samuni, Neumann et al.

for review only

## Requirements

  - `R` (v. >= 4.5.0)
  - `Stan` (v. >= 2.36.0)
  - `cmdstanr` (v. >= 0.9.0.9000)
  - `bamoso` (see `03_panassociations.R`)


## Workflow

We fitted four models and for each we stored results on a local machine but did not include these result files here for storage space reasons.
In order to reproduce the figures in the manuscript, you need to rerun the models on your machine.
We provide additional code to run smaller versions (with data subsets) of the models (that fit substantially faster compared to the models that are presented in the manuscript) if desired.


## Scripts

The following lists the data sets, Stan model files and R scripts to run the analyses.

### *Pan* grooming model

  - `data_files/data_grooming.csv` and `data_files/data_grooming_predictions.csv`
  
  - `pan_grooming.stan`
  
  - `01_grooming.R` and `01_grooming_figures.R`

### *Pan* cooperation model

  - `data_files/pandata.csv`
  
  - `pan_cooperation.stan`
  
  - `02_panmodel.R` and `02_panmodel_figures.R`

### *Pan* association model

  - `data_files/pandata.csv`

  - `03_panassociations.R` and `03_panassociations_figures.R`


### Bayaka cooperation

  - `data_files/human_data.csv`
  
  - `bayaka_model.stan`
  
  - `04_bayakamodel.R` and `04_bayakamodel_figures.R`



## Reproducibility

The sampling algorithms of Stan cannot reliably be reproduced across different machines.
We still added seeds so that we (and anyone else) can at least reproduce results on a specific machine.
