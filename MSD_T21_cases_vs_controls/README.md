# Linear modeling of MSD plasma cytokines across co-occurring conditions

-   [Overview]
-   [Repository contents]
-   [System Requirements]
-   [Data sources]
-   [R Environment Setup and Running Analyses]

## Overview

This analysis project is part of the “Systematic multi‑omic deconvolution of the clinical heterogeneity of Down syndrome” manuscript

This workflow performs linear regression modeling of **MSD plasma cytokine measurements** to evaluate associations with co-occurring clinical conditions in participants from the [Human Trisome Project](https://www.trisome.org/) (HTP) in individuals with Down syndrome (DS). The goal is to identify cytokines significantly associated with the co-occurring conditions while adjusting for confounder covariates.

Please refer to the top-level `README.md` in the `ds-conditions-multiomics/` repository for a full overview of all analyses and general data access instructions.

------------------------------------------------------------------------

## Repository contents

```         
Analysis_1/ 
  ├── Analysis_1.R            # Main analysis script 
  ├── helper_functions.R      # Custom R functions used in analysis 
  ├── data/                   # Input datasets (not included in repository) 
  ├── results/                # Model outputs and summary tables 
  ├── plots/                  # Generated plots 
  ├── rdata                   # Workspace images RDS objects 
  ├── renv.lock               # Reproducible package versions 
  └── README.md               # This README file
```

------------------------------------------------------------------------

## System Requirements

The R packages used in this analysis can be run on any standard computer with enough RAM to support the operations.  

This analysis was originally run on a system with 36 GB RAM running MacOS 14.1 and R version 4.3.2.  

The `renv` package can be used to manage the R environment.  

Exact versions of all R packages can be found in the renv.lock file.

------------------------------------------------------------------------

## Data Sources

Human Trisome Project (HTP) datasets used in this study can be obtained from the associated Synapse repository:

\* [Sample metadata and co-occurring conditions](https://doi.org/10.7303/syn3148878)\
\* [MSD plasma immune markers](https://doi.org/10.7303/syn31475487)\

These datasets originate from the [Human Trisome Project](https://www.trisome.org/).  

Download the required files and place them in the `data/` directory before running the analysis.

No pre-processing beyond what is described in the manuscript is required prior to running the script.

------------------------------------------------------------------------

## R Environment Setup and Running Analyses

1.  Clone the repository.

    ```         
    git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/ds-conditions-multiomics.git
    ```

2.  Change to desired R Project directory and open R project via `.Rproj` file.

3.  Set up reproducible R environment (requires `renv` package to be installed).

    Option A. Restore the R environment.\
    This will install the exact versions of all R packages but requires matching R version.

    ```         
    install.packages("renv")
    renv::restore()
    ```

    Option B. Initialize the R environment.\
    This will install all R packages but will not ensure identical versions.

    ```         
    install.packages("renv")
    renv::init(bioconductor = TRUE)
    ```

### Helper functions

The `helper_functions.R` script contains project-specific functions used throughout the analysis, including:

-   Custom ggplot theme setup for consistent figure formatting
-   Functions to visualize results in the form of a volcano plots

These functions are customized for this project and require no modification for standard execution of the workflow.  

