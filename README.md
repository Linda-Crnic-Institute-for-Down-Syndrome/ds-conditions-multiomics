# Multi‑omic Analysis of Co-occurring Conditions in Down Syndrome

Code and data-processing workflows supporting the manuscript:
“Systematic multi‑omic deconvolution of the clinical heterogeneity of Down syndrome”

## Overview
This repository contains the complete analysis framework used to characterize multi‑omic signatures across clinical conditions in individuals with Down syndrome (DS).  
It includes:

* R scripts and functions
* Data preprocessing workflows
* Statistical modeling and pipelines
* Reproducibility environment (via renv)
* Documentation for running analyses end-to-end

Each analysis workflow is presented as a self-contained R Project. The goal is to provide a fully reproducible, transparent workflow consistent with open‑science practices.  

## Repository Structure

```
ds-conditions-multiomics/
│
├── Analysis_1/            # Self-contained R Project directory for specific analysis workflow
│    ├── Analysis_1.R            # Analysis script
│    ├── helper_functions.R      # Associated R functions
│    ├── data/                   # Raw or external data
│    ├── results/                # Resulting tables, processed data, model outputs
│    ├── figures/                # Generated visualizations and plots
│    ├── rdata/                  # Workspace images RDS objects
│    ├── renv.lock               # R package versions for reproducibility
│    └── README.md
├── .zenodo.json           # Metadata for DOI registration (optional)
├── LICENSE.md             # Software license
└── README.md
```

## Data Availability
Input datasets used in this study can be obtained from the associated Synapse repository:
* [Sample metadata](URL)
* [Co-occurring conditions](URL)
* [MSD](URL)
* [SOMAscan](URL)

Download each dataset to the appropriate `/data/` directories within each R project.  

Alternatively, datasets can be obtained via the [INCLUDE Data Hub](URL):  
DOI   
DOI   

## Software & Dependencies
Key packages include:
* renv
* tidyverse  
* ggplot2  
The renv.lock files provide exact versions.

## R Environment Setup and Running Analyses
1. Clone the repository
   ```
   git clone https://github.com/Linda-Crnic-Institute-for-Down-Syndrome/ds-conditions-multiomics.git
   ``` 
2. Change to desired R Project directory and open R project via `.Rproj` file
3. Set up reproducible R environment (requires `renv` package to be installed)
   a. Restore the R environment  
   This will install the exact versions of all R packages but requires matching R version.
   ```
   install.packages("renv")
   renv::restore()
   ```
   b. Initialize the R environment  
   This will install all R packages but will not ensure identical versions.
   ```
   install.packages("renv")
   renv::init(bioconductor = TRUE)
   ```
4. Follow workflow in analysis script.


## Citation
If you use this code, please cite:  
**Manuscript**
Systematic multi‑omic deconvolution of the clinical heterogeneity of Down syndrome.
Authors, Journal, Year. DOI (insert once available)

**Code**
ZENODO DOI BADGE GOES HERE

## License
This project is licensed under the MIT License – see the LICENSE file for details.
