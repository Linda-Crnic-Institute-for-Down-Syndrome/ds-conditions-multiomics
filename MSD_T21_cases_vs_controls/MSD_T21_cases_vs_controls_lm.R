################################################
# Title: linear regression comparing cases vs. controls (co-occuring conditions) in MSD
# in individuals with Down syndrome
# Author(s):
#   - Srija Chillamcherla
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Linear regression modelling for differential abundance of plasma cytokines and
# immune related factors as measured on by Mesoscale Discovery (MSD) multiplexed
# ELISA.
# See README.md for more details

### Data type(s):
#   A. HTP sample meta data 
#      - HTP_Metadata_v0.5_Synapse.txt
#   B. HTP MSD data
#      - HTP_MSD_Cytokines_Synapse.txt
#

### Workflow:
#   Step 1 - Read in meta data and sample data
#   Step 2 - Linear regression model setup and assessment
#   Step 3 - Model results


### Change Log:
# v1.0
# Initial version
#


# 0 General Setup -----
## initialize renv ----
# renv::init(bioconductor = TRUE) # ~first time only 

## 0.1 Load required libraries ----
library("readxl") # used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse") # data manipulation, visualization, and tidy workflows
library("patchwork") # for assembling plot layouts
library("skimr") # data table summaries
# library("janitor") # data cleaning and other functions
library("rstatix") # for outlier detection 
library("broom.mixed") # for tidying models
library("Hmisc") # statistical summaries and correlation utilities
library("limma")  # linear modeling for high-dimensional omics data
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::count,
  dplyr::rename,
  dplyr::bind_rows
)
# LIBRARY("DEVTOOLS")
library("here") # generates path to current project directory
#

## 0.2 Set required parameters ----
# Input data files
# Make sure you download the required files and name them appropriately 
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_Synapse.txt")
# coc - co-occuring conditions datafile
htp_coc_file <- here("data", "HTP_Cooccuring_conditions_v0.5_Synapse.txt")
# msd - cytokines data file
msd_data_file <- here("data", "HTP_MSD_Cytokines_Synapse.txt")
#
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e", "FALSE" = "#4CAF50", "TRUE" = "#E04B4B")
out_file_prefix <- "MSD_T21_cases_vs_controls_lm.R_" # should match this script title
# End required parameters ###
source(here("helper_functions.R")) # load helper functions
#


# 1.0 Read in data ----
## 1.1 Read in and inspect HTP metadata  ----
meta_data <- htp_meta_data_file %>% 
  read_tsv() %>% 
  filter(Event_name != "Current") %>%  
  rename(Age = Age_at_visit) %>% 
  mutate( # Set factor orders
    Sex = fct_relevel(Sex, "Female"),
    Karyotype = fct_relevel(Karyotype, "Control")
  )
#
# inspect
# meta_data %>% distinct(LabID)
# meta_data %>% count(Karyotype)
#

## 1.2 Read in and inspect HTP co-occurring conditions data ----
coc_data <- htp_coc_file %>% 
  read_tsv() %>% 
  pivot_longer(cols = -c(RecordID:MRAbstractionStatus), names_to = "condition", values_to = "has_cond") %>% 
  mutate(has_cond = as_factor(has_cond)) %>% 
  filter(Event_name != "Current") %>% 
  rename(Age = Age_at_visit) %>% 
  mutate( # Set factor orders
    Sex = fct_relevel(Sex, "Female"),
    Karyotype = fct_relevel(Karyotype, "Control")
  )
#
# inspect
# coc_data %>% distinct(LabID)
# coc_data %>% 
#   inner_join(meta_data) %>% # keeps 244258 rows
#   filter(Karyotype == "T21") %>% # keeps 155828 rows
#   filter(has_cond == TRUE) %>% # keeps 7657 rows
#   count(condition)
#

## 1.3 Read in MSD data ----
msd <- msd_data_file %>% 
  read_tsv() %>% 
  mutate(log2_value = log2(value)) %>% 
  inner_join(meta_data)
#
# inspect
# msd %>% distinct(LabID) # 477 samples
# msd %>% distinct(Analyte) # 54 measurements each for 477 samples = 25758
# 


# 2.0 Linear regression: Data preparation -----
## 2.1 Prep data for linear regressions ----
#
# Filtering is done in two stages/parts
# First stage: the extreme outliers are removed PER FEATURE PER CONDITION 
# Second stage: all the features which do not pass the minimum requirements are removed
coc_prep_dat <- msd %>% 
  inner_join(coc_data %>% filter(!is.na(has_cond))) %>% # keeps 3893130 rows
  filter(Karyotype == "T21") %>% # keeps 2904768 rows
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")), # ensure factor levels in correct order
  ) %>% 
  group_by(Analyte, condition, has_cond) %>% 
  mutate(extrem_analyte_conc = rstatix::is_extreme(log2(value))) %>%
  ungroup() %>% 
  filter(extrem_analyte_conc != TRUE) # keeps 2889084 rows
#
coc_prep_dat_final <- coc_prep_dat %>%
  group_by(condition) %>% 
  # remove conditions which have only TRUE or FALSE
  mutate(has_cond_levels = has_cond %>% fct_drop() %>% levels() %>% length()) %>%
  filter(has_cond_levels > 1) %>%
  # condition × features × has_cond grouping
  group_by(Analyte, condition, has_cond) %>%
  add_count(condition, has_cond) %>%
  rename(min_count = n) %>%
  group_by(Analyte, condition) %>%
  # remove any features which do not have at least 10 counts per TRUE or FALSE in each "feature + condition" combination
  filter(!any(min_count < 10)) %>%
  ungroup() %>%
  # now check balance of categorical covariates
  group_by(Analyte, condition, has_cond) %>%
  mutate(Sex_levels = Sex %>% fct_drop() %>% levels() %>% length()) %>%
  # remove sex-specific conditions (require > 1 level, or lm()/DESeq2 gives error)
  group_by(Analyte, condition) %>%  
  # remove any features from a particular condition if it does not at least 1 sex level in each "feature + condition" combination
  filter(Sex_levels > 1) %>%
  add_count(Sex, name = "sex_count") %>%
  group_by(Analyte, condition) %>%  
  # remove any features from a particular condition if it does not at least 3 samples per sex in each "feature + condition" combination
  filter(!any(sex_count < 3)) %>%
  ungroup() %>%
  # log transform the quantitative measurement
  mutate(log2_analyte_conc = log2(value)) # keeps 1,285,790 rows
#
coc_prep_dat_final %>% 
  distinct(condition)  # 100 conditions remain after all the filters
#
coc_prep_dat_final %>% 
  distinct(LabID) # 346 individuals are being analyzed here
#


# 3.0 For-loop with all the filtered conditions ------
# filter for T21 and run the regression model
# Linear regression (AgeSexSource Fixed model)
#

## 3.1 Run For-loop ----
# Load filtered conditions that pass the filters during the data prep for linear regression (ref to 2.1)
conditions <- coc_prep_dat_final %>% distinct(condition)
#
# initiate the vector
all_cond_analyte_lm_data_final.vec <- character()

for (cond in conditions$condition) {
  cat("Processing condition:", cond, "\n")

  # Step 1: filter data for this condition
  T21_lm_data_final.tmp <- coc_prep_dat_final %>%
    filter(condition == cond) %>%
    select(LabID, RecordID, Age, Sex, Sample_source, Analyte,
           analyte_conc = value, log2_analyte_conc, condition, has_cond) %>%
    nest(data = c(LabID, RecordID, Age, Sex, Sample_source, analyte_conc, log2_analyte_conc, has_cond))

  if (nrow(T21_lm_data_final.tmp) == 0) {
    cat("No data for condition:", cond, "- skipping\n")
    next
  }
  
  # Step 1.5: check categorical variable levels
  dat_check <- T21_lm_data_final.tmp$data[[1]]

  vars_to_check <- c("has_cond", "Sex", "Sample_source")
  invalid_vars <- vars_to_check[sapply(vars_to_check, function(v) nlevels(factor(dat_check[[v]])) < 2)]

  if (length(invalid_vars) > 0) {
    cat("Skipping condition:", cond, "because variables have <2 levels ->", paste(invalid_vars, collapse = ", "), "\n")
    next
  }

  
  # Step 2: fit model if all factors have >=2 levels
  regressions_multi_analyte_AgeSexSource.tmp <- T21_lm_data_final.tmp %>%
    mutate(
      fit = map(data, ~ lm(log2_analyte_conc ~ has_cond + Age + Sex + Sample_source, data = .x)),
      tidied = map(fit, broom::tidy),
      glanced = map(fit, broom::glance),
      augmented = map(fit, broom::augment),
      vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
    )

  # Step 3: extract results
  name.tmp <- cond
  lm_analyte_multi_ageSexSource_results.tmp <- regressions_multi_analyte_AgeSexSource.tmp %>%
    unnest(tidied) %>%
    filter(str_detect(term, "has_cond")) %>%
    select(Analyte, condition, term, log2FoldChange = estimate, pval = p.value) %>%
    mutate(FoldChange = 2^log2FoldChange) %>%
    arrange(pval) %>%
    mutate(BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))) %>%
    select(Analyte, condition, FoldChange, log2FoldChange, pval, BHadj_pval, everything())

  assign(name.tmp, lm_analyte_multi_ageSexSource_results.tmp, pos = 1)
  all_cond_analyte_lm_data_final.vec <- c(all_cond_analyte_lm_data_final.vec, name.tmp)
}
# 

## 3.2 Gather all results ----
all_cond_analytes_lm_data_final_results <- mget(all_cond_analyte_lm_data_final.vec, envir = .GlobalEnv) %>%
  bind_rows(.id = "condition")
#
all_cond_analytes_lm_data_final_results %>%
  distinct(condition) #100 conditions (all filtered conditions should be returned)
#
all_cond_analytes_lm_data_final_results %>% arrange(BHadj_pval) %>%
  filter(BHadj_pval < 0.1) %>%
  count(condition) %>%
  arrange(-n)
# 

## 3.3 Save the file ----
# all_cond_analytes_lm_data_final_results %>%
#   write_tsv(file = here("results", paste0(out_file_prefix, "results_conditions_SexAgeSource", ".txt")))


########## End of Script #############

##### Save workspace    ----
################################################################################
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)

# Reload the workspace
# load(here("rdata", paste0(out_file_prefix, ".RData")))

# session_info ----
date()
sessionInfo()
################################################
