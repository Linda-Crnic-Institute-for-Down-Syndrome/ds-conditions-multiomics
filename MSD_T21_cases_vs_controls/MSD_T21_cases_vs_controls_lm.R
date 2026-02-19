################################################
# Title: Differential abundance analysis of plasma MSD cytokines and immune markers by co-occurring condition status (cases vs. controls) in individuals with Down syndrome (T21)
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
# ELISA from individuals with Down syndrome (T21)
# See README.md for more details

### Data type(s):
# The Human Trisome Project (HTP) is a large-scale longitudinal, multi-omics research initiative 
# generating multi-omics data to advance understanding of Down syndrome (Trisomy 21). 
#   A. HTP sample meta data 
#      - HTP_Metadata_v0.5_Synapse.txt
#   B. HTP co-occuring conditions data
#       - HTP_Cooccuring_conditions_v0.5_Synapse.txt
#   C. HTP MSD data
#      - HTP_MSD_Cytokines_Synapse.txt
#

### Workflow:
#   Step 1 - Read in meta data and sample data
#   Step 2 - Linear regression model setup and assessment
#   Step 3 - Model results
#   Step 4 - Visualize results


### Change Log:
# v1.0
# Initial version
#

# 0 General Setup -----
## Initialize renv ----
# renv::init(bioconductor = TRUE) # for first time only 
#
# To install the exact versions of all R packages but requires matching R version
# renv::restore()

## 0.1 Load required libraries ----
library("readxl") # used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse") # data manipulation, visualization, and tidy workflows
library("skimr") # data table summaries
library("rstatix") # for outlier detection 
library("limma")  # linear modeling for high-dimensional omics data
library("broom") # extract tidy model results
library("ggplot2") # to create plots
library("ggrepel") # required for labelling features
library("ggforce") # required for zooming and sina
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::count,
  dplyr::rename,
  dplyr::bind_rows
)
library("here") # generates path to current project directory
#

## 0.2 Set parameters ----
# HTP datasets used in this study can be obtained from the 
# associated Synapse repository (further details in README.md)
# Download each dataset to the appropriate /data/ directories within R project.
# Input data files
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_Synapse.txt")
#
htp_cooccurring_data_file <- here("data", "HTP_Cooccuring_conditions_v0.5_Synapse.txt")
# msd - cytokines data file
msd_data_file <- here("data", "HTP_MSD_Cytokines_Synapse.txt")
#
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e", "FALSE" = "#4CAF50", "TRUE" = "#E04B4B")
out_file_prefix <- "MSD_T21_cases_vs_controls_lm.R_" # should match this script title
# End required parameters ###
source(here("helper_functions.R")) # load helper functions
#


# 1 Read in data ----
## 1.1 Read in and inspect HTP metadata  ----
meta_data <- htp_meta_data_file %>% 
  read_tsv() %>% 
  filter(Event_name != "Current") %>%  
  mutate( # Set factor orders
    Sample_source_code = as_factor(Sample_source_code), # convert to factor - default is numerical order
    Sex = fct_relevel(Sex, "Female"),
    Karyotype = fct_relevel(Karyotype, "Control")
  ) %>% 
  select(-c("Data_contact", "Date_exported", "Script"))
#

## 1.2 Read in and inspect HTP co-occurring conditions data ----
htp_cooccurring_data <- htp_cooccurring_data_file %>% 
  read_tsv() %>% 
  rename(condition = Condition,
         has_cond = History_of_condition) %>% 
  mutate(has_cond = as_factor(has_cond)) %>% 
  select(-c("Data_contact", "Date_exported", "Script")) %>% 
  inner_join(meta_data) %>% 
  distinct() %>% 
  mutate( # Set factor orders
    Sample_source_code = as_factor(Sample_source_code), # convert to factor - default is numerical order
    Sex = fct_relevel(Sex, "Female"),
    Karyotype = fct_relevel(Karyotype, "Control")
  )
#

## 1.3 Read in MSD data ----
msd <- msd_data_file %>% 
  read_tsv() %>% 
  mutate(log2_value = log2(value)) %>% 
  inner_join(meta_data)
#
# inspect
msd %>% distinct(LabID) # 477 samples
msd %>% distinct(Analyte) # 54 measurements each for 477 samples = 25758

# 2 Linear regression: Data preparation -----
## 2.1 Prep data for linear regressions ----
#
# Filtering is done in two stages/parts
# First stage: the extreme outliers are removed PER FEATURE PER CONDITION 
# Second stage: all the features which do not pass the minimum requirements are removed
regression_prep_dat <- msd %>% 
  inner_join(htp_cooccurring_data %>% filter(!is.na(has_cond))) %>% # keeps 3893130 rows
  filter(Karyotype == "T21") %>% # keeps 2904768 rows
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")), # ensure factor levels in correct order
    has_cond = fct_relevel(has_cond, c("FALSE", "TRUE")) # ensure factor levels in correct order
  ) %>% 
  group_by(Analyte, condition, has_cond) %>% 
  mutate(extrem_analyte_conc = rstatix::is_extreme(log2(value))) %>%
  ungroup() %>% 
  filter(extrem_analyte_conc != TRUE) # keeps 2889084 rows
#
regression_prep_dat_final <- regression_prep_dat %>%
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
regression_prep_dat_final %>% 
  distinct(condition)  # 100 conditions remain after all the filters
#
regression_prep_dat_final %>% 
  distinct(LabID) # 346 individuals are being analyzed here
#


# 3 For-loop with all the filtered conditions ------
# filter for T21 and run the regression model
# Linear regression (AgeSexSource Fixed model)
#
## 3.1 Run For-loop ----
# Load filtered conditions that pass the filters during the data prep for linear regression (ref to 2.1)
conditions <- regression_prep_dat_final %>% distinct(condition)
#
# initiate the vector
all_cond_analyte_lm_data_final.vec <- character()

for (cond in conditions$condition) {
  cat("Processing condition:", cond, "\n")
  
  # Step 1: filter data for this condition
  T21_lm_data_final.tmp <- regression_prep_dat_final %>%
    filter(condition == cond) %>%
    select(LabID, Age, Sex, Sample_source_code, Analyte,
           analyte_conc = value, log2_analyte_conc, condition, has_cond) %>%
    nest(data = c(LabID, Age, Sex, Sample_source_code, analyte_conc, log2_analyte_conc, has_cond))
  
  if (nrow(T21_lm_data_final.tmp) == 0) {
    cat("No data for condition:", cond, "- skipping\n")
    next
  }
  
  # Step 1.5: check categorical variable levels
  dat_check <- T21_lm_data_final.tmp$data[[1]]
  
  vars_to_check <- c("has_cond", "Sex", "Sample_source_code")
  invalid_vars <- vars_to_check[sapply(vars_to_check, function(v) nlevels(factor(dat_check[[v]])) < 2)]
  
  if (length(invalid_vars) > 0) {
    cat("Skipping condition:", cond, "because variables have <2 levels ->", paste(invalid_vars, collapse = ", "), "\n")
    next
  }
  
  
  # Step 2: fit model if all factors have >=2 levels
  regressions_multi_analyte_AgeSexSource.tmp <- T21_lm_data_final.tmp %>%
    mutate(
      fit = map(data, ~ lm(log2_analyte_conc ~ has_cond + Age + Sex + Sample_source_code, data = .x)),
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
    # group_by(condition) %>% 
    mutate(BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))) %>%
    # ungroup() %>% 
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

## 3.3 Export the results ----
all_cond_analytes_lm_data_final_results %>%
  write_tsv(file = here("results", paste0(out_file_prefix, "results_conditions_SexAgeSource", ".txt")))
#


# 4 Plots ----
## 4.1 Volcano plot ----
# Using Obesity as an example here
v_plot <- all_cond_analytes_lm_data_final_results %>%
  filter(condition == "Obesity") %>% # change to any other condition(s)
  volcano_plot_lab_lm(
    title="Diff. abundance of MSD cytokines in Obesity (cases vs. controls) - linear regression",
    subtitle = paste0( "T21s only; AgeSexSource adjusted model\n[Down: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange <1) %>% nrow(),
                       "; Up: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange >1) %>% nrow(), "]")
  ) 
#
v_plot
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "lm_volcano_plot_sig_cond", ".pdf")), device = cairo_pdf, width = 7, height = 6, units = "in")
#

## 4.2 Sina plots of top associations ----
# Adjust the dataset for covariates to align with the linear regression model
#
### 4.2.1 Adjust the MSD dataset ----
msd_sample_data <- msd %>% 
  select(-value) %>% 
  distinct()
msd_unadj_data <- msd %>%
  select(LabID, Analyte, log2_value) %>% 
  spread(key = LabID, value = log2_value) %>% # need to be log2 transformed for batch correction
  column_to_rownames(var = "Analyte")
msd_sex_vec <- msd_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>% 
  enframe(name=NULL, value="LabID") %>% 
  inner_join(
    meta_data %>% select(LabID, Sex)
  ) %>% 
  # distinct() %>% 
  pull(Sex)
msd_age_vec <- msd_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>% 
  enframe(name=NULL, value="LabID") %>% 
  inner_join(
    meta_data %>% select(LabID, Age)
  ) %>% 
  # distinct() %>% 
  pull(Age)
msd_source_vec <- msd_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Sample_source_code)
  ) %>%
  #distinct() %>%
  pull(Sample_source_code)
msd_design <- msd_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Karyotype)
  ) %>%
  distinct(LabID, Karyotype) %>%
  model.matrix(~Karyotype, data = .)
#
msd_SourceSexAge_adj <- msd_unadj_data %>%
  limma::removeBatchEffect(
    batch = msd_source_vec, # categorical as batch
    batch2 = msd_sex_vec, # categorical as batch
    covariates = msd_age_vec, # numeric as co-variate
    design = msd_design
  ) %>%
  as_tibble(rownames = "Analyte") %>%  # convert back to tibble
  gather(-Analyte, key = LabID, value = value) %>%
  mutate(value_adj = 2^value) %>% 
  rename(log2_value_adj = value) %>% # rename
  inner_join(msd_sample_data) 
#

### 4.2.2 sina plots ----
# Using Obesity as an example here
msd_SourceSexAge_adj %>% 
  filter(Karyotype == "T21") %>% # have to be only ind. with DS
  filter(Analyte %in% c(all_cond_analytes_lm_data_final_results %>%
                          filter(condition == "Obesity", BHadj_pval < 0.1) %>%
                          slice_max(order_by = log2FoldChange, n = 5) %>%
                          pull(Analyte))) %>% 
  # to get the correct order of the features based on effect size and significance
  mutate(Analyte = fct_relevel(Analyte, c(all_cond_analytes_lm_data_final_results %>%
                                            filter(condition == "Obesity", BHadj_pval < 0.1) %>%
                                            slice_max(order_by = log2FoldChange, n = 5) %>%
                                            pull(Analyte)))) %>% 
  # filtered to the condition of interest
  inner_join(htp_cooccurring_data %>%  filter(condition == "Obesity")) %>% 
  # remove any participants without any status for the condition
  filter(!is.na(has_cond)) %>% 
  # extreme outliers are calculated per feature and per status of the condition
  group_by(Analyte, has_cond) %>%
  mutate(extreme_score = rstatix::is_extreme(log2_value_adj)) %>%
  ungroup() %>%
  filter(extreme_score == FALSE) %>% 
  # relevling the status of condition 
  mutate(has_cond = fct_relevel(has_cond, c("FALSE", "TRUE"))) %>% 
  ggplot(aes(x = has_cond,  y = log2_value_adj, color = has_cond)) +
  scale_color_manual(values = c("FALSE" = "#4CAF50", "TRUE" = "#E04B4B")) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  labs(title = "Obesity",
       subtitle = paste("Significant associations\nAgeSexSource adj\nextreme outliers are removed"),
       y = "log2_adj_val") +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    aspect.ratio = 1.3,
  ) +
  facet_wrap(~Analyte, ncol = 5, scales = "free_y")
#
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "obesity_sig_sina_top5", ".pdf")), device = cairo_pdf, width = 10, height = 4, units = "in")


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
