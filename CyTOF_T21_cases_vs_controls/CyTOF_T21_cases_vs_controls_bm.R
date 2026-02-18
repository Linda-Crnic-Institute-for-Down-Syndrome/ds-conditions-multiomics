################################################
# Title: Differential abundance analysis of immune cell populations by co-occurring condition status (cases vs. controls) in individuals with Down syndrome (T21) 
# Author(s):
#   - Srija Chillamcherla
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Beta regression modelling for differential abundance of immune cell populations
# measured by the CyTOF FlowSOMv2 from individuals with Down syndrome (T21).
# See README.md for more details
# 

### Data type(s):
# The Human Trisome Project (HTP) is a large-scale longitudinal, multi-omics research initiative 
# generating multi-omics data to advance understanding of Down syndrome (Trisomy 21) 
#   A. HTP sample meta data 
#      - HTP_Metadata_v0.5_Synapse.txt
#   B. HTP co-occuring conditions data
#       - HTP_Cooccuring_conditions_v0.5_Synapse.txt
#   B. Human Trisome Project Mass cytometry v2 (CyTOF)
#      - HTP_CyTOF_CD45posCD66low_FlowSOMv2_cluster_percentage_Synapse.txt
#

### Workflow:
#   Step 1 - Read in meta data and sample data
#   Step 2 - Beta regression model (bm) setup and assessment
#   Step 3 - Model results
#   Step 4 - Visualize bm results


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
## 0.1 Load required libraries ----
library("readxl") # used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse") # data manipulation, visualization, and tidy workflows
library("skimr") # data table summaries
library("rstatix") # for outlier detection 
library("betareg") # to run beta regression
library("limma")  # linear modeling for high-dimensional omics data
library("ggplot2") # to create plots
library("ggrepel") # required for labelling features
library("ggforce") # required for zooming and sina
library("broom") # extract tidy model results
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::count,
  dplyr::rename,
  dplyr::bind_rows
)
library("here") # generates path to current project directory
#

## 0.2 Set file name parameters ----
# HTP datasets used in this study can be obtained from the 
# associated Synapse repository (further details in README.md)
# Download each dataset to the appropriate /data/ directories within R project.
# Input data files
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_Synapse.txt")
#
htp_cooccurring_data_file <- here("data", "HTP_Cooccuring_conditions_v0.5_Synapse.txt")
# CyTOF data file
cytof_data_file <- here("data", "HTP_CyTOF_CD45posCD66low_FlowSOMv2_cluster_percentage_Synapse.txt")
#
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e", "FALSE" = "#4CAF50", "TRUE" = "#E04B4B")
out_file_prefix <- "CyTOF_T21_cases_vs_controls_bm.R_" # should match this script title
# End required parameters ###
source(here("helper_functions.R")) # load helper functions
#


# 1 Read in data ----
## 1.1 Read in and inspect HTP metadata  ----
meta_data <- htp_meta_data_file %>% 
  read_tsv() %>% 
  filter(Event_name != "Current") %>%  
  mutate( # Set factor orders
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
    Sex = fct_relevel(Sex, "Female"),
    Karyotype = fct_relevel(Karyotype, "Control")
  )
#

## 1.3 Read in CyTOF data ----
cytof <- cytof_data_file %>% 
  read_tsv() %>% 
  rename(proportion_subset = Value,
         cluster = Cell_cluster_name)
#
# inspect
cytof %>% distinct(LabID) # 388 samples
cytof %>% distinct(uniqueID) 
# 


# 2 Linear regression: Data preparation -----
## 2.1 Prep data for linear regressions ----
#
# Filtering is done in two stages/parts
# First stage: the extreme outliers are removed PER FEATURE PER CONDITION 
# Second stage: all the features which do not pass the minimum requirements are removed
regression_prep_dat <- cytof %>% 
  inner_join(htp_cooccurring_data %>% filter(!is.na(has_cond))) %>% 
  filter(Karyotype == "T21") %>% 
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")), # ensure factor levels in correct order
    has_cond = fct_relevel(has_cond, c("FALSE", "TRUE")) # ensure factor levels in correct order
    ) %>% 
  # 0 replacement with NA
  mutate(
    proportion_subset = if_else(proportion_subset == 0, NA_real_, proportion_subset)
  ) %>%
  # Outlier removal
  group_by(uniqueID, condition, has_cond) %>%  
  mutate(
    logit_proportion_subset = gtools::logit(proportion_subset),
    extreme_proportion_subset = rstatix::is_extreme(logit_proportion_subset),
  ) %>% 
  filter(extreme_proportion_subset == FALSE 
  ) %>% 
  ungroup()
#
regression_prep_dat_final <- regression_prep_dat %>%  
  filter(!is.na(proportion_subset)) %>% 
  group_by(condition) %>%
  # remove conditions which have only TRUE or FALSE
  mutate(has_cond_levels = has_cond %>% fct_drop() %>% levels() %>% length()) %>% 
  filter(has_cond_levels > 1 ) %>% # keeps 2566695 rows
  ungroup() %>% 
  # condition × features × has_cond grouping
  group_by(uniqueID, condition, has_cond) %>%
  add_count(condition, has_cond) %>%
  rename(min_count = n) %>% 
  group_by(uniqueID, condition) %>% 
  # remove any features which do not have at least 10 counts per TRUE or FALSE in each "feature + condition" combination
  filter(!any(min_count < 10)) %>% # keeps 1275571 rows
  ungroup() %>% 
  # now check balance of categorical covariates
  group_by(uniqueID, condition, has_cond) %>% 
  mutate( 
    Sex_levels = Sex %>% fct_drop() %>% levels() %>% length()
  ) %>%
  group_by(uniqueID, condition) %>% 
  # remove sex-specific conditions (require > 1 level, or lm() gives error)
  filter(Sex_levels > 1 ) %>% 
  add_count(Sex, name = "sex_count") %>% 
  group_by(uniqueID, condition) %>% 
  # remove any features from a particular condition if it does not at least 3 samples per sex in each "feature + condition" combination
  filter(!any(sex_count < 3)) %>% # require at least 3 females and 3 males in each category
  ungroup()
#
regression_prep_dat_final %>% 
  distinct(condition) 
#
regression_prep_dat_final %>% 
  distinct(LabID) 
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
all_cond_cytof_betareg_data_final.vec <- character()

for (cond in conditions$condition) {
  print(paste("Processing condition:", cond))
  
  # Step 1: filter to one cond and nest the data
  T21_betareg_data_final.tmp <- NULL # initializing
  T21_betareg_data_final.tmp <- regression_prep_dat_final %>%
    filter(condition == cond) %>%
    select(LabID, uniqueID, cluster, proportion_subset, Sex, Age, Karyotype, Sample_source_code, condition, has_cond) %>%
    nest(data = c(LabID, proportion_subset, Sex, Age, Karyotype, Sample_source_code, has_cond))
  
  # Debug: Check the condition data
  if (nrow(T21_betareg_data_final.tmp) == 0) {
    print(paste("No data for condition:", cond))
    next
  }
  
  # Step 2: Multi Sex + Age + Source [PREFERRED MODEL]
  regressions_multi_cytof_AgeSexSource.tmp <- NULL # initializing
  regressions_multi_cytof_AgeSexSource.tmp <- T21_betareg_data_final.tmp %>%
    mutate(
      fit = map(data, ~ betareg(proportion_subset ~ has_cond + Sex + Age + Sample_source_code, data = .x, link = "logit")),
      tidied = map(fit, broom::tidy), # see ?tidy.betareg
      glanced = map(fit, broom::glance), # see ?glance.betareg
      augmented = map(fit, broom::augment), # see ?augment.betareg
      vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
    )
  
  # step 3 Extract results for multivar model
  name.tmp <- NULL # initializing
  name.tmp <- paste0(cond)
  betareg_cytof_multi_ageSexSource_results.tmp <- regressions_multi_cytof_AgeSexSource.tmp %>%
    unnest(tidied) %>%
    select(uniqueID, cluster, condition, component, term, estimate, p.value) %>%
    filter(str_detect(term, "has_cond")) |>
    group_by(uniqueID) %>%
    summarize(
      uniqueID = first(uniqueID),
      level = first(level),
      cluster = first(cluster),
      pval = nth(p.value, n = 1),
      FoldChange = exp(nth(estimate, n = 1))
    ) %>%
    ungroup() %>%
    group_by(level) %>%
    arrange(level, pval) %>%
    filter(!str_detect(cluster, "Trash")) %>% # remove excluded clusters
    # group_by(condition) %>% 
    mutate(
      BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))
    ) %>%
    ungroup() %>%
    select(uniqueID, level, cluster, FoldChange, pval, BHadj_pval) %>%
    mutate(Model = "betareg(proportion ~ has_cond + Sex + Age + Sample_source_code)")
  
  betareg_cytof_multi_ageSexSource_results.tmp %>% assign(name.tmp, ., pos=1)
  all_cond_cytof_betareg_data_final.vec <- c(all_cond_cytof_betareg_data_final.vec, name.tmp)
}
#

## 3.2 Gather all results ----
all_cond_cytof_betareg_data_final_results <- mget(all_cond_cytof_betareg_data_final.vec, envir = .GlobalEnv) %>%
  bind_rows(.id = "condition")
#
all_cond_cytof_betareg_data_final_results %>% arrange(BHadj_pval) %>%
  filter(BHadj_pval < 0.1) %>%
  count(condition) %>%
  arrange(-n)
#

## 3.3  Export the results ----
all_cond_cytof_betareg_data_final_results %>%
  write_tsv(file = here("results", paste0(out_file_prefix, "results_conditions_SexAgeSource", ".txt")))
#


# 4 Plots ----
## 4.1 Volcano plot ----
# Using Hidradenitis suppurativa as an example here
v_plot <- all_cond_cytof_betareg_data_final_results %>%
  filter(condition = "Hidradenitis suppurativa") %>% 
  volcano_plot_lab_lm(
    title= "Diff. abundance of cells in Hidradenitis suppurativa (cases vs. controls) - beta regression",
    subtitle = paste0( "T21s only; AgeSexSource adjusted model\n[Down: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange <1) %>% nrow(),
                       "; Up: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange >1) %>% nrow(), "]")
  ) 
#
v_plot
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "bm_volcano_plot_sig_cond", ".pdf")), device = cairo_pdf, width = 7, height = 6, units = "in")
#

## 4.2 Sina plots of top associations ----
# Adjust the dataset for covariates to align with the linear regression model
#
### 4.2.1 Adjust the CyTOF dataset ----
# dataset with only T21:
T21_betareg_data_final <- regression_prep_dat_final %>%
  select(LabID, uniqueID, cluster, proportion_subset, Sex, Age, Karyotype, Sample_source_code, condition, has_cond) %>% 
  nest(data = c(LabID, proportion_subset, Sex, Age, Karyotype, Sample_source_code, has_cond))
#
## Adj the CyTOF data 
T21_betareg_data_final_adj <- T21_betareg_data_final %>% 
  unnest(data) %>% 
  mutate(logit_proportion_subset = gtools::logit(proportion_subset)) %>% 
  # adjust() FUNCTION NO LONGER IN EFFECTSIZE - moved to datawizard
  nest(data = -c(uniqueID, cluster)) %>% 
  mutate(
    adjust_sex = map(data, ~ datawizard::adjust(data = .x, effect = "Sex", select = c("logit_proportion_subset"), keep_intercept = TRUE), multilevel = T),
    adjust_ageSex =  map(adjust_sex, ~ datawizard::adjust(data = .x, effect = "Age", select = c("logit_proportion_subset"), keep_intercept = TRUE)),
    adjust_ageSexSource = map(adjust_ageSex, ~ datawizard::adjust(data = .x, effect = "Sample_source_code", select = c("logit_proportion_subset"), keep_intercept = TRUE, multilevel = T)),
  ) %>% 
  select(uniqueID, cluster, adjust_ageSexSource) %>% 
  unnest(adjust_ageSexSource) %>% 
  mutate(proportion_subset = gtools::inv.logit(logit_proportion_subset)) 
#

### 4.2.2 sina plots ----
# Using Hidradenitis suppurativa as an example here
T21_betareg_data_final_adj %>% 
  filter(Karyotype == "T21") %>% # have to be only ind. with DS
  filter(Aptamer %in% c(all_cond_cytof_betareg_data_final_results %>%
                          filter(condition == "Hidradenitis suppurativa") %>%
                          filter(BHadj_pval < 0.1) %>%
                          mutate(cluster = fct_drop(cluster)) %>%
                          arrange(-FoldChange) %>% 
                          pull(uniqueID))) %>% 
  # to get the correct order of the features based on effect size and significance
  mutate(uniqueID = fct_relevel(uniqueID, c(all_cond_cytof_betareg_data_final_results %>%
                                            filter(condition == "Hidradenitis suppurativa") %>% 
                                            filter(BHadj_pval < 0.1) %>%
                                            mutate(cluster = fct_drop(cluster)) %>%
                                            arrange(-FoldChange) %>% 
                                            pull(uniqueID)))) %>% 
  # filtered to the condition of interest
  inner_join(htp_cooccurring_data %>% filter(condition == "Hidradenitis suppurativa")) %>% 
  # remove any participants without any status for the condition
  filter(!is.na(has_cond)) %>% 
  # extreme outliers are calculated per feature and per status of the condition
  ggroup_by(uniqueID, condition, has_cond) %>%  
  mutate(
    logit_proportion_subset = gtools::logit(proportion_subset),
    extreme_proportion_subset = rstatix::is_extreme(logit_proportion_subset),
  ) %>% 
  filter(extreme_proportion_subset == FALSE) %>% 
  ungroup() %>% 
  # relevel the status of condition
  mutate(has_cond = fct_relevel(has_cond, c("FALSE", "TRUE"))) %>% 
  ggplot(aes(x = has_cond,  y = gtools::logit(proportion_subset), color = has_cond)) +
  scale_color_manual(values = c("FALSE" = "#4CAF50", "TRUE" = "#E04B4B")) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=1), breaks = seq(-12, 12, by = 1)) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  labs(title = "Hidradenitis suppurativa",
       subtitle = paste("Significant associations\nAgeSexSource adj\nextreme outliers are removed"),
       y = "% among subset") +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    aspect.ratio = 1.3,
  ) +
  facet_wrap(~uniqueID, ncol = 5, scales = "free_y")
#
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "Sig_CD45_logit_adj_hiradentitis", ".pdf")), device = cairo_pdf, width = 10, height = 4, units = "in")


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
