################################################################################
# Title: Differential abundance analysis of T21 plasma metabolites vs. BMI
# Author(s):
#   - Zenitha Sundararajan
#   - Matthew Galbraith
# Affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################################################

### Summary:  
# Linear regression modelling for differential abundance of plasma metabolomics and lipodomics via 
# liquid chromatography-mass spectrometry (LC-MS) from individuals with Down syndrome (T21).
# See README.md for more details

### Data type(s):
#   A. Human Trisome Project (HTP) sample meta data 
#      - HTP_Metadata_v0.5_Synapse.txt
#   B. HTP Metabolomics data
#      - HTP_Plasma_Metabolomics_Synapse.txt
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
## initialize renv
# renv::init(bioconductor = TRUE) # run first time only
#
# To install the exact versions of all R packages but requires matching R version
# renv::restore()

## 0.1 Load required libraries ----
library("tidyverse") # data manipulation, visualization, and tidy workflows
library("broom") # for tidying models
library("tictoc") # timer
library("ggplot2") #data visualization
library("ggforce") # required for zooming and sina
library("rstatix") # required for statistical tests
library("conflicted") # force all conflicts to become errors
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::count,
  dplyr::rename,
  dplyr::bind_rows
)
library("here") # generates path to current project directory
#

## 0.2 Set required parameters ----
#
# Input data files
# These files need to be downloaded from ****
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_Synapse.txt")
# metabolomics data file
Metabolomics_data_file <- here("data", "HTP_Plasma_Metabolomics_Synapse.txt")
#
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e")
#
out_file_prefix <- "Metabolomics_T21_BMI_lm.R_" # should match this script title
# End required parameters ###
source(here("helper_functions.R")) # load helper functions
#


# 1 Read in meta data and sample data ----
## 1.1 Read in and inspect HTP metadata  -----
meta_data <- htp_meta_data_file %>% 
  read_tsv() %>% # 1722 obs/rows
  filter(Event_name != "Current") %>% 
  mutate(
    Sex = fct_relevel(Sex, "Female"),
    Karyotype = fct_relevel(Karyotype, "Control"),
    Sample_source = as_factor(Sample_source_code), # convert to factor - default is numerical order
    BMI = round(Weight_kg/((Height_cm/100) ** 2 ),2 )
  ) 
# 

## 1.2 Read in and inspect Metabolomics data  -----
metabolomics_data <- Metabolomics_data_file %>% 
  read_tsv()  
# inspect
metabolomics_data %>% distinct(LabID) # 419 samples
metabolomics_data %>% distinct(Analyte) # 174 measurements each for 419 samples = 72,906
#

## 1.3 Merge Metabolomics with meta data -----
metabolomics_data_meta_data <- metabolomics_data %>%
  inner_join(meta_data,
             by = c("LabID")) %>% # keeps 72,906 rows
  filter(!is.na(Weight_kg), !is.na(Height_cm)) # keeps 70,470 rows
#


# 2 Linear regression model setup and assessment ----
## 2.1 Prepare data for linear regression ----
regressions_data <- metabolomics_data_meta_data %>% # 70,470 rows
  filter(Karyotype == "T21") |> # keeps 54,636 rows
  select(ParticipantID, LabID, BMI, Sex, Age, Sample_source, Analyte, Value) %>% 
  mutate(
    Sex = fct_relevel(Sex, c("Female", "Male")), # ensure factor levels in correct order
  ) %>% 
  group_by(Analyte) %>% # group by analyte
  mutate(extreme_value = rstatix::is_extreme(log2(Value))
  ) %>% 
  ungroup() %>% 
  filter(extreme_value != TRUE) %>% # keeps 54,402 rows (extreme outliers are removed per feature)
  # check here for:
  # 1) a minimum number of samples per group (eg 5) and 
  # 2) that there are >1 levels for categorical variables of interest (prevents errors in regression step)
  group_by(Analyte) %>%  # CHECK CORRECT GROUPING
  mutate( # count number of levels for EACH categorical variable
    Sex_levels = Sex %>% fct_drop() %>% levels() %>% length()
  ) %>%
  filter(Sex_levels > 1) %>% # need to require >1 level for each categorical level or lm() gives error
  ungroup() %>%
  select(ParticipantID, LabID, BMI, Sex, Age, Sample_source, Analyte, Value) |> 
  nest(data = c(ParticipantID, LabID, BMI, Sex, Age, Sample_source, Value)) # nesting allows for easy testing of all features ~ at once
#

## 2.2 Linear regression: Multi-model: BMI + Age + Sex + Sample_source ----
tic("Running linear regressions for multi model with Sex + Age + Sample_source...")
regressions_multi_AgeSexSourceFixed <- regressions_data %>% 
  mutate(
    fit = map(data, ~ lm(log2(Value) ~ BMI +  Age + Sex + Sample_source, data = .x)),
    tidied = map(fit, broom::tidy), # see ?tidy.lm
    glanced = map(fit, broom::glance), # see ?glance.lm
    augmented = map(fit, broom::augment), # see ?augment.lm
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc() # ~ 1.351 sec elapsed; no warnings
#


# 3 Model results ----
## 3.1 Extract results from preferred model: Multi-model: BMI + Age + Sex + Sample_source ----
lm_T21_Analyte_BMI_Age_Sex_Sample_source_results <- regressions_multi_AgeSexSourceFixed |> 
  unnest(tidied) |> 
  filter(str_detect(term, "BMI")) |> # safer than summarize/first approach below
  select(Analyte, term, log2FoldChange = estimate, pval = p.value) |> 
  mutate(FoldChange = 2^log2FoldChange) |> 
  arrange(pval) %>% 
  mutate(BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))) %>% 
  select(
    Analyte,
    FoldChange,
    log2FoldChange,
    pval,
    BHadj_pval,
    everything()
  )
#
lm_T21_Analyte_BMI_Age_Sex_Sample_source_results %>% 
  filter(BHadj_pval < 0.1) # inspect the data
#

## 3.2 Export the results ----
lm_T21_Analyte_BMI_Age_Sex_Sample_source_results %>%
write_tsv(file = here("results", "HTP T21 bmi plasma metabolome linear regressions results.txt"))
#


# 4 Plots ----
## 4.1 Volcano plots -----
v_plot <- lm_T21_Analyte_BMI_Age_Sex_Sample_source_results %>% 
  volcano_plot_lab_lm(
    title = "Diff. abundance of plasma metabolites (T21 vs. BMI) - linear regression",
    subtitle = paste0("T21s only; AgeSexSource adjusted model\n[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]")
  )
v_plot
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "lm_volcano_plot_sig", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#

## 4.2 Scatter plots of top associations ----
# Adjust the dataset for covariates to align with the linear regression model
### 4.2.1 Adjust the Metabolomics dataset ----
metabolomics_sample_data <-  metabolomics_data_meta_data %>% #70,470 rows
  select(-Value) %>%
  distinct()
#
metabolomics_unadj_data <- metabolomics_data_meta_data %>% #70,470 rows
  select(LabID, Analyte, Value) %>%
  mutate(Value = log2(Value)) %>% # need to log2 transform for batch correction
  pivot_wider(names_from = LabID, values_from = Value) %>%
  column_to_rownames(var = "Analyte")
#
metabolomics_source_vec <- metabolomics_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Sample_source)
  ) %>%
  pull(Sample_source)
#
metabolomics_Sex_vec <- metabolomics_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Sex)
  ) %>%
  pull(Sex)
#
metabolomics_Age_vec <- metabolomics_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Age) 
  ) %>%
  pull(Age)
#
metabolomics_design <- metabolomics_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Karyotype)
  ) %>%
  model.matrix(~ Karyotype, data = .)
#
metabolomics_SourceSexAge_adj <- metabolomics_unadj_data %>%
  limma::removeBatchEffect(
    batch = metabolomics_source_vec, # categorical as batch
    batch2 = metabolomics_Sex_vec, # categorical as batch
    covariates = metabolomics_Age_vec,  # numeric as co-variate
    design = metabolomics_design
  ) %>%
  as_tibble(rownames = "Analyte") %>%  # convert back to tibble
  gather(-Analyte, key = LabID, value = adjusted_conc) %>%
  mutate(adjusted_conc = 2^adjusted_conc) %>% # remove log2 transformation
  inner_join(metabolomics_sample_data) %>%
  rename(Abundance_adj = adjusted_conc)
#

### 4.2.2 Adjust the metadata for BMI ----
sample_data_BMI <- metabolomics_data_meta_data %>% # 70,470 rows
  distinct(ParticipantID, LabID, Karyotype, Sex, Age, BMI, Sample_source) %>% # Keeps 405
  select(-BMI)
#
unadj_data_BMI <-  metabolomics_data_meta_data %>% #70,470 rows
  distinct(LabID, BMI) %>% # Keeps 405
  mutate(feature = "BMI",BMI = log2(BMI)) %>% # need to log2 transform for batch correction
  pivot_wider(names_from = LabID, values_from = BMI) %>%
  column_to_rownames(var = "feature") #convert tibble to dataframe
#
source_vec_BMI <- unadj_data_BMI %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>% #build tibble
  inner_join(
    meta_data %>% select(LabID, Sample_source) 
  ) %>%
  pull(Sample_source)
#
Sex_vec_BMI <- unadj_data_BMI %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Sex) 
  ) %>%
  pull(Sex)
#
Age_vec_BMI <- unadj_data_BMI %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Age) 
  ) %>%
  pull(Age)
#
design_BMI <- unadj_data_BMI %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Karyotype)
  ) %>%
  model.matrix(~ Karyotype, data = .)
#
SourceSexAge_adj_BMI <- unadj_data_BMI %>%
  limma::removeBatchEffect(
    batch = source_vec_BMI,  # categorical as batch
    batch2 = Sex_vec_BMI,  # categorical as batch
    covariates = Age_vec_BMI, #numeric as co variate
    design = design_BMI
  ) %>%
  as_tibble(rownames = "BMI") %>%  # convert back to tibble
  gather(-BMI, key = LabID, value = adjusted_bmi) %>%
  mutate(adjusted_bmi = 2^adjusted_bmi) %>% # remove log2 transformation
  inner_join(sample_data_BMI) %>%
  rename(BMI_SourceSexAgeAdj = adjusted_bmi)
#


## 4.2.3 Get interesting/significant features ----
top_signif_by_FC <- bind_rows(
  lm_T21_Analyte_BMI_Age_Sex_Sample_source_results |> 
    filter(BHadj_pval<0.1) %>%  
    arrange(-log2FoldChange) %>%  
    slice_max(order_by = log2FoldChange, n = 5), # Upregulated
  lm_T21_Analyte_BMI_Age_Sex_Sample_source_results |> 
    filter(BHadj_pval<0.1) %>%  
    arrange(log2FoldChange) %>%  
    slice_min(order_by = log2FoldChange, n = 5) # Downregulated
)
#

## 4.2.4 Scatter plots ----
metabolomics_SourceSexAge_adj %>% 
  inner_join(SourceSexAge_adj_BMI %>% 
               select(ParticipantID,LabID,BMI_SourceSexAgeAdj) #Get adjusted bmi
  ) %>% #Keeps 70,470 rows
  filter(Karyotype == "T21") %>% #Keeps 15,834 rows
  filter(Analyte %in% top_signif_by_FC$Analyte) %>%  # filter to features of interest
  mutate(Analyte = fct_relevel(Analyte, top_signif_by_FC$Analyte)) %>%  # control plotting order
  group_by(Analyte) %>%  #Group by only analyte
  mutate(
    extreme_y = rstatix::is_extreme(log2(Abundance_adj))
  ) %>%
  filter(extreme_y == FALSE) %>% # remove extreme outliers
  mutate(density = getDenCols(BMI_SourceSexAgeAdj, Abundance_adj, transform = TRUE)) %>% 
  arrange(density) %>% 
  ungroup() %>%
  ggplot(aes(BMI_SourceSexAgeAdj, log2(Abundance_adj), color = density)) + #
  geom_point() +
  scale_color_viridis_c() +
  geom_smooth(method = "lm") +
  facet_wrap(~ Analyte , scales = "free_y", nrow = 2) +# facet per feature; each feature on it's own scale
  theme(aspect.ratio = 1) + # set fixed aspect ratio
  labs(title = "Top significant metabolites by fold-change: T21 vs. BMI",
       subtitle = paste("Significant associations\nAgeSexSource adj\nextreme outliers are removed"),
       y = "log2_adj_val",
       x = "log2_adj_bmi"
  )
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "scatter_T21_top_signif_by_FC", ".pdf")), device = cairo_pdf, width = 15, height = 5, units = "in")
#


########## End of Script #############

##### Save workspace    ----
################################################################################

save.image(file = here("rdata", paste0(out_file_prefix, ".RData")),
           compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)

# Reload the workspace
#load(here("rdata", paste0(out_file_prefix, ".RData")))

# session_info ----
date()
sessionInfo()
################################################################################