################################################
# Title: Differential abundance analysis of plasma proteins by co-occurring condition status (cases vs. controls) in individuals with Down syndrome (T21)
# Author(s):
#   - Srija Chillamcherla
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# Linear regression modelling for differential abundance of plasma proteins and
# measured by the SOMAscan assay from individuals with Down syndrome (T21)
# See README.md for more details
# 

### Data type(s):
# The Human Trisome Project (HTP) is a large-scale longitudinal, multi-omics research initiative 
# generating multi-omics data to advance understanding of Down syndrome (Trisomy 21) 
#   A. HTP sample meta data 
#      - HTP_Metadata_v0.5_Synapse.txt
#   B. HTP co-occuring conditions data
#       - HTP_Cooccuring_conditions_v0.5_Synapse.txt
#   B. HTP SOMAscan data
#      - HTP_SOMAscan_Proteomics_Synapse.txt
#

### Workflow:
#   Step 1 - Read in meta data and sample data
#   Step 2 - Linear regression model (lm) setup and assessment
#   Step 3 - Model results
#   Step 4 - Visualize lm results
#   Step 5 - Gene set enrichment analysis (GSEA)
#   Step 6 - Visualize GSEA results


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
library("fgsea") # to run GSEA 
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

## 0.2 Set file name parameters ----
# HTP datasets used in this study can be obtained from the 
# associated Synapse repository (further details in README.md)
# Download each dataset to the appropriate /data/ directories within R project.
# Input data files
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_Synapse.txt")
#
htp_cooccurring_data_file <- here("data", "HTP_Cooccuring_conditions_v0.5_Synapse.txt")
# SOMAscan data file
soma_data_file <-here("data", "HTP_SOMAscan_Proteomics_Synapse.txt")
#
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e", "FALSE" = "#4CAF50", "TRUE" = "#E04B4B")
out_file_prefix <- "SOMAscan_T21_cases_vs_controls_lm.R_" # should match this script title
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

## 1.3 Read in SOMAacan data ----
soma <- soma_data_file %>% 
  read_tsv() %>% 
  mutate(log2_value = log2(Value)) 
#
# inspect
soma %>% distinct(LabID) # 419 samples
soma %>% distinct(Aptamer) # 4632 measurements each for 419 samples = 1940808
# 
# check if any samples are not in the metadata
soma %>% anti_join(meta_data)
#
soma %>% 
  inner_join(meta_data) %>% 
  group_by(Karyotype) %>%
  skimr::skim()
#
# subset coc data to have labIDs in soma
htp_cooccurring_data <- htp_cooccurring_data %>% 
  semi_join(soma, by = "LabID") 
#


# 2 Linear regression: Data preparation -----
# All the covariates should be adjusted - Age, Sex and Sample source
# Get the unique conditions (239 in total)
unique_conditions <- htp_cooccurring_data %>% distinct(condition)
#
# Create an empty vector to store the results
all_cond_aptamers_lm_data_final.vec <- character()
#
## For loop for all co-occurring conditions vs SomaScan aptamers ----
for (cond in unique_conditions$condition) {
  print(paste("Processing condition:", cond))
  
  # Filter coc_data for the current condition and remove NAs from has_cond
  condition_data <- htp_cooccurring_data %>%
    filter(condition == cond) %>%
    filter(!is.na(has_cond))
  
  # Join with soma and perform filtering steps
  # Filtering is done in two stages/parts
  # First stage: the extreme outliers are removed PER FEATURE PER CONDITION 
  # Second stage: all the features which do not pass the minimum requirements are removed
  regression_prep_dat <- soma %>%
    inner_join(condition_data) %>%
    filter(Karyotype == "T21") %>%
    mutate(Sex = fct_relevel(Sex, c("Female", "Male")),  # ensure factor levels in correct order
           has_cond = fct_relevel(has_cond, c("FALSE", "TRUE"))) %>%  # ensure factor levels in correct order
    group_by(Aptamer, condition, has_cond) %>%
    mutate(extrem_aptamer_conc = rstatix::is_extreme(log2(Value))) %>%
    ungroup() %>%
    filter(extrem_aptamer_conc != TRUE)
  
  # Apply more filtering
  regression_prep_dat_final <- regression_prep_dat %>% 
    # remove conditions which have only TRUE or FALSE
    mutate(has_cond_levels = has_cond %>% fct_drop() %>% levels() %>% length()) %>%
    filter(has_cond_levels > 1) %>%
    # condition × features × has_cond grouping
    group_by(Aptamer, condition, has_cond) %>%
    add_count(condition, has_cond) %>%
    rename(min_count = n) %>%
    group_by(Aptamer, condition) %>% 
    # remove any features which do not have at least 10 counts per TRUE or FALSE in each "feature + condition" combination
    filter(!any(min_count < 10)) %>%
    ungroup() %>%
    # now check balance of categorical covariates
    group_by(Aptamer, condition, has_cond) %>%
    mutate(Sex_levels = Sex %>% fct_drop() %>% levels() %>% length()) %>%
    # remove sex-specific conditions (require > 1 level, or lm()/DESeq2 gives error)
    group_by(Aptamer, condition) %>%  
    # remove any features from a particular condition if it does not at least 1 sex level in each "feature + condition" combination
    filter(Sex_levels > 1) %>%
    add_count(Sex, name = "sex_count") %>%
    group_by(Aptamer, condition) %>%  
    # remove any features from a particular condition if it does not at least 3 samples per sex in each "feature + condition" combination
    filter(!any(sex_count < 3)) %>%
    ungroup() %>%
    # log transform the quantitative measurement
    mutate(log2_aptamer_conc = log2(Value))
  
  # If coc_prep_dat_final is empty, skip this condition and move to the next
  if (nrow(regression_prep_dat_final) == 0) {
    print(paste("No data for condition:", cond, "- skipping"))
    next
  }
  
  # Step 1: filter to one condition and nest the data
  T21_lm_data_final.tmp <- regression_prep_dat_final %>%
    select(LabID, Age, Sex, Sample_source_code, Aptamer, GeneSymbol, aptamer_conc = Value, log2_aptamer_conc, condition, has_cond) %>%
    nest(data = c(LabID, Age, Sex, Sample_source_code, aptamer_conc, log2_aptamer_conc, has_cond))
  
  # Step 2: Multi Sex + Age + Source (preferred model)
  regressions_multi_aptamer_AgeSexSource.tmp <- T21_lm_data_final.tmp %>%
    mutate(
      fit = map(data, ~ lm(log2_aptamer_conc ~ has_cond + Age + Sex + Sample_source_code, data = .x)),
      tidied = map(fit, broom::tidy),
      glanced = map(fit, broom::glance),
      augmented = map(fit, broom::augment),
      vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
    )
  
  # Step 3: Extract results for multivariable model
  name.tmp <- paste0(cond)
  lm_aptamer_multi_ageSexSource_results.tmp <- regressions_multi_aptamer_AgeSexSource.tmp %>%
    unnest(tidied) %>%
    filter(str_detect(term, "has_cond")) %>%
    select(Aptamer, GeneSymbol, condition, term, log2FoldChange = estimate, pval = p.value) %>%
    mutate(FoldChange = 2^log2FoldChange) %>%
    arrange(pval) %>%
    # group_by(condition) %>% 
    mutate(BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))) %>%
    # ungroup() %>% 
    select(Aptamer, GeneSymbol, condition, FoldChange, log2FoldChange, pval, BHadj_pval, everything())
  
  # Store results and track condition
  lm_aptamer_multi_ageSexSource_results.tmp %>% assign(name.tmp, ., pos = 1)
  all_cond_aptamers_lm_data_final.vec <- c(all_cond_aptamers_lm_data_final.vec, name.tmp)
}
#
# Gather all results
all_cond_aptamers_lm_data_final_results <- mget(all_cond_aptamers_lm_data_final.vec, envir = .GlobalEnv) %>%
  bind_rows(.id = "condition")
#
all_cond_aptamers_lm_data_final_results %>%
  distinct(condition) 
#
all_cond_aptamers_lm_data_final_results %>% arrange(BHadj_pval) %>%
  filter(BHadj_pval < 0.1) %>%
  count(condition) %>%
  arrange(-n)
#

## 3.3 Export the results ----
all_cond_aptamers_lm_data_final_results %>%
  write_tsv(file = here("results", paste0(out_file_prefix, "results_conditions_SexAgeSource", ".txt")))
# 


# 3 Plots ----
## 3.1 Volcano plot ----
# Using Obstructive sleep apnea as an example here
v_plot <- all_cond_aptamers_lm_data_final_results %>%
  filter(condition == "Obstructive sleep apnea") %>% # change to any other condition(s)
  volcano_plot_lab_lm(
    title="Diff. abundance of plasma proteins in Obstructive sleep apnea (cases vs. controls) - linear regression",
    subtitle = paste0( "T21s only; AgeSexSource adjusted model\n[Down: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange <1) %>% nrow(),
                       "; Up: ",(.) %>% filter(BHadj_pval < 0.1 & FoldChange >1) %>% nrow(), "]")
  ) 
#
v_plot
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "lm_volcano_plot_sig_cond", ".pdf")), device = cairo_pdf, width = 7, height = 6, units = "in")
#

## 3.2 Sina plots of top associations ----
# Adjust the dataset for covariates to align with the linear regression model
#
### 3.2.1 Adjust the SOMAscan dataset ----
soma_sample_data <- soma %>%
  select(-Value) %>%
  distinct()
soma_unadj_data <- soma %>%
  select(LabID, Aptamer, log2_value) %>%
  spread(key = LabID, value = log2_value) %>% # need to log2 transform for batch correction
  column_to_rownames(var = "Aptamer")
soma_sex_vec <- soma_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name=NULL, value="LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Sex)
  ) %>%
  # distinct() %>%
  pull(Sex)
soma_age_vec <- soma_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name=NULL, value="LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Age)
  ) %>%
  # distinct() %>%
  pull(Age)
soma_source_vec <- soma_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Sample_source_code)
  ) %>%
  #distinct() %>%
  pull(Sample_source_code)
soma_design <- soma_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "LabID") %>%
  inner_join(
    meta_data %>% select(LabID, Karyotype)
  ) %>%
  distinct(LabID, Karyotype) %>%
  model.matrix(~Karyotype, data = .)
#
soma_SourceSexAge_adj <- soma_unadj_data %>%
  limma::removeBatchEffect(
    batch = soma_source_vec, # categorical as batch
    batch2 = soma_sex_vec, # categorical as batch
    covariates = soma_age_vec, #numeric as co variate
    design = soma_design
  ) %>%
  as_tibble(rownames = "Aptamer") %>%  # convert back to tibble
  gather(-Aptamer, key = LabID, value = Value) %>%
  mutate(value_adj = 2^Value) %>%
  rename(log2_value_adj = Value) %>% # rename
  inner_join(soma_sample_data)
#
# join with meta data
soma_SourceSexAge_adj <- soma_SourceSexAge_adj %>% 
  inner_join(meta_data)
#

### 3.2.2 sina plots ----
# Using Obstructive sleep apnea as an example here
soma_SourceSexAge_adj %>% 
  filter(Karyotype == "T21") %>% # have to be only ind. with DS
  filter(Aptamer %in% c(all_cond_aptamers_lm_data_final_results %>%
                          filter(condition == "Obstructive sleep apnea", BHadj_pval < 0.1) %>%
                          slice_max(order_by = log2FoldChange, n = 5) %>%
                          pull(Aptamer))) %>% 
  # to get the correct order of the features based on effect size and significance
  mutate(Aptamer = fct_relevel(Aptamer, c(all_cond_aptamers_lm_data_final_results %>%
                                            filter(condition == "Obstructive sleep apnea", BHadj_pval < 0.1) %>%
                                            slice_max(order_by = log2FoldChange, n = 5) %>%
                                            pull(Aptamer)))) %>% 
  # filtered to the condition of interest
  inner_join(htp_cooccurring_data %>%  filter(condition == "Obstructive sleep apnea")) %>% 
  # remove any participants without any status for the condition
  filter(!is.na(has_cond)) %>% 
  # extreme outliers are calculated per feature and per status of the condition
  group_by(Aptamer, has_cond) %>%
  mutate(extreme_score = rstatix::is_extreme(log2_value_adj)) %>%
  ungroup() %>%
  filter(extreme_score == FALSE) %>% 
  # relevling the status of condition 
  mutate(has_cond = fct_relevel(has_cond, c("FALSE", "TRUE"))) %>% 
  ggplot(aes(x = has_cond,  y = log2_value_adj, color = has_cond)) +
  scale_color_manual(values = c("FALSE" = "#4CAF50", "TRUE" = "#E04B4B")) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  labs(title = "Obstructive sleep apnea",
       subtitle = paste("Significant associations\nAgeSexSource adj\nextreme outliers are removed"),
       y = "log2_adj_val") +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    aspect.ratio = 1.3,
  ) +
  facet_wrap(~Aptamer, ncol = 5, scales = "free_y")
#
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "Obstructive_sleep_apnea_sig_sina_top5", ".pdf")), device = cairo_pdf, width = 10, height = 4, units = "in")


# 4 GSEA Hallmarks analysis ----
# Load the Hallmarks file
hallmarks <- here("GSEA/h.all.v7.4.symbols.gmt") %>%
  fgsea::gmtPathways(gmt.file = .)
#
# Create a function to run the GSEA on all the aptamers from linear model results
get_rnk <- function(query) {
  query %>% 
    get() %>% 
    arrange(-abs(FoldChange)) %>% # to ensure keeping "strongest" correlations where dup gene names in next step
    distinct(GeneSymbol, .keep_all = TRUE) %>% # where Gene_name is not distinct, this keeps the first row of values ( = strongest correlation from above)
    mutate(rank_metric = log2FoldChange * -log10(pval)) %>% # weighting fold-changes by pval => makes results more stable
    select(ID = GeneSymbol, r = rank_metric) %>%  # weighting fold-changes by pval
    #select(ID = GeneSymbol, r = FoldChange) %>% # use raw fold-changes # BHadj_pval
    arrange(-r) %>%
    tibble::deframe() # convert to named num vector
}
#
ranksList <- c("lm_T21_aptamer_Age_Sex_Sample_source_results")
#
names(ranksList) <- paste0(ranksList, "_ranks") # add names
map(ranksList, get_rnk) %>% # generate rnk vectors for each
  list2env(envir = .GlobalEnv) # save each rnk vector as object using names
unweighted_gsea_res_list <- map(names(ranksList), ~ run_fgsea2(geneset = hallmarks, ranks = get(.x), weighted = FALSE))
names(unweighted_gsea_res_list) <- ranksList
unweighted_gsea_res_tbl <- bind_rows(unweighted_gsea_res_list, .id = "Query") 
#
unweighted_gsea_results <- unweighted_gsea_res_tbl %>% 
  mutate(pathway = pathway %>%  str_remove("^HALLMARK_") %>% str_replace_all("_", " ") %>% str_to_title()) %>% 
  arrange(padj) %>% 
  mutate(leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ", "))) %>% 
  select(pathway, ES, NES, pval, padj, everything()) 
#
unweighted_gsea_results %>% 
  filter(padj < 0.1)
#
## 4.1 Export GSEA results ----
bind_rows(
  list(
    "unweighted" = unweighted_gsea_results 
  ),
  .id = "Comparison"
) %>%
  split(.$Comparison) %>%
  export_excel(filename = "GSEA_Hallmarks")
#

## 4.2 GSEA enrichment plots ----
# Using Other GI history as example here
plotEnrichment2(
  pathway = hallmarks$HALLMARK_ESTROGEN_RESPONSE_LATE,
  stats = `ranks_Obstructive sleep apnea_sourceSexAge_adj`,
  res = (unweighted_gsea_results %>% filter(uniqueID == "Obstructive sleep apnea")),
  title = "Co-occuring conditions TRUE vs. FALSE\nEstrogen response late\nObstructive sleep apnea"
)
#
ggsave(filename = here("plots", paste0(out_file_prefix, "enrichment_estro_Obstructive sleep apnea", ".pdf")), device = cairo_pdf, width = 5, height = 4, units = "in")


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
