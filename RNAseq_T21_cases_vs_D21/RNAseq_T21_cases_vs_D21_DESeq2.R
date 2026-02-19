################################################
# Title: Differential abundance analysis of whole-blood transcriptome in T21 cases per condition vs. D21
# Author(s):
#   - Srija Chillamcherla
#   - Matthew Galbraith
# affiliation(s):
#   - Linda Crnic Institute for Down syndrome
#   - University of Colorado Anschutz
################################################

### Summary:  
# DESeq2 modelling for differential abundance of whole-blood bulk transcriptome genes 
# (RNA-seq) from individuals with (T21) and without (D21) Down syndrome.
# See README.md for more details
# 

### Data type(s):
# The Human Trisome Project (HTP) is a large-scale longitudinal, multi-omics research initiative 
# generating multi-omics data to advance understanding of Down syndrome (Trisomy 21) 
#   A. HTP sample meta data 
#      - HTP_Metadata_v0.5_Synapse.txt
#   B. HTP co-occuring conditions data
#       - HTP_Cooccuring_conditions_v0.5_Synapse.txt
#   C. HTP whole blood transcriptome data
#      - HTP_WholeBlood_RNAseq_Counts_Synapse.txt.gz
#      - HTP_WholeBlood_RNAseq_FPKMs_Synapse.txt.gz
#

### Workflow:
#   Step 1 - Read in counts data for all samples
#   Step 2 - Filter by minimum cpm  
#          - Groups and/or covariates setup
#          - Generate DESeqDataSet(s)
#          - Run DESeq2 analysis
#          - Get summary data for all samples
#   Step 3 - Model results
#   Step 4 - Visualise results
# 

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
library("edgeR") # required for cpm function
library("DESeq2") # Performs differential expression analysis of count data 
library("apeglm") # Implements adaptive shrinkage of effect sizes for improved log2 fold-change estimates
library("RColorBrewer") # Provides visually appealing color palettes for data visualization.
library("ggplot2") # required for current sample-sample distance heatmap
library("readxl") # reading Excel files
library("openxlsx") # required for exporting results as Excel workbooks
library("limma") # adjust for batch effects in the dataset
library("genefilter") # used for function rowVars
library("ggforce") # used for sina plots
library("tidyverse") # required for ggplot2, dplyr etc
library("ggrepel") # required for using geom_text and geom_text_repel() to make sample labels 
library("circlize") # color scale generation
library("BiocParallel") # enables mutli-cpu for some of DEseq2 functions
library("dplyr") # data manipulation: filtering, summarizing, and transforming data frames.
library("ggplot2") # for creating versatile and layered plots
library("ggforce") # required for zooming and sina
library("purrr") # Functional programming tools for iteration and mapping over lists and vectors
library("rstatix") # grid – Provides low-level graphics functionality for fine control of plots in R.
library("grid")
library("conflicted") # force all conflicts to become errors
ncores <- parallel::detectCores() - 1
# register(MulticoreParam(workers = ncores)) # enables mutli-cpu for some of DEseq2 functions; can set number of workers - default is all cores; 
conflicts_prefer( # declare preferences in cases of conflict
  dplyr::filter,
  dplyr::count,
  dplyr::select,
  dplyr::rename,
  base::paste,
  matrixStats::rowVars
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
# whole blood transcriptome data
counts_file <- here("data", "HTP_WholeBlood_RNAseq_Counts_Synapse.txt.gz")
counts_format <- "HTSeq" # featureCounts or HTSeq or custom
rpkms_file <- here("data", "HTP_WholeBlood_RNAseq_FPKMs_Synapse.txt.gz")
gene_anno_file <- here("data", "gene_annotation_gencode.v33.basic.txt")
#
min_cpm <- 0.5 # used for low count filtering; default is 0.5
predictor <- "Karyotype" # set the predictor for the helper function
min_samples <- "auto" # used for low count filtering; use a number, "all", or "auto" (sets to half number of samples)
#
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e", "FALSE" = "#4CAF50", "TRUE" = "#E04B4B")
out_file_prefix <- "RNAseq_T21_cases_vs_D21_DESeq2.R_" # should match this script title
# End required parameters ###
source(here("helper_functions_DESeq.R")) # load helper functions
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
  select(-c("Data_contact", "Date_exported", "Script")) %>% 
  rename(Sampleid = LabID)
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
# The analysis is between D21s (with or without any co-occurring conditions) and
# T21s with co-occurring conditions
# Co-occurring conditions data should be filtered accordingly
htp_cooccurring_data <- htp_cooccurring_data %>% 
  mutate(keep = case_when(
    Karyotype == "T21" & has_cond == TRUE ~ TRUE,
    Karyotype == "T21" & is.na(has_cond) ~ FALSE,
    Karyotype == "Control" ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  filter(keep == TRUE)
#

## 1.3 Read in counts data ----
counts_data_long <- counts_file %>% 
  read_tsv() %>% 
  rename(
    Sampleid = LabID,
    raw_count = Value
  )
# convert the counts data to wide format
counts_data <- counts_data_long %>% 
  pivot_wider(id_cols = EnsemblID, names_from = Sampleid, values_from = raw_count)
#

## 1.4 Read rpkms data ---
rpkms_data <- rpkms_file %>% read_tsv() %>% 
  mutate(log2_value = log2(Value)) %>% 
  rename(Sampleid = LabID)
#
# inspect
counts_data_long %>% distinct(Sampleid) 
counts_data_long %>% distinct(EnsemblID, Gene_name) 
rpkms_data %>% distinct(Sampleid) 
#
#subset meta data to contain only the Sampleid of counts_data
meta_data <- htp_cooccurring_data %>% 
  inner_join(meta_data) %>% 
  inner_join(counts_data_long %>% 
               distinct(Sampleid))
#

## 1.5 Read gene annotation ----
gene_anno <- gene_anno_file %>% 
  read_tsv()
#


# 2 DESeq2 for all the conditions -------
## 2.1 Perform DESeq2 for all the conds using a for-loop ----
# Convert the age into a z-score
meta_data_updated <- meta_data %>% 
  mutate(
    Age_zscore = (Age - mean(Age)) / sd(Age), # convert to Z-scores to scale and center
  ) %>% 
  mutate(Sex = fct_relevel(Sex, c("Female", "Male")))
#
#Get the unique conditions (239 in total)
unique_conditions <- htp_cooccurring_data %>% distinct(condition)
#
# Create an empty vector to store the results
all_cond_transcriptome_DESeq2_data_final.vec <- character()
#
# Some conditions fail this step due to low detection of certain genes.
# Using tryCatch allows us to catch these errors without stopping the loop,
# so processing can continue for the remaining conditions.
# File to store conditions that fail the APEGLM log2FC shrinkage step
#
error_log_file <- "failed_conditions.txt"
#
# For loop to iterate over each condition
for (cond in unique(unique_conditions$condition)) {
  
  # Log the condition being processed
  cat("Running analysis for ", cond, "\n")
  
  # update metadata used in this analysis
  condition_meta_data <- meta_data_updated %>%
    filter(condition == cond) %>%
    mutate(Sex = fct_relevel(Sex, c("Female", "Male"))) %>%
    mutate(Karyotype = fct_relevel(Karyotype, c("Control", "T21")))
  
  counts_data <- counts_data_long %>%
    pivot_wider(names_from = Sampleid, values_from = raw_count) %>%
    rename(Geneid = EnsemblID) %>% 
    select(Geneid, condition_meta_data %>% pull(Sampleid))
  
  # Filter by minimum counts per million
  min_samples = "auto"
  cat("min_samples is set to: ", min_samples, "\n")
  
  if (min_samples == "all") {
    min_samples = ncol(counts_data) - 1
  } else if (min_samples == "auto") {
    min_samples = (ncol(counts_data) - 1) / 2
  }
  
  before <- counts_data %>%
    transmute(Geneid = Geneid, row_sum = rowSums(select(., -Geneid))) %>%
    filter(row_sum > 0) %>%
    nrow()
  
  cpm_data <- counts_data %>%
    column_to_rownames("Geneid") %>%
    cpm()
  
  keep <- cpm_data %>%
    as_tibble(rownames = "Geneid") %>%
    pivot_longer(-Geneid, names_to = "Sampleid", values_to = "cpm") %>%
    mutate(cpm > min_cpm) %>%
    filter(`cpm > min_cpm` == TRUE) %>%
    dplyr::count(Geneid) %>%
    filter(n >= min_samples)
  
  counts_filtered <- counts_data %>%
    filter(Geneid %in% keep$Geneid)
  
  # Groups and/or Covariates setup
  groups <- condition_meta_data %>%
    select(Sampleid, condition, Karyotype, Sex, Age_zscore, Sample_source_code)
  
  # Generate DESeqDataSet object(s)
  multivar_formula <- as.formula(paste0("~", "Karyotype", " + Sex + Age_zscore + Sample_source_code"))
  
  # Multivariable model
  dds_multi <- DESeqDataSetFromMatrix(
    countData = counts_filtered %>%
      select(Geneid, groups %>% pull(Sampleid)) %>%
      column_to_rownames("Geneid"),
    colData = groups,
    design = multivar_formula
  )
  
  # Run multivariable model(s)
  dds_multi <- DESeq(dds_multi, parallel = TRUE)
  
  # TryCatch block to handle potential errors
  tryCatch({
    # Extract results
    
    # look at the function in helper script and modify accordingly
    res_multi <- dds_multi %>%
      get_results_tbl(
        contrast = c("Karyotype", "T21", "Control"), # relevel to ensure the factor
        shrink_type = "apeglm"
      )
    
    # Store results and track condition
    name.tmp <- paste0(cond)
    assign(name.tmp, res_multi, pos = 1)
    all_cond_transcriptome_DESeq2_data_final.vec <- c(all_cond_transcriptome_DESeq2_data_final.vec, name.tmp)
    
  }, error = function(e) {
    # If an error occurs, write the condition to the error log
    cat("Error for condition: ", cond, "\n", file = error_log_file, append = TRUE)
    cat("Error message: ", e$message, "\n", file = error_log_file, append = TRUE)
  })
  
}
#

## 2.3 Gather all results ----
all_cond_transcriptome_DESeq2_data_final_results <- mget(all_cond_transcriptome_DESeq2_data_final.vec, envir = .GlobalEnv) %>%
  bind_rows(.id = "condition")
#
all_cond_transcriptome_DESeq2_data_final_results %>%
  distinct(condition) 
#
all_cond_transcriptome_DESeq2_data_final_results %>% arrange(padj) %>%
  filter(padj < 0.1) %>%
  count(condition) %>%
  arrange(-n)
# 

## 2.4 Export the results ----
all_cond_transcriptome_DESeq2_data_final_results %>%
  write_tsv(file = here("results", paste0(out_file_prefix, "results_conditions_SexAgeSource", ".txt")))

# 3 Plots --------
## 3.1 Volcano plot ----
# Using Cognitive decline as an example here
v_plot <- all_cond_transcriptome_DESeq2_data_final_results %>%
  filter(condition == "Cognitive decline") %>% # change to any other condition(s)
  volcano_plot_lab(
    title="Diff. abundance of genes in Cognitive decline (T21 cases vs. D21) - DESeq2",
    subtitle = paste0( "AgeSexSource adjusted model\n[Down: ",(.) %>% filter(padj < 0.1 & FoldChange <1) %>% nrow(),
                       "; Up: ",(.) %>% filter(padj < 0.1 & FoldChange >1) %>% nrow(), "]")
  ) 
#
v_plot
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "lm_volcano_plot_sig_cond", ".pdf")), device = cairo_pdf, width = 7, height = 6, units = "in")
#

## 3.2 Sina plots of top associations ----
# Use covariates adjusted dataset to align with the DESeq2 
### 3.2.1 Adjust the RPKMS data ----
rpkms_sample_data <- rpkms_data %>%
  filter(Sampleid %in% c(meta_data$Sampleid)) %>% 
  select(-Value, -log2_value) %>%
  distinct()
rpkms_unadj_data <- rpkms_data %>%
  filter(Sampleid %in% c(meta_data$Sampleid)) %>% 
  select(Sampleid, EnsemblID, log2_value) %>%
  spread(key = Sampleid, value = log2_value) %>% # need to log2 transform for batch correction
  column_to_rownames(var = "EnsemblID")
rpkms_sex_vec <- rpkms_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name=NULL, value="Sampleid") %>%
  inner_join(
    meta_data %>% select(Sampleid, Sex)
  ) %>%
  distinct() %>%
  pull(Sex)
rpkms_age_vec <- rpkms_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name=NULL, value="Sampleid") %>%
  inner_join(
    meta_data %>% select(Sampleid, Age)
  ) %>%
  distinct() %>%
  pull(Age)
rpkms_source_vec <- rpkms_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "Sampleid") %>%
  inner_join(
    meta_data %>% select(Sampleid, Sample_source_code)
  ) %>%
  distinct() %>%
  pull(Sample_source_code)
rpkms_design <- rpkms_unadj_data %>% # get batch information, ensuring order matches unadj_data
  colnames() %>%
  enframe(name = NULL, value = "Sampleid") %>%
  inner_join(
    meta_data %>% select(Sampleid, Karyotype)
  ) %>%
  distinct(Sampleid, Karyotype) %>%
  model.matrix(~Karyotype, data = .)
#
rpkms_SourceSexAge_adj <- rpkms_unadj_data %>%
  limma::removeBatchEffect(
    batch = rpkms_source_vec, # categorical as batch
    batch2 = rpkms_sex_vec, # categorical as batch
    covariates = rpkms_age_vec, # numeric as co-variate
    design = rpkms_design
  ) %>%
  as_tibble(rownames = "EnsemblID") %>%  # convert back to tibble
  gather(-EnsemblID, key = Sampleid, value = Value) %>%
  mutate(value_adj = 2^Value) %>%
  rename(log2_value_adj = Value) %>% # rename
  inner_join(rpkms_sample_data)
#
# join it with meta data
rpkms_SourceSexAge_adj <- rpkms_SourceSexAge_adj %>% 
  inner_join(meta_data)
#

# sina plot 
# Use covariates adjusted dataset to align with the DESeq2 
# Using Cognitive decline as an example here
rpkms_SourceSexAge_adj %>% 
  filter(Gene_name %in% c(all_cond_transcriptome_DESeq2_data_final_results %>%
                            filter(condition == "Cognitive decline") %>%
                            filter(padj < 0.1) %>%
                            slice_max(order_by = log2FoldChange_adj, n = 5) %>%
                            pull(Gene_name))) %>%
  mutate(Gene_name = fct_relevel(Gene_name, c(all_cond_transcriptome_DESeq2_data_final_results %>%
                                               filter(condition == "Cognitive decline") %>%
                                               filter(padj < 0.1) %>%
                                               slice_max(order_by = log2FoldChange_adj, n = 5) %>%
                                               pull(Gene_name)))) %>%
  # filtered to the condition of interest
  inner_join(htp_cooccurring_data %>%  filter(condition == "Cognitive decline")) %>% 
  # extreme outliers are calculated per feature and per status of the condition
  group_by(Gene_name, Karyotype) %>%
  mutate(extreme_score = rstatix::is_extreme(log2_value_adj)) %>%
  ungroup() %>%
  filter(extreme_score == FALSE) %>% 
  # relevling the karyotype 
  mutate(Karyotype = fct_relevel(Karyotype, c("Control", "T21"))) %>% 
  ggplot(aes(x = Karyotype,  y = log2_value_adj, color = Karyotype)) +
  scale_color_manual(values = standard_colors) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  labs(title = "Cognitive decline",
       subtitle = paste("Significant associations\nAgeSexSource adj\nextreme outliers are removed"),
       y = "log2(RPKM)") +
  theme(
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    legend.position = "right",
    aspect.ratio = 1.3,
  ) +
  facet_wrap(~Gene_name, ncol = 5, scales = "free_y")
#
# save the plot
ggsave(filename = here("plots", paste0(out_file_prefix, "cognitive_decline_sig_sina_top5", ".pdf")), device = cairo_pdf, width = 10, height = 4, units = "in")



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
