#' ---
#' title: "`r params$title`"
#' author: NULL
#' date: NULL
#' output:
#'   html_document:
#'     toc: true
#'     toc_float:
#'       collapsed: false
#'   fontsize: 11pt
#'   mainfont: arial
#' params:
#'   version: 0.1 # Do not edit except for new version
#'   title: Modelling of Karyotype vs. Frequencies/Proportions for CyTOF FlowSOM clusters - CD45+
#'   project: Human Trisome Project P4C cohort | P4C omics paper   
#'   date: !r strftime(Sys.Date(), format = "%B %d %Y")
#'   author:
#'     - Matthew Galbraith # can list additional authors
#'   email:
#'     - matthew.galbraith@cuanschutz.edu # can list additional emails
#'   affiliation: Department of Pharmacology & Linda Crnic Institute for Down Syndrome, University of Colorado Anschutz Medical Campus
#'   display_code: false #true/false to turn on/off display of R code as default
#'   collapse_output: true #true/false to turn on/off collapse of output as default
#' ---

#' ***
#' ##Project: `r params$project`  
#' Date:  `r params$date`  
#' Report version: `r params$version`   
#' Author(s): `r params$author`  
#' `r params$email`  
#' `r params$affiliation`  
#'   
#' ***  
#' ### Summary  
#' Analyzing differences in immune cell types in T21 vs Control individuals
#'   
#' **Data type(s):**  
#'   
#'   
#' A. CyTOF cluster frequencies (as %) for CD45+  
#'    (see )  
#'    
#'     
#'   
#+ include=FALSE
# **Workflow:**
# 
# 1. ...
# 
# **Comments:**
# 
# *   
#     + ...

#' ***

#+ include=FALSE
# TO PRODUCE REPORT: library("rmarkdown"); render("script_name")
############################################################
# Script original author: Matthew Galbraith                #
# version: 0.1  Date: 02_10_2022                           #
############################################################
# Change Log:
# v0.1
# Initial version (based on P4C_T21_CyTOF_Live_betareg_v0.1.R)
#


#+ knitr_setup, include=FALSE
knitr::opts_chunk$set(echo = params$display_code)
knitr::opts_chunk$set(collapse = params$collapse_output)
knitr::opts_chunk$set(tidy = TRUE)
knitr::opts_chunk$set(fig.align = "center")
knitr::opts_chunk$set(fig.height = 6)
knitr::opts_chunk$set(fig.width = 10)

#+ general_setup, include=FALSE, message=FALSE, warning=FALSE
# Clear workspace
if (exists("params")) { rm(list=ls()[! grepl("params", ls())]) } else { rm(list=ls()) } # This excludes "params" from the list for removal
# Load required libraries
library("knitr")
library("DT") # Used for sortable, searchable datatables in reports
library("readxl") # Used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse")
library("magrittr")
library("ggrepel") # required for labelling genes
library("ggforce") # required for zooming and sina
library("plotly") # required for interactive plots/IFNA1
library("tictoc") # timer
library("skimr")
library("janitor")
library("patchwork")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("betareg")
library("conflicted")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
library("here") # generates path to current project directory
# detach("package:here", unload=TRUE) # run this to reset here()
############################
### EDITABLE PARAMETERS ###
# Make sure these files exist!!
meta_data_file <- here("data", "HTP_CLEANED_02_2021_v0.5_MASTER_RecordID_vs_LabID.Labels.tsv")
#
CD45_meta30_clusters_data_file <- here("data", "P4C_CyTOF_Workflow_v0.1_meta30_sample_summary.txt")
CD45_merging1_clusters_data_file <- here("data", "P4C_CyTOF_Workflow_v0.1_merging1_sample_summary.txt")
CD45_merging2_clusters_data_file <- here("data", "P4C_CyTOF_Workflow_v0.1_merging2_sample_summary.txt")
CD45_merging3_clusters_data_file <- here("data", "P4C_CyTOF_Workflow_v0.1_merging3_sample_summary.txt")
# Bcells_clusters_data_file <- here("data", "")
# Tcells_clusters_data_file <- here("data", "")
#
# meta30_clusters_naming_merging_file <- here("data", "P4C_CyTOF_Workflow_v0.1_meta30_clusters_naming_merging.xlsx")
meta30_clusters_naming_merging_file <- here("data", "P4C_CyTOF_Workflow_v0.1_meta30_clusters_naming_merging_V3.xlsx")
#
CD45_meta30_clusters_summary_file <- here("data", "P4C_CyTOF_Workflow_v0.1_meta30_summary.txt")
CD45_merging1_clusters_summary_file <- here("data", "P4C_CyTOF_Workflow_v0.1_merging1_summary.txt")
CD45_merging2_clusters_summary_file <- here("data", "P4C_CyTOF_Workflow_v0.1_merging2_summary.txt")
CD45_merging3_clusters_summary_file <- here("data", "P4C_CyTOF_Workflow_v0.1_merging3_summary.txt")
#
CD45_cells_tSNE_file <- here("data", "P4C_CyTOF_Workflow_v0.1_tSNE_px440_500_cells_per_sample.txt.gz")
# Bcells_tSNE_file <- here("data", "")
# Tcells_tSNE_file <- here("data", "")
#
standard_colors <- c("Control" = "grey60", "T21" = "#009b4e")
#
out_file_prefix <- "P4C_T21_CyTOF_FlowSOM_betareg_v0.1_"
### END EDITABLE PARAMETERS ###
###############################
source(here("helper_functions.R")) # load helper functions


#' ***
#' **File locations and variables:**
#+ report_variables, echo=FALSE, collapse=TRUE
cat("Working directory:")
here()
cat("Meta data file:")
meta_data_file
#
cat("CyTOF clusters data file(s):")
CD45_meta30_clusters_data_file
CD45_merging1_clusters_data_file
CD45_merging2_clusters_data_file
CD45_merging3_clusters_data_file
#
cat("CyTOF clusters summary file(s):")
CD45_meta30_clusters_summary_file
CD45_merging1_clusters_summary_file
CD45_merging2_clusters_summary_file
CD45_merging3_clusters_summary_file
#
cat("CyTOF tSNE data file(s):")
CD45_cells_tSNE_file
#
cat("\nPrefix for output files:")
out_file_prefix


#' ***
#' ### 1.1 Read in  and inspect HTP metadata  
#'   
#+ read_meta, warning = FALSE, message = FALSE, collapse = FALSE
# 1.1 Read in and inspect HTP metadata  ------
#
meta_data <- meta_data_file %>% 
  read_tsv() %>% 
  filter(Event_name != "Current") %>%  # remove unneeded lines
  rename(Age = Age_at_visit) %>% 
  # select(RecordID, LabID, Event_name, Sex, Karyotype, Age = Age_at_visit) %>% # select relevant variables
  mutate( # Set factor orders
    Event_name = fct_relevel(Event_name, c("Visit 1", "Visit 2", "Visit 3")),
    Sex = fct_relevel(Sex, "Female"),
    Karyotype = fct_relevel(Karyotype, "Control")
  )
#
# cat("HTP metadata summary:")
# meta_data %>% 
#   skimr::skim()
# # cat("HTP comorbidity summary:")
# # htp_comorb %>% 
# #   skimr::skim_to_wide()
# # META DATA SUMMARY
# cat("Expt/Cohort summary information:")
# # meta_data %>% distinct(Sampleid)
# # meta_data %>% distinct(RecordID)
# # meta_data %>% count(Source)
# # meta_data %>% count(Event_name)
# # meta_data %>% count(Sex)
# # meta_data %>% count(Karyotype)
# # meta_data %>% count(Sex, Karyotype)
# meta_data %>% count(Karyotype, Source)
# meta_data %>% 
#   group_by(Sex, Karyotype) %>% 
#   summarize(
#     n = n(),
#     median_Age = median(Age),
#     mean_Age = mean(Age),
#     sd = sd(Age),
#     range = paste0(min(Age), " - ", max(Age))
#   )
# #


#' ***
#' ### 1.2 Read in P4C CyTOF CD45+ clusters frequency data 
#'   
#+ read_freq_data, warning=FALSE, message=FALSE, collapse=FALSE
#
# 1.2 Read in P4C CyTOF CD45+ clusters frequency data ----------
#
# Read in CD45+ clusters frequencies: meta30 ----
CD45_meta30_clusters_data <- CD45_meta30_clusters_data_file %>%
  read_tsv() %>% 
  mutate(meta30_cluster = fct_inorder(as.character(meta30_cluster))) %>% 
  inner_join(meta_data) # all are matched
CD45_meta30_clusters_data
CD45_meta30_clusters_data %>% distinct(meta30_cluster) # 30 unique
CD45_meta30_clusters_data %>% distinct(RecordID) # 388
CD45_meta30_clusters_data %>% distinct(LabID) # 388
#
meta30_summary <- CD45_meta30_clusters_summary_file %>% 
  read_tsv()
meta30_summary
#
meta30_clusters_naming_merging <- meta30_clusters_naming_merging_file %>% 
  read_excel() %>% 
  select(meta30 = meta30_cluster, expert_naming, merging1, merging2) %>% 
  mutate(meta30 = fct_inorder(as.character(meta30))) %>% 
  unite(meta30_name, meta30, expert_naming, sep = ".", remove = FALSE)
meta30_clusters_naming_merging
#
# Sina plot Meta30 cluster proportions -----
CD45_meta30_clusters_data %>% 
  mutate(
    # arrange clusters by overall percentage
    meta30_cluster = fct_relevel(meta30_cluster, meta30_summary %>% arrange(perc_meta30) %>% pull(meta30_cluster) %>% as.character())
  ) %>% 
  ggplot(aes(meta30_cluster, gtools::logit(perc_meta30/100), color = Karyotype)) +
  geom_sina() +
  geom_boxplot(aes(group = paste(meta30_cluster, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  labs(title = "Meta30 % distributions among Live cells by Karyotype", y = "% among live")
CD45_meta30_clusters_data %>% 
  mutate(
    # arrange clusters by overall percentage
    meta30_cluster = fct_relevel(meta30_cluster, meta30_summary %>% arrange(perc_meta30) %>% pull(meta30_cluster) %>% as.character())
  ) %>% 
  ggplot(aes(meta30_cluster, gtools::logit(perc_meta30/100), color = Karyotype)) +
  geom_sina() +
  geom_boxplot(aes(group = paste(meta30_cluster, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  labs(title = "Meta30 % distributions among Live cells by Karyotype", y = "% among live") +
  facet_wrap(~ meta30_cluster, scales = "free") + theme(aspect.ratio = 1.3)
#
# Read in CD45+ clusters frequencies: merging1 ----
CD45_merging1_clusters_data <- CD45_merging1_clusters_data_file %>%
  read_tsv() %>% 
  mutate(merging1 = fct_inorder(as.character(merging1))) %>% 
  inner_join(meta_data) # all are matched
CD45_merging1_clusters_data
CD45_merging1_clusters_data %>% distinct(merging1) # 8 unique
CD45_merging1_clusters_data %>% distinct(RecordID) # 388
CD45_merging1_clusters_data %>% distinct(LabID) # 388
#
merging1_summary <- CD45_merging1_clusters_summary_file %>% 
  read_tsv()
merging1_summary
# Sina plot merging1 cluster proportions -----
CD45_merging1_clusters_data %>% 
  mutate(
    # arrange clusters by overall percentage
    merging1 = fct_relevel(merging1, merging1_summary %>% arrange(perc_merging1) %>% pull(merging1) %>% as.character())
  ) %>% 
  ggplot(aes(merging1, gtools::logit(perc_merging1/100), color = Karyotype)) +
  geom_sina() +
  geom_boxplot(aes(group = paste(merging1, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  labs(title = "Merging1 % distributions among Live cells by Karyotype", y = "% among live")
# CD45_meta30_clusters_data %>% 
#   mutate(
#     # arrange clusters by overall percentage
#     meta30_cluster = fct_relevel(meta30_cluster, meta30_summary %>% arrange(perc_meta30) %>% pull(meta30_cluster) %>% as.character())
#   ) %>% 
#   ggplot(aes(meta30_cluster, gtools::logit(perc_meta30/100), color = Karyotype)) +
#   geom_sina() +
#   geom_boxplot(aes(group = paste(meta30_cluster, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
#   scale_color_manual(values=standard_colors) +
#   scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
#   labs(title = "Meta30 % distributions among Live cells by Karyotype", y = "% among live") +
#   facet_wrap(~ meta30_cluster, scales = "free") + theme(aspect.ratio = 1.3)
#
# Read in CD45+ clusters frequencies: merging2 ----
CD45_merging2_clusters_data <- CD45_merging2_clusters_data_file %>%
  read_tsv() %>% 
  mutate(merging2 = fct_inorder(as.character(merging2))) %>% 
  inner_join(meta_data) # all are matched
CD45_merging2_clusters_data
CD45_merging2_clusters_data %>% distinct(merging2) # 18 unique
CD45_merging2_clusters_data %>% distinct(RecordID) # 388
CD45_merging2_clusters_data %>% distinct(LabID) # 388
#
merging2_summary <- CD45_merging2_clusters_summary_file %>% 
  read_tsv()
merging2_summary
# Sina plot merging2 cluster proportions -----
CD45_merging2_clusters_data %>% 
  mutate(
    # arrange clusters by overall percentage
    merging2 = fct_relevel(merging2, merging2_summary %>% arrange(perc_merging2) %>% pull(merging2) %>% as.character())
  ) %>% 
  ggplot(aes(merging2, gtools::logit(perc_merging2/100), color = Karyotype)) +
  geom_sina() +
  geom_boxplot(aes(group = paste(merging2, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  labs(title = "Merging2 % distributions among Live cells by Karyotype", y = "% among live")
# CD45_meta30_clusters_data %>% 
#   mutate(
#     # arrange clusters by overall percentage
#     meta30_cluster = fct_relevel(meta30_cluster, meta30_summary %>% arrange(perc_meta30) %>% pull(meta30_cluster) %>% as.character())
#   ) %>% 
#   ggplot(aes(meta30_cluster, gtools::logit(perc_meta30/100), color = Karyotype)) +
#   geom_sina() +
#   geom_boxplot(aes(group = paste(meta30_cluster, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
#   scale_color_manual(values=standard_colors) +
#   scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
#   labs(title = "Meta30 % distributions among Live cells by Karyotype", y = "% among live") +
#   facet_wrap(~ meta30_cluster, scales = "free") + theme(aspect.ratio = 1.3)
#
# Read in CD45+ clusters frequencies: merging3 ----
CD45_merging3_clusters_data <- CD45_merging3_clusters_data_file %>%
  read_tsv() %>% 
  mutate(merging3 = fct_inorder(as.character(merging3))) %>% 
  inner_join(meta_data) # all are matched
CD45_merging3_clusters_data
CD45_merging3_clusters_data %>% distinct(merging3) # 24 unique
CD45_merging3_clusters_data %>% distinct(RecordID) # 388
CD45_merging3_clusters_data %>% distinct(LabID) # 388
#
merging3_summary <- CD45_merging3_clusters_summary_file %>% 
  read_tsv() %>% 
  arrange(num) %>% 
  mutate(
    merging3 = fct_inorder(merging3)
  )
merging3_summary
# Sina plot merging3 cluster proportions -----
CD45_merging3_clusters_data %>% 
  mutate(
    # arrange clusters by overall percentage
    merging3 = fct_relevel(merging3, merging3_summary %>% arrange(merging3) %>% pull(merging3) %>% as.character())
  ) %>% 
  ggplot(aes(merging3, gtools::logit(perc_merging3/100), color = Karyotype)) +
  geom_sina() +
  geom_boxplot(aes(group = paste(merging3, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  labs(title = "Merging3 % distributions among Live cells by Karyotype", y = "% among live")
# CD45_meta30_clusters_data %>% 
#   mutate(
#     # arrange clusters by overall percentage
#     meta30_cluster = fct_relevel(meta30_cluster, meta30_summary %>% arrange(perc_meta30) %>% pull(meta30_cluster) %>% as.character())
#   ) %>% 
#   ggplot(aes(meta30_cluster, gtools::logit(perc_meta30/100), color = Karyotype)) +
#   geom_sina() +
#   geom_boxplot(aes(group = paste(meta30_cluster, Karyotype)), notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75, position = position_dodge(width = 0.9)) +
#   scale_color_manual(values=standard_colors) +
#   scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
#   labs(title = "Meta30 % distributions among Live cells by Karyotype", y = "% among live") +
#   facet_wrap(~ meta30_cluster, scales = "free") + theme(aspect.ratio = 1.3)
#
# Make combined frequency table ----
cytof_data_combined <- bind_rows(
  CD45_meta30_clusters_data %>% 
    select(LabID, cluster = meta30_cluster, count = n_meta30, percent = perc_meta30) %>% 
    mutate(level = "meta30"),
  CD45_merging1_clusters_data %>% 
    select(LabID, cluster = merging1, count = n_merging1, percent = perc_merging1) %>% 
    mutate(level = "merging1"),
  CD45_merging2_clusters_data %>% 
    select(LabID, cluster = merging2, count = n_merging2, percent = perc_merging2) %>% 
    mutate(level = "merging2"),
  CD45_merging3_clusters_data %>% 
    select(LabID, cluster = merging3, count = n_merging3, percent = perc_merging3) %>% 
    mutate(level = "merging3")
) %>%
  inner_join(meta_data) %>% 
  mutate(
    Sample_type = "PBMCs",
    Data_type = "CyTOF CD45+ FlowSOM clusters",
  ) %>%
  group_by(level, cluster, Karyotype) %>%
  mutate(Extreme_outlier = rstatix::is_extreme(gtools::logit(percent / 100))) %>%
  ungroup() %>% 
  # select(
  #   RecordID,
  #   LabID,
  #   ExperimentID,
  #   Sample_type,
  #   Data_type,
  #   uniquePopulationName,
  #   percentOf,
  #   Units = units,
  #   Value = value,
  #   Extreme_outlier,
  #   Date_exported,
  #   Data_contact,
  #   Script = script
  # ) %>%
  unite(level, cluster, col = "uniqueID", sep = "|", remove = FALSE)
#
# Export combined frequency table for Explorer ----
# bind_rows(
#   live_cells_data, # NEEDS SOME ADJUSTING TO WORK
#   Bcells_data,
#   Tcells_data
# ) %>% 
#   mutate(
#     Sample_type = "PBMCs",
#     Data_type = "CyTOF",
#   ) %>% 
#   group_by(uniquePopulationName, Karyotype) %>% 
#   mutate(Extreme_outlier = rstatix::is_extreme(gtools::logit(value / 100))) %>% 
#   ungroup() %>% 
#   select(
#     RecordID,
#     LabID,
#     ExperimentID,
#     Sample_type,
#     Data_type,
#     uniquePopulationName,
#     percentOf,
#     Units = units,
#     Value = value,
#     Extreme_outlier,
#     Date_exported,
#     Data_contact,
#     Script = script
#   ) %>% 
#   write_tsv(path = here("data", "P4C_CyTOF_Frequencies_longFormat_TrisomExplorer_080521.txt.gz"))
# #


#' ***
#' ### 2. Beta regression model for proportions: Karyotype
#+ beta_karyotype, warning=FALSE, message=FALSE, collapse=FALSE,include=TRUE
#
# 2 Beta regression model for proportions: Karyotype  --------------------------------------------------------
#
# 2.1 check distributions / 0%s: ------
cytof_data_combined %>% 
  mutate(proportion = percent / 100) %>% 
  ggplot(aes(proportion)) + geom_density() + 
  facet_wrap(~uniqueID, scales = "free")
cytof_data_combined %>% 
  mutate(proportion = percent / 100) %>% 
  ggplot(aes(gtools::logit(proportion))) + geom_density() + 
  facet_wrap(~uniqueID, scales = "free")
# check 0s:
cat("Clusters with most zeroes:")
cytof_data_combined %>% 
  group_by(uniqueID) %>% 
  summarize(
    n_zeros = sum(percent == 0)
  ) %>% 
  arrange(-n_zeros)
cat("Clusters with no zeroes:")
cytof_data_combined %>% 
  group_by(uniqueID) %>% 
  mutate(
    n_zeros = sum(percent == 0)
  ) %>% 
  filter(n_zeros == 0) %>% 
  distinct(uniqueID)
#
# 2.2 check for outliers  --------------------------------------------------------
cytof_data_combined %>% filter(Extreme_outlier == TRUE) %>% count(uniqueID) %>% arrange(-n)
# cat("Clusters with extreme outliers (zeroes replaced with NA):")
# cytof_data_combined %>% 
#   mutate(proportion = percent / 100) %>% 
#   filter(str_detect(uniqueID, "CD4\\+CD8\\+\\|CD3\\+|CD4\\+CD8\\+\\|Live|MDSC\\|Live")) %>% # can be list
#   group_by(uniqueID, Karyotype) %>% 
#   mutate(
#     outlier = rstatix::is_outlier(gtools::logit(proportion)),
#     extreme = rstatix::is_extreme(gtools::logit(proportion))
#   ) %>% 
#   ungroup() %>%
#   ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) +
#   geom_sina(data = . %>% filter(extreme == FALSE & outlier == FALSE)) +
#   geom_boxplot(data = . %>% filter(extreme == FALSE), notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
#   geom_sina(data = . %>% filter(extreme == FALSE & outlier == TRUE), color = "orange") +
#   # geom_sina(data = . %>% filter(extreme == TRUE), color = "red") +
#   geom_point(data = . %>% filter(extreme == TRUE), color = "red") + # geom_sina fails with too few points
#   scale_color_manual(values = standard_colors) +
#   facet_wrap(~ uniqueID, scales = "free", nrow = 1) +
#   labs(
#     title = "Subpopulations with highest number of extreme outliers",
#     subtitle = "Unadjusted; Red = Extreme outliers (beyond Q1/Q3 -/+ 3 x IQR",
#     x = NULL
#   ) +
#   theme(
#     aspect.ratio = 1.3,
#     axis.text.x = element_blank(),
#     legend.position = "bottom"
#   )
# ggsave(filename = here("plots", paste0(out_file_prefix, "sina_outliers", ".png")), width = 10, height = 5, units = "in")
# ggsave(filename = here("plots", paste0(out_file_prefix, "sina_outliers", ".pdf")), device = cairo_pdf, width = 10, height = 5, units = "in")
# #
#
# 2.3 Q-Q Plot(s) --------------------------------------------------------
cytof_data_combined %>% 
  filter(level == "meta30") %>% 
  mutate(proportion = percent / 100) %>% 
  group_by(uniqueID, Karyotype) %>% 
  mutate(
    logit_proportion = gtools::logit(proportion),
    n_zeros = sum(proportion == 0), # no 0s allowed; may want to replace with NA and filter by minimum N
    extreme_proportion = rstatix::is_extreme(logit_proportion)
  ) %>% 
  ungroup() %>%
  # filter(n_zeros == 0) %>%  # no 0s allowed; may want to replace with NA and filter by minimum N
  ggpubr::ggqqplot("logit_proportion", facet.by = "uniqueID", title = "Q-Q: meta30", scales = "free")
#
cytof_data_combined %>% 
  filter(level == "merging1") %>% 
  mutate(proportion = percent / 100) %>% 
  group_by(uniqueID, Karyotype) %>% 
  mutate(
    logit_proportion = gtools::logit(proportion),
    n_zeros = sum(proportion == 0), # no 0s allowed; may want to replace with NA and filter by minimum N
    extreme_proportion = rstatix::is_extreme(logit_proportion)
  ) %>% 
  ungroup() %>%
  # filter(n_zeros == 0) %>%  # no 0s allowed; may want to replace with NA and filter by minimum N
  ggpubr::ggqqplot("logit_proportion", facet.by = "uniqueID", title = "Q-Q: merging1", scales = "free")
#
cytof_data_combined %>% 
  filter(level == "merging2") %>% 
  mutate(proportion = percent / 100) %>% 
  group_by(uniqueID, Karyotype) %>% 
  mutate(
    logit_proportion = gtools::logit(proportion),
    n_zeros = sum(proportion == 0), # no 0s allowed; may want to replace with NA and filter by minimum N
    extreme_proportion = rstatix::is_extreme(logit_proportion)
  ) %>% 
  ungroup() %>%
  # filter(n_zeros == 0) %>%  # no 0s allowed; may want to replace with NA and filter by minimum N
  ggpubr::ggqqplot("logit_proportion", facet.by = "uniqueID", title = "Q-Q: merging2", scales = "free")
#
# 2.4 Set up models for assessment  -----------------------------------------------------------
betareg_data <- cytof_data_combined %>% 
  mutate(proportion = percent / 100) %>% 
  # 0 replacement with NA
  mutate(proportion = if_else(proportion == 0, NA_real_, proportion)) %>%
  group_by(uniqueID, Karyotype) %>% 
  mutate(
    logit_proportion = gtools::logit(proportion),
    extreme_proportion = rstatix::is_extreme(logit_proportion)
  ) %>% 
  ungroup() %>% 
  filter(extreme_proportion == FALSE) %>%
  select(RecordID, LabID, uniqueID, level, cluster, percent, proportion, Sex, Age, Karyotype, Sample_source)
# # filter by minimum N per group
betareg_data_final <- betareg_data %>%
  filter(!is.na(proportion)) %>%
  group_by(uniqueID) %>%  # CHECK CORRECT GROUPING
  add_count(Karyotype, Sex) %>% # count by EACH categorical variable
  filter(n >= 15) %>% # require at least NN samples in each category # CURRENTLY KEEPING 56 of 56 clusters
  mutate( # count number of levels for EACH categorical variable
    Karyotype_levels = Karyotype %>% fct_drop() %>% levels() %>% length(),
    Sex_levels = Sex %>% fct_drop() %>% levels() %>% length()
    ) %>%
  filter(Karyotype_levels > 1 & Sex_levels > 1) %>% # need to require >1 level for each categorical level or lm() gives error
  ungroup() %>%
  select(RecordID, LabID, uniqueID, level, cluster, percent, proportion, Sex, Age, Karyotype, Sample_source) %>% 
  nest(data = c(RecordID, LabID, percent, proportion, Sex, Age, Karyotype, Sample_source))
#
## Nested models
# 2.4.1 Simple  -----------------------------------------------------------
tic("Simple beta regressions")
betareg_simple <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    # vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc() # ~ 10 sec; no warnings
# 2.4.2 Multi Sex+Age  -----------------------------------------------------------
tic("Multivariable beta regressions with Sex")
betareg_multi_Karyotype_Sex <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc() # ~ 6 sec; no warnings
tic("Multivariable beta regressions with Age")
betareg_multi_Karyotype_Age <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Age, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc() # ~ 6 sec; no warnings
tic("Multivariable beta regressions with Sex + Age")
betareg_multi_Karyotype_SexAge <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc() # ~ 6 sec; no warnings
# 2.4.3 Multi Sex+Age+Source  -----------------------------------------------------------
tic("Multivariable beta regressions with Sex + Age + Source")
betareg_multi_Karyotype_SexAgeSource <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc() # ~ 6 sec; no warnings
# 2.4.4 Variable dispersions -----------------------------------------------------------
tic("Multivariable beta regressions with Sex + Age + Source | variable Sex dispersion")
betareg_multi_Karyotype_SexAgeSource_varSex <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source | Sex, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc()
tic("Multivariable beta regressions with Sex + Age + Source | variable Age dispersion")
betareg_multi_Karyotype_SexAgeSource_varAge <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source | Age, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc()
tic("Multivariable beta regressions with Sex + Age + Source | variable Source dispersion")
betareg_multi_Karyotype_SexAgeSource_varSource <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source | Sample_source, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc()
tic("Multivariable beta regressions with Sex + Age + Source | variable Sex + Age + Source dispersion")
betareg_multi_Karyotype_SexAgeSource_varSexAgeSource <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source | Sex + Age + Sample_source, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc()
# 2.4.3 Multi Sex+Age alternative link function(s) -----------------------------------------------------------
tic("Multivariable beta regressions with Sex + Age +Source with loglog link")
betareg_multi_Karyotype_SexAgeSource_loglog <- betareg_data_final %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source, data = .x, link = "loglog")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment), # see ?augment.betareg
    vifs = map(fit, ~car::vif(mod = .x) %>% as_tibble(rownames = "term"))
  )
toc()
#
# 2.5 Check model assumptions using gvlma (TESTING - NOT FINAL + potentially slow) --------------------------------------------------------
# gvlma does not work with betareg models
#
# 2.6 AIC comparison(s) --------------------------------------------------------
# The Akaike Information Critera (AIC) scoring metric for performance of models. It enables us to quantify the model fit while penalizing superfluous predictors. We can look at the difference in AIC between models to determine if one is a better fit. Support for for suboptimal model by ∆AIC: 0 – 2: substantial | 4 – 7: considerably less | >10: essentially none
simple_glance <- betareg_simple %>% unnest(glanced)
multi_Sex_glance <- betareg_multi_Karyotype_Sex %>% unnest(glanced)
multi_Age_glance <- betareg_multi_Karyotype_Age %>% unnest(glanced)
multi_SexAge_glance <- betareg_multi_Karyotype_SexAge %>% unnest(glanced)
multi_SexAgeSource_glance <- betareg_multi_Karyotype_SexAgeSource %>% unnest(glanced)
multi_SexAgeSource_varSex_glance <- betareg_multi_Karyotype_SexAgeSource_varSex %>% unnest(glanced)
multi_SexAgeSource_varAge_glance <- betareg_multi_Karyotype_SexAgeSource_varAge %>% unnest(glanced)
multi_SexAgeSource_varSource_glance <- betareg_multi_Karyotype_SexAgeSource_varSource %>% unnest(glanced)
multi_SexAgeSource_varSexAgeSource_glance <- betareg_multi_Karyotype_SexAgeSource_varSexAgeSource %>% unnest(glanced)
multi_SexAgeSource_loglog_glance <- betareg_multi_Karyotype_SexAgeSource_loglog %>% unnest(glanced)
#
a1 <- simple_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_Sex_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
        ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniquePopulationName), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniquePopulationName), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "simple vs. Sex", x = NULL)
#
a2 <- simple_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_Age_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "simple vs. Age", x = NULL)
#
a3 <- simple_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAge_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "simple vs. Sex+Age", x = NULL)
#
a4 <- simple_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "simple vs. Sex+Age+Source", x = NULL)
#
a5 <- multi_Sex_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAge_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "Sex vs. Sex+Age", x = NULL)
#
a6 <- multi_SexAge_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "Sex+Age vs. Sex+Age+Source", x = NULL)
#
av1 <- multi_SexAgeSource_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_varSex_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "Sex+Age+Source vs. \nSex+Age+Source | varSex", x = NULL)
#
av2 <- multi_SexAgeSource_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_varAge_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "Sex+Age+Source vs. \nSex+Age+Source | varAge", x = NULL)
#
av3 <- multi_SexAgeSource_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_varSource_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "Sex+Age+Source vs. \nSex+Age+Source | varSource", x = NULL)
#
av4 <- multi_SexAgeSource_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_varSexAgeSource_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "Sex+Age+Source vs. \nSex+Age+Source | varSexAgeSource", x = NULL)
#
a7 <- multi_SexAgeSource_glance %>% select(uniqueID, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_loglog_glance %>% select(uniqueID, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(uniqueID = fct_inorder(uniqueID)) %>% 
  ggplot(aes(uniqueID, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme() + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    # axis.text.x = element_blank(),  # turn off x-axis ticks with too many
    # axis.ticks.x = element_blank() # turn off x-axis labels with too many
  ) + 
  # geom_text_repel(data = . %>%  slice_max(AIC_diff, n = 5), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = uniqueID), xlim = c(NA, NA), size = 3) + # labels if too many x
  # geom_text_repel(aes(label = uniquePopulationName), xlim = c(NA, NA), size = 3) + # labels if too many x
  labs(title = "Sex+Age+Source vs. \nSex+Age+Source | loglog", x = NULL)
#
library(patchwork)
(a1 + a2) / (a3 + a4) + plot_annotation(title = "Beta regression: AIC comparisons") # plot_spacer()
(a5 + a6) / (plot_spacer() + plot_spacer()) + plot_annotation(title = "Beta regression: AIC comparisons") # plot_spacer()
(av1 + av2) / (av3 + av4) + plot_annotation(title = "Beta regression w/var: AIC comparisons") # plot_spacer()
a7
cat("Preferred model(s) based on AIC comparisons:
    proportion ~ Karyotype + Sex + Age + Source
    proportion ~ Karyotype + Sex + Age + Source | varSource?
    ") # 
#
# 2.7 Check max VIFs -----
betareg_multi_Karyotype_Sex %>% unnest(vifs) %>% select(uniqueID, term, VIF = value) %>% group_by(term) %>% slice_max(VIF, n = 1)
betareg_multi_Karyotype_Age %>% unnest(vifs) %>% select(uniqueID, term, VIF = value) %>% group_by(term) %>% slice_max(VIF, n = 1)
betareg_multi_Karyotype_SexAge %>% unnest(vifs) %>% select(uniqueID, term, VIF = value) %>% group_by(term) %>% slice_max(VIF, n = 1)
betareg_multi_Karyotype_SexAgeSource %>% unnest(vifs) %>% select(uniqueID, term, GVIF) %>% group_by(term) %>% slice_max(GVIF, n = 1)


#' ***
#' ### 3. Model results
#+ model_results, warning=FALSE, message=FALSE, collapse=FALSE,include=TRUE
#
# 3.1 Extract and plot results for simple model(s) ------------------------------------------------------------
betareg_simple_results <- betareg_simple %>% 
  unnest(tidied) %>% 
  select(uniqueID, level, cluster, component, term, estimate, p.value) %>% 
  group_by(uniqueID) %>% 
  summarize(
    uniqueID = first(uniqueID),
    level = first(level),
    cluster = first(cluster),
    logit_Control = nth(estimate, n = 1),
    logit_T21 = nth(estimate, n = 2) + logit_Control,
    pval = nth(p.value, n = 2),
    Control = gtools::inv.logit(logit_Control), # = Checked
    T21 = gtools::inv.logit(logit_T21), # = Checked
    FoldChange = exp(nth(estimate, n = 2)) # = Checked
  ) %>% 
  group_by(level) %>% 
  arrange(level, pval) %>% 
  filter(!str_detect(cluster, "EXCLUDE")) %>% # remove excluded clusters from merging3
  mutate(
    FoldChange = T21 / Control, # similar to exp(nth(estimate, n = 1))
    BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))
  ) %>% 
  ungroup() %>% 
  select(uniqueID, level, cluster, mean_prop_Control = Control, mean_prop_T21 = T21, FoldChange, pval, BHadj_pval) %>% 
  mutate(Model = "betareg(proportion ~ Karyotype)")
# manual check of transformation and FoldChange extraction for betareg models -----
# # get means for 2 groups and calc FoldChange
# betareg_data_final[1,] %>% 
#   unnest(data) %>% 
#   group_by(Karyotype) %>% 
#   summarize(mean = mean(proportion, na.rm = TRUE)) %>% 
#   pivot_wider(names_from = Karyotype, values_from = mean) %>% 
#   mutate(
#     FoldChange = T21 / Control, # not sure why this is different
#     diff = T21 - Control,
#     gtools::logit(T21) - gtools::logit(Control),
#     # exp(gtools::logit(T21) - gtools::logit(Control)),
#     # gtools::inv.logit(gtools::logit(T21) - gtools::logit(Control))
#     )
# # exponentiate model coefs and calc FoldChange
# betareg_simple[1,] %>% 
#   unnest(tidied) %>% 
#   select(uniquePopulationName, component, term, estimate) %>% 
#   filter(component != "precision") %>% 
#   pivot_wider(names_from = term, values_from = estimate) %>% 
#   mutate(
#     FoldChange = exp((`(Intercept)` + KaryotypeT21)) / exp(`(Intercept)`), # 0.5496285 = same as exp below
#     m1 = exp(`(Intercept)`)/(1+exp(`(Intercept)`)),
#     m2 = exp(`(Intercept)` + KaryotypeT21)/(1+exp(`(Intercept)` + KaryotypeT21)),
#     ratio = m2 / m1 # 0.5591801 = same as inv.logit on coefs
#     )
# # inv.logit model coefs and calc FoldChange
# betareg_simple[1,] %>% 
#   unnest(tidied) %>% 
#   select(uniquePopulationName, component, term, estimate) %>% 
#   filter(component != "precision") %>% 
#   pivot_wider(names_from = term, values_from = estimate) %>% 
#   mutate(
#     m1 = gtools::inv.logit(`(Intercept)`),
#     m2 = gtools::inv.logit((`(Intercept)` + KaryotypeT21)),
#     FoldChange = gtools::inv.logit((`(Intercept)` + KaryotypeT21)) / gtools::inv.logit(`(Intercept)`) # 0.5591801
#     )
# # by direct transformation of estimates
# betareg_simple[1,] %>% 
#   unnest(tidied) %>% 
#   select(uniquePopulationName, component, term, estimate) %>% 
#   mutate(
#     exp = exp(estimate), # 0.5496285 = ALMOST SAME AS ABOVE!
#     invlogit = gtools::inv.logit(estimate) # DEFINITELY WRONG
#     )
# # # from https://stats.stackexchange.com/questions/120472/how-to-interpret-coefficients-of-a-beta-regression-model-with-logit-link
# # betareg_simple[1,] %>% 
# #   unnest(tidied) %>% 
# #   select(uniquePopulationName, component, term, estimate) %>% 
# #   filter(component != "precision") %>% 
# #   pivot_wider(names_from = term, values_from = estimate) %>% 
# #   mutate(
# #     Diff = exp(`(Intercept)` + KaryotypeT21) / (1 + exp(`(Intercept)` + KaryotypeT21))
# #     )
# # See also https://stats.stackexchange.com/questions/297659/interpretation-of-betareg-coef
#
# Export long to text
betareg_simple_results %>% 
  left_join(meta30_clusters_naming_merging, by = c("cluster" = "meta30")) %>% # add names to meta30 clusters
  write_tsv(file = here("results", paste0(out_file_prefix, "betareg_results_simple_Karyotype", ".txt")))
# Export to Excel ---------
betareg_simple_results %>% 
  left_join(meta30_clusters_naming_merging, by = c("cluster" = "meta30")) %>% # add names to meta30 clusters
  split(.$level) %>% 
  export_excel(filename = "betareg_results_simple_Karyotype")
# Volcano simple ----
meta30_clusters_naming_merging %>% 
  rename(cluster = meta30) %>% 
  inner_join(
    betareg_simple_results %>% 
  filter(level == "meta30")
  )
meta30_clusters_naming_merging %>% 
  rename(cluster = meta30) %>% 
  inner_join(
    betareg_simple_results %>% 
      filter(level == "meta30")
  ) %>% 
  mutate(cluster = meta30_name) %>% 
  volcano_plot_lab(
    title="Meta30: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_meta30", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_meta30", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_simple_results %>% 
  filter(level == "merging1")
betareg_simple_results %>% 
  filter(level == "merging1") %>% 
  volcano_plot_lab(
    title="merging1: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_merging1", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_merging1", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_simple_results %>% 
  filter(level == "merging2")
betareg_simple_results %>% 
  filter(level == "merging2") %>% 
  volcano_plot_lab(
    title="merging2: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_merging2", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_merging2", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_simple_results %>% 
  filter(level == "merging3")
betareg_simple_results %>% 
  filter(level == "merging3") %>% 
  volcano_plot_lab(
    title="merging3: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_merging3", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple_Karyotype_merging3", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
# 3.2 Extract and plot results for multivar model: proportion ~ Karyotype + Sex + Age + Source (PREFERRED MODEL) ------------------------------------------------------------
betareg_multi_Karyotype_SexAgeSource_results <- betareg_multi_Karyotype_SexAgeSource %>% 
  unnest(tidied) %>% 
  select(uniqueID, level, cluster, component, term, estimate, p.value) %>% 
  group_by(uniqueID) %>% 
  summarize(
    uniqueID = first(uniqueID),
    level = first(level),
    cluster = first(cluster),
    logit_Control = nth(estimate, n = 1),
    logit_T21 = nth(estimate, n = 2) + logit_Control,
    pval = nth(p.value, n = 2),
    Control = gtools::inv.logit(logit_Control), # = Checked
    T21 = gtools::inv.logit(logit_T21), # = Checked
    FoldChange = exp(nth(estimate, n = 2)) # = Checked
  ) %>% 
  group_by(level) %>% 
  arrange(level, pval) %>% 
  filter(!str_detect(cluster, "EXCLUDE")) %>% # remove excluded clusters from merging3
  mutate(
    FoldChange = T21 / Control, # similar to exp(nth(estimate, n = 1))
    BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))
  ) %>% 
  ungroup() %>% 
  select(uniqueID, level, cluster, mean_prop_Control = Control, mean_prop_T21 = T21, FoldChange, pval, BHadj_pval) %>% 
  mutate(Model = "betareg(proportion ~ Karyotype + Sex + Age + Sample_source)")
#
# Export long to text
betareg_multi_Karyotype_SexAgeSource_results %>% 
  left_join(meta30_clusters_naming_merging, by = c("cluster" = "meta30")) %>% # add names to meta30 clusters
  write_tsv(file = here("results", paste0(out_file_prefix, "betareg_results_multi_Karyotype_SexAgeSource", ".txt")))
# Export to Excel ---------
betareg_multi_Karyotype_SexAgeSource_results %>% 
  left_join(meta30_clusters_naming_merging, by = c("cluster" = "meta30")) %>% # add names to meta30 clusters
  split(.$level) %>% 
  export_excel(filename = "betareg_results_multi_Karyotype_SexAgeSource")
# Volcano ----
meta30_clusters_naming_merging %>% 
  rename(cluster = meta30) %>% 
  inner_join(
    betareg_multi_Karyotype_SexAgeSource_results %>% 
      filter(level == "meta30")
  )
meta30_clusters_naming_merging %>% 
  rename(cluster = meta30) %>% 
  inner_join(
    betareg_multi_Karyotype_SexAgeSource_results %>% 
      filter(level == "meta30")
  ) %>% 
  mutate(cluster = meta30_name) %>% 
  volcano_plot_lab(
    title="Meta30: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_meta30", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_meta30", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_multi_Karyotype_SexAgeSource_results %>% 
  filter(level == "merging1")
betareg_multi_Karyotype_SexAgeSource_results %>% 
  filter(level == "merging1") %>% 
  volcano_plot_lab(
    title="merging1: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_merging1", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_merging1", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_multi_Karyotype_SexAgeSource_results %>% 
  filter(level == "merging2")
betareg_multi_Karyotype_SexAgeSource_results %>% 
  filter(level == "merging2") %>% 
  volcano_plot_lab(
    title="merging2: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_merging2", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_merging2", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_multi_Karyotype_SexAgeSource_results %>% 
  filter(level == "merging3")
betareg_multi_Karyotype_SexAgeSource_results %>% 
  filter(level == "merging3") %>% 
  volcano_plot_lab(
    title="merging3: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source ","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]", "\ndown                                                   up")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_merging3", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_merging3", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
# 3.3 Extract and plot results for multivar model: proportion ~ Karyotype + Sex + Age + Source | variable source dispersion (PREFERRED MODEL?) ------------------------------------------------------------
betareg_multi_Karyotype_SexAgeSource_varSource_results <- betareg_multi_Karyotype_SexAgeSource_varSource %>% 
  unnest(tidied) %>% 
  select(uniqueID, level, cluster, component, term, estimate, p.value) %>% 
  group_by(uniqueID) %>% 
  summarize(
    uniqueID = first(uniqueID),
    level = first(level),
    cluster = first(cluster),
    logit_Control = nth(estimate, n = 1),
    logit_T21 = nth(estimate, n = 2) + logit_Control,
    pval = nth(p.value, n = 2),
    Control = gtools::inv.logit(logit_Control), # = Checked
    T21 = gtools::inv.logit(logit_T21), # = Checked
    FoldChange = exp(nth(estimate, n = 2)) # = Checked
  ) %>% 
  group_by(level) %>% 
  arrange(level, pval) %>% 
  filter(!str_detect(cluster, "EXCLUDE")) %>% # remove excluded clusters from merging3
  mutate(
    FoldChange = T21 / Control, # similar to exp(nth(estimate, n = 1))
    BHadj_pval = p.adjust(pval, method = "BH", n = length(pval))
  ) %>% 
  ungroup() %>% 
  select(uniqueID, level, cluster, mean_prop_Control = Control, mean_prop_T21 = T21, FoldChange, pval, BHadj_pval) %>% 
  mutate(Model = "betareg(proportion ~ Karyotype + Sex + Age + Sample_source | Sample_source)")
#
# Export long to text
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  left_join(meta30_clusters_naming_merging, by = c("cluster" = "meta30")) %>% # add names to meta30 clusters
  write_tsv(file = here("results", paste0(out_file_prefix, "betareg_multi_Karyotype_SexAgeSource_varSource", ".txt")))
# Export to Excel ---------
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  left_join(meta30_clusters_naming_merging, by = c("cluster" = "meta30")) %>% # add names to meta30 clusters
  split(.$level) %>% 
  export_excel(filename = "betareg_multi_Karyotype_SexAgeSource_varSource")
# Volcano ----
meta30_clusters_naming_merging %>% 
  rename(cluster = meta30) %>% 
  inner_join(
    betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
      filter(level == "meta30")
  )
meta30_clusters_naming_merging %>% 
  rename(cluster = meta30) %>% 
  inner_join(
    betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
      filter(level == "meta30")
  ) %>% 
  mutate(cluster = meta30_name) %>% 
  volcano_plot_lab(
    title="Meta30: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source|Source\n","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_meta30", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_meta30", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging1")
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging1") %>% 
  volcano_plot_lab(
    title="merging1: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source|Source\n","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_merging1", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_merging1", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging2")
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging2") %>% 
  volcano_plot_lab(
    title="merging2: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source|Source\n","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_merging2", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_merging2", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging3")
betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging3") %>% 
  volcano_plot_lab(
    title="merging3: Diff. proportion in T21", 
    subtitle = paste0("betareg:proportion~Karyotype+Sex+Age+Source|Source\n","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]")
  )
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_merging3", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi_Karyotype_SexAgeSource_varSource_merging3", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
#
# Individual plots meta30 clusters ------------------------------------------------
sig_meta30_up <- betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "meta30") %>% 
  filter(BHadj_pval<0.1) %>%
  # arrange(pval) %>% 
  # slice_min(BHadj_pval, n = 4) %>%
  # .[1:5,] %>% 
  arrange(-FoldChange) %>% 
  slice_max(order_by = FoldChange, n = 5)
sig_meta30_dn <- betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "meta30") %>% 
  filter(BHadj_pval<0.1) %>%
  arrange(FoldChange) %>% 
  slice_min(order_by = FoldChange, n = 5)
# 
betareg_data_final %>% 
  unnest(data) %>% 
  filter(uniqueID %in% sig_meta30_up$uniqueID) %>%
  mutate(cluster = fct_relevel(cluster, sig_meta30_up$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, proportion, color = Karyotype)) +
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
    ) +
  labs(title="Meta30: cluster proportions by Karyotype",
       subtitle="Top upreg. by betareg; unadjusted",
       x = "ratio among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopUpreg_sina_unadj", ".png")), width=12, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopUpreg_sina_unadj", ".pdf")), device=cairo_pdf, width=12, height=3, units="in")
#
# unajusted sina plots ----
# with logit transformation/scale but % labels
betareg_data_final %>% 
  unnest(data) %>% 
  filter(uniqueID %in% sig_meta30_up$uniqueID) %>%
  mutate(cluster = fct_relevel(cluster, sig_meta30_up$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) + 
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
  ) +
  labs(title="Meta30: cluster % by Karyotype",
       subtitle="Top upreg. by betareg; unadjusted; logit transform with % labels",
       x = NULL,
       y = "% among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopUpreg_sina_unadj_logit", ".png")), width=12, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopUpreg_sina_unadj_logit", ".pdf")), device=cairo_pdf, width=12, height=3, units="in")
#
betareg_data_final %>% 
  unnest(data) %>% 
  filter(uniqueID %in% sig_meta30_dn$uniqueID) %>%
  mutate(cluster = fct_relevel(cluster, sig_meta30_dn$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) + 
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
  ) +
  labs(title="Meta30: cluster % by Karyotype",
       subtitle="Top downreg. by betareg; unadjusted; logit transform with % labels",
       x = NULL,
       y = "% among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopDnreg_sina_unadj_logit", ".png")), width=12, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopDnreg_sina_unadj_logit", ".pdf")), device=cairo_pdf, width=12, height=3, units="in")
#
# adjusted sina plots ------
betareg_data_final %>% 
  filter(level == "meta30") %>% 
  unnest(data) %>%
  mutate(logit_proportion = gtools::logit(proportion)) %>%
  effectsize::adjust(effect = "Sex", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Age", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Sample_source", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  as_tibble() %>%
  mutate(proportion = gtools::inv.logit(logit_proportion)) %>% 
  filter(uniqueID %in% sig_meta30_up$uniqueID) %>%
  mutate(cluster = fct_relevel(cluster, sig_meta30_up$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) + 
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
  ) +
  labs(title="Meta30: cluster % by Karyotype",
       subtitle="Top upreg. by betareg; Age/Sex/Source-adjusted; logit transform with % labels",
       x = NULL,
       y = "% among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopUpreg_sina_SexAgeSourceadj_logit", ".png")), width=12, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopUpreg_sina_SexAgeSourceadj_logit", ".pdf")), device=cairo_pdf, width=12, height=3, units="in")
#
betareg_data_final %>% 
  filter(level == "meta30") %>% 
  unnest(data) %>%
  mutate(logit_proportion = gtools::logit(proportion)) %>%
  effectsize::adjust(effect = "Sex", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Age", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Sample_source", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  as_tibble() %>%
  mutate(proportion = gtools::inv.logit(logit_proportion)) %>% 
  filter(uniqueID %in% sig_meta30_dn$uniqueID) %>%
  mutate(cluster = fct_relevel(cluster, sig_meta30_dn$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) + 
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
  ) +
  labs(title="Meta30: cluster % by Karyotype",
       subtitle="Top downreg. by betareg; Age/Sex/Source-adjusted; logit transform with % labels",
       x = NULL,
       y = "% among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopDnreg_sina_SexAgeSourceadj_logit", ".png")), width=12, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "Meta30_TopDnreg_sina_SexAgeSourceadj_logit", ".pdf")), device=cairo_pdf, width=12, height=3, units="in")
#
# Individual plots merging3 clusters ------------------------------------------------
sig_merging3_up <- betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging3") %>% 
  filter(BHadj_pval<0.1) %>%
  # arrange(pval) %>% 
  # slice_min(BHadj_pval, n = 4) %>%
  # .[1:5,] %>% 
  arrange(-FoldChange) %>% 
  slice_max(order_by = FoldChange, n = 5)
sig_merging3_dn <- betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging3") %>% 
  filter(BHadj_pval<0.1) %>%
  arrange(FoldChange) %>% 
  slice_min(order_by = FoldChange, n = 5)
# 
# adjusted sina plots merging3 ------
betareg_data_final %>% 
  filter(level == "merging3") %>% 
  unnest(data) %>%
  mutate(logit_proportion = gtools::logit(proportion)) %>%
  effectsize::adjust(effect = "Sex", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Age", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Sample_source", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  as_tibble() %>%
  mutate(proportion = gtools::inv.logit(logit_proportion)) %>% 
  filter(uniqueID %in% sig_merging3_up$uniqueID) %>%
  mutate(cluster = fct_relevel(cluster, sig_merging3_up$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) + 
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
  ) +
  labs(title="merging3: cluster % by Karyotype",
       subtitle="Top upreg. by betareg; Age/Sex/Source-adjusted; logit transform with % labels",
       x = NULL,
       y = "% among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "merging3_TopUpreg_sina_SexAgeSourceadj_logit", ".png")), width=12, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "merging3_TopUpreg_sina_SexAgeSourceadj_logit", ".pdf")), device=cairo_pdf, width=12, height=3, units="in")
#
betareg_data_final %>% 
  filter(level == "merging3") %>% 
  unnest(data) %>%
  mutate(logit_proportion = gtools::logit(proportion)) %>%
  effectsize::adjust(effect = "Sex", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Age", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Sample_source", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  as_tibble() %>%
  mutate(proportion = gtools::inv.logit(logit_proportion)) %>% 
  filter(uniqueID %in% sig_merging3_dn$uniqueID) %>%
  mutate(cluster = fct_relevel(cluster, sig_merging3_dn$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) + 
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
  ) +
  labs(title="merging3: cluster % by Karyotype",
       subtitle="Top downreg. by betareg; Age/Sex/Source-adjusted; logit transform with % labels",
       x = NULL,
       y = "% among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "merging3_TopDnreg_sina_SexAgeSourceadj_logit", ".png")), width=12, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "merging3_TopDnreg_sina_SexAgeSourceadj_logit", ".pdf")), device=cairo_pdf, width=12, height=3, units="in")
#
# adjusted sina plots for Fig. 5 ------
betareg_data_final %>% 
  filter(level == "merging3") %>% 
  unnest(data) %>%
  mutate(logit_proportion = gtools::logit(proportion)) %>%
  effectsize::adjust(effect = "Sex", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Age", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  effectsize::adjust(effect = "Sample_source", select = c("logit_proportion"), keep_intercept = TRUE) %>%
  as_tibble() %>%
  mutate(proportion = gtools::inv.logit(logit_proportion)) %>% 
  # Signif vs both Karyotype AND IFN score:
  filter(uniqueID %in% c("merging3|naïve CD8+ T", "merging3|CD27+ B", "merging3|naïve CD4+ T",  "merging3|non-classical monocytes",  "merging3|CD4+ TCM", "merging3|CD8+ TEM")) %>%
  # mutate(cluster = fct_relevel(cluster, sig_meta30_up$cluster %>% as.character())) %>% # control plotting order
  arrange(cluster) %>% 
  ggplot(aes(Karyotype, gtools::logit(proportion), color = Karyotype)) + 
  geom_sina() +
  geom_boxplot(notch=TRUE, varwidth=FALSE, outlier.shape=NA, coef=FALSE, width=0.3, color="black", fill="transparent", size=0.75) +
  scale_color_manual(values=standard_colors) +
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  facet_wrap(~ cluster, scales = "free", nrow=1) +
  theme(
    axis.text.x = element_blank(),
    aspect.ratio = 1.3,
    legend.position="bottom"
  ) +
  labs(title="merging3: cluster % by Karyotype",
       subtitle="Signif vs both Karyotype AND IFN score; Age/Sex/Source-adjusted; logit transform with % labels",
       x = NULL,
       y = "% among CD45+"
  )
ggsave(filename=here("plots", paste0(out_file_prefix, "merging3_Fig5_sina_SexAgeSourceadj_logit", ".png")), width=20, height=3, units="in")
ggsave(filename=here("plots", paste0(out_file_prefix, "merging3_Fig5_sina_SexAgeSourceadj_logit", ".pdf")), device=cairo_pdf, width=20, height=3, units="in")
#



#' ***
#' ### 4. tSNE plots
#+ tSNE_plots, warning=FALSE, message=FALSE, collapse=FALSE,include=TRUE
#
# 4.1 read in data for tSNE plots --------------------------------------------------------
tSNE_data_px440_500_cells <- CD45_cells_tSNE_file %>% 
  read_tsv() %>% 
  rename(meta30 = meta30_cluster)
#
meta30_cluster_labels <- tSNE_data_px440_500_cells %>% 
  mutate(meta30 = fct_inorder(as.character(meta30))) %>% 
  group_by(meta30) %>% 
  summarize(
    TSNE1 = median(TSNE1),
    TSNE2 = median(TSNE2)
  ) %>% 
  inner_join(meta30_clusters_naming_merging)
#
meta30_colors <- c("#9f4c78","#62ba3d","#b05dd3","#51c367","#d4469b","#aeb933","#626edd","#79a43f","#844ea3","#388336","#da3e5c","#58cca2","#d04f30","#47c4ce","#dd862f","#6071b6","#cda744","#5f9fd8","#5d781e","#cc8dd1","#68b274","#e47d93","#31754b","#a84949","#37967d","#e2926d","#626c30","#9b632d","#a8b26d","#897a2e")
#
merging3_labels <- tSNE_data_px440_500_cells %>% 
  group_by(merging3) %>% 
  summarize(TSNE1 = median(TSNE1), TSNE2 = median(TSNE2)) %>% 
  inner_join(merging3_summary) %>% 
  filter(!str_detect(merging3, "EXCLUDE"))
#
merging3_colors <- c("#636ad8","#60b646","#af5dcb","#a9b539","#d14396","#5bbe7d","#d34258","#4bb7a7","#d0522d","#5faadc","#d69b3a","#5778be","#757327","#9b72b8","#447d41","#de87bb","#b2ad6a","#9f4866","#a26431","#de8373", "grey", "grey", "grey", "grey")
names(merging3_colors) <- merging3_summary %>% arrange(merging3) %>% pull(merging3_label)
#
# 4.2 tSNE plots --------------------------------------------------------
#
# tSNE px440 plot by cluster: meta30 ------
t440_meta30 <- tSNE_data_px440_500_cells %>% 
  # filter(Karyotype == "T21") %>% # not for this comparison
  group_by(sample_id) %>%
  slice_sample(n = 100) %>% # Subsample to reduce overplotting
  ungroup() %>%
  mutate(meta30 = as_factor(meta30)) %>% # may want to sort to put less abundant clusters on top
  ggplot(aes(TSNE1, TSNE2, color = meta30)) +
  geom_point(size = 0.05) +
  scale_color_manual(values = meta30_colors) +
  geom_text(data = meta30_cluster_labels, aes(label = meta30), color = "black") +
  guides(colour = guide_legend(override.aes = list(size=2), ncol = 2)) + # note ncol = 2
  theme(
    aspect.ratio = 1,
    legend.key=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank()
  ) +
  labs(
    title = "Meta30 Clusters",
    subtitle = "px = 440; 100 cells per sample"
  )
t440_meta30
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_by_meta30_clusters", ".png")), width = 7, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_by_meta30_clusters", ".pdf")), device = cairo_pdf, width = 7, height = 5, units = "in")
#
# tSNE px440 plot by cluster: merging1 ------
t440_merging1 <- tSNE_data_px440_500_cells %>% 
  # filter(Karyotype == "T21") %>% # not for this comparison
  group_by(sample_id) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  # mutate(meta12_cluster = as_factor(meta12_cluster)) %>% # may want to sort to put less abundant clusters on top
  ggplot(aes(TSNE1, TSNE2, color = merging1)) +
  geom_point(size = 0.05) +
  # geom_blank() + # use to get PDF version but no legend - see below to get legend
  guides(colour = guide_legend(override.aes = list(size=2))) +
  geom_text(data = merging3_labels, aes(label = num), color = "black") + # turn off for PNG
  theme(
    aspect.ratio = 1,
    legend.key=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank()
  ) +
  labs(
    title = "Major Lineages",
    subtitle = "Merging3 cluster labels; px = 440; 100 cells per sample"
  )
t440_merging1
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_by_merging1_clusters", ".png")), width = 7, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_by_merging1_clusters", ".pdf")), device = cairo_pdf, width = 7, height = 5, units = "in")
# get just the legend
cowplot::get_legend(t440_merging1) %>% cowplot::ggdraw()
ggsave(filename = here("plots", paste0(out_file_prefix, "LEGEND_tSNE_px440_by_merging1_clusters", ".pdf")), device = cairo_pdf, width = 7, height = 5, units = "in")
#
# tSNE px440 plot by betareg fold-change: meta30 ------
limits <- betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "meta30") %>% 
  pull(FoldChange) %>% 
  log2() %>% 
  max(abs(.)) %>% 
  plyr::round_any(0.01, ceiling) * c(-1, 1) # ENSURE CENTER ON 0
#
t440_meta30_FC <- tSNE_data_px440_500_cells %>% 
  mutate(meta30 = fct_inorder(as.character(meta30))) %>% 
  # filter(Karyotype == "T21") %>% # not for this comparison
  group_by(sample_id) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  inner_join(
    betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
      filter(level == "meta30") %>% 
      select(meta30 = cluster, FoldChange, BHadj_pval)
  ) %>% 
  ggplot(aes(TSNE1, TSNE2, color = log2(FoldChange))) +
  geom_point(size = 0.05) +
  scale_color_distiller(palette = "RdBu", limit = limits, oob = scales::squish) + # ENSURE CENTER ON 0
  geom_text(
    data = meta30_cluster_labels %>% 
      inner_join(betareg_multi_Karyotype_SexAgeSource_varSource_results %>% filter(level == "meta30") %>% select(meta30 = cluster, FoldChange, BHadj_pval)) %>% 
      filter(BHadj_pval < 0.1), 
    aes(label = meta30), color = "black") +
  theme(
    aspect.ratio = 1,
    legend.key=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank()
  ) +
  labs(
    title = "Meta30 Clusters: Fold-change T21/Control",
    subtitle = "betaregression; px = 440; 100 cells per sample"
  )
t440_meta30_FC
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_meta30_by_betareg", ".png")), width = 7, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_meta30_by_betareg", ".pdf")), device = cairo_pdf, width = 7, height = 5, units = "in")
#
# tSNE px440 plot by betareg fold-change: merging3 ------
limits <- betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
  filter(level == "merging3") %>% 
  pull(FoldChange) %>% 
  log2() %>% 
  max(abs(.)) %>% 
  plyr::round_any(0.01, ceiling) * c(-1, 1) # ENSURE CENTER ON 0
#
t440_merging3_FC <- tSNE_data_px440_500_cells %>% 
  # mutate(meta30 = fct_inorder(as.character(meta30))) %>% 
  # filter(Karyotype == "T21") %>% # not for this comparison
  group_by(sample_id) %>%
  slice_sample(n = 100) %>%
  ungroup() %>%
  inner_join(
    betareg_multi_Karyotype_SexAgeSource_varSource_results %>% 
      filter(level == "merging3") %>% 
      select(merging3 = cluster, FoldChange, BHadj_pval)
  ) %>% 
  ggplot(aes(TSNE1, TSNE2, color = log2(FoldChange))) +
  geom_point(size = 0.05) +
  # geom_blank() + # use to get PDF version but no legend - see below to get legend
  scale_color_distiller(palette = "RdBu", limit = limits, oob = scales::squish) + # ENSURE CENTER ON 0
  geom_text( # turn off for PNG
    data = merging3_labels %>%
      inner_join(betareg_multi_Karyotype_SexAgeSource_varSource_results %>% filter(level == "merging3") %>% select(merging3 = cluster, FoldChange, BHadj_pval)) %>%
      filter(BHadj_pval < 0.1),
    aes(label = num), color = "black") +
  theme(
    aspect.ratio = 1,
    legend.key=element_blank(),
    axis.text.x=element_blank(), axis.text.y=element_blank()
  ) +
  labs(
    title = "merging3 Clusters: Fold-change T21/Control",
    subtitle = "betaregression; px = 440; 100 cells per sample"
  )
t440_merging3_FC
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_merging3_by_betareg", ".png")), width = 7, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "tSNE_px440_merging3_by_betareg", ".pdf")), device = cairo_pdf, width = 7, height = 5, units = "in")
# get just the legend
cowplot::get_legend(t440_merging3_FC) %>% cowplot::ggdraw()
ggsave(filename = here("plots", paste0(out_file_prefix, "LEGEND_tSNE_px440_merging3_by_betareg", ".pdf")), device = cairo_pdf, width = 7, height = 5, units = "in")
#




# Correlation matrix plot -----
CD45_meta30_rcorr <- CD45_meta30_clusters_data %>% 
  select(LabID, cluster = meta30_cluster, percent = perc_meta30) %>% 
  pivot_wider(names_from = cluster, values_from = percent) %>% 
  column_to_rownames(var = "LabID") %>% 
  as.matrix() %>% 
  Hmisc::rcorr(type = "spearman")
#
CD45_meta30_rcorr %>% 
  .$r %>% 
  corrplot::corrplot(method = "color", order = "hclust", addrect = 5)
#
CD45_meta30_hclust <- CD45_meta30_rcorr %>% 
  .$r %>% 
  dist() %>% 
  hclust()
CD45_meta30_hclust %>% str()
CD45_meta30_hclust$order
#
# this does not work yet:
CD45_meta30_rcorr %>% 
  .$r %>% 
  as_tibble(rownames = "meta30") %>% 
  pivot_longer(-meta30, names_to = "cluster", values_to = "rho") %>% 
  mutate(
    meta30 = fct_relevel(meta30, CD45_meta30_hclust$order %>% as.character()),
    cluster = fct_relevel(meta30, CD45_meta30_hclust$order %>% as.character())
  ) %>% 
  inner_join(
    "/Users/mattgalbraith/R_testing/cytof_Workflow/results/P4C_CyTOF_Workflow_v0.1_meta30_clusters_naming_merging.xlsx" %>% 
      read_excel() %>% 
      select(meta30 = meta30_cluster, expert_naming, merging1, merging2) %>% 
      mutate(meta30 = as.character(meta30))
  ) %>% 
  tidyHeatmap::heatmap(meta30, cluster, rho)


#+ save_workspace, include = FALSE
# ##################
# # Save workspace #
# ##################
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)
# ################
# # RESTART HERE #
# ################
# load(here("rdata", paste0(out_file_prefix, ".RData")))


#' 
#' ***
#' ### Session Info
#+ session_info, collapse=TRUE
# Report generated at:
date()

sessionInfo()