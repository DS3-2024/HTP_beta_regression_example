################################################
# Title: beta regression template
# Project: Beta regression analysis of mass cytometry data in R
# Author(s):
#   - author 1
# email(s):
#   - author1@institution.edu
# affiliation(s):
#   - Department of XXX
#   - University of XXX
################################################
# Script original author: Matthew Galbraith
# version: 0.1  Date: 05_23_2024
################################################
# Change Log:
# v0.1
# Initial version
#

### Summary:  
# Beta regression modelling for differential abundance of immune cell clusters/subpopulations as measured by mass cytometry and clustered using FlowSOM.
# See PMID 37379383.
#  

### Data type(s):
#   A. HTP sample meta data
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#   B. HTP CD45+CD66lo-gated immune cell percentage composition data
#      Where/who did this data come from?
#      What is the source of the original data and where is it stored?
#

### Workflow:
#   Step 1 - Read in and inspect sample meta data + immune cell percentage composition data
#   Step 2 - Data exploration
#   Step 3 - Beta regression model setup and assessment
#   Step 4 - Model results
#   Step 5 - Plot individual features
#  

## Comments:
#  Any further relevant details?  
#  


# 0 General Setup -----
# RUN FIRST TIME
# renv::init()
## 0.1 Load required libraries ----
library("readxl") # used to read .xlsx files
library("openxlsx") # used for data export as Excel workbooks
library("tidyverse") # data wrangling and ggplot2
library("ggrepel") # for labelling features
library("ggforce") # for sina plots
library("tictoc") # timer
library("skimr") # data summary
library("janitor") # data cleaning
library("patchwork") # assembling multiple plots
library("gtools") # for logit() and inv.logit() functions
library("betareg") # beta regression
library("conflicted")
conflict_prefer("filter", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("count", "dplyr")
library("here") # generates path to current project directory
# detach("package:here", unload=TRUE) # run this to reset here()
source(here("helper_functions.R")) # load helper functions
#

## 0.2 renv setup ----
# see vignette("renv")
# The general workflow when working with renv is:
#  1 Call renv::init() to initialize a new project-local environment with a private R library,
#  2 Work in the project as normal, installing and removing new R packages as they are needed in the project,
#    recommend using renv::install() as this will automatically update lockfile
#  3 Call renv::snapshot() to save the state of the project library to the lockfile (called renv.lock),
#  4 Continue working on your project, installing and updating R packages as needed.
#  5 Call renv::snapshot() again to save the state of your project library if
# your attempts to update R packages were successful, or call renv::restore() to
# revert to the previous state as encoded in the lockfile if your attempts to
# update packages introduced some new problems.
#

## 0.3 Set required parameters ----
# Input data files
htp_meta_data_file <- here("data", "HTP_Metadata_v0.5_Synapse.txt") # comments/notes on this file?
htp_CD45_FlowSOM_percentage_data_file <- here("data", "HTP_CyTOF_CD45posCD66low_FlowSOM_cluster_percentage_Synapse.txt") # comments/notes on this file?
# Other parameters
# standard_colors <- c("Group1" = "#F8766D", "Group2" = "#00BFC4")
standard_colors <- c("Control" = "gray60", "T21" = "#009b4e")
out_file_prefix <- "DS3_CyTOF_betareg_v0.1_"
# End required parameters ###


# 1 Read in and inspect data ----
## 1.1 Read in meta data ----
htp_meta_data <- htp_meta_data_file |> 
  read_tsv() |> 
  mutate(
    Karyotype = fct_relevel(Karyotype, c("Control", "T21")), # convert to factor and set order
    Sex = fct_relevel(Sex, "Female"), # convert to factor and set order
    Sample_source_code = as_factor(Sample_source_code) # convert to factor - default is numerical order
  )
# inspect
htp_meta_data
htp_meta_data |> skimr::skim()
#
here("data/HTP_Metadata_v0.5_dictionary.txt") |> read_tsv()
#


## 1.2 Read in immune cell percentage composition data ----
htp_CD45_FlowSOM_percentage_data <- htp_CD45_FlowSOM_percentage_data_file |> 
  read_tsv()
# inspect
htp_CD45_FlowSOM_percentage_data # 7,756 rows
htp_CD45_FlowSOM_percentage_data |> skimr::skim()
htp_CD45_FlowSOM_percentage_data |> distinct(Cell_cluster_name) # 20 Cell clusters
htp_CD45_FlowSOM_percentage_data |> distinct(LabID) # 388 LabIDs
#
here("data/HTP_CyTOF_CD45posCD66low_FlowSOM_cluster_dictionary.txt") |> read_tsv()
#

## 1.3 Join meta data with data type 1 and data type 2 ----
htp_meta_CD45_FlowSOM_percentage_data <- htp_CD45_FlowSOM_percentage_data |> 
  inner_join(htp_meta_data, join_by(LabID))
# check number of rows returned !!!


# 2 Data exploration  ----
## 2.1 basic check of data distribution(s) ----
htp_meta_CD45_FlowSOM_percentage_data |> 
  mutate(proportion = Value / 100) |> 
  ggplot(aes(proportion)) + geom_density() + 
  facet_wrap(~Cell_cluster_name, scales = "free") +
  labs(title = "Distributions of cell cluster percentages (untransformed)")
htp_meta_CD45_FlowSOM_percentage_data |> 
  mutate(proportion = Value / 100) |> 
  ggplot(aes(gtools::logit(proportion))) + geom_density() + 
  facet_wrap(~Cell_cluster_name, scales = "free") +
  labs(title = "Distributions of cell cluster percentages (logit-transformed)")
#
htp_meta_CD45_FlowSOM_percentage_data |> 
  mutate(
    proportion = Value / 100,
    logit_proportion = gtools::logit(proportion)
  ) |> 
  ggplot(aes(Karyotype, proportion, color = Karyotype)) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  facet_wrap(~ Cell_cluster_name, scales = "free", nrow = 3) +
  theme(
    aspect.ratio = 1.3,
    strip.text.x = element_text(hjust = 0) # left-align facet labels
  ) +
  labs(
    title = "Distributions of cell cluster percentages (untransformed)"
  )
#
htp_meta_CD45_FlowSOM_percentage_data |> 
  mutate(
    proportion = Value / 100,
    logit_proportion = gtools::logit(proportion)
  ) |> 
  ggplot(aes(Karyotype, logit_proportion, color = Karyotype)) +
  geom_sina() +
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  scale_color_manual(values = standard_colors) +
  facet_wrap(~ Cell_cluster_name, scales = "free", nrow = 3) +
  theme(
    aspect.ratio = 1.3,
    strip.text.x = element_text(hjust = 0) # left-align facet labels
  ) +
  labs(
    title = "Distributions of cell cluster percentages (logit-transformed)"
  )
#

## 2.2 Check for zero values ----
# (problem for beta regression)
htp_meta_CD45_FlowSOM_percentage_data %>% 
  group_by(Cell_cluster_name) %>% 
  summarize(
    n_zeros = sum(Value == 0)
  ) %>% 
  arrange(-n_zeros)
#

## 2.3 Check for extreme outliers ----
# Features with greatest numbers of extreme outliers
htp_meta_CD45_FlowSOM_percentage_data |>  
  mutate(
    proportion = Value / 100,
    logit_proportion = gtools::logit(proportion)
    ) |> 
  group_by(Cell_cluster_name, Karyotype) |>  # important to think about appropriate grouping here
  mutate(extreme = rstatix::is_extreme(logit_proportion)) |>  
  ungroup() |> 
  filter(extreme == TRUE) |>  
  count(Cell_cluster_name, name = "n_extreme") |>  
  arrange(-n_extreme)
# Plot with outliers highlighted
htp_meta_CD45_FlowSOM_percentage_data |> 
  filter(Cell_cluster_name %in% c("non-classical monocytes", "CD56hi NK", "Classical monocytes and M-MDSCs")) |> 
  mutate(
    proportion = Value / 100,
    logit_proportion = gtools::logit(proportion)
  ) |> 
  group_by(Cell_cluster_name, Karyotype) |>  # important to think about appropriate grouping here
  mutate(extreme = rstatix::is_extreme(logit_proportion)) |>  
  ungroup() |> 
  ggplot(aes(Karyotype, logit_proportion)) +
  geom_sina(
    data = . %>% filter(extreme == FALSE),
    aes(color = Karyotype)
  ) +
  geom_sina(
    data = . %>% filter(extreme == TRUE),
    aes(color = "extreme")
  ) +
  geom_boxplot(
    data = . %>% filter(extreme == FALSE),
    notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75
  ) +
  scale_color_manual(values = c(standard_colors, "extreme" = "red")) +
  facet_wrap(~ Cell_cluster_name, scales = "free") +
  theme(aspect.ratio = 1.3) +
  labs(
    title = "Cell clusters with most 'extreme' outliers"
  )
#

## 2.3 Shapiro-Wilk test ----
# very stringent but worth comparing 'best' and 'worst' features
# untransformed
htp_meta_CD45_FlowSOM_percentage_data |>  
  mutate(
    proportion = Value / 100,
    logit_proportion = gtools::logit(proportion)
  ) |> 
  group_by(Cell_cluster_name, Karyotype) |>  # important to think about appropriate grouping here
  mutate(extreme = rstatix::is_extreme(logit_proportion)) |>  
  ungroup() |> 
  filter(extreme == FALSE) |>  
  group_by(Cell_cluster_name) |> 
  rstatix::shapiro_test(proportion) |> 
  arrange(p)
# logit transformed
htp_meta_CD45_FlowSOM_percentage_data |>  
  mutate(
    proportion = Value / 100,
    logit_proportion = gtools::logit(proportion)
  ) |> 
  group_by(Cell_cluster_name, Karyotype) |>  # important to think about appropriate grouping here
  mutate(extreme = rstatix::is_extreme(logit_proportion)) |>  
  ungroup() |> 
  filter(extreme == FALSE) |>  
  group_by(Cell_cluster_name) |> 
  rstatix::shapiro_test(logit_proportion) |> 
  arrange(p)
#
# Q-Q plots (OPTIONAL)
htp_meta_CD45_FlowSOM_percentage_data |> 
  filter(Cell_cluster_name %in% c("naïve CD8+ T", "CD4+ TEM", "CD8a+ γδ T", "CD27+ B")) |> 
  mutate(
    Cell_cluster_name = fct_relevel(Cell_cluster_name, c("naïve CD8+ T", "CD4+ TEM", "CD8a+ γδ T", "CD27+ B")),
    proportion = Value / 100,
    logit_proportion = gtools::logit(proportion)
  ) |> 
  group_by(Cell_cluster_name, Karyotype) |>  # important to think about appropriate grouping here
  mutate(extreme = rstatix::is_extreme(logit_proportion)) |>  
  ungroup() |> 
  ggpubr::ggqqplot("logit_proportion", facet.by = "Cell_cluster_name", title = "Q-Q plot(s)")
#


# 3 Beta regression modelling ----
## 3.1 Set up models for assessment ----
betareg_data <- htp_meta_CD45_FlowSOM_percentage_data |>  
  mutate(
    proportion = Value / 100,
    proportion = if_else(proportion == 0, NA_real_, proportion), # 0 replacement with NA
    logit_proportion = gtools::logit(proportion)
  ) |> 
  group_by(Cell_cluster_name, Karyotype) |>  # important to think about appropriate grouping here
  mutate(extreme = rstatix::is_extreme(logit_proportion)) |>  
  ungroup() |> 
  filter(extreme == FALSE) |>  # remove extreme outliers
  filter(!is.na(proportion)) |> # remove NA / zero values
  # filter by minimum N per group (otherwise map() may fail unhelpfully):
  group_by(Cell_cluster_name) |> 
  add_count(Karyotype, Sex) |> # count by EACH categorical variable
  filter(n >= 15) |> # require at least NN samples in each category
  mutate( # count number of levels for EACH categorical variable
    Karyotype_levels = Karyotype |> fct_drop() |> levels() |> length(),
    Sex_levels = Sex |> fct_drop() |> levels() |> length()
  ) |> 
  ungroup() |> 
  # need to have >1 level for each categorical level or lm() gives error
  filter(Karyotype_levels > 1 & Sex_levels > 1) |> 
  select(LabID, Cell_cluster_name, Value, proportion, logit_proportion, Sex, Age, Karyotype, Sample_source_code) |> 
  nest(data = -Cell_cluster_name)
betareg_data
#
## 3.2 Simple model ----
tic("Running beta regressions for simple model...")
betareg_simple <- betareg_data %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment) # see ?augment.betareg
  )
toc() # Check warnings
betareg_simple
#
# Inspect a model object
betareg_simple |> pluck("fit", 1) # row 1 = Classical monocytes and M-MDSCs
betareg_simple |> pluck("fit", 1) |> class()
betareg_simple |> pluck("fit", 1) |> summary()
betareg_simple |> pluck("fit", 1) |> broom::tidy()
#
# 3.3 Multi-variate model(s) ----
tic("Running beta regressions with Sex + Age + Source")
betareg_multi_SexAgeSource <- betareg_data %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source_code, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment) # see ?augment.betareg
  ) # Check warnings
toc()
betareg_multi_SexAgeSource
#
# 3.4 Variable dispersions ----
tic("Running beta regressions with Sex + Age + Source | variable Source")
betareg_multi_SexAgeSource_varSource <- betareg_data %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source_code | Sample_source_code, data = .x, link = "logit")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment) # see ?augment.betareg
  )
toc() # Check warnings
betareg_multi_SexAgeSource_varSource
#
# 3.5 Alternative link function(s) ----
tic("Running beta regressions with Sex + Age + Source with loglog link")
betareg_multi_SexAgeSource_loglog <- betareg_data %>% 
  mutate(
    fit = map(data, ~ betareg(proportion ~ Karyotype + Sex + Age + Sample_source_code, data = .x, link = "loglog")),
    tidied = map(fit, broom::tidy), # see ?tidy.betareg
    glanced = map(fit, broom::glance), # see ?glance.betareg
    augmented = map(fit, broom::augment) # see ?augment.betareg
  )
toc() # Check warnings
betareg_multi_SexAgeSource_loglog
#


## 3.6 Compare models using AIC ----
# Akaike Information Criteria (AIC): lower values are 'better'.
# Incorporates log-likelihood statistics and Maximum Likelihood Estimation and penalizes complicated models with more covariates.
# i.e. a model with more covariates must 'overcome' the additional complexity to be considered 'better'.
# Often will need some compromise to choose single preferred model specification across all features.
simple_glance <- betareg_simple %>% unnest(glanced)
multi_SexAgeSource_glance <- betareg_multi_SexAgeSource %>% unnest(glanced)
multi_SexAgeSource_glance_varSource <- betareg_multi_SexAgeSource_varSource %>% unnest(glanced)
multi_SexAgeSource_glance_loglog <- betareg_multi_SexAgeSource_loglog %>% unnest(glanced)
#
a1 <- simple_glance %>% select(Cell_cluster_name, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_glance %>% select(Cell_cluster_name, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(Cell_cluster_name = fct_inorder(Cell_cluster_name)) %>% 
  ggplot(aes(Cell_cluster_name, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # turn off with too many
  labs(title = "simple  vs.\n multi_SexAgeSource model") +
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = Cell_cluster_name), size = 3, nudge_x = 2)
#
a2 <- multi_SexAgeSource_glance %>% select(Cell_cluster_name, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_glance_varSource %>% select(Cell_cluster_name, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(Cell_cluster_name = fct_inorder(Cell_cluster_name)) %>% 
  ggplot(aes(Cell_cluster_name, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # turn off with too many
  labs(title = "multi_SexAgeSource  vs.\n multi_SexAgeSource|varSource model") +
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = Cell_cluster_name), size = 3, nudge_x = 2)
#
a3 <- multi_SexAgeSource_glance %>% select(Cell_cluster_name, AIC1 = AIC) %>% 
  inner_join(multi_SexAgeSource_glance_loglog %>% select(Cell_cluster_name, AIC2 = AIC)) %>% 
  mutate(AIC_diff = AIC1 - AIC2) %>% 
  arrange(-AIC_diff) %>% 
  mutate(Cell_cluster_name = fct_inorder(Cell_cluster_name)) %>% 
  ggplot(aes(Cell_cluster_name, AIC_diff)) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 10, linetype = 2) +
  geom_hline(yintercept = -10, linetype = 2) +
  geom_point() +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + # turn off with too many
  labs(title = "multi_SexAgeSource  vs.\n multi_SexAgeSource|loglog model") +
  geom_text_repel(data = . %>%  slice_max(abs(AIC_diff), n = 6), aes(label = Cell_cluster_name), size = 3, nudge_x = 2)
#
a1 + a2 + a3 + plot_spacer() + plot_annotation(title = "Beta regression: AIC comparisons")
#
a1 + expand_limits(y = c(NA, 157)) +
a2 + expand_limits(y = c(NA, 157)) +
a3 + expand_limits(y = c(NA, 157)) +
  plot_spacer() + plot_annotation(title = "Beta regression: AIC comparisons")
#

# 4 Model results ----
## 4.1 Extract model results for simple model ----
betareg_simple_results <- betareg_simple |> 
  unnest(tidied) |> 
  filter(str_detect(term, "Karyotype")) |> 
  transmute(
    Cell_cluster_name, 
    FoldChange = estimate |> exp(), # check for transformation and adjust accordingly
    pval = p.value,
    BHadj_pval = p.adjust(pval, method = "BH"),
    Model = "betareg(proportion ~ Karyotype)" # Update accordingly
  ) |> 
  arrange(pval)
# Volcano plot
v1 <- betareg_simple_results %>% 
  volcano_plot_lab_lm(
    title = "Cell clusters: T21 vs. Control",
    subtitle = paste0("betareg: proportion ~ Karyotype\n","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]")
  )
v1
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_simple", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
# ggrastr::rasterize(., layers='Point', dpi = 600, dev = "ragg_png")
#
# Export results as tab-delimited text
betareg_simple_results |> 
  write_tsv(file = here("results", paste0(out_file_prefix, "betareg_simple_results", ".txt")))
# Export results as Excel (input must be named list)
list(
  "betareg_simple" = betareg_simple_results
) |> 
  export_excel(filename = "betareg_simple_results")
#

## 4.2 Extract model results for multi model ----
betareg_multi_SexAgeSource_results <- betareg_multi_SexAgeSource |> 
  unnest(tidied) |> 
  filter(str_detect(term, "Karyotype")) |> 
  transmute(
    Cell_cluster_name, 
    FoldChange = estimate |> exp(), # check for transformation and adjust accordingly
    pval = p.value,
    BHadj_pval = p.adjust(pval, method = "BH"),
    Model = "betareg(proportion ~ Karyotype+Sex+Age+Source)" # Update accordingly
  ) |> 
  arrange(pval)
# Volcano plot
v2 <- betareg_multi_SexAgeSource_results %>% 
  volcano_plot_lab_lm(
    title = "Cell clusters: T21 vs. Control",
    subtitle = paste0("betareg: proportion ~ Karyotype+Sex+Age+Source\n","[Down: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange < 1) %>% nrow(), "; Up: ", (.) %>% filter(BHadj_pval < 0.1 & FoldChange > 1) %>% nrow(), "]")
  )
v2
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi", ".png")), width = 5, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_multi", ".pdf")), device = cairo_pdf, width = 5, height = 5, units = "in")
# ggrastr::rasterize(., layers='Point', dpi = 600, dev = "ragg_png")
#
# Export results as tab-delimited text
betareg_multi_SexAgeSource_results |> 
  write_tsv(file = here("results", paste0(out_file_prefix, "betareg_multi_SexAgeSource_results", ".txt")))
# Export results as Excel (input must be named list)
list(
  "betareg_multi_SexAgeSource" = betareg_multi_SexAgeSource_results
) |> 
  export_excel(filename = "betareg_multi_SexAgeSource_results")
#


## 4.3 Arrange multiple plots using patchwork ----
v1 + v2 + plot_layout(guides = "collect", nrow = 1) # different y scales
#
v1 + expand_limits(y = c(NA, 52)) + # same y scale
  v2 + expand_limits(y = c(NA, 52)) +
  plot_layout(guides = "collect", nrow = 1)
ggsave(filename = here("plots", paste0(out_file_prefix, "volcano_combined", ".pdf")), device = cairo_pdf, width = 12, height = 5, units = "in")
#


# 5 Plot individual features ----
## 5.1 Get interesting/significant features ----
top_signif_by_FC <- bind_rows(
  betareg_multi_SexAgeSource_results |> 
    filter(BHadj_pval<0.1) |> 
    arrange(-FoldChange) |> 
    slice_max(order_by = FoldChange, n = 5), # Upregulated
  betareg_multi_SexAgeSource_results |> 
    filter(BHadj_pval<0.1) |> 
    arrange(FoldChange) |> 
    slice_min(order_by = FoldChange, n = 5) # Downregulated
)
#
## 5.2 Sina plots ----
# sina plots use horizontal jitter of data points to show density
# 
# With logit tramsformed data
s1 <- betareg_data |> 
  unnest(data) |> # extreme outliers already removed
  filter(Cell_cluster_name %in% top_signif_by_FC$Cell_cluster_name) |> # filter to features of interest
  mutate(Cell_cluster_name = fct_relevel(Cell_cluster_name, top_signif_by_FC$Cell_cluster_name)) |> # control plotting order
  ggplot(aes(Karyotype, logit_proportion, color = Karyotype)) +
  geom_sina() + 
  geom_boxplot(notch = TRUE, varwidth = FALSE, outlier.shape = NA, coef = FALSE, width = 0.3, color = "black", fill = "transparent", size = 0.75) +
  facet_wrap(~ Cell_cluster_name, scales = "free_y", nrow = 2) + # facet per feature; each feature on it's own scale
  scale_color_manual(values = standard_colors) + # use standardized colors
  theme(aspect.ratio = 1.3) + # set fixed aspect ratio
  labs(
    title = "Top significant cell clusters by fold-change: T21 vs. Control",
    subtitle = "Unadjusted data; extreme outliers removed"
  )
s1
#
# ggrastr::rasterize(., layers='Point', dpi = 600, dev = "ragg_png")
#
# With transformation of labels back to original percentage scale
s2 <- s1+
  scale_y_continuous(labels=function(x)round(gtools::inv.logit(x) * 100, digits=3), breaks = seq(-12, 12, by = 1)) +
  labs(
    y = "% among CD45+CD66lo" # need to change y-axis label to match
  )
s2
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_top_signif_by_FC", ".png")), width = 15, height = 5, units = "in")
ggsave(filename = here("plots", paste0(out_file_prefix, "sina_top_signif_by_FC", ".pdf")), device = cairo_pdf, width = 15, height = 5, units = "in")
#




################################################
# save workspace ----
save.image(file = here("rdata", paste0(out_file_prefix, ".RData")), compress = TRUE, safe = TRUE) # saves entire workspace (can be slow)
# To reload previously saved workspace:
# load(here("rdata", paste0(out_file_prefix, ".RData")))

# session_info ----
date()
sessionInfo()
################################################

