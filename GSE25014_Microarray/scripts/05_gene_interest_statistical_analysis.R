# ============================================================
# Lesson 05: Statistical Analysis of Genes of Interest
# ============================================================
#
# Dataset:
#   GSE25014
#
# Platform:
#   Affymetrix Human Genome U133 Plus 2.0
#
# Objective:
#
#   Evaluate the expression of selected genes of interest
#   using the gene-level expression matrix generated in the
#   previous lessons and the original sample metadata
#   downloaded from GEO.
#
#
# Experimental design:
#
#   Cell types:
#
#       1. pulmonary artery endothelial cells
#       2. pulmonary microvascular endothelial cells
#
#
#   Treatments:
#
#       control
#       heme
#
#
# IMPORTANT:
#
#   The two cell types will be analyzed SEPARATELY.
#
#   Comparison 1:
#
#       PAEC:
#       control vs heme
#
#   Comparison 2:
#
#       PMVEC:
#       control vs heme
#
#
# Statistical workflow:
#
#   Gene expression
#          ↓
#   GEO metadata
#          ↓
#   Sample identification
#          ↓
#   Cell type and treatment
#          ↓
#   Genes of interest
#          ↓
#   Descriptive statistics
#          ↓
#   Normality test
#          ↓
#   Statistical test
#          ↓
#   Raw p-value
#          ↓
#   Benjamini-Hochberg correction
#          ↓
#   FDR
#          ↓
#   Visualization
#
#
# Statistical tests:
#
#   Both groups approximately normal:
#       Welch's t-test
#
#   At least one group not approximately normal:
#       Wilcoxon rank-sum test
#
#
# Multiple testing correction:
#
#   Benjamini-Hochberg (BH)
#
# ============================================================


# ============================================================
# Block 0: Load Required Packages
# ============================================================
library(tidyverse)
library(ggpubr)
library(rstatix)

# ============================================================
# Block 1: Load Gene-Level Expression Data
# ============================================================

# Gene-level expression matrix generated in previous lessons.
gene_expression <- readRDS("data/processed/GSE25014_expression_gene_level.rds")

# ------------------------------------------------------------
# Inspect expression object
# ------------------------------------------------------------
print(class(gene_expression))
print(dim(gene_expression))
print(head(gene_expression[, 1:5]))

# ------------------------------------------------------------
# Number of genes and samples
# ------------------------------------------------------------
print(paste("Number of genes:", nrow(gene_expression)))
print(paste("Number of samples:", ncol(gene_expression)))

# ------------------------------------------------------------
# Inspect gene identifiers
# ------------------------------------------------------------
print("First gene identifiers:")
print(head(rownames(gene_expression)))

# ------------------------------------------------------------
# Check duplicated gene identifiers
# ------------------------------------------------------------
duplicated_genes <- unique(
  rownames(gene_expression)[duplicated(rownames(gene_expression))])

print(paste("Duplicated gene identifiers:", length(duplicated_genes)))

if(length(duplicated_genes) > 0){
  print(duplicated_genes)}

# ============================================================
# Block 2: Load Original GEO Metadata
# ============================================================

# The experimental information comes directly from GEO.
#
# We will NOT manually create the experimental design.
#
# Cell type and treatment are extracted from the original
# GEO metadata.
metadata <- read.csv("data/metadata/GSE25014_metadata.csv", check.names = FALSE)

# ------------------------------------------------------------
# Inspect metadata
# ------------------------------------------------------------
print(class(metadata))
print(dim(metadata))
print(colnames(metadata))
print(head(metadata))

# ============================================================
# Block 3: Inspect Cell Types Available in GEO
# ============================================================
print("Cell types available in GEO:")
print(unique(metadata$`cell type:ch1`))

# ============================================================
# Block 4: Inspect Treatments Available in GEO
# ============================================================
print("Treatments available in GEO:")
print(unique(metadata$`treatment:ch1`))

# ============================================================
# Block 5: Inspect Complete Experimental Design
# ============================================================

# This table allows us to verify how many samples belong
# to each combination of cell type and treatment.
experimental_design <- table(metadata$`cell type:ch1`, metadata$`treatment:ch1`)
print("Experimental design obtained from GEO:")
print(experimental_design)

# ============================================================
# Block 6: Create Clean Metadata
# ============================================================

# geo_accession is the official GEO sample identifier.
#
# Example:
#
#   GSM614501
#   GSM614502
#   GSM614503
#
# These identifiers will be used to connect the expression
# matrix with the GEO metadata.

metadata_clean <- metadata %>%
  select(
    SampleID = geo_accession,
    Title = title,
    CellType = `cell type:ch1`,
    Treatment = `treatment:ch1`) %>%
  mutate(
    SampleID = str_trim(SampleID),
    Title = str_trim(Title),
    CellType = str_trim(CellType),
    Treatment = str_trim(Treatment))

# ------------------------------------------------------------
# Inspect cleaned metadata
# ------------------------------------------------------------
print(metadata_clean)

# ============================================================
# Block 7: Standardize Expression Sample IDs
# ============================================================

# The expression matrix contains names such as:
#
#   GSM614501_Sample_1.CEL
#
# while GEO metadata contains:
#
#   GSM614501
#
# Therefore, we remove the file-specific suffix.
expression_samples <- colnames(gene_expression)
print("Original expression sample names:")
print(expression_samples)

# ------------------------------------------------------------
# Remove "_Sample_X.CEL"
# ------------------------------------------------------------
expression_samples_clean <- sub("_Sample_[0-9]+\\.CEL$", "", expression_samples)
print("Standardized expression sample IDs:")
print(expression_samples_clean)

# ============================================================
# Block 8: Create Sample Mapping Table
# ============================================================
sample_map <- tibble(ExpressionSample = expression_samples,
                     SampleID = expression_samples_clean)
print(sample_map)

# ============================================================
# Block 9: Validate Sample IDs
# ============================================================

# ------------------------------------------------------------
# Check that all expression samples exist in GEO metadata
# ------------------------------------------------------------
all_expression_in_metadata <- all(expression_samples_clean %in%
                                    metadata_clean$SampleID)

print(paste("All expression samples found in GEO metadata:",
            all_expression_in_metadata))

# ------------------------------------------------------------
# Identify missing samples
# ------------------------------------------------------------
missing_metadata <- expression_samples_clean[!expression_samples_clean %in%
                                               metadata_clean$SampleID]

if(length(missing_metadata) > 0){
  print("Expression samples missing from GEO metadata:")
  print(missing_metadata)}

# ------------------------------------------------------------
# Check duplicated GEO IDs
# ------------------------------------------------------------
duplicated_metadata_ids <- metadata_clean$SampleID[duplicated(
  metadata_clean$SampleID)]

print(paste("Duplicated GEO sample IDs:", length(duplicated_metadata_ids)))

if(length(duplicated_metadata_ids) > 0){
  print(duplicated_metadata_ids)}

# ============================================================
# Block 10: Reorder Metadata According to Expression Matrix
# ============================================================

# The order of the metadata must correspond exactly to
# the order of the expression matrix columns.
metadata_analysis <- metadata_clean %>%
  slice(match(expression_samples_clean, SampleID))

# ------------------------------------------------------------
# Inspect reordered metadata
# ------------------------------------------------------------
print(metadata_analysis)

# ------------------------------------------------------------
# Verify sample order
# ------------------------------------------------------------
same_sample_order <- all(expression_samples_clean == metadata_analysis$SampleID)
print(paste("Expression and metadata have the same sample order:",
            same_sample_order))

# ------------------------------------------------------------
# Stop if samples do not match
# ------------------------------------------------------------
if(!same_sample_order){
  stop("ERROR: Expression matrix and GEO metadata do not match.")}

# ============================================================
# Block 11: Inspect Experimental Design After Matching
# ============================================================
print("Experimental design used in the analysis:")
print(table(metadata_analysis$CellType, metadata_analysis$Treatment))

# ============================================================
# Block 12: Select the Two Cell Types of Interest
# ============================================================

# We explicitly select the two endothelial cell types
# analyzed in this lesson.
cell_types_interest <- c("pulmonary artery endothelial cells",
                         "pulmonary microvascular endothelial cells")

print("Cell types selected for analysis:")
print(cell_types_interest)

# ------------------------------------------------------------
# Verify that the cell types exist
# ------------------------------------------------------------
cell_types_missing <- cell_types_interest[!cell_types_interest %in%
                                            metadata_analysis$CellType]

if(length(cell_types_missing) > 0){
  stop(paste("ERROR: Cell type not found:", 
             paste(cell_types_missing, collapse = ", ")))}

# ============================================================
# Block 13: Select Treatments
# ============================================================

# According to the GEO metadata, the comparison of interest
# is:
#
#       control vs heme
treatments_interest <- c("control", "heme")
print("Treatments selected:")
print(treatments_interest)

# ------------------------------------------------------------
# Verify treatments
# ------------------------------------------------------------
treatments_missing <- treatments_interest[!treatments_interest %in%
                                            metadata_analysis$Treatment]

if(length(treatments_missing) > 0){
  stop(paste("ERROR: Treatment not found:", 
             paste(treatments_missing, collapse = ", ")))}

# ============================================================
# Block 14: Define Genes of Interest
# ============================================================
genes_interest <- c("CDH5", "TJP1", "CTTN", "PTK2", "HMOX1", "NFE2L2", "GPX4",
                    "SLC7A11", "ACSL4", "TFRC", "FTH1", "NCOA4", "AIFM2",
                    "SOD2", "NQO1", "AKT1", "MAPK1", "MAPK3", "BMPR2", "TGFB1")

print(paste("Number of genes of interest:", length(genes_interest)))

# ============================================================
# Block 15: Check Gene Availability
# ============================================================
genes_found <- genes_interest[genes_interest %in% rownames(gene_expression)]
genes_missing <- genes_interest[!genes_interest %in% rownames(gene_expression)]

# ------------------------------------------------------------
# Genes found
# ------------------------------------------------------------
print("Genes found in the expression matrix:")
print(genes_found)

# ------------------------------------------------------------
# Genes missing
# ------------------------------------------------------------
print("Genes missing from the expression matrix:")
print(genes_missing)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
print(paste("Genes found:", length(genes_found),
            "of", length(genes_interest)))

# ============================================================
# Block 16: Create Gene-Level Analysis Matrix
# ============================================================
expression_interest <- gene_expression[genes_found, ,drop = FALSE]
print(dim(expression_interest))
print(head(expression_interest))

# ============================================================
# Block 17: Convert Expression Matrix to Long Format
# ============================================================

# Original matrix:
#
#              GSM1   GSM2   GSM3
#   CDH5       ...    ...    ...
#   HMOX1      ...    ...    ...
#
#
# Long format:
#
#   Gene   SampleID   Expression
#
# This format is easier to use with ggplot2 and statistical
# analysis.
expression_long <- expression_interest %>%
  as.data.frame() %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene,
               names_to = "ExpressionSample",
               values_to = "Expression") %>%
  mutate(SampleID = sub("_Sample_[0-9]+\\.CEL$", "", ExpressionSample))

# ============================================================
# Block 18: Add GEO Metadata
# ============================================================
expression_long <- expression_long %>%
  left_join(metadata_analysis, by = "SampleID")

# ------------------------------------------------------------
# Inspect final analysis table
# ------------------------------------------------------------
print(head(expression_long))
print(dim(expression_long))

# ============================================================
# Block 19: Filter Cell Types of Interest
# ============================================================

# This is the key step that ensures PAECs and PMVECs
# are analyzed separately.
expression_long <- expression_long %>%
  filter(CellType %in% cell_types_interest)

# ------------------------------------------------------------
# Check number of samples per cell type and treatment
# ------------------------------------------------------------
print(expression_long %>%
        distinct(SampleID, CellType, Treatment) %>%
        count(CellType, Treatment))

# ============================================================
# Block 20: Check for Missing Metadata
# ============================================================
missing_metadata_rows <- expression_long %>%
  filter(is.na(CellType) | is.na(Treatment))

print(paste("Rows with missing experimental information:",
            nrow(missing_metadata_rows)))

if(nrow(missing_metadata_rows) > 0){
  print(missing_metadata_rows)
  stop("ERROR: Some expression samples do not have GEO metadata.")}

# ============================================================
# Block 21: Convert Treatment to Factor
# ============================================================
expression_long <- expression_long %>%
  mutate(Treatment = factor(Treatment, levels = c("control", "heme")))

print(levels(expression_long$Treatment))

# ============================================================
# Block 22: Expression Summary
# ============================================================

# Descriptive statistics are calculated independently for:
#
#   Gene
#   Cell type
#   Treatment

expression_summary <- expression_long %>%
  group_by(Gene, CellType, Treatment) %>%
  summarise(n = sum(!is.na(Expression)),
            Mean = mean(Expression, na.rm = TRUE),
            Median = median(Expression, na.rm = TRUE),
            SD = sd(Expression, na.rm = TRUE),
            Min = min(Expression, na.rm = TRUE),
            Max = max(Expression, na.rm = TRUE),
            .groups = "drop")

print(expression_summary)

# ============================================================
# Block 23: Test Normality
# ============================================================

# Shapiro-Wilk test
#
# H0:
#   The data follow a normal distribution.
#
# H1:
#   The data do not follow a normal distribution.
#
#
# p > 0.05:
#   No statistical evidence against normality.
#
#
# p <= 0.05:
#   Evidence against normality.
#
#
# IMPORTANT:
#
# The sample size per group is small.
#
# Therefore, the Shapiro-Wilk test should NOT be interpreted
# alone.
#
# Q-Q plots should also be examined.
normality_results <- expression_long %>%
  group_by(Gene, CellType, Treatment) %>%
  shapiro_test(Expression) %>%
  ungroup()

print(normality_results)

# ============================================================
# Block 24: Classify Normality
# ============================================================
normality_results <- normality_results %>%
  mutate(Normality = case_when(p > 0.05 ~ "Approximately normal",
                               p <= 0.05 ~ "Evidence against normality",
                               TRUE ~ NA_character_))

print(normality_results)

# ============================================================
# Block 25: Q-Q Plot for PAECs
# ============================================================
qq_paec <- expression_long %>%
  filter(CellType == "pulmonary artery endothelial cells") %>%
  filter(Gene == "CDH5") %>%
  ggplot(aes(sample = Expression)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Treatment) +
  labs(title = "CDH5 - PAECs - Q-Q Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_classic()

print(qq_paec)

# ============================================================
# Block 26: Q-Q Plot for PMVECs
# ============================================================
qq_pmv_ec <- expression_long %>%
  filter(CellType == "pulmonary microvascular endothelial cells") %>%
  filter(Gene == "CDH5") %>%
  ggplot(aes(sample = Expression)) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ Treatment) +
  labs(title = "CDH5 - PMVECs - Q-Q Plot",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_classic()

print(qq_pmv_ec)

# ============================================================
# Block 27: Function to Select Statistical Test
# ============================================================

# The comparison is performed separately for each cell type.
#
#
# PAEC:
#
#   control vs heme
#
#
# PMVEC:
#
#   control vs heme
#
#
# Decision:
#
#   Both groups approximately normal
#             ↓
#       Welch's t-test
#
#
#   At least one group not approximately normal
#             ↓
#       Wilcoxon rank-sum test


run_gene_test <- function(data){
  # ----------------------------------------------------------
  # Get treatment groups
  # ----------------------------------------------------------
  treatments <- unique(data$Treatment[!is.na(data$Treatment)])
  
  # ----------------------------------------------------------
  # Check that exactly two treatments exist
  # ----------------------------------------------------------
  if(length(treatments) != 2){
    return(tibble(Treatment1 = NA_character_,
                  Treatment2 = NA_character_,
                  N1 = NA_integer_,
                  N2 = NA_integer_,
                  Normality_p1 = NA_real_,
                  Normality_p2 = NA_real_,
                  Test = NA_character_,
                  p = NA_real_))}
  
  # ----------------------------------------------------------
  # Extract expression values
  # ----------------------------------------------------------
  group1 <- data %>%
    filter(Treatment == treatments[1]) %>%
    pull(Expression) %>%
    na.omit()
  
  group2 <- data %>%
    filter(Treatment == treatments[2]) %>%
    pull(Expression) %>%
    na.omit()
  
  # ----------------------------------------------------------
  # Check sample size
  # ----------------------------------------------------------
  if(length(group1) < 3 | length(group2) < 3){
    return(tibble(Treatment1 = treatments[1],
                  Treatment2 = treatments[2],
                  N1 = length(group1),
                  N2 = length(group2),
                  Normality_p1 = NA_real_,
                  Normality_p2 = NA_real_,
                  Test = "Insufficient observations",
                  p = NA_real_))}
  
  # ----------------------------------------------------------
  # Shapiro-Wilk test
  # ----------------------------------------------------------
  shapiro1 <- shapiro.test(group1)
  shapiro2 <- shapiro.test(group2)
  
  # ----------------------------------------------------------
  # Select statistical test
  # ----------------------------------------------------------
  if(shapiro1$p.value > 0.05 & shapiro2$p.value > 0.05){
    
    # --------------------------------------------------------
    # Welch's t-test
    # --------------------------------------------------------
    test_result <- t.test(group1, group2, paired = FALSE, var.equal = FALSE)
    test_name <- "Welch t-test"
    } else {
    
    # --------------------------------------------------------
    # Wilcoxon rank-sum test
    # --------------------------------------------------------
    test_result <- wilcox.test(group1, group2, paired = FALSE, exact = FALSE)
    test_name <- "Wilcoxon"}
  
  # ----------------------------------------------------------
  # Return results
  # ----------------------------------------------------------
  tibble(Treatment1 = treatments[1],
         Treatment2 = treatments[2],
         N1 = length(group1),
         N2 = length(group2),
         Normality_p1 = shapiro1$p.value,
         Normality_p2 = shapiro2$p.value,
         Test = test_name,
         p = test_result$p.value)}

# ============================================================
# Block 28: Run Statistical Analysis
# ============================================================

# IMPORTANT:
#
# group_by(Gene, CellType)
#
# ensures that each cell type is analyzed independently.
#
# Therefore:
#
#   CDH5 + PAEC
#       → control vs heme
#
#   CDH5 + PMVEC
#       → control vs heme
#
# are two independent statistical comparisons.
statistical_results <- expression_long %>%
  group_by(Gene, CellType) %>%
  group_modify(~ run_gene_test(.x)) %>%
  ungroup()

print(statistical_results)

# ============================================================
# Block 29: Multiple Testing Correction
# ============================================================

# We are performing:
#
#   20 genes
#       ×
#   2 cell types
#
# = 40 statistical comparisons.
#
#
# The raw p-values are therefore corrected using
# Benjamini-Hochberg.
#
# The correction is applied to the complete set of
# comparisons.
statistical_results <- statistical_results %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

print(statistical_results)

# ============================================================
# Block 30: Classify Statistical Significance
# ============================================================
statistical_results <- statistical_results %>%
  mutate(Significance = case_when(p_adj < 0.001 ~ "***",
                                  p_adj < 0.01 ~ "**",
                                  p_adj < 0.05 ~ "*",
                                  TRUE ~ "ns"))

print(statistical_results)

# ============================================================
# Block 31: Sort Statistical Results
# ============================================================
statistical_results_sorted <- statistical_results %>%
  arrange(p_adj)

print(statistical_results_sorted)

# ============================================================
# Block 32: Display Results by Cell Type
# ============================================================

# ------------------------------------------------------------
# PAEC results
# ------------------------------------------------------------
paec_results <- statistical_results %>%
  filter(CellType == "pulmonary artery endothelial cells") %>%
  arrange(p_adj)

print("==============================================")
print("PAEC - RESULTS")

print("==============================================")
print(paec_results)

# ------------------------------------------------------------
# PMVEC results
# ------------------------------------------------------------
pmvec_results <- statistical_results %>%
  filter(CellType == "pulmonary microvascular endothelial cells") %>%
  arrange(p_adj)

print("==============================================")
print("PMVEC - RESULTS")

print("==============================================")

print(pmvec_results)

# ============================================================
# Block 33: Export Statistical Results
# ============================================================
write.csv(statistical_results,
          "results/GSE25014_gene_interest_statistics.csv",
          row.names = FALSE)

# ------------------------------------------------------------
# Export PAEC results
# ------------------------------------------------------------
write.csv(paec_results,
  "results/GSE25014_PAEC_statistics.csv",
  row.names = FALSE)

# ------------------------------------------------------------
# Export PMVEC results
# ------------------------------------------------------------
write.csv(
  pmvec_results,
  "results/GSE25014_PMVEC_statistics.csv",
  row.names = FALSE)

# ============================================================
# Block 34: Prepare Data for Visualization
# ============================================================
plot_data <- expression_long %>%
  mutate(Treatment = factor(Treatment, levels = c("control", "heme")))

# ============================================================
# Block 35: Function to Generate Gene Plot
# ============================================================
gene_plot <- function(gene_name, cell_type){
  # ----------------------------------------------------------
  # Select expression data
  # ----------------------------------------------------------
  plot_df <- plot_data %>%
    filter(Gene == gene_name,
           CellType == cell_type)
  
  # ----------------------------------------------------------
  # Select statistical result
  # ----------------------------------------------------------
  stat_df <- statistical_results %>%
    filter(Gene == gene_name,
           CellType == cell_type)
  
  # ----------------------------------------------------------
  # Check data
  # ----------------------------------------------------------
  if(nrow(plot_df) == 0){
    stop("Gene or cell type not found.")}
  
  if(nrow(stat_df) == 0){
    stop("Statistical result not found.")}
  
  # ----------------------------------------------------------
  # Create statistical label
  # ----------------------------------------------------------
  p_label <- paste0("Test: ", stat_df$Test, "\n", "p = ",
                    format.pval(stat_df$p,
                                digits = 3,
                                eps = 0.001),
                    "\n", "FDR = ", 
                    format.pval(stat_df$p_adj,
                                digits = 3,
                                eps = 0.001))
  
  # ----------------------------------------------------------
  # Determine annotation position
  # ----------------------------------------------------------
  y_max <- max(plot_df$Expression, na.rm = TRUE)
  y_min <- min(plot_df$Expression, na.rm = TRUE)
  y_range <- y_max - y_min
  
  if(y_range == 0){
    y_range <- 1
  }
  
  y_position <- y_max + 0.20 * y_range
  
  # ----------------------------------------------------------
  # Generate plot
  # ----------------------------------------------------------
  p <- ggplot(plot_df, aes(x = Treatment, y = Expression)) +
    geom_boxplot(width = 0.6, outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 2) +
    annotate("text", x = 1.5, y = y_position, label = p_label, size = 4) +
    labs(title = paste(gene_name, "-", cell_type),
         x = "Treatment",
         y = "Gene expression") +
    theme_classic()
  return(p)
}

# ============================================================
# Block 36: Example Plot - CDH5 in PAEC
# ============================================================
cdh5_paec <- gene_plot(gene_name = "CDH5",
                       cell_type = "pulmonary artery endothelial cells")

print(cdh5_paec)

# ============================================================
# Block 37: Example Plot - CDH5 in PMVEC
# ============================================================
cdh5_pmvec <- gene_plot(gene_name = "CDH5",
                        cell_type = "pulmonary microvascular endothelial cells")


print(cdh5_pmvec)

# ============================================================
# Block 38: Example Plot - HMOX1 in PAEC
# ============================================================
hmox1_paec <- gene_plot(gene_name = "HMOX1",
                        cell_type = "pulmonary artery endothelial cells")

print(hmox1_paec)

# ============================================================
# Block 39: Example Plot - HMOX1 in PMVEC
# ============================================================
hmox1_pmvec <- gene_plot(gene_name = "HMOX1",
                         cell_type = "pulmonary microvascular endothelial cells")

print(hmox1_pmvec)


# ============================================================
# Block 40: Save Example Figures
# ============================================================
ggsave(filename = "results/GSE25014_CDH5_PAEC.png",
       plot = cdh5_paec,
       width = 7,
       height = 5,
       dpi = 300)

ggsave(filename = "results/GSE25014_CDH5_PMVEC.png",
       plot = cdh5_pmvec,
       width = 7,
       height = 5,
       dpi = 300)

ggsave(filename = "results/GSE25014_HMOX1_PAEC.png",
       plot = hmox1_paec,
       width = 7,
       height = 5,
       dpi = 300)

ggsave(filename = "results/GSE25014_HMOX1_PMVEC.png",
       plot = hmox1_pmvec,
       width = 7,
       height = 5,
       dpi = 300)

# ============================================================
# Block 41: Create Directory for All Gene Plots
# ============================================================
plot_directory <- "results/GSE25014_gene_plots"

if(!dir.exists(plot_directory)){
  dir.create(plot_directory, recursive = TRUE)}

# ============================================================
# Block 42: Generate Plots for All Genes and Cell Types
# ============================================================

# One figure will be generated for each:
#
#   Gene × Cell Type
#
#
# Therefore:
#
#   20 genes × 2 cell types
#
# = 40 figures.
for(gene_name in genes_found){
  for(cell_type in cell_types_interest){
    
    # --------------------------------------------------------
    # Generate plot
    # --------------------------------------------------------
    plot_current <- gene_plot(gene_name = gene_name,
                              cell_type = cell_type)
    
    # --------------------------------------------------------
    # Create safe file names
    # --------------------------------------------------------
    safe_gene <- str_replace_all(gene_name, "[^A-Za-z0-9]", "_")
    safe_cell <- str_replace_all(cell_type, "[^A-Za-z0-9]", "_")
    
    # --------------------------------------------------------
    # Create file name
    # --------------------------------------------------------
    file_name <- paste0(plot_directory, "/", safe_gene, "_", safe_cell, ".png")
    
    # --------------------------------------------------------
    # Save plot
    # --------------------------------------------------------
    ggsave(filename = file_name,
           plot = plot_current,
           width = 7,
           height = 5,
           dpi = 300)}}

# ============================================================
# Block 43: Identify Significant Results
# ============================================================

significant_results <- statistical_results %>%
  filter(p_adj < 0.05) %>%
  arrange(p_adj)

print("==============================================")
print("SIGNIFICANT RESULTS")
print("==============================================")

print(significant_results)

# ============================================================
# Block 44: Significant Results by Cell Type
# ============================================================

# ------------------------------------------------------------
# PAEC
# ------------------------------------------------------------
significant_paec <- statistical_results %>%
  filter(CellType == "pulmonary artery endothelial cells", p_adj < 0.05) %>%
  arrange(p_adj)

print("Significant genes in PAEC:")
print(significant_paec)

# ------------------------------------------------------------
# PMVEC
# ------------------------------------------------------------
significant_pmvec <- statistical_results %>%
  filter(CellType == "pulmonary microvascular endothelial cells", p_adj < 0.05) %>%
  arrange(p_adj)

print("Significant genes in PMVEC:")

print(significant_pmvec)

# ============================================================
# Block 45: Save Significant Results
# ============================================================
write.csv(significant_results,
          "results/GSE25014_significant_gene_results.csv",
          row.names = FALSE)

# ============================================================
# Block 46: Save Final Analysis Objects
# ============================================================
saveRDS(expression_long,
        "data/processed/GSE25014_gene_interest_expression_long.rds")

saveRDS(statistical_results,
        "data/processed/GSE25014_gene_interest_statistics.rds")

saveRDS(normality_results,
        "data/processed/GSE25014_gene_interest_normality.rds")

# ============================================================
# End of Lesson 05
# ============================================================