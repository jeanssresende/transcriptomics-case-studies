# ============================================================
# Transcriptomics Case Studies
# Case Study: GSE25014
#
# Lesson 04: Probe-to-Gene Processing
# ============================================================
#
# Dataset:
#   GSE25014 - Gene expression data of endothelium exposed to heme
#
# Platform:
#   GPL570 - Affymetrix Human Genome U133 Plus 2.0 Array
#
# Objective:
#   Convert the RMA probe-set expression matrix into a
#   gene-level expression matrix.
#
# Strategy:
#
#   1. Load the RMA expression matrix
#   2. Load GPL570 probe annotation
#   3. Remove probes without gene annotation
#   4. Remove probes mapping to multiple genes
#   5. Keep probes mapping to exactly one gene
#   6. Group multiple probes representing the same gene
#   7. Summarize probes using the median expression
#   8. Generate a gene × sample matrix
#
# Important:
#
#   We do NOT select the probe with the highest CV.
#
#   Multiple probes for the same gene will be summarized
#   using the median expression value.
#
# ============================================================


# ------------------------------------------------------------
# 1. Install required packages
# ------------------------------------------------------------

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}


# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

library(dplyr)
library(readr)


# ------------------------------------------------------------
# 3. Define directories
# ------------------------------------------------------------

case_study_dir <- "GSE25014_Microarray"

if (dir.exists(case_study_dir)) {
  base_dir <- case_study_dir
} else {
  base_dir <- "."
}


processed_dir <- file.path(
  base_dir,
  "data",
  "processed"
)


metadata_dir <- file.path(
  base_dir,
  "data",
  "metadata"
)


results_dir <- file.path(
  base_dir,
  "results"
)


# ------------------------------------------------------------
# 4. Create results directory
# ------------------------------------------------------------

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Load RMA expression matrix
# ------------------------------------------------------------

expr_matrix <- read.csv(
  file.path(
    processed_dir,
    "GSE25014_expression_RMA.csv"
  ),
  row.names = 1,
  check.names = FALSE
)


# ------------------------------------------------------------
# 6. Inspect expression matrix
# ------------------------------------------------------------

dim(expr_matrix)

head(expr_matrix)

summary(expr_matrix)


# ------------------------------------------------------------
# 7. Load GPL570 annotation
# ------------------------------------------------------------

annotations <- read.csv(
  file.path(
    metadata_dir,
    "GSE25014_probe_annotations.csv"
  ),
  check.names = FALSE
)


# ------------------------------------------------------------
# 8. Inspect annotation
# ------------------------------------------------------------

dim(annotations)

colnames(annotations)

head(annotations)


# ------------------------------------------------------------
# 9. Standardize Gene Symbol column
# ------------------------------------------------------------

annotations <- annotations %>%
  mutate(
    GeneSymbol = trimws(
      as.character(`Gene symbol`)
    )
  )


# ------------------------------------------------------------
# 10. Count genes associated with each probe
# ------------------------------------------------------------

annotations <- annotations %>%
  mutate(
    Number_of_genes = case_when(
      
      is.na(GeneSymbol) |
        GeneSymbol == "" ~ 0L,
      
      grepl(
        "///",
        GeneSymbol,
        fixed = TRUE
      ) ~ lengths(
        strsplit(
          GeneSymbol,
          split = "///",
          fixed = TRUE
        )
      ),
      
      TRUE ~ 1L
    )
  )


# ------------------------------------------------------------
# 11. Check annotation categories
# ------------------------------------------------------------

table(
  annotations$Number_of_genes
)


# ------------------------------------------------------------
# 12. Identify probes without gene annotation
# ------------------------------------------------------------

probes_without_gene <- annotations %>%
  filter(
    Number_of_genes == 0
  )


nrow(probes_without_gene)


# ------------------------------------------------------------
# 13. Identify probes mapping to multiple genes
# ------------------------------------------------------------

ambiguous_probes <- annotations %>%
  filter(
    Number_of_genes > 1
  )


nrow(ambiguous_probes)


# ------------------------------------------------------------
# 14. Keep only probes mapping to exactly one gene
# ------------------------------------------------------------

annotations_clean <- annotations %>%
  filter(
    Number_of_genes == 1
  )


nrow(annotations_clean)


# ------------------------------------------------------------
# 15. Check whether probe IDs are present in expression matrix
# ------------------------------------------------------------

sum(
  annotations_clean$ProbeID %in%
    rownames(expr_matrix)
)


# ------------------------------------------------------------
# 16. Keep only probes present in both datasets
# ------------------------------------------------------------

annotations_clean <- annotations_clean %>%
  filter(
    ProbeID %in% rownames(expr_matrix)
  )


# ------------------------------------------------------------
# 17. Match annotation order to expression matrix
# ------------------------------------------------------------

annotations_clean <- annotations_clean[
  match(
    rownames(expr_matrix),
    annotations_clean$ProbeID
  ),
]


# Remove probes without annotation after matching

annotations_clean <- annotations_clean[
  !is.na(annotations_clean$ProbeID),
]


# ------------------------------------------------------------
# 18. Check dimensions after filtering
# ------------------------------------------------------------

dim(annotations_clean)


# ------------------------------------------------------------
# 19. Extract expression values for valid probes
# ------------------------------------------------------------

expr_clean <- expr_matrix[
  annotations_clean$ProbeID,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 20. Confirm row alignment
# ------------------------------------------------------------

all(
  rownames(expr_clean) ==
    annotations_clean$ProbeID
)


# ------------------------------------------------------------
# 21. Add Gene Symbol to expression data
# ------------------------------------------------------------

expr_clean_df <- as.data.frame(
  expr_clean
)

expr_clean_df$GeneSymbol <-
  annotations_clean$GeneSymbol


# ------------------------------------------------------------
# 22. Check how many probes represent each gene
# ------------------------------------------------------------

probe_count_per_gene <- expr_clean_df %>%
  count(
    GeneSymbol,
    name = "Number_of_probes"
  ) %>%
  arrange(
    desc(Number_of_probes)
  )


head(
  probe_count_per_gene,
  30
)


# ------------------------------------------------------------
# 23. Count genes represented by multiple probes
# ------------------------------------------------------------

genes_multiple_probes <- probe_count_per_gene %>%
  filter(
    Number_of_probes > 1
  )


nrow(genes_multiple_probes)


# ------------------------------------------------------------
# 24. Prepare expression matrix for summarization
# ------------------------------------------------------------

expr_values <- expr_clean_df %>%
  select(
    -GeneSymbol
  )


# ------------------------------------------------------------
# 25. Summarize multiple probes per gene
# ------------------------------------------------------------
#
# For genes represented by multiple probes,
# we calculate the median expression value
# across their probes.
#
# ------------------------------------------------------------

gene_expression <- expr_clean_df %>%
  group_by(
    GeneSymbol
  ) %>%
  summarise(
    across(
      everything(),
      median,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 26. Convert Gene Symbol into row names
# ------------------------------------------------------------

gene_expression <- as.data.frame(
  gene_expression
)

rownames(gene_expression) <-
  gene_expression$GeneSymbol

gene_expression$GeneSymbol <- NULL


# ------------------------------------------------------------
# 27. Check final gene-level matrix
# ------------------------------------------------------------

dim(gene_expression)

head(gene_expression)

summary(gene_expression)


# ------------------------------------------------------------
# 28. Check for duplicated gene names
# ------------------------------------------------------------

sum(
  duplicated(
    rownames(gene_expression)
  )
)


# ------------------------------------------------------------
# 29. Investigate genes of interest
# ------------------------------------------------------------

genes_of_interest <- c(
  "CDH5",
  "MAPK1",
  "HMOX1",
  "NQO1",
  "NFE2L2"
)


genes_present <- intersect(
  genes_of_interest,
  rownames(gene_expression)
)


genes_present


gene_expression[
  genes_present,
  ,
  drop = FALSE
]


# ------------------------------------------------------------
# 30. Save cleaned annotation
# ------------------------------------------------------------

write.csv(
  annotations_clean,
  file = file.path(
    results_dir,
    "GSE25014_clean_probe_annotation.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 31. Save ambiguous probes
# ------------------------------------------------------------

write.csv(
  ambiguous_probes,
  file = file.path(
    results_dir,
    "GSE25014_ambiguous_probes.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 32. Save probe counts per gene
# ------------------------------------------------------------

write.csv(
  probe_count_per_gene,
  file = file.path(
    results_dir,
    "GSE25014_probe_count_per_gene.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 33. Save final gene-level expression matrix
# ------------------------------------------------------------

write.csv(
  gene_expression,
  file = file.path(
    processed_dir,
    "GSE25014_expression_gene_level.csv"
  ),
  row.names = TRUE
)


# ------------------------------------------------------------
# 34. Save final matrix as RDS
# ------------------------------------------------------------

saveRDS(
  gene_expression,
  file = file.path(
    processed_dir,
    "GSE25014_expression_gene_level.rds"
  )
)


# ------------------------------------------------------------
# 35. Final statistics
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("GSE25014 - Lesson 04\n")
cat("============================================\n")

cat(
  "Original probe sets:",
  nrow(expr_matrix),
  "\n"
)

cat(
  "Probes without gene annotation:",
  nrow(probes_without_gene),
  "\n"
)

cat(
  "Probes mapping to multiple genes:",
  nrow(ambiguous_probes),
  "\n"
)

cat(
  "Probes retained:",
  nrow(expr_clean),
  "\n"
)

cat(
  "Final number of genes:",
  nrow(gene_expression),
  "\n"
)

cat(
  "Number of samples:",
  ncol(gene_expression),
  "\n"
)

cat(
  "Duplicated genes:",
  sum(
    duplicated(
      rownames(gene_expression)
    )
  ),
  "\n"
)


# ------------------------------------------------------------
# 36. List output files
# ------------------------------------------------------------

print(
  list.files(
    processed_dir
  )
)

print(
  list.files(
    results_dir
  )
)


# ------------------------------------------------------------
# 37. Finish
# ------------------------------------------------------------

print("Lesson 04 complete.")

print("Probe-to-gene processing was successfully performed.")

print("The final gene-level expression matrix is ready for downstream analysis.")