# ============================================================
# Transcriptomics Case Studies
# Case Study: GSE25014
#
# Lesson 03: RMA Normalization
# ============================================================
#
# Dataset:
#   GSE25014 - Gene expression data of endothelium exposed to heme
#
# Platform:
#   GPL570 - Affymetrix Human Genome U133 Plus 2.0 Array
#
# Samples:
#   24
#
# Technology:
#   Microarray
#
# Objective:
#   Import the raw CEL files and perform RMA normalization.
#
# RMA includes:
#
#   1. Background correction
#   2. Quantile normalization
#   3. Probe-set summarization
#
# Important:
#
#   The probe-to-gene treatment will NOT be performed
#   in this lesson.
#
#   This will be performed in Lesson 04.
#
# ============================================================


# ------------------------------------------------------------
# 1. Install required packages
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("affy", quietly = TRUE)) {
  BiocManager::install(
    "affy",
    ask = FALSE,
    update = FALSE
  )
}


# ------------------------------------------------------------
# 2. Load package
# ------------------------------------------------------------

library(affy)


# ------------------------------------------------------------
# 3. Define directories
# ------------------------------------------------------------

case_study_dir <- "GSE25014_Microarray"

if (dir.exists(case_study_dir)) {
  base_dir <- case_study_dir
} else {
  base_dir <- "."
}


cel_dir <- file.path(
  base_dir,
  "data",
  "raw",
  "GSE25014",
  "CEL"
)


processed_dir <- file.path(
  base_dir,
  "data",
  "processed"
)


results_dir <- file.path(
  base_dir,
  "results"
)


# ------------------------------------------------------------
# 4. Create output directories
# ------------------------------------------------------------

dir.create(
  processed_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 5. Find CEL files
# ------------------------------------------------------------

cel_files <- list.files(
  cel_dir,
  pattern = "\\.CEL$",
  full.names = TRUE,
  ignore.case = TRUE
)


# Number of CEL files

length(cel_files)


# Show first files

head(cel_files)


# ------------------------------------------------------------
# 6. Check number of samples
# ------------------------------------------------------------

if (length(cel_files) != 24) {
  
  warning(
    paste(
      "Expected 24 CEL files, but found",
      length(cel_files)
    )
  )
  
}


# ------------------------------------------------------------
# 7. Import CEL files
# ------------------------------------------------------------

raw_data <- ReadAffy(
  celfile.path = cel_dir
)


# Show object

raw_data


# ------------------------------------------------------------
# 8. Check sample names
# ------------------------------------------------------------

sampleNames(raw_data)


# ------------------------------------------------------------
# 9. Check dimensions
# ------------------------------------------------------------

dim(raw_data)


# ------------------------------------------------------------
# 10. Extract raw expression matrix
# ------------------------------------------------------------
#
# raw_data is an AffyBatch object.
#
# exprs() extracts the expression matrix from the object.
#
# Rows    = probe sets
# Columns = samples
#
# ------------------------------------------------------------

raw_expression <- exprs(raw_data)


# Check dimensions

dim(raw_expression)


# Summary of raw expression values

summary(raw_expression)


# First rows

head(raw_expression)


# ------------------------------------------------------------
# 11. Boxplot before RMA
# ------------------------------------------------------------
#
# IMPORTANT:
#
# We use exprs(raw_data) because raw_data itself is an
# AffyBatch object.
#
# boxplot() needs the numerical expression matrix.
#
# ------------------------------------------------------------

#pdf(
#  file.path(
#    results_dir,
#    "GSE25014_boxplot_before_RMA.pdf"
#  ),
#  width = 12,
#  height = 6
#)

#boxplot(
#  exprs(raw_data),
#  main = "GSE25014 - Before RMA",
#  las = 2,
#  cex.axis = 0.5
#)

#dev.off()


# ------------------------------------------------------------
# 12. Apply RMA normalization
# ------------------------------------------------------------

rma_data <- rma(
  raw_data
)


# Show RMA object

rma_data


# ------------------------------------------------------------
# 13. Extract RMA expression matrix
# ------------------------------------------------------------
#
# rma_data is an ExpressionSet object.
#
# exprs() extracts the normalized expression matrix.
#
# ------------------------------------------------------------

expr_matrix <- exprs(rma_data)


# Check dimensions

dim(expr_matrix)


# Summary of normalized values

summary(expr_matrix)


# First rows

head(expr_matrix)


# ------------------------------------------------------------
# 14. Check expression scale
# ------------------------------------------------------------

range(expr_matrix)

summary(expr_matrix)


# ------------------------------------------------------------
# 15. Boxplot after RMA
# ------------------------------------------------------------
#
# IMPORTANT:
#
# Use exprs(rma_data), NOT rma_data directly.
#
# rma_data is an ExpressionSet.
#
# exprs(rma_data) is the numerical expression matrix.
#
# ------------------------------------------------------------

#pdf(
#  file.path(
#    results_dir,
#    "GSE25014_boxplot_after_RMA.pdf"
#  ),
#  width = 12,
#  height = 6
#)

#boxplot(
#  exprs(rma_data),
#  main = "GSE25014 - After RMA",
#  las = 2,
#  cex.axis = 0.5
#)

#dev.off()


# ------------------------------------------------------------
# 16. Compare before and after RMA
# ------------------------------------------------------------
#
# Again, we use exprs() to extract the numerical matrices.
#
# ------------------------------------------------------------

#pdf(
#  file.path(
#    results_dir,
#    "GSE25014_boxplot_before_after_RMA.pdf"
#  ),
#  width = 14,
#  height = 6
#)


#par(mfrow = c(1, 2))


# Before RMA

#boxplot(
#  exprs(raw_data),
#  main = "Before RMA",
#  las = 2,
#  cex.axis = 0.5
#)


# After RMA

#boxplot(
#  exprs(rma_data),
#  main = "After RMA",
#  las = 2,
#  cex.axis = 0.5
#)


#par(mfrow = c(1, 1))


#dev.off()


# ------------------------------------------------------------
# 17. Save RMA expression matrix
# ------------------------------------------------------------

write.csv(
  expr_matrix,
  file = file.path(
    processed_dir,
    "GSE25014_expression_RMA.csv"
  )
)


# ------------------------------------------------------------
# 18. Save RMA ExpressionSet
# ------------------------------------------------------------

saveRDS(
  rma_data,
  file = file.path(
    processed_dir,
    "GSE25014_RMA_ExpressionSet.rds"
  )
)


# ------------------------------------------------------------
# 19. Save sample information
# ------------------------------------------------------------

sample_info <- data.frame(
  SampleID = sampleNames(rma_data),
  stringsAsFactors = FALSE
)


write.csv(
  sample_info,
  file = file.path(
    processed_dir,
    "GSE25014_RMA_samples.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------
# 20. Final checks
# ------------------------------------------------------------

cat("\n")
cat("============================================\n")
cat("GSE25014 - Lesson 03\n")
cat("============================================\n")

cat(
  "Number of CEL files:",
  length(cel_files),
  "\n"
)

cat(
  "Number of samples:",
  ncol(expr_matrix),
  "\n"
)

cat(
  "Number of probe sets:",
  nrow(expr_matrix),
  "\n"
)

cat(
  "Minimum expression value:",
  min(expr_matrix),
  "\n"
)

cat(
  "Maximum expression value:",
  max(expr_matrix),
  "\n"
)


# ------------------------------------------------------------
# 21. Check output files
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
# 22. Finish
# ------------------------------------------------------------

print("Lesson 03 complete.")

print("RMA normalization was successfully performed.")

print("The expression matrix is ready for probe-to-gene treatment in Lesson 04.")

