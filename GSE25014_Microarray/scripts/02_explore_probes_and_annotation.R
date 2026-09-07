# ============================================================
# Transcriptomics Case Studies
# Case Study: GSE25014
#
# Lesson 02: Exploring Probes and Annotation
# ============================================================
#
# Dataset:
#   GSE25014 - Gene expression data of endothelium exposed to heme
#
# Platform:
#   GPL570 - Affymetrix Human Genome U133 Plus 2.0 Array
#
# Objective:
#   Explore the relationship between microarray probes and
#   gene annotations.
#
# In this lesson we will:
#
#   1. Explore the GPL570 annotation
#   2. Identify probes without gene annotation
#   3. Identify probes mapping to multiple genes
#   4. Identify genes represented by multiple probes
#   5. Investigate genes of interest
#
# Important:
#
#   No probes will be removed permanently in this lesson.
#   Probe filtering and probe-to-gene processing will be
#   performed after RMA in a later lesson.
#
# ============================================================


# ============================================================
# 1. Load packages
# ============================================================

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}

library(dplyr)
library(readr)


# ============================================================
# 2. Define directories
# ============================================================

case_study_dir <- "GSE25014_Microarray"

if (dir.exists(case_study_dir)) {
  base_dir <- case_study_dir
} else {
  base_dir <- "."
}

metadata_dir <- file.path(
  base_dir,
  "data",
  "metadata"
)

results_dir <- file.path(
  base_dir,
  "results"
)

dir.create(
  results_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Import probe annotation
# ============================================================

annotations <- read.csv(
  file.path(
    metadata_dir,
    "GSE25014_probe_annotations.csv"
  ),
  check.names = FALSE
)


# ============================================================
# 4. Explore the annotation table
# ============================================================

dim(annotations)

colnames(annotations)

head(annotations)


# ============================================================
# 5. Explore the main annotation fields
# ============================================================

annotations %>%
  select(
    ProbeID,
    `Gene symbol`,
    `Gene title`,
    ID,
    `Gene ID`
  ) %>%
  head(20)


# ============================================================
# 6. Check duplicated ProbeIDs
# ============================================================
#
# Each probe should normally have one entry in the annotation
# table used by this analysis.
#
# ============================================================

sum(
  duplicated(
    annotations$ProbeID
  )
)


# ============================================================
# 7. Create a standardized GeneSymbol column
# ============================================================

annotations <- annotations %>%
  mutate(
    GeneSymbol = trimws(
      as.character(`Gene symbol`)
    )
  )


# ============================================================
# 8. Identify probes without gene annotation
# ============================================================

probes_without_gene <- annotations %>%
  filter(
    is.na(GeneSymbol) |
      GeneSymbol == ""
  )


nrow(probes_without_gene)


# ============================================================
# 9. Identify probes with gene annotation
# ============================================================

probes_with_gene <- annotations %>%
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != ""
  )


nrow(probes_with_gene)


# ============================================================
# 10. Summarize annotation completeness
# ============================================================

annotation_summary <- tibble(
  Category = c(
    "With gene annotation",
    "Without gene annotation"
  ),
  
  Probes = c(
    nrow(probes_with_gene),
    nrow(probes_without_gene)
  )
)

annotation_summary


# ============================================================
# 11. Identify probes mapping to multiple genes
# ============================================================
#
# The GPL570 annotation can contain multiple gene symbols
# associated with a single probe.
#
# Example:
#
#   MIR4640///DDR1
#
# The separator used in this annotation is "///".
#
# ============================================================

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


# ============================================================
# 12. Extract ambiguous probes
# ============================================================

ambiguous_probes <- annotations %>%
  filter(
    Number_of_genes > 1
  )


nrow(ambiguous_probes)


# ============================================================
# 13. Inspect ambiguous probes
# ============================================================

ambiguous_probes %>%
  select(
    ProbeID,
    GeneSymbol,
    `Gene title`,
    ID,
    `Gene ID`,
    Number_of_genes
  ) %>%
  head(20)


# ============================================================
# 14. Examples of ambiguous probes
# ============================================================

annotations %>%
  filter(
    Number_of_genes > 1
  ) %>%
  select(
    ProbeID,
    GeneSymbol,
    `Gene ID`
  ) %>%
  head(10)


# ============================================================
# 15. Identify genes represented by multiple probes
# ============================================================
#
# The relationship can also occur in the opposite direction:
#
#   Probe A ─┐
#   Probe B ─┼──> Gene X
#   Probe C ─┘
#
# One gene may therefore be represented by several probes.
#
# ============================================================

probe_count_per_gene <- annotations %>%
  
  filter(
    !is.na(GeneSymbol),
    GeneSymbol != "",
    Number_of_genes == 1
  ) %>%
  
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


# ============================================================
# 16. Genes with more than one probe
# ============================================================

genes_multiple_probes <- probe_count_per_gene %>%
  filter(
    Number_of_probes > 1
  )


nrow(genes_multiple_probes)


# ============================================================
# 17. Distribution of probes per gene
# ============================================================

table(
  probe_count_per_gene$Number_of_probes
)


# ============================================================
# 18. Investigate genes of interest
# ============================================================
#
# We will investigate some genes that will be relevant
# throughout the transcriptomic analysis.
#
# ============================================================

genes_of_interest <- c(
  "CDH5",
  "MAPK1",
  "HMOX1",
  "NQO1",
  "NFE2L2"
)


genes_of_interest_annotation <- annotations %>%
  filter(
    GeneSymbol %in% genes_of_interest
  ) %>%
  select(
    ProbeID,
    GeneSymbol,
    `Gene title`,
    ID,
    `Gene ID`,
    Number_of_genes
  ) %>%
  arrange(
    GeneSymbol
  )


genes_of_interest_annotation


# ============================================================
# 19. Count probes for genes of interest
# ============================================================

genes_of_interest_annotation %>%
  count(
    GeneSymbol,
    name = "Number_of_probes"
  )


# ============================================================
# 20. Investigate CDH5
# ============================================================

annotations %>%
  filter(
    GeneSymbol == "CDH5"
  ) %>%
  select(
    ProbeID,
    GeneSymbol,
    `Gene title`,
    `Gene ID`
  )


# ============================================================
# 21. Investigate MAPK1
# ============================================================

annotations %>%
  filter(
    GeneSymbol == "MAPK1"
  ) %>%
  select(
    ProbeID,
    GeneSymbol,
    `Gene title`,
    `Gene ID`
  )


# ============================================================
# 22. Create a probe-level annotation summary
# ============================================================

probe_summary <- annotations %>%
  mutate(
    Annotation_status = case_when(
      
      Number_of_genes == 0 ~
        "No gene annotation",
      
      Number_of_genes == 1 ~
        "One gene",
      
      Number_of_genes > 1 ~
        "Multiple genes"
    )
  ) %>%
  count(
    Annotation_status,
    name = "Number_of_probes"
  )


probe_summary


# ============================================================
# 23. Save annotation summary
# ============================================================

write.csv(
  probe_summary,
  file = file.path(
    results_dir,
    "probe_annotation_summary.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 24. Save ambiguous probes
# ============================================================

write.csv(
  ambiguous_probes,
  file = file.path(
    results_dir,
    "ambiguous_probes.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 25. Save probe count per gene
# ============================================================

write.csv(
  probe_count_per_gene,
  file = file.path(
    results_dir,
    "probe_count_per_gene.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 26. Save genes of interest annotation
# ============================================================

write.csv(
  genes_of_interest_annotation,
  file = file.path(
    results_dir,
    "genes_of_interest_annotation.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 27. Final exploration
# ============================================================

dim(annotations)

nrow(probes_without_gene)

nrow(ambiguous_probes)

nrow(genes_multiple_probes)

n_distinct(
  annotations$GeneSymbol[
    annotations$Number_of_genes == 1
  ]
)


# ============================================================
# 28. Lesson complete
# ============================================================

print("Lesson 02 complete.")

print("Probe annotation has been explored.")

print("Probe filtering will be performed after RMA.")
