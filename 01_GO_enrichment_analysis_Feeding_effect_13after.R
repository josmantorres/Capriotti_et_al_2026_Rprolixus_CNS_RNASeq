# ============================================================
# GO enrichment analysis using ermineR
# ============================================================
#
# Description:
# Gene Ontology (GO) enrichment analysis of RNA-seq differential
# expression results from the nervous system of Rhodnius prolixus.
#
# The analysis evaluates transcriptional changes associated with
# feeding status 13 hours after feeding by comparing:
#
#   - Fed: 13 hours after feeding
#   - Unfed: fasting condition
#
# GO enrichment is evaluated separately for:
#   - Biological Process (BP)
#   - Molecular Function (MF)
#   - Cellular Component (CC)
#
# The Gene Score Resampling (GSR) method implemented in ermineR
# is used for the enrichment analysis.
#
# ============================================================


# ------------------------------------------------------------
# 1. Install required packages
# ------------------------------------------------------------

# Install ermineR from GitHub if necessary:
# install.packages("remotes")
# remotes::install_github("PavlidisLab/ermineR")

# rJava installation on Ubuntu may require additional
# system-level dependencies. See:
# https://www.r-bloggers.com/2018/02/installing-rjava-on-ubuntu/


# ------------------------------------------------------------
# 2. Load required packages
# ------------------------------------------------------------

library(ermineR)
library(rJava)
library(httr)


# Increase the timeout for downloading annotation resources
# when internet connection speed is slow.
options(timeout = max(10000, getOption("timeout")))


# ------------------------------------------------------------
# 3. Input and output directories
# ------------------------------------------------------------

# Replace these paths with the directories used on your system.

input_dir <- "path/to/input/files"
output_dir <- "path/to/output/files"


# ------------------------------------------------------------
# 4. Load GO annotation data
# ------------------------------------------------------------

# The GO annotation file must contain:
#   - a gene identifier column ("names")
#   - a GO term column ("Goterm")
#
# Each row represents a gene-GO term association.

goAnnots <- read.table(
  file.path(input_dir, "GO_terms.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Prepare GO annotation list for ermineR
# ------------------------------------------------------------

# Create a list of GO terms associated with each gene.
GOList <- split(
  goAnnots$Goterm,
  goAnnots$names
)

# Convert the GO annotation list to the format required by ermineR.
annotationList <- makeAnnotation(
  GOList,
  return = TRUE
)


# ------------------------------------------------------------
# 6. Load differential expression scores
# ------------------------------------------------------------

# The score file contains gene-level p-values obtained from
# the differential expression analysis comparing Fed vs Unfed
# animals 13 hours after feeding.
#
# Because ermineR is configured with logTrans = TRUE and
# bigIsBetter = FALSE, smaller p-values correspond to stronger
# evidence for differential expression.

scoreList_Fed_vs_Unfed <- read.table(
  file.path(input_dir, "IDs_vs_pvalue_Fed_vs_Unfed.txt"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 7. GO enrichment analysis
# ------------------------------------------------------------
#
# ermineR parameters:
#
#   aspects:
#       "B" = Biological Process
#       "M" = Molecular Function
#       "C" = Cellular Component
#
#   test = "GSR":
#       Gene Score Resampling
#
#   pAdjust = "FDR":
#       False Discovery Rate correction
#
#   iterations = 200000:
#       Number of resampling iterations
#
#   geneReplicates = "mean":
#       Mean score used for replicated gene identifiers
#
#   stats = "mean":
#       Mean statistic used by the GSR analysis
#
#   logTrans = TRUE:
#       Log-transform the input scores
#
#   bigIsBetter = FALSE:
#       Smaller p-values correspond to stronger evidence
#       for differential expression.
#
#   minClassSize = 10:
#       Minimum number of genes in a GO category
#
#   maxClassSize = 200:
#       Maximum number of genes in a GO category
#
# ------------------------------------------------------------


# ============================================================
# Biological Process
# ============================================================

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "B",
  scores = scoreList_Fed_vs_Unfed,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_BP_Fed_vs_Unfed.txt"
  )
)


# ============================================================
# Molecular Function
# ============================================================

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "M",
  scores = scoreList_Fed_vs_Unfed,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_MF_Fed_vs_Unfed.txt"
  )
)


# ============================================================
# Cellular Component
# ============================================================

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "C",
  scores = scoreList_Fed_vs_Unfed,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_CC_Fed_vs_Unfed.txt"
  )
)


# ============================================================
# End of analysis
# ============================================================
