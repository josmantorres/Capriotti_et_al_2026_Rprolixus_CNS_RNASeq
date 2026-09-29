# =========================================================
# Differential expression analysis: Fed vs Unfed
# =========================================================
# Tools: R, DESeq2, apeglm
#
# Contrast: Fed vs Unfed (reference level = Unfed)
#   log2FC > 0  -> higher expression in Fed
#   log2FC < 0  -> higher expression in Unfed
#
# Differential expression criterion:
#   apeglm s-value < 0.05, evaluated against |log2FC| > 1
#   (lfcThreshold = 1)
#
# Settings: independentFiltering = TRUE, cooksCutoff = TRUE
# No gene-filtering step is applied before the analysis; genes with
# zero counts across all samples are excluded automatically by DESeq2.
#
# Input files (tab-delimited, placed in `input_dir`):
#   - counts_table.txt : first column = gene IDs; remaining columns =
#                        raw (unnormalized) integer counts, one per sample
#   - sample_metadata.txt : one row per sample, in the same order as the
#                        count columns; must contain a column named
#                        "treatment" with values "Fed" and "Unfed"
#
# Outputs (written to `output_dir`):
#   - DESeq2_results_all_genes.txt
#   - DEGs_Fed_vs_Unfed.txt
#   - MAplot_Fed_vs_Unfed.tif
#   - rlog_values.txt
#   - Sample_distance_heatmap.tif
#   - PCA_rlog.tif
#   - sessionInfo.txt
# =========================================================

library(DESeq2)
library(apeglm)
library(pheatmap)
library(RColorBrewer)

# ---------------------------------------------------------
# 0. Parameters
# ---------------------------------------------------------

input_dir  <- "input"    # folder with the input files
output_dir <- "output"   # folder for the results
alpha_s    <- 0.05       # s-value threshold (also alpha for results/summary)
lfc_thr    <- 1          # |log2FC| threshold (results and lfcShrink)
dec_out    <- "."        # decimal separator for exported tables

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------
# 1. Import count matrix and sample metadata
# ---------------------------------------------------------

counts <- read.delim(file.path(input_dir, "counts_table.txt"),
                     sep = "\t", header = TRUE)

count_matrix <- data.matrix(counts[, 2:ncol(counts)])
rownames(count_matrix) <- counts[, 1]

sample_info <- read.delim(file.path(input_dir, "sample_metadata.txt"),
                          sep = "\t", header = TRUE)

# Consistency checks
stopifnot(ncol(count_matrix) == nrow(sample_info))
stopifnot("treatment" %in% colnames(sample_info))
# If the metadata contains a sample-name column, also verify the order:
# stopifnot(all(colnames(count_matrix) == sample_info$Sample))

# Reference level (Unfed) defined before creating the DESeq2 object
sample_infotreatment<-factor(sampleinfotreatment,
                                levels = c("Unfed", "Fed"))

# ---------------------------------------------------------
# 2. Build the DESeq2 object and fit the model
# ---------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData   = sample_info,
  design    = ~ treatment
)

dds <- DESeq(dds, test = "Wald")

resultsNames(dds)
# "Intercept" "treatment_Fed_vs_Unfed"

coef_name <- "treatment_Fed_vs_Unfed"

# ---------------------------------------------------------
# 3. Wald test against the |log2FC| > lfc_thr threshold
# ---------------------------------------------------------
# With lfcThreshold > 0, pvalue/padj test H0: |LFC| <= lfc_thr.
# Genes flagged by Cook's distance receive pvalue = NA.
# Independent filtering only sets padj to NA for low-mean genes.

res <- results(
  dds,
  name = coef_name,
  lfcThreshold = lfc_thr,
  altHypothesis = "greaterAbs",
  alpha = alpha_s,
  independentFiltering = TRUE,
  cooksCutoff = TRUE
)

# ---------------------------------------------------------
# 4. Log2 fold-change shrinkage with apeglm (s-values)
# ---------------------------------------------------------
# With lfcThreshold > 0, the s-value is the probability that the sign of
# the effect is incorrect or that its magnitude is below the threshold.

res_shrunk <- lfcShrink(
  dds,
  coef = coef_name,
  res = res,
  type = "apeglm",
  lfcThreshold = lfc_thr,
  svalue = TRUE
)

# apeglm does not apply Cook's cutoff: s-values are set to NA for genes
# that results() flagged as outliers (pvalue = NA) or with all-zero counts.
res_shrunksvalue[is.na(respvalue)] <- NA

# Wald statistics (against the LFC threshold) added to the final table
res_shrunkpvalueWald<-respvalue
res_shrunkpadjWald  <-respadj

summary(res_shrunk, alpha = alpha_s)

# ---------------------------------------------------------
# 5. Export results
# ---------------------------------------------------------

res_table <- as.data.frame(res_shrunk)

write.table(res_table,
            file = file.path(output_dir, "DESeq2_results_all_genes.txt"),
            sep = "\t", dec = dec_out,
            row.names = TRUE, col.names = NA, quote = FALSE)

# Differentially expressed genes: s-value < alpha_s
degs <- res_table[!is.na(res_tablesvalue)&restablesvalue < alpha_s, ]
degsRegulation<-ifelse(degslog2FoldChange > 0, "Up_in_Fed", "Up_in_Unfed")
degs <- degs[order(degs$svalue), ]

table(degs$Regulation)

write.table(degs,
            file = file.path(output_dir, "DEGs_Fed_vs_Unfed.txt"),
            sep = "\t", dec = dec_out,
            row.names = TRUE, col.names = NA, quote = FALSE)

# ---------------------------------------------------------
# 6. MA plot
# ---------------------------------------------------------

tiff(file.path(output_dir, "MAplot_Fed_vs_Unfed.tif"),
     width = 15, height = 15, units = "cm", res = 300,
     pointsize = 10, compression = "lzw")
plotMA(res_shrunk, alpha = alpha_s, ylim = c(-10, 10),
       main = "Fed vs Unfed", cex = 0.5)
abline(h = c(-lfc_thr, lfc_thr), col = "blue")
dev.off()

# ---------------------------------------------------------
# 7. Exploratory analyses: rlog, sample distances, PCA
# ---------------------------------------------------------
# blind = TRUE: transformation independent of the experimental design,
# so that the quality-control plots are unbiased.

rld <- rlogTransformation(dds, blind = TRUE, fitType = "parametric")

write.table(data.frame(assay(rld)),
            file = file.path(output_dir, "rlog_values.txt"),
            sep = "\t", dec = dec_out,
            row.names = TRUE, col.names = NA, quote = FALSE)

sample_dists <- dist(t(assay(rld)))            # Euclidean distance
dist_matrix  <- as.matrix(sample_dists)
rownames(dist_matrix) <- colnames(dds)
colnames(dist_matrix) <- colnames(dds)

heat_colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

tiff(file.path(output_dir, "Sample_distance_heatmap.tif"),
     width = 10, height = 10, units = "cm", res = 300,
     pointsize = 10, compression = "lzw")
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = heat_colors)
dev.off()

# PCA on the 500 most variable genes (plotPCA default: ntop = 500)
tiff(file.path(output_dir, "PCA_rlog.tif"),
     width = 15, height = 15, units = "cm", res = 300,
     pointsize = 10, compression = "lzw")
print(plotPCA(rld, intgroup = "treatment"))
dev.off()

# ---------------------------------------------------------
# 8. Record software versions for reproducibility
# ---------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(output_dir, "sessionInfo.txt"))


