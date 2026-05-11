############################
# 0) Install packages
############################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

packages <- c(
  "airway",
  "DESeq2",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "SummarizedExperiment",
  "gplots",
  "RColorBrewer",
  "ggplot2"
)

missing_packages <- packages[
  !vapply(packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  BiocManager::install(missing_packages, ask = FALSE, update = FALSE)
}


############################
# 2) Load packages
############################

library(airway)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(SummarizedExperiment)
library(gplots)
library(RColorBrewer)
library(ggplot2)

############################
# 3) Load the airway dataset
############################

data("airway")
se <- airway

############################
# 4) Inspect the data
############################

# Count matrix
head(assay(se))

# Total counts per sample
colSums(assay(se))

# Sample metadata
colData(se)

# Genomic coordinates
rowRanges(se)

############################
# 5) Convert to DESeqDataSet
############################

# Design:
# cell controls for cell-line differences
# dex tests the treatment effect

dds <- DESeqDataSet(se, design = ~ cell + dex)

# Make untreated the reference condition
dds$dex <- relevel(dds$dex, ref = "untrt")

# Estimate size factors for normalized counts
dds <- estimateSizeFactors(dds)

############################
# 6) rlog transformation for visualization
############################

# rlog stabilizes variance for visualization tasks such as PCA and clustering.
# Use raw counts with DESeq() for actual differential expression testing.

rld <- rlog(dds, blind = FALSE)

# Inspect transformed expression values
head(assay(rld))

############################
# 7) Compare ordinary log transform vs rlog
############################

png(
  filename = "02_output_files/figures/normalized_vs_rlog_counts.png",
  width = 1200,
  height = 600,
  res = 150
)

# Save current plotting layout
opar <- par(mfrow = c(1, 2))

# Left plot:
# log2(normalized counts + 1)
plot(
  log2(1 + counts(dds, normalized = TRUE)[, 1:2]),
  col = rgb(0, 0, 0, 0.2),
  pch = 16,
  cex = 0.3,
  main = "log2(normalized counts + 1)",
  xlab = "Sample 1",
  ylab = "Sample 2"
)

# Right plot:
# rlog-transformed values
plot(
  assay(rld)[, 1:2],
  col = rgb(0, 0, 0, 0.2),
  pch = 16,
  cex = 0.3,
  main = "rlog-transformed counts",
  xlab = "Sample 1",
  ylab = "Sample 2"
)

# Restore plotting layout
par(opar)

# Save PNG
dev.off()

############################
# 8) Sample-to-sample distances
############################

# dist() computes distances between rows.
# Because samples are columns in assay(rld), transpose first.

sampleDists <- dist(t(assay(rld)))
sampleDists

# Convert distance object to matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Rename rows for easier interpretation
rownames(sampleDistMatrix) <- paste(rld$dex, rld$cell, sep = "-")
colnames(sampleDistMatrix) <- paste(rld$dex, rld$cell, sep = "-")

# Build blue color palette
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# Hierarchical clustering
hc <- hclust(sampleDists)

# Save heatmap
png(
  filename = "02_output_files/figures/sample_distance_heatmap.png",
  width = 1000,
  height = 900,
  res = 150
)

heatmap.2(
  sampleDistMatrix,
  Rowv = as.dendrogram(hc),
  Colv = as.dendrogram(hc),
  symm = TRUE,
  trace = "none",
  col = colors,
  margins = c(10, 10),
  main = "Sample-to-sample distances"
)

dev.off()

############################
# 9) PCA plot
############################

# Built-in PCA plot
png(
  filename = "02_output_files/figures/pca_builtin.png",
  width = 900,
  height = 700,
  res = 150
)

plotPCA(rld, intgroup = c("dex", "cell"))

dev.off()

# Extract PCA data for customized ggplot
pca_data <- plotPCA(
  rld,
  intgroup = c("dex", "cell"),
  returnData = TRUE
)

# Percentage of variance explained by PC1 and PC2
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Customized PCA plot
p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = dex, shape = cell)) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of rlog-transformed counts") +
  theme_bw()

# Display plot
print(p_pca)

# Save plot
ggsave(
  filename = "02_output_files/figures/pca_customized.png",
  plot = p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

############################
# 10) Run DESeq2 differential expression
############################

# DESeq() performs:
# size factor estimation
# dispersion estimation
# model fitting
# hypothesis testing

dds <- DESeq(dds)

############################
# 11) Extract results
############################

# Treatment effect:
# trt compared with untrt

res <- results(dds, contrast = c("dex", "trt", "untrt"))

# Inspect result table
head(res)

# Metadata describing result columns
mcols(res, use.names = TRUE)

# Summary of significant genes and filtering
summary(res)

############################
# 12) Understand result columns
############################

# baseMean:
# mean normalized count across all samples
#
# log2FoldChange:
# estimated effect size on log2 scale
# positive = higher in treated
# negative = lower in treated
#
# lfcSE:
# standard error of log2 fold change
#
# stat:
# Wald test statistic
#
# pvalue:
# raw p-value
#
# padj:
# multiple-testing adjusted p-value using BH/FDR

############################
# 13) Other comparisons
############################

# Example: compare one cell line against another
# numerator = N061011
# denominator = N61311

res_cell <- results(dds, contrast = c("cell", "N061011", "N61311"))
head(res_cell)

############################
# 14) Multiple testing summary
############################

# How many genes have raw p < 0.05?
sum(res$pvalue < 0.05, na.rm = TRUE)

# How many genes have non-missing p-values?
sum(!is.na(res$pvalue))

# How many genes pass FDR < 0.1?
sum(res$padj < 0.1, na.rm = TRUE)

# Subset significant genes
resSig <- subset(res, padj < 0.1)

# Strongest down-regulated significant genes
head(resSig[order(resSig$log2FoldChange), ])

# Strongest up-regulated significant genes
head(resSig[order(resSig$log2FoldChange, decreasing = TRUE), ])

############################
# 15) Plot counts for the top gene
############################

# Find the gene with the smallest adjusted p-value
topGene <- rownames(res)[which.min(ifelse(is.na(res$padj), Inf, res$padj))]

# Get plotting data from plotCounts()
count_data <- plotCounts(
  dds,
  gene = topGene,
  intgroup = "dex",
  returnData = TRUE
)

# Plot counts on log scale
p_counts <- ggplot(count_data, aes(x = dex, y = count, fill = dex)) +
  scale_y_log10() +
  geom_dotplot(binaxis = "y", stackdir = "center") +
  ggtitle(paste("Counts for top gene:", topGene)) +
  xlab("Treatment condition") +
  ylab("Normalized count") +
  theme_bw()

# Display plot
print(p_counts)

# Save plot
ggsave(
  filename = "02_output_files/figures/top_gene_counts.png",
  plot = p_counts,
  width = 6,
  height = 5,
  dpi = 300
)

############################
# 16) MA plot
############################

# MA plot:
# x-axis = average expression
# y-axis = log2 fold change
# red points are significant genes by default threshold

png(
  filename = "02_output_files/figures/ma_plot.png",
  width = 900,
  height = 700,
  res = 150
)

plotMA(
  res,
  ylim = c(-5, 5),
  main = "MA plot: treated vs untreated"
)

points(
  res[topGene, ]$baseMean,
  res[topGene, ]$log2FoldChange,
  col = "dodgerblue",
  cex = 2,
  lwd = 2
)

text(
  res[topGene, ]$baseMean,
  res[topGene, ]$log2FoldChange,
  labels = topGene,
  pos = 2,
  col = "dodgerblue"
)

dev.off()

############################
# 17) Dispersion plot
############################

png(
  filename = "02_output_files/figures/dispersion_estimates.png",
  width = 900,
  height = 700,
  res = 150
)

plotDispEsts(dds)

dev.off()

############################
# 18) Add gene annotation
############################

# The row names are Ensembl IDs.
# Add HGNC symbols and Entrez IDs for easier interpretation.

res$hgnc_symbol <- unname(
  mapIds(
    org.Hs.eg.db,
    keys = rownames(res),
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
)

res$entrezgene <- unname(
  mapIds(
    org.Hs.eg.db,
    keys = rownames(res),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
)

# Sort by p-value
resOrdered <- res[order(res$pvalue), ]

# Show top results
head(resOrdered)

############################
# 19) Export results
############################

# Export all ordered DESeq2 results
write.csv(
  as.data.frame(resOrdered),
  file = "02_output_files/results/results_ordered_by_pvalue.csv",
  row.names = TRUE
)

# Export significant genes only
write.csv(
  as.data.frame(resSig),
  file = "02_output_files/results/significant_genes_FDR_0.1.csv",
  row.names = TRUE
)

# Export cell-line comparison results
write.csv(
  as.data.frame(res_cell),
  file = "02_output_files/results/cell_comparison_N061011_vs_N61311.csv",
  row.names = TRUE
)

############################
# 20) Save session information
############################

sink("02_output_files/results/sessionInfo.txt")
sessionInfo()
sink()

############################
# 21) Final message
############################

message("Analysis complete.")
message("Figures saved in: 02_output_files/figures")
message("Results saved in: 02_output_files/results")

