# Install packages:
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "airway",         # Example RNA-seq dataset used in the tutorial
  "DESeq2",         # Differential expression analysis
  "org.Hs.eg.db",   # Human gene annotation
  "AnnotationDbi",  # ID mapping utilities such as mapIds()
  "SummarizedExperiment",
  "gplots",         # heatmap.2()
  "RColorBrewer",   # heatmap colors
  "ggplot2"         # PCA plot and count plot
))


# Load data:
library(airway)
data("airway")
se <- airway

head(assay(se))

# Assessing differential expression:
library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex) # transform the data format for analysis
rld <- rlog(dds) # rlog transformation

sampleDists <- dist(t(assay(rld))) # rlog transformation calculate the Euclidean distance between samples treated vs. untreated

dds$dex <- relevel(dds$dex, "untrt") # making sure the default log2 fold changes are calculated as treated over untreated

dds <- DESeq(dds) # run the differential expression pipeline

(res <- results(dds)) # building the result table

summary(res)

resSig <- subset(res, padj < 0.1) # considering all genes with an adjusted p value below 10 as significant
top_DEG <- resSig[ order( resSig$log2FoldChange ), ]

library(org.Hs.eg.db)
top_DEG$hgnc_symbol <- unname(mapIds(org.Hs.eg.db, rownames(top_DEG), "SYMBOL", "ENSEMBL")) # adding gene names
write.csv(top_DEG,
          file = "../02_output_files/results/top_DEG_dex_vs_untreated.csv",
          row.names = top_DEG$hgnc_symbol)


# Visualizations:
library("gplots")
library("RColorBrewer")

## Heatmap:
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$dex, rld$cell, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
png("../02_output_files/figures/heatmap_dex_untreated.png", width = 1200, height = 1200)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )
dev.off()

png("../02_output_files/figures/pca_dex_untreated.png", width = 600, height = 600)
plotPCA(rld, intgroup = c("dex", "cell"))
dev.off()




