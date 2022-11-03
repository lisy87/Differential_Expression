## DIFFERENTIAL GENE EXPRESSION ANALYSIS USING DESeq2
# http://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#quick-start
# http://master.bioconductor.org/packages/release/bioc/html/DESeq2.html

# Pipeline based on:
# https://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf


# Installing DESeq2

if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("DESeq2")

library(DESeq2)
library (magrittr)

# GETTING THE DATA!!

setwd("PATH")

# Defining tumor name
tumor <- "name.of.the.tumor_" 

# Get the table of read counts obtained with featureCounts
readcounts <- read.table(paste0(tumor,"gene_count.txt"), header = TRUE)
row.names(readcounts) <- readcounts$Geneid
readcounts <- readcounts[ , -c (1:6)]

# To assign the sample names
orig_names <- names(readcounts)
names(readcounts) <-c("sample1", "sample2", "sample3", ...) # Please, obtain the names from the gene_count.txt to kept the correct order

# Check the data
str(readcounts)
head(readcounts, n=3)

# Generate the sample condition file
# This file will determine the comparison levels. So, we have to change it if we want to change the comparison.
# We will generate the most general file with our samples using the sample_names.txt obtained previously.
sample_info <- read.delim(paste0(tumor,"sample_condition.txt"))
sample_info

# Generate the DESeqDataSet
DESeq.ds <- DESeqDataSetFromMatrix(countData = readcounts , colData = sample_info , design = ~ condition)

# Check
colData(DESeq.ds) %>% head
assay(DESeq.ds, "counts") %>% head
rowData(DESeq.ds) %>% head

# Test counts
counts(DESeq.ds) %>% str

# Remove genes without any counts
DESeq.ds <- DESeq.ds[rowSums(counts(DESeq.ds)) > 0, ]
# Test counts again
counts(DESeq.ds) %>% str

# NORMALIZING READS COUNTS!!

# Calculate the size factor and add it to the data set
DESeq.ds <- estimateSizeFactors(DESeq.ds)
sizeFactors(DESeq.ds)

# Check colData and corroborate that the sizeFactors have been add to it
colData(DESeq.ds)

# Retrieve the _normalized_ read counts
counts.sf_normalized <- counts(DESeq.ds, normalized = TRUE )

# Log2 transformation of read counts
log.norm.counts <- log2(counts.sf_normalized + 1)

#Saving normalized counts list
write.table(log.norm.counts, file = "log.norm.counts.txt")
write.table(log.norm.counts, file = paste0(tumor, 'log.norm.counts.txt'))

# Graph untransformed and log2-transformed read counts
pdf(paste0(tumor,"boxplot_raw_VS_log2_readcounts.pdf"), width = 8, height = 11)
par(mfrow=c(2,1))
par(mar=c(4,3,3,1), las = 2)
boxplot(counts.sf_normalized, notch = TRUE, main = "untransformed read counts", ylab = "read counts", cex.axis = 0.6)
boxplot(log.norm.counts, notch = TRUE, main = " log2 transformed read counts ", ylab = "log2(read counts)", cex.axis = 0.6)
dev.off()

# To graph only log2-transformed read counts
pdf(paste0(tumor,"boxplot_log2transf_readcounts.pdf"), width = 8, height = 11)
par(mfrow=c(1,1))
par(mar=c(6,3,3,1), las = 2)
boxplot(log.norm.counts, notch = TRUE, main = " log2 transformed read counts ", ylab = "log2(read counts)", cex.axis = 0.6)
dev.off()

# Exploring normalized read counts
par(mfrow=c(1,1))
plot(log.norm.counts[ ,1:2], cex =0.1, main = "Normalized log2 ( read counts )")

# Check for heteroskedasticity plotting the mean vs. the standard deviation:
BiocManager::install(c("vsn"), lib = "/path/")
install.packages("hexbin")
library(vsn)
library(ggplot2)
library(hexbin)

# mean-sd plot
msd_plot1 <- meanSdPlot(log.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot1$gg + ggtitle ("sequencing depth normalized log2 (read counts)") + ylab ("standard deviation")

# Obtaining regularized log - transformed values
DESeq.rlog <- rlog(DESeq.ds, blind = TRUE)
rlog.norm.counts <- assay(DESeq.rlog)

# mean-sd plot for rlog-transformed data
library(vsn)
library(ggplot2)
msd_plot2 <- meanSdPlot(rlog.norm.counts, ranks = FALSE, plot = FALSE)
msd_plot2$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation")

# Plotting the graphs

pdf(paste0(tumor,"pairwise_comparisons_replicate_samples.pdf"), width = 8, height = 11)
par(mfrow=c(2,1))
par(mar=c(3,2,3,1))
plot(log.norm.counts[ ,1:2], cex =0.1, main = " Normalized log2 ( read counts )")
plot(rlog.norm.counts[ ,1:2], cex =0.1, main = " Normalized r-log2 ( read counts )")
dev.off()

pdf(paste0(tumor,"mean-sd-plots.pdf"), width = 7, height = 7)
par(mfrow=c(2,1))
msd_plot1$gg + ggtitle ("sequencing depth normalized log2 (read counts)") + ylab ("standard deviation")
msd_plot2$gg + ggtitle("rlog-transformed read counts") + ylab("standard deviation")
dev.off()

# Exploring global read count patterns

# Calculate the correlation between columns of a matrix
distance.m_rlog <-as.dist(1 - cor(rlog.norm.counts, method = "pearson"))

# Plotting Dendrogram
pdf(paste0(tumor,"Dendrogram.pdf"), width = 7, height = 7)
plot(hclust(distance.m_rlog), labels = colnames(rlog.norm.counts), main = "rlog-transformed read counts\ndistance : Pearson correlation ", cex.axis =0.6, cex.lab= 0.8, cex=0.7, cex.main=0.8)
dev.off()

# Principal Components Analysis (PCA)

library (ggplot2)

P <- plotPCA(DESeq.rlog)
P <- P + theme_bw() + ggtitle("Rlog transformed counts")
P

# Add names to the points
nudge <- position_nudge(y = 1)
Q <- P + theme_bw() + ggtitle("Rlog transformed counts") + geom_text(aes(label = name), position = nudge, size = 2, color ="black")
Q

#Saving plot
pdf(paste0(tumor,"PCA.pdf"), width = 7, height = 7)
Q
dev.off()

# DIFFERENTIAL GENE EXPRESSION ANALYSIS!!

# We would like to compare the resistant lines versus parental lines samples
# with the parental values used as the denominator for the fold change calculation.

# DESeq2 uses the levels of the condition to determine the order of the comparison
str(colData(DESeq.ds)$condition)
colData(DESeq.ds)$condition <- relevel(colData(DESeq.ds)$condition, "parental")

# running the DGE analysis
DESeq.ds <- DESeq(DESeq.ds)
DGE.results <- results(DESeq.ds, independentFiltering = TRUE , alpha = 0.05)
summary(DGE.results)

#Saving results in files
sink(paste0(tumor,'summary_DGE.results.txt'))
cat ("*************************************\n")
cat ("summary DGE results\n")
cat ("*************************************\n")
summary(DGE.results)
cat ("*************************************\n")
cat ("DGE p < 0.05\n")
cat ("*************************************\n")
table(DGE.results$padj < 0.05)
cat ("*************************************\n")
cat ("All DE genes")
rownames(subset(DGE.results, padj < 0.05))
sink()

upregulated <- rownames(subset(DGE.results, log2FoldChange > 0 & padj < 0.05))
downregulated <- rownames(subset(DGE.results, log2FoldChange < 0 & padj < 0.05))

write(downregulated, file = "downregulated_DGE.txt")
write(upregulated, file = "upregulated_DGE.txt")

# Exploratory plots following DGE analysis

# Histograms
pdf(paste0(tumor,"pvalues_freq.pdf"), width = 7, height = 7)
par(mar=c(5,4,3,1))
hist(DGE.results$pvalue, col = "grey", border = "white", xlab = " ", ylab = " ", main = "frequencies of p-values")
dev.off()

# MA plot 
pdf(paste0(tumor,"MA.pdf"), width = 7, height = 7)
par(mar=c(5,4,3,1), las = 1)
plotMA(DGE.results, alpha = 0.05 , main = " Parental vs Resiatant", ylim = c(-4 ,4))
dev.off()

# Heatmaps
install.packages("NMF")
library(NMF)

DGE.results.sorted <- DGE.results[order(DGE.results$padj) , ]
DGEgenes <- rownames(subset(DGE.results.sorted, padj < 0.05))
hm.mat_DGEgenes <- log.norm.counts[DGEgenes, ]

# Saving results to file
sink(paste0(tumor,'DGEgenes_normreadcounts.txt'))
cat ("normalized read counts for DE genes\n")
cat ("*************************************\n")
hm.mat_DGEgenes
sink()

# Plot the normalized read counts of DE genes sorted by the adjusted p-value
pdf(paste0(tumor,"heatmap_pvalue.pdf"), width = 7, height = 7, onefile = FALSE)
aheatmap(hm.mat_DGEgenes, Rowv = NA , Colv = NA)
dev.off()

# Combine the heatmap with hierarchical clustering
pdf(paste0(tumor,"heatmap_clustering.pdf"), width = 7, height = 7, onefile = FALSE)
aheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average")
dev.off()

# Scale the read counts per gene to emphasize the sample-type-specific differences
pdf(paste0(tumor,"heatmap_zscore.pdf"), width = 7, height = 7, onefile = FALSE)
aheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average", scale = "row")
dev.off()


# Ordered heatmap
install.packages('pheatmap') #https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap  
library(pheatmap)
pheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average", scale = "row", show_rownames = F, angle_col = 45)

pdf(paste0(tumor,"heatmap_zscore_order.pdf"), width = 7, height = 7, onefile = FALSE)
pheatmap(hm.mat_DGEgenes, Rowv = TRUE, Colv = TRUE, distfun = "euclidean", hclustfun = "average", scale = "row", show_rownames = F, angle_col = 45)
dev.off()

# Volcano plots

par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DGE.results)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)

# Add labels to DGE genes if you want
# with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))

pdf(paste0(tumor,"volcanoplot.pdf"), width = 7, height = 7, onefile = FALSE)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(DGE.results)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)
dev.off()

# EXPLORING THE DIFFERENTIALLY EXPRESSED GENES

# Plot read counts for a single gene
ibrary(grDevices)
#plotCounts(dds = DESeq.ds, gene = "ENSG00000111885.7", normalized = TRUE, transform = FALSE, main = "DGE genes expression")

# Plot read counts for multiple genes
# https://hbctraining.github.io/DGE_workshop/lessons/06_DGE_visualizing_results.html

# GENE ONTOLOGY ANALYSIS

# Installing org.Hs.eg.db library
BiocManager :: install("org.Hs.eg.db", lib = "/path/")
library(org.Hs.eg.db, lib.loc = "/path/")

# list the types of keywords that are available to query the annotation database
keytypes(org.Hs.eg.db)

# list columns that can be retrieved from the annotation data base
columns(org.Hs.eg.db)

# Convert DGEgenes to eliminate version information from the ENSEMBL code
DGEgenes
s <- DGEgenes
DGEgenes2 <- sapply(strsplit(s, split='.', fixed=TRUE), function(x) (x[1]))
DGEgenes2

#Changing the names of the genes in DGE.results for downstream analysis:
DGE.results_names <- rownames(DGE.results)
t <- DGE.results_names
DGE.results_names2 <- sapply(strsplit(t, split='.', fixed=TRUE), function(x) (x[1]))
DGE.results2 <- DGE.results
rownames(DGE.results2) <- DGE.results_names2

# Obtaining down and up regulated genes with new names
upregulated_newnames <- rownames(subset(DGE.results2, log2FoldChange > 0 & padj < 0.05))
downregulated_newnames <- rownames(subset(DGE.results2, log2FoldChange < 0 & padj < 0.05))

write(downregulated_newnames, file = paste0(tumor,"downregulated_DGE_newnames.txt"))
write(upregulated_newnames, file =  paste0(tumor,"upregulated_DGE_newnames.txt"))

# Make a batch retrieval for all DE genes
anno <- select(org.Hs.eg.db, keys = DGEgenes2, keytype = "ENSEMBL", columns = c("ENSEMBL","GENENAME","CHR"))
write.csv(anno, file = paste0(tumor,'DGEgenes_ENSEMBL.csv'))

# Make a batch retrieval for up DE genes
anno_up <- select(org.Hs.eg.db, keys = upregulated_newnames, keytype = "ENSEMBL", columns = c("ENSEMBL","GENENAME","CHR"))
write.csv(anno_up, file = paste0(tumor,'upregulated_DEgenes_ENSEMBL.csv'))

# Make a batch retrieval for down DE genes
anno_down <- select(org.Hs.eg.db, keys = downregulated_newnames, keytype = "ENSEMBL", columns = c("ENSEMBL","GENENAME","CHR"))
write.csv(anno_down, file = paste0(tumor,'downregulated_DEgenes_ENSEMBL.csv'))

# Merge the information of the DGE analysis with the information about the genes
out.df <- merge (as.data.frame(DGE.results2), anno, by.x = "row.names", by.y = "ENSEMBL")
# Saving info to file
write.csv(out.df, file = paste0(tumor,'DGEgenes_ENSEMBL_all_info.csv'))

# GO TERM ENRICHMENT ANALYSIS OF DEG 
# Based on https://github.com/friedue/course_RNA-seq2019/blob/master/Day04/ReadCountsToDGE_03_GOtermEnrichment.Rmd
BiocManager::install("goseq", lib = "/path/") # package for GO term enrichment

library(goseq) 
library(DESeq2)
library(magrittr)
library(ggplot2)

gene.vector2 <- row.names(DGE.results2) %in% DGEgenes2 %>% as.integer
names(gene.vector2) <- row.names(DGE.results2)

# Quantifying the length bias (= weight for each gene)
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", lib = "/path/")
BiocManager::install("GenomicFeatures", lib = "/path/")
library(GenomicFeatures, lib = "/path/")
pwf <- nullp(gene.vector2, "hg38", "ensGene")
head(pwf)

# Test for enrichment of GO terms
GO.wall <- goseq(pwf, "hg38", "ensGene", use_genes_without_cat=TRUE)
head(GO.wall)

# Retrieving the GO categories assigned to each gene
go_gns <- getgo( rownames(DGE.results2), 'hg38', 'ensGene') %>% stack
write.csv(go_gns, file = paste0(tumor,'GO_genes.csv'))

# Adding gene information to the goseq results:
merge(GO.wall, go_gns, by.x = "category", by.y = "values") %>% dim

# Using REVIGO to summarize results:

# Generating REVIGO's input
# Export a list of over-represented GO terms and the corresponding p-value
# which is the input that REVIGO needs
subset(GO.wall, over_represented_pvalue < 0.01, 
       select = c("category","over_represented_pvalue")) %>%
  write.table(.,
              file = paste0(tumor, "Enriched_GOterms_goseq.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

# You can use the output of previous command in the REVIGO web: http://revigo.irb.hr/
# REVIGO generates multiple plots; the easiest way to obtain them in high quality is to download the R scripts that it offers and to run those yourself.


## FILTERING DEG!!

# Selecting the threshold

#In a loop testing different thresholds:
for (i in seq(0,4, by = 0.1)) {
  object <- paste("DGE.results.filt_",i, sep="")
  print(object)
  x <- results(DESeq.ds, independentFiltering = TRUE , alpha = 0.05 , lfcThreshold = i, altHypothesis="greaterAbs")
  summary(x)
}

#Saving results to file
sink(paste0(tumor,'test_lfcThreshold.txt'))
for (i in seq(0,4, by = 0.1)) {
  object <- paste("DGE.results.filt_",i, sep="")
  print(object)
  x <- results(DESeq.ds, independentFiltering = TRUE , alpha = 0.05 , lfcThreshold = i, altHypothesis="greaterAbs")
  summary(x)
}
sink() 

# Select the values to do the table in terminal (it's more easy ;))
# printf Regulation'\t'Threshold'\t'Count'\n' > "tumor"_table_threshold.tsv 
# capture down: cat "tumor"_test_lfcThreshold.txt | grep "LFC <" | cut -d ":" -f 2 | cut -d "," -f 1 >> "tumor"_table_threshold.tsv
# capture up: cat "tumor"_test_lfcThreshold.txt | grep "LFC >" | cut -d ":" -f 2 | cut -d "," -f 1 >> "tumor"_table_threshold.tsv
# Create the Threshold column
for (i in seq(0,4, by = 0.1)) {
  print(i)
}
# Create the tsv file editing the file "tumor"_table_threshold.tsv

# Import csv file
table <- read.table(file = paste0(tumor,"table_threshold.tsv"), header = TRUE, sep ="\t", dec = ".")
table
 
# Create the graph
library(ggplot2)
ggplot(table, aes(x = factor(Threshold), y = factor(Count), colour = Regulation, group = Regulation)) + geom_line()

# Saving graph
pdf(paste0(tumor,"Threshold_graph.pdf"), width = 14, height = 10)
ggplot(table, aes(x = factor(Threshold), y = factor(Count), colour = Regulation, group = Regulation)) + geom_line()
dev.off()

# Select the threshold and create the variables with them

thresh_1 <- "selected value 1"
thresh_2 <- "selected value 2"

# Create a new object

DGE.results.filt_thresh_1 <- results(DESeq.ds, independentFiltering = TRUE , alpha = 0.05 , lfcThreshold = thresh_1, altHypothesis="greaterAbs")
summary(DGE.results.filt_thresh_1)

# And run all the protocol using the new object DGE.results.filt_thresh_1  

