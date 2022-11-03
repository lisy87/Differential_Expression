# Differential Gene Expression in resistant tumoral lines
Differential expression analysis of cellular lines resistant to cisplatin in different head and neck squamous cell carcinoma (HNSCC)

# General workflow

The Fastq raw data files were analysed with FastQC using default parameters, and raw reads were filtered with Trimmomatic to remove low-quality reads and Illumina adapters. The filtered reads were mapped against the GRCh38 reference genome (primary assembly and annotation) with STAR. The BAM files were filtered with samtools to retain mapped reads with quality scores higher than 40. Read counts were obtained for all samples from the BAM files using the featureCounts tool with the default options for paired-end reads. We performed differential expression analysis using the DESeq2 and the parental condition as the denominator of the comparison. Significance was set at a false discovery rate (FDR)-adjusted p-value < 0.05 and log2-fold change (log2FC) values of 1. Batch retrieval of differentially expressed genes was performed using the org.Hs.e.g.db package. Gene Ontology enrichment analysis of differentially expressed gene identification was carried out using the goseq library.

Citation:
