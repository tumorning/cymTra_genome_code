library(DESeq2)
library(pheatmap)
library(ggplot2)
library(ggVolcano)
library(RColorBrewer)


rm(list=ls())

setwd('/home/liunyw/project/orchid/DEG')

stages <- c('T2', 'T3', 'T4', 'T5', 'T6')

for ( stage in stages) {

    ##################################################################################################################################
    print('------------------------------------------------------')
    print(paste("T1 vs", stage, sep=" "))
    print('------------------------------------------------------')

    all_counts_df <- read.table(paste("/home/liunyw/project/orchid/DEG/T1vs", stage, ".counts.txt", sep=""), sep = "\t", quote="", header = T, check.names = F, row.name = 1)

    # head(all_counts_df)
    df_all_counts <- all_counts_df[rowSums(all_counts_df) != 0, ]

    df_meta <- read.table(paste("/home/liunyw/project/orchid/DEG/T1vs", stage, ".meta.txt", sep=""), sep = "\t", quote="", header = T, check.names = F, row.name = 1, stringsAsFactors = TRUE)
    #print(df_meta)


    ##################################################################################################################################

    ##################################################################################################################################
    #构建DESeqDataSet
    dds <- DESeqDataSetFromMatrix(countData = df_all_counts,
                                colData = df_meta,
                                design = ~ con)


    dds <- DESeq(dds)
    res <- results(dds)
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized=TRUE)
    write.table(normalized_counts, quote=FALSE, file=paste("/home/liunyw/project/orchid/DEG/deseq2.normalized.T1vs", stage, ".txt", sep="."), sep="\t")
    ##################################################################################################################################



    ##################################################################################################################################
    #Extracting transformed values
    rld <- rlog(dds, blind=FALSE)
    vsd <- vst(dds, blind=FALSE)

    #Heatmap of the sample-to-sample distances
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)
    pdf(paste("/home/liunyw/project/orchid/DEG/deseq2.pheatmap.T1vs", stage, ".pdf", sep="."))
    Heapmap <- pheatmap(rld_cor, annotation = df_meta)
    dev.off()

    #PCA
    pdf(paste("/home/liunyw/project/orchid/DEG/deseq2.pca.T1vs", stage, ".pdf", sep="."))
    PCA <- plotPCA(rld, intgroup = c("con"))
    print(PCA)
    dev.off()
    ##################################################################################################################################


    res <- results(dds, contrast=c("con", "control", "treated"))

    res <- res[order(res$pvalue),]

    write.csv(res, quote=FALSE, file=paste("/home/liunyw/project/orchid/DEG/deseq2.T1vs", stage, "all.gene.txt", sep="."))

    diff_gene <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
    write.csv(diff_gene, quote = FALSE, file=paste("/home/liunyw/project/orchid/DEG/deseq2.T1vs", stage, "different_gene.txt", sep="."))


    up_diff_gene <- subset(diff_gene, log2FoldChange > 1)
    write.csv(up_diff_gene, quote = FALSE, file = paste("/home/liunyw/project/orchid/DEG/deseq2.T1vs", stage, "different_up_gene.txt", sep="."))


    down_diff_gene <- subset(diff_gene, log2FoldChange < -1)
    write.csv(down_diff_gene, quote = FALSE, file = paste("/home/liunyw/project/orchid/DEG/deseq2.T1vs", stage, "different_downs_gene.txt", sep="."))
}