#! /home/liunyw/miniforge3/envs/eggnog/bin/R

rm(list=ls())
setwd('~/cafe/enrichment')
library('clusterProfiler')

#GO
decrease_geneid_pvalue <- read.table('Decrease.geneid.pvalue.txt')
enrichgene <- as.factor(decrease_geneid_pvalue$V1)
term2gene <- read.csv(file = 'GO2gene.txt', header=F, sep=",")
term2name <- read.table(file = 'go2term.txt', header=F, sep="\t", quote = "")
x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
#x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
write.table(x, file = 'Decrease.geneid.pvalue.GOenrichment.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(file = 'Decrease.geneid.pvalue.GOenrichment.bar.pdf', height = 12, width = 16)
barplot(x,showCategory=21,drop=T)
dev.off()

#KEGG
decrease_geneid_pvalue <- read.table('Decrease.geneid.pvalue.txt')
enrichgene <- as.factor(decrease_geneid_pvalue$V1)
term2gene <- read.csv(file = 'ko2gene.txt', header=F, sep="\t")
term2name <- read.table(file = 'ko2name.txt', header=F, sep="\t", quote = "")
# x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
write.table(x, file = 'Decrease.geneid.pvalue.KEGGenrichment.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(file = 'Decrease.geneid.pvalue.KEGGenrichment.bar.pdf', height = 12, width = 16)
barplot(x,showCategory=30,drop=T)
dev.off()



#GO 非显著扩张
decrease_geneid <- read.table('Decrease.geneid.txt')
enrichgene <- as.factor(decrease_geneid$V1)
term2gene <- read.csv(file = 'GO2gene.txt', header=F, sep=",")
term2name <- read.table(file = 'go2term.txt', header=F, sep="\t", quote = "")
x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
#x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
write.table(x, file = 'Decrease.geneid.GOenrichment.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(file = 'Decrease.geneid.GOenrichment.bar.pdf', height = 12, width = 16)
barplot(x,showCategory=21,drop=T)
dev.off()

#KEGG
decrease_geneid <- read.table('Decrease.geneid.txt')
enrichgene <- as.factor(decrease_geneid$V1)
term2gene <- read.csv(file = 'ko2gene.txt', header=F, sep="\t")
term2name <- read.table(file = 'ko2name.txt', header=F, sep="\t", quote = "")
# x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
write.table(x, file = 'Decrease.geneid.KEGGenrichment.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(file = 'Decrease.geneid.KEGGenrichment.bar.pdf', height = 12, width = 16)
barplot(x,showCategory=30,drop=T)
dev.off()

