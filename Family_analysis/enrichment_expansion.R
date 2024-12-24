#! /home/liunyw/miniforge3/envs/eggnog/bin/R

rm(list=ls())
setwd('~/cafe/enrichment')
library('clusterProfiler')

#GO
increase_geneid_pvalue <- read.table('Increase.geneid.pvalue.txt')
enrichgene <- as.factor(increase_geneid_pvalue$V1)
term2gene <- read.csv(file = 'GO2gene.txt', header=F, sep=",")
term2name <- read.table(file = 'go2term.txt', header=F, sep="\t", quote = "")
x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
write.table(x, file = 'Increase.geneid.pvalue.enrichment.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(file = 'Increase.geneid.pvalue.enrichment.bar.pdf', height = 12, width = 16)
barplot(x,showCategory=20,drop=T)
dev.off()

#KEGG
increase_geneid_pvalue <- read.table('Increase.geneid.pvalue.txt')
enrichgene <- as.factor(increase_geneid_pvalue$V1)
term2gene <- read.csv(file = 'ko2gene.txt', header=F, sep="\t")
term2name <- read.table(file = 'ko2name.txt', header=F, sep="\t", quote = "")
x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
write.table(x, file = 'Increase.geneid.pvalue.KEGGenrichment.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(file = 'Increase.geneid.pvalue.KEGGenrichment.bar.pdf', height = 12, width = 16)
barplot(x,showCategory=20,drop=T)
dev.off()




#GO 非显著扩张
increase_geneid <- read.table('Increase.geneid.txt')
enrichgene <- as.factor(increase_geneid_pvalue$V1)
term2gene <- read.csv(file = 'GO2gene.txt', header=F, sep=",")
term2name <- read.table(file = 'go2term.txt', header=F, sep="\t", quote = "")
x <- enricher(enrichgene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.01, pAdjustMethod = "BH", qvalueCutoff = 0.05)
write.table(x, file = 'Increase.geneid.enrichment.txt', sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
pdf(file = 'Increase.geneid.enrichment.bar.pdf', height = 12, width = 16)
barplot(x,showCategory=20,drop=T)
dev.off()



