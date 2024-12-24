
rm(list=ls())
setwd("/home/liunyw/project/tml/0506_wgcna/")

library(WGCNA)
library(reshape2)
library(stringr)

options(stringsAsFactors = FALSE)

enableWGCNAThreads()

#read files
exprMat <- "GOODFPKM.txt"
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
robustY = ifelse(corType=="pearson",T,F)

datExpr <- read.table(exprMat, sep='\t', row.names=1, header=T,
                      quote="", comment="", check.names=F)
dim(datExpr)
traitMat <- "NEWtraits1102.csv"
datTrait <- read.table(traitMat, sep=",", stringsAsFactors = T, row.names=1, header=T, quote="", comment="", check.names=F)
datTraits = data.frame(datTrait)


m.mad <- apply(datExpr,1,mad)
datExprVar <- datExpr[which(m.mad >
                  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
dim(datExprVar)

datExpr <- as.data.frame(t(datExprVar))

gsg = goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK){
if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr)#
head(datExpr)[,1:8]

###
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

###
sampleTree = hclust(dist(datExpr), method = "average")
pdf(file="0506_sample_clustering.pdf",width = 14,height = 10)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
dev.off()

###
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed=FALSE);
# Plot the sample dendrogram and the colors underneath.
pdf(file="0506_Sample_dendrogram_and_trait_heatmap.pdf",width = 14,height = 10)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

###
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector=powers, 
                        networkType=type, verbose=5)
power=sft$powerEstimate
power #
# Plot the results:
pdf(file="0506_Scale_independence_and_Mean_connectivity.pdf",width = 14,height = 10)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
       
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h       
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###
pdf(file="0506_K.pdf",width = 14,height = 10)
k <- softConnectivity(datExpr,power=power)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k,main="Check Scale free topology\n")
dev.off()

###
#one-step
net = blockwiseModules(datExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, "0506.tom"),
                       verbose = 3)
table(net$colors)
###
0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15 
3179 5155 4032 2955 1762 1317  870  800  660  638  621  591  461  402  380  316 
  16   17   18   19   20   21   22   23   24 
 224  212  193  178  108  107   85   59   59
# open a graphics window
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
MEs = orderMEs(net$MEs);
geneTree = net$dendrograms[[1]];
pdf(file="0506_Module_colors.pdf",width = 14,height = 10)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

###
unmergedColors = labels2colors(net$unmergedColors)
mergedColors   = labels2colors(net$colors)
pdf(file="0506_Mergedcolors.pdf",width = 14,height = 10)
plotDendroAndColors(net$dendrograms[[1]],
 cbind(unmergedColors[net$blockGenes[[1]]], mergedColors[net$blockGenes[[1]]]),
 c("Dynamic Tree Cut" , "Merged colors"),
 dendroLabels = FALSE,
 hang = 0.03,
 addGuide = TRUE,
 guideHang = 0.05
)
dev.off()

###
pdf(file="0516_Eigengene_adjacency_heatmap.pdf",width = 14,height = 10)
plotEigengeneNetworks(MEs, 
                      "Eigengene adjacency heatmap", 
                      marHeatmap = c(3,4,2,2), 
                      plotDendrograms = FALSE, 
                      xLabelsAngle = 90)
dev.off()





###
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

traitMat <- "NEWtraits0506.csv"
datTrait <- read.table(traitMat, sep=',', stringsAsFactors = T, row.names=1, header=T, quote="", comment="", check.names=F)
datTraits = data.frame(datTrait)
trait <- "NEWtraits0506.csv"
if(trait != "") {
  traitData <- read.table(file=trait, sep=",", header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(datExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
# Will display correlations and their p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
# Display the correlation values within a heatmap plot这一步是画最关键的性状和模块关联的图
png(file="0506_NEW_Module_and_trait_pearson.png",width = 1400,height = 1000)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData),
                 yLabels = colnames(MEs_col),
                 cex.lab = 0.5,
                 ySymbols = colnames(MEs_col), colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix, setStdMargins = FALSE,
                 cex.text = 0.7, zlim = c(-1,1),
                 main = paste("Module-trait relationships"))
dev.off()
#
which.trait <- "Caryophyllene"
modTraitCor[, which.trait]



###
moduleColors = labels2colors(net$colors)
unique(moduleColors)
head(names(datExpr)[moduleColors=="black"])#change color
###
pdf(file="0411_BlackHEATMAP.pdf",width = 14,height = 10)
which.module="black"
plotMat(
t(scale(datExpr[,moduleColors==which.module ]) ),
nrgcols=30,
rlabels=F,
rcols=which.module,
main=which.module,
cex.main=2
)
dev.off()

### Isolate traits from the clinical traits，Caryphyllene as example：
Caryophyllene = as.data.frame(datTraits$Caryophyllene);
names(Caryophyllene) = "Caryophyllene"
MET = orderMEs(cbind(MEs, Caryophyllene))
pdf("0519_MMvsGS_pink.pdf");
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
= 90)
dev.off()

###MEblack and Caryophyllene 
pdf("0411_brown_module_heatmap_and_the_eigengene.pdf")
which.module="brown"
ME=MEs_col[,paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="MPP")
dev.off()


###Gene_significance_across_modules
png(file="0411_Gene_significance_across_modules_Caryophyllene.png",width = 1400,height = 1000)
moduleColors = labels2colors(net$colors)
which.trait <- "Caryophyllene"
y <- datTraits[, which.trait]
GS <- as.numeric(cor(y ,datExpr, use="p"))
GeneSignificance <-  abs(GS)
ModuleSignificance <- tapply(
GeneSignificance,
moduleColors, mean, na.rm=T)
plotModuleSignificance(GeneSignificance, moduleColors)
dev.off()
###
datKME = signedKME(
datExpr,
MEs,
outputColumnName="MM.")
FilterGenes= abs(GS)> .2 & abs(datKME$MM.6)>.8
###result as:
# [1] "blue"          "cyan"          "brown"         "yellow"       
# [5] "salmon"        "black"         "magenta"       "turquoise"    
# [9] "red"           "purple"        "tan"           "pink"         
# [13] "grey"          "green"         "lightyellow"   "grey60"       
# [17] "lightgreen"    "lightcyan"     "darkturquoise" "midnightblue" 
# [21] "darkgrey"      "greenyellow"   "darkgreen"     "darkred"      
# [25] "royalblue"
####
dimnames(data.frame(datExpr))[[2]][FilterGenes]
trait_hubGenes<-colnames(datExpr)[FilterGenes]
###
pdf("0517_hub_unsigned_correlations.pdf")
plotNetworkHeatmap(datExpr,plotGenes = trait_hubGenes,networkType = "unsigned",useTOM = TRUE,power=power,main="unsigned correlations")
dev.off()
###
write.table(dimnames(data.frame(datExpr))[[2]][FilterGenes],"0108_black_filtergenes1.txt",quote=F)

###
TOM = TOMsimilarityFromExpr(datExpr, power = 18)
modules = c("brown")
genes = colnames(datExpr)
inModule = is.finite(match(moduleColors,modules))
inModule = (moduleColors==module);
modgenes = genes[inModule]
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modgenes, modgenes)
cyt = exportNetworkToCytoscape(
                                 modTOM,
                                 edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), "-module_edges.txt", sep=""),
                                 nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), "-module_nodes.txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.02,
                                 nodeNames = modgenes,
                                 nodeAttr = moduleColors[inModule]
);

save.image("0411WGCNA.RData")