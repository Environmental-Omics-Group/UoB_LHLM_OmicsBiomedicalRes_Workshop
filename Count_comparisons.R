##### source for installing packages #####
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# install packages eg:-
# BiocManager::install("plotly") 
# install.packages("VennDiagram")

##### load installed libraries to use in project #####

library(plotly)
library(VennDiagram)

##### Read input file #####
df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)
dim(df)
View(df)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])

##### Count and subset of Differential Expressed Genes (DEGs) #####
### Number of significant DEG IDs by P-value and P-adjusted (Q) cutoffs ###

cat(paste0('significant p-values = ', sum(AluVsCon$pvalue <= 0.05, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(AluVsCon$padj <= 0.05, na.rm = T), '\n'))

cat(paste0('significant p-values = ', sum(IndVsCon$pvalue <= 0.05, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(IndVsCon$padj <= 0.05, na.rm = T), '\n'))

cat(paste0('significant p-values = ', sum(AluIndVsCon$pvalue <= 0.05, na.rm = T), '\n'))
cat(paste0('significant p-adjust = ', sum(AluIndVsCon$padj <= 0.05, na.rm = T), '\n'))

##### Set cutoff to filter #####
padjCutoff <- 0.05

##### Subset significant DEG #####
AluVsConSig <- subset(AluVsCon, AluVsCon$padj <= padjCutoff)
IndVsConSig <- subset(IndVsCon, IndVsCon$padj <= padjCutoff)
AluIndVsConSig <- subset(AluIndVsCon, AluIndVsCon$padj <= padjCutoff)

##### Subset up regulated genes #####
AluVsConUp <- subset(AluVsCon, AluVsCon$log2FoldChange > 0 & AluVsCon$padj <= padjCutoff)
mean(AluVsConUp$log2FoldChange)

IndVsConUp <- subset(IndVsCon, IndVsCon$log2FoldChange > 0 & IndVsCon$padj <= padjCutoff)
mean(IndVsConUp$log2FoldChange)

AluIndVsConUp <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange > 0 & AluIndVsCon$padj <= padjCutoff)
mean(AluIndVsConUp$log2FoldChange)

##### Subset down regulated genes #####
AluVsConDown <- subset(AluVsCon, AluVsCon$log2FoldChange < 0 & AluVsCon$padj <= padjCutoff)
mean(AluVsConDown$log2FoldChange)

IndVsConDown <- subset(IndVsCon, IndVsCon$log2FoldChange < 0 & IndVsCon$padj <= padjCutoff)
mean(IndVsConDown$log2FoldChange)

AluIndVsConDown <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange < 0 & AluIndVsCon$padj <= padjCutoff)
mean(AluIndVsConDown$log2FoldChange)
