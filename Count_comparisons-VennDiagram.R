##### Venn Diagram for significant genes #####
### Subset significant Gene IDs by each treatment ###

AluVsConGenes <- row.names(AluVsConSig)
IndVsConGenes <- row.names(IndVsConSig)
AluIndVsConGenes <- row.names(AluIndVsConSig)

vennPlot <- venn.diagram(list(AluIndVsConGenes, AluVsConGenes, IndVsConGenes), NULL, fill=c("red", "green", "blue"), alpha=c(0.5,0.5,0.5), cex=3, cat.fonface=4, category.names = c("AluIndVsCon","AluVsCon", "IndVsCon"))
grid.draw(vennPlot)

dev.off() ## Export / Save the image before running this line.

##### Subset up / down regulated genes #####
### Up regulated genes ###
AluVsConUp <- subset(AluVsCon, AluVsCon$log2FoldChange > 0 & AluVsCon$padj <= padjCutoff)
mean(AluVsConUp$log2FoldChange)

IndVsConUp <- subset(IndVsCon, IndVsCon$log2FoldChange > 0 & IndVsCon$padj <= padjCutoff)
mean(IndVsConUp$log2FoldChange)

AluIndVsConUp <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange > 0 & AluIndVsCon$padj <= padjCutoff)
mean(AluIndVsConUp$log2FoldChange)

##### Venn Diagram for Up regulated genes #####
AluVsConUpGenes <- row.names(AluVsConUp)
IndVsConUpGenes <- row.names(IndVsConUp)
AluIndVsConUpGenes <- row.names(AluIndVsConUp)

vennPlot <- venn.diagram(list(AluIndVsConUpGenes, AluVsConUpGenes, IndVsConUpGenes), NULL, fill=c("red", "green", "blue"), alpha=c(0.5,0.5,0.5), cex=3, cat.fonface=4, category.names = c("AluIndVsCon","AluVsCon", "IndVsCon"))
grid.draw(vennPlot)

dev.off() ## Export / Save the image before running this line.

### Down regulated genes ###
AluVsConDown <- subset(AluVsCon, AluVsCon$log2FoldChange < 0 & AluVsCon$padj <= padjCutoff)
mean(AluVsConDown$log2FoldChange)

IndVsConDown <- subset(IndVsCon, IndVsCon$log2FoldChange < 0 & IndVsCon$padj <= padjCutoff)
mean(IndVsConDown$log2FoldChange)

AluIndVsConDown <- subset(AluIndVsCon, AluIndVsCon$log2FoldChange < 0 & AluIndVsCon$padj <= padjCutoff)
mean(AluIndVsConDown$log2FoldChange)

##### Venn Diagram for Down regulated genes #####
AluVsConDownGenes <- row.names(AluVsConDown)
IndVsConDownGenes <- row.names(IndVsConDown)
AluIndVsConDownGenes <- row.names(AluIndVsConDown)

vennPlot <- venn.diagram(list(AluVsConDownGenes, IndVsConDownGenes, AluIndVsConDownGenes), NULL, fill=c("red", "green", "blue"), alpha=c(0.5,0.5,0.5), cex=3, cat.fonface=4, category.names = c("AluIndVsCon","AluVsCon", "IndVsCon"))
grid.draw(vennPlot)

## venn.diagram(list(AluVsConDownGenes, IndVsConDownGenes, AluIndVsConDownGenes), "ex1.png", fill=c("red", "green", "blue"), alpha=c(0.5,0.5,0.5), cex=3, cat.fonface=4, category.names = c("AluIndVsCon","AluVsCon", "IndVsCon")) ##

dev.off() ## Export / Save the image before running this line.

##### Extract significant gene IDs #####

write.table(AluVsConGenes,"Significant_AL_Genes.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(IndVsConGenes,"Significant_IN_Genes.txt",row.names=F,col.names=F,quote=F,sep="\t")
write.table(AluIndVsConGenes,"Significant_ALIN_Genes.txt",row.names=F,col.names=F,quote=F,sep="\t")

SharedGeneID <- intersect(AluVsConGenes,IndVsConGenes)
write.table(SharedGeneID,"AL_and_IN_Shared_Genes.txt",row.names=F,col.names=F,quote=F,sep="\t")

SharedGeneID <- intersect(AluVsConUpGenes,IndVsConUpGenes)
write.table(SharedGeneID,"AL_and_IN_Shared_Up_regulated_Genes.txt",row.names=F,col.names=F,quote=F,sep="\t")

SharedGeneID <- intersect(AluVsConDownGenes,IndVsConDownGenes)
write.table(SharedGeneID,"AL_and_IN_Shared_Down_regulated_Genes.txt",row.names=F,col.names=F,quote=F,sep="\t")