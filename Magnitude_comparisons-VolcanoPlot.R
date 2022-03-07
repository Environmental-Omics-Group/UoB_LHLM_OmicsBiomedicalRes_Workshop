##### Volcano plot ######
##### Load installed libraries to use in project #####

library(plotly)

##### Read input file #####
df <- read.csv("DESeq2_results_all.csv", header = T, row.names = 1, check.names = F)
dim(df)

AluIndVsCon <- as.data.frame(df[,1:3])
AluVsCon <- as.data.frame(df[,4:6])
IndVsCon <- as.data.frame(df[,7:9])

##### Volcano plot ######
### Alu vs Con ###
valid_index <- which(!is.na(AluVsCon$padj))
df <- data.frame(pvals = AluVsCon$padj[valid_index], lfcs = AluVsCon$log2FoldChange[valid_index], gene = rownames(AluVsCon)[valid_index])
df$significant = as.factor(df$pvals <= 0.05 & abs(df$lfcs) >= 2)

p <- ggplot(df, aes(x = lfcs, y = -log10(pvals), colour = significant, gene = gene)) +
  geom_point(alpha = 0.8, size = 1) +
  labs(title = 'Aluminium Volcano Plot',
       x = 'log2 fold change',
       y = '-log10 p-value') +
  scale_color_manual(values = c('#619CFF', '#F564E3'))

ggplotly(p)

dev.off() ## Export / Save the image before running this line.

### Ind vs Con ###
valid_index <- which(!is.na(IndVsCon$padj))
df <- data.frame(pvals = IndVsCon$padj[valid_index], lfcs = IndVsCon$log2FoldChange[valid_index], gene = rownames(IndVsCon)[valid_index])
df$significant = as.factor(df$pvals <= 0.05 & abs(df$lfcs) >= 2)

p <- ggplot(df, aes(x = lfcs, y = -log10(pvals), colour = significant, gene = gene)) +
  geom_point(alpha = 0.8, size = 1) +
  labs(title = 'Indium Volcano Plot',
       x = 'log2 fold change',
       y = '-log10 p-value') +
  scale_color_manual(values = c('#619CFF', '#F564E3'))

ggplotly(p)

dev.off() ## Export / Save the image before running this line.

### AluInd vs Con ###
valid_index <- which(!is.na(AluIndVsCon$padj))
df <- data.frame(pvals = AluIndVsCon$padj[valid_index], lfcs = AluIndVsCon$log2FoldChange[valid_index], gene = rownames(AluIndVsCon)[valid_index])
df$significant = as.factor(df$pvals <= 0.05 & abs(df$lfcs) >= 2)

p <- ggplot(df, aes(x = lfcs, y = -log10(pvals), colour = significant, gene = gene)) +
  geom_point(alpha = 0.8, size = 1) +
  labs(title = 'AluInd Volcano Plot',
       x = 'log2 fold change',
       y = '-log10 p-value') +
  scale_color_manual(values = c('#619CFF', '#F564E3'))

ggplotly(p)

dev.off() ## Export / Save the image before running this line.