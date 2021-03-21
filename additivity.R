rm(list = ls())

# load libraries
library(ggplot2)
library(plotly)


# -----------------------------------------------------
# 1. load DESeq2 results
# -----------------------------------------------------
deseq.df <- read.table('DESeq2_results_all.csv', header = T, sep = ',',  row.names = 1, check.names = F)
  
aluminium_indium_dres <- deseq.df[,1:3]
aluminium_dres <- deseq.df[,4:6]
indium_dres <- deseq.df[,7:9]

# check the loaded values
head(aluminium_indium_dres)
dim(aluminium_indium_dres)

head(aluminium_dres)
dim(aluminium_dres)

head(indium_dres)
dim(indium_dres)


# -----------------------------------------------------
# 2. calculate observed and estimated log2FC for DEGs 
# -----------------------------------------------------
aluminium_indium_lfc <- abs(aluminium_indium_dres$log2FoldChange)
aluminium_lfc <- abs(aluminium_dres$log2FoldChange)
indium_lfc <- abs(indium_dres$log2FoldChange)

observed_lfc <- aluminium_indium_lfc
expected_lfc <- (aluminium_lfc + indium_lfc)/2


# -----------------------------------------------------
# 3. find the overlapping DEGs
# -----------------------------------------------------
# DEGs on all three comparisons
degs <- (aluminium_indium_dres$padj < 0.05) & (aluminium_dres$padj < 0.05) & (indium_dres$padj < 0.05)
# OR, DEGs on any two comparisons but not in all three comparisons
degs <- (
  ((aluminium_indium_dres$padj < 0.05) & (aluminium_dres$padj < 0.05)) | 
  ((aluminium_indium_dres$padj < 0.05) & (indium_dres$padj < 0.05)) |
  ((aluminium_dres$padj < 0.05) & (indium_dres$padj < 0.05))
) & !((aluminium_indium_dres$padj < 0.05) & (aluminium_dres$padj < 0.05) & (indium_dres$padj < 0.05))
# OR, DEGs on any two comparisons including those in all three comparisons
degs <- ((aluminium_indium_dres$padj < 0.05) & (aluminium_dres$padj < 0.05)) | 
        ((aluminium_indium_dres$padj < 0.05) & (indium_dres$padj < 0.05)) |
        ((aluminium_dres$padj < 0.05) & (indium_dres$padj < 0.05))

# remove na and outliers
degs[is.na(degs)] <- F
degs[(aluminium_indium_lfc < quantile(aluminium_indium_lfc, 0.001)) | (aluminium_indium_lfc > quantile(aluminium_indium_lfc, 0.999))] <- F
degs[(aluminium_lfc < quantile(aluminium_lfc, 0.001)) | (aluminium_lfc > quantile(aluminium_lfc, 0.999))] <- F
degs[(indium_lfc < quantile(indium_lfc, 0.001)) | (indium_lfc > quantile(indium_lfc, 0.999))] <- F

# show the number of DEGs
sum(degs)


# -----------------------------------------------------
# 4. plot additivity and fit linear model 
# -----------------------------------------------------
plot.data <- data.frame(Observed = observed_lfc[degs], Expected = expected_lfc[degs])

p <- ggplot(plot.data, aes(x = Observed, y = Expected)) +
     geom_point() +
     geom_smooth(method = 'lm', formula = y~x) +
     geom_line(aes(y = Observed), color = 'red')
p
#ggplotly(p)
