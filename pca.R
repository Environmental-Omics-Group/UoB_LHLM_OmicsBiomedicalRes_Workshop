rm(list = ls())

# load libraries
library(ggplot2)
library(plotly)

# set work dir
setwd("~/MScTox_BIOM33_Workshop/Session 2")


# -----------------------------------------------------
# 1. load sample sheet metadata
# -----------------------------------------------------
sample_sheet_file <- 'sample_sheet.csv'
sample_sheet <- read.table(sample_sheet_file, header = T, sep = ',', row.names = 1)

# take a subset if needed
# =====================================
# QUESTION 1: Please try to used different condition subsets to plot your PCA:
#   a. Control, Indium, and Aluminium_Indium 
#   b. Control, Aluminium, and Aluminium_Indium 
#   c. Control, Indium, and Aluminium
#   d. Indium, Aluminium, and Aluminium_Indium
# How does the PCA plot change?
# =====================================
sample_sheet <- sample_sheet[sample_sheet$Condition %in% c('Control','Indium','Aluminium_Indium'),]

# check the loaded metadata
head(sample_sheet)
dim(sample_sheet)


# -----------------------------------------------------
# 2. load read counts
# -----------------------------------------------------
# =====================================
# QUESTION 2: Please try to used different input matrices to plot your PCA:
#   a. raw counts 
#   b. normalised counts
#   c. VST counts
# How does the PCA plot change?
# =====================================
read_counts_file <- 'gene_norm_counts.csv'
read_counts <- read.table(read_counts_file, header = T, sep = ',',  row.names = 1, check.names = F)

# convert to a matrix (and round to integer)
read_counts <- as.matrix(read_counts)

# check the loaded values
head(read_counts)
dim(read_counts)


# -----------------------------------------------------
# 3. match the order of sample sheet and read counts
# -----------------------------------------------------
# match sample names
match_index <- match(sample_sheet$Sample_Name, colnames(read_counts))
read_counts <- read_counts[,match_index]

# check if there are any missing samples
stopifnot(sample_sheet$Sample_Name == colnames(read_counts))
stopifnot(!any(is.na(sample_sheet$Sample_Name)))
stopifnot(!any(is.na(colnames(read_counts))))

# check the matched values
head(read_counts)
dim(read_counts)


# -----------------------------------------------------
# 4. Define PCA plot
# -----------------------------------------------------
pca.plot <- function(read.counts, classes, 
                     comps = c(1, 2), ntop = min(500, nrow(read.counts)), standard = T,
                     col = c('lightblue', 'orange', 'MediumVioletRed', 'SpringGreen')){
  top_index <- order(apply(read.counts, 1, var), decreasing = TRUE)[1:ntop]
  pca <- prcomp(scale(t(read.counts[top_index,]), center = standard, scale = standard))
  prp <- pca$sdev^2 * 100 / sum(pca$sdev^2)
  
  pca_comps <- pca$x[,comps]
  prp_comps <- round(prp[comps], 2)
  
  df <- data.frame(pc1 = pca_comps[,1], pc2 = pca_comps[,2], condition = classes)
  p  <- ggplot(df, aes(x = pc1, y = pc2, color = condition)) + 
    geom_point(size = 3) + 
    labs(title = paste0('Principal Component Analysis - Axes ', comps[1] , ' and ', comps[2]), 
         x = paste0('PC', comps[1], ' (', prp_comps[1], '%)'), 
         y = paste0('PC', comps[2], ' (', prp_comps[2], '%)')) + 
    geom_text(label = colnames(read.counts), vjust = 0, nudge_y = 1) +
    scale_color_manual(values = col)
  return(p)
}


# -----------------------------------------------------
# 5. Plot PCA on the read counts
# -----------------------------------------------------
groups <- sample_sheet$Condition

# =====================================
# QUESTION 3: Please try to used different number of top genes (with the largest variants) to plot your PCA:
#   a. a small number, e.g., 100
#   b. half of the genes
#   c. all the genes
# How does the PCA plot change?
# =====================================
p <- pca.plot(read_counts, groups, comps = c(1,2), ntop = 2000)
p
# ggplotly(p)


library(ggplot2)
library(stringr)
library(dplyr)


rna_norm_numeric <- select(rna_norm_counts, -X)

normrna_pca <- prcomp(t(rna_norm_numeric), center = TRUE, scale = TRUE)
pca_summary <- summary(normrna_pca)
print(pca_summary)

pca_df <- data.frame(normrna_pca$x, cond = sample_sheet[row.names(normrna_pca$x),]$REF) #, make = stringr::word(rownames(rna_norm_numeric), 1))

ggplot(pca_df, aes(x = PC1, y = PC2, col = cond)) +
  geom_point(size = 3) +
  labs(x = "PC1", y= "PC2", title = "PCA for normalised RNA data" ) + theme(legend.position = "bottom")
