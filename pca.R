rm(list = ls())

# load libraries
library(ggplot2)
library(plotly)


# -----------------------------------------------------
# 1. load sample sheet metadata
# -----------------------------------------------------
sample_sheet_file <- 'sample_sheet.csv'
sample_sheet <- read.table(sample_sheet_file, header = T, sep = ',', row.names = 1)

# take a subset if needed
sample_sheet <- sample_sheet[sample_sheet$Condition %in% c('Control','Indium','Aluminium_Indium'),]

# check the loaded metadata
head(sample_sheet)
dim(sample_sheet)


# -----------------------------------------------------
# 2. load read counts
# -----------------------------------------------------
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
p <- pca.plot(read_counts, groups, comps = c(1,2), ntop = 2000)
p
ggplotly(p)
