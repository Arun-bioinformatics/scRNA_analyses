# Setup and main plots.

# Run the following instalation if packages "scran" and "scater" are not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scater", version = "3.8")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scran", version = "3.8")

install.packages("pheatmap")

# Load libraries onto Rstudio
library(scater) 
library(scran)
library(pheatmap)

# Loading libraries and data into the RStudio programme: our data file was saved with the name 
# "ft_data.rds", however the name "sce" has been assigned as a simplified way to refer to our 
# datafile

sce = readRDS("ft_data.rds") 

# Establishing preservation methods as levels of our data
preservation_method = cellInfo$source 
colData(sce)$source = factor(colData(sce)$source, levels = c("Fresh", 
                                                             "cryopreserved", 
                                                             "O.N. cultured",
                                                             "2-day cultured",
                                                             "6-day cultured"))
levels(colData(sce)$source)

# 1) T-SNE plot distinguishing cells by preservation
set.seed(100)
sce = runTSNE(sce)
plotTSNE(sce, colour_by = "source")

# 2) T-SNE plot distinguishing cells by batch origin
set.seed(100)
sce = runTSNE(sce)
plotTSNE(sce, colour_by = "batch")

# 3) Cell population principal component analysis 
sce = runPCA(sce, ncomponents = 2)
plotPCA(sce, colour_by = "source", ncomponents = 2)

# 4) Heatmap of distinguishingly expressed genes across each cell preservation techniques:
logCount = logcounts(sce)
logCount_var = matrixStats::rowVars(logcounts(sce))
# We now want to find the 50 most variable genes.
gene_idx = order(logCount_var, decreasing=TRUE)
gene_idx = gene_idx[1:50]
# We want to sort the cells by preservation method
col_idx = order(colData(sce)$source) 
# We now create a data frame containing the preservation method information by copying a reordered version of the source cell information column:
my_sample_col <- data.frame(method = colData(sce)$source[col_idx])  
row.names(my_sample_col) <- colData(sce)$Sample[col_idx]
# We will therefore scale all the measurements to ensure that different genes are visually comparable:
logCount_scaled = t(scale( t(logCount)))
# Heatmap code
pheatmap(logCount_scaled[gene_idx, col_idx], 
         cluster_rows=FALSE, cluster_cols=FALSE, 
         annotation_col = my_sample_col)


