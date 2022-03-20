# Comparison of logarithmic distribution across different preservation techniques
# The establised cell preservation techniques implemented for the study are the following:
# "Fresh", "cryopreserved", "O.N. cultured","2-day cultured","6-day cultured"
# In this case, the script comparers Fresh cells to cryopreserved cells
# find entries (cells) which are either Fresh OR cryopreserved

subset_idx = which( sce$source == "Fresh" | sce$source == "cryopreserved" )

# do findMarkers only on these cells
out = findMarkers(sce[, subset_idx], sce$source[subset_idx])
out[["Fresh"]]

#generate a text
write.table(out[["Fresh"]], file="Fresh_to_cryopreserved.txt", sep="\t", quote=FALSE)

# Expression plots
plotExpression(sce, features=c("GAPDH", "EPCAM"), x="source", colour_by="source")

# FTE cell marker genes expression plots
plotExpression(sce, features=c("KRT7", "PAX8", "CCDC17", "CCDC78"), x="source", colour_by="source")

# Wilcoxon text code example
wilcox.test(logCount["KRT7", which(sce$source=="Fresh")], logCount["KRT7", which(sce$source=="Fresh")])

# Correlation plot example
geneA = "KRT7"
geneB = "PAX8"

logCounts = logcounts(sce)

df = data.frame( geneA = logCounts[geneA,], geneB = logCounts[geneB,], sources = sce$source)
ggplot(df, aes(x=geneA, y=geneB, colour=sce$source)) + geom_point() + xlab("GeneA") + ylab("GeneB")

# Spearman correlation using the same df
install.packages("Hmisc")
library("Hmisc")

idx = which(df$geneA > 0 & df$geneB > 0)
# calculating the correlation using rcorr() function
rcorr(df$geneA[idx], df$geneB[idx], type = "spearman")

