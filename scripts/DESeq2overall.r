#----Create DESeq2 object----
library(DESeq2)

dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~condition+batch)

# Filter low counts and run DEG analysis----
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$condition <- relevel(dds$condition, ref = 'IRE1_WT') # define control

dds <- DESeq(dds)

dir.create(paste0(mainDir, 'rds/'), recursive = F)
saveRDS(dds, file = './rds/dds.rds')

# Normalise dds counts with varianceStabilizingTransformation----
vsd <- varianceStabilizingTransformation(dds, blind = F) # if more than 100 samples use vst() function

# Remove batch effects----
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch = vsd$batch)

saveRDS(vsd, file = './rds/vsd.rds')

#----Visualisations----
dir.create(paste0(mainDir, 'plots/'), recursive = F)

library(ggplot2)

# PCA plot----
pdf('./plots/pca.pdf')
plotPCA(vsd, intgroup = 'condition') + theme(aspect.ratio = 1, panel.background = NULL)
dev.off()

# Top variable genes heat map----

topVarGenes <- head(order(-rowVars(assay(vsd))), 500)
mat <- assay(vsd)[topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[ , 'condition'])
colnames(df) <- 'Treatment'
rownames(df) <- colnames(vsd)

pdf(paste0(mainDir, './plots/heatmap.pdf'))
pheatmap(mat, cluster_rows = T, show_rownames = F, cluster_cols = T,
         annotation_col = df, angle_col = 45, cutree_cols = T, border_color = NA,
         scale = 'row')
dev.off()

# Sample distances----
library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

pdf(paste0(mainDir, './plots/sampleDistances.pdf'))
pheatmap(sampleDistMatrix, cluster_rows = F, cluster_cols = F, angle_col = 45,
         annotation_col = df, border_color = NA, cellwidth = 15, cellheight = 15)
dev.off()

# Cluster tree---
library(matrixStats)
rv <- rowVars(assay(vsd))
o <- order(rv, decreasing = T)
dists <- dist(t(assay(vsd)[head(o, 500), ]))
hc <- hclust(dists)

pdf('./plots/dendrogram.pdf')
plot(hc, labels = vsd$run)
dev.off()

#----session----
sessionInfo()
