#----Differentiall Gene Expression Analysis----
resultsNames(dds)

condition <- 'IRE1_KO'
control <- 'IRE1_WT'
comparison <- paste0('condition_', condition, '_vs_', control)

resultsNames(dds)[2] == comparison

threshold <- log2(1.5)

dir.create(paste0(mainDir, './plots/', comparison), recursive = F)

# Run DEG analysis----
res <- results(dds, name = comparison, lfcThreshold = threshold)
summary(res)
sum(res$padj < 0.1, na.rm = T)

# Add more info to results: Symbols, GO, ENTREZ, etc.
library(org.Hs.eg.db)

res$symbol <- mapIds(org.Hs.eg.db, keys = row.names(res), column = 'SYMBOL',
                     keytype = 'ENSEMBL', multiVals = 'first')

res$name <- mapIds(org.Hs.eg.db, keys = row.names(res), column = 'GENENAME',
                     keytype = 'ENSEMBL', multiVals = 'first')

res$entrez <- mapIds(org.Hs.eg.db, keys = row.names(res), column = 'ENTREZID',
                     keytype = 'ENSEMBL', multiVals = 'first')

res$go <- mapIds(org.Hs.eg.db, keys = row.names(res), column = 'GO',
                     keytype = 'ENSEMBL', multiVals = 'first')

# Sort the most DEGs, filter and save
res <- res[order(res$pvalue), ]
res.sig <- subset(res, padj < 0.1)

dir.create(paste0(mainDir, './results/'), recursive = F)
write.csv(as.data.frame(res), file = paste0(mainDir, './results/', comparison, '.csv'))
write.csv(as.data.frame(res.sig), file = paste0(mainDir, './results/', comparison, '.significant.csv'))

#----Volcano plot----
library(EnhancedVolcano)
pdf(paste0('./plots/', comparison, '/volcano.', comparison, '.pdf'))
EnhancedVolcano(res, lab = res$symbol, pCutoff = 0.1, x = 'log2FoldChange',
                title = comparison, y = 'padj', legendPosition = 'right',
                FCcutoff = threshold, drawConnectors = F, boxedLabels = F,
                subtitle = '', caption = '', labSize = 2, pointSize = 1,
                colAlpha = 0.5) + theme(aspect.ratio = 1, panel.background = element_rect(linetype = NULL))
dev.off()

#----session----
sessionInfo()
