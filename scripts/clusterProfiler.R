library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(dplyr)

#----Read data and create a directory for plots----
dir.create(paste0(mainDir, './plots/', comparison, '/Gene_Ontology'), recursive = F)

#----Gene Ontology----
# Biological Pathways----
ont <- 'BP' # replace BP with MF for Molecular Functions or CC for Cell Components

d <- read.csv(paste0('./results/', comparison, '.csv'))
geneList <- d[ , 3] # extract Log2FoldChange
names(geneList) <- as.character(d[ , 10]) # extract entrez ids
geneList <- na.omit(geneList)
geneList <- geneList[unique(names(geneList))]
geneList <- sort(geneList, decreasing = T)

gene <- names(geneList)[abs(geneList) > 2]
gene <- na.omit(gene)
gene <- gene[!duplicated(gene)]

setwd(paste0('./plots/', comparison, '/Gene_Ontology'))

# Over-representation analysis----
ego <- enrichGO(gene = gene, universe = names(geneList),
                OrgDb = org.Hs.eg.db, ont = ont, readable = T,
                pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                pAdjustMethod = 'BH', minGSSize = 10, maxGSSize = 500)

filter(ego, p.adjust < 0.05, qvalue < 0.05)
mutate(ego, geneRatio = parse_ratio(GeneRatio)) %>% arrange(desc(geneRatio))
select(ego, -geneID) %>% head(n = 30)
ego_for_plot <- mutate(ego, richFactor = Count / as.numeric(sub('/\\d+', '', BgRatio)))

write.csv(ego_for_plot@result, file = paste0(mainDir, './results/enrichGO.', ont, '.', comparison, '.csv'))

pdf(paste0('./enrichGO.', ont, '.', comparison, '.pdf'))
dotplot(ego_for_plot, showCategory = 20, font.size = 6, orderBy = 'x') +
  xlab('Enrichment factor') + ylab(NULL) +
  ggtitle('Enriched Gene Ontology') + theme(aspect.ratio = 2) +
  scale_color_gradient(high = 'whitesmoke', low = 'indianred2')
dev.off()

# Gene Set Enrichment Analysis----
gse <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, ont = ont,
             minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0)

gse_for_plot <- arrange(gse, abs(NES)) %>% group_by(sign(NES))

write.csv(gse_for_plot@result, file = paste0(mainDir, './results/gseGO.', ont, '.', comparison, '.csv'))

pdf(paste0('./gseGO.', ont, '.', comparison, '.pdf'))
dotplot(gse_for_plot, showCategory = 10, font.size = 6, orderBy = 'x', split = '.sign') +
  facet_grid(.~.sign) + ylab(NULL) +
  ggtitle('Gene Set Enrichment') + theme(aspect.ratio = 2) +
  scale_color_gradient(high = 'whitesmoke', low = 'indianred2')
dev.off()

# Separate GSEA plots, if necessary----
figure <- list()
index <- 1:20

gsea <- function(x) {
  figure <- gseaplot2(gse, geneSetID = x, title = gse@result$Description[[x]])
}

figures <- lapply(index, gsea)

for (i in 1:length(figures)) {
  pdf(paste0('./gsea.', ont, '.', comparison, '_', gse@result$Description[[i]], '.pdf'),
      width = 4, height = 4.5); print(figures[i]); dev.off()
}

setwd(mainDir)

#----session----
sessionInfo()
