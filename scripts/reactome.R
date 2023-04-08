#----Reactome Pathway Analysis----
library(ReactomePA)

dir.create(paste0(mainDir, './plots/', comparison, '/Reactome'), recursive = F)
setwd(paste0('./plots/', comparison, '/Reactome'))

# Reactome over-representation analysis----
x <- enrichPathway(gene = gene, pvalueCutoff = 0.05, readable = T, organism = 'human')

filter(x, p.adjust < 0.05, qvalue < 0.1)
mutate(x, geneRatio = parse_ratio(GeneRatio)) %>% arrange(desc(geneRatio))
select(x, -geneID) %>% head
y <- mutate(x, richFactor = Count / as.numeric(sub('/\\d+', '', BgRatio)))

write.csv(y@result, file = paste0(mainDir, './results/enrichPathway.', comparison, '.csv'))

pdf(paste0('./enrichPathway.', comparison, '.pdf'))
dotplot(y, showCategory = 20, font.size = 6, orderBy = 'x') +
  xlab('Enrichment factor') + ylab(NULL) +
  ggtitle('Enriched Reactome') + theme(aspect.ratio = 2) +
  scale_color_gradient(high = 'whitesmoke', low = 'indianred2')
dev.off()

# Reactome GSEA----
g <- gsePathway(geneList, pvalueCutoff = 0.2, pAdjustMethod = 'BH', eps = 0,
                organism = 'human')
t <- arrange(g, abs(NES)) %>% group_by(sign(NES))

write.csv(t@result, file = paste0(mainDir, './results/gsePathway.', comparison, '.csv'))

pdf(paste0('./gsePathway.', comparison, '.pdf'))
dotplot(t, showCategory = 10, font.size = 6, orderBy = 'x', split = '.sign') +
  facet_grid(.~.sign) + ylab(NULL) +
  ggtitle('Reactome GSEA') + theme(aspect.ratio = 2) +
  scale_color_gradient(high = 'whitesmoke', low = 'indianred2')
dev.off()

setwd(mainDir)

#----session----
sessionInfo()
