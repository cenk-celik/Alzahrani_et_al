#----KEGG analysis----
library(pathview)

dir.create(paste0(mainDir, './plots/', comparison, '/KEGG'), recursive = F)
setwd(paste0('./plots/', comparison, '/KEGG'))

gsek <- gseKEGG(geneList, organism = 'human', minGSSize = 10, pvalueCutoff = 0.05, eps = 0)
k <- arrange(gsek, abs(NES)) %>% group_by(sign(NES))

pdf(paste0('./gseKEGG.', comparison, '.pdf'))
dotplot(k, showCategory = 10, font.size = 6, orderBy = 'x', split = '.sign') +
  facet_grid(.~.sign) + ylab(NULL) +
  ggtitle('KEGG') + theme(aspect.ratio = 2) +
  scale_color_gradient(high = 'whitesmoke', low = 'indianred2')
dev.off()

setwd(mainDir)

#----session----
sessionInfo()