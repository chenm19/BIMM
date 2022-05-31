
#' ------------------------------------------------
#'  REACTOME enrichment analysis 
#' ------------------------------------------------
library(clusterProfiler); library(ReactomePA)

genes = as.character(na.omit(unique(c(edges[,1], edges[,2]))))
enrichs <- enrichPathway(gene=genes, pvalueCutoff = 0.05, readable=TRUE)
reactom <- as.data.frame(enrichs$reactome)
temp <- factor(reactom$ID, levels = unique(reactom$ID))
reactom$Name <- temp

library(scales); cols <- viridis_pal()(10)
#' show_col(viridis_pal()(10))
ggplot(reactom, aes(x=Name,y=score,fill=score))+
    geom_bar(stat="identity", width=0.39)+
                 scale_x_discrete(name ='Name') +
                           ggtitle('Reactome enrichment')







