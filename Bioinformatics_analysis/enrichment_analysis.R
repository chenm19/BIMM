
edge <- read.csv('../tf_files_edges.csv',stringsAsFactors=F,row.names=1)
node <- read.csv(file='../tf_files_nodes.csv',stringsAsFactors=F,row.names=1)

edge[,1] <- node[match(edge[,1],node[,1]),2]
edge1 <- edge[-which(is.na(edge[,1])),]

new.gs = unique(c(edge1[,1], edge1[,2]))
    
#' ------------------------------------------------
#'  draw GOplot: notice "scale (-log); zscore"    |
#' ------------------------------------------------
library(GOplot)
go.enrich <- GO.enrich(geneset=new.gs)

try <- do.call(rbind,go.enrich)

try$name <- try$ID
try$category <- as.character(sapply(rownames(try),function(i){substr(i,start=1,stop=2)}))
rownames(try) <- NULL
colnames(try) <- c('ID','term','GeneRatio','BgRatio','pvalue','adj_pval','qvalue','genes','Count','score','size','name','category')

b <- data.frame(ID=new.gs,logFC=1,stringsAsFactors=F)
circ <- circle_dat(terms=try,genes=b)

circ$zscore <- as.numeric(sapply(try$genes,function(i){length(strsplit(i,'[/]')[[1]])})/sqrt(163))

circ$count <- try$Count
circ$logFC <- sample(seq(0,1,0.1),size=109,replace=T); circ$zscore <- (-log10(try$score))
reduced_circ <- reduce_overlap(circ, overlap = 0.85)

cols <- c("chartreuse4", "brown2", "cornflowerblue")

pdf('GO_enrichment.pdf',width=20,height=6.6)
print (GOBubble(reduced_circ,labels=10))
dev.off()

#' -----------------------------------

enrichs <- kegg.enrich(geneset=new.gs)

reactom <- as.data.frame(enrichs$reactome)

#' write.csv(reactom,file='REACTOME_enrich_results.csv',quote=F)
#' reactom  <- read.csv(file='REACTOME_enrich_results.csv',stringsAsFactors=F)

try <- factor(reactom$ID, levels = unique(reactom$ID))
reactom$Name <- try

library(scales); cols <- viridis_pal()(10)
#' show_col(viridis_pal()(10))

im1 <- ggplot(reactom, aes(x=Name,y=score,fill=score))+
    geom_bar(stat="identity", width=0.39)+
             scale_fill_gradient2(low=cols[1], mid=cols[5], high=cols[10],
                                  midpoint=as.numeric(summary(reactom$score)[3])) +
                 scale_x_discrete(name ='Name') + coord_flip()+ theme_bw(base_size=12)+
                     #' theme_Publication(base_family='Arial',base_size=7)+
                       theme(axis.text.x = element_text(color='black',size=12),
                             axis.text.y =element_text(color='black',size=12),axis.title=element_blank())+
                           ggtitle('Reactome enrichment')
                           #' geom_point(aes(y=Count),stat="identity",alpha=.5,size=1)


pdf('reactome_new.pdf',height=4.5,width=10); print (im1); dev.off()





