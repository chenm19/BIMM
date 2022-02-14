
degs <- read.csv(file='Up_and_Down_Genes3vs0h.csv',stringsAsFactors=F)
genes <- read.csv(file='A549_single_cell_data.csv',stringsAsFactors=F)[,1]

max(degs$adjPval)
ps <- which(abs(degs[,2])>0.5)

#' -----  TF analysis ------
library(RcisTarget)
gene <- degs[ps,1]

#' Select motif database to use (i.e. organism and distance around TSS)
data(motifAnnotations_hgnc)

motifRankings <- importRankings("./cisTarget_database/hg19-tss-centered-10kb-7species.mc9nr.feather")

motifAnnotations_hgnc[(directAnnotation==TRUE),]
# Motif enrichment analysis:
motifEnrichmentTable_wGenes <- cisTarget(gene, motifRankings,motifAnnot=motifAnnotations_hgnc)

# Calculate AUC
motifs_AUC <- calcAUC(gene, motifRankings, nCores=6)

# Select significant motifs, add TF annotation & format as table
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,motifAnnot=motifAnnotations_hgnc)

xs <- as.character(as.data.frame(motifEnrichmentTable[,5])[,1])
try <- lapply(xs,function(x){
                  try <- strsplit(x,'[ (]')[[1]][1]
                  if (grepl(';',try)){ strsplit(try,';')[[1]][1] } else {try}})

motifEnrichmentTable$TFs <- unlist(try)

#' @param motifEnrichmentTable; res.tf

res.tf <- motifEnrichmentTable
tf.list <- unique(na.omit(res.tf$TFs))

signifMotifNames <- motifEnrichmentTable$motif

incidenceMatrix <- getSignificantGenes(gene,
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000-20, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

#' construct the network edges
library(reshape2)
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")

library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))

#' get TF labels
motifs1 <- motifEnrichmentTable[match(motifs,motifEnrichmentTable$motif),]$TFs

motifEnrichmentTable[which(is.na(motifs1)),]

#' construct the network nodes
nodes <- data.frame(id=c(motifs, genes),   
      label=c(motifs1, genes),    
      title=c(motifs1, genes), # tooltip 
      color=c(rep("purple", length(motifs1)), rep("skyblue", length(genes))))

ps <- which(is.na(nodes[,2]))[-1]
nodes[1,2] <- 'ERK'; nodes[1,3] <- 'ERK'

temp <- sapply(c(1:78),function(i){
    fo <- strsplit(as.character(motifEnrichmentTable[i,6]),';')[[1]]
    fo1 <- gsub(' ','',fo)[-length(fo)]
    return (fo1)})
temp1 <- unique(setdiff(unlist(temp),unique(nodes[,2])))

nodes[ps,2] <- temp1[1:20]; nodes[ps,3] <- temp1[1:20]

#' ------- plot tf network --------
library("igraph"); 

library(ggsci); col1 <- c(pal_material('red')(10)[3],pal_material('pink')(10)[3], pal_material('deep-orange')(10)[3],pal_material('red')(10)[4])

#' -----------------
#'  simplify edges |
#' -----------------

identical(nodes$id[match(as.character(edges[,1]),nodes$id)],as.character(edges[,1]))
edges[,1] <- as.factor(nodes$label[match(as.character(edges[,1]),nodes$id)])
new.edges <- unique(edges)

eds <- lapply(unique(edges[,1]),function(ed){
    edge1 <- edges[edges[,1]==ed,]
    if (nrow(edge1)>15){edge2=edge1[1:15,]} else {edge2 = edge1}
    return (edge2)
})
eds1 = do.call(rbind,eds)
new.edge  <- eds1
library(igraph)


library(RColorBrewer)
n <- 337
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
cell.col=sample(col_vector, n, replace=T)

names(cell.col) <- unique(c(new.edge[,1],new.edge[,2]))

colnames(new.edge) <- c('cell_from','cell_to')

library(dplyr)

pdf('vis.pdf',height=8,width=8)
NetView(new.edge,col=cell.col,vertex.label.cex=0.1,vertex.size=1.5, arrow.width=0,edge.max.width=3,
        label=FALSE)
graphics.off()

#' ------ visualization by RedgeR -------
library(RedeR)
nodes1 <- nodes[,1:3]; colnames(nodes1) <- c('id','Symbol','nodeAlias')

net <- graph_from_data_frame(d=new.edges, vertices=nodes1, directed=T)
gt.me  <- subg(g=net, dat=nodes, refcol=1)
gt.me  <- att.setv(g=gt.me, from="Symbol", to="nodeAlias")
gt.me <- att.setv(gt.me, from="weight", to="nodeColor", breaks=seq(-1,1,0.2), pal=2)    
gt.me <- att.setv(gt.me, from="size", to="nodeSize", nquant=10, isrev=TRUE, xlim=c(0,20,1))

#' saveRDS(gt.me,file='Network_data.RDS')

gt.me <- readRDS(file='Network_data.RDS')
rdp <- RedPort()
calld(rdp)
addGraph(rdp,gt.me)

resetd(rdp)

N0 <- addGraph(rdp, gt.me, gcoord=c(70,55), gscale=80, isNest=TRUE, theme='tm1', gzoom=30)
## N1 <- addGraph(rdp, gt6, gcoord=c(20,70), gscale=50, isNest=TRUE, theme='tm1', gzoom=30)
## N2 <- addGraph(rdp, gt12, gcoord=c(70,55), gscale=80, isNest=TRUE, theme='tm1', gzoom=30)

N3 <- nestNodes(rdp, nodes=V(gt.me)$name, parent=N1, theme='tm2')

mergeOutEdges(rdp)

relax(rdp,50,400)

scl <- gt3$legNodeColor$scale
leg <- gt3$legNodeColor$legend 
addLegend.color(rdp, colvec=scl, labvec=leg, title="node color (logFC)")

selectNodes(rdp,"RET")

addGraph(rdp,sg)



