
#' ----------------  
#' TF analysis  
#' ----------------  
library(RcisTarget)

#' take deg0 as an example
#gene = deg0
#database = '../cisTarget_database/hg19-tss-centered-10kb-7species.mc9nr.feather'

data(motifAnnotations_hgnc)
motifRankings <- importRankings(database)
motifAnnotations_hgnc[(directAnnotation==TRUE),]
motifEnrichmentTable_wGenes <- cisTarget(gene, motifRankings,motifAnnot=motifAnnotations_hgnc)
motifs_AUC <- calcAUC(gene, motifRankings, nCores=6)
motifEnrichmentTable <- addMotifAnnotation(motifs_AUC,motifAnnot=motifAnnotations_hgnc)

temp <- as.character(as.data.frame(motifEnrichmentTable[,5])[,1])
try <- lapply(temp,function(x){
  try <- strsplit(x,'[ (]')[[1]][1]
  if (grepl(';',try)){ strsplit(try,';')[[1]][1] } else {try}})
motifEnrichmentTable$TFs <- unlist(try)

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





