
args = commandArgs(trailingOnly=TRUE)
print (length(args))

if (length(args)==0){ DEG <- TRUE }
if (length(args)>=1){ DEG <- args[1] }
library(Seurat)

if (DEG == TRUE){
  load('data_A549.rda')
  sc.data <- as.matrix(data_A549[[1]]$RNA)
  labels <- data_A549$labels
  colnames(sc.data) <- paste0('cell',1:ncol(sc.data))
  sc.label <- labels
  rownames(sc.label) <- paste0('cell',1:ncol(sc.data))
  sc.label[,1] <- paste0(sc.label[,1],'h')
  colnames(sc.label) <- 'condition'

  seurat.rna <- CreateSeuratObject(counts = sc.data,assay = 'RNA',
                                 project = 'RNA',min.cells = 0,
                                 meta.data = sc.label, verbose=FALSE)

  Idents(seurat.rna) <- seurat.rna@meta.data$condition 
  seurat.rna <- NormalizeData(seurat.rna, normalization.method = "LogNormalize", scale.factor = 10000,
                              verbose=FALSE)
  seurat.rna <- FindVariableFeatures(seurat.rna, selection.method = "vst", nfeatures = 2000,
                                     verbose=FALSE)
  seurat.rna <- ScaleData(seurat.rna,features = rownames(seurat.rna),verbose=FALSE)
  seurat.rna <- RunPCA(seurat.rna, features = VariableFeatures(object = seurat.rna),verbose=FALSE)

  seurat.rna <- FindNeighbors(seurat.rna, dims = 1:30,verbose=FALSE)
  seurat.rna <- FindClusters(seurat.rna, resolution = 0.5,verbose=FALSE)
  seurat.rna <- RunUMAP(seurat.rna, dims = 1:25,verbose=FALSE)

  col1 = RColorBrewer::brewer.pal(8, "Set1")[1:3]
  seurat.rna@meta.data$time = factor(seurat.rna@meta.data$condition, levels=c('0h','1h','3h'))

  pdf('UMAP-embedding.pdf')
  print (DimPlot(seurat.rna, reduction = "umap",group.by = 'time',cols=col1))
  graphics.off()

  Idents(seurat.rna) <- seurat.rna@meta.data$time
  deg.0 <- FindMarkers(seurat.rna,ident.1='0h',verbose=FALSE)
  deg.1 <- FindMarkers(seurat.rna,ident.1='1h',verbose=FALSE)
  deg.3 <- FindMarkers(seurat.rna,ident.1='3h',verbose=FALSE)

  #; DEG list of differnt conditions
  deg0 = unique(rownames(deg.0[abs(deg.0[,2])>=0.5 & deg.0[,5]<0.05,]))
  deg1 = unique(rownames(deg.1[abs(deg.1[,2])>=0.5 & deg.1[,5]<0.05,]))
  deg3 = unique(rownames(deg.3[abs(deg.3[,2])>=0.5 & deg.3[,5]<0.05,]))

  library(xlsx)
  write.xlsx(deg0,file='A549_DEGs.xlsx',sheetName='condition_0h',
           col.names=TRUE, row.names=TRUE, append=TRUE)
  write.xlsx(deg1,file='A549_DEGs.xlsx',sheetName='condition_1h',
           col.names=TRUE, row.names=TRUE, append=TRUE)
  write.xlsx(deg3,file='A549_DEGs.xlsx',sheetName='condition_3h',
           col.names=TRUE, row.names=TRUE, append=TRUE)
}

if (length(args)>=2){ 
  TF <- args[2] 
  if (TF){ source('TF_analysis.R')}
}

if (length(args)>=3){ 
  enrich <- args[3] 
  if (enrich){ source('enrich_analysis.R') }
}

if (length(args)>=4){ 
  surv <- args[4] 
  if (surv) {source('clinical_analysis.R')}
}
print ('FINISH')