
####################################################
#Input Data
####################################################
library(Signac); library(Seurat)

sc.data <- as.matrix(x[[1]])
colnames(sc.data) <- paste0('cell',1:ncol(sc.data))
sc.label <- labels
rownames(sc.label) <- paste0('cell',1:ncol(sc.data))
sc.label[,1] <- paste0(sc.label[,1],'h')
colnames(sc.label) <- 'condition'
write.csv(sc.data,file='A549_single_cell_data.csv',quote=F)
write.csv(sc.label,file='A549_single_cell_label.csv',quote=F)

seurat.rna@meta.data$celltype <- paste0(seurat.rna@meta.data$Time,'h')
Idents(seurat.rna) <- seurat.rna@meta.data$celltype 

deg.0 <- FindMarkers(seurat.rna,ident.1='0h')
deg.1 <- FindMarkers(seurat.rna,ident.1='1h')
deg.3 <- FindMarkers(seurat.rna,ident.1='3h')

deg.03 <- FindMarkers(seurat.rna,ident.1='3h',ident.2='0h')

library(xlsx)

write.xlsx(deg.0,file='Up_and_Down_Regulated_Genes.xlsx',sheetName='condition_0h',
           col.names=TRUE, row.names=TRUE, append=TRUE)
write.xlsx(deg.1,file='Up_and_Down_Regulated_Genes.xlsx',sheetName='condition_1h',
           col.names=TRUE, row.names=TRUE, append=TRUE)
write.xlsx(deg.3,file='Up_and_Down_Regulated_Genes.xlsx',sheetName='condition_3h',
           col.names=TRUE, row.names=TRUE, append=TRUE)
write.xlsx(deg.03,file='Up_and_Down_Regulated_Genes.xlsx',sheetName='condition_3h_vs_0h',
           col.names=TRUE, row.names=TRUE, append=TRUE)

seurat.rna <- NormalizeData(seurat.rna, normalization.method = "LogNormalize", scale.factor = 10000)
seurat.rna <- FindVariableFeatures(seurat.rna, selection.method = "vst", nfeatures = 2000)

seurat.rna <- ScaleData(seurat.rna,features = rownames(seurat.rna))

seurat.rna <- RunPCA(seurat.rna, features = VariableFeatures(object = seurat.rna))

seurat.rna <- FindNeighbors(seurat.rna, dims = 1:10)
seurat.rna <- FindClusters(seurat.rna, resolution = 0.5)
seurat.rna <- RunUMAP(seurat.rna, dims = 1:10)

seurat.rna@meta.data$celltype <- paste0(seurat.rna@meta.data$Time,'h')
Idents(seurat.rna) <- seurat.rna@meta.data$celltype 

pdf('try1.pdf')
print (DoHeatmap(object = seurat.rna))
print (DimPlot(seurat.rna, reduction = "umap"))
graphics.off()

pdf('try2.pdf')
print (FeaturePlot(seurat.rna, features = c("TGFBR1","TGFB2","SMAD3","FOXO3")))
graphics.off()

pdf('try3.pdf')
VlnPlot(seurat.rna, features = c("TGFBR1","TGFB2","SMAD3","FOXO3"))
graphics.off()

temp <- as.matrix(seurat.rna@assays$RNA@scale.data)
temp <- as.matrix(seurat.rna@assays$RNA@counts)

ab <- temp[rownames(temp)%in%c("TGFBR1","TGFB2","SMAD3","FOXO3"),]

ps1 <- which(seurat.rna@meta.data$celltype=='0h')
ps2 <- which(seurat.rna@meta.data$celltype=='1h')
ps3 <- which(seurat.rna@meta.data$celltype=='3h')

foo <- data.frame(rbind(rowMeans(ab[,ps1]), rowMeans(ab[,ps2]), rowMeans(ab[,ps3])))
foo$time <- c('0h','1h','3h')

#' ------- plot heatmap -------
var.features = unique(c(rownames(deg.0[abs(deg.0[,2])>1,]),'TGFBR1','FOXO3','TGFB2','SMAD3',rownames(deg.3[abs(deg.3[,2])>0.7,])))
                #,c('BTBD11','AKAP13','ASPH','CNTN1','KYNU'))

rownames(deg.1[abs(deg.1[,2])>0.4,]),

pdf('heatmap_minghan.pdf')
DoHeatmap(object = seurat.rna, features = var.features)
graphics.off()

#' ------- plot barplot -------
library(reshape2); library(ggplot2); library(dplyr)
com.p <- melt(foo)
colnames(com.p) <- c('method','data','value')

tels1 <- lapply(c("TGFBR1","TGFB2","SMAD3","FOXO3"),function(d){
    tel <- com.p[com.p[,2]==d,]
    tel1 <- tel[match(c('0h','1h','3h'),tel$time),]
    tel1$time <- factor(tel1$time,levels=c('0h','1h','3h'))
    return (tel1) })

for ( i in 1:4){
    d <- c("TGFBR1","TGFB2","SMAD3","FOXO3")[i]
    pdf(paste0('boxplot1_',d,'.pdf'),width=6,height=3.9)
    p <- ggplot(tels1[[i]], aes(x=time, y=value, fill=time)) + 
        geom_bar(stat="identity", color="black", width=1,position=position_dodge()) +
        scale_fill_brewer(palette='Set1')+ theme_bw() +
        theme(axis.text.x = element_text(color='black',size=10,angle=60,vjust=0.8),
              axis.text.y =element_text(color='black',size=12),
              legend.text=element_text(size=14),
              text = element_text(size = 11),
              axis.title = element_text(face="bold"),
              panel.background = element_rect(colour = NA),
              plot.background = element_rect(colour = NA),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
    print (p)
    graphics.off()
}

