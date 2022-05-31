
library(survival); library(survivalROC)

#' --------------------------
#'  RNA-seq data organize   |
#' --------------------------
rnaseq <- read.delim(file='TCGA_firehose_LUAD_rnaseqv2_RSEM_normalized_data.txt',stringsAsFactors=F,row.names=1,check.names=F,skip=1)
rnaseq.colname <- read.delim(file='TCGA_firehose_LUAD_rnaseqv2_RSEM_normalized_data.txt',stringsAsFactors=F,row.names=1,check.names=F,skip=0)
colnames(rnaseq) <- colnames(rnaseq.colname)

Names <- sapply(rownames(rnaseq),function(i){strsplit(i,'[|]')[[1]][1]})
rnaseq$name <- Names
rnaseq1 <- rnaseq[-which(duplicated(Names)),]
rownames(rnaseq1) <- Names[-which(duplicated(Names))]

#' -------------------------
#'  Survival data organize |
#' -------------------------
dd.surv <- read.delim('GDC_LUAD_xena_phenotype.tsv',header=T,sep='\t',stringsAsFactors=F,check.names=F)
rownames(dd.surv) <- dd.surv[,1]

dd.surv1 <- dd.surv[,colnames(dd.surv)%in%c('vital_status.diagnoses','days_to_death.diagnoses','days_to_last_follow_up.diagnoses','tumor_stage.diagnoses','sample_type.samples')]
dd.surv2 <- subset(dd.surv1,select=c(4,1,2,3,5))
colnames(dd.surv2) <- c('status','time','followup','stage','type')

ps1 <- which(dd.surv2[,1]=='')
dd.surv3 <- dd.surv2[-ps1,]

dd.surv3[,1] <- as.numeric(as.factor(dd.surv3[,1]))-1
dd.surv3[which(dd.surv3$status==0),2] <- dd.surv3[which(dd.surv3$status==0),3]

head(dd.surv3)
ps2 <- which(is.na(dd.surv3[,2]))
dd.surv4 <- dd.surv3[-which(is.na(dd.surv3[,2])),]
data.surv <- dd.surv4

temp <- as.character(sapply(colnames(rnaseq1),function(i){
    return (paste0(strsplit(i,'-')[[1]][1:4],collapse='-'))}))
colnames(rnaseq1) <- temp

#' @param ovalapping patients with both clinical data and RNAseq data
com.ss=intersect(rownames(data.surv),colnames(rnaseq1))         

ps3 <- as.numeric(na.omit(match(com.ss,rownames(data.surv))))
ps4 <- as.numeric(na.omit(match(com.ss,colnames(rnaseq1))))

com.surv <- data.surv[ps3,]; com.seq <- rnaseq1[,ps4]
identical(rownames(com.surv),colnames(com.seq))

Com.surv <- com.surv; 
all(identical(colnames(com.seq),rownames(Com.surv))); #' check
all(identical(colnames(com.seq),rownames(com.surv))); #' check

tem1 <- rownames(Com.surv)[Com.surv$type=='Primary Tumor']
tumor.seq <- com.seq[,tem1]; tumor.surv <- Com.surv[tem1,]

tem2 <- rownames(Com.surv)[Com.surv$type=='Solid Tissue Normal']
norm.seq <- com.seq[,tem2]; norm.surv <- Com.surv[tem2,]

#' --------------------
#'  survival analysis 
#' --------------------

#' @param genes
library(robustbase)
temp <- colMedians(as.matrix(tumor.seq[rownames(tumor.seq)%in%genes,]))

ft1 <- which(temp > quantile(temp,probs=0.7))
ft2 <- which(temp < quantile(temp,probs=0.3))

foo1 <- tumor.surv[ft1,]; foo1$category <- 'high'
foo2 <- tumor.surv[ft2,]; foo2$category <- 'low'
foot <- rbind(foo1,foo2)

library(survminer)
fit <- survfit(Surv(time,status)~category,data=foot)

pdf('survival_plot.pdf')
ggsurvplot(fit,
          pval = TRUE, 
          risk.table = TRUE, 
          risk.table.col = "strata",
          linetype = "strata",
          surv.median.line = "hv", 
          ggtheme = theme_bw(),
          palette = c("#E41A1C","#377EB8","#4DAF4A"))
graphics.off()

#' ------------------------
#'  survROC for ROC curves
#' ------------------------
temp <- as.matrix(tumor.seq[rownames(tumor.seq)%in%genes,])
temps <- cbind(tumor.surv,t(temp))

cutoff <- 365/3
roc1 = survivalROC(Stime=temps$time,  
                     status=temps$status,
                     marker = temps$val1,
                     predict.time = cutoff, method="KM")

pdf('roc_curve.pdf',heigh=3.6,width=3.6)
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "", "AUC = ",round(roc1$AUC,3)), 
     ylab="TP",main="Method = KM,  Year = 1", col ="#E7B800") 
abline(0,1)
graphics.off()


