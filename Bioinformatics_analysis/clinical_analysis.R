
library(survival); library(readxl)

read_xlsx <- function(filename) {
    sheets <- readxl::excel_sheets(filename)
    x <-    lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
    names(x) <- sheets; return (x) }

#' @param  working path
#' --------------------------
#'  RNA-seq data organize   |
#' --------------------------

#' @param input LUAD rnaseq data
rnaseq <- read.delim(file='TCGA_firehose_LUAD_rnaseqv2_RSEM_normalized_data.txt',stringsAsFactors=F,row.names=1,check.names=F,skip=1)
rnaseq.colname <- read.delim(file='TCGA_firehose_LUAD_rnaseqv2_RSEM_normalized_data.txt',stringsAsFactors=F,row.names=1,check.names=F,skip=0)
colnames(rnaseq) <- colnames(rnaseq.colname)

Names <- sapply(rownames(rnaseq),function(i){strsplit(i,'[|]')[[1]][1]})
rnaseq$name <- Names
rnaseq1 <- rnaseq[-which(duplicated(Names)),]
rownames(rnaseq1) <- Names[-which(duplicated(Names))]
# rnaseq1 <- rnaseq2

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

dat.surv <- dd.surv4

cc.stages <- dat.surv$stage
cc.stages[which(cc.stages%in%c('stage i','stage ia','stage ib'))] <- 'Stage-I'
cc.stages[which(cc.stages%in%c('stage ii','stage iia','stage iib'))] <- 'Stage-II'

cc.stages[-grep('-',cc.stages)] <- 'Stage-II+'
dat.surv$stage <- cc.stages

data.surv <- dat.surv

#' --------------------------------------------
#'  survival data and RNA-seq data in common  |
#' --------------------------------------------

temp <- as.character(sapply(colnames(rnaseq1),function(i){
    return (paste0(strsplit(i,'-')[[1]][1:4],collapse='-'))}))
colnames(rnaseq1) <- temp

#' @param ovalapping patients with both clinical data and RNAseq data
com.ss=intersect(rownames(data.surv),colnames(rnaseq1))         

ps3 <- as.numeric(na.omit(match(com.ss,rownames(data.surv))))
ps4 <- as.numeric(na.omit(match(com.ss,colnames(rnaseq1))))

com.surv <- data.surv[ps3,]; com.seq <- rnaseq1[,ps4]
identical(rownames(com.surv),colnames(com.seq))

#' @param com.surv survival data with ovalap samples; including tumor and adjacent normal
#' @param com.seq RNA-seq data with overlap samples

#' ---------------------------------------------------------
#'  add the category of early or late stage information    |
#' ---------------------------------------------------------
tmp1 <- which(com.surv$stage=='Stage-II+')

Com.surv <- com.surv; Com.surv$stage <- 'early'
Com.surv$stage[tmp1] <- 'late'

#' @param Com.surv survival data with ovalap samples; with stage information

all(identical(colnames(com.seq),rownames(Com.surv))); #' check
all(identical(colnames(com.seq),rownames(com.surv))); #' check

tem1 <- rownames(Com.surv)[Com.surv$type=='Primary Tumor']
tumor.seq <- com.seq[,tem1]; tumor.surv <- Com.surv[tem1,]

tem2 <- rownames(Com.surv)[Com.surv$type=='Solid Tissue Normal']
norm.seq <- com.seq[,tem2]; norm.surv <- Com.surv[tem2,]

#' @param tumor.seq tumor tissue with RNA-seq data
#' @param tumor.surv tumor tissue with clinical data 

#' ------------------------------
#'  tumor vs normal; survival   |
#' ------------------------------

edge <- read.csv('../tf_files_edges.csv',stringsAsFactors=F,row.names=1)
node <- read.csv(file='../tf_files_nodes.csv',stringsAsFactors=F,row.names=1)

edge[,1] <- node[match(edge[,1],node[,1]),2]
edge1 <- edge[-which(is.na(edge[,1])),]

#' for (i in 1:10){
#' new.gs <- c(edge1[which(edge1[,2]==sigs[i]),1],edge1[which(edge1[,1]==sigs[i]),2])
new.gs <- sigs
library(robustbase)
temp <- colMedians(as.matrix(tumor.seq[rownames(tumor.seq)%in%new.gs,]))

ft1 <- which(temp > quantile(temp,probs=0.7))
ft2 <- which(temp < quantile(temp,probs=0.3))

foo1 <- tumor.surv[ft1,]; foo1$category <- 'high'
foo2 <- tumor.surv[ft2,]; foo2$category <- 'low'
foot <- rbind(foo1,foo2)

library(survminer)
print (i)
t <- format(as.numeric(summary(coxph(Surv(time,status)~category,
                                     data=foot))$logtest[3]),digits=2)

print (t)
foot1 <- foot[foot$time<3000,]
fit <- survfit(Surv(time,status)~category,data=foot)

pdf('surv_sigs.pdf',width=6.9,height=6.8)
ggsurvplot(fit,
          pval = TRUE, #conf.int = TRUE,
          risk.table = TRUE, # Add risk table
          risk.table.col = "strata", # Change risk table color by groups
          linetype = "strata", # Change line type by groups
          surv.median.line = "hv", # Specify median survival
          ggtheme = theme_bw(), # Change ggplot2 theme
          palette = c("#E41A1C","#377EB8","#4DAF4A"))
graphics.off()

library(RColorBrewer)
brewer.pal(3,'Set1')

#' ------ roc curves ------
temp <- as.matrix(tumor.seq[rownames(tumor.seq)%in%new.gs,])
#; tumor.surv$value <- temp

temps <- cbind(tumor.surv,t(temp))

temps$val1 <- 1.0047*(1.321e-04*temps$EGF +
              # -1.775e-05*temps$ERBB2 +
              # -3.040e-06*temps$FOXO3
              + 4.515e-05*temps$MAPK1 +
              1.320e-04*temps$MYC + #2.137e-04*temps$SMAD3 +
              1.720e-04*temps$TGFB2) #+ 1.235e-05*temps$TGFBR1

coxph(Surv(time,status)~val1+stage,data=temps)
      
cutoff <- 365/3
Mayo4.2 = survivalROC(Stime=temps$time,  
                     status=temps$status,
                     marker = temps$val1,
                     predict.time = cutoff, method="KM")

pdf('auc_curve_3.pdf',heigh=3.6,width=3.6)
plot(Mayo4.2$FP, Mayo4.2$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "", "AUC = ",round(Mayo4.2$AUC,3)), 
     ylab="TP",main="Method = KM,  Year = 1", col ="#E7B800") # , "#2E9FDF"))
abline(0,1)
graphics.off()

cutoff <- 365/2
pdf('auc_curve_2.pdf',heigh=3.6,width=3.6)
plot(Mayo4.2$FP, Mayo4.2$TP, type="l", xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "", "AUC = ",round(Mayo4.2$AUC,3)), 
     ylab="TP",main="Method = KM,  Year = 1", col ="#2E9FDF")
abline(0,1)
graphics.off()

## Define a function
fun_survivalROC <- function(lp, t) {
    res <- with(foot,
                survivalROC(Stime        = time,
                            status       = status,
                            marker       = get(category),
                            predict.time = t,
                            method       = "KM"))       # KM method without smoothing

    ## Plot ROCs
    with(res, plot(TP ~ FP, type = "l", main = sprintf("t = %.0f, AUC = %.2f", t, AUC)))
    abline(a = 0, b = 1, lty = 2)

    res
}

## 2 x 5 layout
layout(matrix(1:10, byrow = T, ncol = 5))

## Model with age and sex
res_survivalROC <- lapply(1:10 * 365.25, function(t) {
    fun_survivalROC(lp = "lp.age.sex", t)
})

#' -----------------------------------
#'  analysis of downregulated genes  | 
#' -----------------------------------

foos <- lapply(1:9,function(i){
    temp1 <- as.data.frame(as.numeric(tumor.seq[rownames(tumor.seq)==sigs[i],]))
    temp2 <- as.data.frame(as.numeric(norm.seq[rownames(norm.seq)==sigs[i],]))
    colnames(temp1) <- 'data'; colnames(temp2) <- 'data'
    temp1$type <- 'tumor'; temp2$type <- 'normal'
    if (sum(is.na(temp))<506){
        temps <- rbind(temp1,temp2)
        temps$sig <- sigs[i]
        temps$type <- factor(temps$type,levels=c('tumor','normal'))
    }
    return (temps)
})

#fools <- do.call(rbind,foos)
        
sigs <- c("TGFBR1","TGFB2","SMAD3","FOXO3",'EGF','ERBB2','MAPK1','MYC','TGFB2')

for ( i in 1:9){
    r1 <- summary(foos[[i]]$data)[1]
    r2 <- summary(foos[[i]]$data)[5]
    pdf(paste0('TvsN_plot_',sigs[i],'.pdf'),width=3.9,height=3.9)
    p <- ggplot(foos[[i]], aes(x=type, y=data, fill=type)) +
        geom_boxplot(outlier.shape=NA,
                     notch=TRUE)+ ylim(c(r1,r2))+
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

for ( i in 1:9){
    t1 <- foos[[i]]
    t1.1 <- t1[t1[,2]=='tumor',1]; t1.2 <- t1[t1[,2]=='normal',1]
    print (paste(sigs[i],'pval=',format(wilcox.test(t1.1,t1.2)$p.value,digits=2)))
}
