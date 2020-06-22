#Eddy October 2016
# application of edgeR models to identify differentially expressed transcripts over time and between host races
#   using RNAseq data

#base code for main tables

#source("https://bioconductor.org/biocLite.R")
#biocLite("edgeR")
#biocLite("locfit")

library(edgeR)
library(stringr)
library(locfit)
library(statmod)
library(ggplot2)
library(reshape)
library(heatmap3)

setwd("~/Documents/Cerasi/pomonella/RSEMresults/")

#######create table of counts#######

file_list<-list.files(pattern='*genes.results')
file_list
table<-c()

for (file in file_list){
  sampleID<-str_match(file,"(^[A-Za-z0-9]*_[A-Za-z0-9]*_[A-Za-z0-9]*_[0-9]*)")[,1]
  print(sampleID)
  df<-read.table(file,header=T,sep="\t",row.names=NULL)
  df2<-df[,c(1,5)]
  colnames(df2)[colnames(df2)=="expected_count"] <- paste(sampleID,"_exp_count",sep="")
  if (length(table) ==0){table<-df2
  } else {
  table<-merge(table,df2,by="gene_id") }
}  


#annotation file (must set quotes)
#annotation file has been reduced down to gene_id hence the duplication in places
annotationFile<-read.table("ref_Rhagoletis_zephyria_1.0_top_level.gtf.unique.info.mergedgeneid",sep="\t",quote = "",header=TRUE,row.names=NULL)

table<-merge(annotationFile,table, by="gene_id",all=TRUE)

#flybase ids from top hit blastn
#with gg codes flybase !!!this is reduced to unique gg codes per gene (thus have lost unique pp codes per gene)!!!
#is top hit over e-6
blastn.flybase<-read.table("../../GenomeFiles/Rhag_zephyria_transcripts.blastxFlybase.tab.onehit_evaluee6_geneid_singlehit",sep="\t",quote = "",header=TRUE,row.names=NULL)
table<-merge(blastn.flybase,table, by="gene_id",all=TRUE)

#flybase ids from blastp top 6 considered those e-value e-6
blastp.flybase<-read.table("Rhag_zephyria_tblastnProteinfa.tab.fivehit.columnsuseful.evaluee6.reduced",sep="\t",header=TRUE,row.names=NULL,quote = "")
blastp.flybase$protein_id<-NULL
table<-merge(blastp.flybase,table, by="gene_id",all=TRUE)

#swissport annotations
# e-6
blastp.swissport<-read.table("rhago_transcriptome_STARextract_swissport_e6_reduce_geneid",sep="\t",header=TRUE,row.names=NULL,quote = "")
table<-merge(blastp.swissport,table,  by="gene_id",all=TRUE)


#########greg's filtering function represented in >= 50% of samples by at least one count#######
filterMinCount<- function(x) {
  pres<-x >=1
  out=F
  if ((sum(pres)/length(pres)) >= 0.5) {out=T}
  return(out)
}

sum(table$Rp_A_2M_01_exp_count)# = 15304315 (24536873)
sum(table$Rp_A_2M_02_exp_count)# = 30540959 (50950016)
sum(table$Rp_A_2M_02_exp_count)
colnames(table)
#remove two weirdo's
table$Rp_H_4M_04_exp_count <- NULL
table$Rp_H_4M_02_exp_count <- NULL

#with flybase
filterInd.flybase<-apply(table[,(-1:-29)],1,filterMinCount) #you will need to change this to represent the annotations you have brought in
table <- table[filterInd.flybase,]
#17274 passing filter

########create DGElist object#######
colnames(table)
#two duds out flybase
table.dge<-DGEList(counts=table[,30:77],genes=table[,(1:29)])
table.dge$samples
names(table)

#normalisation
?calcNormFactors
table.dge<- calcNormFactors(table.dge) #default is TMM
table.dge$samples

#########basic MDS plot ########
plotMDS(table.dge)
#coloured/labeled figure:
plotMDS(table.dge, top=500,pch = c(rep(15,5),rep(2,5),rep(4,5),rep(19,5),rep(7,5),rep(15,5),rep(2,5),rep(4,3),rep(19,5),rep(7,5)),col=c(rep("blue",25), rep("red",23)))
legend(0.8,1.3,bty = "n",legend=(c("2 mo.","3 mo.","4 mo.","5 mo.","6 mo.","Apple","Haw")),pch = c(15,2,4,19,7,18,18),col=c("black","black","black","black","black","blue","red"))

plotMDS(table.dge, top=1000,pch = c(rep(15,5),rep(2,5),rep(4,5),rep(19,5),rep(7,5),rep(15,5),rep(2,5),rep(4,3),rep(19,5),rep(7,5)),col=c(rep("blue",25), rep("red",23)))
legend(0.9,0.1,bty = "n",legend=(c("2 mo.","3 mo.","4 mo.","5 mo.","6 mo.","Apple","Haw")),pch = c(15,2,4,19,7,18,18),col=c("black","black","black","black","black","blue","red"))



########pomonella design 1########
#month * fruit
month<-factor(c(rep(2,5),rep(3,5),rep(4,5),rep(5,5),rep(6,5),rep(2,5),rep(3,5),rep(4,3),rep(5,5),rep(6,5)))

fruit<-factor(c(rep("Apple",25), rep("Haw",23)))

design <- model.matrix(~month*fruit)
rownames(design) <- colnames(table.dge)
design

#full model
#                         (Intercept)  month3 month4 month5 month6    fruitApple        month3:fruitApple month4:fruitApple month5:fruitApple month6:fruitApple
#Rp_A_2M_01_exp_count           1      0      0      0      0          1                 0                 0                 0                 0
#Rp_A_2M_02_exp_count           1      0      0      0      0          1                 0                 0                 0                 0
#Rp_A_2M_03_exp_count           1      0      0      0      0          1                 0                 0                 0                 0
#Rp_A_2M_04_exp_count           1      0      0      0      0          1                 0                 0                 0                 0
#Rp_A_2M_05_exp_count           1      0      0      0      0          1                 0                 0                 0                 0
#Rp_A_3M_01_exp_count           1      1      0      0      0          1                 1                 0                 0                 0
#Rp_A_3M_02_exp_count           1      1      0      0      0          1                 1                 0                 0                 0
#Rp_A_3M_03_exp_count           1      1      0      0      0          1                 1                 0                 0                 0
#Rp_A_3M_04_exp_count           1      1      0      0      0          1                 1                 0                 0                 0
#Rp_A_3M_05_exp_count           1      1      0      0      0          1                 1                 0                 0                 0
#Rp_A_4M_01_exp_count           1      0      1      0      0          1                 0                 1                 0                 0
#Rp_A_4M_02_exp_count           1      0      1      0      0          1                 0                 1                 0                 0
#Rp_A_4M_03_exp_count           1      0      1      0      0          1                 0                 1                 0                 0
#Rp_A_4M_04_exp_count           1      0      1      0      0          1                 0                 1                 0                 0
#Rp_A_4M_05_exp_count           1      0      1      0      0          1                 0                 1                 0                 0
#Rp_A_5M_01_exp_count           1      0      0      1      0          1                 0                 0                 1                 0
#Rp_A_5M_02_exp_count           1      0      0      1      0          1                 0                 0                 1                 0
#Rp_A_5M_03_exp_count           1      0      0      1      0          1                 0                 0                 1                 0
#Rp_A_5M_04_exp_count           1      0      0      1      0          1                 0                 0                 1                 0
#Rp_A_5M_05_exp_count           1      0      0      1      0          1                 0                 0                 1                 0
#Rp_A_6M_01_exp_count           1      0      0      0      1          1                 0                 0                 0                 1
#Rp_A_6M_02_exp_count           1      0      0      0      1          1                 0                 0                 0                 1
#Rp_A_6M_03_exp_count           1      0      0      0      1          1                 0                 0                 0                 1
#Rp_A_6M_04_exp_count           1      0      0      0      1          1                 0                 0                 0                 1
#Rp_A_6M_05_exp_count           1      0      0      0      1          1                 0                 0                 0                 1
#Rp_H_2M_01_exp_count           1      0      0      0      0          0                 0                 0                 0                 0
#Rp_H_2M_02_exp_count           1      0      0      0      0          0                 0                 0                 0                 0
#Rp_H_2M_03_exp_count           1      0      0      0      0          0                 0                 0                 0                 0
#Rp_H_2M_04_exp_count           1      0      0      0      0          0                 0                 0                 0                 0
#Rp_H_2M_05_exp_count           1      0      0      0      0          0                 0                 0                 0                 0
#Rp_H_3M_01_exp_count           1      1      0      0      0          0                 0                 0                 0                 0
#Rp_H_3M_02_exp_count           1      1      0      0      0          0                 0                 0                 0                 0
#Rp_H_3M_03_exp_count           1      1      0      0      0          0                 0                 0                 0                 0
#Rp_H_3M_04_exp_count           1      1      0      0      0          0                 0                 0                 0                 0
#Rp_H_3M_05_exp_count           1      1      0      0      0          0                 0                 0                 0                 0
#Rp_H_4M_01_exp_count           1      0      1      0      0          0                 0                 0                 0                 0
#Rp_H_4M_03_exp_count           1      0      1      0      0          0                 0                 0                 0                 0
#Rp_H_4M_05_exp_count           1      0      1      0      0          0                 0                 0                 0                 0
#Rp_H_5M_01_exp_count           1      0      0      1      0          0                 0                 0                 0                 0
#Rp_H_5M_02_exp_count           1      0      0      1      0          0                 0                 0                 0                 0
#Rp_H_5M_03_exp_count           1      0      0      1      0          0                 0                 0                 0                 0
#Rp_H_5M_04_exp_count           1      0      0      1      0          0                 0                 0                 0                 0
#Rp_H_5M_05_exp_count           1      0      0      1      0          0                 0                 0                 0                 0
#Rp_H_6M_01_exp_count           1      0      0      0      1          0                 0                 0                 0                 0
#Rp_H_6M_02_exp_count           1      0      0      0      1          0                 0                 0                 0                 0
#Rp_H_6M_03_exp_count           1      0      0      0      1          0                 0                 0                 0                 0
#Rp_H_6M_04_exp_count           1      0      0      0      1          0                 0                 0                 0                 0
#Rp_H_6M_05_exp_count           1      0      0      0      1          0                 0                 0                 0                 0

#########estimate disperson squareroot=coefficient of variation of biological variation#######
table.dge<- estimateDisp(table.dge, design, robust=TRUE)
table.dge$common.dispersion
#0.08347977
sqrt(0.06750284)
#0.26 coefficent of variation

#plot
plotBCV(table.dge)

##############contrasts##################
fit <- glmFit(table.dge, design)

colnames(design)
#[1] "(Intercept)"       "month3"            "month4"            "month5"            "month6"            "fruitApple"        "month3:fruitApple" "month4:fruitApple"
#[9] "month5:fruitApple" "month6:fruitApple"

#########WITHIN APPLE########
#apple three to two months
lrt3monthvs2Apple <- glmLRT(fit,contrast=c(0,1,0,0,0,0,1,0,0,0))
Apple3monthvs2 <- data.frame(topTags(lrt3monthvs2Apple,n=nrow(table),sort="none"))
#create table for 2M->3M->4M->5M->6M
FDR.table.apple<-Apple3monthvs2[,c(1,34)]
colnames(FDR.table.apple)[colnames(FDR.table.apple)=="FDR"] <- "month3vs2_apple_FDR"
LogFC.table.apple<-Apple3monthvs2[,c(1,30)]
colnames(LogFC.table.apple)[colnames(LogFC.table.apple)=="logFC"] <- "month3vs2_apple_logFC"

#2nd table for greg November 2016
FDR.table.from2M.greg.apple<-Apple3monthvs2[,c(1,16,34)]
colnames(FDR.table.from2M.greg.apple)[colnames(FDR.table.from2M.greg.apple)=="FDR"] <- "month3vs2_apple_FDR"
LogFC.table.from2M.greg.apple<-Apple3monthvs2[,c(1,16,30)]
colnames(LogFC.table.from2M.greg.apple)[colnames(LogFC.table.from2M.greg.apple)=="logFC"] <- "month3vs2_apple_logFC"

#seperate table everything from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’
FDR.table.from2M.apple<-Apple3monthvs2[,c(1,34)]
colnames(FDR.table.from2M.apple)[colnames(FDR.table.from2M.apple)=="FDR"] <- "month3vs2_apple_FDR"
LogFC.table.from2M.apple<-Apple3monthvs2[,c(1,30)]
colnames(LogFC.table.from2M.apple)[colnames(LogFC.table.from2M.apple)=="logFC"] <- "month3vs2_apple_logFC"

#apple four to two months
lrt4monthvs2Apple <- glmLRT(fit,contrast=c(0,0,1,0,0,0,0,1,0,0))
Apple4monthvs2 <- data.frame(topTags(lrt4monthvs2Apple,n=nrow(table),sort="none"))
Apple4monthvs2.FDR05 <- subset(Apple4monthvs2, FDR < 0.05)
FDR.table.from2M.apple$month4vs2_apple_FDR<-Apple4monthvs2[,34]
LogFC.table.from2M.apple$month4vs2_apple_logFC<-Apple4monthvs2[,30]

#table for greg
FDR.table.from2M.greg.apple$month4vs2_apple_FDR<-Apple4monthvs2[,34]
LogFC.table.from2M.greg.apple$month4vs2_apple_logFC<-Apple4monthvs2[,30]

#apple five to two months
lrt5monthvs2Apple <- glmLRT(fit,contrast=c(0,0,0,1,0,0,0,0,1,0))
Apple5monthvs2 <- data.frame(topTags(lrt5monthvs2Apple,n=nrow(table),sort="none"))
Apple5monthvs2.FDR05 <- subset(Apple5monthvs2, FDR < 0.05)
FDR.table.from2M.apple$month5vs2_apple_FDR<-Apple5monthvs2[,34]
LogFC.table.from2M.apple$month5vs2_apple_logFC<-Apple5monthvs2[,30]

#table for greg
FDR.table.from2M.greg.apple$month5vs2_apple_FDR<-Apple5monthvs2[,34]
LogFC.table.from2M.greg.apple$month5vs2_apple_logFC<-Apple5monthvs2[,30]

#apple six to two months
lrt6monthvs2Apple <- glmLRT(fit,contrast=c(0,0,0,0,1,0,0,0,0,1))
Apple6monthvs2 <- data.frame(topTags(lrt6monthvs2Apple,n=nrow(table),sort="none"))
Apple6monthvs2.FDR05 <- subset(Apple6monthvs2, FDR < 0.05)
FDR.table.from2M.apple$month6vs2_apple_FDR<-Apple6monthvs2[,34]
LogFC.table.from2M.apple$month6vs2_apple_logFC<-Apple6monthvs2[,30]

#table for greg
FDR.table.from2M.greg.apple$month6vs2_apple_FDR<-Apple6monthvs2[,34]
LogFC.table.from2M.greg.apple$month6vs2_apple_logFC<-Apple6monthvs2[,30]

#apple four to three months
lrt4monthvs3Apple <- glmLRT(fit,contrast=c(0,-1,1,0,0,0,-1,1,0,0))
Apple4monthvs3 <- data.frame(topTags(lrt4monthvs3Apple,n=nrow(table),sort="none"))
Apple4monthvs3.FDR05 <- subset(Apple4monthvs3, FDR < 0.05)
#
compare row names same or not
all.equal(Apple3monthvs2$gene_id,Apple4monthvs3$gene_id)
test<-Apple3monthvs2$gene_id==Apple4monthvs3$gene_id
length(test[test==TRUE])
length(test[test==FALSE])
FDR.table.apple$month4vs3_apple_FDR<-Apple4monthvs3[,34]
LogFC.table.apple$month4vs3_apple_logFC<-Apple4monthvs3[,30]

#apple five to four months
lrt5monthvs4Apple <- glmLRT(fit,contrast=c(0,0,-1,1,0,0,0,-1,1,0))
Apple5monthvs4 <- data.frame(topTags(lrt5monthvs4Apple,n=nrow(table),sort="none"))
Apple5monthvs4.FDR05 <- subset(Apple5monthvs4, FDR < 0.05)
FDR.table.apple$month5vs4_apple_FDR<-Apple5monthvs4[,34]
LogFC.table.apple$month5vs4_apple_logFC<-Apple5monthvs4[,30]

#apple six to five months
lrt6monthvs5Apple <- glmLRT(fit,contrast=c(0,0,0,-1,1,0,0,0,-1,1))
Apple6monthvs5 <- data.frame(topTags(lrt6monthvs5Apple,n=nrow(table),sort="none"))
Apple6monthvs5.FDR05 <- subset(Apple6monthvs5, FDR < 0.05)
FDR.table.apple$month6vs5_apple_FDR<-Apple6monthvs5[,34]
LogFC.table.apple$month6vs5_apple_logFC<-Apple6monthvs5[,30]

LogFC.table.apple %>% filter(.,gene_id=='gene17502')

#######APPLE VS HAW##########
#two months haw vs apple
lrt2monthApplevsHaw <- glmLRT(fit, contrast = c(0,0,0,0,0,1,0,0,0,0))
ApplevsHaw2Month <- data.frame(topTags(lrt2monthApplevsHaw,n=nrow(table),sort="none"))
ApplevsHaw2Month.FDR05 <- subset(ApplevsHaw2Month, FDR < 0.05)
##create table for fruit month interaction
FDR.table.fruitInteraction<-ApplevsHaw2Month[,c(1,34)]
colnames(FDR.table.fruitInteraction)[colnames(FDR.table.fruitInteraction)=="FDR"] <- "ApplevsHaw2Month_FDR"
LogFC.table.fruitInteraction<-ApplevsHaw2Month[,c(1,30)]
colnames(LogFC.table.fruitInteraction)[colnames(LogFC.table.fruitInteraction)=="logFC"] <- "ApplevsHaw2Month_LogFC"

#three months apple vs haw
lrt3monthApplevsHaw <- glmLRT(fit, contrast = c(0,0,0,0,0,1,1,0,0,0))
ApplevsHaw3Month <- data.frame(topTags(lrt3monthApplevsHaw,n=nrow(table),sort="none"))
ApplevsHaw3Month.FDR05 <- subset(ApplevsHaw3Month, FDR < 0.05)
FDR.table.fruitInteraction$ApplevsHaw3Month_FDR<-ApplevsHaw3Month[,34]
LogFC.table.fruitInteraction$ApplevsHaw3Month_LogFC<-ApplevsHaw3Month[,30]

#four months apple vs haw
lrt4monthApplevsHaw <- glmLRT(fit, contrast = c(0,0,0,0,0,1,0,1,0,0))
ApplevsHaw4Month <- data.frame(topTags(lrt4monthApplevsHaw,n=nrow(table),sort="none"))
ApplevsHaw4Month.FDR05 <- subset(ApplevsHaw4Month, FDR < 0.05)
FDR.table.fruitInteraction$ApplevsHaw4Month_FDR<-ApplevsHaw4Month[,34]
LogFC.table.fruitInteraction$ApplevsHaw4Month_LogFC<-ApplevsHaw4Month[,30]

#five months apple vs haw
lrt5monthApplevsHaw <- glmLRT(fit, contrast = c(0,0,0,0,0,1,0,0,1,0))
ApplevsHaw5Month <- data.frame(topTags(lrt5monthApplevsHaw,n=nrow(table),sort="none"))
ApplevsHaw5Month.FDR05 <- subset(ApplevsHaw5Month, FDR < 0.05)
FDR.table.fruitInteraction$ApplevsHaw5Month_FDR<-ApplevsHaw5Month[,34]
LogFC.table.fruitInteraction$ApplevsHaw5Month_LogFC<-ApplevsHaw5Month[,30]

#six months apple vs haw
lrt6monthApplevsHaw <- glmLRT(fit, contrast = c(0,0,0,0,0,1,0,0,0,1))
ApplevsHaw6Month <- data.frame(topTags(lrt6monthApplevsHaw,n=nrow(table),sort="none"))
ApplevsHaw6Month.FDR05 <- subset(ApplevsHaw6Month, FDR < 0.05)
FDR.table.fruitInteraction$ApplevsHaw6Month_FDR<-ApplevsHaw6Month[,34]
LogFC.table.fruitInteraction$ApplevsHaw6Month_LogFC<-ApplevsHaw6Month[,30]

#apple 5M vs Haw 6M
lrt5Mhaw6Mapple <- glmLRT(fit,contrast=c(0,0,0,-1,1,1,0,0,0,1))
Haw5MvsApple6Month <- data.frame(topTags(lrt5Mhaw6Mapple,n=nrow(table),sort="none"))
Haw5MvsApple6Month.FDR05 <- subset(Haw5MvsApple6Month, FDR < 0.05)
length(Haw5MvsApple6Month.FDR05$gene_id) #518


#######WITHIN HAW########
#haw 3M ->2M
lrt3monthvs2Haw <- glmLRT(fit, contrast = c(0,1,0,0,0,0,0,0,0,0))
Haw3monthvs2 <- data.frame(topTags(lrt3monthvs2Haw,n=nrow(table),sort="none"))
Haw3monthvs2.FDR05 <- subset(Haw3monthvs2, FDR < 0.05)

#create table for 2M->3M->4M->5M->6M
FDR.table.haw<-Haw3monthvs2[,c(1,34)]
colnames(FDR.table.haw)[colnames(FDR.table.haw)=="FDR"] <- "month3vs2_haw_FDR"
LogFC.table.haw<-Haw3monthvs2[,c(1,30)]
colnames(LogFC.table.haw)[colnames(LogFC.table.haw)=="logFC"] <- "month3vs2_haw_logFC"

#seperate table everything from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’
FDR.table.from2M.haw<-Haw3monthvs2[,c(1,34)]
colnames(FDR.table.from2M.haw)[colnames(FDR.table.from2M.haw)=="FDR"] <- "month3vs2_haw_FDR"
LogFC.table.from2M.haw<-Haw3monthvs2[,c(1,30)]
colnames(LogFC.table.from2M.haw)[colnames(LogFC.table.from2M.haw)=="logFC"] <- "month3vs2_haw_logFC"

#haw 4M ->2M
lrt4monthvs2Haw <- glmLRT(fit, contrast = c(0,0,1,0,0,0,0,0,0,0))
Haw4monthvs2 <- data.frame(topTags(lrt4monthvs2Haw,n=nrow(table),sort="none"))
Haw4monthvs2.FDR05 <- subset(Haw4monthvs2, FDR < 0.05)
FDR.table.from2M.haw$month4vs2_haw_FDR<-Haw4monthvs2[,34]
LogFC.table.from2M.haw$month4vs2_haw_logFC<-Haw4monthvs2[,30]

#haw 5M ->2M
lrt5monthvs2Haw <- glmLRT(fit, contrast = c(0,0,0,1,0,0,0,0,0,0))
Haw5monthvs2 <- data.frame(topTags(lrt5monthvs2Haw,n=nrow(table),sort="none"))
Haw5monthvs2.FDR05 <- subset(Haw5monthvs2, FDR < 0.05)
FDR.table.from2M.haw$month5vs2_haw_FDR<-Haw5monthvs2[,34]
LogFC.table.from2M.haw$month5vs2_haw_logFC<-Haw5monthvs2[,30]

#haw 6M ->2M
lrt6monthvs2Haw <- glmLRT(fit, contrast = c(0,0,0,0,1,0,0,0,0,0))
Haw6monthvs2 <- data.frame(topTags(lrt6monthvs2Haw,n=nrow(table),sort="none"))
Haw6monthvs2.FDR05 <- subset(Haw6monthvs2, FDR < 0.05)
FDR.table.from2M.haw$month6vs2_haw_FDR<-Haw6monthvs2[,34]
LogFC.table.from2M.haw$month6vs2_haw_logFC<-Haw6monthvs2[,30]

#haw 4M->3M
lrt4monthvs3Haw <- glmLRT(fit, contrast = c(0,-1,1,0,0,0,0,0,0,0))
Haw4monthvs3 <- data.frame(topTags(lrt4monthvs3Haw,n=nrow(table),sort="none"))
Haw4monthvs3.FDR05 <- subset(Haw4monthvs3, FDR < 0.05)
FDR.table.haw$month4vs3_haw_FDR<-Haw4monthvs3[,34]
LogFC.table.haw$month4vs3_haw_logFC<-Haw4monthvs3[,30]

#haw 5M->4M
lrt5monthvs4Haw <- glmLRT(fit, contrast = c(0,0,-1,1,0,0,0,0,0,0))
Haw5monthvs4 <- data.frame(topTags(lrt5monthvs4Haw,n=nrow(table),sort="none"))
Haw5monthvs4.FDR05 <- subset(Haw5monthvs4, FDR < 0.05)
FDR.table.haw$month5vs4_haw_FDR<-Haw5monthvs4[,34]
LogFC.table.haw$month5vs4_haw_logFC<-Haw5monthvs4[,30]

#haw 6M->5M
lrt6monthvs5Haw <- glmLRT(fit, contrast = c(0,0,0,-1,1,0,0,0,0,0))
Haw6monthvs5 <- data.frame(topTags(lrt6monthvs5Haw,n=nrow(table),sort="none"))
Haw6monthvs5.FDR05 <- subset(Haw6monthvs5, FDR < 0.05)
FDR.table.haw$month6vs5_haw_FDR<-Haw6monthvs5[,34]
LogFC.table.haw$month6vs5_haw_logFC<-Haw6monthvs5[,30]

#Apple from baseline of Haw 2M#
#apple 2M -> 2M (haw)
lrt2monthvs2AppleToHaw2M <- glmLRT(fit, contrast = c(0,0,0,0,0,1,0,0,0,0))
Apple2monthvs2MHaw <- data.frame(topTags(lrt2monthvs2AppleToHaw2M,n=nrow(table),sort="none"))
Apple2monthvs2MHaw.FDR05 <- subset(Apple2monthvs2MHaw, FDR < 0.05)
#seperate table everything from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’ using apple baseline
FDR.table.from2Mhaw.apple<-Apple2monthvs2MHaw[,c(1,34)]
colnames(FDR.table.from2Mhaw.apple)[colnames(FDR.table.from2Mhaw.apple)=="FDR"] <- "Apple2monthvs2Haw_FDR"
LogFC.table.from2Mhaw.apple<-Apple2monthvs2MHaw[,c(1,30)]
colnames(LogFC.table.from2Mhaw.apple)[colnames(LogFC.table.from2Mhaw.apple)=="logFC"] <- "Apple2monthvs2Haw_logFC"

#table for greg November 2016
FDR.table.from2Mhaw.greg.apple<-Apple2monthvs2MHaw[,c(1,16,34)]
colnames(FDR.table.from2Mhaw.greg.apple)[colnames(FDR.table.from2Mhaw.greg.apple)=="FDR"] <- "Apple2monthvs2Haw_FDR"
LogFC.table.from2Mhaw.greg.apple<-Apple2monthvs2MHaw[,c(1,16,30)]
colnames(LogFC.table.from2Mhaw.greg.apple)[colnames(LogFC.table.from2Mhaw.greg.apple)=="logFC"] <- "Apple2monthvs2Haw_logFC"

#apple 3M ->2M (haw)
lrt3monthvs2AppleToHaw2M <- glmLRT(fit, contrast = c(0,1,0,0,0,1,1,0,0,0))
Apple3monthvs2MHaw <- data.frame(topTags(lrt3monthvs2AppleToHaw2M,n=nrow(table),sort="none"))
Apple3monthvs2MHaw.FDR05 <- subset(Apple3monthvs2MHaw, FDR < 0.05)
FDR.table.from2Mhaw.apple$Apple3monthvs2Haw_FDR<-Apple3monthvs2MHaw[,34]
LogFC.table.from2Mhaw.apple$Apple3monthvs2Haw_logFC<-Apple3monthvs2MHaw[,30]

FDR.table.from2Mhaw.greg.apple$Apple3monthvs2Haw_FDR<-Apple3monthvs2MHaw[,34]
LogFC.table.from2Mhaw.greg.apple$Apple3monthvs2Haw_logFC<-Apple3monthvs2MHaw[,30]

#apple 4M ->2M (haw)
lrt4monthvs2AppleToHaw2M <- glmLRT(fit, contrast = c(0,0,1,0,0,1,0,1,0,0))
Apple4monthvs2MHaw <- data.frame(topTags(lrt4monthvs2AppleToHaw2M,n=nrow(table),sort="none"))
Apple4monthvs2MHaw.FDR05 <- subset(Apple4monthvs2MHaw, FDR < 0.05)
FDR.table.from2Mhaw.apple$Apple4monthvs2Haw_FDR<-Apple4monthvs2MHaw[,34]
LogFC.table.from2Mhaw.apple$Apple4monthvs2Haw_logFC<-Apple4monthvs2MHaw[,30]

FDR.table.from2Mhaw.greg.apple$Apple4monthvs2Haw_FDR<-Apple4monthvs2MHaw[,34]
LogFC.table.from2Mhaw.greg.apple$Apple4monthvs2Haw_logFC<-Apple4monthvs2MHaw[,30]

#apple 5M ->2M (haw)
lrt5monthvs2AppleToHaw2M <- glmLRT(fit, contrast = c(0,0,0,1,0,1,0,0,1,0))
Apple5monthvs2MHaw <- data.frame(topTags(lrt5monthvs2AppleToHaw2M,n=nrow(table),sort="none"))
Apple5monthvs2MHaw.FDR05 <- subset(Apple5monthvs2MHaw, FDR < 0.05)
FDR.table.from2Mhaw.apple$Apple5monthvs2Haw_FDR<-Apple5monthvs2MHaw[,34]
LogFC.table.from2Mhaw.apple$Apple5monthvs2Haw_logFC<-Apple5monthvs2MHaw[,30]

FDR.table.from2Mhaw.greg.apple$Apple5monthvs2Haw_FDR<-Apple5monthvs2MHaw[,34]
LogFC.table.from2Mhaw.greg.apple$Apple5monthvs2Haw_logFC<-Apple5monthvs2MHaw[,30]

#apple 6M ->2M (haw)
lrt6monthvs2AppleToHaw2M <- glmLRT(fit, contrast = c(0,0,0,0,1,1,0,0,0,1))
Apple6monthvs2MHaw <- data.frame(topTags(lrt6monthvs2AppleToHaw2M,n=nrow(table),sort="none"))
Apple6monthvs2MHaw.FDR05 <- subset(Apple6monthvs2MHaw, FDR < 0.05)
FDR.table.from2Mhaw.apple$Apple6monthvs2Haw_FDR<-Apple6monthvs2MHaw[,34]
LogFC.table.from2Mhaw.apple$Apple6monthvs2Haw_logFC<-Apple6monthvs2MHaw[,30]

FDR.table.from2Mhaw.greg.apple$Apple6monthvs2Haw_FDR<-Apple6monthvs2MHaw[,34]
LogFC.table.from2Mhaw.greg.apple$Apple6monthvs2Haw_logFC<-Apple6monthvs2MHaw[,30]

#######fruit time interaction#########
colnames(design)
fruitTimeinteraction<-glmLRT(fit, coef=7:10)
fruitTimeinteraction.tags <- data.frame(topTags(fruitTimeinteraction,n=nrow(table),sort="none"))
fruitTimeinteraction.tags.FDR05 <- subset(fruitTimeinteraction.tags, FDR < 0.05)
length(fruitTimeinteraction.tags.FDR05$gene_id) #625

##########month################
Month<-glmLRT(fit, coef=2:5)
Month.tags <- data.frame(topTags(Month,n=nrow(table),sort="none"))
Month.tags.FDR05 <- subset(Month.tags, FDR < 0.05)
length(Month.tags.FDR05$gene_id) #3785

########fruit###############
Fruit<-glmLRT(fit, coef=6)
Fruit.tags <- data.frame(topTags(Fruit,n=nrow(table),sort="none"))
Fruit.tags.FDR05 <- subset(Fruit.tags, FDR < 0.05)
length(Fruit.tags.FDR05$gene_id) #1134

##intersect fruit:time month##
fruittimeinteraction.gene_id<-fruitTimeinteraction.tags.FDR05$gene_id
length(fruittimeinteraction.gene_id) #625

Month.tags.FDR05.withFruitTimeinteraction <- Month.tags.FDR05[Month.tags.FDR05$gene_id %in% fruittimeinteraction.gene_id, ]
length(Month.tags.FDR05.withFruitTimeinteraction$gene_id) #351 

Month.tags.FDR05.NOFruitTimeinteraction <- Month.tags.FDR05[!Month.tags.FDR05$gene_id %in% fruittimeinteraction.gene_id, ]
length(Month.tags.FDR05.NOFruitTimeinteraction$gene_id) #3434 

Fruit.tags.FDR05.withFruitTimeinteraction <- Fruit.tags.FDR05[Fruit.tags.FDR05$gene_id %in% fruittimeinteraction.gene_id, ]
length(Fruit.tags.FDR05.withFruitTimeinteraction$gene_id) #242 

Fruit.tags.FDR05.NOFruitTimeinteraction <- Fruit.tags.FDR05[!Fruit.tags.FDR05$gene_id %in% fruittimeinteraction.gene_id, ]
length(Fruit.tags.FDR05.NOFruitTimeinteraction$gene_id) #892 


##############time series 2M->3M->4M->5M->6M###############
#graph all genes
library(ggplot2)

#APPLE
#graph all genes with a significant interaction
FDR.table.apple.sigrows<-FDR.table.apple[apply(FDR.table.apple[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid<-FDR.table.apple.sigrows$gene_id
LogFC.table.sigFDR.apple <- LogFC.table.apple[LogFC.table.apple$gene_id %in% sig_geneid, ]

library(reshape)
df <- melt(LogFC.table.sigFDR.apple, "gene_id")
ggplot(df, aes(variable, value,group = gene_id)) +
  geom_line()

#do significant from 2M->3M->4M->5M->6M but log values from from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’

LogFC.table.from2M.sigFDR.apple <- LogFC.table.from2M.apple[LogFC.table.from2M.apple$gene_id %in% sig_geneid, ]

df_from2M <- melt(LogFC.table.from2M.sigFDR.apple, "gene_id")
ggplot(df_from2M, aes(variable, value,group = gene_id)) +
  geom_line()

#significant from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’

FDR.table.from2M.apple.sigrows<-FDR.table.from2M.apple[apply(FDR.table.from2M.apple[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid_from2M<-FDR.table.from2M.apple.sigrows$gene_id
LogFC.table.from2M.sigFDR.apple.from2M <- LogFC.table.from2M.apple[LogFC.table.from2M.apple$gene_id %in% sig_geneid_from2M, ]

library(reshape)
df_from2M <- melt(LogFC.table.from2M.sigFDR.apple.from2M, "gene_id")
ggplot(df_from2M, aes(variable, value,group = gene_id)) +
  geom_line()

#HAW
FDR.table.haw.sigrows<-FDR.table.haw[apply(FDR.table.haw[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid<-FDR.table.haw.sigrows$gene_id
LogFC.table.sigFDR.haw <- LogFC.table.haw[LogFC.table.haw$gene_id %in% sig_geneid, ]

library(reshape)
df <- melt(LogFC.table.sigFDR.haw, "gene_id")
ggplot(df, aes(variable, value,group = gene_id)) +
  geom_line()

#do significant from 2M->3M->4M->5M->6M but log values from from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’

LogFC.table.from2M.sigFDR.haw <- LogFC.table.from2M.haw[LogFC.table.from2M.haw$gene_id %in% sig_geneid, ]

df_from2M <- melt(LogFC.table.from2M.sigFDR.haw, "gene_id")
ggplot(df_from2M, aes(variable, value,group = gene_id)) +
  geom_line()

#significant from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’

FDR.table.from2M.haw.sigrows<-FDR.table.from2M.haw[apply(FDR.table.from2M.haw[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid_from2M<-FDR.table.from2M.haw.sigrows$gene_id
LogFC.table.from2M.sigFDR.haw.from2M <- LogFC.table.from2M.haw[LogFC.table.from2M.haw$gene_id %in% sig_geneid_from2M, ]

library(reshape)
df_from2M <- melt(LogFC.table.from2M.sigFDR.haw.from2M, "gene_id")
ggplot(df_from2M, aes(variable, value,group = gene_id)) +
  geom_line()

#significant from 2M ‘2M->3M’ 2M->4M’ ‘2M->5M’ ‘2M->6M’ haw 2m baseline

FDR.table.from2Mhaw.apple.sigrows<-FDR.table.from2Mhaw.apple[apply(FDR.table.from2Mhaw.apple[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig_geneid_from2Mhaw.apple<-FDR.table.from2Mhaw.apple.sigrows$gene_id
LogFC.table.from2M.sigFDR.apple.from2Mhaw <- LogFC.table.from2Mhaw.apple[LogFC.table.from2Mhaw.apple$gene_id %in% sig_geneid_from2Mhaw.apple, ]

library(reshape)
df_from2M <- melt(LogFC.table.from2M.sigFDR.apple.from2Mhaw, "gene_id")
ggplot(df_from2M, aes(variable, value,group = gene_id)) +
  geom_line()

#how many siginificant genes are shared?#

test<-LogFC.table.sigFDR.apple$gene_id %in% LogFC.table.sigFDR.haw$gene_id
sum(test)
#185 shared genes of significance in 2M->3M->4M->5M->6M
length(LogFC.table.sigFDR.apple$gene_id) #2083
length(LogFC.table.sigFDR.haw$gene_id) #398

test2<-LogFC.table.from2M.sigFDR.apple.from2M$gene_id %in% LogFC.table.from2M.sigFDR.haw.from2M$gene_id
sum(test2)
#2578 shared genes of significance in 2M->3M 2M->4M 2M->5M 2M->6M
length(LogFC.table.from2M.sigFDR.apple.from2M$gene_id) #2578
length(LogFC.table.from2M.sigFDR.haw.from2M$gene_id) #5133

#using 2M haw baseline for apple
test6<-LogFC.table.from2M.sigFDR.apple.from2Mhaw$gene_id %in% LogFC.table.from2M.sigFDR.haw.from2M$gene_id
sum(test6)
#2696 shared genes of significance in 2M->3M 2M->4M 2M->5M 2M->6M
length(LogFC.table.from2M.sigFDR.apple.from2Mhaw$gene_id) #4648
length(LogFC.table.from2M.sigFDR.haw.from2M$gene_id) #4597
colnames(LogFC.table.from2M.sigFDR.haw.from2M)

#genes in series compared to genes sig from 2M
test3<-LogFC.table.from2M.sigFDR.apple.from2M$gene_id %in% LogFC.table.sigFDR.apple$gene_id
sum(test3) #1644 shared out of 2083

#genes in series compared to genes sig from 2M
test4<-LogFC.table.from2M.sigFDR.haw.from2M$gene_id %in% LogFC.table.sigFDR.haw$gene_id
sum(test4) #278 out of a posible 398

#output top flybase hit for gene ontology analysis
LogFC.table.from2M.sigFDR.apple.from2M.flybase<-table[table$gene_id %in% LogFC.table.from2M.sigFDR.apple.from2M$gene_id,]
head(LogFC.table.from2M.sigFDR.apple.from2M.flybase)
length(LogFC.table.from2M.sigFDR.apple.from2M.flybase$gene_id)

LogFC.table.from2M.sigFDR.apple.from2M.flybase<-LogFC.table.from2M.sigFDR.apple.from2M.flybase[c(1,16)]
head(LogFC.table.from2M.sigFDR.apple.from2M.flybase)
length(LogFC.table.from2M.sigFDR.apple.from2M.flybase$gene_id)

LogFC.table.from2M.sigFDR.haw.from2M.flybase<-table[table$gene_id %in% LogFC.table.from2M.sigFDR.haw.from2M$gene_id,]
head(LogFC.table.from2M.sigFDR.haw.from2M.flybase)
length(LogFC.table.from2M.sigFDR.haw.from2M.flybase$gene_id)

LogFC.table.from2M.sigFDR.haw.from2M.flybase<-LogFC.table.from2M.sigFDR.haw.from2M.flybase[c(1,16)]
head(LogFC.table.from2M.sigFDR.haw.from2M.flybase)
length(LogFC.table.from2M.sigFDR.haw.from2M.flybase$gene_id)

LogFC.table.from2M.sigFDR.apple.from2Mhaw.flybase<-table[table$gene_id %in% LogFC.table.from2M.sigFDR.apple.from2Mhaw$gene_id,]
head(LogFC.table.from2M.sigFDR.apple.from2Mhaw.flybase)
length(LogFC.table.from2M.sigFDR.apple.from2Mhaw.flybase$gene_id)

LogFC.table.from2M.sigFDR.apple.from2Mhaw.flybase<-LogFC.table.from2M.sigFDR.apple.from2Mhaw.flybase[c(1,16)]
head(LogFC.table.from2M.sigFDR.apple.from2Mhaw.flybase)
length(LogFC.table.from2M.sigFDR.apple.from2Mhaw.flybase$gene_id)

LogFC.table.from2M.haw.flybase<-table[c(1,16)]
head(LogFC.table.from2M.haw.flybase)
length(LogFC.table.from2M.haw.flybase$gene_id)
LogFC.table.from2M.apple.flybase<-LogFC.table.from2M.haw.flybase

###################fruit time interactions#########################
length(fruitTimeinteraction.tags.FDR05$gene_id) #625
fruittimeinteraction.gene_id<-fruitTimeinteraction.tags.FDR05$gene_id

test<-LogFC.table.from2M.sigFDR.apple.from2M$gene_id %in% fruitTimeinteraction.tags.FDR05$gene_id
sum(test) #488

test<-LogFC.table.from2M.sigFDR.haw.from2M$gene_id %in% fruitTimeinteraction.tags.FDR05$gene_id
sum(test) #342

test<-LogFC.table.from2M.sigFDR.apple.from2Mhaw$gene_id %in% fruitTimeinteraction.tags.FDR05$gene_id
sum(test) #468

#genes with time 
#genes in apple/haw time series from 2M that have/have not fruit:time interaction
#apple
length(LogFC.table.from2M.sigFDR.apple.from2M$gene_id) #5133
LogFC.table.from2M.sigFDR.apple.from2M.withFruitTimeinteraction <- LogFC.table.from2M.sigFDR.apple.from2M[LogFC.table.from2M.sigFDR.apple.from2M$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.from2M.sigFDR.apple.from2M.withFruitTimeinteraction$gene_id) #488

LogFC.table.from2M.sigFDR.apple.from2M.NOFruitTimeinteraction <- LogFC.table.from2M.sigFDR.apple.from2M[!LogFC.table.from2M.sigFDR.apple.from2M$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.from2M.sigFDR.apple.from2M.NOFruitTimeinteraction$gene_id) #4645

#haw
length(LogFC.table.from2M.sigFDR.haw.from2M$gene_id) #4597
LogFC.table.from2M.sigFDR.haw.from2M.withFruitTimeinteraction <- LogFC.table.from2M.sigFDR.haw.from2M[LogFC.table.from2M.sigFDR.haw.from2M$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.from2M.sigFDR.haw.from2M.withFruitTimeinteraction$gene_id) #342

LogFC.table.from2M.sigFDR.haw.from2M.NOFruitTimeinteraction <- LogFC.table.from2M.sigFDR.haw.from2M[!LogFC.table.from2M.sigFDR.haw.from2M$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.from2M.sigFDR.haw.from2M.NOFruitTimeinteraction$gene_id) #4255

#apple from haw 2M
length(LogFC.table.from2M.sigFDR.apple.from2Mhaw$gene_id) #4648
LogFC.table.from2M.sigFDR.apple.from2Mhaw.withFruitTimeinteraction <- LogFC.table.from2M.sigFDR.apple.from2Mhaw[LogFC.table.from2M.sigFDR.apple.from2Mhaw$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.from2M.sigFDR.apple.from2Mhaw.withFruitTimeinteraction$gene_id) #468

LogFC.table.from2M.sigFDR.apple.from2Mhaw.NOFruitTimeinteraction <- LogFC.table.from2M.sigFDR.apple.from2Mhaw[!LogFC.table.from2M.sigFDR.apple.from2Mhaw$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.from2M.sigFDR.apple.from2Mhaw.NOFruitTimeinteraction$gene_id) #4180

#genes with fruit 
# genes in 2Mapple->2Mhaw, 3Mapple->3Mhaw, 4Mapple->4Mhaw, 5Mapple->5Mhaw, 6Mapple->6Mhaw that are/are not in fruit:month interactions
length(fruittimeinteraction.gene_id) #625

FDR.table.fruitInteraction.sigrows<-FDR.table.fruitInteraction[apply(FDR.table.fruitInteraction[, -1], MARGIN = 1, function(x) any(x < 0.05)), ]
sig.fruitInteraction<-FDR.table.fruitInteraction.sigrows$gene_id
length(sig.fruitInteraction) #2902
LogFC.table.fruitInteraction.sigFDR <- LogFC.table.fruitInteraction[LogFC.table.fruitInteraction$gene_id %in% sig.fruitInteraction, ]
length(LogFC.table.fruitInteraction.sigFDR$gene_id) #2902

LogFC.table.fruitInteraction.withFruitTimeinteraction<-LogFC.table.fruitInteraction.sigFDR[LogFC.table.fruitInteraction.sigFDR$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.fruitInteraction.withFruitTimeinteraction$gene_id) #563

LogFC.table.fruitInteraction.NOFruitTimeinteraction<-LogFC.table.fruitInteraction.sigFDR[!LogFC.table.fruitInteraction.sigFDR$gene_id %in% fruittimeinteraction.gene_id, ]
length(LogFC.table.fruitInteraction.NOFruitTimeinteraction$gene_id) #2339

#figures of within haw/apple with fruit:time interactions#
#apple time series from apple2M with fruit:time interaction
head(LogFC.table.from2M.sigFDR.apple.from2Mhaw.withFruitTimeinteraction)
library(ggplot2)
library(reshape)
df <- melt(LogFC.table.from2M.sigFDR.apple.from2Mhaw.withFruitTimeinteraction, "gene_id")
ggplot(df, aes(variable, value,group = gene_id)) +
  geom_line()

#haw time series from apple2M with fruit:time interaction
head(LogFC.table.from2M.sigFDR.haw.from2M.withFruitTimeinteraction)
df <- melt(LogFC.table.from2M.sigFDR.haw.from2M.withFruitTimeinteraction, "gene_id")
ggplot(df, aes(variable, value,group = gene_id)) +
  geom_line()#+

###########generating table for greg################
head(LogFC.table.from2M.greg.apple)
head(LogFC.table.from2M.haw)

head(FDR.table.from2Mhaw.greg.apple)
head(FDR.table.from2M.haw)

FDR.table.from2M.greg.apple$Lowest_FDR_from2M.TimeSeries.apple<-apply(FDR.table.from2M.greg.apple[,3:6],1,min)
FDR.table.from2Mhaw.greg.apple$Lowest_FDR_from2Mhaw.TimeSeries.apple<-apply(FDR.table.from2Mhaw.greg.apple[,3:7],1,min)
FDR.table.from2M.haw$Lowest_FDR_from.TimeSeries.haw<-apply(FDR.table.from2M.haw[,2:5],1,min)

Table.from2M.greg.apple<-merge(LogFC.table.from2M.greg.apple,FDR.table.from2M.greg.apple,by="gene_id")
Table.from2Mhaw.greg.apple<-merge(LogFC.table.from2Mhaw.greg.apple,FDR.table.from2Mhaw.greg.apple,by="gene_id")
Table.from2M.haw<-merge(LogFC.table.from2M.haw,FDR.table.from2M.haw,by="gene_id")

Table.from2M.apple<-merge(Table.from2M.greg.apple,Table.from2Mhaw.greg.apple,by="gene_id")

colnames(Table.from2M.greg.apple)
colnames(Table.from2Mhaw.greg.apple)
colnames(Table.from2M.haw)
colnames(Table.from2M.apple)
Table.from2M.apple$FlyBase_FBgn_tophit.y.x<-NULL
Table.from2M.apple$FlyBase_FBgn_tophit.y.y<-NULL
Table.from2M.apple$FlyBase_FBgn_tophit.x.y<-NULL

fruitTimeinteraction.greg<-fruitTimeinteraction.tags[c(1,37)]
colnames(fruitTimeinteraction.greg)[colnames(fruitTimeinteraction.greg)=="FDR" ]<- "FDR_FruitTimeInteraction"
Table.from2M.apple<-merge(Table.from2M.apple,fruitTimeinteraction.greg,by="gene_id")
Table.from2M.haw<-merge(Table.from2M.haw,fruitTimeinteraction.greg,by="gene_id")

Fruit.greg<-Fruit.tags[c(1,34)]
colnames(Fruit.greg)[colnames(Fruit.greg)=="FDR" ]<- "FDR_Fruit"
Table.from2M.apple<-merge(Table.from2M.apple,Fruit.greg,by="gene_id")
Table.from2M.haw<-merge(Table.from2M.haw,Fruit.greg,by="gene_id")

Month.greg<-Month.tags[c(1,37)]
colnames(Month.greg)[colnames(Month.greg)=="FDR" ]<- "FDR_Month"
Table.from2M.apple<-merge(Table.from2M.apple,Month.greg,by="gene_id")
Table.from2M.haw<-merge(Table.from2M.haw,Month.greg,by="gene_id")

colnames(Table.from2M.haw)
colnames(Table.from2M.apple)

#merging in 2->3->4->5->6

colnames(FDR.table.apple)
FDR.table.apple$Lowest_FDR_consTimeSeries<-apply(FDR.table.apple[,2:5],1,min)
timeseries.cons.apple<-merge(LogFC.table.apple,FDR.table.apple,by="gene_id")
colnames(FDR.table.haw)
FDR.table.haw$Lowest_FDR_consTimeSeries<-apply(FDR.table.haw[,2:5],1,min)
timeseries.cons.haw<-merge(LogFC.table.haw,FDR.table.haw,by="gene_id")

Table.from2M.apple<-merge(Table.from2M.apple, timeseries.cons.apple,by="gene_id")
Table.from2M.haw<-merge(Table.from2M.haw,timeseries.cons.haw,by="gene_id")

#with apple->haw 2M, apple->haw 3M, apple->haw 4M, apple->haw 5M, apple->haw 6M in 
head(FDR.table.fruitInteraction)
head(LogFC.table.fruitInteraction)

colnames(FDR.table.fruitInteraction)
FDR.table.fruitInteraction$Lowest_FDR_FruitBetweenMonths<-apply(FDR.table.fruitInteraction[,2:6],1,min)
apple.haw.month.contrasts<-merge(LogFC.table.fruitInteraction,FDR.table.fruitInteraction,by='gene_id')

Table.from2M.apple<-merge(Table.from2M.apple, apple.haw.month.contrasts,by="gene_id")
Table.from2M.haw<-merge(Table.from2M.haw,apple.haw.month.contrasts,by="gene_id")

colnames(Table.from2M.haw)
colnames(Table.from2M.apple)

flybase_annot<-Table.from2M.apple[c(1:2)]
head(flybase_annot)

Table.from2M.haw<-merge(Table.from2M.haw,flybase_annot,by="gene_id")

write.table(Table.from2M.apple,"AppleBackToHaw2MremoveBadHaw4M/Table.apple.JuneEJD",quote=F,sep="\t",row.names=FALSE)
write.table(Table.from2M.haw,"AppleBackToHaw2MremoveBadHaw4M/Table.haw.JuneEJD",quote=F,sep="\t",row.names=FALSE)

colnames(Table.from2M.apple)
colnames(Table.from2M.haw)
#with haw haw2M->3M haw2M->4M haw2M->5M haw2>->6M

head(LogFC.table.from2M.haw)
head(FDR.table.from2M.haw)

colnames(FDR.table.from2M.haw)
FDR.table.from2M.haw$Lowest_FDR_timeSeries2Mhaw<-apply(FDR.table.from2M.haw[,2:5],1,min)

haw.fromHAW2M<-merge(LogFC.table.from2M.haw,FDR.table.from2M.haw,by="gene_id")

Table.from2Mapple.greg.haw.cons.monthcontrasts.haw2M<-merge(Table.from2Mapple.greg.haw.cons.monthcontrasts,haw.fromHAW2M,by="gene_id")
write.table(Table.from2Mapple.greg.haw.cons.monthcontrasts.haw2M,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mapple.greg.haw.cons.monthcontrasts.haw2M",quote=F,sep="\t",row.names=FALSE)





#insulin tor wnt signalling
insulin.flybase<-read.table("insulinSignalingFlybase.txt",sep="\t",quote = "",header=FALSE,row.names=NULL)
colnames(insulin.flybase)[colnames(insulin.flybase)=="V1" ]<- "FlyBase_FBgn_tophit.x"

tor.flybase<-read.table("torSignalingFlybase.txt",sep="\t",quote = "",header=FALSE,row.names=NULL)
colnames(tor.flybase)[colnames(tor.flybase)=="V1" ]<- "FlyBase_FBgn_tophit.x"

wnt.flybase<-read.table("wntSignalingFlybase.txt",sep="\t",quote = "",header=FALSE,row.names=NULL,fill=TRUE)
colnames(wnt.flybase)[colnames(wnt.flybase)=="V1" ]<- "FlyBase_FBgn_tophit.x"

Table.from2Mapple.greg.haw$host<-rep("haw",nrow(Table.from2Mapple.greg.haw))
Table.from2M.greg.apple$host<-rep("apple",nrow(Table.from2M.greg.apple))
Table.from2M.greg.apple$M2_LogFC<-rep(0,nrow(Table.from2M.greg.apple))
Table.from2M.greg.apple$M2_FDR<-rep("NA",nrow(Table.from2M.greg.apple))

colnames(Table.from2Mapple.greg.haw)
Table.from2Mapple.greg.haw<-Table.from2Mapple.greg.haw[,c(17,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)]
Table.from2M.greg.apple<-Table.from2M.greg.apple[,c(15,1,2,16,3,4,5,6,17,7,8,9,10,11,12,13,14)]


Table.from2M.greg.apple.insulin <- Table.from2M.greg.apple[grep(paste(insulin.flybase$FlyBase_FBgn_tophit.x,collapse="|"), Table.from2M.greg.apple$FlyBase_FBgn_tophit.x),]
Table.from2Mapple.greg.haw.insulin <- Table.from2Mapple.greg.haw[grep(paste(insulin.flybase$FlyBase_FBgn_tophit.x,collapse="|"), Table.from2Mapple.greg.haw$FlyBase_FBgn_tophit.x),]

Table.from2M.greg.apple.tor <- Table.from2M.greg.apple[grep(paste(tor.flybase$FlyBase_FBgn_tophit.x,collapse="|"), Table.from2M.greg.apple$FlyBase_FBgn_tophit.x),]
Table.from2Mapple.greg.haw.tor <- Table.from2Mapple.greg.haw[grep(paste(tor.flybase$FlyBase_FBgn_tophit.x,collapse="|"), Table.from2Mapple.greg.haw$FlyBase_FBgn_tophit.x),]

Table.from2M.greg.apple.wnt <- Table.from2M.greg.apple[grep(paste(wnt.flybase$FlyBase_FBgn_tophit.x,collapse="|"), Table.from2M.greg.apple$FlyBase_FBgn_tophit.x),]
Table.from2Mapple.greg.haw.wnt <- Table.from2Mapple.greg.haw[grep(paste(wnt.flybase$FlyBase_FBgn_tophit.x,collapse="|"), Table.from2Mapple.greg.haw$FlyBase_FBgn_tophit.x),]

Table.from2M.greg.apple.insulin<-merge(Table.from2M.greg.apple.insulin,insulin.flybase,by="FlyBase_FBgn_tophit.x")
Table.from2Mapple.greg.haw.insulin<-merge(Table.from2Mapple.greg.haw.insulin,insulin.flybase,by="FlyBase_FBgn_tophit.x")

#will need to be fixed up post R to deal with doubles in flybse gg column
Table.from2M.greg.apple.tor<-merge(Table.from2M.greg.apple.tor,tor.flybase,by="FlyBase_FBgn_tophit.x",all=TRUE)
Table.from2Mapple.greg.haw.tor<-merge(Table.from2Mapple.greg.haw.tor,tor.flybase,by="FlyBase_FBgn_tophit.x",all=TRUE)

Table.from2M.greg.apple.wnt<-merge(Table.from2M.greg.apple.wnt,wnt.flybase,by="FlyBase_FBgn_tophit.x",all=TRUE)
Table.from2Mapple.greg.haw.wnt<-merge(Table.from2Mapple.greg.haw.wnt,wnt.flybase,by="FlyBase_FBgn_tophit.x",all=TRUE)

write.table(Table.from2M.greg.apple.insulin,"AppleBackToHaw2MremoveBadHaw4M/Table.from2M.greg.apple.insulin",quote=F,sep="\t",row.names=FALSE)
write.table(Table.from2M.greg.apple.tor,"AppleBackToHaw2MremoveBadHaw4M/Table.from2M.greg.apple.tor",quote=F,sep="\t",row.names=FALSE)
write.table(Table.from2M.greg.apple.wnt,"AppleBackToHaw2MremoveBadHaw4M/Table.from2M.greg.apple.wnt",quote=F,sep="\t",row.names=FALSE)
write.table(Table.from2Mapple.greg.haw.insulin,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mapple.greg.haw.insulin",quote=F,sep="\t",row.names=FALSE)
write.table(Table.from2Mapple.greg.haw.tor,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mapple.greg.haw.tor",quote=F,sep="\t",row.names=FALSE)
write.table(Table.from2Mapple.greg.haw.wnt,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mapple.greg.haw.wnt",quote=F,sep="\t",row.names=FALSE)

write.table(Table.from2M.greg.apple,"AppleBackToHaw2MremoveBadHaw4M/Table.from2M.greg.apple",quote=F,sep="\t",row.names=FALSE)
write.table(Table.from2Mapple.greg.haw,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mapple.greg.haw",quote=F,sep="\t",row.names=FALSE)

#do haw from 2M, insulin/tor/wnt
haw.fromHAW2M.flybase<-merge(haw.fromHAW2M,flybase.annotations,by="gene_id")

haw.fromHAW2M.flybase.insulin <- haw.fromHAW2M.flybase[grep(paste(insulin.flybase$FlyBase_FBgn_tophit.x,collapse="|"), haw.fromHAW2M.flybase$FlyBase_FBgn_tophit),]
haw.fromHAW2M.flybase.tor <- haw.fromHAW2M.flybase[grep(paste(tor.flybase$FlyBase_FBgn_tophit.x,collapse="|"), haw.fromHAW2M.flybase$FlyBase_FBgn_tophit),]
haw.fromHAW2M.flybase.wnt <- haw.fromHAW2M.flybase[grep(paste(wnt.flybase$FlyBase_FBgn_tophit.x,collapse="|"), haw.fromHAW2M.flybase$FlyBase_FBgn_tophit),]

colnames(insulin.flybase)[colnames(insulin.flybase)=="FlyBase_FBgn_tophit.x" ]<- "FlyBase_FBgn_tophit"
colnames(wnt.flybase)[colnames(wnt.flybase)=="FlyBase_FBgn_tophit.x" ]<- "FlyBase_FBgn_tophit"
colnames(tor.flybase)[colnames(tor.flybase)=="FlyBase_FBgn_tophit.x" ]<- "FlyBase_FBgn_tophit"

Table.haw.fromHAW2M.flybase.insulin<-merge(haw.fromHAW2M.flybase.insulin,insulin.flybase,by="FlyBase_FBgn_tophit")

#will need to be fixed up post R to deal with doubles in flybse gg column
Table.haw.fromHAW2M.flybase.tor<-merge(haw.fromHAW2M.flybase.tor,tor.flybase,by="FlyBase_FBgn_tophit",all=TRUE)
Table.haw.fromHAW2M.flybase.wnt<-merge(haw.fromHAW2M.flybase.wnt,wnt.flybase,by="FlyBase_FBgn_tophit",all=TRUE)

write.table(Table.haw.fromHAW2M.flybase.insulin,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mhaw.greg.haw.insulin",quote=F,sep="\t",row.names=FALSE)
write.table(Table.haw.fromHAW2M.flybase.tor,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mhaw.greg.haw.tor",quote=F,sep="\t",row.names=FALSE)
write.table(Table.haw.fromHAW2M.flybase.wnt,"AppleBackToHaw2MremoveBadHaw4M/Table.from2Mhaw.greg.haw.wnt",quote=F,sep="\t",row.names=FALSE)

#end
