#R
#GJR 6/12/2020
#perform sliding window analyses of 500bp, non overlapping windows using absolute allele frequency difference


###### 8/6/2019, modified 6/12/2020 GJR ##########
#### sliding window analyses, no overlap #########
##################################################

#use 500bp windows, run in series for host, eclApple, and eclHaw
#need to re-load data each time to cut down environment and memory for parallel execution

# host differences
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
load('UrbanPoolseq.v2.Rdat')
freqCol<-34 #pdifHost
snps<-snps[,c(2:3,freqCol)]
rm(list=setdiff(ls(), c("snps")))
gc()
window<-500
processors=10
scaffolds<-unique(snps$scaffold)
library(doParallel)
cl <- makeCluster(processors)
registerDoParallel(cl)
out<-foreach(scaf=scaffolds) %dopar% {
    out <- data.frame(scaffold=character(),start=numeric(),absFdif=numeric(),stringsAsFactors=F)
    ind<-snps$scaffold==scaf
    srange<-snps$position[ind]
    fdifs<-snps[ind,3]
    start<-1
    while (start <= max(srange)) {
        ind<-srange >= start & srange < (start+window)
        if (sum(ind) >= 1) {
            out<-rbind(out,data.frame(scaffold=scaf,start=start,absFdif=mean(abs(fdifs[ind])),stringsAsFactors=F))
        }
        start<-start+window
    }
    return(out)
}
stopCluster(cl)
library(dplyr)
windows<-bind_rows(out, .id = "column_label")
rm(list=setdiff(ls(), c("windows")))
gc()
save.image('slideHostDiff.Rdat')

# apple Ecl
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
load('UrbanPoolseq.v2.Rdat')
freqCol<-36 #pdifEclApple
snps<-snps[,c(2:3,freqCol)]
rm(list=setdiff(ls(), c("snps")))
gc()
window<-500
processors=10
scaffolds<-unique(snps$scaffold)
library(doParallel)
cl <- makeCluster(processors)
registerDoParallel(cl)
out<-foreach(scaf=scaffolds) %dopar% {
    out <- data.frame(scaffold=character(),start=numeric(),absFdif=numeric(),stringsAsFactors=F)
    ind<-snps$scaffold==scaf
    srange<-snps$position[ind]
    fdifs<-snps[ind,3]
    start<-1
    while (start <= max(srange)) {
        ind<-srange >= start & srange < (start+window)
        if (sum(ind) >= 1) {
            out<-rbind(out,data.frame(scaffold=scaf,start=start,absFdif=mean(abs(fdifs[ind])),stringsAsFactors=F))
        }
        start<-start+window
    }
    return(out)
}
stopCluster(cl)
library(dplyr)
windows<-bind_rows(out, .id = "column_label")
rm(list=setdiff(ls(), c("windows")))
gc()
save.image('slideAppleEclDiff.Rdat')



# haw Ecl
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
load('UrbanPoolseq.v2.Rdat')
freqCol<-35 #pdifEclHaw
snps<-snps[,c(2:3,freqCol)]
rm(list=setdiff(ls(), c("snps")))
gc()
window<-500
processors=10
scaffolds<-unique(snps$scaffold)
library(doParallel)
cl <- makeCluster(processors)
registerDoParallel(cl)
out<-foreach(scaf=scaffolds) %dopar% {
    out <- data.frame(scaffold=character(),start=numeric(),absFdif=numeric(),stringsAsFactors=F)
    ind<-snps$scaffold==scaf
    srange<-snps$position[ind]
    fdifs<-snps[ind,3]
    start<-1
    while (start <= max(srange)) {
        ind<-srange >= start & srange < (start+window)
        if (sum(ind) >= 1) {
            out<-rbind(out,data.frame(scaffold=scaf,start=start,absFdif=mean(abs(fdifs[ind])),stringsAsFactors=F))
        }
        start<-start+window
    }
    return(out)
}
stopCluster(cl)
library(dplyr)
windows<-bind_rows(out, .id = "column_label")
rm(list=setdiff(ls(), c("windows")))
gc()
save.image('slideHawEclDiff.Rdat')


#load results, create indices for loci that are 'outliers', defined as top 97.5 percentile
load('slideHostDiff.Rdat')
slideHost<-windows
rm(windows)
load('slideAppleEclDiff.Rdat')
slideApple<-windows
rm(windows)
load('slideHawEclDiff.Rdat')
slideHaw<-windows
rm(windows)

#create combined data frame
slide<-slideHost[,-1]
names(slide)[3]<-'absFdifHost'
slide$absFdifApple<-slideApple$absFdif
slide$absFdifHaw<-slideHaw$absFdif
threshQuant=0.975
indHost<-slide$absFdifHost > ( quantile(slide$absFdifHost,threshQuant) )
indApple<-slide$absFdifApple > ( quantile(slide$absFdifApple,threshQuant) )
indHaw<-slide$absFdifHaw > ( quantile(slide$absFdifHaw,threshQuant) )

#plot distributions of abs allele freq difs for single snps and for 500bp windows
pdf('AlleleFreqDistsSlideAndSnps.pdf')
par(mfcol = c(3,3))
hist(abs(snps$pdifEclHaw),las=1,xlab='absolute value(frequency difference)',main=paste('median',median(abs(snps$pdifEclHaw))))
hist(slide$absFdifHaw,las=1,xlab='absolute value(frequency difference)',main=paste('median',median(slide$absFdifHaw)))

hist(abs(snps$pdifEclApple),las=1,xlab='absolute value(frequency difference)',main=paste('median',median(abs(snps$pdifEclApple))))
hist(slide$absFdifApple,las=1,xlab='absolute value(frequency difference)',main=paste('median',median(slide$absFdifApple)))

hist(abs(snps$pdifHost),las=1,xlab='absolute value(frequency difference)',main=paste('median',median(abs(snps$pdifHost))))
hist(slide$absFdifHost,las=1,xlab='absolute value(frequency difference)',main=paste('median',median(slide$absFdifHost)))
dev.off()


#correlation for AFD b/t emergence pools for apple vs. haw flies when including all data even including non-associated loci
cor.test(slide$absFdifApple,slide$absFdifHaw)

#	Pearson's product-moment correlation

#data:  slide$absFdifApple and slide$absFdifHaw
#t = 567.07, df = 1524108, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.4160933 0.4187153
#sample estimates:
#      cor 
#0.4174052 

#correlation for AFD b/t emergence pools for apple vs. haw flies including only outlier loci in both apple and haw flies
indBoth<-indApple & indHaw
cor.test(slide$absFdifApple[indBoth],slide$absFdifHaw[indBoth])
#	Pearson's product-moment correlation

#data:  slide$absFdifApple[indBoth] and slide$absFdifHaw[indBoth]
#t = 20.607, df = 8976, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.1926955 0.2321987
#sample estimates:
#      cor 
#0.2125339 


#test overlap in outliers associated with ecl in haw vs. apple flies

#                   sigEclHaw nSigEclHaw
#       sigEclApple       a       b
#       nSigEclApple      c       d
a<-sum(indHaw & indApple)
b<-sum(indHaw & !indApple)
c<-sum(!indHaw & indApple)
d<-nrow(slide)-a-b-c
fisher.test(rbind(c(a,b),c(c,d)))
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 15.01338 15.83525
#sample estimates:
#odds ratio 
#  15.42111


#test overlap in outliers associated with ecl in haw vs. host differences

#                   sigEclHaw nSigEclHaw
#       sigHost       a       b
#       nSigHost      c       d
a<-sum(indHaw & indHost)
b<-sum(indHaw & !indHost)
c<-sum(!indHaw & indHost)
d<-nrow(slide)-a-b-c
fisher.test(rbind(c(a,b),c(c,d)))
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 6.295544 6.713867
#sample estimates:
#odds ratio 
#  6.501782 


#test overlap in outliers associated with ecl in apple vs. host differences

#                   sigEclApple nSigEclApple
#       sigHost       a       b
#       nSigHost      c       d
a<-sum(indApple & indHost)
b<-sum(indApple & !indHost)
c<-sum(!indApple & indHost)
d<-nrow(slide)-a-b-c
fisher.test(rbind(c(a,b),c(c,d)))
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 6.278446 6.694548
#sample estimates:
#odds ratio 
#  6.482914 
