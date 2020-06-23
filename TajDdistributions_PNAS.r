#R
#GJR 10/29/2019
# load and analyze results of sliding window (500bp) calculations of Tajima's D, caclulated using popoolation

setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
load('UrbanPoolseq.v2.Rdat')

source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')
ids<-"null"
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select* from urbanapoolTajD"
tajAll<-queryDb(ids,queryString)

tajAll$Min<-sapply(tajAll$window_center,function(x) x - 499)
tajAll$Max<-sapply(tajAll$window_center,function(x) x + 500)

scaf<-snps$scaffold[1]
pos<-snps$position[1]

ind<-tajAll$scaffold==scaf & pos >= tajAll$Min & pos <= tajAll$Max
dat<-tajAll[ind,]
dif<-abs(dat$window_center - pos)
ind<-which(dif==min(dif))
tajda<-dat$Urbana_apple_TajD[ind]
tajdh<-dat$Urbana_haw_TajD[ind]

sigInd<-which(snps$fdrHost < 0.05 & abs(snps$pdifHost) > 0.8)
tajSig<-data.frame(cbind(rep(0,length(sigInd)), rep(0,length(sigInd)) ) )
names(tajSig)<-c('tajdApple','tajdHaw')
j=1
for (i in sigInd) {
    scaf<-snps$scaffold[i]
    pos<-snps$position[i]
    ind<-tajAll$scaffold==scaf & pos >= tajAll$Min & pos <= tajAll$Max
    dat<-tajAll[ind,]
    dif<-abs(dat$window_center - pos)
    ind<-which(dif==min(dif))
    if (length(ind) > 1) {
        ind<-ind[1]
    }
    tajSig[j,1]<-dat$Urbana_apple_TajD[ind]
    tajSig[j,2]<-dat$Urbana_haw_TajD[ind]
    j=j+1
}


pdf('TajDAllvsSigHostFDR0.05pdifGT0.8.pdf')
boxplot(tajAll$Urbana_apple_TajD,tajAll$Urbana_haw_TajD,tajSig$tajdApple,tajSig$tajdHaw,las=1)
dev.off()

plotDat<-data.frame(cat=rep('allApple',length((tajAll$Urbana_apple_TajD))),taj=tajAll$Urbana_apple_TajD)
plotDat<-rbind(plotDat, data.frame(cat=rep('allHaw',length((tajAll$Urbana_haw_TajD))),taj=tajAll$Urbana_haw_TajD) )
plotDat<-rbind(plotDat, data.frame(cat=rep('sigApple',length((tajSig$tajdApple))),taj=tajSig$tajdApple) )
plotDat<-rbind(plotDat, data.frame(cat=rep('sigHaw',length((tajSig$tajdHaw))),taj=tajSig$tajdHaw) )
plotDat<-plotDat[!is.na(plotDat$taj),]

library(ggplot2)

pdf('TajDAllvsSigHostFDR0.05pdifGT0.8.violin.pdf')
p<-ggplot(plotDat, aes(x=cat, y=taj)) + 
  geom_violin(trim=FALSE)
p + geom_boxplot(width=0.1,outlier.shape=NA) + theme_classic()
dev.off()


