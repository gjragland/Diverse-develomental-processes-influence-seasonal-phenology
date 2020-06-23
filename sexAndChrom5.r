#R
#GJR 6/13/2020

##### GJR 6/13/2019 -- check correlation for allele freq diffs excluding chrom 5
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
load('UrbanPoolseq.v2.Rdat')


#load previously generated dataframe with chromosome assignments and LD groups (have to be w/in 200bp of mapped RAD snp) for SNPs
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/')
snpsWChrom<-read.table('snpsWChrom.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
snpsWld<-read.table('snpsWld.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
ind<-snpsWChrom$snpId %in% snps$snpId
snpsWChrom<-snpsWChrom[ind,]
ind<-snpsWld$snpId %in% snps$snpId
snpsWld<-snpsWld[ind,]



ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrEclHaw < 0.05
cor.test(snpsWChrom$pdifEclHaw[ind],snpsWChrom$pdifEclApple[ind])
#      cor 
#0.9570023 -- 21,105 SNPs
ind<-ind & snpsWChrom$chromosome != 5
cor.test(snpsWChrom$pdifEclHaw[ind],snpsWChrom$pdifEclApple[ind])
#     cor 
#0.961435 -- 20,072 SNPs
#correlation for ecl-association persists when chromosme 5 is removed

ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrHost < 0.05
cor.test(snpsWChrom$pdifEclApple[ind],snpsWChrom$pdifHost[ind])
#       cor 
#0.04953328 --242 SNPs
ind<-ind & snpsWChrom$chromosome != 5
cor.test(snpsWChrom$pdifEclApple[ind],snpsWChrom$pdifHost[ind])
#       cor 
#0.04707222 --211 SNPs

ind<-snpsWChrom$fdrEclHaw < 0.05 & snpsWChrom$fdrHost < 0.05
cor.test(snpsWChrom$pdifEclHaw[ind],snpsWChrom$pdifHost[ind])
#      cor 
#-0.169298 --203 SNPs
ind<-ind & snpsWChrom$chromosome != 5
cor.test(snpsWChrom$pdifEclHaw[ind],snpsWChrom$pdifHost[ind])
#       cor 
#-0.1867622 --195 SNPs


#also check to see whether sex-determining factors are highly differentiated b/t the bulks
ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$chromosome == 5
a<-snpsWChrom$pdifEclApple[ind]
ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$chromosome != 5
b<-snpsWChrom$pdifEclApple[ind]
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/SexRatioAnalyses')
pdf('AlleleFreqDifsEclosionBulksChrom5VsAll.apple.pdf')
boxplot(abs(a),abs(b),las=1,xaxt='n',ylab='Allele Frequency Difference (absolute value)', main='Differences b/t eclosion bulks, Apple Flies')
axis(1, at=1:2, labels=c('Chromosome V (38,804 loci)','All Other Chroms. (177,190 loci)'))
dev.off()
ind<-snpsWChrom$fdrEclHaw < 0.05 & snpsWChrom$chromosome == 5
a<-snpsWChrom$pdifEclHaw[ind]
ind<-snpsWChrom$fdrEclHaw < 0.05 & snpsWChrom$chromosome != 5
b<-snpsWChrom$pdifEclHaw[ind]
pdf('AlleleFreqDifsEclosionBulksChrom5VsAll.haw.pdf')
boxplot(abs(a),abs(b),las=1,xaxt='n',ylab='Allele Frequency Difference (absolute value)', main='Differences b/t eclosion bulks, Haw Flies')
axis(1, at=1:2, labels=c('Chromosome V (19,443 loci)','All Other Chroms. (152,222)'))
dev.off()

ind<-snpsWChrom$fdrHost < 0.05 & snpsWChrom$chromosome == 5
a<-snpsWChrom$pdifHost[ind]
ind<-snpsWChrom$fdrHost < 0.05 & snpsWChrom$chromosome != 5
b<-snpsWChrom$pdifHost[ind]
pdf('AlleleFreqDifsHostChrom5VsAll.pdf')
boxplot(abs(a),abs(b),las=1,xaxt='n',ylab='Allele Frequency Difference (absolute value)', main='Differences b/t Apple and Haw Flies')
axis(1, at=1:2, labels=c('Chromosome V (893 loci)','All Other Chroms. (3,437 loci)'))
dev.off()


ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrEclHaw < 0.05 & abs(snpsWChrom$pdifEclApple) > 0.4 &  abs(snpsWChrom$pdifEclHaw) > 0.4
cor.test(snpsWChrom$pdifEclHaw[ind],snpsWChrom$pdifEclApple[ind])
#      cor 
#0.9662702 -- 15,013 SNPs 
# So, if alternatively we focus on SNPs differentiated well above the expectation for sex ratio differences, we still get a very strong, positive correlation
