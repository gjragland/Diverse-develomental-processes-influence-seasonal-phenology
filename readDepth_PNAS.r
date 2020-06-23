#R
#GJR 5/25/2020
#check robustness of correlation results to coverage depth

##################################################################
###### 5/25/2020 Re-analyses for PNAS                   ##########
##################################################################

############## Correlation and read depth ########################
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
load('UrbanPoolseq.v2.Rdat')


cor.test(snps$pdifEclHaw,snps$pdifEclApple) #28,133,855 loci
##overall correlation is actually quite high for including a bunch of neutral/low coverage loci
#	Pearson's product-moment correlation

#data:  snps$pdifEclHaw and snps$pdifEclApple
#t = 1504.7, df = 28133853, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2725796 0.2732636
#sample estimates:
#      cor 
#0.2729216 

depth<-apply(snps[10:17],1,sum)
mean(depth) #119

depthHaw<-apply(snps[10:13],1,sum)
depthApple<-apply(snps[14:17],1,sum)


#there is a modest negative correlation between read depth and allele freq diff between bulks
cor.test(abs(snps$pdifEclHaw),depthHaw) #r = -0.1313835
cor.test(abs(snps$pdifEclApple),depthApple) #r = -0.1233404


#however, filtering to 60x and 100x min coverage per pop and no filtering for AFD -> higher correlation than unfiltered data

ind<-depthHaw >= 60 & depthApple >=60
cor.test(snps$pdifEclHaw[ind],snps$pdifEclApple[ind]) # r = 0.3412434, 7673103 loci

ind<-depthHaw >= 100 & depthApple >=100
cor.test(snps$pdifEclHaw[ind],snps$pdifEclApple[ind]) # r = 0.323521, 1153338 loci


#similarly, filtering for coverage and for significant loci -> very high correlations

ind<-depthHaw >= 60 & depthApple >= 60 & snps$fdrEclHaw < 0.05 & snps$fdrEclApple < 0.05
cor.test(snps$pdifEclHaw[ind],snps$pdifEclApple[ind]) # r = 0.9607694, 26829 loci


ind<-depthHaw >= 100 & depthApple >= 100 & snps$fdrEclHaw < 0.05 & snps$fdrEclApple < 0.05
cor.test(snps$pdifEclHaw[ind],snps$pdifEclApple[ind]) # r = 0.9588505, 4048 loci

ind<-snps$fdrEclHaw < 0.05 & snps$fdrEclApple < 0.05
