#R
#GJR 6/22/2020

#annotated analyses of poolseq data to accompany PNAS ms.
#organized by order of results in the manuscript

################################################################################### 
############################ INITIAL DATA MANIPULATION ############################ 
################################################################################### 

#import, deterimine which rows are indels, and creates 2 new data sets, snps and indels
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
#flat file of SQL database with SNP and INDEL allele frequencies and results of Fisher's Exact tests
data<-read.table('poolseqSnpsUrbana.txt',stringsAsFactors=F,row.names=NULL,header=T)
gt1<-function(x) {
    out<-F
    if (nchar(x) > 1) {out<-T}
    return(out)
}
isIndel<-sapply(data$ref,gt1) | sapply(data$alt,gt1)
snps<-data[!isIndel,]
indels<-data[isIndel,]
rm(data)
gc()

#fix problem with character columns
a<-apply(snps[,6:23],2,as.numeric)
snps[,6:23]<-a
rm(a)
gc()

#filter to min 10x per pop coverage
ind<-rowSums(snps[,6:7]) >= 10 & rowSums(snps[,8:9]) >= 10 
snps<-snps[ind,]



#new, 6/26/2018
#import pvals from LR tests implemented on summit using LRfreqTestl.r; see github for details
hrpval<-read.table('LRpvalsHostRace.txt',header=T,row.names=NULL)
aepval<-read.table('LRpvalsAppleEcl.txt',header=T,row.names=NULL)
hepval<-read.table('LRpvalsHawEcl.txt',header=T,row.names=NULL)

snps<-data.frame(snps,LRhostPval=hrpval,LReclApple=aepval,LReclHaw=hepval)
names(snps)[25:27]<-c('LRhostPval','LReclApplePval','LReclHawPval')
rm(hrpval,aepval,hepval)
gc()

#filter very high coverage loci (probably repetitive regions)
coverage<-rowSums(snps[,6:7])
highVal<-quantile(coverage,0.99)
snps<-snps[coverage < highVal,]

#filter two other screwy loci
ind<-is.na(snps$LReclHawPval)
snps<-snps[!ind,]

#estimate fdr for snps based on fisher test
snps$fdrHost<-p.adjust(snps$urbana_appleave_hawave_fisher_pvalue,method='BH')
snps$fdrEclApple<-p.adjust(snps$urbana_appleearly_applelate_fisher_pvalue,method='BH')
snps$fdrEclHaw<-p.adjust(snps$urbana_hawearly_hawlate_fisher_pvalue,method='BH')


##estimate fdr based on LR test p-values;
# note that the LR calculation above suffers from underflow, yielding NAs at high coverage. Used log-sum-exp trick, but it's an approximation
snps$fdrHostLR<-p.adjust(snps$LRhostPval,method='BH')
snps$fdrEclAppleLR<-p.adjust(snps$LReclApplePval,method='BH')
snps$fdrEclHawLR<-p.adjust(snps$LReclHawPval,method='BH')

#fisher and lr test qualitatively similar (results correlated at r = 0.9)
#cor.test(snps$urbana_appleave_hawave_fisher_pvalue,snps$LRhostPval)
##note, fisher test actually seems to be more conservative, yields a larger p-value in most cases
#e.g., ...
a<-snps$urbana_appleave_hawave_fisher_pvalue - snps$LRhostPval
hist(a)

#thus, we move forward using Fisher's exact test results, a la 'popoolation'


################################################################################### 
################################################################################### 
################################################################################### 


############################################################################################################## 
######### Snp overrepresentation stats and allele frequencies, focusing on eclosion associations ############## 
############################################################################################################## 

#nSnps
nrow(snps) # 28133855
#nIndels
nrow(indels)  # 1667189

#significant association, apple
sum(snps$fdrEclApple < 0.05) #586460

#significant association, haw
sum(snps$fdrEclHaw < 0.05) #459563

#significant association, both
sum(snps$fdrEclHaw < 0.05 & snps$fdrEclApple < 0.05) #53335


#Test of overrepresentation of overlapping SNPSs
a<-sum( snps$fdrEclHaw < 0.05 & snps$fdrEclApple < 0.05 )
b<-sum( !(snps$fdrEclHaw < 0.05) & snps$fdrEclApple < 0.05 )
c<-sum( snps$fdrEclHaw < 0.05 & !(snps$fdrEclApple < 0.05) )
d<-sum( !(snps$fdrEclHaw < 0.05) & !(snps$fdrEclApple < 0.05) )

#                   sigEclHaw nSigEclHaw
#       sigEclApple       a       b
#       nSigEclApple      c       d

fisher.test(rbind(c(a,b),c(c,d)))
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 6.62091 6.74868
#sample estimates:
#odds ratio 
#  6.683671 


#Calculate allele frequency diffs for snps, apple - haw and early - late
snps$pdifHost<-(snps$urbana_appleave_maj/(snps$urbana_appleave_maj+snps$urbana_appleave_min)) -
    (snps$urbana_hawave_maj/(snps$urbana_hawave_maj+snps$urbana_hawave_min))
snps$pdifEclHaw<-(snps$urbana_hawearly_maj/(snps$urbana_hawearly_maj+snps$urbana_hawearly_min)) -
    (snps$urbana_hawlate_maj/(snps$urbana_hawlate_maj+snps$urbana_hawlate_min))
snps$pdifEclApple<-(snps$urbana_appleearly_maj/(snps$urbana_appleearly_maj+snps$urbana_appleearly_min)) -
    (snps$urbana_applelate_maj/(snps$urbana_applelate_maj+snps$urbana_applelate_min))



# r = 0.96 correlation b/t allele frequency difs for haw ecl vs apple ecl
#tested using permutation test, see github 'permCor.r'
#result: p < 0.0001 (10k permutations, no values as extreme as observed)
#non-permuted approximation...
ind<-snps$fdrEclApple < 0.05 & snps$fdrEclHaw < 0.05
cor.test(snps$pdifEclApple[ind],snps$pdifEclHaw[ind])

#	Pearson's product-moment correlation

#data:  snps$pdifEclApple[ind] and snps$pdifEclHaw[ind]
#t = 747.24, df = 53333, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.9546656 0.9561456
#sample estimates:
#      cor 
#0.9554116 


#For tests of robustness to diffrent read depths, see 'readDepth.r'
#For sliding window analyses see 'slide500.r'
#For analysis of the putative effects of sex on pool differences, see 'sexAndChrom5.r' 


######################################################################################################## 
######### various analyses including indels, still focusing on ecl associations  ####################### 
######################################################################################################## 


#fix problem with character columns
a<-apply(indels[,18:23],2,as.numeric)
indels[,18:23]<-a
rm(a)
gc()

#filter to min 10x per pop coverage
getCounts<-function(x) {sum(as.numeric(unlist(strsplit(x,','))))}
ind<-( sapply(indels$urbana_appleave_maj,getCounts) + sapply(indels$urbana_appleave_min,getCounts) ) >= 10 &
    ( sapply(indels$urbana_hawave_maj,getCounts) + sapply(indels$urbana_hawave_min,getCounts) ) >= 10
a<-ind
a[is.na(a)]<-F
indels<-indels[a,]

coverage<-sapply(indels$urbana_appleave_maj,getCounts) + sapply(indels$urbana_appleave_min,getCounts)
highVal<-quantile(coverage,0.99)
indels<-indels[coverage < highVal,]

#fix problem with character columns
a<-apply(indels[,18:23],2,as.numeric)
indels[,18:23]<-a
rm(a)
gc()

#filter to min 10x per pop coverage
getCounts<-function(x) {sum(as.numeric(unlist(strsplit(x,','))))}
ind<-( sapply(indels$urbana_appleave_maj,getCounts) + sapply(indels$urbana_appleave_min,getCounts) ) >= 10 &
    ( sapply(indels$urbana_hawave_maj,getCounts) + sapply(indels$urbana_hawave_min,getCounts) ) >= 10
a<-ind
a[is.na(a)]<-F
indels<-indels[a,]

coverage<-sapply(indels$urbana_appleave_maj,getCounts) + sapply(indels$urbana_appleave_min,getCounts)
highVal<-quantile(coverage,0.99)
indels<-indels[coverage < highVal,]


#perform fisher exact tests (replaces existing columns that were previously miscalculated)
fishTest<-function(x) {
    x<-as.character(x)
    if (sum(grepl("NULL",x)) > 0) {pval<-NA} else {
    maj1<-as.numeric(unlist(strsplit(x[1],',')))
    min1<-as.numeric(unlist(strsplit(x[2],',')))
    maj2<-as.numeric(unlist(strsplit(x[3],',')))
    min2<-as.numeric(unlist(strsplit(x[4],',')))
    mat<-cbind(c(maj1,min1),c(maj2,min2))
    stat<-fisher.test(mat)
    pval<-stat$p }
    return(pval)
}

library(parallel)
nCores=8
cl <- makeCluster(nCores)
dat<-indels[,6:9]
clusterExport(cl=cl, varlist=c("dat","fishTest"),envir=environment())
out<-unlist(parApply(cl, dat, 1, fishTest ))
stopCluster(cl)
indels$urbana_appleave_hawave_fisher_pvalue<-out

cl <- makeCluster(nCores)
dat<-indels[,10:13]
clusterExport(cl=cl, varlist=c("dat","fishTest"),envir=environment())
out<-unlist(parApply(cl, dat, 1, fishTest ))
stopCluster(cl)
indels$urbana_appleearly_applelate_fisher_pvalue<-out

cl <- makeCluster(nCores)
dat<-indels[,14:17]
clusterExport(cl=cl, varlist=c("dat","fishTest"),envir=environment())
out<-unlist(parApply(cl, dat, 1, fishTest ))
stopCluster(cl)
indels$urbana_hawearly_hawlate_fisher_pvalue<-out


noNA<-function(x) {
    x[!is.na(x)]
}

#estimate fdr
indels$fdrHost<-p.adjust(indels$urbana_appleave_hawave_fisher_pvalue,method='BH')
indels$fdrEclApple<-p.adjust(indels$urbana_appleearly_applelate_fisher_pvalue,method='BH')
indels$fdrEclHaw<-p.adjust(indels$urbana_hawearly_hawlate_fisher_pvalue,method='BH')

indels$fdrHost[is.na(indels$fdrHost)]<-1
indels$fdrEclApple[is.na(indels$fdrEclApple)]<-1
indels$fdrEclHaw[is.na(indels$fdrEclHaw)]<-1


#number sig apple ecl
sum(indels$fdrEclApple < 0.05) #85523

#number sig haw ecl
sum(indels$fdrEclHaw < 0.05) #69720

#number sig both
sum(indels$fdrEclHaw < 0.05 & indels$fdrEclApple < 0.05) #9382


######################################################################################################## 
######### genomic distribution of variants, still focusing on ecl associations  ####################### 
########################################################################################################


#load previously generated dataframe with chromosome assignments and LD groups (have to be w/in 200bp of mapped RAD snp) for SNPs
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/')
snpsWChrom<-read.table('snpsWChrom.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
snpsWld<-read.table('snpsWld.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
ind<-snpsWChrom$snpId %in% snps$snpId
snpsWChrom<-snpsWChrom[ind,]
ind<-snpsWld$snpId %in% snps$snpId
snpsWld<-snpsWld[ind,]


#snps associated with eclosion highly overrepresented on chroms 1 and 3, p < 2.2e-16, df = 4
ind<-snpsWChrom$fdrEclHaw < 0.01 & snpsWChrom$fdrEclApple < 0.01

a<-table(snpsWChrom$chromosome[ind])
b<-aggregate(position ~ chromosome,data=snpsWChrom,length)[,2]
probs<-b/sum(b)
expected<-probs*sum(a)
#chi squared goodness of fit test, predicted probs are even across chroms (null hypothesis)
chisq.test(a,p=probs)
#snps associated with eclosion highly overrepresented on chroms 1 and 3, p < 2.2e-16, df = 4
# note, expected values are the proportion of total snps on each chrom X total number of snps

#	Chi-squared test for given probabilities

#data:  a
#X-squared = 4587.8, df = 4, p-value < 2.2e-16


######################################################################################################## 
######### overlap with Loci from Ragland et al. 2017 ################################################### 
########################################################################################################


#see 'CompareToFenvilleRAD_PNAS.r'


######################################################################################################## 
######### Lists of flybase annotations for gene models with significant SNPS and indels ################ 
########################################################################################################

#indels AND snps
fdrHost<-p.adjust(c(indels$urbana_appleave_hawave_fisher_pvalue,snps$urbana_appleave_hawave_fisher_pvalue),method='BH')
fdrEclApple<-p.adjust(c(indels$urbana_appleearly_applelate_fisher_pvalue,snps$urbana_appleearly_applelate_fisher_pvalue),method='BH')
fdrEclHaw<-p.adjust(c(indels$urbana_hawearly_hawlate_fisher_pvalue,snps$urbana_hawearly_hawlate_fisher_pvalue),method='BH')

fdrHost[is.na(fdrHost)]<-1
fdrEclApple[is.na(fdrEclApple)]<-1
fdrEclHaw[is.na(fdrEclHaw)]<-1

uniScafs<-unique(snps$scaffold)

a<-sum( fdrEclHaw < 0.01 & fdrEclApple < 0.01 ) #9822 snps and indels sig assoc. with ecl at fdr < 0.01
b<-sum( !(fdrEclHaw) < 0.01 & fdrEclApple < 0.01 )
c<-sum( fdrEclHaw < 0.01 & !(fdrEclApple < 0.01) )
d<-sum( !(fdrEclHaw) < 0.01 & !(fdrEclApple < 0.01) )

fisher.test(rbind(c(a,b),c(c,d)))
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
#  9.787423 10.207326
#sample estimates:
#odds ratio 
#  9.994476 

#pull annotations and write to file

#requires function to query local SQL database
source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')

ind<-fdrEclHaw < 0.01 & fdrEclApple < 0.01
snpIds<-c(indels$snpId,snps$snpId)
ids<-snpIds[ind]
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)
#candidates affecting eclosion timing
#write.table(annoList,'SnpsAndIndelsInGenesFdrEclHawLt0.1FdrEclAppleLt0.01.txt',sep="\t",quote=F,row.names=F)

#Generate lists to check for enrichment of similar funcitonal categories on high LD regions of chroms 1 and 3
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/')
#snpsWChrom<-read.table('snpsWChrom.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
snpsWld<-read.table('snpsWld.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)

ind<-( snpsWld$LDgr == 'H' | snpsWld$LDgr == 'M' ) & ( snpsWld$chrom == 1 | snpsWld$chrom == 3 )
ids<-snpsWld$snpId[ind]
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
locs<-locs[!duplicated(locs)]
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)

setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
write.table(annoList,'snpsAnnoHighAndMedLDChrom1And3.txt',row.names=F,quote=F,sep='\t')


######################################################################################################## 
######### Check whether sig SNPs occur in/near (proximity to TSS) DE transcripts ####################### 
########################################################################################################

#See 'TSS_PNAS.r'


######################################################################################################## 
######### Test for overlap of ecl-associated SNPs with Burke et al. #################################### 
########################################################################################################

#See 'CompareToDrosophila_PNAS.r'


######################################################################################################## 
######### Analyses of host race differences ############################################################ 
########################################################################################################

#number SNPs and Indes associated with host race divergence
sum(fdrHost < 0.05) #14650
#number SNPs and Indes associated with host race divergence and apple ecl
sum(fdrHost < 0.05 & fdrEclApple < 0.05) #698
#number SNPs and Indes associated with host race divergence and haw ecl
sum(fdrHost < 0.05 & fdrEclHaw < 0.05) #596

a<-sum( fdrEclHaw < 0.05 & fdrHost < 0.05 ) 
b<-sum( !(fdrEclHaw) < 0.05 & fdrHost < 0.05 )
c<-sum( fdrEclHaw < 0.05 & !(fdrHost < 0.05) )
d<-sum( !(fdrEclHaw) < 0.05 & !(fdrHost < 0.05) )

#                   sigEclHaw nSigEclHaw
#       sigHost       a       b
#       nSigHost      c       d

#significant overlap b/t variants associated with ecl in haw and those associated with host divergence
fisher.test(rbind(c(a,b),c(c,d)))
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 2.196912 2.593120
#sample estimates:
#odds ratio 
#  2.388621 

a<-sum( fdrEclApple < 0.05 & fdrHost < 0.05 ) 
b<-sum( !(fdrEclApple) < 0.05 & fdrHost < 0.05 )
c<-sum( fdrEclApple < 0.05 & !(fdrHost < 0.05) )
d<-sum( !(fdrEclApple) < 0.05 & !(fdrHost < 0.05) )

#                   sigEclApple nSigEclApple
#       sigHost       a       b
#       nSigHost      c       d

#significant overlap b/t variants associated with ecl in haw and those associated with host divergence
fisher.test(rbind(c(a,b),c(c,d)))
#	Fisher's Exact Test for Count Data

#data:  rbind(c(a, b), c(c, d))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 2.026043 2.362565
#sample estimates:
#odds ratio 
#  2.189408 


#write lists of annotations for gene models containing sig SNPs and indels
ind<-(fdrEclHaw < 0.05 | fdrEclApple < 0.05) & fdrHost < 0.05)
snpIds<-c(indels$snpId,snps$snpId)
ids<-snpIds[ind]
#query<-sprintf("select * from annotation where snpId = %s", id)
queryString<-"select * from annotation where snpId = %s"
geneList<-queryDb(ids,queryString)
locs<-geneList$loc[geneList$effect!='intergenic_region']
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(locs,queryString)
#candidates affecting eclosion timing and host-differentiated
write.table(annoList,'SnpsAndIndelsInGenesFdrEclHawLt0.05OrFdrEclAppleLt0.05AndFdrHostLt0.05.txt',sep="\t",quote=F,row.names=F)



#Analysis of Tajima's D: see 'TajDdistributions_PNAS.r'


