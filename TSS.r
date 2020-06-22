#R
#GJR 5/23/2019
# test for overrepresentation of SNPs near TSS of differentially expressed transcripts

library(GenomicFeatures)

data<-makeTxDbFromGFF('/media/raglandlab/ExtraDrive1/genomes/Rzeph/ref_Rhagoletis_zephyria_1.0_top_level.gff3')

pr<-promoters(data,upstream=1e6, downstream=1e6,columns=c('GENEID','EXONCHROM'))


chroms<-as.vector(seqnames(pr))
geneids<-unlist(as.vector(pr$GENEID))
pranges<-data.frame(ranges(pr))
pranges<-data.frame(pranges,geneids,stringsAsFactors=F)
locs<-unique(pranges$geneids)

startCoord<-data.frame(start=vector(length=length(locs)),end=vector(length=length(locs)))
j=1;
for (i in locs) {
    subrange<-pranges[pranges$geneids==i,]
   startCoord[j,1]<-min(subrange$start)
    startCoord[j,2]<-max(subrange$end)
    j=j+1
}

startCoord$start[startCoord$start < 0] <- 0
startCoord<-data.frame(startCoord,locid=locs)

locsToChroms<-data.frame(chroms,locid=geneids)
a<-match(startCoord$locid,locsToChroms$locid)
coords<-data.frame(startCoord,chrom=locsToChroms$chroms[a])

setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')

write.table(coords,'CoordinatesForCisReg.txt',sep="\t",row.names=F,quote=F)

setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
options(stringsAsFactors=F)
coords<-read.table('CoordinatesForCisReg.txt',header=T,row.names=NULL)


source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')

queryString<-"select * from feature_alias where loc ='%s'"
locs<-coords$locid[1:5]
annoList<-queryDb(locs,queryString)


######### start here ##############


load('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/UrbanPoolseq.v2.Rdat')

load('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/clus.Rdat')
names(appledat)[11]<-'Lowest_FDR_from2M.apple.TimeSeries'
pomAll<-merge(appledat[,c(1:6,11)],hawdat[,c(1:5,10)])
rm(list=setdiff(ls(), c("snps","pomAll")))
gc()

library(GenomicFeatures)
library(GenomicRanges)
data<-makeTxDbFromGFF('/media/raglandlab/ExtraDrive1/genomes/Rzeph/ref_Rhagoletis_zephyria_1.0_top_level.gff3')
gr=GRanges(seqnames=c("NW_016157085.1"),
           ranges=IRanges(start=c(100),end=c(100))
           )





geneids<-pomAll$gene_id[pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 | pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05]

source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')
queryString<-"select loc from feature_alias where gene_id ='%s'"
annoList<-queryDb(geneids,queryString)


pr<-promoters(data,upstream=5e3, downstream=5e3,columns=c('GENEID','EXONCHROM'))
deGr<-subset(pr,GENEID %in% annoList$loc)

nSigBoth<- sum(snps$fdrEclApple < 0.05 & snps$fdrEclHaw < 0.05) 
ind<-sample(1:nrow(snps), nSigBoth )
gr=GRanges(seqnames=snps$scaffold[ind],
           ranges=IRanges(snps$position[ind],end=snps$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)


# for snps sig assoc with ecl time
#point estimate
ind<-snps$fdrEclApple < 0.05 & snps$fdrEclHaw < 0.05
gr=GRanges(seqnames=snps$scaffold[ind],
           ranges=IRanges(snps$position[ind],end=snps$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)
#[1] 1880 snps overlap within 5kb of TSS of differentially expressed genes (DE in BOTH apple and Haw)
# 4661  snps overlap within 5kb of TSS of differentially expressed genes (DE in apple OR Haw)

#added 12/3/2019 -- get set of SNPs in overlap
a<-intersect(gr,deGr,ignore.strand=T)
snpSet<-as.data.frame(a)
names(snpSet)[1:2]<-c('scaffold','position')

#load previously generated dataframe with chromosome assignments and LD groups (have to be w/in 200bp of mapped RAD snp) for SNPs
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/')
snpsWChrom<-read.table('snpsWChrom.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)

snpSet<-merge(snpsWChrom[,1:4],snpSet)

s<-Sys.time()
nRep=10000
out<-vector(length=nRep)
for (i in 1:nRep) {
    ind<-sample(1:nrow(snps), nSigBoth )
    gr=GRanges(seqnames=snps$scaffold[ind],
               ranges=IRanges(snps$position[ind],end=snps$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    out[i]<-sum(a > 0)
}
Sys.time() - s
#(1-ecdf(out)(1880)) # = 0
#so p < 0.0001, not one in 10,000 random sets, though expectation was
#median(out) # = 1685 # for DE BOTH

(1-ecdf(out)(4661)) # = 0, with expectation
median(out) # = 3976 for DE in apple OR Haw



setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
pdf('SigEclAndDETimeDistAndPointEstimate.pdf')
boxplot(out,ylim=c(1500,1900),las=1)
points(1,1880,pch=20)
dev.off()

ind<-snps$fdrEclApple < 0.05 & snps$fdrEclHaw < 0.05
gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
           ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
           )
deAndSnp<-findOverlaps(gr,deGr)
deInd<-to(deAndSnp)
geneids<-unique(unlist(as.vector(deGr$GENEID[deInd])))

source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(geneids,queryString)
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/Candidates')
write.table(annoList,'SnpEclAppleAndHawFdrLt0.05AndAFDGt0.5.noChromFilter.txt',sep="\t",row.names=F,quote=F)




load('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/clus.Rdat')
geneids<-appledat$gene_id[appledat$Lowest_FDR_FruitBetweenMonths < 0.05]

source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')
queryString<-"select loc from feature_alias where gene_id ='%s'"
annoList<-queryDb(geneids,queryString)
deGr<-subset(pr,GENEID %in% annoList$loc)


#point estimate
ind<-(snps$fdrEclApple < 0.05 | snps$fdrEclHaw < 0.05) & snps$fdrHost < 0.05
gr=GRanges(seqnames=snps$scaffold[ind],
           ranges=IRanges(snps$position[ind],end=snps$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)
# 35 snps are associated with eclosion and host race and DE between host races

nSigBoth<-sum( (snps$fdrEclApple < 0.05 | snps$fdrEclHaw < 0.05) & snps$fdrHost < 0.05) 
s<-Sys.time()
nRep=10000
out<-vector(length=nRep)
for (i in 1:nRep) {
    ind<-sample(1:nrow(snps), nSigBoth )
    gr=GRanges(seqnames=snps$scaffold[ind],
               ranges=IRanges(snps$position[ind],end=snps$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    out[i]<-sum(a > 0)
}
Sys.time() - s
(1-ecdf(out)(35)) # = 0.2324
#so p = 0.2324, no clear signal of snps lining up with host differences in diapause expression
median(out) # = 1685


#point estimate
ind<-snps$fdrHost < 0.05
gr=GRanges(seqnames=snps$scaffold[ind],
           ranges=IRanges(snps$position[ind],end=snps$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)
# 476 snps are associated with eclosion and host race and DE between host races

nSigBoth<-sum( (snps$fdrHost < 0.05) )
s<-Sys.time()
nRep=10000
out<-vector(length=nRep)
for (i in 1:nRep) {
    ind<-sample(1:nrow(snps), nSigBoth )
    gr=GRanges(seqnames=snps$scaffold[ind],
               ranges=IRanges(snps$position[ind],end=snps$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    out[i]<-sum(a > 0)
}
Sys.time() - s
(1-ecdf(out)(476)) # = 0
#so p < 0.0001, here there is acctually a signal of host race snps being near DE genes
median(out) #386



processors=2

nRep=50
library(doParallel)
cl <- makeCluster(processors)
registerDoParallel(cl)
out<-foreach(i=1:nRep,.combine='c') %dopar% {
    ind<-sample(1:nrow(snps), nSigBoth )
    gr=GRanges(seqnames=snps$scaffold[ind],
               ranges=IRanges(snps$position[ind],end=snps$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    #out[i]<-sum(a > 0)
    return( sum(a > 0) )
}
stopCluster(cl)
return(1-ecdf(out)(medVal))

processors=2
nRep=50
library(doParallel)
cl <- makeCluster(processors)
registerDoParallel(cl)
out<-foreach(i=1:nRep,.combine='c') %dopar% {
    b<-snps$position[i]
    return( sum(a > 0) )
}
stopCluster(cl)



#point estimate
ind<-snps$fdrEclApple < 0.05 & snps$fdrEclHaw < 0.05
gr=GRanges(seqnames=snps$scaffold[ind],
           ranges=IRanges(snps$position[ind],end=snps$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)


##### 6/13/2019 -- Thresholding to 0.4 freq diff to account for sex ratio diffs among bulks -- note, need to load 'snpsWChrom'

#load previously generated dataframe with chromosome assignments and LD groups (have to be w/in 200bp of mapped RAD snp) for SNPs
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/')
snpsWChrom<-read.table('snpsWChrom.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
snpsWld<-read.table('snpsWld.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
ind<-snpsWChrom$snpId %in% snps$snpId
snpsWChrom<-snpsWChrom[ind,]


#point estimate
ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrEclHaw < 0.05 & abs(snpsWChrom$pdifEclApple) > 0.4 &  abs(snpsWChrom$pdifEclHaw) > 0.4
gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
           ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)
#[1] 740 snps overlap within 5kb of TSS of differentially expressed genes

nSigBoth<- sum(snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrEclHaw < 0.05 & abs(snpsWChrom$pdifEclApple) > 0.4 &  abs(snpsWChrom$pdifEclHaw) > 0.4) 

s<-Sys.time()
nRep=10000
out<-vector(length=nRep)
for (i in 1:nRep) {
    ind<-sample(1:nrow(snpsWChrom), nSigBoth )
    gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
               ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    out[i]<-sum(a > 0)
}
Sys.time() - s
(1-ecdf(out)(740)) # = 0.0012
#so p = 0.0012
median(out) # = 665






#point estimate
ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrEclHaw < 0.05 & abs(snpsWChrom$pdifEclApple) > 0.5 &  abs(snpsWChrom$pdifEclHaw) > 0.5
gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
           ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)
#[1] 157 snps overlap within 5kb of TSS of differentially expressed genes

nSigBoth<- sum(snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrEclHaw < 0.05 & abs(snpsWChrom$pdifEclApple) > 0.5 &  abs(snpsWChrom$pdifEclHaw) > 0.5) 

s<-Sys.time()
nRep=10000
out<-vector(length=nRep)
for (i in 1:nRep) {
    ind<-sample(1:nrow(snpsWChrom), nSigBoth )
    gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
               ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    out[i]<-sum(a > 0)
}
Sys.time() - s
(1-ecdf(out)(157)) # = 0
#so p < 0.0001
median(out) # = 24
#So, at a threshold of 0.5 allele freq dif, there's actually a reasonably strong signal of cis regulation, so those 157 could be good candidates


ind<-snpsWChrom$fdrEclApple < 0.05 & snpsWChrom$fdrEclHaw < 0.05 & abs(snpsWChrom$pdifEclApple) > 0.5 &  abs(snpsWChrom$pdifEclHaw) > 0.5
gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
           ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
           )
deAndSnp<-findOverlaps(gr,deGr)
deInd<-to(deAndSnp)
geneids<-unique(unlist(as.vector(deGr$GENEID[deInd])))

source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')
queryString<-"select * from feature_alias where loc ='%s'"
annoList<-queryDb(geneids,queryString)
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/Candidates')
write.table(annoList,'SnpEclAppleAndHawFdrLt0.05AndAFDGt0.5.txt',sep="\t",row.names=F,quote=F)



########## check genes DE b/t hosts with genes with AFD for host and ecl ##########3

geneids<-appledat$gene_id[appledat$Lowest_FDR_FruitBetweenMonths < 0.05]

load('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved/clus.Rdat')
source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')
queryString<-"select loc from feature_alias where gene_id ='%s'"
annoList<-queryDb(geneids,queryString)
# need to remake 'data' object for granges, gets overwritten in the load statement

pr<-promoters(data,upstream=5e3, downstream=5e3,columns=c('GENEID','EXONCHROM'))
deGr<-subset(pr,GENEID %in% annoList$loc)





#point estimate
ind<-snpsWChrom$fdrEclApple < 0.05  & snpsWChrom$fdrHost < 0.05 & abs(snpsWChrom$pdifEclApple) > 0.4
ind<-ind | ( snpsWChrom$fdrEclHaw < 0.05  & snpsWChrom$fdrHost < 0.05 & abs(snpsWChrom$pdifEclHaw) > 0.4 )
gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
           ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)
#[1] 10 snps overlap within 5kb of TSS of differentially expressed genes

nSigBoth<-sum(ind)

s<-Sys.time()
nRep=10000
out<-vector(length=nRep)
for (i in 1:nRep) {
    ind<-sample(1:nrow(snpsWChrom), nSigBoth )
    gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
               ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    out[i]<-sum(a > 0)
}
Sys.time() - s
(1-ecdf(out)(10)) # = 0.7342
#so p = 0.7342
median(out) # = 13

#point estimate
ind<-snpsWChrom$fdrEclApple < 0.05  & snpsWChrom$fdrHost < 0.05 & abs(snpsWChrom$pdifHost) > 0.4
ind<-ind | ( snpsWChrom$fdrEclHaw < 0.05  & snpsWChrom$fdrHost < 0.05 & abs(snpsWChrom$pdifEclHaw) > 0.4 )
gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
           ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
           )
suppressWarnings(a<-countOverlaps(gr,deGr))
sum(a > 0)
#[1] 10 snps overlap within 5kb of TSS of differentially expressed genes

nSigBoth<-sum(ind)

s<-Sys.time()
nRep=10000
out<-vector(length=nRep)
for (i in 1:nRep) {
    ind<-sample(1:nrow(snpsWChrom), nSigBoth )
    gr=GRanges(seqnames=snpsWChrom$scaffold[ind],
               ranges=IRanges(snpsWChrom$position[ind],end=snpsWChrom$position[ind])
               )
    suppressWarnings(a<-countOverlaps(gr,deGr))
    out[i]<-sum(a > 0)
}
Sys.time() - s
(1-ecdf(out)(10)) # = 0.7342
#so p = 0.7342
median(out) # = 13
