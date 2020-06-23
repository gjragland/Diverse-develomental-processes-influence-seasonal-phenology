#R
#GJR 5/23/2019
# Test for overrepresentation of SNPs near TSS of DE transcripts

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


#permute 10k times to generate null distribution
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





