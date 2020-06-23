#R
#GJR 6/14/2019
#Compare urbana poolseq results to drosophila development time experiment from:
# Burke, M.K., Dunham, J.P., Shahrestani, P., Thornton, K.R., Rose, M.R. and Long, A.D., 2010. Genome-wide analysis of a long-term evolution experiment with Drosophila. Nature, 467(7315), p.587.

options(stringsAsFactors=F)
droSnp<-read.table('/media/raglandlab/ExtraDrive4/dros_stuff/mysqlTables/dros.majmin.counts.fisher.snpId',sep="\t",header=T,row.names=NULL)

ind<-rowSums(droSnp[,5:6]) > 10 & rowSums(droSnp[,7:8]) > 10
droSnp<-droSnp[ind,]
droSnp$Fdr<-p.adjust(droSnp$fisher_pval,method='BH')
pdif<-function(x) {
    ( x[1]/( x[1]+x[2] ) ) - ( x[3]/( x[3]+x[4] ) )
}
droSnp$pdif<-apply(droSnp[,5:8],1,pdif)
#ind<-droSnp$Fdr < 0.05 & abs(droSnp$pdif) > 0.65
ind<-droSnp$Fdr < 0.01 
hist(abs(droSnp$pdif[ind]))
#association strength fairly normally distributed, heavily right-skewed, though (alt allele more frequently selected for, apparently)

source('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq/gitrepo/UrbanaPoolSeq/queryDb.r')

queryString<-"select * from annotation where snpId ='%s'"
snpIds<-droSnp$snpId[ind]
#snpIds<-snpIds[1:1000]
annoList<-queryDb(snpIds,queryString,'test_db')
ind<-annoList$effect != "upstream_gene_variant" & annoList$effect != "downstream_gene_variant" & annoList$effect != "intergenic_region" 
annoList<-annoList[ind,] # genome is too compact, need to get rid up up/downstream variants and intergenic
setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
write.table(annoList,'GenesFromDrosDevoPoolseq.txt',sep="\t",row.names=F,quote=F)
#needed to do some id conversion to get to flybase numbers, stored here:
flybase<-read.table('GenesFromDrosDevoPoolseq.flybase.txt',sep="\t",row.names=NULL,header=T,stringsAsFactors=F)
flybase<-unique(flybase$flybase)

queryString<-"select Flybase_FBgn from feature_alias"
allFlybase<-queryDb('null',queryString)
allFlybase<-unique(allFlybase[,1])
flybase<-flybase[flybase %in% allFlybase]

setwd('/media/raglandlab/ExtraDrive4/Urbana_PoolSeq')
pomSnps<-read.table('SnpsAndIndelsInGenesFdrEclHawLt0.01FdrEclAppleLt0.01.txt',sep="\t",row.names=NULL,header=T,stringsAsFactors=F)
#pomSnps<-read.table('SnpsAndIndelsInGenesFdrEclHawLt0.005FdrEclAppleLt0.005.txt',,sep="\t",row.names=NULL,header=T,stringsAsFactors=F)
pomFlybase<-unique(pomSnps$Flybase_FBgn)
pomFlybase<-pomFlybase[-which(is.na(pomFlybase))]




a<-sum(pomFlybase %in% flybase)
b<-sum(!(pomFlybase %in% flybase))
c<-sum(!(flybase %in% pomFlybase))
d<-length(allFlybase)-a-b-c

#                    Snp   noSnp
# in Dros list        a     c
# not in Dros list    b     d

fisher.test(cbind(c(a,c),c(b,d)))
#	Fisher's Exact Test for Count Data for FDR < 0.005 both host races

#data:  cbind(c(a, c), c(b, d))
#p-value = 0.000226
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 1.230239 1.947560
#sample estimates:
#odds ratio 
#  1.554115 

