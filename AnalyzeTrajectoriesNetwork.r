
#load up differential expression data
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved')
load('clus.Rdat')
#remove a custom function, that interferes with a WGCNA function used below
rm(collapseRows)

#set up the matrix of gene DE over time in both host races
names(appledat)[11]<-'Lowest_FDR_from2M.apple.TimeSeries'
pomAll<-merge(appledat[,c(1:6,11)],hawdat[,c(1:5,10)])
names(pomAll)[2]<-'flyid'
pomAll$flyid[is.na(pomAll$flyid)]<-'noHit'
ind<-pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 & pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05
ind<-ind & appledat$Lowest_FDR_FruitBetweenMonths < 0.05


doNetStats<-function(ind) {

    mat<-as.matrix(cbind(pomAll[ind,3:6],pomAll[ind,8:11]))
    # for fixing fragged gene models, if desired
    #rbind(
    #    c('FBgn0283451' 'gene14043'), #broad
    #    c('FBgn0000546' 'gene9007'), #Ecr
    #    c('FBgn0004647','gene18485',)) #notch

    # get rid of redundant transcripts that map to the same flybase gene model
    # somewhat conservative, will probably collapse some true gene duplications
    library(WGCNA)
    rowID<-pomAll$gene_id[ind]
    rowGroup<-pomAll$flyid[ind]
    rowGroup[rowGroup=='noHit']<-rowID[rowGroup=='noHit']
    rownames(mat)<-rowID
    a<-collapseRows(mat,rowGroup=rowGroup,rowID=rowID,connectivityBasedCollapsing=T)
    mat<-a$datETcollapsed

    library(coexnet)


    #makeSummarizedExperimentFromDataFrame(mat)

    #findThreshold(mat, method='correlation', plotting = FALSE)
    #[1] 0.94


    thresh=0.94
    softPower = 1;
    adj = adjacency(t(mat), power = softPower);
    diag(adj)<-0
    sum( apply( adj, 2, function(x) sum( x >= thresh ) ) )/2
    #[1] 23755 retained connections
    library(igraph)
    removeEdges<-function(x) {
        x[x < thresh]<-0
        return(x)
    }

    threshAdj<-apply( adj, 1, removeEdges  )
    gr<-graph_from_adjacency_matrix(threshAdj, mode = "undirected",weighted=T)
    degreeNodes<-degree(gr)
    trNodes<-transitivity(gr)
    hubNodes<-hub_score(gr)$vector

    out<-data.frame(id=names(degreeNodes),degree=degreeNodes,hubScores=hubNodes)

    coords<-layout.fruchterman.reingold(gr)
    plot(gr,coords,vertex.size=3,vertex.label=NA,vertex.color='blue')
    return(out)
    
}

#pdf('NetPlotSigTimeAll.pdf')
ind<-pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 & pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05
out<-doNetStats(ind)
#write.table(out,'DEtimeBothRaces.NetworkStats.txt',quote=T,sep="\t",row.names=F)
#dev.off()
alldat<-data.frame(out,cat='all')

#pdf('NetPlotSigTimeSigHost.pdf')
ind<-pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 & pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05
ind<-ind & appledat$Lowest_FDR_FruitBetweenMonths < 0.05
out<-doNetStats(ind)
#dev.off()
alldat<-rbind(alldat,data.frame(out,cat='sigHost'))

#pdf('NetPlotSigTimeRandSamp774.pdf')
set.seed(123)
ind<-sample(1:nrow(pomAll),sum(ind))
out<-doNetStats(ind)
#dev.off()
alldat<-rbind(alldat,data.frame(out,cat='rand'))

library(ggplot2)
p <- ggplot(alldat, aes(x=cat, y=hubScores)) + 
  geom_violin() + ylim(0,0.001)

ggplot(alldat[alldat$cat=='sigHost',], aes(x=hubScores)) +
  geom_histogram(aes(y=..count../sum(..count..))) + ylim(0,0.8)


ggplot(df, aes(height)) + stat_ecdf(geom = "step")

pdf('ecdfNetworkAllvsSigHost.pdf')
ggplot(alldat[alldat$cat!='rand',], aes(hubScores,color=cat)) + 
    stat_ecdf(geom = "step") + theme_classic()
dev.off()

# for fixing fragged gene models, if desired
rbind(
    c('FBgn0283451' 'gene14043'), #broad
    c('FBgn0000546' 'gene9007'), #Ecr
    c('FBgn0004647','gene18485',)) #notch

# get rid of redundant transcripts that map to the same flybase gene model
# somewhat conservative, will probably collapse some true gene duplications
library(WGCNA)
rowID<-pomAll$gene_id[ind]
rowGroup<-pomAll$flyid[ind]
rowGroup[rowGroup=='noHit']<-rowID[rowGroup=='noHit']
rownames(mat)<-rowID
a<-collapseRows(mat,rowGroup=rowGroup,rowID=rowID,connectivityBasedCollapsing=T)
mat<-a$datETcollapsed

library(coexnet)


#makeSummarizedExperimentFromDataFrame(mat)

#findThreshold(mat, method='correlation', plotting = FALSE)
#[1] 0.94


thresh=0.94
softPower = 1;
adj = adjacency(t(mat), power = softPower);
diag(adj)<-0
sum( apply( adj, 2, function(x) sum( x >= thresh ) ) )/2
#[1] 23755 retained connections
library(igraph)
removeEdges<-function(x) {
    x[x < thresh]<-0
    return(x)
}

threshAdj<-apply( adj, 1, removeEdges  )
gr<-graph_from_adjacency_matrix(threshAdj, mode = "undirected",weighted=T)
degreeNodes<-degree(gr)
trNodes<-transitivity(gr)
hubNodes<-hub_score(gr)$vector

out<-data.frame(id=names(degreeNodes),degree=degreeNodes,hubScores=hubNodes)

coords<-layout.fruchterman.reingold(gr)
plot(gr,coords,vertex.size=3,vertex.label=NA)
#write.table(out,'DEtimeBothRaces.NetworkStats.txt',quote=T,sep="\t",row.names=F)

gr<-graph_from_adjacency_matrix(adj, mode = "undirected",weighted=T)
strengthNodes<-strength(gr)

  a b c
a 0 1 1
b 1 0 0
c 1 0 0
mat<-cbind(c(0,1,1),c(1,0,0),c(1,0,0))

colnames(mat)<-c('a','b','c')
rownames(mat)<-c('a','b','c')


library(igraph)

gr<-graph_from_adjacency_matrix(mat, mode = "undirected")

mat<-cbind(c(0,0.8,0.5),c(0.8,0,0),c(0.5,0,0))
colnames(mat)<-c('a','b','c')
rownames(mat)<-c('a','b','c')

gr<-graph_from_adjacency_matrix(mat, mode = "undirected",weighted=T)

plot(gr,edge.width=10*E(gr)$weight)

data(iris)

mat<-abs(cor(iris[,1:4]))
diag(mat)<-0
gr<-graph_from_adjacency_matrix(mat, mode = "undirected",weighted=T)
plot(gr,edge.width=10*E(gr)$weight,layout=layout.fruchterman.reingold)

#Summing up the edge weights of the adjacent edges for each vertex
strength(gr)

biocLite("coexnet")
findThreshold(expData, method, plotting = FALSE)


######## 7/16/2019 -- Analysis of Drosophila genetic interaction networks, but first re-analyze networks form Pomonella expression data #############


#### Reanalyze pomonella networ

ind<-pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 & pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05
out<-doNetStats(ind)
ind2<-appledat$Lowest_FDR_FruitBetweenMonths < 0.05 & ind
ids<-pomAll$flyid
ids[ids=='noHit']<-pomAll$gene_id[ids=='noHit']
geneids<-ids[ind2]

deHostInd<-out$id %in% geneids
deTimeNotHostInd<-out$id %in% ids[ind & !ind2]
median(out$hubScores[deHostInd])
# 0.0002357644
median(out$hubScores[deTimeNotHostInd])
# 0.0002357644
ks.test(out$hubScores[deHostInd],out$hubScores[deTimeNotHostInd])
#	Two-sample Kolmogorov-Smirnov test

#data:  out$hubScores[deHostInd] and out$hubScores[deTimeNotHostInd]
#D = 0.13681, p-value = 1.103e-08
#alternative hypothesis: two-sided
# Note:Essentially same result if you use degree


#### Analyze drosophila network -- seems there's no relationship betwen connectedness estimated from drosophila and divergent transcription

options(stringsAsFactors=F)
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/NetworkAnalysis')
ppi<-read.table('flybase_ppi.txt',header=T,row.names=NULL,sep="\t")




calcDegree<-function(x,y) {
    outid<-character()
    outDegree<-vector()
    all<-unique(c(x,y))
    for (id in all) {
        outid<-c(outid,id)
        outDegree<-c(outDegree, sum(x==id | y==id) )
    }
    return(data.frame(id=outid,degree=outDegree))
}

ppiDegree<-calcDegree(ppi[,1],ppi[,2])
#need to assign degree = 0 to all transcripts de over time in pomonella but not in the dros interaction list
deTimeIds<-pomAll$flyid[ind]
deTimeIds<-deTimeIds[!(deTimeIds=='noHit')]
deHostIds<-pomAll$flyid[ind2]
deHostIds<-deHostIds[!(deHostIds=='noHit')]
zeroIds<-deTimeIds[ !( deTimeIds %in% ppiDegree$id ) ]
zeroIds<-unique(zeroIds)

tempDat<-data.frame(id=zeroIds,degree=0)
ppiDegree<-rbind(ppiDegree,tempDat)
mean(ppiDegree$degree[ ppiDegree$id %in% deHostIds ])
#6.872059
mean(ppiDegree$degree[ !(ppiDegree$id %in% deHostIds) & ppiDegree$id %in% deTimeIds ])
#5.439655
ks.test( ppiDegree$degree[ ppiDegree$id %in% deHostIds ] , ppiDegree$degree[ !(ppiDegree$id %in% deHostIds) & ppiDegree$id %in% deTimeIds ] )
#	Two-sample Kolmogorov-Smirnov test

#data:  ppiDegree$degree[ppiDegree$id %in% deHostIds] and ppiDegree$degree[!(ppiDegree$id %in% deHostIds) & ppiDegree$id %in% ppiDegree$degree[ppiDegree$id %in% deHostIds] and     deTimeIds]
#D = 0.038835, p-value = 0.5153
#alternative hypothesis: two-sided

setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/NetworkAnalysis')
ppi<-read.table('fly_genetic_interactions.txt',header=T,row.names=NULL,sep="\t")
ppiDegree<-calcDegree(ppi[,1],ppi[,2])
#need to assign degree = 0 to all transcripts de over time in pomonella but not in the dros interaction list
deTimeIds<-pomAll$flyid[ind]
deTimeIds<-deTimeIds[!(deTimeIds=='noHit')]
deHostIds<-pomAll$flyid[ind2]
deHostIds<-deHostIds[!(deHostIds=='noHit')]
zeroIds<-deTimeIds[ !( deTimeIds %in% ppiDegree$id ) ]
zeroIds<-unique(zeroIds)

tempDat<-data.frame(id=zeroIds,degree=0)
ppiDegree<-rbind(ppiDegree,tempDat)
mean(ppiDegree$degree[ ppiDegree$id %in% deHostIds ])
#6.872059
mean(ppiDegree$degree[ !(ppiDegree$id %in% deHostIds) & ppiDegree$id %in% deTimeIds ])
#5.439655
ks.test( ppiDegree$degree[ ppiDegree$id %in% deHostIds ] , ppiDegree$degree[ !(ppiDegree$id %in% deHostIds) & ppiDegree$id %in% deTimeIds ] )
#same result as above

# analyze Drosophila transcription factor (tf) data -- also nothing here

setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/NetworkAnalysis')
tf<-read.table('tf_gene.txt',header=T,sep="\t",row.names=NULL,stringsAsFactors=F,quote="\"",comment.char="")

id<-character()
degree<-vector()
for (i in unique(tf[,1])) {
    id<-c(id,i)
    degree<-c(degree, sum( tf[,1] == i ) )
}
inDegree<-data.frame(id,degree)

id<-character()
degree<-vector()
for (i in unique(tf[,2])) {
    id<-c(id,i)
    degree<-c(degree, sum( tf[,2] == i ) )
}
outDegree<-data.frame(id,degree)

ind2<-ind2<-appledat$Lowest_FDR_FruitBetweenMonths < 0.05

allIds<-unique(pomAll$flyid)
allIds<-allIds[!(allIds=='noHit')]
deHostIds<-unique(pomAll$flyid[ind2])
deHostIds<-deHostIds[!(deHostIds=='noHit')]

median( inDegree$degree[ inDegree$id %in% deHostIds ] )
median( inDegree$degree[ !(inDegree$id %in% deHostIds) & inDegree$id %in% allIds ] )

ks.test(inDegree$degree[ inDegree$id %in% deHostIds ],inDegree$degree[ !(inDegree$id %in% deHostIds) & inDegree$id %in% allIds ])
#	Two-sample Kolmogorov-Smirnov test

#data:  inDegree$degree[inDegree$id %in% deHostIds] and inDegree$degree[!(inDegree$id %in% deHostIds) & inDegree$id %in% inDegree$degree[inDegree$id %in% deHostIds] and     allIds]
#D = 0.12452, p-value = 0.9691
