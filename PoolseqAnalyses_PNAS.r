#R
#GJR 12/8/2017
# Analyze RNAseq data
#start by loading results, includes 'appledat' and 'hawdat', which should contain all relevant edgeR results
#also should load the function 'chop'
#may be other functions in there, but we will re-code here to allow modifications

############################################################################################################
################## load data, libraries, themes, etc. ######################################################
############################################################################################################


#load DE output files of edgeR GLMs (see pomonellaEdgeRForPNAS.r)
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/prelimResults')

options(stringsAsFactors=F)
appledat<-read.table('Table.apple.JuneEJD',header=T,row.names=NULL)
hawdat<-read.table('Table.haw.JuneEJD',header=T,row.names=NULL)


#Change to working directory for results
setwd('/media/raglandlab/ExtraDrive1/RpomDiapauseRNAseqTraj_GJR/pomResultsOutliersRemoved')




#load key libraries
library(gplots)
library(RColorBrewer)
library(ggplot2)
library(WGCNA)
#set color palette
#pcolors<-brewer.pal(10,"PiYG")
pcolors<-brewer.pal(10,"BrBG")
options(stringsAsFactors=F)

#set black background theme for presentation figures
theme_black = function() {
 
  theme_classic() %+replace%
 
    theme(
      # Specify axis options
      axis.line = element_line(color="white"),
      panel.background = element_rect(fill = "black"),
      plot.background = element_rect(fill = "black"),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      legend.background = element_rect(color = NA, fill = "black"),
      legend.text = element_text(color = "white")
      #guides(colour = guide_legend(override.aes = list(col="white")))
    )
}

#chopping function for heatmap display
chop<-function(x,thresh=2) {
  elchop<-function(y) {
    if (is.na(y)==FALSE & abs(y) > thresh) {
      if (y > thresh) {y<-thresh}
      if (y < -thresh) {y<--thresh}
    }
    return(y)
  }
  sapply(x,elchop)
}


#create merged expression set for each host race relative to it's own 2-month values
names(appledat)[11]<-'Lowest_FDR_from2M.apple.TimeSeries'
pomAll<-merge(appledat[,c(1:6,11)],hawdat[,c(1:5,10)])
names(pomAll)[2]<-'flyid'
pomAll$flyid[is.na(pomAll$flyid)]<-'noHit'

#7152 transcripts DE over time in apple OR haw flies
sum(pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 | pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05)
#5133 transcripts DE over time in apple 
sum(pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05)
#4597 transcripts DE over time in haw 
sum(pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05)
#2578 DE in both
sum(pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05 & pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05)

ind<-pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05 & pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05
#2578 total transcripts DE across time in both apple and haw
inda<-pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05
indh<-pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05
sum(inda)
sum(inda & ind)
#5133 transcripts DE in apple, 2578 shared with haw
sum(indh)
sum(indh & ind)
#4597 transcripts DE in haw, 2578 shared with apple

#                     sigApple                      nsigApple

#  sigHaw       sum(ind) = 2578               sum(!inda & indh) = 2019 
#  nsigHaw      sum(!indh & inda) = 2555      sum(!indh & !inda)  = 10122 

fisher.test(cbind(c(2578,2555),c(2019,10122)))
#	Fisher's Exact Test for Count Data

#data:  cbind(c(2578, 2555), c(2019, 10122))
#p-value < 2.2e-16
#alternative hypothesis: true odds ratio is not equal to 1
#95 percent confidence interval:
# 4.700695 5.443903
#sample estimates:
#odds ratio 
#  5.057962 


meanAp<-rowMeans(pomAll[3:6])
meanHw<-rowMeans(pomAll[8:11])
ind2<-meanAp*meanHw > 0
sum(ind & ind2)
#2456 of 2578 (95%) regulated over time in same direction in apple and haw
sum(inda & meanAp > 0)/sum(inda)
#54% upregulated in apple
sum(indh & meanHw > 0)/sum(indh)
#53% upregulated in Haw


############################################################################################################
################## Common expression trajectories in apple and haw #########################################
############################################################################################################


#analyze the common trajectories, excluding trancripts likely to have a sig host fruit effect
#indSameTraj<-ind & appledat$Lowest_FDR_FruitBetweenMonths > 0.1
#can also run all same except on all genes w/sig in both (time effect), regardless of host effect
indSameTraj<-ind

#creat heatmap of common trajectory transcripts
#kernal density plot on the colorbar shows values are heavily right skewed; the strongest response is for increased rather than decreased expression
#probably reflects an increase in developmental activity from a developmental nadir reached sometime in the first month or so

mat<-as.matrix(cbind(pomAll[indSameTraj,3:6],pomAll[indSameTraj,8:11]))
choppedData<-apply(mat,2,chop,thresh=2)
clus<-heatmap.2(mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F)

#SigTimeNoHost...
#pdf('heatMapCommonTrajectories.pdf')
#heatmap.2(choppedData,col=pcolors,trace="none",Colv=F,Rowv=clus$rowDendrogram,density.info='density')
#dev.off()

#SigTime...
#pdf('heatMapCommonTrajectoriesSigTime.pdf')
#heatmap.2(choppedData,col=pcolors,trace="none",Colv=F,Rowv=clus$rowDendrogram,density.info='density')
#dev.off()



#### Run WGCNA without the 'W', i.e., no soft thresholding

softPower = 1;
adjacency = adjacency(t(mat), power = softPower);


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#     labels = FALSE, hang = 0.04);


minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)


#pdf('wgcnaCluster.pdf')
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Gene dendrogram and module colors")
#dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(t(mat), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
#plot(METree, main = "Clustering of module eigengenes",
#     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(t(mat), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


modNums<-length(unique(mergedColors))
upGenes<-matrix(nrow=ncol(mat),ncol=modNums)
downGenes<-matrix(nrow=ncol(mat),ncol=modNums)
cn=1
for (j in unique(mergedColors)) {
    for (i in 1:ncol(mat)) {
        submat<-mat[mergedColors==j,i]
        upGenes[i,cn]<-mean(submat[submat > 0])
        downGenes[i,cn]<-mean(submat[submat < 0])
    }
    cn<-cn+1
}
colnames(upGenes)<-unique(mergedColors)

plotTraj<-function(x,y,main) {
    range<-c(min(y,na.rm=T),max(x,na.rm=T))
    #plot(2:6,c(0,x[1:4]),ylim=range,col='green',type='l',main=main) #apple
    plot(2:6,c(0,x[1:4]),ylim=c(-1.1,1.2),col='green',type='l',main=main,las=1) #apple, alt version with set y limits
    lines(2:6,c(0,(y[1:4])),col='green' )
    lines(2:6,c(0,x[5:8]),col='red') #haw
    lines(2:6,c(0,(y[5:8])),col='red')
    lines(2:6,rep(0,5),lty=3)
}
#pdf('WGCNA_modulesSigTimeBothNoHostEffect_MeanTraj.pdf')
#par(mfrow = c(2,3))
#for (i in 1:ncol(upGenes)) {
#    plotTraj(upGenes[,i],downGenes[,i],paste(colnames(upGenes)[i],"_n=",sum(mergedColors==colnames(upGenes)[i]),sep=""))
#}
#dev.off()


library(tidyr)
long <-data.frame(mat,mergedColors) %>% gather(treat, expr, -mergedColors )
long$fruit<-'haw'
long$fruit[grepl('apple',long$treat)]<-'apple'
long$month<-3
long$month[grepl('month4',long$treat)]<-4
long$month[grepl('month5',long$treat)]<-5
long$month[grepl('month6',long$treat)]<-6
long$month<-as.numeric(long$month)
long$expr<-as.numeric(long$expr)
long$direction<-'up'
long$direction[long$expr < 0]<-'down'
for (i in unique(mergedColors)) {
    long<-rbind(long,c(i, 'month2', 0,  'haw',   2,  'up'))
    long<-rbind(long,c(i, 'month2', 0,  'haw',   2,  'down'))
    long<-rbind(long,c(i, 'month2', 0,  'apple',   2,  'up'))
    long<-rbind(long,c(i, 'month2', 0,  'apple',   2,  'down'))
}
long$month<-as.numeric(long$month)
long$expr<-as.numeric(long$expr)

sePlus<-function(x) {
    if (length(x)==1) {x<-c(x,x)}
    mean(x)+2*sd(x)/sqrt(length(x))
}
seMinus<-function(x) {
    if (length(x)==1) {x<-c(x,x)}
    mean(x)-2*sd(x)/sqrt(length(x))
}

#pdf('temp.pdf')
#ggplot(long[long$mergedColors=='black',], aes(x=month,y=expr, color=interaction(fruit,direction))) +
#  stat_summary(geom="ribbon", fun.ymin="seMinus", fun.ymax="sePlus", aes(fill=interaction(fruit,direction)),  alpha=0.3) +
#    theme_classic()

p<-list()
for (i in 1:length(unique(long$mergedColors))) {
    modCol<-unique(long$mergedColors)[i]
    p[[i]]<-ggplot(long[long$mergedColors==modCol,], aes(x=month,y=expr, color=interaction(fruit,direction))) +
  stat_summary(geom="ribbon", fun.ymin="seMinus", fun.ymax="sePlus", aes(fill=interaction(fruit,direction)),  alpha=0.3) +
    theme_classic()+ggtitle(modCol)
}


#pdf('temp.pdf',width=8,height=5)
#grid.arrange(grobs = p)
#dev.off()




#for all genes sig in both                                        
#pdf('WGCNA_modulesSigTimeBoth_MeanTraj.pdf')
#par(mfrow = c(3,3))
#for (i in 1:ncol(upGenes)) {
#    plotTraj(upGenes[,i],downGenes[,i],paste(colnames(upGenes)[i],"_n=",sum(mergedColors==colnames(upGenes)[i]),sep=""))
#}
#dev.off()



out<-pomAll[indSameTraj,]
out$module<-mergedColors
#write.table(out,"WGCNA_modulesSigTimeBothNoHostEffect.txt",sep="\t",row.names=F,quote=F)
###for all genes sig in both  
write.table(out,"WGCNA_modulesSigTimeBoth.txt",sep="\t",row.names=F,quote=F)





###### Cell Fate genes demonstrating pattern of differentiation over time ##################
##plot prospero, numb, brat, mir, and pon
gnames<-c("brat","numb","prospero","mir","pon")
genes<-c("FBgn0010300","FBgn0002973","FBgn0004595","FBgn0021776","FBgn0025739")
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),dir=character(),stringsAsFactors=F)
for (i in 1:length(genes)) {
    y<-colMeans(rbind(as.numeric(pomAll[pomAll$flyid==genes[i],3:6]),as.numeric(pomAll[pomAll$flyid==genes[i],8:11],na.rm=T)))
    dir<-'up'
    if (y[4] < 0) {dir<-'down'}
    subframe<-data.frame(x=x,y=c(0,y),gene=rep(gnames[i],length(x)),dir=dir)
    plotMat<-rbind(plotMat,subframe)
}


#ppcolors<-c(pcolors[length(pcolors)-2],pcolors[3])
#pdf('CellFateGenesCommonTraj.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=gene,pch=gene,col=dir)) + theme_classic() + geom_line(lwd=2) + scale_shape_manual(values=c(15:18,25)) +
    geom_point(size=7,bg=ppcolors[2]) + scale_colour_manual(values=ppcolors) + geom_hline(yintercept = 0,lty=3) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=16)) + labs(x='Weeks Overwinter',y='Relative Expression')
#dev.off()

#### added 8/23/2019 to add ftz-f1 and remove pon, which is not a primary specifier
##plot prospero, numb, brat, mir, and pon
gnames<-c("brat","numb","prospero","mir","ftz-f1")
genes<-c("FBgn0010300","FBgn0002973","FBgn0004595","FBgn0021776","FBgn0001078")
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),dir=character(),stringsAsFactors=F)
for (i in 1:length(genes)) {
    y<-colMeans(rbind(as.numeric(pomAll[pomAll$flyid==genes[i],3:6]),as.numeric(pomAll[pomAll$flyid==genes[i],8:11],na.rm=T)))
    dir<-'up'
    if (y[4] < 0) {dir<-'down'}
    subframe<-data.frame(x=x,y=c(0,y),gene=rep(gnames[i],length(x)),dir=dir)
    plotMat<-rbind(plotMat,subframe)
}


ppcolors<-c(pcolors[length(pcolors)-2],pcolors[3])
pdf('CellFateGenesAndFtzf1CommonTraj.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=gene,pch=gene,col=dir)) + theme_classic() + geom_line(lwd=2) + scale_shape_manual(values=c(15:18,25)) +
    geom_point(size=7,bg=ppcolors[2]) + scale_colour_manual(values=ppcolors) + geom_hline(yintercept = 0,lty=3) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=16)) + labs(x='Weeks Overwinter',y='Relative Expression')
dev.off()




############# Heat Map of Insulin Signaling ###################
insulin<-read.table('/media/raglandlab/ExtraDrive1/flybaseFxnlGeneLists/insulinSignalingFlybase.txt',stringsAsFactors=F,header=T,row.names=NULL,sep="\t")
#chico spans two fragments in assembly, so two transcripts that are really from one gene model (chrom 5)
#will average b/t fragments

mat<-as.matrix(cbind(pomAll[indSameTraj & pomAll$flyid %in% insulin$flyid,3:6],pomAll[indSameTraj & pomAll$flyid %in% insulin$flyid,8:11]))
flyids<-pomAll$flyid[indSameTraj & pomAll$flyid %in% insulin$flyid]
rownames(mat)<-insulin$SYMBOL[match(flyids,insulin$flyid)]
choppedData<-apply(mat,2,chop,thresh=2)
clus<-heatmap.2(mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F)

#pdf('heatMapInsulinSigTime.pdf')
#heatmap.2(choppedData,col=pcolors,trace="none",Colv=F,Rowv=clus$rowDendrogram,density.info='density')
#dev.off()

                                        # + scale_colour_gradientn(colors = brewer.pal(10,"PiYG"))

######### Trajectories of Wnt signalling, filter for sig time, average b/t appl/haw, then average over multiple gene models #########


meanTraj<-function(x) {
    y<-rbind(x[1:4],x[5:8])
    colMeans(y)
}

traj<-t(apply(mat,1,meanTraj))
maxFold<-apply(traj,1,function(x) abs(max(x)))
traj<-traj[maxFold > 1,]
#average multi-models for genes
ind<-rownames(traj)=='nej'
a<-colMeans(traj[ind,])
traj<-traj[!ind,]
traj<-rbind(traj,a)
rownames(traj)[nrow(traj)]<-'nej'

aveRows<-function(x,traj) {
    ind<-rownames(traj)==x
    a<-colMeans(traj[ind,])
    traj<-traj[!ind,]
    traj<-rbind(traj,a)
    rownames(traj)[nrow(traj)]<-x
    return(traj)
}
traj<-aveRows('nej',traj)
traj<-aveRows('Apc',traj)
traj<-aveRows('fz2',traj)

traj<-data.frame(traj)
traj$M0<-0
names(traj)<-c('M3','M4','M5','M6','M2')

library(tidyr)
traj$gname<-rownames(traj)
plotDat<-gather(traj, key = month, value = measurement,
       M2, M3, M4, M5, M6)
plotDat$month.num<-sapply(plotDat$month,function(x) as.numeric(unlist(strsplit(x,'M')[1])[2]))

#pdf('WntGenesCommonTraj.pdf')
#ppcolors<-c(pcolors[length(pcolors)-2],pcolors[3])
#ggplot(data=plotDat,aes(x=month.num,y=measurement,group=gname,pch=gname)) + theme_classic() + geom_line(lwd=2) + scale_shape_manual(values=c(0:2,5:6,15:18,25)) +
#    geom_point(size=7,bg=ppcolors[2]) + scale_colour_manual(values=ppcolors) + geom_hline(yintercept = 0,lty=3) +
#    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=16)) + labs(x='Months Overwinter',y='Relative Expression')
#dev.off()

############# Heat Map of Wnt Signaling ###################
wnt<-read.table('/media/raglandlab/ExtraDrive1/flybaseFxnlGeneLists/wntSignalingFlybase.txt',stringsAsFactors=F,header=T,row.names=NULL,sep="\t")
mat<-as.matrix(cbind(pomAll[indSameTraj & pomAll$flyid %in% wnt$flyid,3:6],pomAll[indSameTraj & pomAll$flyid %in% wnt$flyid,8:11]))
flyids<-pomAll$flyid[indSameTraj & pomAll$flyid %in% wnt$flyid]
rownames(mat)<-wnt$SYMBOL[match(flyids,wnt$flyid)]
choppedData<-apply(mat,2,chop,thresh=2)
clus<-heatmap.2(mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F)

#pdf('heatMapWntSigTime.pdf')
#heatmap.2(choppedData,col=pcolors,trace="none",Colv=F,Rowv=clus$rowDendrogram,density.info='density')
#dev.off()



##################### distribution of hub scores, as estimated in 'AnalyzeTrajectoriesNetwork_PNAS.r' ##############
# NOTE: All main analyses of network connectivity/hubness in 'AnalyzeTrajectoriesNetwork_PNAS.r'

netStats<-read.table('DEtimeBothRaces.NetworkStats.txt',sep="\t",row.names=NULL,stringsAsFactors=F,header=T)
dat<-unique(appledat$FlyBase_FBgn_tophit.x.x[indSameTraj])
dat<-data.frame(id=dat,hubscore=0)
dat<-dat[!is.na(dat$id),]
dat<-merge( netStats,dat,all=T)
dat$hubScore[is.na(dat$hubScore)]<-0
#distribution of hub scores is roughly exponential
pdf('HistDETimeHubscores.pdf')
hist(dat$hubScore,xlab = 'Hub Score', main='Distribution of Hub Scores for all trancripts DE across time')
dev.off()




############################################################################################################
################## Host Race difference in trajectories ####################################################
############################################################################################################


inda<-pomAll$Lowest_FDR_from2M.apple.TimeSeries < 0.05
indh<-pomAll$Lowest_FDR_from.TimeSeries.haw < 0.05

fdrsTime<-pomAll$Lowest_FDR_from2M.apple.TimeSeries
     
#create merged expression set for each host race relative to HAW 2-month values
names(appledat)[11]<-'Lowest_FDR_from2M.apple.TimeSeries'
pomAll<-merge(appledat[,c(1,12:16,45)],hawdat[,c(1:5,34)])
ind<-pomAll$Lowest_FDR_FruitBetweenMonths < 0.05
names(pomAll)[12]<-'flyid'
pomAll$flyid[is.na(pomAll$flyid)]<-'noHit'


ind<-pomAll$Lowest_FDR_FruitBetweenMonths < 0.05
#2902 transcripts DE between host races (at at least one time point)

sum(appledat$ApplevsHaw2Month_FDR < 0.05)
#1134 gene DE at 2 months
sum(appledat$ApplevsHaw3Month_FDR < 0.05)
#45 gene DE at 3 months
sum(appledat$ApplevsHaw4Month_FDR < 0.05)
#6 gene DE at 4 months
sum(appledat$ApplevsHaw5Month_FDR < 0.05)
#1755 gene DE at 5 months
sum(appledat$ApplevsHaw6Month_FDR < 0.05)
#329 gene DE at 6 months

sum(ind & (inda | indh))
#2160 of the 2902 (74%) are significantly DE over time



#NADH ubiquinone oxidoreductase
a<-cbind(
c('FBgn0031684','FBgn0047038','FBgn0031228','FBgn0031021','FBgn0035046','FBgn0030718'),
c('ND-13A','ND-13B','ND-15','ND-18','ND-19','ND-20'))


#Succinate dehydrogenase
b<-cbind(c('FBgn0037873','FBgn0039112'),
c('SdhC','SdhD'))

#Ubiquinol-cytochrome c oxidoreductase
c<-cbind(c('FBgn0034245','FBgn0030733','FBgn0035600','FBgn0038271','FBgn0036728'),
c('UQCR-6.4','UQCR-14','Cyt-c1','UQCR-C1','UQCR-Q'))

#Cytochrome c oxidase
d<-cbind(c('FBgn0034877','FBgn0015031','FBgn0040773'),
c('levy','cype','COX7C'))

#F0/F1 ATP synthase
e<-cbind(c('FBgn0020235','FBgn0028342','FBgn0019644','FBgn0016120','FBgn0038224','FBgn0035032','FBgn0010612'),
c('ATPsynγ','ATPsynδ','ATPsynB','ATPsynD','ATPsynE','ATPsynF','ATPsynG'))

oxphos<-data.frame(rbind(a,b,c,d,e),stringsAsFactors=F)
names(oxphos)<-c('flyid','SYMBOL')


#set apple and haw colors for plots (for black background)
ppcolors<-c('#00A29F','#B84025')

##plot NADH ubiquinone oxidoreductase
gnames<-c('ND-13A','ND-13B','ND-15','ND-18','ND-19','ND-20')
genes<-c('FBgn0031684','FBgn0047038','FBgn0031228','FBgn0031021','FBgn0035046','FBgn0030718')
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:(length(genes)-1)) { #exclude last gene for visual purposes
    ya<-as.numeric(pomAll[pomAll$flyid==genes[i],2:6])
    yh<-as.numeric(pomAll[pomAll$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}


#ppcolors<-c(rgb(red=80,green=86,blue=113,maxColorValue=255),rgb(red=189,green=96,blue=130,maxColorValue=255))
#pdf('NADHubiquinoneOxidoreductase.blackback.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),pch=gene,col=pop)) + theme_black() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="white") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month')
#dev.off()

#Succinate dehydrogenase
gnames<-c('SdhC','SdhD')
genes<-c('FBgn0037873','FBgn0039112')
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:length(genes)) { 
    ya<-as.numeric(pomAll[pomAll$flyid==genes[i],2:6])
    yh<-as.numeric(pomAll[pomAll$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}

#ppcolors<-c(rgb(red=80,green=86,blue=113,maxColorValue=255),rgb(red=189,green=96,blue=130,maxColorValue=255))
#pdf('SuccinateDehydrogenase.blackback.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),pch=gene,col=pop)) + theme_black() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="white") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month')
#dev.off()


#Ubiquinol-cytochrome c oxidoreductase
gnames<-c('UQCR-6.4','UQCR-14','Cyt-c1','UQCR-C1','UQCR-Q')
genes<-c('FBgn0034245','FBgn0030733','FBgn0035600','FBgn0038271','FBgn0036728')
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:length(genes)) { 
    ya<-as.numeric(pomAll[pomAll$flyid==genes[i],2:6])
    yh<-as.numeric(pomAll[pomAll$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}

ppcolors<-c(rgb(red=80,green=86,blue=113,maxColorValue=255),rgb(red=189,green=96,blue=130,maxColorValue=255))
#pdf('Ubiquinol-cytochromeCoxidoreductase.blackback.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),pch=gene,col=pop)) + theme_black() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="white") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month')
#dev.off()
#pdf('Ubiquinol-cytochromeCoxidoreductase.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),pch=gene,col=pop)) + theme_classic() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="black") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month') + ylim(-.8,.35)
#dev.off()

#Cytochrome c oxidase
gnames<-c('levy','cype','COX7C')
genes<-c('FBgn0034877','FBgn0015031','FBgn0040773')
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:length(genes)) { 
    ya<-as.numeric(pomAll[pomAll$flyid==genes[i],2:6])
    yh<-as.numeric(pomAll[pomAll$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}

#ppcolors<-c(rgb(red=80,green=86,blue=113,maxColorValue=255),rgb(red=189,green=96,blue=130,maxColorValue=255))
#ppcolors<-c(rgb(red=184,green=65,blue=38,maxColorValue=255),rgb(red=51,green=51,blue=51,maxColorValue=255))
#ppcolors<-c('#B84025','#00A29F')
#pdf('CytochromeCoxidase.blackback.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),pch=gene,col=pop)) + theme_black() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="white") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month')
#dev.off()

#F0/F1 ATP synthase
gnames<-c('ATPsynγ','ATPsynδ','ATPsynB','ATPsynD','ATPsynE','ATPsynF','ATPsynG')
genes<-c('FBgn0020235','FBgn0028342','FBgn0019644','FBgn0016120','FBgn0038224','FBgn0035032','FBgn0010612')
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:length(genes)) { #exclude first two for graphical purposes 
    ya<-as.numeric(pomAll[pomAll$flyid==genes[i],2:6])
    yh<-as.numeric(pomAll[pomAll$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}


#pdf('F0-F1_ATP_synthase.blackback.pdf')
#ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),pch=gene,col=pop)) + theme_black() + geom_line(lwd=2) + 
#    geom_hline(yintercept = 0,lty=3,col="white") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
#    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month')
#dev.off()
pdf('F0-F1_ATP_synthase.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),pch=gene,col=pop)) + theme_classic() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="black") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25,0:1)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month') + ylim(-.8,.35)
dev.off()
####### Rpls ###########

rpls<-read.table('/media/raglandlab/ExtraDrive1/flybaseFxnlGeneLists/RplsDEbtRaces.txt',header=F,row.names=NULL)
dat<-pomAll[ind,]
dups<-dat$flyid[duplicated(dat$flyid)]
dat<-dat[!(dat$flyid %in% dups),]

gnames<-as.vector(rpls[,2])
genes<-as.vector(rpls[,1])
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:length(genes)) { 
    ya<-as.numeric(dat[dat$flyid==genes[i],2:6])
    yh<-as.numeric(dat[dat$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}


#pdf('RibosomalProteins.blackback.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),col=pop)) + theme_black() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="white") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month')
#dev.off()



#separate plots for large and small subunit
gnames<-as.vector(rpls[,2])[grepl('RpL',rpls[,2])]
genes<-as.vector(rpls[,1])[grepl('RpL',rpls[,2])]
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:length(genes)) { 
    ya<-as.numeric(dat[dat$flyid==genes[i],2:6])
    yh<-as.numeric(dat[dat$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}
large<-plotMat

#pdf('RibosomalProteinsLargeSubunit.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),col=pop)) + theme_classic() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="black") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month') + ylim(-.9,.2)
#dev.off()

gnames<-as.vector(rpls[,2])[grepl('RpS',rpls[,2])]
genes<-as.vector(rpls[,1])[grepl('RpS',rpls[,2])]
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),stringsAsFactors=F)
for (i in 1:length(genes)) { 
    ya<-as.numeric(dat[dat$flyid==genes[i],2:6])
    yh<-as.numeric(dat[dat$flyid==genes[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)))
    plotMat<-rbind(plotMat,subframe)
}
small<-plotMat
large<-data.frame(large,unit='large')
small<-data.frame(small,unit='small')



#ppcolors<-c(rgb(red=80,green=86,blue=113,maxColorValue=255),rgb(red=189,green=96,blue=130,maxColorValue=255))
#pdf('RibosomalProteinsSmallSubunit.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),col=pop)) + theme_classic() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="black") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month') + ylim(-.9,.2)
#dev.off()


plotMat<-rbind(large,small)
#ppcolors<-c(rgb(red=80,green=86,blue=113,maxColorValue=255),rgb(red=189,green=96,blue=130,maxColorValue=255))
pdf('RibosomalProteinsAll.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gene,pop),col=pop,pch=unit)) + theme_classic() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="black") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month') 
dev.off()





###### Histone proteins ##################
his<-read.table('allFlyIds.withNames.txt',header=T,sep="\t",row.names=NULL,stringsAsFactors=F,quote="\"",comment.char="")
#greping for His2,His3, and His4 reveals that the rhago histone genes all annotate to only one drosophila copy of each:
#FBgn0053865 His2A:CG33865
#FBgn0053910 His2B:CG33910
#FBgn0053866 His3:CG33866
#FBgn0053909 His4:CG33909
#https://onlinelibrary.wiley.com/doi/full/10.1002/prot.21720

gnames<-c('His2A','His2A','His2A',
          'His2B','His2B','His2B',
          'His3','His3',
          'His4','His4','His4')
ids<-c('gene1175','gene1178','gene19717',
       'gene1176','gene1179','gene19716',
       'gene1174','gene7382',
       'gene1173','gene1177','gene9918')
genes<-c('His2A.1','His2A.2','His2A.3',
          'His2B.1','His2B.2','His2B.3',
          'His3.1','His3.2',
          'His4.1','His4.2','His4.3')
x<-c(2:6)
plotMat<-data.frame(x=numeric(),y=numeric(),gene=character(),pop=character(),gname=character(),stringsAsFactors=F)
for (i in 1:length(ids)) {  
    ya<-as.numeric(pomAll[pomAll$gene_id==ids[i],2:6])
    yh<-as.numeric(pomAll[pomAll$gene_id==ids[i],8:11])
    yh<-c(0,yh)
    subframe<-data.frame(x=x,y=c(ya,yh),gene=rep(gnames[i],length(x)),pop=c(rep('apple',5),rep('haw',5)),gname=genes[i])
    plotMat<-rbind(plotMat,subframe)
}
pdf('HistoneProteins.pdf')
ggplot(data=plotMat,aes(x=x,y=y,group=interaction(gname,pop),pch=gene,col=pop)) + theme_classic() + geom_line(lwd=2) + 
    geom_hline(yintercept = 0,lty=3,col="black") +  geom_point(size=7) + scale_colour_manual(values=ppcolors) + scale_shape_manual(values=c(15:18,25)) +
    theme(legend.position=c(.1, .7)) + theme(axis.text = element_text(size=28)) + labs(x='Weeks Overwinter',y='Expression relative to Haw 2-month')
dev.off()

ind<-pomAll$Lowest_FDR_FruitBetweenMonths < 0.05
mat<-as.matrix(cbind(pomAll[ind,2:6],pomAll[ind,8:11]))
choppedData<-apply(mat,2,chop,thresh=2)
clus<-heatmap.2(mat,col=brewer.pal(10,"PiYG"),trace="none",Colv=F)

#Trajectories for significant host differences in at least one month
pdf('heatMapDEHostTrajectories.pdf')
heatmap.2(choppedData,col=pcolors,trace="none",Colv=F,Rowv=clus$rowDendrogram,density.info='density')
dev.off()


#sublists
i2<-appledat$ApplevsHaw2Month_FDR < 0.05
i=nrow(appledat)
demonth<-matrix(rep(0,(i*5)),nrow=i,ncol=5)
for (i in 40:44) {
    j=i-39
    deind<-appledat[,i] < 0.05
    demonth[deind,j]<-j+1
}



#### Run WGCNA without the 'W', i.e., no soft thresholding

softPower = 1;
adjacency = adjacency(t(mat), power = softPower);


# Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
#plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#     labels = FALSE, hang = 0.04);


minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                deepSplit = 2, pamRespectsDendro = FALSE,
                minClusterSize = minModuleSize);
table(dynamicMods)


#pdf('wgcnaCluster.pdf')
# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
#plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
#                    dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05,
#                    main = "Gene dendrogram and module colors")
#dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(t(mat), colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
#plot(METree, main = "Clustering of module eigengenes",
#     xlab = "", sub = "")

MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(t(mat), dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


modNums<-length(unique(mergedColors))
upGenes<-matrix(nrow=ncol(mat),ncol=modNums)
downGenes<-matrix(nrow=ncol(mat),ncol=modNums)
cn=1
for (j in unique(mergedColors)) {
    for (i in 1:ncol(mat)) {
        submat<-mat[mergedColors==j,i]
        upGenes[i,cn]<-mean(submat[submat > 0])
        downGenes[i,cn]<-mean(submat[submat < 0])
    }
    cn<-cn+1
}
colnames(upGenes)<-unique(mergedColors)

plotTraj<-function(x,y,main) {
    range<-c(min(y,na.rm=T),max(x,na.rm=T))
    #plot(2:6,c(0,x[1:4]),ylim=range,col='green',type='l',main=main) #apple
    plot(2:6,x[1:5],ylim=c(-1.1,1.2),col='green',type='l',main=main,las=1) #apple, alt version with set y limits
    lines(2:6,y[1:5],col='green' )
    lines(2:6,c(0,x[6:9]),col='red') #haw
    lines(2:6,c(0,(y[6:9])),col='red')
    lines(2:6,rep(0,5),lty=3)
}
#pdf('WGCNA_modulesSigHost_MeanTraj.pdf')
par(mfrow = c(3,3))
for (i in 1:ncol(upGenes)) {
    plotTraj(upGenes[,i],downGenes[,i],paste(colnames(upGenes)[i],"_n=",sum(mergedColors==colnames(upGenes)[i]),sep=""))
}
#dev.off()


library(tidyr)
long <-data.frame(mat,mergedColors) %>% gather(treat, expr, -mergedColors )
long$fruit<-'haw'
long$fruit[grepl('Apple',long$treat)]<-'apple'
long$month<-3
long$month[grepl('month4',long$treat)]<-4
long$month[grepl('month5',long$treat)]<-5
long$month[grepl('month6',long$treat)]<-6
long$month[grepl('2month',long$treat)]<-2
long$month[grepl('4month',long$treat)]<-4
long$month[grepl('5month',long$treat)]<-5
long$month[grepl('6month',long$treat)]<-6
long$month<-as.numeric(long$month)
long$expr<-as.numeric(long$expr)
long$direction<-'up'
long$direction[long$expr < 0]<-'down'
for (i in unique(mergedColors)) {
    long<-rbind(long,c(i, 'month2', 0,  'haw',   2,  'up'))
    long<-rbind(long,c(i, 'month2', 0,  'haw',   2,  'down'))
}
long$month<-as.numeric(long$month)
long$expr<-as.numeric(long$expr)

sePlus<-function(x) {
    if (length(x)==1) {x<-c(x,x)}
    mean(x)+2*sd(x)/sqrt(length(x))
}
seMinus<-function(x) {
    if (length(x)==1) {x<-c(x,x)}
    mean(x)-2*sd(x)/sqrt(length(x))
}

#pdf('temp.pdf')
#ggplot(long[long$mergedColors=='black',], aes(x=month,y=expr, color=interaction(fruit,direction))) +
#  stat_summary(geom="ribbon", fun.ymin="seMinus", fun.ymax="sePlus", aes(fill=interaction(fruit,direction)),  alpha=0.3) +
#    theme_classic()

p<-list()
for (i in 1:length(unique(long$mergedColors))) {
    modCol<-unique(long$mergedColors)[i]
    p[[i]]<-ggplot(long[long$mergedColors==modCol,], aes(x=month,y=expr, color=interaction(fruit,direction))) +
  stat_summary(geom="ribbon", fun.ymin="seMinus", fun.ymax="sePlus", aes(fill=interaction(fruit,direction)),  alpha=0.3) +
    theme_classic()+ggtitle(modCol)+ylim(-1,1)
}
library(gridExtra)


pdf('WGCNA_modulesSigHost_RibbonTraj.pdf',width=15.5,height=8)
grid.arrange(grobs = p)
dev.off()



out<-pomAll
out<-data.frame(out,demonth)
out<-out[ind,]
out$module<-mergedColors
names(out)[13:17]<-c('deM2','deM3','deM4','deM5','deM6')
write.table(out,"WGCNA_modulesSigHostEffectPlusDeMonths.txt",sep="\t",row.names=F,quote=F)


}


