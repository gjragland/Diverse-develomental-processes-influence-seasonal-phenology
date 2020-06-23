#R
# 6/18/2019
#Correlations of allele freqs from RAD data at fennville
#plus overlap of fenville rad data with Urbana pooseq



#### Enrichment tests for ecl loci between Urbana (pool-seq) and Fennville (RAD-seq)
### MBC 6/17/2019

### get the list of RAD SNPs passing HWE at Fennville for apple and haw  

FenEcl <- read.table("/media/raglandlab/ExtraDrive4/mySQL_RAD/Ecl_hwe_RADtoRzeph.txt", header=T)


##This file only has apple and haw freq difs, to get the averaged, we'll query these 
#positions against a larger flat file with the averages 

Ecl_all <- read.table("/media/raglandlab/ExtraDrive4/mySQL_RAD/Ecl_flat.txt", header=T)
Ecl_loc <- matrix(unlist(strsplit(as.character(Ecl_all$RAD_loc),"_")), nrow = nrow(Ecl_all), ncol=2, byrow = T)
Ecl_avg <- data.frame("RAD_contig" = Ecl_loc[,1], "RAD_pos" = Ecl_loc[,2], Ecl_all[,14:15])
Ecl_all <- merge(FenEcl, Ecl_avg, sort=F)

#Now query the Zepheryia positions to get IDs and then get the fisher test results. 

library(RMySQL)

mydb = dbConnect(MySQL(), user='raglandlab', password='pomonella', dbname='PomUrbanaGrant')

for (i in 1:nrow(Ecl_all)){
  dat<-dbGetQuery(mydb, paste("select * from main where scaffold = '",Ecl_all$Zeph_scaf[i],"' and position = '",Ecl_all$Zeph_pos[i],"'",sep=""))
  #dat<-dbGetQuery(mydb, paste("select * from annotation where snpId = '",dat$snpId,"'",sep=""))
  if (i==1) {SNPids<-dat} else {SNPids=rbind(SNPids,dat)}
}

for (i in 1:nrow(SNPids)){
  dat<-dbGetQuery(mydb, paste("select * from snpFisherurbana where snpId = '",SNPids$snpId[i],"'",sep=""))
  #dat<-dbGetQuery(mydb, paste("select * from annotation where snpId = '",dat$snpId,"'",sep=""))
  if (i==1) {UrbFish<-dat} else {UrbFish=rbind(UrbFish,dat)}
}

#get rid of any duplicates (there are occasional positions with multiple SNPids)
ind <- duplicated(UrbFish[,-c(1,27)])
UrbFish <- UrbFish[!ind,]

#Get rid of indels 
gt1<-function(x) {
  out<-F
  if (nchar(x) > 1) {out<-T}
  return(out)
}
isIndel<-sapply(UrbFish$ref,gt1) | sapply(UrbFish$alt,gt1)
UrbFish<-UrbFish[!isIndel,]
#indels<-data[isIndel,]
rm(isIndel)


#convert to numeric 
a<-apply(UrbFish[,6:17],2,as.numeric)
UrbFish[,6:17] <- a
rm(a)

#Get rid of anything with low coverage in the Pool seq 
ind<-rowSums(UrbFish[,6:7]) >= 10 & rowSums(UrbFish[,8:9]) >= 10 & rowSums(UrbFish[,10:11]) >= 10 & 
  rowSums(UrbFish[,12:13]) >= 10 & rowSums(UrbFish[,14:15]) >= 10 & rowSums(UrbFish[,16:17]) >= 10
nrow(UrbFish)
#[1] 8489
sum(ind)
#[1] 8489
#Nothing needs to be filtered on converage. 

#Next, we'll match up the RAD seq and PoolSeq tables so we can run the enrichment tests 

colnames(Ecl_all)[5:6] <- c("scaffold", "position")

Urb_Fen_comp <- merge(Ecl_all, UrbFish, by = c("scaffold", "position"), sort=F)

#blueprint for enrichment tests in Urbana Vs. Fennville 

#         sig Pseq   Nsig Pseq
# sig RAD       a     b
# Nsig RAD      c     d

##Apple 

a <- sum(Urb_Fen_comp$AppleEclPval < 0.05 & Urb_Fen_comp$urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05)
b <- sum(Urb_Fen_comp$AppleEclPval < 0.05 & Urb_Fen_comp$urbana_appleearly_applelate_fisher_pvalue_adjust >= 0.05)
c <- sum(Urb_Fen_comp$AppleEclPval > 0.05 & Urb_Fen_comp$urbana_appleearly_applelate_fisher_pvalue_adjust < 0.05)
d <- sum(Urb_Fen_comp$AppleEclPval > 0.05 & Urb_Fen_comp$urbana_appleearly_applelate_fisher_pvalue_adjust >= 0.05)
apple_fish <- fisher.test(cbind(c(a,c),c(b,d)))

# Fisher's Exact Test for Count Data
# 
# data:  cbind(c(a, c), c(b, d))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  4.113051 6.924410
# sample estimates:
# odds ratio 
#   5.340971 

#Haw Urbana vs. Fennville

#         sig Pseq   Nsig Pseq
# sig RAD       a     b
# Nsig RAD      c     d

a <- sum(Urb_Fen_comp$HawEclPval < 0.05 & Urb_Fen_comp$urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05)
b <- sum(Urb_Fen_comp$HawEclPval < 0.05 & Urb_Fen_comp$urbana_hawearly_hawlate_fisher_pvalue_adjust >= 0.05)
c <- sum(Urb_Fen_comp$HawEclPval > 0.05 & Urb_Fen_comp$urbana_hawearly_hawlate_fisher_pvalue_adjust < 0.05)
d <- sum(Urb_Fen_comp$HawEclPval > 0.05 & Urb_Fen_comp$urbana_hawearly_hawlate_fisher_pvalue_adjust >= 0.05)
haw_fish <- fisher.test(cbind(c(a,c),c(b,d)))

# Fisher's Exact Test for Count Data
# 
# data:  cbind(c(a, c), c(b, d))
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   6.059012 11.001940
# sample estimates:
# odds ratio 
#   8.137841 
