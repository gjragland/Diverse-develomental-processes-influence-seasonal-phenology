#R
#GJR 6/1/2018
#implement LR test of Lynch et al 2014:
# Lynch, Michael, et al. "Population-genetic inference from pooled-sequencing data." Genome biology and evolution 6.5 (2014): 1210-1218.
#modified version to deal with underflow
#uses the log-sum-exp trick, see https://hips.seas.harvard.edu/blog/2013/01/09/computing-log-sum-exp/
#   and http://wittawat.com/posts/log-sum_exp_underflow.html
# this is an approximation, though, and can be worse with smaller numbers
# range of likelhoods is large, so that the constant choice renders really small values to essentially zero

LRfreqTestl<-function(x,pHatE=0.01) {
    nM1<-x[1]
    nm1<-x[2]
    nM2<-x[3]
    nm2<-x[4]
    N<-x[5]
    #default pHatE is Illumina error rate
    #nM1=25 # count of major in pop1
    #nm1=15 # count of minor in pop1

    #nM2=16 # count of major in pop2
    #nm2=25 # count of minor in pop2

    #N=48 # number of individuals sampled in each bulk

    #pHatE<-0.01 #error rate (illumina), from equation 3a
    epHat<-(3 * pHatE)/2 #epsilon hat from 3a


    pHatM1<-nM1/(nM1+nm1) #frequency of major allele based on counts of major and minor
    pHatM2<-nM2/(nM2+nm2) #frequency of major allele based on counts of major and minor
    pHatMpooled<-(nM1+nM2) / ( nM1+nM2+nm1+nm2 ) #pooled estimate

    
    pHatMLE<-function(pHatM,epHat) {
        ( pHatM*( 1 - ( (2*epHat)/3 ) ) - ( epHat/3 ) ) / ( 1 - ( (4*epHat)/3 ) )
    }

    #get MLEs for p in the two populations, and the pooled estimate
    pHatM1MLE<-pHatMLE(pHatM1,epHat)
    pHatM2MLE<-pHatMLE(pHatM2,epHat)
    pHatMpooledMLE<-pHatMLE(pHatMpooled,epHat)



    phiM<-function(p,e) { p*( 1 - (4*e)/3 ) + (e/3) } #probability of observing major allele given p (freq of major allele) and error e
    phim<-function(p,e) { p*( ((4*e)/3) - 1 ) + ( 1 - e ) } #probability of observing minor allele given p (freq of major allele) and error e


    i<-0:(2*N)

    loglikeFun<-function(i,pHatM,nM,nm) {
        #choose((2*N),i) * (pHatM^i) * ( (1 - pHatM)^( (2*N) - i ) ) * ( phiM(( i/(2*N) ),epHat)^nM ) * ( phim(( i/(2*N) ),epHat)^nm )
        #changing function to log space to avoid small numbers
        log(choose((2*N),i)) + log(pHatM^i) + ( (2*N) - i )*log( (1 - pHatM) ) + nM*log( phiM(( i/(2*N) ),epHat)) + nm*log(phim(( i/(2*N) ),epHat))
 
    }

    
    
    sumLikes<-function(i,pHatM,nM,nm) {
        a<-sapply(i,loglikeFun,pHatM=pHatM,nM=nM,nm=nm)
        b<-max(a)
        logsum<-b + sum( sapply( a,function(x) exp(x-b) ) )
    }
    sumLikes(i,pHatM,nM,nm)
    
    LpooledPop1<-sumLikes(i,pHatM=pHatMpooledMLE,nM=nM1,nm=nm1) 
    LpooledPop2<-sumLikes(i,pHatM=pHatMpooledMLE,nM=nM2,nm=nm2)

    LPop1<-sumLikes(i,pHatM=pHatM1MLE,nM=nM1,nm=nm1) 
    LPop2<-sumLikes(i,pHatM=pHatM2MLE,nM=nM2,nm=nm2) 

    LR<-2*( ( LPop1 + LPop2 ) - ( LpooledPop1 + LpooledPop2 ) )
    pchisq(LR,1,lower.tail=F)
}
