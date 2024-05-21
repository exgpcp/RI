#!/usr/bin/Rscript
#####################weighted MAP##############################################
###################for the 300 datasets in the first synethetic example######################
args = commandArgs(trailingOnly=TRUE)
library(parallel)
library(DEoptim)
library(pracma)
library(TruncatedNormal)
library(matrixcalc) # is.positive.definite
#library(matlib)
library(Matrix)
library(base)
library(mvtnorm)
library(truncnorm)
#library(mlegp)

expo_quad_kernel<-function(theta00,theta11,xn,xm){ # 1,0.1
  return(theta00*exp(-theta11/2*sum((xn - xm)**2)))
}

expo_quad_kernel2<-function(theta00,theta11,xn,T){ # 1,0.1
  return( sqrt(pi/2/theta11)*theta00*(erf(sqrt(theta11/2)*(T-xn))+erf(sqrt(theta11/2)*xn)))
}

expo_quad_kernel3<-function(theta00,theta11,T){ # 1,0.1
  return( 2*theta00/theta11*(sqrt(pi*theta11/2)*T*erf(sqrt(theta11/2)*T)+ exp(-theta11/2*(T**2)) -1))
}

inten1<-function(x){
  return(2*exp(-x/15)+exp(-((x-25)/10)**2))
}

inten2<-function(x){
  return(10+x-x)
}

inten3<-function(x){
  if(x<=25){return(2+1/25*x)}
  if(x<=50){return(5-2/25*x)}
  if(x<=75){return(3/50*x-2)}
  return(1/50*x+1)
  
}

#measure_sup1=3 #supreme
T1=50
#bin_num1=1000
#x1=seq(0,T1,length.out=bin_num1)
#intensity1=inten1(x1)

#measure_sup2=11 #supreme
T2=5
#bin_num2=1000
#x2=seq(0,T2,length.out=bin_num2)
#intensity2=inten2(x2)

#measure_sup3=3.1 #supreme
T3=100
#bin_num3=1000
#x3=seq(0,T3,length.out=bin_num3)
#intensity3=sapply(x3,inten3)


nam <- paste0("data_",args[1])
points_inhomo=get(nam)
noise_var=1e-5
ccc=4
N=length(points_inhomo)
T=T1
results=matrix(0,10,3)
for (mm in 1:10){
  delta_m=T/(mm-1)
  t=seq(0,T,length.out=mm)
  negloglikelihood31<-function(samps,par){
    theta00=par[1]
    theta11=par[2]
    m=length(samps)
    t=seq(0,T,length.out=m)
    delta_m=T/(m-1)
    c_vec=rep(delta_m,m)
    c_vec[1]=delta_m/2
    c_vec[m]=delta_m/2
    comp1=-sum(c_vec*samps)
    samp_event=rep(0,N)
    fin<-function(x){
      index_L=floor(x/delta_m)+1
      index_R=ceiling(x/delta_m)+1
      if (x-t[index_L]<=delta_m/2){
        fin=index_L
      }
      else{
        fin=index_R
      }
      return(samps[fin])
    }
    samp_event=sapply(points_inhomo,fin)
    comp_total=comp1+sum(log(samp_event))
    mmm=m+1
    cov=matrix(0,mmm,mmm)
    sampinteg=c(samps,-comp1)
    for (i in 1:mmm)
    {
      for (j in i:mmm){
        if (j<mmm & i<mmm){
          cov[i,j]=expo_quad_kernel(theta00,theta11,t[i],t[j])
        }
        if(j==mmm&i<mmm){
          cov[i,j]=expo_quad_kernel2(theta00,theta11,t[i],T)
        }
        if(j==mmm & i==mmm){
          cov[i,j]=expo_quad_kernel3(theta00,theta11,T)  
        }
        if(j!=i){
          cov[j,i]=cov[i,j]
        }
      }
    }
    cov_noise=cov+diag(mmm)*1e-5
    finalvalue=dmvnorm(sampinteg, mean=rep(0,mmm),sigma = cov_noise,log=TRUE)-log(pmvnorm(lower=rep(0,mmm),
                                                                                          upper=rep(Inf,mmm),
                                                                                          mean=rep(0,mmm),
                                                                                          sigma = cov_noise)[1])
    
    if(is.positive.definite(cov_noise)==FALSE){finalvalue=-1e10}
    return(-comp_total-finalvalue/ccc)
  }



hypeasy<-function(par){
  return(negloglikelihood31(par[-c(1,2)],par[c(1,2)]))
}

outnew <- DEoptim(fn=hypeasy, lower = rep(0,2+mm), upper = c(30,1,rep(10,mm)),control=list(trace=TRUE,itermax=1000))#,fnMap=Mapfun)
results[mm,3]=outnew$optim$bestval
results[mm,1:2]=outnew$optim$bestmem[1:2]
}

###############################################################################
index=which.min(results[,3])
a=results[index,1:2]
b=c(a,index)
save(b, file =paste0( "/MAP/func1_priorest",args[1],".rda"))








