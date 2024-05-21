#################################routine MLE estimate hyperparameters
###################for the first 100 datasets in the first synthetic example

#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

library(parallel)
library(DEoptim)
library(pracma)
library(TruncatedNormal)
library(matrixcalc) # is.positive.definite
library(matlib)
library(Matrix)
library(base)
library(mvtnorm)
library(truncnorm)

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

measure_sup1=3 #supreme
T1=50
#bin_num1=1000
#x1=seq(0,T1,length.out=bin_num1)
#intensity1=inten1(x1)

measure_sup2=11 #supreme
T2=5
#bin_num2=1000
#x2=seq(0,T2,length.out=bin_num2)
#intensity2=inten2(x2)

measure_sup3=3.1 #supreme
T3=100
#bin_num3=1000
#x3=seq(0,T3,length.out=bin_num3)
#intensity3=sapply(x3,inten3)


#########################################################################################
nam <- paste0("data_",args[1])#args[1] range from 1 to 100
points_inhomo=get(nam)
N=length(points_inhomo)


groundest<-function(c,func,T,u1,u2,noise_var){
  fv=c*sapply(points_inhomo,func)
  kernelloglik<-function(par){
    theta00=par[1]
    theta11=par[2]
    cov=matrix(0,N,N)
    for (i in 1:N)
    {
      for (j in i:N){
        cov[i,j]=expo_quad_kernel(theta00,theta11,points_inhomo[i],points_inhomo[j])
        cov[j,i]=cov[i,j]
      }
    }
    
    cov_noise=cov+diag(N)*noise_var
    finalvalue=dmvnorm(fv, mean=rep(0,N),sigma = cov_noise,log=TRUE)-log(pmvnorm(lower=rep(0,N),
                                                                                 upper=rep(Inf,N),
                                                                                 mean=rep(0,N),
                                                                                 sigma = cov_noise)[1])
    return(-finalvalue)
  }
  
  out <- DEoptim(fn=kernelloglik, lower = c(0,0), upper = c(u1 ,u2),control=list(trace=TRUE,itermax=3000))
  return(out$optim$bestmem)
}

cc=1
a=groundest(cc,inten1,T1,20,1,1e-5)
save(a, file =paste0( "/groundest/func1_s1_groundest",args[1],".rda"))









