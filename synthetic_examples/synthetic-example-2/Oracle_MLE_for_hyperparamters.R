#################################Oracle MLE estimate hyperparameters
###################for the first 100 datasets in the second synthetic example
#!/usr/bin/Rscript
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

#nam <- paste0("data_",args[1])#args[1] range from 301 to 400
#points_inhomo=get(nam)
nam <- paste0("syndata/data_",args[1],".rda")#args[1] range from 301 to 400
load(nam)
points_inhomo=dataa
N=length(points_inhomo)
T2=5
cc=1 #cc=2 corresponds to the second 100 datasets in the second synthetic example

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

a=groundest(cc,inten2,T2,1000,10,1e-5)
save(a, file =paste0( "/groundest/func2_s1_groundest",args[1],".rda"))
