#!/usr/bin/Rscript


args = commandArgs(trailingOnly=TRUE)
#print(paste0('Hello! I am a task number: ', args[1]) )
#install.packages("DEoptim")

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


#########################generate events from the groundtruth intensity function

inhomo_simulation<-function(measure_sup,T,func,c){
  t=0
  points_homo=c()
  points_inhomo=c()
  points_thinned=c()
  while (t<T){
    interval=rexp(1,rate=measure_sup*c)
    points_homo=append(points_homo,interval)
    t=t+interval
    if(t>=T) {break}
    u=runif(1)
    if(u*measure_sup<func(t)){
      points_inhomo=append(points_inhomo,t)
    }
    else{
      points_thinned=append(points_thinned,t)}
    
  }
  return(points_inhomo)
  
}

set.seed(1234)
#########function1 300 datasets


measure_sup1=3 #supreme
T1=50
bin_num1=1000
x1=seq(0,T1,length.out=bin_num1)
intensity1=inten1(x1)

for(c in 1:3){
  for(k in 1:100){
    nam <- paste("points_inhomo_f1_size", c,"data",k, sep = "_")
    assign(nam,inhomo_simulation(measure_sup1,T1, inten1,c))
    
  }
}






#################function2: 300 datasets
measure_sup2=11 #supreme
T2=5
bin_num2=1000
x2=seq(0,T2,length.out=bin_num2)
intensity2=inten2(x2)

for(c in 1:3){
  for(k in 1:100){
    nam <- paste("points_inhomo_f2_size", c,"data",k, sep = "_")
    assign(nam,inhomo_simulation(measure_sup2,T2, inten2,c))
    
  }
}



########################function3 300 datasets
measure_sup3=3.1 #supreme

T3=100
bin_num3=1000
x3=seq(0,T3,length.out=bin_num3)
intensity3=sapply(x3,inten3)

for(c in 1:3){
  for(k in 1:100){
    nam <- paste("points_inhomo_f3_size", c,"data",k, sep = "_")
    assign(nam,inhomo_simulation(measure_sup3,T3, inten3,c))
    
  }
}






#################################estimate hyperparameter


#########################################################################################


nam <- paste0("points_inhomo_f2_size_1_data_",args[1])
points_inhomo=get(nam)
N=length(points_inhomo)





groundest<-function(c,func,T,u1,u2,noise_var){
  #mm is # of grid, usually 100
  #func=inten2
  #mm=10
  #T=T2
  #t=seq(0,T,length.out=mm)
  fv=c*sapply(points_inhomo,func)
  #plot(t,fv,type='l')
  
  
  #fit1=mlegp(t,fv)    
  #summary(fit1)
  #theta1_copy=2*fit1$beta
  #theta0_copy=fit1$sig2
  #print(theta1_copy)
  #print(theta0_copy)
  
  
  
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
    #print(cov_noise)
    #print(dmvnorm(fv, mean=rep(0,mm),sigma = cov_noise,log=TRUE))
    #finalvalue=dtmvnorm(fv, mu=rep(0,mm), sigma=cov_noise, lb=rep(0,mm), ub=rep(Inf,mm), type = c("mc", "qmc"), log = TRUE, B = 1e4)
    finalvalue=dmvnorm(fv, mean=rep(0,N),sigma = cov_noise,log=TRUE)-log(pmvnorm(lower=rep(0,N),
                                                                                 upper=rep(Inf,N),
                                                                                 mean=rep(0,N),
                                                                                 sigma = cov_noise)[1])
    
    return(-finalvalue)
  }
  
  
  out <- DEoptim(fn=kernelloglik, lower = c(0,0), upper = c(u1 ,u2),control=list(trace=TRUE,itermax=3000))
  return(out$optim$bestmem)
}
a=groundest(1,inten2,T2,1000,10,1e-5)
save(a, file =paste0( "/scratch/groups/juliapr/output_Bingjing/groundest/func2_s1_groundest",args[1],".rda"))

#1.088164e+02 1.547089e-10






