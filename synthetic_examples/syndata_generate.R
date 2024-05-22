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
library(mlegp)

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

#########################generate events from the ground truth intensity function
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
#########first synthetic example: generate 300 datasets from orginal, double and triple of inten1
measure_sup1=3 #supreme
T1=50
bin_num1=1000
x1=seq(0,T1,length.out=bin_num1)
intensity1=inten1(x1)
mmm=0
for(c in 1:3){
  for(k in 1:100){
    mmm=mmm+1
    nam <- paste("data_", mmm, sep = "")
    data=inhomo_simulation(measure_sup1,T1, inten1,c)
    save(data, file=paste(nam,'.rda',sep=''))
  }
}

#########second synthetic example: generate 300 datasets from from orginal, double and triple of inten2
measure_sup2=11 #supreme
T2=5
bin_num2=1000
x2=seq(0,T2,length.out=bin_num2)
intensity2=inten2(x2)

for(c in 1:3){
  for(k in 1:100){
    mmm=mmm+1
    nam <- paste("data_", mmm, sep = "")
    data=inhomo_simulation(measure_sup2,T2, inten2,c)
    save(data, file=paste(nam,'.rda',sep=''))
  }
}

