#INLA COX Process
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
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


##############################for synthetic example 1
coverage1_1=rep(0,300)
coverage1_2=rep(0,300)
l2_dist1_1=rep(0,300)
l2_dist1_2=rep(0,300)
width1_1=rep(0,300)
width1_2=rep(0,300)
time_1=rep(0,300)

for (jjj in 1:300){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  ccc=ceil(jjj/100)
  T=50
  Ngrid=100
  hyper <- list(prec = list(param = c(.1, .1)))
  grid<-seq(0,T,T/Ngrid)
  #grid<-seq(0,5,by=1)
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  m <- length(dataa)
  grid_trimmed <- setdiff(x = grid, y = dataa)
  sorting <- sort(c(grid_trimmed, dataa), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE)
  time <- field[time_index]
  event_out <- c(rep(0, length(grid_trimmed)), rep(1,m))[ordering]
  E <- diff(sgrid)
  E_log = log(E)
  
  data_run <-data.frame(y =  event_out[-1], time = time, E_log = E_log, E= E)
  formula <- y ~  f(time, model="rw1", hyper = hyper, constr = FALSE, scale.model = TRUE)-1
  family <- "poisson"
  mod <- INLA::inla(formula, family = "poisson", data = data_run, offset = data_run$E_log,
                    control.predictor = list(compute=TRUE))
  
  ff<-function(x){
    return(exp(mod$summary.random$time$`0.5quant`)[1+sum(x>=0.25+mod$summary.random$time$ID)])
  }
  
  ff2<-function(x){
    return(exp(mod$summary.random$time$`0.025quant`)[1+sum(x>=0.25+mod$summary.random$time$ID)])
  }
  
  ff3<-function(x){
    return(exp(mod$summary.random$time$`0.975quant`)[1+sum(x>=0.25+mod$summary.random$time$ID)])
  }
  bins=100
  #x1=seq(0,T1,T1/bins)[1:bins]
  x1=dataa
  intensity1=sapply(x1,inten1)*ccc
  med=sapply(x1,ff)
  low=sapply(x1,ff2)
  high=sapply(x1,ff3)
  width1=sum(high-low)/length(dataa)
  
  coverage=sum((intensity1>=low)*(intensity1<=high))/length(dataa)
  l2dist_med=sum((intensity1-med)**2)
  coverage1_1[jjj]=coverage
  l2_dist1_1[jjj]=l2dist_med
  width1_1[jjj]=width1
  
  bins=100
  x1=seq(0,T1,T1/bins)[1:bins]
  intensity1=sapply(x1,inten1)*ccc
  med=sapply(x1,ff)
  low=sapply(x1,ff2)
  high=sapply(x1,ff3)
  width2=sum(high-low)/bins
  
  coverage=sum((intensity1>=low)*(intensity1<=high))/bins
  l2dist_med=sum((intensity1-med)**2)
  coverage1_2[jjj]=coverage
  l2_dist1_2[jjj]=l2dist_med
  width1_2[jjj]=width2
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_1[jjj]=time.taken
}


##############################for synthetic example 2
coverage2_1=rep(0,300)
coverage2_2=rep(0,300)
l2_dist2_1=rep(0,300)
l2_dist2_2=rep(0,300)
width2_1=rep(0,300)
width2_2=rep(0,300)
time_2=rep(0,300)

for (jjj in 1:300){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',(jjj+300),'.rda',sep=''))
  ccc=ceil(jjj/100)
  T=5
  Ngrid=100
  hyper <- list(prec = list(param = c(.1, .1)))
  grid<-seq(0,T,T/Ngrid)
  #grid<-seq(0,5,by=1)
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  m <- length(dataa)
  grid_trimmed <- setdiff(x = grid, y = dataa)
  sorting <- sort(c(grid_trimmed, dataa), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE)
  time <- field[time_index]
  event_out <- c(rep(0, length(grid_trimmed)), rep(1,m))[ordering]
  E <- diff(sgrid)
  E_log = log(E)
  data_run <-data.frame(y =  event_out[-1], time = time, E_log = E_log, E= E)
  formula <- y ~  f(time, model="rw1", hyper = hyper, constr = FALSE, scale.model = TRUE)-1
  family <- "poisson"
  mod <- INLA::inla(formula, family = "poisson", data = data_run, offset = data_run$E_log,
                    control.predictor = list(compute=TRUE))
  dd=T/Ngrid/2
  ff<-function(x){
    return(exp(mod$summary.random$time$`0.5quant`)[1+sum(x>=dd+mod$summary.random$time$ID)])
  }
  
  ff2<-function(x){
    return(exp(mod$summary.random$time$`0.025quant`)[1+sum(x>=dd+mod$summary.random$time$ID)])
  }
  
  ff3<-function(x){
    return(exp(mod$summary.random$time$`0.975quant`)[1+sum(x>=dd+mod$summary.random$time$ID)])
  }
  bins=100
  #x1=seq(0,T1,T1/bins)[1:bins]
  x1=dataa
  intensity1=sapply(x1,inten2)*ccc
  med=sapply(x1,ff)
  low=sapply(x1,ff2)
  high=sapply(x1,ff3)
  width1=sum(high-low)/length(dataa)
  
  coverage=sum((intensity1>=low)*(intensity1<=high))/length(dataa)
  l2dist_med=sum((intensity1-med)**2)
  coverage2_1[jjj]=coverage
  l2_dist2_1[jjj]=l2dist_med
  width2_1[jjj]=width1
  
  
  bins=100
  x1=seq(0,T,T/bins)[1:bins]
  intensity1=sapply(x1,inten2)*ccc
  med=sapply(x1,ff)
  low=sapply(x1,ff2)
  high=sapply(x1,ff3)
  width2=sum(high-low)/bins
  
  coverage=sum((intensity1>=low)*(intensity1<=high))/bins
  l2dist_med=sum((intensity1-med)**2)
  coverage2_2[jjj]=coverage
  l2_dist2_2[jjj]=l2dist_med
  width2_2[jjj]=width2
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_2[jjj]=time.taken
}

##############################for synthetic example 3
coverage3_1=rep(0,300)
coverage3_2=rep(0,300)
l2_dist3_1=rep(0,300)
l2_dist3_2=rep(0,300)
width3_1=rep(0,300)
width3_2=rep(0,300)
time_3=rep(0,300)

for (jjj in 1:300){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',(jjj+600),'.rda',sep=''))
  ccc=ceil(jjj/100)
  T=100
  Ngrid=100
  hyper <- list(prec = list(param = c(.1, .1)))
  grid<-seq(0,T,T/Ngrid)
  lengthout <- length(grid) - 1
  field <- grid[-1] - diff(grid)/2
  m <- length(dataa)
  grid_trimmed <- setdiff(x = grid, y = dataa)
  sorting <- sort(c(grid_trimmed, dataa), index.return=TRUE)
  sgrid <- sorting$x
  ordering <- sorting$ix
  time_index <- cut(x = sgrid[-1], breaks = grid, labels = FALSE)
  time <- field[time_index]
  event_out <- c(rep(0, length(grid_trimmed)), rep(1,m))[ordering]
  E <- diff(sgrid)
  E_log = log(E)
  
  data_run <-data.frame(y =  event_out[-1], time = time, E_log = E_log, E= E)
  formula <- y ~  f(time, model="rw1", hyper = hyper, constr = FALSE, scale.model = TRUE)-1
  family <- "poisson"
  mod <- INLA::inla(formula, family = "poisson", data = data_run, offset = data_run$E_log,
                    control.predictor = list(compute=TRUE))
  
  dd=T/Ngrid/2
  ff<-function(x){
    return(exp(mod$summary.random$time$`0.5quant`)[1+sum(x>=dd+mod$summary.random$time$ID)])
  }
  
  ff2<-function(x){
    return(exp(mod$summary.random$time$`0.025quant`)[1+sum(x>=dd+mod$summary.random$time$ID)])
  }
  
  ff3<-function(x){
    return(exp(mod$summary.random$time$`0.975quant`)[1+sum(x>=dd+mod$summary.random$time$ID)])
  }
  bins=100
  x1=dataa
  intensity1=sapply(x1,inten3)*ccc
  med=sapply(x1,ff)
  low=sapply(x1,ff2)
  high=sapply(x1,ff3)
  width1=sum(high-low)/length(dataa)
  
  coverage=sum((intensity1>=low)*(intensity1<=high))/length(dataa)
  l2dist_med=sum((intensity1-med)**2)
  coverage3_1[jjj]=coverage
  l2_dist3_1[jjj]=l2dist_med
  width3_1[jjj]=width1
  
  
  bins=100
  x1=seq(0,T,T/bins)[1:bins]
  intensity1=sapply(x1,inten3)*ccc
  med=sapply(x1,ff)
  low=sapply(x1,ff2)
  high=sapply(x1,ff3)
  width2=sum(high-low)/bins
  
  coverage=sum((intensity1>=low)*(intensity1<=high))/bins
  l2dist_med=sum((intensity1-med)**2)
  coverage3_2[jjj]=coverage
  l2_dist3_2[jjj]=l2dist_med
  width3_2[jjj]=width2
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_3[jjj]=time.taken
}

save(coverage3_1,coverage3_2,l2_dist3_1,l2_dist3_2,time_3, width3_1,width3_2,
     coverage2_1,coverage2_2,l2_dist2_1,l2_dist2_2,time_2,width2_1,width2_2,
     coverage1_1,coverage1_2,l2_dist1_1,l2_dist1_2,time_1,width1_1,width1_2,
     file = "analysis.RData")



#######################################print outputs shown in tables##############################
###################synthetic example 1
round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage1_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width1_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist1_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage1_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width1_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)

round(quantile(l2_dist1_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage1_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width1_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist1_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage1_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width1_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_1[100:200]),2)
round(sd(time_1[100:200]),2)


round(quantile(l2_dist1_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage1_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width1_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist1_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage1_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width1_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_1[200:300]),2)
round(sd(time_1[200:300]),2)



###################synthetic example 2
round(quantile(l2_dist2_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width2_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist2_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage2_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width2_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_2[1:100]),2)
round(sd(time_2[1:100]),2)

round(quantile(l2_dist2_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width2_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist2_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage2_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width2_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_2[100:200]),2)
round(sd(time_2[100:200]),2)


round(quantile(l2_dist2_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width2_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist2_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage2_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width2_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_2[200:300]),2)
round(sd(time_2[200:300]),2)



###################synthetic example 3
round(quantile(l2_dist3_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage3_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width3_1[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist3_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage3_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width3_2[1:100], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_3[1:100]),2)
round(sd(time_3[1:100]),2)

round(quantile(l2_dist3_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage3_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width3_1[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist3_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage3_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width3_2[100:200], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_3[100:200]),2)
round(sd(time_3[100:200]),2)


round(quantile(l2_dist3_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage3_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width3_1[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(l2_dist3_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(quantile(100*coverage3_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),0)
round(quantile(width3_2[200:300], probs =c(0.5,0.025,0.25,0.75,0.975) , na.rm = FALSE),2)
round(mean(time_3[200:300]),2)
round(sd(time_3[200:300]),2)

