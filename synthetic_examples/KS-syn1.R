
#library(parallel)
#library(pracma)
#library(TruncatedNormal)
#library(matrixcalc) # is.positive.definite
#library(matlib)
#library(Matrix)
#library(base)
#library(mvtnorm)
#library(truncnorm)
library(spatstat.geom)     # owin, ppp, lookup.linim
library(spatstat.explore)  # bw.CvL (optional)


setwd("C:/Users/bingj/OneDrive - Stanford/syndata_rev")
inten1<-function(x){
  return(2*exp(-x/15)+exp(-((x-25)/10)**2))
}

###############################KS-PPL##################################

l2_dist1_1=rep(0,100)
time_1=rep(0,100)

for (jjj in 1:100){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  
  Tmax   <- 50
  Ngrid  <- 100
  xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]   # (0, 50], length 100
  
  ## --- embed 1D times on a thin horizontal strip and build a planar ppp ---
  # Use height = 1 so per-area intensity matches per-x intensity along y=0
  W   <- owin(xrange = c(0, Tmax), yrange = c(-0.5, 0.5))
  Xpp <- ppp(x = dataa, y = rep(0, length(dataa)), window = W)
  
  h_ppl   <- bw.ppl(Xpp)       # Poisson likelihood CV
  lam_img <- density.ppp(Xpp, sigma = h_ppl, at = "pixels", edge = TRUE)

  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img, x = xxx, y = rep(0, Ngrid))
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  
  #compute SSE
  intensity1=sapply(xxx,inten1)
  l2dist=sum((intensity1-lambda_xxx)**2)
  l2_dist1_1[jjj]=l2dist
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_1[jjj]=time.taken
}


round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)



###############################KS-SCOTT##################################
l2_dist1_1=rep(0,100)
time_1=rep(0,100)

for (jjj in 1:100){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  Tmax   <- 50
  Ngrid  <- 100
  xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]   # (0, 50], length 100
  
  ## --- embed 1D times on a thin horizontal strip and build a planar ppp ---
  # Use height = 1 so per-area intensity matches per-x intensity along y=0
  W   <- owin(xrange = c(0, Tmax), yrange = c(-0.5, 0.5))
  Xpp <- ppp(x = dataa, y = rep(0, length(dataa)), window = W)
  h_scott   <- bw.scott(Xpp)      
  lam_img <- density.ppp(Xpp, sigma = h_scott[1], at = "pixels", edge = TRUE)
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img, x = xxx, y = rep(0, Ngrid))
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  intensity1=sapply(xxx,inten1)
  l2dist=sum((intensity1-lambda_xxx)**2)
  l2_dist1_1[jjj]=l2dist
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_1[jjj]=time.taken
}

round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)


######################KS-fixed value Tmax/2#####################
l2_dist1_1=rep(0,100)
time_1=rep(0,100)

for (jjj in 1:100){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  Tmax   <- 50
  Ngrid  <- 100
  xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]   # (0, 50], length 100
  
  ## --- embed 1D times on a thin horizontal strip and build a planar ppp ---
  # Use height = 1 so per-area intensity matches per-x intensity along y=0
  W   <- owin(xrange = c(0, Tmax), yrange = c(-0.5, 0.5))
  Xpp <- ppp(x = dataa, y = rep(0, length(dataa)), window = W)
  
  ## --- choose a fixed bandwidth (Tmax/2) ---
  lam_img <- density.ppp(Xpp, sigma = Tmax/2, at = "pixels", edge = TRUE)
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img, x = xxx, y = rep(0, Ngrid))
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  intensity1=sapply(xxx,inten1)
  l2dist=sum((intensity1-lambda_xxx)**2)
  l2_dist1_1[jjj]=l2dist
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_1[jjj]=time.taken
}

round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)

