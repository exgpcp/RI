library(spatstat.geom)     # owin, ppp, lookup.linim
library(spatstat.explore)  # bw.CvL (optional)
setwd("C:/Users/bingj/OneDrive - Stanford/syndata_rev")

inten2<-function(x){
  return(10+x-x)
}



##############################KS-PPL################################################
l2_dist2_1=rep(0,100)
time_2=rep(0,100)

for (jjj in 301:400){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  Tmax   <- 5
  Ngrid  <- 100
  xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]  
  
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
  intensity2=sapply(xxx,inten2)
  l2dist=sum((intensity2-lambda_xxx)**2)
  l2_dist2_1[jjj-300]=l2dist
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_2[jjj-300]=time.taken
}

round(quantile(l2_dist2_1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2),2)
round(sd(time_2),2)



#######################################KS-SCOTT########################################
l2_dist2_1=rep(0,100)
time_2=rep(0,100)
for (jjj in 301:400){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  Tmax   <- 5
  Ngrid  <- 100
  xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]  
  ## --- embed 1D times on a thin horizontal strip and build a planar ppp ---
  # Use height = 1 so per-area intensity matches per-x intensity along y=0
  W   <- owin(xrange = c(0, Tmax), yrange = c(-0.5, 0.5))
  Xpp <- ppp(x = dataa, y = rep(0, length(dataa)), window = W)
  
  ## --- choose a fixed bandwidth (reviewer-suggested selectors) ---
  h_scott   <- bw.scott(Xpp)       # Poisson likelihood CV
  lam_img <- density.ppp(Xpp, sigma = h_scott[1], at = "pixels", edge = TRUE)
  
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img, x = xxx, y = rep(0, Ngrid))
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  
  intensity2=sapply(xxx,inten2)
  l2dist=sum((intensity2-lambda_xxx)**2)
  l2_dist2_1[jjj-300]=l2dist
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_2[jjj-300]=time.taken
}

round(quantile(l2_dist2_1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2),2)
round(sd(time_2),2)


#####################################KS-FIXED value Tmax/2
l2_dist2_1=rep(0,100)
time_2=rep(0,100)
for (jjj in 301:400){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  Tmax   <- 5
  Ngrid  <- 100
  xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]  
  ## --- embed 1D times on a thin horizontal strip and build a planar ppp ---
  # Use height = 1 so per-area intensity matches per-x intensity along y=0
  W   <- owin(xrange = c(0, Tmax), yrange = c(-0.5, 0.5))
  Xpp <- ppp(x = dataa, y = rep(0, length(dataa)), window = W)
  ## --- choose a fixed bandwidth (reviewer-suggested selectors) ---
  lam_img <- density.ppp(Xpp, sigma =Tmax/2, at = "pixels", edge = TRUE)
  
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img, x = xxx, y = rep(0, Ngrid))
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  intensity2=sapply(xxx,inten2)
  l2dist=sum((intensity2-lambda_xxx)**2)
  l2_dist2_1[jjj-300]=l2dist
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time_2[jjj-300]=time.taken
}

round(quantile(l2_dist2_1, probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2),2)
round(sd(time_2),2)



