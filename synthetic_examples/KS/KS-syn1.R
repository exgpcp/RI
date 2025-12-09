library(spatstat.geom)     # owin, ppp, lookup.linim
library(spatstat.explore)  # bw.CvL (optional)
library(sparr)
library(parallel)


setwd("C:/Users/bingj/OneDrive - Stanford/syndata_rev")
inten1<-function(x){
  return(2*exp(-x/15)+exp(-((x-25)/10)**2))
}

###############################KS-PPL##################################
coverage1_1=rep(0,100)
width1_1=rep(0,100)
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
  lam_img <- suppressWarnings(
    density.ppp(Xpp, sigma = h_ppl, at = "pixels", edge = TRUE, se = TRUE)
  )

  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img$estimate, x = xxx, y = rep(0, Ngrid))
  SEfix <- lam_img$SE
  SEfix$v[is.na(SEfix$v) | SEfix$v < 0] <- 0
  se_xxx <- lookup.im(SEfix, x = xxx, y = rep(0, Ngrid))
  a=1.96
  low=-se_xxx*a+lambda_xxx
  high=se_xxx*a+lambda_xxx
  #compute SSE
  intensity1=sapply(xxx,inten1)
  l2dist=sum((intensity1-lambda_xxx)**2)
  width2=mean(2*a*se_xxx)
  coverage=sum((intensity1>=low)*(intensity1<=high))/length(xxx)
  l2_dist1_1[jjj]=l2dist
  coverage1_1[jjj]=coverage
  width1_1[jjj]=width2
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  time_1[jjj]=time.taken
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  lines(xxx, high, col = "red", lwd = 2, lty = 2)
  lines(xxx, low, col = "red", lwd = 2, lty = 2)
  lines(xxx, intensity1, col = "blue", lwd = 2)
}

round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)

###############################KS-SCOTT##################################
coverage1_1=rep(0,100)
width1_1=rep(0,100)
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
  lam_img <-density.ppp(Xpp, sigma =  h_scott[1], at = "pixels", edge = TRUE, se = TRUE)
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img$estimate, x = xxx, y = rep(0, Ngrid))
  SEfix <- lam_img$SE
  SEfix$v[is.na(SEfix$v) | SEfix$v < 0] <- 0
  se_xxx <- lookup.im(SEfix, x = xxx, y = rep(0, Ngrid))
  a=1.96
  low=-se_xxx*a+lambda_xxx
  high=se_xxx*a+lambda_xxx
  #compute SSE
  intensity1=sapply(xxx,inten1)
  l2dist=sum((intensity1-lambda_xxx)**2)
  width2=mean(2*a*se_xxx)
  coverage=sum((intensity1>=low)*(intensity1<=high))/length(xxx)
  l2_dist1_1[jjj]=l2dist
  coverage1_1[jjj]=coverage
  width1_1[jjj]=width2
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  time_1[jjj]=time.taken
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  lines(xxx, high, col = "red", lwd = 2, lty = 2)
  lines(xxx, low, col = "red", lwd = 2, lty = 2)
  lines(xxx, intensity1, col = "blue", lwd = 2)
}

round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)

###############################KS-adaptive##################################
coverage1_1=rep(0,100)
width1_1=rep(0,100)
l2_dist1_1=rep(0,100)
time_1=rep(0,100)
B = 500     # number of bootstrap samples
CORES = max(1, detectCores() - 1)
for (jjj in 1:100){
  start.time <- Sys.time()
  print(jjj)
  load(paste('data_',jjj,'.rda',sep=''))
  Tmax   <- 50
  Ngrid  <- 100
  xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]
  W   <- owin(xrange = c(0, Tmax), yrange = c(-0.5, 0.5))
  Xpp <- ppp(x = dataa, y = rep(0, length(dataa)), window = W)
  hp <- BOOT.density(Xpp)
  h0 <- OS(Xpp)
  lam_img <- bivariate.density(Xpp, h0=h0, hp=hp, adapt=TRUE, davies.baddeley=0.025)
  lambda_xxx <- lookup.im(lam_img$z, x=xxx, y=rep(0,Ngrid)) * length(dataa)
  intensity1 = sapply(xxx, inten1)
  l2_dist1_1[jjj] = sum((intensity1 - lambda_xxx)^2)
  
  ######################################
  ### PARALLEL BOOTSTRAP BEGINS HERE ###
  ######################################
  cl <- makeCluster(CORES)
  clusterExport(cl, c("dataa", "W", "xxx", "Ngrid",
                      "BOOT.density", "OS",
                      "bivariate.density", "ppp", "lookup.im"))
  boot_est_list <- parLapply(cl, 1:B, function(b){
    data_boot <- sample(dataa, replace=TRUE)
    Xboot <- ppp(x=data_boot, y=rep(0,length(data_boot)), window=W)
    hp_b <- BOOT.density(Xboot)
    h0_b <- OS(Xboot, nstar = "npoints")
    
    lam_boot <- bivariate.density(Xboot, h0=h0_b, hp=hp_b,
                                  adapt=TRUE, davies.baddeley=0.025)
    lookup.im(lam_boot$z, x=xxx, y=rep(0,Ngrid)) * length(data_boot)
  })
  stopCluster(cl)
  ## Convert bootstrap results into matrix
  boot_mat = do.call(cbind, boot_est_list)
  
  ## Bootstrap 95% CI
  low95  <- apply(boot_mat, 1, quantile, 0.025)
  high95 <- apply(boot_mat, 1, quantile, 0.975)
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  time_1[jjj] <- time.taken
  coverage1_1[jjj] = mean(intensity1 >= low95 & intensity1 <= high95)
  width1_1[jjj] = mean(high95 - low95)
  
  plot(xxx, lambda_xxx, type="l", ylim=c(0,15))
  lines(xxx, intensity1, col="blue", lwd=2)
  lines(xxx, low95, col="red", lty=2)
  lines(xxx, high95, col="red", lty=2)
  rug(dataa)
}

round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)

###############################KS-DIGGLE##################################
coverage1_1=rep(0,100)
width1_1=rep(0,100)
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
  h_dig <- bw.diggle(Xpp, hmax = 10)     
  lam_img <- suppressWarnings(density.ppp(Xpp, sigma = h_dig, at = "pixels", edge = TRUE, se = TRUE))
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img$estimate, x = xxx, y = rep(0, Ngrid))
  SEfix <- lam_img$SE
  SEfix$v[is.na(SEfix$v) | SEfix$v < 0] <- 0
  se_xxx <- lookup.im(SEfix, x = xxx, y = rep(0, Ngrid))
  a=1.96
  low=-se_xxx*a+lambda_xxx
  high=se_xxx*a+lambda_xxx
  #compute SSE
  intensity1=sapply(xxx,inten1)
  l2dist=sum((intensity1-lambda_xxx)**2)
  width2=mean(2*a*se_xxx)
  coverage=sum((intensity1>=low)*(intensity1<=high))/length(xxx)
  l2_dist1_1[jjj]=l2dist
  coverage1_1[jjj]=coverage
  width1_1[jjj]=width2
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  time_1[jjj]=time.taken
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  lines(xxx, high, col = "red", lwd = 2, lty = 2)
  lines(xxx, low, col = "red", lwd = 2, lty = 2)
  lines(xxx, intensity1, col = "blue", lwd = 2)
}

round(quantile(l2_dist1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width1_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_1[1:100]),2)
round(sd(time_1[1:100]),2)

