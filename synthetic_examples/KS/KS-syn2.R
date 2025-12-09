library(spatstat.geom)     # owin, ppp, lookup.linim
library(spatstat.explore)  # bw.CvL (optional)
library(sparr)
library(parallel)
setwd("C:/Users/bingj/OneDrive - Stanford/syndata_rev")

inten2<-function(x){
  return(10+x-x)
}


##############################KS-PPL################################################
coverage2_1=rep(0,100)
width2_1=rep(0,100)
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
  #lam_img <- density.ppp(Xpp, sigma = h_ppl, at = "pixels", edge = TRUE)
  lam_img <- suppressWarnings(density.ppp(Xpp, sigma = h_ppl, at = "pixels", edge = TRUE,se=TRUE))
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img$estimate, x = xxx, y = rep(0, Ngrid))
  SEfix <- lam_img$SE
  SEfix$v[is.na(SEfix$v) | SEfix$v < 0] <- 0
  se_xxx <- lookup.im(SEfix, x = xxx, y = rep(0, Ngrid))
  a=1.96
  low=-se_xxx*a+lambda_xxx
  high=se_xxx*a+lambda_xxx
  
  #compute SSE
  intensity2=sapply(xxx,inten2)
  l2dist=sum((intensity2-lambda_xxx)**2)
  width2=mean(2*a*se_xxx)
  coverage=sum((intensity2>=low)*(intensity2<=high))/length(xxx)
  l2_dist2_1[jjj-300]=l2dist
  coverage2_1[jjj-300]=coverage
  width2_1[jjj-300]=width2
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  time_2[jjj-300]=time.taken

  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  lines(xxx, high, col = "red", lwd = 2, lty = 2)
  lines(xxx, low, col = "red", lwd = 2, lty = 2)
  lines(xxx, intensity2, col = "blue", lwd = 2)
}

round(quantile(l2_dist2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2[1:100]),2)
round(sd(time_2[1:100]),2)

#######################################KS-SCOTT########################################
coverage2_1=rep(0,100)
width2_1=rep(0,100)
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
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lam_img <- suppressWarnings(density.ppp(Xpp, sigma =h_scott[1], at = "pixels", edge = TRUE,se=TRUE))
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img$estimate, x = xxx, y = rep(0, Ngrid))
  SEfix <- lam_img$SE
  SEfix$v[is.na(SEfix$v) | SEfix$v < 0] <- 0
  se_xxx <- lookup.im(SEfix, x = xxx, y = rep(0, Ngrid))
  a=1.96
  low=-se_xxx*a+lambda_xxx
  high=se_xxx*a+lambda_xxx
  #compute SSE
  intensity2=sapply(xxx,inten2)
  l2dist=sum((intensity2-lambda_xxx)**2)
  width2=mean(2*a*se_xxx)
  coverage=sum((intensity2>=low)*(intensity2<=high))/length(xxx)
  l2_dist2_1[jjj-300]=l2dist
  coverage2_1[jjj-300]=coverage
  width2_1[jjj-300]=width2
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  time_2[jjj-300]=time.taken
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  lines(xxx, high, col = "red", lwd = 2, lty = 2)
  lines(xxx, low, col = "red", lwd = 2, lty = 2)
  lines(xxx, intensity2, col = "blue", lwd = 2)
}

round(quantile(l2_dist2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2[1:100]),2)
round(sd(time_2[1:100]),2)

##############################KS-DIGGLE################################################
coverage2_1=rep(0,100)
width2_1=rep(0,100)
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
  h_dig <- bw.diggle(Xpp,hmax=2)    
  #lam_img <- density.ppp(Xpp, sigma = h_ppl, at = "pixels", edge = TRUE)
  lam_img <- suppressWarnings(density.ppp(Xpp, sigma = h_dig, at = "pixels", edge = TRUE,se=TRUE))
  ## --- evaluate λ(t) at grid points (y = 0) ---
  lambda_xxx <- lookup.im(lam_img$estimate, x = xxx, y = rep(0, Ngrid))
  SEfix <- lam_img$SE
  SEfix$v[is.na(SEfix$v) | SEfix$v < 0] <- 0
  se_xxx <- lookup.im(SEfix, x = xxx, y = rep(0, Ngrid))
  a=1.96
  low=-se_xxx*a+lambda_xxx
  high=se_xxx*a+lambda_xxx
  #compute SSE
  intensity2=sapply(xxx,inten2)
  l2dist=sum((intensity2-lambda_xxx)**2)
  width2=mean(2*a*se_xxx)
  coverage=sum((intensity2>=low)*(intensity2<=high))/length(xxx)
  l2_dist2_1[jjj-300]=l2dist
  coverage2_1[jjj-300]=coverage
  width2_1[jjj-300]=width2
  end.time <- Sys.time()
  time.taken <- as.numeric(end.time - start.time, units = "secs")
  time_2[jjj-300]=time.taken
  # quick plot
  plot(xxx, lambda_xxx , type = "l", xlab = "t", ylab = expression(hat(lambda)(t)),
       main = "1D intensity via density.ppp (fixed bandwidth)")
  rug(dataa)
  lines(xxx, high, col = "red", lwd = 2, lty = 2)
  lines(xxx, low, col = "red", lwd = 2, lty = 2)
  lines(xxx, intensity2, col = "blue", lwd = 2)
}


round(quantile(l2_dist2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2[1:100]),2)
round(sd(time_2[1:100]),2)


##############################KS-ADAPTIVE################################################
coverage2_1=rep(0,100)
width2_1=rep(0,100)
l2_dist2_1=rep(0,100)
time_2=rep(0,100)
B = 200     # number of bootstrap samples
CORES = max(1, detectCores() - 1)

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
  ## Bandwidths
  hp <- tryCatch(bw.SJ(dataa), error = function(e) bw.nrd0(dataa))
  h0 <- 2.5 * hp
  ## Original estimator
  lam_img <- bivariate.density(Xpp, h0=h0, hp=hp, adapt=TRUE, davies.baddeley=0.025)
  lambda_xxx <- lookup.im(lam_img$z, x=xxx, y=rep(0,Ngrid)) * length(dataa)
  ## True λ
  intensity2 = sapply(xxx, inten2)
  ## L2 error
  l2_dist2_1[jjj] = sum((intensity2 - lambda_xxx)^2)
  
  ######################################
  ### PARALLEL BOOTSTRAP BEGINS HERE ###
  ######################################
  cl <- makeCluster(CORES)
  clusterExport(cl, c("dataa", "W", "xxx", "Ngrid", 
                      "bw.SJ", "bw.nrd0", "bivariate.density", "ppp", "lookup.im"))
  boot_est_list <- parLapply(cl, 1:B, function(b){
    data_boot <- sample(dataa, replace=TRUE)
    hp_b <- tryCatch(bw.SJ(data_boot), error=function(e) bw.nrd0(data_boot))
    h0_b <- 2.5 * hp_b
    Xboot <- ppp(x=data_boot, y=rep(0,length(data_boot)), window=W)
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
  
  end.time<-Sys.time()
  time_2[jjj-300] <- as.numeric(end.time - start.time, units = "secs")
  coverage2_1[jjj-300] = mean(intensity2 >= low95 & intensity2 <= high95)
  width2_1[jjj-300] = mean(high95 - low95)
  plot(xxx, lambda_xxx, type="l",ylim=c(5,20))
  lines(xxx, intensity2, col="blue", lwd=2)
  lines(xxx, low95, col="red", lty=2)
  lines(xxx, high95, col="red", lty=2)
  rug(dataa)
}

round(quantile(l2_dist2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2[1:100]),2)
round(sd(time_2[1:100]),2)

