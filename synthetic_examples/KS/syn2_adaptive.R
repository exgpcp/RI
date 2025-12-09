#################################KS
###################for the first 100 datasets in the first synthetic example
#!/usr/bin/Rscript
args = commandArgs(trailingOnly=TRUE)
library(spatstat)
library(parallel)
library(spatstat.geom)     # owin, ppp, lookup.linim
library(spatstat.explore)
library(sparr)


nam <- paste0("/home/groups/juliapr/Bingjing/code/currentcode/syndata/data_",args[1],".rda")#args[1] range from 1 to 100
load(nam)
inten2<-function(x){
  return(10+x-x)
}


B = 200     # number of bootstrap samples
CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
start.time <- Sys.time()
Tmax   <- 5
Ngrid  <- 100
xxx    <- seq(0, Tmax, length.out = Ngrid + 1)[-1]
W   <- owin(xrange = c(0, Tmax), yrange = c(-0.5, 0.5))
Xpp <- ppp(x = dataa, y = rep(0, length(dataa)), window = W)
hp <- BOOT.density(Xpp)
h0 <- OS(Xpp)
lam_img <- bivariate.density(Xpp, h0=h0, hp=hp, adapt=TRUE, davies.baddeley=0.025)
lambda_xxx <- lookup.im(lam_img$z, x=xxx, y=rep(0,Ngrid)) * length(dataa)
intensity2 = sapply(xxx, inten2)
l2_dist= sum((intensity2 - lambda_xxx)^2)

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
  h0_b <- OS(Xboot)
  
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
time<- time.taken
coverage = mean(intensity2 >= low95 & intensity2 <= high95)
width = mean(high95 - low95)

outfile <- paste0("/home/groups/juliapr/Bingjing/code/currentcode/KS/syn2data_",
                  args[1],
                  ".rda")

save(l2_dist, width, coverage, time,
     file = outfile)
