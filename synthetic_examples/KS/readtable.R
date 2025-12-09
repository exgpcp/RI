coverage2_1 = rep(0,100)
width2_1    = rep(0,100)
l2_dist2_1  = rep(0,100)
time_2      = rep(0,100)

for (i in 1:100) {
  
  infile <- paste0("/home/groups/juliapr/Bingjing/code/currentcode/KS/syn2data_",
                   i+300, ".rda")
  
  load(infile)   # loads: l2_dist, width, coverage, time
  
  l2_dist2_1[i] = l2_dist
  width2_1[i]   = width
  coverage2_1[i]= coverage
  time_2[i]     = time
}


round(quantile(l2_dist2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(quantile(100*coverage2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),0)
round(quantile(width2_1[1:100], probs =c(0.5,0.25,0.75) , na.rm = FALSE),2)
round(mean(time_2[1:100]),2)
round(sd(time_2[1:100]),2)
