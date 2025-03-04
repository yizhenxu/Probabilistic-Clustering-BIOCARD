# data simulated from data_processing_BMsim11Y.R
# stan model specified in BMsim11Y_stan.R
# rho_theta ~ Unif(0, 1.5)
## run yizhenxu/simparallel.sh for parallelizing k posterior draws: BMsim11Y.R 
# cd ./AD
# qsub bash1.sh
rm(list=ls())
taskid = commandArgs(trailingOnly=TRUE); print(taskid)

#rm(list=ls()) 
# coerce the value to an integer
#taskid = as.numeric(Sys.getenv('SGE_TASK_ID'))# 1:3000

source("BMsim11Y_setting.R")
catlen = c(2,4,3)

sfn = paste0("/fastscratch/myscratch/yxu2/simu11_",taskid,".RData")
oldfn = paste0("/users/yxu2/AD/simu11_",taskid,".RData")

if(!(file.exists(sfn) || file.exists(oldfn))){
  
  book = expand.grid(1:3, 1:1000) 
  book$id = 1:nrow(book) 
  modid = book[taskid,1]
  dataid = book[taskid,2]
  
  if(modid == 2){ # only run 2-cluster model
    g1 = gamma1 # same as that in BMsim_setting.R
    
    library(data.table)
    library(rstan)
    
    #setwd("./AD") # working dir already set at qsub level
    
    load("simdata11.RData")
    
    #fn = paste0("mod",modid,".stan")
    if(modid == 1) fn = "mod1.stan"
    if(modid > 1) fn = "mod2.stan"
    
    d = simd[dataID == dataid,]
    
    x = cbind(1, as.matrix(d[, c("X1", "X2")]))
    z = as.matrix( d[, "age"])
    y = as.matrix(d[, paste0("Y",1:K), with=F])
    
    nit = 10000; nwup = 1000; nc = 1; seednum = 8
    
    if(modid == 1){
      staninput = list(Nsubs = length(unique(d[, fID])),
                       Nobs = nrow(d),
                       K = ncol(y),
                       Npreds = ncol(x),
                       NZ = ncol(z),
                       X = x,
                       Z = z,
                       time = as.numeric(d$age),
                       lsub = d[,.N, by = fID]$N,
                       y = y,
                       gamma1 = g1,
                       catlen = catlen,
                       catsum = cumsum(c(0,catlen))[1:3])
      
      parvec = c("gamma2","gamma3",
                 "betat","betaX","Sig",
                 "alpha",
                 "lp__", "rho_theta", "d") 
      
    }
    
    
    if(modid %in% c(2,3)){
      
      staninput = list(Nsubs = length(unique(d[,fID])),
                       Nobs = nrow(d),
                       K = ncol(y),
                       Npreds = ncol(x),
                       NZ = ncol(z),
                       X = x,
                       Z = z,
                       time = as.numeric(d$age),
                       lsub = d[,.N, by=fID]$N,
                       endidx = cumsum( d[,.N, by=fID]$N ),
                       y = y,
                       gamma1 = g1,
                       L = modid,
                       catlen = catlen,
                       catsum = cumsum(c(0,catlen))[1:3])
      
      parvec = c("gamma2","gamma3",
                 "betat","betaX","Sig",
                 "alpha0","alpha",
                 "lp__","rho_theta","d", "lambda")
    }
    
    
    options(mc.cores = nc) #parallel::detectCores()
    fit = stan(file = fn, data = staninput,
               pars = parvec , include = T,
               chains = nc, iter = nit, warmup = nwup, thin = 1, seed = seednum)
    
    res = extract(fit, permuted =F,inc_warmup=F)
    save(res, file = sfn)  # simu updated to resctrict the range of rho_theta
  } # modid
} # check file.exists

