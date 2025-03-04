#qrsh -l mem_free=10G,h_vmem=10G,h_fsize=20G
#module load conda_R/4.3.x

#do this:
if(0){
  fidx = 1 # test data being the fidx_th fold
  source("./AD/BM2_opt_Feb2023_AccuracyCV.R")
}

setwd("./AD")

#route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/"
#setwd(route)

a = 3;b=1

type = 'B'; ageopt = 2 # 

source("data_processing_Feb2023_noWMH.R")

# trim data to be the training set for the fidx_th fold

eachct = d[, .N, by=Study.ID]
idlist = eachct[N > 1, Study.ID]
Nid = length(idlist)
nfold = 5
fsize = c(rep(round(Nid/nfold), nfold - 1), Nid - round(Nid/nfold)*(nfold-1))

if(fidx == 1){
  testidx = 1:fsize[1]
} else {
  testidx = (sum(fsize[1:(fidx-1)])+1):sum(fsize[1:fidx])
}
trainidx = (1:Nid)[-testidx]

set.seed(1)
idlist_forsampling = sample(idlist, Nid)
trainid = idlist_forsampling[trainidx]
testid = idlist_forsampling[testidx]

d = d[Study.ID %in% trainid,]

# create data list for stan 

agedrt = d[,age[.N]-age[1], by = get(id)]$V1
#min and max go into l and u of prior tuning
agemin = as.numeric(d[, age - c(NA, age[-.N]), by=get(id)]$V1) #here

library(invgamma)
lu = c(min(agemin,na.rm=T), max(agedrt)) #here
foo = function(par){
  a = pinvgamma(lu[1], shape = par[1], scale = par[2])
  b = 1 - pinvgamma(lu[2], shape = par[1], scale = par[2])
  return((a-0.01)^2+(b-0.01)^2)
}
rhoab = optim(c((sum(lu)/(lu[2]-lu[1]))^2+2, sum(lu)/2*(sum(lu)/(lu[2]-lu[1]))^2+1),foo)$par

catlen = c(length(cog), length(mri), length(csf))

xname = c("apoe", "SEX", "education")
x <- cbind(1, as.matrix(d[, xname,with=F]))

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

y <- d[, all,with=F]
yfill0 = y
yfill0[is.na(yfill0)] = 0
subjvec = rep(1:length(unique(d[,get(id)])),d[,.N,by=get(id)]$N)


staninput = list(Nsubs = length(unique(d[,get(id)])),
                 Nobs = nrow(d),
                 K = ncol(y),
                 Npreds = ncol(x),
                 NZ = ncol(z),
                 X = x,
                 Z = z,
                 time = as.numeric(d$age),
                 lsub = d[,.N, by=get(id)]$N,
                 endidx = cumsum( d[,.N, by=get(id)]$N ),
                 subjvec = subjvec,
                 y_observed = 1*(!is.na(y)),
                 y = yfill0,
                 gamma1 = gamma1,
                 rho_lower = rhoab[1],
                 rho_upper = rhoab[2],
                 L = 2,
                 y1on0 = as.numeric(1*(!is.na(y[,1]) & y[,1] == min(y[,1],na.rm=T)) ), 
                 a = a,b=b,
                 catlen = catlen,
                 catsum = cumsum(c(0,catlen))[1:3])


nit = 1000; nwup = 600; nc = 10; seednum = 8

fn <- paste0("mod",fidx,".stan")
if (file.exists(fn))  file.remove(fn)
source("scode.R")
writeLines(scode1, fn)
rstan:::rstudio_stanc(fn)

options(mc.cores = nc) #parallel::detectCores()
fit = stan(file = fn, data = staninput,
           pars = c("gamma2","gamma3",
                    "betat","betaX","Sig",
                    "alpha0","alpha","rho_theta1","rho_theta2",
                    "lp__","log_lik", "logpyz", "lambda","d","gamma3mu") , include = T,
           chains = nc, iter = nit, warmup = nwup, thin = 1, seed = seednum)

save(fit, file = paste0("BM2_opt2_Feb2023_noWMH_Accuracy_fold",fidx,".RData"))  #here
#
