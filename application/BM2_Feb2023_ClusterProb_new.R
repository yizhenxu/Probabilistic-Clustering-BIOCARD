rm(list = ls())

local = 0

library(MASS) # mvrnom

type = 'B'; agesopt = 2 # 

if(local == 1){
  route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/"
  source(paste0(route, "PlotFunctions.R"))
  fn = paste0("/Users/yizhenxu/BM2_opt2_Feb2023_noWMH.RData")
  #fn = "BM2_opt2_Feb2023_noWMH.RData"
} else {
  setwd("./AD")
  #route = "/users/yxu2/AD/"
  route = ""
  source("PlotFunctions.R")
  fn = "BM2_opt2_Feb2023_noWMH.RData" # old syntax with all data
  fn = "BM2_opt2_Feb2023_noWMH_Accuracy_newSyntax.RData" # new syntax excluding last observation
  
}


library(rstan)

type = "B"; ageopt = 2
source(paste0(route, "data_processing_Feb2023_noWMH.R"))

if(file.exists(fn)) load(fn)

od = 1; L = 2;  nY = length(all)
res = extract(fit, permuted =F,inc_warmup=F)
dn = dimnames(res)
lp = res[,,which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)


zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

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

datlist = list(Nsubs = length(unique(d[,get(id)])),
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
               y1on0 = 1*(!is.na(y[,1]) & y[,1] == min(y[,1],na.rm=T)), 
               catlen = catlen,
               catsum = cumsum(c(0,catlen))[1:3])

npost = 400

betaX = array(res[, odch[od], grep("betaX",dn$parameters )], dim = c(npost, datlist$K, 4))
betat = res[, odch[od], grep("betat",dn$parameters )]
gamma2 = res[, odch[od], grep("gamma2",dn$parameters )]
gamma3 = res[, odch[od], grep("gamma3\\[",dn$parameters )]
alpha0 = res[, odch[od], grep("alpha0",dn$parameters )]
alpha = res[, odch[od], grep("alpha\\[",dn$parameters )]
rho = res[, odch[od], grep("rho",dn$parameters )]
theta = array(res[, odch[od], grep("d\\[",dn$parameters )], dim = c(npost, datlist$Nobs, datlist$L))
sigma = array(res[, odch[od], grep("Sig",dn$parameters )], dim = c(npost, datlist$L, datlist$K) )
lambda = res[, odch[od], grep("lambda\\[",dn$parameters )]

posttheta = res[, odch[od], grep("d\\[",dn$parameters )]

#########################################################################
#########################################################################
# to compile, (shift command .) to show hidden file, 
# download gfortran from https://mac.r-project.org/tools/
# locate under /usr/local/gfortran/lib the dylib's reported missing in the errors
# create new folder /usr/local/lib/ and move those dylib's into this folder
if(local){
  library("Rcpp")
  Rcpp::sourceCpp(paste0(route,"cfun_AD.cpp")) 
} else {
  library("Rcpp")
  Rcpp::sourceCpp("cfun_AD.cpp")
}

id = 1; s = 1; cluster = 1
idx = which( datlist$subjvec == id ); Ni = length(idx)
cc  = 3
hidx = idx[1:cc]
x = rnorm(cc)
tmp = ntheta_posterior_c(x, cluster, 
                   matrix(datlist$X[hidx,],nrow=cc), 
                   matrix(datlist$Z[hidx,],nrow=cc), datlist$time[hidx], 
                   matrix(datlist$y[hidx,],nrow=cc), matrix(datlist$y_observed[hidx,],nrow=cc),datlist$y1on0[hidx], 
                   betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                   gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )


tmp = ntheta_gradient_c(x,  cluster, 
                        matrix(datlist$X[hidx,],nrow=cc), 
                        matrix(datlist$Z[hidx,],nrow=cc), datlist$time[hidx], 
                        matrix(datlist$y[hidx,],nrow=cc),matrix(datlist$y_observed[hidx,],nrow=cc), datlist$y1on0[hidx], 
                        betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                        gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )

#compare optim(R) and optsol(rcpp)
ptm_mclapply <- proc.time()  
tmp = optsol( cluster, matrix(datlist$X[hidx,],nrow=cc), 
                  matrix(datlist$Z[hidx,],nrow=cc), datlist$time[hidx], 
                  matrix(datlist$y[hidx,],nrow=cc),matrix(datlist$y_observed[hidx,],nrow=cc), datlist$y1on0[hidx], 
                  betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                  gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
proc.time() - ptm_mclapply # 

tmp = calc_prob_laplace_persondraw(idx, s = 1,
                           datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                           betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho,lambda)

getall_probClst(idx, length(idx), 100,
                datlist$X, datlist$Z, datlist$time, datlist$y,  datlist$y_observed, datlist$y1on0,
                betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho, lambda)

check = vector("list", length(unique(datlist$subjvec)))
see = vector("list", length(unique(datlist$subjvec)))

s = 1 # only look at first posterior draw
 ct = 1
for( id in unique(datlist$subjvec)){
  idx = which( datlist$subjvec == id )
  check[[ct]] = calc_prob_laplace_persondraw(idx-1, s = 1,
                             datlist$X, datlist$Z, datlist$time, datlist$y,  datlist$y_observed, datlist$y1on0,
                             betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho,lambda)
  see[[ct]] = check[[ct]][2]
  ct = ct + 1
}

 # define allidx and eachct
 idlist1 = unique(d$Study.ID)
 allidx1 = which(d[,Study.ID] %in% idlist1)
 
 eachct1 =  d[Study.ID %in% idlist1,.N, by=Study.ID]$N
 
 ptm_mclapply <- proc.time()  
 getPs = getall_probClst(allidx1, eachct1,  400,
                 datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                 betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho, lambda)
 proc.time() - ptm_mclapply # 
 
 IDshow = rep(idlist1, eachct1)
 show = cbind(IDshow, apply(getPs, 1, function(x) mean(x,na.rm=T)),round(getPs[,1:400],2))
 
 save(getPs, IDshow, show, file = "BM2_Feb2023_Projection_P.RData")
 
 
 
 
 
 
 
 
 
 if(0){ # tested in R
   
   # theta is vector of length Ni
   theta_posterior = function(theta, Ni, cluster, X, Z, time, y, K,
                              betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho){
     
     di = - ifelse(cluster == 1, 0, alpha0)  + Z %*% alpha + theta
     #mu[person time, outcome] 
     mui = X %*% t(betaX) + 
       matrix(time ,ncol=1)%*% matrix(betat, nrow=1) + 
       mapply(function(j) gamma1[j]/(1+exp(-gamma2[j]*(di-gamma3[j]))), 1:K)
     #gamma1/(1 + exp(-gamma2[s,]*(di-gamma3[s,]) ) );
     
     yi = matrix(y,nrow = Ni)
     obs_ptj = matrix(mapply(function(j) -(1/(2*sigma[cluster,j]))*(yi[,j] - mui[,j])^2, 1:K), nrow = Ni)
     
     
     Dmat = exp(-(outer(time, time, "-")^2)/(2*rho[cluster]^2)) + diag(1e-10,Ni)
     logpj = sum(apply(obs_ptj, 1, function(x) sum(x,na.rm=T))) - t(theta) %*% solve(Dmat) %*% theta/2 # check
     return(logpj)
   }
   
   theta_gradient = function(theta, Ni, cluster, X, Z, time, y, K,
                             betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho){
     
     yi = matrix(y,nrow = Ni)
     
     di = - ifelse(cluster == 1, 0, alpha0)  + Z %*% alpha + theta
     tmp = matrix(mapply(function(j) exp(-gamma2[j]*(di-gamma3[j])), 1:K),nrow=Ni)
     #mu[person time, outcome] 
     mui = X %*% t(betaX) + 
       matrix(time,ncol=1)%*% matrix(betat, nrow=1) + 
       mapply(function(j) gamma1[j]/(1+tmp[,j]), 1:K)
     
     r = matrix(mapply(function(j) (1/sigma[cluster,j])*(yi[,j] - mui[,j])*
                         gamma1[j]*tmp[,j]*gamma2[j]/(1+tmp[,j])^2, 1:K), nrow=Ni)
     
     Dmat = exp(-(outer(time, time, "-")^2)/(2*rho[cluster]^2)) + diag(1e-10,Ni)
     drvt = apply(r, 1,  function(x) sum(x,na.rm=T)) - solve(Dmat) %*% theta
     return(drvt)
   }
   
   theta_V = function(time, rho, cluster, Ni){
     
     Dmat = exp(-(outer(time, time, "-")^2)/(2*rho[cluster]^2)) + diag(1e-10,Ni)
     return( Dmat )
   }
   
   #########################################################################
   # for each person id and each posterior index s
   
   get_theta_all = function(id, s, datlist, 
                            betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho){
     idx = which( datlist$subjvec == id ); Ni = length(idx)
     theta_all = rnorm(datlist$L)
     for(hl in 1:Ni){
       hidx = idx[1:hl]
       thetapar = thetaV =  vector(mode='list', length = datlist$L)
       for(cluster in 1:2){
         optsol = optim(rep(0,hl),  theta_posterior, theta_gradient, "method" = "BFGS", control = list(fnscale = -1),
                        Ni = hl, cluster = cluster, 
                        X = datlist$X[hidx,], Z = datlist$Z[hidx,], time = datlist$time[hidx], y = datlist$y[hidx,], K = datlist$K,
                        betaX = betaX[s,,], betat = betat[s,], sigma = sigma[s,,], 
                        gamma1 = gamma1, gamma2 = gamma2[s,], gamma3 = gamma3[s,], 
                        alpha0 = alpha0[s], alpha = alpha[s,], rho = rho[s,])
         thetapar[[cluster]] = optsol$par
         thetaV[[cluster]] = theta_V( datlist$time[hidx], rho[s,], cluster, length(hidx))
       }
       
       thetasamp = c()
       for(cluster in 1:2){
         thetasamp = cbind(thetasamp, mvrnorm(n = 1, thetapar[[cluster]], thetaV[[cluster]]) )
       }
       theta_all = rbind(theta_all, thetasamp[hl,])
     } # hl
     return(theta_all)
   }
   
   
   # for each person id and each posterior index s
   # simulate m theta samples for the cluster probabilities under posterior index s
   calc_prob_persondraw_R = function(m, id, s, datlist, 
                                     betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho,lambda){
     idx = which( datlist$subjvec == id ); Ni = length(idx)
     prob_all = c()
     for(hl in 1:Ni){
       hidx = idx[1:hl]
       thetapar = thetaV =  vector(mode='list', length = datlist$L)
       prob = c()
       for(cluster in 1:2){
         optsoll = optim(rep(0,hl),  theta_posterior, theta_gradient, "method" = "BFGS", control = list(fnscale = -1),
                         Ni = length(hidx), cluster = cluster, 
                         X = datlist$X[hidx,], Z = datlist$Z[hidx,], time = datlist$time[hidx], y = datlist$y[hidx,], K = datlist$K,
                         betaX = betaX[s,,], betat = betat[s,], sigma = sigma[s,,], 
                         gamma1 = gamma1, gamma2 = gamma2[s,], gamma3 = gamma3[s,], 
                         alpha0 = alpha0[s], alpha = alpha[s,], rho = rho[s,])
         thetapar[[cluster]] = optsoll$par
         thetaV[[cluster]] = theta_V( datlist$time[hidx], rho[s,], cluster, length(hidx))
         thetasamp = mvrnorm(n = m, thetapar[[cluster]], thetaV[[cluster]])
         denom = apply(thetasamp, 1, function(x) dmvnorm(x, thetapar[[cluster]], thetaV[[cluster]], log=F))
         numer = apply(thetasamp, 1, function(x) exp(theta_posterior(x, length(hidx),  cluster, 
                                                                     datlist$X[hidx,], datlist$Z[hidx,], datlist$time[hidx], datlist$y[hidx,], datlist$K, 
                                                                     betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                                                                     gamma2[s,], gamma3[s,], alpha0[s], alpha[s,], rho[s,])) )
         prob = c(prob, mean(numer/denom))
       }
       prob_all = rbind(prob_all, lambda[s,]*prob)
       
     } # hl
     
     return(prob_all)
   }
   
   #########################################################################
   # apply across posterior draws
   library(parallel) # mclapply
   library(mvtnorm) # dmvnorm
   id = 1
   get_clst_prob(m=100, id, s=1, datlist, 
                 betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho)
   ptm_mclapply <- proc.time()  
   prob_all_draws = mclapply(1:npost,mc.cores = parallel::detectCores() / 2, 
                             function(s) get_clst_prob(m=100, id, s, datlist, 
                                                       betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho))
   proc.time() - ptm_mclapply # 
   
   
   probs = lapply(prob_all_draws, function(mat) t(apply(mat ,1, function(x) x/sum(x) ) ) )
   probs_array = simplify2array(probs)
   apply(probs_array,c(1,2),summ)[,1]
   
   
   
   ptm_mclapply <- proc.time()  
   optim(rep(0,cc),  theta_posterior, theta_gradient, "method" = "BFGS", control = list(fnscale = -1),
         cluster = cluster, 
         X = datlist$X[hidx,], Z = datlist$Z[hidx,], time = datlist$time[hidx], y = datlist$y[hidx,],yobs = datlist$y_observed[hidx,], y1on0 = datlist$y1on0[hidx],
         betaX = betaX[s,,], betat = betat[s,], sigma = sigma[s,,], 
         gamma1 = gamma1, gamma2 = gamma2[s,], gamma3 = gamma3[s,], 
         alpha0 = alpha0[s], alpha = alpha[s,], rho = rho[s,])[[1]]
   theta_V( datlist$time[hidx], rho[s,], cluster, length(hidx))
   proc.time() - ptm_mclapply # 
   
   #tmp = calc_prob_persondraw(m=100, idx, s = 1,
   #                 datlist$X, datlist$Z, datlist$time, datlist$y, datlist$K,
   #                 betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho, lambda)
   
   #Rcpp::sourceCpp("cfun_AD.cpp")
   #tmp = calc_prob_persondraw(m=10, idx, s = 1,
   #                     datlist$X, datlist$Z, datlist$time, datlist$y, datlist$K,
   #                     betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho,lambda)
   
 }
