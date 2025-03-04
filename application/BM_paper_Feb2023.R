### 10/19/2024
### add calculation of DIC


type = 'B'; ageopt = 2 # 

source("data_processing_Feb2023_noWMH.R")

a = 3;b=1
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
                 y1on0 = 1*(!is.na(y[,1]) & y[,1] == min(y[,1],na.rm=T)), 
                 a = a,b=b,
                 catlen = catlen,
                 catsum = cumsum(c(0,catlen))[1:3])


### dic for 10 chains


calc_likelihood = function(staninput, L,
                           betaXj, betatj, gamma1, gamma2j, gamma3j, dscorej, Sigj, lambdaj){
  
  K = staninput$K; Nobs = staninput$Nobs; N = staninput$Nsubs
  mu1 = cbind(staninput$X, staninput$time) %*% rbind(betaXj, betatj) # Nobs x K
  
  mu_tmp = lapply(1:K, function(k) gamma1[k]/(1 + exp( -gamma2j[k] * (dscorej - gamma3j[k])) ))
  mu = array( unlist( mu_tmp ) , dim = c( Nobs, L, K ))
  
  yloglik = array(NA, dim = c( Nobs, L, K ))
  person_loglik = matrix(NA, nrow = N, ncol = L)
  for(ell in 1:L){
    mu[,ell,] = mu[,ell,] + mu1
    yloglik[,ell,] =  dnorm(y, mean = mu[,ell,], 
                            sd = matrix(sqrt(Sigj[ell,]), nrow=1) %x% matrix(1, nrow = Nobs, ncol=1), 
                            log = TRUE) # Nobs x K 
    yloglik[,ell,][staninput$y_observed == 0] = 0
    tmp = pnorm(y[,1], mean = mu[,ell,1], 
                sd = matrix(sqrt(Sigj[ell,1]), nrow=1) %x% matrix(1, nrow = Nobs, ncol=1), 
                log = TRUE) # Nobs x 1 
    yloglik[staninput$y1on0 == 1,ell,1] = tmp[staninput$y1on0 == 1]
    
    library(data.table)
    dttmp = cbind(yloglik[,ell,], ID = rep(1:N, staninput$lsub) ); 
    colnames(dttmp)[1:K] = paste0("Y", 1:K); dttmp = as.data.frame(dttmp)
    setDT(dttmp)
    if(L==1){
      person_loglik[,ell] = dttmp[, sum(.SD), by=ID ]$V1
    } else {
      person_loglik[,ell] = dttmp[, sum(.SD), by=ID ]$V1+ log(lambdaj[ell])
    }
    
  }
  
  if(L==1){
    loglik = sum(person_loglik)
  } else {
    loglik = sum(log(apply(exp(person_loglik), 1, sum)))
  }
  
  return(loglik)
}


calc_dic = function(pidx, gamma1, res, L, staninput){
  dn = dimnames(res)
  
  K = staninput$K; Nobs = staninput$Nobs; N = staninput$Nsubs
  
  if(L==1)
    loglikvec = lapply(pidx, function(pj) calc_likelihood(staninput, L,
                                                          t(res$betaX[pj,,]), 
                                                          matrix(res$betat[pj,], nrow = 1), gamma1, 
                                                          res$gamma2[pj,], res$gamma3[pj,], 
                                                          matrix(res$d[pj,], ncol = L), 
                                                          matrix(res$Sig[pj,],nrow=1), res$lambda[pj,] ))
  if(L>1)
    loglikvec = lapply(pidx, function(pj) calc_likelihood(staninput, L,
                                                          t(res$betaX[pj,,]), 
                                                          matrix(res$betat[pj,], nrow = 1), gamma1, 
                                                          res$gamma2[pj,], res$gamma3[pj,], 
                                                          res$d[pj,,], 
                                                          res$Sig[pj,,] , res$lambda[pj,] ))
  loglikvec = unlist(loglikvec)
  dic = (-2) * mean(loglikvec) + var(loglikvec)/2
  return(dic)
}

dic = rep(NA,3)
for(L in 1:3){
  fn = paste0(paste0("BM",L,"_opt2_Feb2023_noWMH.RData")) # BM2_opt2_Feb2023.R
  if(file.exists(fn)) load(fn)
  
  library(rstan)
  res = extract(fit, permuted =T,inc_warmup=F)
  
  dic[L] =  calc_dic(pidx = 1:4000, gamma1, res, L, staninput)
  rm(fit)
} 
#> dic
#[1]  73301.43 118813.61 141365.66

### add calculation of diagnostic score
### PSRF does not work
# library(rstan)
# pidx = 1:400
# for(L in 1:3){
#   fn = paste0(paste0("BM",L,"_opt2_Feb2023_noWMH.RData")) # BM2_opt2_Feb2023.R
#   if(file.exists(fn)) load(fn)
#   
#   if(L == 1){
#     allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha", "rho_theta", "d", "lp__")
#   } else {
#     allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha0", "alpha", "rho_theta", "d", "lp__", "lambda")
#   }
#   
#   psrf_list = vector("list", length(allnames))
#   names(psrf_list) = allnames
#   
#   res = extract(fit, permuted =F,inc_warmup=F)
#   dn = dimnames(res)
#   
#   vn = c("gamma2","gamma3","betaX","betat")
#   for(j in 1:length(vn)){
#     ch = res[pidx, , grep(vn[j],dn$parameters )]
#     library(EDISON)
#     samples = vector("list", 10)
#     for(kk in 1:10){
#       samples[[kk]] = t(ch[,kk,]) 
#     }
#   
#     psrf_list[[j]] = rbind(psrf_list[[j]],   psrf(samples))
#   } # j
#   save(psrf_list, file = paste0("postsumm_PSRF_", L,".RData"))
#  
#   rm(fit)
# } 

### Geweke for the posterior mode chain?
setwd("./AD")
library(coda) #geweke
library(rstan)
#pidx = 1:400

mygeweke = function(chain, f1, f2){
  if(is.null(ncol(chain))) chain = matrix(chain, ncol =1)
  rr = tryCatch(geweke.diag(chain, frac1=f1, frac2=f2)$z, error = function(e){return(NA)} )
  if(sum(is.na(rr)) > 0 & ncol(chain) > 1){
    rr  = rep(NA, ncol(chain))
    for(v in 1:ncol(chain)){
      rr[v]  = tryCatch(geweke.diag(chain[,v], frac1=f1, frac2=f2)$z, error = function(e){return(NA)} )
    }
  }
  return(rr)
}


p1 = 0.45; p2 = 0.05

for(L in 1:3){
  fn = paste0(paste0("BM",L,"_opt2_Feb2023_noWMH.RData")) # BM2_opt2_Feb2023.R
  if(file.exists(fn)) load(fn)

  res = extract(fit, permuted =T,inc_warmup=F) # odch[od] = 1
  dn = dimnames(res)
  
  pidx = 1:4000
  
  conv2 = c()#matrix(NA, ncol = K, nrow = 6+L+1+L-1+L+L)
  cname = c()
  
  vn = "gamma2"
  ch = res$gamma2[pidx, ]
  conv2 = rbind(conv2, mygeweke(ch, p1, p2))
  cname = c(cname, vn)
  
  vn = "gamma3"
  ch = res$gamma3[pidx, ]
  conv2 = rbind(conv2, mygeweke(ch, p1, p2))
  cname = c(cname, vn)
  
  vn = "betaX"
  ch = res$betaX[pidx, ,]
  for(kk in 1:4){
    conv2 = rbind(conv2, mygeweke(ch[,,kk], p1, p2) )
  }
  cname = c(cname, c("beta0","betaApoe","betaSex","betaEdu"))
  
  vn = "betat"
  ch = res$betat[pidx, ]
  conv2 = rbind(conv2, mygeweke(ch, p1, p2))
  cname = c(cname, vn)
  
  vn = "Sig"
  if(L==1){
    ch = res$Sig[pidx, ]
    conv2 = rbind(conv2, mygeweke(ch, p1, p2))
  } else {
    ch = res$Sig[pidx,, ]
    for(kk in 1:L){
      conv2 = rbind(conv2, mygeweke(ch[,kk,], p1, p2) )
    }
  }
  cname = c(cname, paste0("Sig", 1:L))
  
  rownames(conv2) = cname
  
  
  conv2_latent = c()
  cname = c()
  
  # alpha0
  if(L > 1){
    if(L==2) ch = as.numeric(res$alpha0[pidx])
    if(L==3) ch = res$alpha0[pidx,]
    conv2_latent = c(conv2_latent, mygeweke(ch, p1, p2))
    cname = c(cname, paste0("alpha0_", 2:L))
  }
  
  # alpha (apoe, age, apoe x age)
  ch = res$alpha[pidx,]
  conv2_latent = c(conv2_latent, mygeweke(ch, p1, p2))
  cname = c(cname, "alphaApoe", "alphaAge", "alphaIntrtn")
  
  # rho
  if(L==1) ch = as.numeric(res$rho_theta[pidx])
  if(L==2) ch = cbind(res$rho_theta1[pidx], res$rho_theta2[pidx])
  if(L==3) ch = res$rho_theta[pidx,]  
  conv2_latent = c(conv2_latent, mygeweke(ch, p1, p2))
  cname = c(cname, paste0("rho",1:L))
  
  # lambda
  if(L > 1){
    ch = res$lambda[pidx,]
    conv2_latent = c(conv2_latent, mygeweke(ch, p1, p2))
    cname = c(cname, paste0("lambda",1:L))
  }
  names(conv2_latent) = cname 
  
  conv2; conv2_latent
  
  
  rm(fit, res)
  save(conv2, conv2_latent, file = paste0("ApplicationGeweke_", L,".RData"))
}# L

# create tables for replying reviewers

L=2
load(paste0("ApplicationGeweke_", L,".RData"))
colnames(conv2) = paste0("Y", 1:11)
round(conv2, 2)

L=1
load(paste0("ApplicationGeweke_", L,".RData"))
colnames(conv2) = paste0("Y", 1:11)
round(conv2, 2)

L=3
load(paste0("ApplicationGeweke_", L,".RData"))
colnames(conv2) = paste0("Y", 1:11)
round(conv2, 2)

### visualize biomarker residual correlation

calc_residual = function(staninput, L,
                           betaXj, betatj, gamma1, gamma2j, gamma3j, dscorej, Sigj, lambdaj){
  
  K = staninput$K; Nobs = staninput$Nobs; N = staninput$Nsubs
  mu1 = cbind(staninput$X, staninput$time) %*% rbind(betaXj, betatj) # Nobs x K
  
  mu_tmp = lapply(1:K, function(k) gamma1[k]/(1 + exp( -gamma2j[k] * (dscorej - gamma3j[k])) ))
  mu = array( unlist( mu_tmp ) , dim = c( Nobs, L, K ))
  
  resid = array(NA, dim = c( Nobs, L, K ))
  
  for(ell in 1:L){
    mu[,ell,] = mu[,ell,] + mu1
    resid[,ell,] = y - mu[,ell,]
    resid[staninput$y1on0 == 1, ell, 1] = NA # those with ceiling effect
    resid[,ell,][staninput$y_observed == 0] = NA # missing Y
  }
  return(resid)
}


summ = function(x){
  rdnum = 2
  x = x[!is.na(x)]
  sx = sort(x)
  n = length(x)
  r = c(round(mean(x),rdnum), round(sx[ceiling(n*0.025)],rdnum), round(sx[floor(n*0.975)],rdnum))
  return(r)
}

calc_resid_corr = function(pidx, gamma1, res, L, cidx, staninput){

  K = staninput$K; Nobs = staninput$Nobs; N = staninput$Nsubs
  
  rr = lapply(pidx, function(pj) calc_residual(staninput, L,
                                            t(res$betaX[pj,,]), 
                                            matrix(res$betat[pj,], nrow = 1), gamma1, 
                                            res$gamma2[pj,], res$gamma3[pj,], 
                                            res$d[pj,,], 
                                            res$Sig[pj,,] , res$lambda[pj,] ))
  
  rrarray = simplify2array(rr)#[1:3290, 1:2, 1:11, 1:4000]; NxLxKxdraws
  
  library(WGCNA)
  show1 = lapply(pidx, function(pj){
    cmat = corAndPvalue(rrarray[,1,,pj], method='spearman')$cor
    nn = ncol(cmat)
    #idx = unlist(lapply(1:(nn-1), function(j) j*nn+1:j)) # upper matrix 
    idx = unlist(lapply(1:(nn-1), function(j) (j-1)*nn+(j+1):nn)) # lower matrix
    return(cmat[idx])
  } )
  show1 = show
  
  show2 = lapply(pidx, function(pj){
    cmat = corAndPvalue(rrarray[,2,,pj], method='spearman')$cor
    nn = ncol(cmat)
    #idx = unlist(lapply(1:(nn-1), function(j) j*nn+1:j))
    idx = unlist(lapply(1:(nn-1), function(j) (j-1)*nn+(j+1):nn))
    return(cmat[idx])
  } )
  
  show1array = simplify2array(show1); show2array = simplify2array(show2)
  summlist = list(apply(show1array, 1, summ), apply(show2array, 1, summ))
  return(summlist)
}

# make scatterplot

L=2
fn = paste0(paste0("BM",L,"_opt2_Feb2023_noWMH.RData")) # BM2_opt2_Feb2023.R
if(file.exists(fn)) load(fn)

res = extract(fit, permuted =T,inc_warmup=F) # odch[od] = 1
dn = dimnames(res)

pidx = 1:4000

K = staninput$K; Nobs = staninput$Nobs; N = staninput$Nsubs

rr = lapply(pidx, function(pj) calc_residual(staninput, L,
                                             t(res$betaX[pj,,]), 
                                             matrix(res$betat[pj,], nrow = 1), gamma1, 
                                             res$gamma2[pj,], res$gamma3[pj,], 
                                             res$d[pj,,], 
                                             res$Sig[pj,,] , res$lambda[pj,] ))

rrarray = simplify2array(rr)#[1:3290, 1:2, 1:11, 1:4000]; NxLxKxdraws
rrmean = apply(rrarray, c(1,2,3), function(x) mean(x, na.rm=T))

for(l in 1:2){
  pdl = rrmean[,l,]
  #paperall = c("MMSE", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", TeX('P-tau$_{181p}$'), TeX("A$\\beta_{42}$/A$\\beta_{40}$") )
  paperall = c("MMSE", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", "P-tau", "AbetaRatio" )
  pdl = as.data.frame(pdl); colnames(pdl) = paperall
  pdf(paste0("resid_scatter_",l,".pdf"))
  pairs(pdl, pch = 19,cex = 0.2)
  dev.off()
}



# visualize

library("ggplot2")

pd = as.data.frame(cbind(1:ncol(summlist[[1]]), t(summlist[[1]])) )
colnames(pd) = c("x","y","lower","upper")
p1 = ggplot(pd, aes(x, y)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Outcome Correlations in Cluster 1") + xlab("")

pd = as.data.frame(cbind(1:ncol(summlist[[2]]), t(summlist[[2]])) )
colnames(pd) = c("x","y","lower","upper")
p2 = ggplot(pd, aes(x, y)) +        # ggplot2 plot with confidence intervals
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper))+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab("Outcome Correlations in Cluster 2") + xlab("")

pdf("resid_corr.pdf")
require(gridExtra)
grid.arrange(p1, p2, nrow=2)
dev.off()

# DW test
x = summlist[[1]][1,]
sum((x[2:length(x)] - x[1:(length(x)-1)])^2)/sum(x^2) 
#0.85 # by upper matrix
#0.85 # by lower matrix
x = summlist[[2]][1,]
sum((x[2:length(x)] - x[1:(length(x)-1)])^2)/sum(x^2) 
#0.87 # by upper matrix
#1.05 # by lower matrix

### visualize qqplot of each transformed biomarker

type = 'B'; ageopt = 2 # 

route = ""
source(paste0(route, "PlotFunctions.R"))

source(paste0(route, "data_processing_Feb2023_noWMH.R"))

# ratio of on-boundary LOD Y1
tmp = 1*(!is.na(y[,1]) & y[,1] == min(y[,1],na.rm=T))
round(table(tmp)/length(tmp),2) # 0 0.44 1 0.56

# qqplots
require(qqplotr)
ytmp = y
ytmp[tmp,1] = NA
pd = data.frame(obs = c(ytmp), index = rep(1:ncol(y), each = nrow(y)) )
pd$index = as.factor(pd$index)
gg <- ggplot(data = pd, mapping = aes(sample = obs, color = index, fill = index)) +
  stat_qq_band(alpha=0.5) +
  stat_qq_line() +
  stat_qq_point() +
  facet_wrap(~ index) +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles")
pdf("qqploty.pdf")
gg
dev.off()

pdf("tmp.pdf") # Y1 not normal univariately (excluding celings)
qqnorm(y[tmp==0,1])
dev.off()

library(rstatix)
tmp = 1*(!is.na(y[,1]) & y[,1] == min(y[,1],na.rm=T))
round(table(tmp)/length(tmp),2) # 0 0.44 1 0.56
ytmp = y
ytmp[tmp,1] = NA
mshapiro_test(ytmp)

###########################################################

type = 'B'; ageopt = 2 # 

route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/"
source(paste0(route, "PlotFunctions.R"))

source(paste0(route, "data_processing_Feb2023_noWMH.R"))

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

### data summary of BIOCARD

# age at recruitment
range(as.numeric(d[, ageori[1], by=Study.ID]$V1)) # 20, 86

# education
range(as.numeric(d[, eduori[1], by=Study.ID]$V1)) # 12, 20

# visit gap
summary(as.numeric(d[, ageori - c(NA, ageori[-.N]), by=get(id)]$V1)) # IQR 1 (0.99, 1.15)

# baseline missingness
tmp = d[, .SD[1], by = Study.ID]
yinit = tmp[, all,with=F]
apply(yinit,2, function(x) mean(is.na(x)))

# conversion time: last diagnosis ascertained by last available diagnosis
d1 = d[!is.na(dx),] # length(unique(d1$RID)) == length(unique(d$RID))
fstS = d1[,dx[1],by=Study.ID]$V1
lstS =  d1[,dx[.N],by=Study.ID]$V1
fstage = d1[, ageori[1], by=Study.ID]$V1

firsttime = function(d1, tag){
  dd = copy(d1)
  dd[, tmp := cumsum(!is.na(dx) & dx == tag), by=Study.ID]
  dd[, target := min(1000, ageori[tmp==1]),by=Study.ID]
  rt = dd[, target[1], by=Study.ID]$V1
  rt[rt==1000] = NA
  return(rt)
}
MCIage = firsttime(d1, "MCI")
ADage = firsttime(d1, "AD")

covtype = 0
covtype[fstS == "NORMAL" & lstS == "MCI"] = "NM-MCI"
covtype[fstS == "NORMAL" & lstS == "AD"] = "NM-AD"
covtype[fstS == "NORMAL" & lstS == "NORMAL"] = "NM-NM"
table(covtype)
NM1 = 1*(fstS == "NORMAL" & lstS == "MCI"); NA1 = 1*(fstS == "NORMAL" & lstS == "AD"); 
NMt1 = MCIage - fstage; NMt1[NM1==0] = NA
NAt1 = ADage - fstage; NAt1[NA1==0] = NA

# age among those converted to MCI/AD
agevec = as.numeric(d[,ageori[1], Study.ID]$V1)
summ(agevec[covtype == "NM-NM"])
summ(agevec[covtype == "NM-MCI"])
summ(agevec[covtype == "NM-AD"])

# length of follow-up, baseline characteristics (age, apoe, SEX, education) and baseline biomarkers for all, NM-NM, NM-MCI, and NM-AD
tab = matrix(NA, nrow = 16, ncol = 4)
rownames(tab) = c("conversion time", "baseline age", "apoe","SEX","education",all)
colnames(tab) = c("all", "NM-NM", "NM-MCI", "NM-AD")
intab = function(v){
  if(length(unique(v)) == 2){
    round(c(mean(v,na.rm=T),mean(v[covtype == "NM-NM"],na.rm=T),mean(v[covtype == "NM-MCI"],na.rm=T),mean(v[covtype == "NM-AD"],na.rm=T)),2)
  } else {
    c(summ(v), summ(v[covtype == "NM-NM"]), summ(v[covtype == "NM-MCI"]), summ(v[covtype == "NM-AD"]))
  }
}
# conversion time
tab[1, 3] = summ(NMt1)# 12.61 (4.25, 23.12)
tab[1, 4] = summ(NAt1)# 13.45 (3.01, 19.17)
# baseline age
v = as.numeric(d[, ageori[1], Study.ID]$V1)
tab[2,]  = intab(v)
# baseline covariates
for(j in 3:16){
  v = as.numeric(d[, get(rownames(tab)[j])[1], Study.ID]$V1)
  tab[j,]  = intab(v)
}
noquote(tab)
#d[d$Study.ID %in% unique(d$Study.ID)[which(v==0)],]

##################################################################
##################################################################
##################################################################
#BIC = -2 * LL + log(N) * k
#Where LL is the log-likelihood of the model, N is the number of examples in the training dataset, and k is the number of parameters in the model.

BIC = ll = q = qre = qfix = rep(NA,3)
od = 1; idod = 1
type = "lp__"
n = nrow(d)
N = length(unique(d$Study.ID))
for(L in 1:3){
  fn = paste0(paste0("/Users/yizhenxu/BM",L,"_opt2_Feb2023_noWMH.RData")) # BM2_opt2_Feb2023.R
  if(file.exists(fn)) load(fn)
  
  qfix[L] = (ncol(x) + 1 + L + 2) * length(all)  # betaX, betat, sig, gamma2&3
  qre[L] = ncol(z) + (L-1) + L # alpha, alpha0, rho
  q[L] = (L-1) + qfix[L] + qre[L] # lambda, fix params, re params 
  
  res = extract(fit, permuted =F,inc_warmup=F)
  dn = dimnames(res)
  lp = res[,,which(dn$parameters == type)]
  odch = order(apply(lp, 2, mean), decreasing = T)
  maxidx = order(lp[, odch[od]],decreasing=T)[idod]
  ll[L] = res[maxidx, odch[od], grep(type,dn$parameters )]
 
  rm(fit)
} 

for(L in 1:3)
  BIC[L] =  -2 * ll[L]+ log(N) * (q[L])
  #BIC[L] =  -2 * ll[L]+  log(N) * qre[L] + log(n) * ((L-1)+qfix[L])

# 56442.14 55012.24 57183.97

# no WMH, updated to add DSBACK
#57418.10 55911.54 58246.37

BIC


##################################################################
##################################################################
##################################################################
### estimated table for application (L = 2)
source(paste0(route, "PlotFunctions.R"))
library(rstan)

type = "B"; ageopt = 2
source(paste0(route, "data_processing_Feb2023_noWMH.R"))

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

fn = paste0("/Users/yizhenxu/BM2_opt2_Feb2023_noWMH.RData")
if(file.exists(fn)) load(fn)

od = 1; L = 2;  nY = length(all)
res = extract(fit, permuted =F,inc_warmup=F)
dn = dimnames(res)
lp = res[,,which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)


require(latex2exp)
paperall = c("MMSE", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", TeX('P-tau$_{181p}$'), TeX("A$\\beta_{42}$/A$\\beta_{40}$") )

pdf("p1_meeting.pdf", height = 5.86, width = 5.68)
#showf_chain_AM2_paper_v12_BIOCARD_meeting("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1) # output CI in the order of Daisy's work
showf_chain_AM2_paper_v12_BIOCARD("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1)
dev.off()

### tables
ch = cbind(res[, odch[od], grep("betat",dn$parameters )],
           res[, odch[od], grep("betaX",dn$parameters )],
           res[, odch[od], grep(paste0("Sig\\[", 1),dn$parameters )],
           res[, odch[od], grep(paste0("Sig\\[", 2),dn$parameters )],
           res[, odch[od], grep("gamma2",dn$parameters )],
           res[, odch[od], grep("gamma3\\[",dn$parameters )])
xname = c("beta_t", "intercept", "beta_apoe","beta_sex","beta_edu",
          "Sig1","Sig2","gamma2", "gamma3")
nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab =  matrix(tmp[1:(nY*nx)], ncol =nY, byrow=T)
colnames(tab) =  c("MMSE Score", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", "P-tau", "Abetaratio")

rownames(tab) = xname

tabr1 = tab[c(2,1,3:5),]
tabr2 = tab[6:9,]
noquote(t(tabr1))
noquote(cbind(t(tabr2),""))

library(xtable)
xtable(t(tabr1))
xtable(cbind(t(tabr2),""))

ch = cbind(res[, odch[od], grep("alpha0",dn$parameters )],
           res[, odch[od], grep("alpha\\[",dn$parameters )],
           res[, odch[od], grep("rho_theta",dn$parameters )],
           res[, odch[od], grep("lambda\\[",dn$parameters )])
xname = c( "alpha0(2)", "alpha_apoe", "alpha_age", "alpha_apoexage","rho(1)", "rho(2)","lambda(1)", "lambda(2)")

nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab = as.data.frame(matrix(tmp,nrow=1) )
colnames(tab) =  xname
tab
xtable(tab)
##################################################################
##################################################################
##################################################################
### L = 1 results

source(paste0(route, "PlotFunctions.R"))
library(rstan)

type = "B"; ageopt = 2
source(paste0(route, "data_processing_Feb2023_noWMH.R"))

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

fn = paste0("/Users/yizhenxu/BM1_opt2.RData")
if(file.exists(fn)) load(fn)

od = 1; L = 1; titleshow = "lp__"; nY = length(all); g1 = NA; ptitle=""
res = extract(fit, permuted =F,inc_warmup=F)
dn = dimnames(res)
lp = res[,,which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)

gamma2 = res[, odch[od], grep("gamma2",dn$parameters )]
gamma3 = res[, odch[od], grep("gamma3",dn$parameters )]

npost = nrow(gamma2)

rds = res[, odch[od], grep("d\\[",dn$parameters , value=T)]

rg = range(rds)

# start plot

ds = seq(rg[1], rg[2],0.01)
f = array(NA, dim=c(length(ds), npost, nY)) # plot points index n, posterior index l, Y index k

if(all(is.na(g1))){
  gamma1plot = rep(1,nY)
} else {
  gamma1plot = g1
}

for(l in 1:npost){
  for(k in 1:nY){
    for(n in 1:length(ds)){
      f[n,l,k] = gamma1plot[k]/(1+exp(-gamma2[l,k]*(ds[n]-gamma3[l,k])))
    }
  }
}

tmp = apply(f, c(1,3), summ_sep)
loglik = res[, odch[od], grep(titleshow,dn$parameters )]
nMRI = nY - 6
colv = rep(c("black","red","green"), c(3,nMRI,3)); ltyv = c(1:3,1:nMRI,1:3)

###

pdf("p2.pdf", height = 4.5, width = 7.1)

par(mfrow=c(2,1))

# plot curve
par(mar=c(0,5,3,3))
if(titleshow == "log_lik"){
  plot(ds, rep(NA,length(ds)), ylim = range(f), xlim = rg,
       main = paste0(ptitle, paste0("log likelihood is ",summ(loglik) )),ylab = "f", xaxt="n",las=1)
}
if(titleshow == "lp__"){
  plot(ds, rep(NA,length(ds)), ylim = range(f),ylab = "", xaxt="n",las=1, xlab = "latent score")
}
for(k in 1:nY){
  lines(ds, tmp[1,,k],col=colv[k], lty = ltyv[k])
  polygon(c(rev(ds), ds), c(rev(tmp[2, ,k]), tmp[3,,k]), border = NA,col = adjustcolor(colv[k],alpha.f=0.1) )
  ## intervals
  lines(ds, tmp[2, ,k], lty = ltyv[k], col = adjustcolor(colv[k],alpha.f=0.3))
  lines(ds, tmp[3, ,k], lty = ltyv[k], col = adjustcolor(colv[k],alpha.f=0.3))
}  
legend("topleft", legend=all, col= colv, lty= ltyv, cex=0.5)

# plot histograms
par(mar=c(5,5,0.2,3))
hlist = c()
for(i in 1:npost){
  a = hist(rds[i,], breaks = unique(c(rg[1],seq(rg[1], rg[2], 0.1),rg[2])), plot=F)
  hlist = cbind(hlist, a$density)
}
hall =  apply(hlist, 1, summ_sep) # rows being mean, l and u for 95% CB, length(breaks) cols
clist = list(rgb(173,216,230,max = 255, alpha = 100, names = "lt.blue"),
             rgb(255,192,203, max = 255, alpha = 100, names = "lt.pink"),
             rgb(204,204,255, max = 255, alpha = 100, names = "lt.purple"))
cvec = c("blue", "red", "purple")
# empty plot 
plot(a$breaks, c(hall[3,],0), ylim = rev(range(hall)), type="n", las=1 , main="", ylab="", xlab = "latent score")  
tmp1 = rbind(processband(x = a$breaks,y=c(hall[2,],0)),
            processband(x = a$breaks,y=c(hall[3,],0))[(length(a$breaks)*2+1):2,])
polygon(tmp1[,1], tmp1[,2], border=NA, col=clist[[1]])
lines(a$breaks, c(hall[1,],0),type="s",col=cvec[1],lty=1)  
dev.off()

### tables
ch = cbind(res[, odch[od], grep("betat",dn$parameters )],
           res[, odch[od], grep("betaX",dn$parameters )],
           res[, odch[od], grep("Sig",dn$parameters )],
           res[, odch[od], grep("gamma2",dn$parameters )],
           res[, odch[od], grep("gamma3\\[",dn$parameters )])
xname = c("beta_t", "intercept", "beta_apoe","beta_sex","beta_edu",
          "Sig","gamma2", "gamma3")
nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab =  matrix(tmp[1:(nY*nx)], ncol =nY, byrow=T)
colnames(tab) =  all
rownames(tab) = xname
noquote(t(tab))


ch = cbind(res[, odch[od], grep("alpha\\[",dn$parameters )],
           res[, odch[od], grep("rho_theta",dn$parameters )])
xname = c( "alpha_apoe", "alpha_age", "alpha_apoexage","rho")

nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab = as.data.frame(matrix(tmp,nrow=1) )
colnames(tab) =  xname

##################################################################
##################################################################
##################################################################
### L = 3 results

source(paste0(route, "PlotFunctions.R"))
library(rstan)

type = "B"; ageopt = 2
source(paste0(route, "data_processing_Feb2023_noWMH.R"))

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

fn = paste0("/Users/yizhenxu/BM3_opt2.RData")
if(file.exists(fn)) load(fn)

od = 1; L = 3; titleshow = "lp__"; nY = length(all); g1 = NA; ptitle=""
res = extract(fit, permuted =F,inc_warmup=F)
dn = dimnames(res)
lp = res[,,which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)

gamma2 = res[, odch[od], grep("gamma2",dn$parameters )]
gamma3 = res[, odch[od], grep("gamma3",dn$parameters )]

npost = nrow(gamma2)
# membership by draws
postzch = lapply(1:dim(res)[1], 
                 function(j){
                   logpyz = matrix(res[j, odch[od], grep("logpyz",dn$parameters )],ncol=L)
                   postz = matrix(NA, nrow = nrow(logpyz), ncol = L)
                   for (n in 1:nrow(logpyz)) {
                     logdenom = log(sum(exp(logpyz[n,])))
                     postz[n,] = exp(logpyz[n,] - logdenom)
                   }
                   return(postz)
                 })
postzch = simplify2array(postzch)
membch = apply(postzch, c(1,3), function(x) which.max(x))
membchall = lapply(1:ncol(membch), function(i) rep(membch[,i],d[,.N,by=Study.ID]$N))
membchall = matrix(unlist(membchall), ncol = npost)

rds = res[, odch[od], grep("d\\[",dn$parameters , value=T)]
rds = simplify2array(lapply(1:nrow(rds), function(i) matrix(rds[i,], ncol=L)))

rg = range(rds)

# start plot

ds = seq(rg[1], rg[2],0.01)
f = array(NA, dim=c(length(ds), npost, nY)) # plot points index n, posterior index l, Y index k

if(all(is.na(g1))){
  gamma1plot = rep(1,nY)
} else {
  gamma1plot = g1
}

for(l in 1:npost){
  for(k in 1:nY){
    for(n in 1:length(ds)){
      f[n,l,k] = gamma1plot[k]/(1+exp(-gamma2[l,k]*(ds[n]-gamma3[l,k])))
    }
  }
}

tmp = apply(f, c(1,3), summ_sep)
loglik = res[, odch[od], grep(titleshow,dn$parameters )]
nMRI = nY - 6
colv = rep(c("black","red","green"), c(3,nMRI,3)); ltyv = c(1:3,1:nMRI,1:3)

###

pdf("p3.pdf", height = 4.5, width = 7.1)
par(mfrow=c(2,1))

# plot curve
par(mar=c(0,5,3,3))
if(titleshow == "log_lik"){
  plot(ds, rep(NA,length(ds)), ylim = range(f), xlim = rg,
       main = paste0(ptitle, paste0("log likelihood is ",summ(loglik) )),ylab = "f", xaxt="n",las=1)
}
if(titleshow == "lp__"){
  plot(ds, rep(NA,length(ds)), ylim = range(f),ylab = "", xaxt="n",las=1, xlab = "latent score")
}
for(k in 1:nY){
  lines(ds, tmp[1,,k],col=colv[k], lty = ltyv[k])
  polygon(c(rev(ds), ds), c(rev(tmp[2, ,k]), tmp[3,,k]), border = NA,col = adjustcolor(colv[k],alpha.f=0.1) )
  ## intervals
  lines(ds, tmp[2, ,k], lty = ltyv[k], col = adjustcolor(colv[k],alpha.f=0.3))
  lines(ds, tmp[3, ,k], lty = ltyv[k], col = adjustcolor(colv[k],alpha.f=0.3))
}  
legend("topleft", legend=all, col= colv, lty= ltyv, cex=0.5)

# plot histograms
par(mar=c(5,5,0.2,3))
hlist = vector("list",L)
for(i in 1:npost){
  memb = membchall[,i]
  for(lev in 1:L){
    a = hist(rds[memb==lev, lev, i], breaks = unique(c(rg[1],seq(rg[1], rg[2], 0.1),rg[2])), plot=F)
    hlist[[lev]] = cbind(hlist[[lev]], a$density)
  }
}
hall = lapply(hlist, function(x) apply(x, 1, summ_sep)) # lev_th element rows being mean, l and u for 95% CB
clist = list(rgb(173,216,230,max = 255, alpha = 100, names = "lt.blue"),
             rgb(255,192,203, max = 255, alpha = 100, names = "lt.pink"),
             rgb(204,204,255, max = 255, alpha = 100, names = "lt.purple"))
cvec = c("blue", "red", "purple")
# empty plot 
plot(a$breaks, c(hall[[1]][3,],0), ylim = rev(range(unlist(hall))), type="n", las=1 , main="", ylab="", xlab = "latent score")  
for(lev in 1:L){
  tmp1 = rbind(processband(x = a$breaks,y=c(hall[[lev]][2,],0)),
               processband(x = a$breaks,y=c(hall[[lev]][3,],0))[(length(a$breaks)*2+1):2,])
  polygon(tmp1[,1], tmp1[,2], border=NA, col=clist[[lev]])
  lines(a$breaks, c(hall[[lev]][1,],0),type="s",col=cvec[lev],lty=1)  
}
dev.off()

### tables
ch = cbind(res[, odch[od], grep("betat",dn$parameters )],
           res[, odch[od], grep("betaX",dn$parameters )],
           res[, odch[od], grep(paste0("Sig\\[", 1),dn$parameters )],
           res[, odch[od], grep(paste0("Sig\\[", 2),dn$parameters )],
           res[, odch[od], grep(paste0("Sig\\[", 3),dn$parameters )],
           res[, odch[od], grep("gamma2",dn$parameters )],
           res[, odch[od], grep("gamma3\\[",dn$parameters )])
xname = c("beta_t", "intercept", "beta_apoe","beta_sex","beta_edu",
          "Sig1","Sig2","Sig3","gamma2", "gamma3")
nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab =  matrix(tmp[1:(nY*nx)], ncol =nY, byrow=T)
colnames(tab) =  all
rownames(tab) = xname
noquote(t(tab))


ch = cbind(res[, odch[od], grep("alpha0",dn$parameters )],
           res[, odch[od], grep("alpha\\[",dn$parameters )],
           res[, odch[od], grep("rho_theta",dn$parameters )],
           res[, odch[od], grep("lambda\\[",dn$parameters )])
xname = c( "alpha0(2)","alpha0(3)", "alpha_apoe", "alpha_age", "alpha_apoexage","rho(1)", "rho(2)","rho(3)","lambda(1)", "lambda(2)","lambda(3)")

nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab = as.data.frame(matrix(tmp,nrow=1) )
colnames(tab) =  xname


##################################################################
##################################################################
##################################################################
### conversion time tables (L=2)

source(paste0(route, "PlotFunctions.R"))
library(rstan)

type = "B"; ageopt = 2
source(paste0(route, "data_processing_Feb2023_noWMH.R"))

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

fn = paste0("/Users/yizhenxu/BM2_opt2_Feb2023_noWMH.RData")
if(file.exists(fn)) load(fn)

od = 1; L = 2; titleshow = "lp__"; nY = length(all); g1 = NA; ptitle=""
res = extract(fit, permuted =F,inc_warmup=F)
dn = dimnames(res)
lp = res[,,which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)

gamma2 = res[, odch[od], grep("gamma2",dn$parameters )]
gamma3 = res[, odch[od], grep("gamma3",dn$parameters )]

npost = nrow(gamma2)
# membership by draws
postzch = lapply(1:dim(res)[1], 
                 function(j){
                   logpyz = matrix(res[j, odch[od], grep("logpyz",dn$parameters )],ncol=L)
                   postz = matrix(NA, nrow = nrow(logpyz), ncol = L)
                   for (n in 1:nrow(logpyz)) {
                     logdenom = log(sum(exp(logpyz[n,])))
                     postz[n,] = exp(logpyz[n,] - logdenom)
                   }
                   return(postz)
                 })
postzch = simplify2array(postzch)
membch = apply(postzch, c(1,3), function(x) which.max(x))
membchall = lapply(1:ncol(membch), function(i) rep(membch[,i],d[,.N,by=Study.ID]$N))
membchall = matrix(unlist(membchall), ncol = npost)

d1 = d[!is.na(dx),] # length(unique(d1$RID)) == length(unique(d$RID))
fstS = d1[,dx[1],by=Study.ID]$V1
lstS =  d1[,dx[.N],by=Study.ID]$V1

fstage = d1[, ageori[1], by=Study.ID]$V1
lstage = d1[, ageori[.N], by=Study.ID]$V1

firsttime = function(d1, tag){
  dd = copy(d1)
  dd[, tmp := cumsum(!is.na(dx) & dx == tag), by=Study.ID]
  dd[, target := min(1000, ageori[tmp==1]),by=Study.ID]
  rt = dd[, target[1], by=Study.ID]$V1
  rt[rt==1000] = NA
  return(rt)
}
MCIage = firsttime(d1, "MCI")
MCIage_end = MCIage; MCIage_end[lstS == "NORMAL"] = NA # those with confirmed end diagnosis MCI/AD
ADage = firsttime(d1, "AD")
#CNage = lstage - fstage; CNage[ !is.na(MCIage) + !is.na(ADage) > 0] = NA

# type 1: last diagnosis ascertained by last available diagnosis

NM1 = 1*(fstS == "NORMAL" & lstS == "MCI"); NA1 = 1*(fstS == "NORMAL" & lstS == "AD"); 
NMt1 = MCIage - fstage; NMt1[NM1==0] = NA
NAt1 = ADage - fstage; NAt1[NA1==0] = NA

covtype = 0
covtype[fstS == "NORMAL" & lstS == "MCI"] = "NM"
covtype[fstS == "NORMAL" & lstS == "AD"] = "NA"
covtype[fstS == "NORMAL" & lstS == "NORMAL"] = "NN"
table(covtype)

### Conversion by initial to final diagnosis

# transition time posterior summary
NMNAch = lapply(1:dim(res)[1], 
                function(j){
                  c(sum(covtype ==  "NN" & membch[,j] == 1) / sum(membch[,j] == 1),
                    sum(covtype ==  "NM" & membch[,j] == 1) / sum(membch[,j] == 1),
                    sum(covtype ==  "NA" & membch[,j] == 1) / sum(membch[,j] == 1),
                    sum(covtype ==  "NN" & membch[,j] == 2) / sum(membch[,j] == 2),
                    sum(covtype ==  "NM" & membch[,j] == 2) / sum(membch[,j] == 2),
                    sum(covtype ==  "NA" & membch[,j] == 2) / sum(membch[,j] == 2),
                    NA,
                    mean(NMt1[membch[,j]==1],na.rm=T),
                    mean(NAt1[membch[,j]==1],na.rm=T),
                    NA,
                    mean(NMt1[membch[,j]==2],na.rm=T),
                    mean(NAt1[membch[,j]==2],na.rm=T))
                })
tab = apply(matrix(unlist(NMNAch), nrow = 12),1,summ)
ct = lapply(1:dim(res)[1], 
            function(j){
              c(sum(membch[,j] == 1),sum(membch[,j] == 2))
            })
ct = apply(matrix(unlist(ct), nrow = 2),1,summ); ctvec = rep("",6); ctvec[c(1,4)] = ct
out = cbind(c("1","","", "2","",""),ctvec, rep(c("CN-CN","CN-MCI","CN-AD"),2), matrix(tab, ncol = 2))
colnames(out) = c("Cluster Index" ,"Subject Count", "Transition type", "Percentage", "Duration")

noquote(out)

### Conversion from the first <status 1> to the first <status 2>
L = 2
covtype = c("CN-MCI","CN-MCIp", "MCI-AD"); lc = length(covtype)
out = matrix(NA, nrow = L*lc, ncol = 4)
colnames(out) = c("Cluster Index", "Conversion Type", "Percentage", "Duration")
out[,1] = rep("", nrow(out)); out[(0:(L-1))*lc+1,1] = 1:L
out[,2] = rep(covtype,L)
fillout = function(idx){
  output = c()
  for(tag in covtype){
    if(tag == "CN-MCI"){a = fstage[idx]; b = MCIage[idx]} 
    if(tag == "CN-MCIp"){a = fstage[idx]; b = MCIage_end[idx]} 
    if(tag == "MCI-AD"){a = MCIage[idx]; b = ADage[idx]} 
    d = b - a
    pc = sum(!is.na(d))/ sum(!is.na(a))
    dr = mean(d, na.rm = T)
    output = c(output, c(pc, dr))
  }
  return(output)
}
NMNAch = lapply(1:dim(res)[1], 
                function(j){
                  c(fillout(which(membch[,j]==1)), fillout(which(membch[,j]==2)))
                })
tab = apply(matrix(unlist(NMNAch), nrow = 6*L),1,summ)

out[,3:4] = matrix(tab, byrow = T, ncol = 2) 
noquote(out)
##################################################################
##################################################################
##################################################################
### conversion time tables (L=3)

source("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/code/PlotFunctions.R")
library(rstan)

type = "B"; ageopt = 2
source("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/code/data_processing_Feb2023_noWMH.R")

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term

fn = paste0("/Users/yizhenxu/BM3_opt2.RData")
if(file.exists(fn)) load(fn)

od = 1; L = 3; titleshow = "lp__"; nY = length(all); g1 = NA; ptitle=""
res = extract(fit, permuted =F,inc_warmup=F)
dn = dimnames(res)
lp = res[,,which(dn$parameters == "lp__")]
odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)

gamma2 = res[, odch[od], grep("gamma2",dn$parameters )]
gamma3 = res[, odch[od], grep("gamma3",dn$parameters )]

npost = nrow(gamma2)
# membership by draws
postzch = lapply(1:dim(res)[1], 
                 function(j){
                   logpyz = matrix(res[j, odch[od], grep("logpyz",dn$parameters )],ncol=L)
                   postz = matrix(NA, nrow = nrow(logpyz), ncol = L)
                   for (n in 1:nrow(logpyz)) {
                     logdenom = log(sum(exp(logpyz[n,])))
                     postz[n,] = exp(logpyz[n,] - logdenom)
                   }
                   return(postz)
                 })
postzch = simplify2array(postzch)
membch = apply(postzch, c(1,3), function(x) which.max(x))
membchall = lapply(1:ncol(membch), function(i) rep(membch[,i],d[,.N,by=Study.ID]$N))
membchall = matrix(unlist(membchall), ncol = npost)

d1 = d[!is.na(dx),] # length(unique(d1$RID)) == length(unique(d$RID))
fstS = d1[,dx[1],by=Study.ID]$V1
lstS =  d1[,dx[.N],by=Study.ID]$V1

fstage = d1[, ageori[1], by=Study.ID]$V1
lstage = d1[, ageori[.N], by=Study.ID]$V1

firsttime = function(d1, tag){
  dd = copy(d1)
  dd[, tmp := cumsum(!is.na(dx) & dx == tag), by=Study.ID]
  dd[, target := min(1000, ageori[tmp==1]),by=Study.ID]
  rt = dd[, target[1], by=Study.ID]$V1
  rt[rt==1000] = NA
  return(rt)
}
MCIage = firsttime(d1, "MCI")
ADage = firsttime(d1, "AD")
#CNage = lstage - fstage; CNage[ !is.na(MCIage) + !is.na(ADage) > 0] = NA

# type 1: last diagnosis ascertained by last available diagnosis

NM1 = 1*(fstS == "NORMAL" & lstS == "MCI"); NA1 = 1*(fstS == "NORMAL" & lstS == "AD"); 
NMt1 = MCIage - fstage; NMt1[NM1==0] = NA
NAt1 = ADage - fstage; NAt1[NA1==0] = NA

covtype = 0
covtype[fstS == "NORMAL" & lstS == "MCI"] = "NM"
covtype[fstS == "NORMAL" & lstS == "AD"] = "NA"
covtype[fstS == "NORMAL" & lstS == "NORMAL"] = "NN"
table(covtype)

# transition time posterior summary
NMNAch = lapply(1:dim(res)[1], 
                function(j){
                  c(sum(covtype ==  "NN" & membch[,j] == 1) / sum(membch[,j] == 1),
                    sum(covtype ==  "NM" & membch[,j] == 1) / sum(membch[,j] == 1),
                    sum(covtype ==  "NA" & membch[,j] == 1) / sum(membch[,j] == 1),
                    sum(covtype ==  "NN" & membch[,j] == 2) / sum(membch[,j] == 2),
                    sum(covtype ==  "NM" & membch[,j] == 2) / sum(membch[,j] == 2),
                    sum(covtype ==  "NA" & membch[,j] == 2) / sum(membch[,j] == 2),
                    sum(covtype ==  "NN" & membch[,j] == 3) / sum(membch[,j] == 3),
                    sum(covtype ==  "NM" & membch[,j] == 3) / sum(membch[,j] == 3),
                    sum(covtype ==  "NA" & membch[,j] == 3) / sum(membch[,j] == 3),
                    NA,
                    mean(NMt1[membch[,j]==1],na.rm=T),
                    mean(NAt1[membch[,j]==1],na.rm=T),
                    NA,
                    mean(NMt1[membch[,j]==2],na.rm=T),
                    mean(NAt1[membch[,j]==2],na.rm=T),
                    NA,
                    mean(NMt1[membch[,j]==3],na.rm=T),
                    mean(NAt1[membch[,j]==3],na.rm=T))
                })
tab = apply(matrix(unlist(NMNAch), nrow = 18),1,summ)
ct = lapply(1:dim(res)[1], 
            function(j){
              c(sum(membch[,j] == 1),sum(membch[,j] == 2),sum(membch[,j] == 3))
            })
ct = apply(matrix(unlist(ct), nrow = 3),1,summ); ctvec = rep("",9); ctvec[c(1,4,7)] = ct
out = cbind(c("1","","", "2","","","3","",""), ctvec, rep(c("CN-CN","CN-MCI","CN-AD"),3), matrix(tab, ncol = 2))
colnames(out) = c("Cluster Index" ,"Subject Count", "Transition type", "Percentage", "Duration")
noquote(out)

### Conversion from the first <status 1> to the first <status 2>
L = 3
out = matrix(NA, nrow = L*2, ncol = 4)
colnames(out) = c("Cluster Index", "Conversion Type", "Percentage", "Duration")
out[,1] = rep("", nrow(out)); out[(0:(L-1))*2+1,1] = 1:L
out[,2] = rep(c("CN-MCI", "MCI-AD"),L)
fillout = function(idx){
  output = c()
  for(tag in c("CN-MCI","MCI-AD")){
    if(tag == "CN-MCI"){a = fstage[idx]; b = MCIage[idx]} 
    if(tag == "MCI-AD"){a = MCIage[idx]; b = ADage[idx]} 
    d = b - a
    pc = sum(!is.na(d))/ sum(!is.na(a))
    dr = mean(d, na.rm = T)
    output = c(output, c(pc, dr))
  }
  return(output)
}
NMNAch = lapply(1:dim(res)[1], 
                function(j){
                  c(fillout(which(membch[,j]==1)), fillout(which(membch[,j]==2)), fillout(which(membch[,j]==3)))
                })
tab = apply(matrix(unlist(NMNAch), nrow = 4*L),1,summ)

out[,3:4] = matrix(tab, byrow = T, ncol = 2) 
noquote(out)
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################









npost = nrow(gamma2)
# membership by draws
postzch = lapply(1:dim(res)[1], 
                 function(j){
                   logpyz = matrix(res[j, odch[od], grep("logpyz",dn$parameters )],ncol=L)
                   postz = matrix(NA, nrow = nrow(logpyz), ncol = L)
                   for (n in 1:nrow(logpyz)) {
                     logdenom = log(sum(exp(logpyz[n,])))
                     postz[n,] = exp(logpyz[n,] - logdenom)
                   }
                   return(postz)
                 })
postzch = simplify2array(postzch)
membch = apply(postzch, c(1,3), function(x) which.max(x))
membchall = lapply(1:ncol(membch), function(i) rep(membch[,i],d[,.N,by=Study.ID]$N))
membchall = matrix(unlist(membchall), ncol = npost)

rds = res[, odch[od], grep("d\\[",dn$parameters , value=T)]
rds = simplify2array(lapply(1:nrow(rds), function(i) matrix(rds[i,], ncol=L)))

rg = range(rds)

# start plot

ds = seq(rg[1], rg[2],0.01)
f = array(NA, dim=c(length(ds), npost, nY)) # plot points index n, posterior index l, Y index k

if(all(is.na(g1))){
  gamma1plot = rep(1,nY)
} else {
  gamma1plot = g1
}

for(l in 1:npost){
  for(k in 1:nY){
    for(n in 1:length(ds)){
      f[n,l,k] = gamma1plot[k]/(1+exp(-gamma2[l,k]*(ds[n]-gamma3[l,k])))
    }
  }
}

tmp = apply(f, c(1,3), summ_sep)
loglik = res[, odch[od], grep(titleshow,dn$parameters )]
nMRI = nY - 6
colv = rep(c("black","red","green"), c(3,nMRI,3)); ltyv = c(1:3,1:nMRI,1:3)

###

pdf("p1.pdf", height = 4.5, width = 7.1)
par(mfrow=c(2,1))

# plot curve
par(mar=c(0,5,3,3))
if(titleshow == "log_lik"){
  plot(ds, rep(NA,length(ds)), ylim = range(f), xlim = rg,
       main = paste0(ptitle, paste0("log likelihood is ",summ(loglik) )),ylab = "f", xaxt="n",las=1)
}
if(titleshow == "lp__"){
  plot(ds, rep(NA,length(ds)), ylim = range(f),ylab = "", xaxt="n",las=1, xlab = "latent score")
}
for(k in 1:nY){
  lines(ds, tmp[1,,k],col=colv[k], lty = ltyv[k])
  polygon(c(rev(ds), ds), c(rev(tmp[2, ,k]), tmp[3,,k]), border = NA,col = adjustcolor(colv[k],alpha.f=0.1) )
  ## intervals
  lines(ds, tmp[2, ,k], lty = ltyv[k], col = adjustcolor(colv[k],alpha.f=0.3))
  lines(ds, tmp[3, ,k], lty = ltyv[k], col = adjustcolor(colv[k],alpha.f=0.3))
}  
legend("topleft", legend=all, col= colv, lty= ltyv, cex=0.5)

# plot histograms
par(mar=c(5,5,0.2,3))
hlist = vector("list",L)
for(i in 1:npost){
  memb = membchall[,i]
  for(lev in 1:L){
    a = hist(rds[memb==lev, lev, i], breaks = unique(c(rg[1],seq(rg[1], rg[2], 0.1),rg[2])), plot=F)
    hlist[[lev]] = cbind(hlist[[lev]], a$density)
  }
}
hall = lapply(hlist, function(x) apply(x, 1, summ_sep)) # lev_th element rows being mean, l and u for 95% CB
clist = list(rgb(173,216,230,max = 255, alpha = 100, names = "lt.blue"),
             rgb(255,192,203, max = 255, alpha = 100, names = "lt.pink"),
             rgb(204,204,255, max = 255, alpha = 100, names = "lt.purple"))
cvec = c("blue", "red", "purple")
# empty plot 
plot(a$breaks, c(hall[[1]][3,],0), ylim = rev(range(unlist(hall))), type="n", las=1 , main="", ylab="", xlab = "latent score")  
for(lev in 1:L){
  tmp1 = rbind(processband(x = a$breaks,y=c(hall[[lev]][2,],0)),
               processband(x = a$breaks,y=c(hall[[lev]][3,],0))[(length(a$breaks)*2+1):2,])
  polygon(tmp1[,1], tmp1[,2], border=NA, col=clist[[lev]])
  lines(a$breaks, c(hall[[lev]][1,],0),type="s",col=cvec[lev],lty=1)  
}
dev.off()

### tables
ch = cbind(res[, odch[od], grep("betat",dn$parameters )],
           res[, odch[od], grep("betaX",dn$parameters )],
           res[, odch[od], grep(paste0("Sig\\[", 1),dn$parameters )],
           res[, odch[od], grep(paste0("Sig\\[", 2),dn$parameters )],
           res[, odch[od], grep("gamma2",dn$parameters )],
           res[, odch[od], grep("gamma3\\[",dn$parameters )])
xname = c("beta_t", "intercept", "beta_apoe","beta_sex","beta_edu",
          "Sig1","Sig2","gamma2", "gamma3")
nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab =  matrix(tmp[1:(nY*nx)], ncol =nY, byrow=T)
colnames(tab) =  all
rownames(tab) = xname
noquote(t(tab))


ch = cbind(res[, odch[od], grep("alpha0",dn$parameters )],
           res[, odch[od], grep("alpha\\[",dn$parameters )],
           res[, odch[od], grep("rho_theta",dn$parameters )],
           res[, odch[od], grep("lambda\\[",dn$parameters )])
xname = c( "alpha0(2)", "alpha_apoe", "alpha_age", "alpha_apoexage","rho(1)", "rho(2)","lambda(1)", "lambda(2)")

nx = length(xname)
tmp = noquote(apply(ch, 2, summ))
tab = as.data.frame(matrix(tmp,nrow=1) )
colnames(tab) =  xname



##################################################################
##################################################################
##################################################################

# 1. use the posterior mean of membership probability in each cluster; postp = matrix of n ppl as rows, L clusters as cols
# 2. obtain the mode as membership; clustering = apply(postp, 1, which.max)
# 3. define distance matrix dis, between each pair of rows/subjects, of outcomes, defined as Euclidean distance
#    between (Ys1,...,YsK) and (Yi1,...,YiK) for person s and person i
# 4. sil = silhouette(clustering,dis); plot(sil)

drawmemb = function(fit, L){
  res = extract(fit, permuted =F,inc_warmup=F)
  dn = dimnames(res)
  lp = res[,,which(dn$parameters == "lp__")]
  odch = order(apply(lp, 2, mean), decreasing = T)
  od = 1
  
  postzch = lapply(1:dim(res)[1], 
                   function(j){
                     logpyz = matrix(res[j, odch[od], grep("logpyz",dn$parameters )],ncol=L)
                     postz = matrix(NA, nrow = nrow(logpyz), ncol = L)
                     for (n in 1:nrow(logpyz)) {
                       logdenom = log(sum(exp(logpyz[n,])))
                       postz[n,] = exp(logpyz[n,] - logdenom)
                     }
                     return(postz)
                   })
  postzch = simplify2array(postzch)
  pmean = apply(postzch, c(1,2), mean)
  memb = apply(pmean, 1, which.max)

  return(memb)
}


type = "B"; ageopt = 2
source("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/code/data_processing_Feb2023_noWMH.R")


clustering = matrix(NA, nrow = 297, ncol = 2)

for(L in 2:3){
  fn = paste0(paste0("/Users/yizhenxu/BM",L,"_opt2.RData"))
  if(file.exists(fn)) load(fn)
  clustering[,L-1] = drawmemb(fit, L)
} 

#dis = distfun(y)
d1 = d[!is.na(dx),] # length(unique(d1$RID)) == length(unique(d$RID))
fstS = d1[,dx[1],by=Study.ID]$V1
lstS =  d1[,dx[.N],by=Study.ID]$V1
fstage = d1[, ageori[1], by=Study.ID]$V1

firsttime = function(d1, tag){
  dd = copy(d1)
  dd[, tmp := cumsum(!is.na(dx) & dx == tag), by=Study.ID]
  dd[, target := min(1000, ageori[tmp==1]),by=Study.ID]
  rt = dd[, target[1], by=Study.ID]$V1
  rt[rt==1000] = NA
  return(rt)
}
MCIage = firsttime(d1, "MCI")
ADage = firsttime(d1, "AD")

covtype = 0
covtype[fstS == "NORMAL" & lstS == "MCI"] = "NM-MCI"
covtype[fstS == "NORMAL" & lstS == "AD"] = "NM-AD"
covtype[fstS == "NORMAL" & lstS == "NORMAL"] = "NM-NM"
table(covtype)
NM1 = 1*(fstS == "NORMAL" & lstS == "MCI"); NA1 = 1*(fstS == "NORMAL" & lstS == "AD"); 
NMt1 = MCIage - fstage; NMt1[NM1==0] = NA
NAt1 = ADage - fstage; NAt1[NA1==0] = NA

covtime = NAt1; covtime[!is.na(NMt1)] = NMt1[!is.na(NMt1)] 
covtime = as.numeric(covtime)
lfu = d[, as.numeric(ageori[.N] - ageori[1]), by = Study.ID]$V1
event = 1*(!is.na(covtime))
covtime[is.na(covtime)] = lfu[is.na(covtime)]

library(survival)
plot(survfit(Surv(covtime, event) ~ clustering[,2]), 
     xlab = "Days", 
     ylab = "Overall survival probability")

s1 = survfit(Surv(covtime, event) ~ 1)
yy = s1$surv[2:276] - s1$surv[1:275]
plot(s1$time[-1], yy)

### using y to calculate sillouette score failed
### because cannot define distance between multivariate longitudinal y between person

  library(Rcpp)
  
  cppFunction('NumericVector distfun(NumericMatrix mat){
 int n = mat.nrow();
 int ncol = mat.ncol();
 int nd = (pow(n,2) - n)/2;
 NumericVector out(nd);
 int ct = 0;
   for (int r1 = 0; r1 < n-1; r1++) {
    for (int r2 = r1+1; r2 < n; r2++) {
      double total = 0;
      for (int c = 0; c < ncol; c++) {
        if (mat(r1, c)<1000 & mat(r2, c)<1000) {
          total += pow(mat(r1, c) - mat(r2, c), 2);
        }
      }
        out(ct) = sqrt(total);
        ct += 1;
    }
  }
  return out;
}')
  
  # take within-person average for each outcome
  ytmp = d[, c(all,"Study.ID"), with=F]
  r = c()
  for(j in 1:length(all)){
    ytmp1 = copy(ytmp); ytmp1[is.na(ytmp1)] = 0
    r = cbind(r, ytmp1[, get(all[j])[.N] - get(all[j])[1], by=Study.ID]$V1)
  }
  # calculate distance between subjects' average outcomes
  ytmp = r; ytmp[is.na(ytmp)] = 2000
  dis0 = dist(ytmp)
  dis = distfun(ytmp)
  attributes(dis) = attributes(dis0)
  
  par(mfrow = c(2,1))
  sil = silhouette(clustering[,1],dis); plot(sil)
  sil = silhouette(clustering[,2],dis); plot(sil)
  

  
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ##################################################################
  ### meeting 05/16/2023 cognitive scores have cluster-specific gamma2 and gamma3
  
  ### when d contains alpha0
  source(paste0(route, "PlotFunctions.R"))
  library(rstan)
  
  type = "B"; ageopt = 2
  source(paste0(route, "data_processing_Feb2023_noWMH.R"))
  
  zname = c("apoe", "age")
  z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
  z = cbind(z, z[,1]*z[,2]) # interaction term
  
  fn = paste0("/Users/yizhenxu/BM2_opt2_Feb2023_noWMH_cog_alpha0_sdch20.RData")
  fn = "/Users/yizhenxu/Documents/BM2_opt2_Feb2023_noWMH_cog_alpha0_sdch20.RData"
  if(file.exists(fn)) load(fn)
  
  od = 1; L = 2;  nY = length(all)
  res = extract(fit, permuted =F,inc_warmup=F)
  dn = dimnames(res)
  lp = res[,,which(dn$parameters == "lp__")]
  odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)
  
  print(traceplot(fit, pars = "lp__", inc_warmup = F))
  
  require(latex2exp)
  paperall = c("MMSE", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", TeX('P-tau$_{181p}$'), TeX("A$\\beta_{42}$/A$\\beta_{40}$") )
  
  png("/Users/yizhenxu/Desktop/d33.png", width = 798, height = 349)
  #showf_chain_AM2_paper_v12_BIOCARD_meeting("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1) # output CI in the order of Daisy's work
  showf_chain_AM2_paper_v13_BIOCARD("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1)
  dev.off()
  
  for(j in 2:10){
    png(paste0("/Users/yizhenxu/Desktop/d33_",j,".png"), width = 798, height = 349)
    #showf_chain_AM2_paper_v12_BIOCARD_meeting("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1) # output CI in the order of Daisy's work
    showf_chain_AM2_paper_v13_BIOCARD("", d, paperall, all,fit, od=j, L = 2,g1=gamma1, odmean = 1)
    dev.off()
  }
  
  ### when d does not contain alpha0
  source(paste0(route, "PlotFunctions.R"))
  library(rstan)
  
  type = "B"; ageopt = 2
  source(paste0(route, "data_processing_Feb2023_noWMH.R"))
  
  zname = c("apoe", "age")
  z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
  z = cbind(z, z[,1]*z[,2]) # interaction term
  
  fn = paste0("/Users/yizhenxu/BM2_opt2_Feb2023_noWMH_cog_sdch20.RData")
  fn = "/Users/yizhenxu/Documents/BM2_opt2_Feb2023_noWMH_cog_sdch20.RData"
  if(file.exists(fn)) load(fn)
  
  od = 1; L = 2;  nY = length(all)
  res = rstan::extract(fit, permuted =F,inc_warmup=F)
  dn = dimnames(res)
  lp = res[,,which(dn$parameters == "lp__")]
  odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)
  
  print(traceplot(fit, pars = "lp__", inc_warmup = F))
  
  require(latex2exp)
  paperall = c("MMSE", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", TeX('P-tau$_{181p}$'), TeX("A$\\beta_{42}$/A$\\beta_{40}$") )
  
  png("/Users/yizhenxu/Desktop/d34.png", width = 798, height = 349)
  #showf_chain_AM2_paper_v12_BIOCARD_meeting("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1) # output CI in the order of Daisy's work
  showf_chain_AM2_paper_v14_BIOCARD("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1,plotj=NA)
  # add dsback
  #showf_chain_AM2_paper_v14_BIOCARD("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1, plotj = 3)
  dev.off()
  
  for(j in 2:10){
    png(paste0("/Users/yizhenxu/Desktop/d34_",j,".png"), width = 798, height = 349)
    showf_chain_AM2_paper_v14_BIOCARD("", d, paperall, all,fit, od=j, L = 2,g1=gamma1, odmean = 1,plotj=NA)
    dev.off()
  }
  
  ### when d does not contain alpha0 AND 
  ### gamma2 for cognitive is ordered for cognitive score for identifiability
  source(paste0(route, "PlotFunctions.R"))
  library(rstan)
  
  type = "B"; ageopt = 2
  source(paste0(route, "data_processing_Feb2023_noWMH.R"))
  
  zname = c("apoe", "age")
  z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
  z = cbind(z, z[,1]*z[,2]) # interaction term
  
  fn = paste0("/Users/yizhenxu/BM2_opt2_Feb2023_noWMH_cog_odg2_sdch20.RData")
  fn = "/Users/yizhenxu/Documents/BM2_opt2_Feb2023_noWMH_cog_odg2_sdch20.RData"
  fn = "BM2_opt2_Feb2023_noWMH_cog_odg2_sdch20.RData"
  if(file.exists(fn)) load(fn)
  
  od = 1; L = 2;  nY = length(all)
  res = rstan::extract(fit, permuted =F,inc_warmup=F)
  dn = dimnames(res)
  lp = res[,,which(dn$parameters == "lp__")]
  odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)
  
  pdf("tmp.pdf")
  print(traceplot(fit, pars = "lp__", inc_warmup = F))
  dev.off()
  
  require(latex2exp)
  paperall = c("MMSE", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", TeX('P-tau$_{181p}$'), TeX("A$\\beta_{42}$/A$\\beta_{40}$") )
  
  png("/Users/yizhenxu/Desktop/d35.png", width = 798, height = 349)
  #showf_chain_AM2_paper_v12_BIOCARD_meeting("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1) # output CI in the order of Daisy's work
  showf_chain_AM2_paper_v14_BIOCARD("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1,plotj=NA)
  # add dsback
  #showf_chain_AM2_paper_v14_BIOCARD("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1, plotj = 3)
  dev.off()
  
  for(j in 2:10){
    png(paste0("/Users/yizhenxu/Desktop/d35_",j,".png"), width = 798, height = 349)
    showf_chain_AM2_paper_v14_BIOCARD("", d, paperall, all,fit, od=j, L = 2,g1=gamma1, odmean = 1,plotj=NA)
    dev.off()
  }
  
  ### when d contains alpha0 and restrict gamma2 to be the same, %%%
  ### onnly gamma3 differ for cognitive scores
  source(paste0(route, "PlotFunctions.R"))
  library(rstan)
  
  type = "B"; ageopt = 2
  source(paste0(route, "data_processing_Feb2023_noWMH.R"))
  
  zname = c("apoe", "age")
  z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
  z = cbind(z, z[,1]*z[,2]) # interaction term
  
  fn = paste0("BM2_opt2_Feb2023_noWMH_cog_alpha0_sameg2_sdch20_1.RData")
  fn = "/Users/yizhenxu/Documents/BM2_opt2_Feb2023_noWMH_cog_alpha0_sameg2_sdch20_1.RData"
  if(file.exists(fn)) load(fn)
  
  od = 1; L = 2;  nY = length(all)
  res = extract(fit, permuted =F,inc_warmup=F)
  dn = dimnames(res)
  lp = res[,,which(dn$parameters == "lp__")]
  odch = order(apply(lp, 2, function(x) mean(x)), decreasing = T)
  
  pdf("tmp.pdf")
  print(traceplot(fit, pars = "lp__", inc_warmup = F))
  dev.off()
  
  require(latex2exp)
  paperall = c("MMSE", "Log Mem", "DSBACK", "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD","T-tau", TeX('P-tau$_{181p}$'), TeX("A$\\beta_{42}$/A$\\beta_{40}$") )
  
  #png("/Users/yizhenxu/Desktop/d36.png", width = 798, height = 349)
  png("d36.png", width = 798, height = 349)
  #showf_chain_AM2_paper_v12_BIOCARD_meeting("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1) # output CI in the order of Daisy's work
  showf_chain_AM2_paper_v15_BIOCARD("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1)
  dev.off()
  
  for(j in 2:10){
    #png(paste0("/Users/yizhenxu/Desktop/d36_",j,".png"), width = 798, height = 349)
    png(paste0("d36_",j,".png"), width = 798, height = 349)
    #showf_chain_AM2_paper_v12_BIOCARD_meeting("", d, paperall, all,fit, od=1, L = 2,g1=gamma1, odmean = 1) # output CI in the order of Daisy's work
    showf_chain_AM2_paper_v15_BIOCARD("", d, paperall, all,fit, od=j, L = 2,g1=gamma1, odmean = 1)
    dev.off()
  }
  