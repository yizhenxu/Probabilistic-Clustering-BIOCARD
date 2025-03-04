
setwd("./AD")

library(abind)
library(data.table)
# modid is the number of clusters specified in the model estimation

#source("BMsim_setting2.R")
source("BMsim11Y_setting.R")

load("simdata11.RData")

#sig1 = sqsig1^2; sig2 = sqsig2^2
sig = c(rbind(sig1, sig2))
betaX = c(cbind(beta0,betaX1,betaX2))
Nobs = N * Ni

book = expand.grid(1:3, 1:1000) 
book$id = 1:nrow(book) 
#modid = book[taskid,1]
#dataid = book[taskid,2]


CPfun = function(x,val,perc = 0.95){
  x = x[!is.na(x)]
  sx = sort(x)
  n = length(x)
  Lp = (1-perc)/2
  Up = 1 - Lp
  if(Lp == 0) Lp = 1/n
  r = c(sx[ceiling(n*Lp)], sx[floor(n*Up)])
  cp = ifelse(val >= r[1] & val <= r[2], 1, 0)
  return(cp)
}

IQR = function(x){
  r = quantile(x, probs = c(0.25,0.75))
  return(as.numeric(r))
}

#library(mvtnorm) #dmvnorm

calc_likelihood = function(simd, dataid = book[taskid,2], L,
                           betaXj, betatj, gamma1, gamma2j, gamma3j, dscorej, Sigj, lambdaj){
  d = simd[dataID == dataid,]
  
  X = cbind(1, as.matrix(d[, c("X1", "X2")]))
  y = as.matrix(d[, paste0("Y",1:K), with=F]) # Nobs x K
  time = as.numeric(d$age)
  
  mu1 = cbind(X, time) %*% rbind(betaXj, betatj) # Nobs x K
  
  mu_tmp = lapply(1:K, function(k) gamma1[k]/(1 + exp( -gamma2j[k] * (dscorej - gamma3j[k])) ))
  mu = array( unlist( mu_tmp ) , dim = c( Nobs, L, K ))

  yloglik = array(NA, dim = c( Nobs, L, K )) 
  person_loglik = matrix(NA, nrow = N, ncol = L)
  for(ell in 1:L){
    mu[,ell,] = mu[,ell,] + mu1
    yloglik[,ell,] =  dnorm(y, mean = mu[,ell,], 
                            sd = matrix(sqrt(Sigj[ell,]), nrow=1) %x% matrix(1, nrow = Nobs, ncol=1), 
                            log = TRUE) # Nobs x K 
    person_loglik[,ell] = rowSums(matrix(t(yloglik[,ell,]), ncol = Ni*K, byrow = T), dims = 1) + 
      log(lambdaj[ell])
  }
  loglik = sum(log(apply(exp(person_loglik), 1, sum)))

  return(loglik)
}


calc_dic = function(pidx, gamma1, res, L, cidx, simd, taskid){
  dn = dimnames(res)
  betaXch = res[pidx, cidx, grep("betaX",dn$parameters )] # p x K
  betatch = res[pidx, cidx, grep("betat",dn$parameters )] # 1 x K
  gamma2ch = res[pidx, cidx, grep("gamma2",dn$parameters )] # 1 x K
  gamma3ch = res[pidx, cidx, grep("gamma3",dn$parameters )] # 1 x K
  dscorech = res[pidx, cidx, grep("d\\[",dn$parameters )] # n x L
  Sigch =  res[pidx, cidx, grep("Sig",dn$parameters )] # L x K
  if(L == 1){
    lambdach = matrix(rep(1, length(pidx)), ncol = 1)
  } else {
    lambdach =  res[pidx, cidx, grep("lambda",dn$parameters )] # L
  }
  
  loglikvec = lapply(pidx, function(pj) calc_likelihood(simd, dataid = book[taskid,2], L,
                                                        t(matrix(betaXch[pj,], nrow = K)), 
                                                        matrix(betatch[pj,], nrow = 1), gamma1, 
                                                        gamma2ch[pj,], gamma3ch[pj,], 
                                                        matrix(dscorech[pj,], nrow = Nobs), 
                                                        matrix(Sigch[pj,], nrow = L) , lambdach[pj,] ))
  loglikvec = unlist(loglikvec)
  dic = (-2) * mean(loglikvec) + var(loglikvec)/2
  return(dic)
}

pidx = 1:9000

for(modid in 1:3){
  
  if(modid == 1){
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha", "rho_theta", "d", "lp__")
  } else {
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha0", "alpha", "rho_theta", "d", "lp__", "lambda")
  }
  
  taskvec = which(book[,1]==modid)
  
  ### CP
  for(CPpc in c(seq(0.75, 0.95, 0.05), 0.99)){
    
    CP = vector("list", length(allnames))
    names(CP) = allnames
    
    for(taskid in taskvec){
      print(CPpc); print(taskid)
      
      sfn = paste0("/fastscratch/myscratch/yxu2/simu11_",taskid,".RData")
      oldfn = paste0("/users/yxu2/AD/simu11_",taskid,".RData")
      
      if(file.exists(sfn)){
        load(sfn)
      } else if(file.exists(oldfn)){
        load(oldfn)
      } else {
        next
      }
      
      dn = dimnames(res)
      
      doit = 1
      
      if(dim(res)[2] == 1){
        cidx = 1
        sm = apply(res[pidx, 1, grep("Sig",dn$parameters )],2,  mean)
        if(abs(max(sm)) >= 100 && modid <=2) doit = 0
        if(modid == 2){
          sm = mean(res[pidx, 1, grep("alpha0",dn$parameters )])
          if(sm < 0.5) doit = 0
        }
        
      } else {
        maxsig = rep(NA,2)
        malpha0 = rep(NA,2)
        for(s in 1:2){
          sm = apply(res[pidx, s, grep("Sig",dn$parameters )],2,  mean)
          maxsig[s] = abs(max(sm)) 
          if(modid == 2) malpha0[s] = mean(res[pidx, s, grep("alpha0",dn$parameters )])
        }
        
        cidxvec = which(maxsig < 100)
        if(length(cidxvec) == 0) doit = 0
        if(length(cidxvec) == 1){
          cidx = cidxvec
          if(modid == 2 & malpha0[cidx] < 0.5) doit = 0
        }
        if(length(cidxvec) == 2) {
          cidx = which.max(c(mean(res[pidx, 1, grep("lp__",dn$parameters )]), 
                             mean(res[pidx, 2, grep("lp__",dn$parameters )])))
          if(modid == 2){
            aidxvec = which(malpha0 > 0.5)
            if(length(aidxvec) == 1) 
              cidx = aidxvec
            if(length(aidxvec) == 0)
              doit = 0
          }
          
        }
      }
      
      if(doit){
        
        vn = c("gamma2","gamma3","betaX","betat")
        for(j in 1:length(vn)){
          ch = res[pidx, cidx, grep(vn[j],dn$parameters )]
          CP[[j]] = rbind(CP[[j]], 
                          unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], get(vn[j])[l], perc = CPpc))))
        }
        
        # Sig
        ch = res[pidx, cidx, grep("Sig",dn$parameters )]
        if(modid == 2) 
          CP$Sig = rbind(CP$Sig , 
                         unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], sig[l], perc = CPpc))))
        if(modid == 1) 
          CP$Sig = rbind(CP$Sig , 
                         c(unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], sig1[l], perc = CPpc))),
                           unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], sig2[l], perc = CPpc)))) )
        if(modid == 3){
          imat = matrix(1:27, nrow = 3)
          CP$Sig = rbind(CP$Sig , 
                         c(unlist(lapply(1:9,  function(l) CPfun(ch[,imat[1,l]], sig1[l], perc = CPpc))),
                           unlist(lapply(1:9,  function(l) CPfun(ch[,imat[2,l]], sig1[l], perc = CPpc))),
                           unlist(lapply(1:9,  function(l) CPfun(ch[,imat[3,l]], sig1[l], perc = CPpc))),
                           unlist(lapply(1:9,  function(l) CPfun(ch[,imat[1,l]], sig2[l], perc = CPpc))),
                           unlist(lapply(1:9,  function(l) CPfun(ch[,imat[2,l]], sig2[l], perc = CPpc))),
                           unlist(lapply(1:9,  function(l) CPfun(ch[,imat[3,l]], sig2[l], perc = CPpc)))) )
        }
        
        
        # alphat
        ch = res[pidx, cidx, grep("alpha\\[",dn$parameters )]
        CP$alpha = rbind(CP$alpha, CPfun(ch, alphat, perc = CPpc))
        
        # alpha0
        if(modid > 0){
          ch = res[pidx, cidx,  grep("alpha0",dn$parameters )]
          if(modid == 2){
            CP$alpha0 = rbind(CP$alpha0, CPfun(ch, alpha0, perc = CPpc))
          }
          if(modid == 3){
            CP$alpha0 = rbind(CP$alpha0, 
                              c(CPfun(ch[,1], alpha0), CPfun(ch[,2], alpha0, perc = CPpc)))
          }
        }
        
        # rho
        ch = res[pidx, cidx, grep("rho_theta",dn$parameters )]
        if(modid == 2) 
          CP$rho_theta = rbind(CP$rho_theta, 
                               unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], rho[l], perc = CPpc))))
        if(modid == 1) 
          CP$rho_theta = rbind(CP$rho_theta , 
                               c(CPfun(ch, rho[1]),CPfun(ch, rho[2], perc = CPpc)) )
        if(modid == 3){
          CP$rho_theta = rbind(CP$rho_theta , 
                               c(CPfun(ch[,1], rho[1], perc = CPpc),
                                 CPfun(ch[,2], rho[1], perc = CPpc),
                                 CPfun(ch[,3], rho[1], perc = CPpc),
                                 CPfun(ch[,1], rho[2], perc = CPpc),
                                 CPfun(ch[,2], rho[2], perc = CPpc),
                                 CPfun(ch[,3], rho[2], perc = CPpc) ) )
        }
        
        # lambda
        ch = res[pidx, cidx, grep("lambda",dn$parameters )]
        if(modid == 2) 
          CP$lambda = rbind(CP$lambda, 
                            unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], lambda[l], perc = CPpc))))
        
        if(modid == 3){
          CP$lambda = rbind(CP$lambda , 
                            c(CPfun(ch[,1], lambda[1], perc = CPpc),
                              CPfun(ch[,2], lambda[1], perc = CPpc),
                              CPfun(ch[,3], lambda[1], perc = CPpc),
                              CPfun(ch[,1], lambda[2], perc = CPpc),
                              CPfun(ch[,2], lambda[2], perc = CPpc),
                              CPfun(ch[,3], lambda[2], perc = CPpc) ) )
        }
 
      } # doit
    } # taskid
    assign(paste0("CP", CPpc*100), CP)
  } # CPpc
  save(CP75, CP80, CP85, CP90, CP95, CP99, file = paste0("postsummCP", modid,".RData"))
  
} # modid

### postm, BIC, DIC

#library(EDISON) # psrf requires multiple chains
library(coda) # use geweke diagnostic for stationarity of single chain

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

for(modid in 1:3){
  
  if(modid == 1){
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha", "rho_theta", "d", "lp__")
  } else {
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha0", "alpha", "rho_theta", "d", "lp__", "lambda")
  }
  
  taskvec = which(book[,1]==modid)
  
  postm = conv = vector("list", length(allnames))
  names(postm) = names(conv) = allnames
  
  #lp = RMSE = 
  BIC = DIC = c()
  
  for(taskid in taskvec){
    print(taskid)
    
    sfn = paste0("/fastscratch/myscratch/yxu2/simu11_",taskid,".RData")
    oldfn = paste0("/users/yxu2/AD/simu11_",taskid,".RData")
    
    if(file.exists(sfn)){
      load(sfn)
    } else if(file.exists(oldfn)){
      load(oldfn)
    } else {
      next
    }
    
    dn = dimnames(res)
    
    doit = 1
    
    if(dim(res)[2] == 1){
      cidx = 1
      sm = apply(res[pidx, 1, grep("Sig",dn$parameters )],2,  mean)
      if(abs(max(sm)) >= 100 && modid <=2) doit = 0
      if(modid == 2){
        sm = mean(res[pidx, 1, grep("alpha0",dn$parameters )])
        if(sm < 0.5) doit = 0
      }
      
    } else {
      maxsig = rep(NA,2)
      malpha0 = rep(NA,2)
      for(s in 1:2){
        sm = apply(res[pidx, s, grep("Sig",dn$parameters )],2,  mean)
        maxsig[s] = abs(max(sm)) 
        if(modid == 2) malpha0[s] = mean(res[pidx, s, grep("alpha0",dn$parameters )])
      }
      
      cidxvec = which(maxsig < 100)
      if(length(cidxvec) == 0) doit = 0
      if(length(cidxvec) == 1){
        cidx = cidxvec
        if(modid == 2 & malpha0[cidx] < 0.5) doit = 0
      }
      if(length(cidxvec) == 2) {
        cidx = which.max(c(mean(res[pidx, 1, grep("lp__",dn$parameters )]), 
                           mean(res[pidx, 2, grep("lp__",dn$parameters )])))
        if(modid == 2){
          aidxvec = which(malpha0 > 0.5)
          if(length(aidxvec) == 1) 
            cidx = aidxvec
          if(length(aidxvec) == 0)
            doit = 0
        }
        
      }
    }
    
    if(doit){
      
      ### postm
      vn = c("gamma2","gamma3","betaX","betat")
      for(j in 1:length(vn)){
        ch = res[pidx, cidx, grep(vn[j],dn$parameters )]
        postm[[j]] = rbind(postm[[j]], apply(ch, 2, mean) )
        conv[[j]] = rbind(conv[[j]], mygeweke(ch, 0.45, 0.05))
      }
      
      # Sig
      ch = res[pidx, cidx, grep("Sig",dn$parameters )]
      postm$Sig = rbind(postm$Sig , apply(ch, 2, mean) )
      conv$Sig = rbind(conv$Sig, mygeweke(ch, 0.45, 0.05))
      
      # alphat
      ch = res[pidx, cidx, grep("alpha\\[",dn$parameters )]
      postm$alpha = rbind(postm$alpha, mean(ch) )
      conv$alpha = c(conv$alpha, mygeweke(ch, 0.45, 0.05))
      
      # alpha0
      ch = res[pidx, cidx,  grep("alpha0",dn$parameters )]
      if(modid == 2){
        postm$alpha0 = rbind(postm$alpha0, mean(ch) )
        conv$alpha0 = c(conv$alpha0, mygeweke(ch, 0.45, 0.05))
      }
      if(modid == 3){
        postm$alpha0 = rbind(postm$alpha0, apply(ch, 2, mean) )
        conv$alpha0 = rbind(conv$alpha0, mygeweke(ch, 0.45, 0.05) )
      }
      
      # rho
      ch = res[pidx, cidx, grep("rho_theta",dn$parameters )]
      if(modid==1){
        postm$rho_theta = rbind(postm$rho_theta, mean(ch) )
        conv$rho_theta = c(conv$rho_theta, mygeweke(ch, 0.45, 0.05))
      } else {
        postm$rho_theta = rbind(postm$rho_theta, apply(ch, 2, mean) )
        conv$rho_theta = rbind(conv$rho_theta, mygeweke(ch, 0.45, 0.05) )
      }
      
      # lambda
      if(modid > 1){
        ch = res[pidx, cidx, grep("lambda",dn$parameters )]
        postm$lambda = rbind(postm$lambda, apply(ch,2,mean))
        conv$lambda = rbind(conv$lambda, mygeweke(ch, 0.45, 0.05) )
      }
      
      ### calculate BIC
      L = modid
      qfix = (3 + 1 + L + 2) * K  # beta0, betaX1, betaX2, betat, sig, gamma2&3
      qre = 1 + (L-1) + L # alphat, alpha0, rho
      q = (L-1) + qfix + qre # lambda, fix params, re params 
      
      type = "lp__"; idod = 1
      lp = res[, cidx, which(dn$parameters == type)]
      maxidx = order(lp, decreasing=T)[idod]
      ll = lp[maxidx]
      BIC =  c(BIC, -2 * ll+ log(N) * q )
      
      ### calculate DIC
      DIC = c(DIC,  calc_dic(pidx, gamma1, res, L, cidx, simd, taskid))
      
      #lp = c(lp, mean(res[pidx, cidx, grep("lp__",dn$parameters )]))
      rm(res); gc()
    } # doit
  } # taskid
  
  save(postm, BIC, DIC, conv, file = paste0("postsumm", modid,".RData"))
  
} # modid


# new geweke fraction
pidx = 4001:9000
for(modid in 1:3){
  
  if(modid == 1){
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha", "rho_theta", "d", "lp__")
  } else {
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha0", "alpha", "rho_theta", "d", "lp__", "lambda")
  }
  
  taskvec = which(book[,1]==modid)
  
  # RMSE = vector("list", length(allnames))
  # names(RMSE) = allnames
  
  conv1 = vector("list", length(allnames))
  names(conv1) = allnames
  
  for(taskid in taskvec){
    print(taskid)
    
    sfn = paste0("/fastscratch/myscratch/yxu2/simu11_",taskid,".RData")
    oldfn = paste0("/users/yxu2/AD/simu11_",taskid,".RData")
    
    if(file.exists(sfn)){
      load(sfn)
    } else if(file.exists(oldfn)){
      load(oldfn)
    } else {
      next
    }
    
    dn = dimnames(res)
    
    doit = 1
    
    if(dim(res)[2] == 1){
      cidx = 1
      sm = apply(res[pidx, 1, grep("Sig",dn$parameters )],2,  mean)
      if(abs(max(sm)) >= 100 && modid <=2) doit = 0
      if(modid == 2){
        sm = mean(res[pidx, 1, grep("alpha0",dn$parameters )])
        if(sm < 0.5) doit = 0
      }
      
    } else {
      maxsig = rep(NA,2)
      malpha0 = rep(NA,2)
      for(s in 1:2){
        sm = apply(res[pidx, s, grep("Sig",dn$parameters )],2,  mean)
        maxsig[s] = abs(max(sm)) 
        if(modid == 2) malpha0[s] = mean(res[pidx, s, grep("alpha0",dn$parameters )])
      }
      
      cidxvec = which(maxsig < 100)
      if(length(cidxvec) == 0) doit = 0
      if(length(cidxvec) == 1){
        cidx = cidxvec
        if(modid == 2 & malpha0[cidx] < 0.5) doit = 0
      }
      if(length(cidxvec) == 2) {
        cidx = which.max(c(mean(res[pidx, 1, grep("lp__",dn$parameters )]), 
                           mean(res[pidx, 2, grep("lp__",dn$parameters )])))
        if(modid == 2){
          aidxvec = which(malpha0 > 0.5)
          if(length(aidxvec) == 1) 
            cidx = aidxvec
          if(length(aidxvec) == 0)
            doit = 0
        }
        
      }
    }
    
    if(doit){
      
      ### another geweke
      vn = c("gamma2","gamma3","betaX","betat")
      for(j in 1:length(vn)){
        ch = res[pidx, cidx, grep(vn[j],dn$parameters )]
        conv1[[j]] = rbind(conv1[[j]], mygeweke(ch, 0.1, 0.5))
      }
      
      # Sig
      ch = res[pidx, cidx, grep("Sig",dn$parameters )]
      conv1$Sig = rbind(conv1$Sig, mygeweke(ch, 0.1, 0.5))
      
      # alphat
      ch = res[pidx, cidx, grep("alpha\\[",dn$parameters )]
      conv1$alpha = c(conv1$alpha, mygeweke(ch, 0.1, 0.5))
      
      # alpha0
      ch = res[pidx, cidx,  grep("alpha0",dn$parameters )]
      if(modid == 2){
        conv1$alpha0 = c(conv1$alpha0, mygeweke(ch, 0.1, 0.5))
      }
      if(modid == 3){
        conv1$alpha0 = rbind(conv1$alpha0, mygeweke(ch, 0.1, 0.5) )
      }
      
      # rho
      ch = res[pidx, cidx, grep("rho_theta",dn$parameters )]
      if(modid==1){
        conv1$rho_theta = c(conv1$rho_theta, mygeweke(ch, 0.1, 0.5))
      } else {
        conv1$rho_theta = rbind(conv1$rho_theta, mygeweke(ch, 0.1, 0.5) )
      }
      
      # lambda
      if(modid > 1){
        ch = res[pidx, cidx, grep("lambda",dn$parameters )]
        conv1$lambda = rbind(conv1$lambda, mygeweke(ch, 0.1, 0.5) )
      }
      
      #lp = c(lp, mean(res[pidx, cidx, grep("lp__",dn$parameters )]))
      rm(res); gc()
    } # doit
  } # taskid
  
  save( conv1, file = paste0("postsumm_conv1_", modid,".RData"))
  
} # modid

### calculate RRMSE

ufun = function(x){# for root MSE
  rr = sqrt(mean(x^2))
  return(rr)
}
#RMSE$betaX = apply(t((t(postm$betaX) - betaX)/betaX),2, ufun) 
rrmsefun = function(x, truth){
  num = apply((t(t(x) - truth))^2, 2, mean)
  if(is.null(ncol(x))){
    dem = sum(x^2)
  } else {
    dem = apply(x, 2, function(u) sum(u^2))
  }
  return(sqrt(num/dem))
}
source("BMsim11Y_setting.R")
# ## when we treat posterior mean (for each replicate) as the estimate
# for(modid in 1:3){ 
#   if(modid == 1){
#     allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha", "rho_theta", "d", "lp__")
#   } else {
#     allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha0", "alpha", "rho_theta", "d", "lp__", "lambda")
#   }
#   RRMSE = vector("list", length(allnames))
#   names(RRMSE) = allnames
#   load(paste0("postsumm", modid,".RData"))
#   
#   RRMSE$gamma2 = rrmsefun( postm$gamma2, gamma2)
#   RRMSE$gamma3 = rrmsefun( postm$gamma3, gamma3)
#   betaX = c(cbind(beta0, betaX1, betaX2))
#   RRMSE$betaX = rrmsefun( postm$betaX, betaX)
#   RRMSE$betat = rrmsefun( postm$betat, betat)
#   RRMSE$Sig = rrmsefun( postm$Sig, sig)
#   RRMSE$alpha = rrmsefun(postm$alpha, alphat)
#   if(modid==1){
#     RRMSE$rho_theta = c(rrmsefun(postm$rho_theta, rho[1]), rrmsefun(postm$rho_theta, rho[2]) )
#   }
#   if(modid==2){
#     RRMSE$alpha0 = rrmsefun(postm$alpha0, alpha0)
#     RRMSE$rho_theta = rrmsefun(postm$rho_theta, rho)
#     RRMSE$lambda = rrmsefun(postm$lambda, lambda)
#   }
#   if(modid==3){
#     RRMSE$alpha0 = rrmsefun(postm$alpha0, c(alpha0, alpha0))
#     RRMSE$rho_theta = c(rrmsefun(postm$rho_theta, rep(rho[1],modid)), 
#                         rrmsefun(postm$rho_theta, rep(rho[2],modid)) )
#     RRMSE$lambda = c(rrmsefun(postm$lambda, rep(lambda[1],modid)),
#                      rrmsefun(postm$lambda, rep(lambda[2],modid)) )
#   }
#   save(RRMSE, file = paste0("postsumm_RRMSE", modid,".RData"))
# }

## when we treat posterior draw as an estimate, summarize for each replicate,
## and summarize mean (95% CI) across replicates
pidx = 1:9000
for(modid in 1:3){
  
  if(modid == 1){
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha", "rho_theta", "d", "lp__")
  } else {
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha0", "alpha", "rho_theta", "d", "lp__", "lambda")
  }
  
  taskvec = which(book[,1]==modid)
  
  # RMSE = vector("list", length(allnames))
  # names(RMSE) = allnames
  
  RRMSE = vector("list", length(allnames))
  names(RRMSE) = allnames
  
  for(taskid in taskvec){
    print(taskid)
    
    sfn = paste0("/fastscratch/myscratch/yxu2/simu11_",taskid,".RData")
    oldfn = paste0("/users/yxu2/AD/simu11_",taskid,".RData")
    
    if(file.exists(sfn)){
      load(sfn)
    } else if(file.exists(oldfn)){
      load(oldfn)
    } else {
      next
    }
    
    dn = dimnames(res)
    
    doit = 1
    
    if(dim(res)[2] == 1){
      cidx = 1
      sm = apply(res[pidx, 1, grep("Sig",dn$parameters )],2,  mean)
      if(abs(max(sm)) >= 100 && modid <=2) doit = 0
      if(modid == 2){
        sm = mean(res[pidx, 1, grep("alpha0",dn$parameters )])
        if(sm < 0.5) doit = 0
      }
      
    } else {
      maxsig = rep(NA,2)
      malpha0 = rep(NA,2)
      for(s in 1:2){
        sm = apply(res[pidx, s, grep("Sig",dn$parameters )],2,  mean)
        maxsig[s] = abs(max(sm)) 
        if(modid == 2) malpha0[s] = mean(res[pidx, s, grep("alpha0",dn$parameters )])
      }
      
      cidxvec = which(maxsig < 100)
      if(length(cidxvec) == 0) doit = 0
      if(length(cidxvec) == 1){
        cidx = cidxvec
        if(modid == 2 & malpha0[cidx] < 0.5) doit = 0
      }
      if(length(cidxvec) == 2) {
        cidx = which.max(c(mean(res[pidx, 1, grep("lp__",dn$parameters )]), 
                           mean(res[pidx, 2, grep("lp__",dn$parameters )])))
        if(modid == 2){
          aidxvec = which(malpha0 > 0.5)
          if(length(aidxvec) == 1) 
            cidx = aidxvec
          if(length(aidxvec) == 0)
            doit = 0
        }
        
      }
    }
    
    if(doit){
      
      # parameters shared across latent classes
      betaX = c(cbind(beta0, betaX1, betaX2))
      truthlist = list(gamma2, gamma3, betaX, betat)
      vn = c("gamma2","gamma3","betaX","betat")
      for(j in 1:length(vn)){
        ch = res[pidx, cidx, grep(vn[j],dn$parameters )]
        RRMSE[[j]] = rbind(RRMSE[[j]], rrmsefun( ch, truthlist[[j]]))
      }
      
      # Sig
      ch = res[pidx, cidx, grep("Sig",dn$parameters )]
      if(modid == 1){
        RRMSE$Sig = rbind(RRMSE$Sig,  
                          c(rrmsefun( ch, sig1), rrmsefun( ch, sig2)) )
      }
      if(modid == 2){
        RRMSE$Sig = rbind(RRMSE$Sig,  rrmsefun( ch, sig))
      }
      if(modid == 3){
        RRMSE$Sig = rbind(RRMSE$Sig,  
                          c(rrmsefun( ch, rep(sig1, each = modid)), rrmsefun( ch, rep(sig2, each = modid))) )
      }
      
      # alphat
      ch = res[pidx, cidx, grep("alpha\\[",dn$parameters )]
      RRMSE$alpha = c(RRMSE$alpha, rrmsefun( ch, alphat))
      
      # alpha0
      ch = res[pidx, cidx,  grep("alpha0",dn$parameters )]
      if(modid == 2){
        RRMSE$alpha0 = c(RRMSE$alpha0, rrmsefun( ch, alpha0))
      }
      if(modid == 3){
        RRMSE$alpha0 = rbind(RRMSE$alpha0, rrmsefun( ch, c(alpha0, alpha0)) )
      }
      
      # rho
      ch = res[pidx, cidx, grep("rho_theta",dn$parameters )]
      if(modid==1){
        RRMSE$rho_theta = rbind(RRMSE$rho_theta,
                                c(rrmsefun( ch, rho[1]), rrmsefun( ch, rho[2]) ))
      }
      if(modid==2){
        RRMSE$rho_theta = rbind(RRMSE$rho_theta, rrmsefun( ch, rho) )
      }
      if(modid==3){
        RRMSE$rho_theta = rbind(RRMSE$rho_theta,
                                c(rrmsefun( ch, rep(rho[1],modid)), 
                                  rrmsefun( ch, rep(rho[2],modid)) ))
      }
      # lambda
      if(modid > 1){
        ch = res[pidx, cidx, grep("lambda",dn$parameters )]
        if(modid==2){
          RRMSE$lambda = rbind(RRMSE$lambda, rrmsefun(ch, lambda) )
        }
        if(modid==3)
        RRMSE$lambda = rbind(RRMSE$lambda, 
                             c(rrmsefun(ch, rep(lambda[1],modid)),
                               rrmsefun(ch, rep(lambda[2],modid)) ))
      }
      
      #lp = c(lp, mean(res[pidx, cidx, grep("lp__",dn$parameters )]))
      rm(res); gc()
    } # doit
  } # taskid
  
  save( RRMSE, file = paste0("postsumm_RRMSE_", modid,".RData"))
  
} # modid


####################################################################



