# By p9.pdf, we see that L=2 has optimal and suboptimal solutions, with 
# average posteriors  -30792.54 and -30920.29, respectively. The suboptimal
# solution has alpha0 < 0.5, as in p10.pdf.
# Hence, the criteria of including a chain from a simulation replicate are
# 1. parameters converged, e.g. sigma^2 for all outcomes < 100
# 2. if L = 2, exclude suboptimal solution, e.g. alpha0 < 0.5
# 3. if both 1 and 2 satisfied, choose the chain (out of the two chains) with the larger posterior

# For revision, write functions for calculating BIC, DIC, potential scale reduction factor (PSRF), and 
# coverage probabilities at 75%, 80%, 85%, 90%, 95%, and 99%. 
# Also provide summaries for lambda as results for the classification of clusters.

setwd("./AD")

library(abind)
# modid is the number of clusters specified in the model estimation

#source("BMsim_setting2.R")
source("BMsim11Y_setting.R")

#sig1 = sqsig1^2; sig2 = sqsig2^2
sig = c(rbind(sig1, sig2))
betaX = c(cbind(beta0,betaX1,betaX2))


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

pidx = 1:9000

for(modid in 1:3){
  
  if(modid == 1){
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha", "rho_theta", "d", "lp__")
  } else {
    allnames = c("gamma2","gamma3","betaX","betat", "Sig", "alpha0", "alpha", "rho_theta", "d", "lp__", "lambda")
  }
  postm = vector("list", length(allnames))
  names(postm) = allnames
  CP = vector("list", length(allnames))
  names(CP) = names(postm)
  
  taskvec = which(book[,1]==modid)
  lp = c()
  
  for(taskid in taskvec){
    
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
        postm[[j]] = rbind(postm[[j]], apply(ch, 2, mean) )
        CP[[j]] = rbind(CP[[j]], 
                        unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], get(vn[j])[l]))))
      }
      
      # Sig
      
      ch = res[pidx, cidx, grep("Sig",dn$parameters )]
      postm$Sig = rbind(postm$Sig , apply(ch, 2, mean) )
      if(modid == 2) 
        CP$Sig = rbind(CP$Sig , 
                        unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], sig[l]))))
      if(modid == 1) 
        CP$Sig = rbind(CP$Sig , 
                       c(unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], sig1[l]))),
                         unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], sig2[l])))) )
      if(modid == 3){
        imat = matrix(1:27, nrow = 3)
        CP$Sig = rbind(CP$Sig , 
                       c(unlist(lapply(1:9,  function(l) CPfun(ch[,imat[1,l]], sig1[l]))),
                         unlist(lapply(1:9,  function(l) CPfun(ch[,imat[2,l]], sig1[l]))),
                         unlist(lapply(1:9,  function(l) CPfun(ch[,imat[3,l]], sig1[l]))),
                         unlist(lapply(1:9,  function(l) CPfun(ch[,imat[1,l]], sig2[l]))),
                         unlist(lapply(1:9,  function(l) CPfun(ch[,imat[2,l]], sig2[l]))),
                         unlist(lapply(1:9,  function(l) CPfun(ch[,imat[3,l]], sig2[l])))) )
      }
        
            
      # alphat
      
      ch = res[pidx, cidx, grep("alpha\\[",dn$parameters )]
      idx = which(names(postm)=="alpha")
      postm$alpha = rbind(postm$alpha, mean(ch) )
      CP$alpha = rbind(CP$alpha, CPfun(ch, alphat))
      
      # alpha0
      
      if(modid > 0){
        ch = res[pidx, cidx,  grep("alpha0",dn$parameters )]
        if(modid == 2){
          postm$alpha0 = rbind(postm$alpha0, mean(ch) )
          CP$alpha0 = rbind(CP$alpha0, CPfun(ch, alpha0))
        }
        if(modid == 3){
          postm$alpha0 = rbind(postm$alpha0, apply(ch, 2, mean) )
          CP$alpha0 = rbind(CP$alpha0, 
                            c(CPfun(ch[,1], alpha0), CPfun(ch[,2], alpha0)))
        }
      }
      
      # rho
      
      ch = res[pidx, cidx, grep("rho_theta",dn$parameters )]
      idx = which(names(postm)=="rho_theta")
      if(modid==1){
        postm$rho_theta = rbind(postm$rho_theta, mean(ch) )
      } else {
        postm$rho_theta = rbind(postm$rho_theta, apply(ch, 2, mean) )
      }
      if(modid == 2) 
        CP$rho_theta = rbind(CP$rho_theta, 
                          unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], rho[l]))))
      if(modid == 1) 
        CP$rho_theta = rbind(CP$rho_theta , 
                       c(CPfun(ch, rho[1]),CPfun(ch, rho[2])) )
      if(modid == 3){
        CP$rho_theta = rbind(CP$rho_theta , 
                       c(CPfun(ch[,1], rho[1]),
                         CPfun(ch[,2], rho[1]),
                         CPfun(ch[,3], rho[1]),
                         CPfun(ch[,1], rho[2]),
                         CPfun(ch[,2], rho[2]),
                         CPfun(ch[,3], rho[2]) ) )
      }
  
      ch = res[pidx, cidx, grep("lambda",dn$parameters )]
      postm$lambda = rbind(postm$lambda, apply(ch,2,mean))
      if(modid == 2) 
        CP$lambda = rbind(CP$lambda, 
                       unlist(lapply(1:ncol(ch),  function(l) CPfun(ch[,l], lambda[l]))))
     
      if(modid == 3){
        CP$lambda = rbind(CP$lambda , 
                       c(CPfun(ch[,1], lambda[1]),
                         CPfun(ch[,2], lambda[1]),
                         CPfun(ch[,3], lambda[1]),
                         CPfun(ch[,1], lambda[2]),
                         CPfun(ch[,2], lambda[2]),
                         CPfun(ch[,3], lambda[2]) ) )
      }
      
      lp = c(lp, mean(res[pidx, cidx, grep("lp__",dn$parameters )]))
    } # doit
  } # taskid
  save(postm, CP, lp, file = paste0("postsumm", modid,".RData"))
  
  
} # modid


####################################################################
