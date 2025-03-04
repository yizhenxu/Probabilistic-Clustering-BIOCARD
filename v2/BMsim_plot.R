# source("BMsim_table_final.R")

library(ggplot2)
library(ggthemes) 
library(latex2exp)
library("gridExtra")

source("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/code/BMsim11Y_setting.R")
sig = c(rbind(sig1, sig2))
betaX = c(cbind(beta0,betaX1,betaX2))

rdnum = 3

summ = function(x){
  x = x[!is.na(x)]
  sx = sort(x)
  n = length(x)
  r = c(round(mean(x),rdnum), round(sx[ceiling(n*0.025)],rdnum), round(sx[floor(n*0.975)],rdnum))
  return(r)
}


g_legend = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}


summtab = vector("list",3)
lpall = vector("list",3)
for(modid in 1:3){
  
  load(paste0("/Users/yizhenxu/postsumm", modid,".RData"))
  ttype = c("Bias","BiasL","BiasU", "CP")
  alltab = list(Ypar = matrix(NA, ncol = length(ttype)*K, nrow = 6),
                sigpar = matrix(NA,  ncol = length(ttype)*K, nrow = ifelse(modid==2,2,modid*2)),
                rhopar = matrix(NA,  ncol = length(ttype), nrow = ifelse(modid==2,2,modid*2)),
                lambdapar = matrix(NA,  ncol = length(ttype), nrow = ifelse(modid==2,2,modid*2)),
                poppar = matrix(NA, ncol =length(ttype) , nrow = (modid-1)+1))
  rownames(alltab$Ypar) = c("gamma2","gamma3",paste0("betaX",0:2),"betat" )
  colnames(alltab$Ypar) = paste0(ttype ,rep(1:K,each=length(ttype)))
  if(modid == 2){
    rownames(alltab$sigpar) = c("Sig1", "Sig2")
    rownames(alltab$rhopar) = c("rho1", "rho2")
    rownames(alltab$lambdapar) = c("lambda1", "lambda2")
  } else {
    rownames(alltab$sigpar) = paste0(rep(c("Sig1", "Sig2"), each = modid), rep(paste0("sighat",1:modid),2))
    rownames(alltab$rhopar) = paste0(rep(c("rho1", "rho2"), each = modid), rep(paste0("rhohat",1:modid),2))
    rownames(alltab$lambdapar) = paste0(rep(c("lambda1", "lambda2"), each = modid), rep(paste0("lambdahat",1:modid),2))
  } 
  
  colnames(alltab$sigpar) =  paste0(ttype ,rep(1:K,each=length(ttype)))
  colnames(alltab$rhopar) = colnames(alltab$lambdapar) =  ttype
  if(modid == 1){
    rownames(alltab$poppar) = c("alphat")
  } else {
    rownames(alltab$poppar) = c(paste0("alpha0_", 2:modid), "alphat")
  }
  colnames(alltab$poppar) =  ttype
  
  rn = rownames(alltab$Ypar)
  
  cpvec = 4*(0:(K-1))+4
  
  # gamma2 bias
  for(k in 1:K){
    bvec = 4*(k-1)+1:3
    alltab$Ypar[1, bvec] = summ(abs(postm$gamma2[,k]-gamma2[k]))
  }
  # gamma2 CP
  alltab$Ypar[1, cpvec] = round(apply(CP$gamma2,2,mean),rdnum)
  
  # gamma3 bias
  for(k in 1:K){
    bvec = 4*(k-1)+1:3
    alltab$Ypar[2, bvec] = summ(abs(postm$gamma3[,k]-gamma3[k]))
  }
  # gamma3 CP
  alltab$Ypar[2, cpvec] = round(apply(CP$gamma3,2,mean),rdnum)
  
  # betaX bias
  ct = 1
  for(j in 3:5){
    for(k in 1:K){
      bvec = 4*(k-1)+1:3
      alltab$Ypar[j, bvec] = summ(abs(postm$betaX[,ct]-betaX[ct]))
      ct = ct + 1
    }
  }
  # betaX CP
  alltab$Ypar[3:5, cpvec] = t(matrix(round(apply(CP$betaX,2,mean),rdnum),3,3))
  
  # betat bias
  for(k in 1:K){
    bvec = 4*(k-1)+1:3
    alltab$Ypar[6, bvec] = summ(abs(postm$betat[,k]-betat[k]))
  }
  # betat CP
  alltab$Ypar[6, cpvec] = round(apply(CP$betat,2,mean),rdnum)
  
  
  if(modid == 1){
    # sig bias
    for(k in 1:K){
      bvec = 4*(k-1)+1:3
      alltab$sigpar[1, bvec] = summ(abs(postm$Sig[,k]-sig1[k]))
      alltab$sigpar[2, bvec] = summ(abs(postm$Sig[,k]-sig2[k]))
    }
    # sig CP
    alltab$sigpar[, cpvec] = matrix(round(apply(CP$Sig,2,mean),rdnum),byrow=T, nrow=2)
  }
  
  if(modid == 2){
    # sig bias
    ct = 1
    for(k in 1:K){
      for(j in 1:2){
        bvec = 4*(k-1)+1:3
        alltab$sigpar[j, bvec] = summ(abs(postm$Sig[,ct]-sig[ct]))
        ct = ct + 1
      }
    }
    # sig CP
    alltab$sigpar[, cpvec] = matrix(round(apply(CP$Sig,2,mean),rdnum),byrow=T, ncol=K)
  }
  
  if(modid == 3){ 
    imat2 = matrix(1:18, nrow  =  2)
    imat3 = matrix(1:27, nrow  =  3)
    # sig bias
    for(k in 1:K){
      for(j in 1:3){
        bvec = 4*(k-1)+1:3
        alltab$sigpar[j, bvec] = summ(abs(postm$Sig[,imat3[j,k]]-sig[imat2[1,k]]))
        alltab$sigpar[j+3, bvec] = summ(abs(postm$Sig[,imat3[j,k]]-sig[imat2[2,k]]))
      }
    }
    # sig CP
    alltab$sigpar[, cpvec] = matrix(round(apply(CP$Sig,2,mean),rdnum),byrow=T, ncol = K)
  }
  
  if(modid == 2){
    # alpha0 bias
    alltab$poppar[1, 1:3] = round(summ(abs(postm$alpha0-alpha0)),rdnum)
    # alpha0 CP
    alltab$poppar[1, 4] = round(mean(CP$alpha0),rdnum)
  }
  
  if(modid == 3){
    # alpha0 bias
    alltab$poppar[1, 1:3] = round(summ(abs(postm$alpha0[,1]-alpha0)),rdnum)
    alltab$poppar[2, 1:3] = round(summ(abs(postm$alpha0[,2]-alpha0)),rdnum)
    # alpha0 CP
    alltab$poppar[1, 4] = round(mean(CP$alpha0[,1]),rdnum)
    alltab$poppar[2, 4] = round(mean(CP$alpha0[,2]),rdnum)
  }
  
  # alphat bias
  alltab$poppar[modid, 1:3] = round(summ(abs(postm$alphat-alphat)),rdnum)
  # alphat CP
  alltab$poppar[modid, 4] = round(mean(CP$alphat),rdnum)
  
  if(modid == 1){
    # rho bias
    alltab$rhopar[1, 1:3] = round(summ(abs(postm$rho-rho[1])),rdnum)
    alltab$rhopar[2, 1:3] = round(summ(abs(postm$rho-rho[2])),rdnum)
    alltab$rhopar[, 4] = round(apply(CP$rho,2,mean),rdnum)
  }
  if(modid == 2 ){
    # rho bias
    rl = lapply(1:modid, function(j) round(summ(abs(postm$rho[,j]-rho[j])),rdnum))
    for(j in 1:modid){
      alltab$rhopar[j, 1:3] = rl[[j]]
    }
    # rho CP
    alltab$rhopar[1:modid, 4] = round(apply(CP$rho,2,mean),rdnum)
  }
  
  if(modid == 3){
    # rho bias
    for(j in 1:3){
      alltab$rhopar[j, 1:3] = round(summ(abs(postm$rho[,j]-rho[1])),rdnum)
      alltab$rhopar[j+3, 1:3] = round(summ(abs(postm$rho[,j]-rho[2])),rdnum)
    }
    
    alltab$rhopar[, 4] = round(apply(CP$rho,2,mean),rdnum)
  }
  
  
  if(modid == 2 ){
    # lambda bias
    rl = lapply(1:modid, function(j) round(summ(abs(postm$lambda[,j]-lambda[j])),rdnum))
    for(j in 1:modid){
      alltab$lambdapar[j, 1:3] = rl[[j]]
    }
    # lambda CP
    alltab$lambdapar[1:modid, 4] = round(apply(CP$lambda,2,mean),rdnum)
  }
  
  if(modid == 3){
    # lambda bias
    for(j in 1:3){
      alltab$lambdapar[j, 1:3] = round(summ(abs(postm$lambda[,j]-lambda[1])),rdnum)
      alltab$lambdapar[j+3, 1:3] = round(summ(abs(postm$lambda[,j]-lambda[2])),rdnum)
    }
    
    alltab$lambdapar[, 4] = round(apply(CP$lambda,2,mean),rdnum)
  }
  summtab[[modid]] = alltab
  lpall[[modid]] = lp
}


if(0){ # posteriors, gamma3 and alpha0 histograms before excluding suboptimal solution of L=2
  
  ### plot of gamma3 
  dat = c()
  for(modid in 1:3){
    load(paste0("postsumm", modid,".RData"))
    
    for(j in 1:3){
      dat = rbind(dat, cbind(c(postm$gamma3), rep(1:K, each = nrow(postm$gamma3)), modid) ) 
    }
  }
  dat = as.data.frame(dat); colnames(dat) = c("gamma3", "Yindex", "modid")
  
  p = vector("list", K)
  for(j in 1:K){
    pd = dat[dat$Yindex == j,]
    p[[j]] = ggplot( data = pd,aes(x=gamma3, fill=as.factor(modid))) +
      geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
      scale_fill_manual(name  ="Number of Clusters", values=c("#69b3a2", "#404080","#E69F00")) +
      labs(fill="") + xlab(TeX("$\\gamma_3")) +
      theme(legend.position = "bottom")
  }
  
  g_legend = function(a.gplot){
    tmp = ggplot_gtable(ggplot_build(a.gplot))
    leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend = tmp$grobs[[leg]]
    return(legend)
  }
  
  lg = g_legend(p[[1]])
  pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p9.pdf",
      width = 8, height=7)
  grid.arrange(arrangeGrob(p[[1]] + theme(legend.position = "none"), 
                           p[[2]] + theme(legend.position = "none"), 
                           p[[3]] + theme(legend.position = "none"), 
                           p[[4]] + theme(legend.position = "none"), 
                           p[[5]] + theme(legend.position = "none"), 
                           p[[6]] + theme(legend.position = "none"),
                           p[[7]] + theme(legend.position = "none"), 
                           p[[8]] + theme(legend.position = "none"), 
                           p[[9]] + theme(legend.position = "none"), ncol=3),
               lg, nrow = 2, heights = c(10,1))
  dev.off()
  
  
  ### plot of alpha0 
  dat = c()
  for(modid in 2:3){
    load(paste0("postsumm", modid,".RData"))
    
    for(j in 2:3){
      dat = rbind(dat, cbind(c(postm$alpha0), rep(1:(modid-1), each = nrow(postm$alpha0)), modid) ) 
    }
  }
  dat = as.data.frame(dat); colnames(dat) = c("alpha0", "aindex", "modid")
  dat$type = dat$aindex + dat$modid
  
  pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p10.pdf",
      width = 8, height=7)
  ggplot( data = dat,aes(x=alpha0, fill = as.factor(type))) +
    geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=20) +
    scale_fill_manual(name  ="Type", values=c("#69b3a2", "#404080","#E69F00"), labels = c("L=2","L=3 alpha0(2)", "L=3 alpha0(3)")) +
    labs(fill="") + xlab(TeX("$\\alpha_0")) +
    theme(legend.position = "bottom")
  
  dev.off()
  
  ### compare average lp__ of the optimal and subopt solutions under L=2
  load(paste0("postsumm", 2,".RData"))
  idx =which(postm$gamma3[,4]>0.5)
  mean(lpall[[2]][idx]) # -30792.5
  
  idx =which(postm$gamma3[,4]<=0.5)
  mean(lpall[[2]][idx]) # -30922.9
  
} # if(0)

### plot of betas and gammas

vnvec = c("gamma2", "gamma3", "betaX0", "betaX1", "betaX2", "betat" )

dat = c()

for(modid in 1:3){
  add = c()
  infotab = summtab[[modid]]$Ypar
  S = nrow(infotab)
  
  for(j in 1:S){
    vn = vnvec[j]
    x = 1:9 - (modid)/4
    add = rbind(add, cbind(x, matrix(infotab[j,], byrow = T, ncol = 4), modid) )
  }
  
  add = cbind(1:9, add); add = as.data.frame(add)
  colnames(add) = c("Yindex","x", "mean", "lower", "upper", "CP", "modid")
  add$var = rep(vnvec, each = 9)
  add$varid = rep(1:S, each = 9)
  dat = rbind(dat, add)
}

#create forest plot
p = vector("list", length(vnvec))
titlevec = c("gamma_{2}","gamma_{3}","beta_{0}","beta_{X_1}","beta_{X_2}","beta_{t}")
for(j in 1:length(vnvec)){
  pd = dat[dat$var == vnvec[j],]
  pd$shape = pd$modid + 14
  pd$col = pd$modid
  ptitle = TeX(paste0("Bias in $\\",titlevec[j]))
  p[[j]] = ggplot(data=pd, aes(y=x, x=mean, xmin=lower, xmax=upper)) +
    geom_point(data=pd, aes(color = as.factor(col)), size = 0.9)+
    #geom_point(shape = pd$modid + 14) + 
    geom_errorbarh(height=0.05, aes(color = as.factor(col))) +
    #scale_x_discrete(breaks=1:9, labels=1:9)
    scale_y_continuous(breaks=1:nrow(pd), labels=pd$Yindex) +
    labs(x = ptitle, y = 'Outcome Index') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    theme_hc()  + scale_color_discrete(name  ="Number of Clusters",
                                           labels = c(1,2,3))+
    theme(legend.position = "bottom")
    
}

g_legend = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

lg = g_legend(p[[1]])

pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p4.pdf", 
    width = 6, height=8)
grid.arrange(arrangeGrob(p[[1]] + theme(legend.position = "none"), 
                         p[[2]] + theme(legend.position = "none"), 
                         p[[3]] + theme(legend.position = "none"), 
                         p[[4]] + theme(legend.position = "none"), 
                         p[[5]] + theme(legend.position = "none"), 
                         p[[6]] + theme(legend.position = "none"), ncol=2),
             lg, nrow = 2, heights = c(10,1))
dev.off()

### plot of sigma

dat = c()

# 1st col truth 2nd col estimate
idxl = list(matrix(c(1,2,1,1), ncol=2),
            matrix(c(1,2,1,2), ncol=2),
            matrix(c(rep(1:2,each = 3), rep(1:3,2)), ncol=2))
for(modid in 1:3){
  add = c()
  infotab = summtab[[modid]]$sigpar
  S = nrow(infotab)
  
  for(j in 1:S){
    add = rbind(add, cbind(matrix(infotab[j,], byrow = T, ncol = 4), modid, 
                           idxl[[modid]][j,1],idxl[[modid]][j,2]) )
  }
  
  add = cbind(1:9, add); add = as.data.frame(add)
  colnames(add) = c("Yindex", "mean", "lower", "upper", "CP", "modid", "true", "est")
  dat = rbind(dat, add)
}
library(data.table)
dat = as.data.table(dat)
dat$tp = 1; dat[(modid == 3 & est == 1 & true == 2)|(modid == 3 & est %in% c(2,3)  & true == 1), tp := 0.3] 
library(ggplot2)

#create forest plot
p = vector("list", 2)
titlevec = c("sigma^2_1","sigma^2_2")
for(j in 1:2){
  pd = dat[dat$true == j,]
  num = nrow(unique(cbind(pd$modid,pd$est)))
  pd$x = pd$Yindex - rep((1:(num))/(num+1), each = K)
  pd$col = pd$modid
  pd$shape = pd$est + 14
  ptitle = TeX(paste0("Bias in $\\",titlevec[j]))
  p[[j]] = ggplot(data=pd, aes(y=x, x=mean, xmin=lower, xmax=upper)) +
    geom_point(data=pd, aes(shape = as.factor(shape), color = as.factor(col),alpha=tp), size = 1.7)+
    geom_errorbarh(height=0.05, aes(color = as.factor(col)), alpha=pd$tp )+
    scale_y_continuous(breaks=1:nrow(pd), labels=pd$Yindex) +
    labs(x = ptitle, y = 'Outcome Index') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    theme_hc() + scale_shape_discrete(name  ="Cluster Index",
                                           labels = c(1,2,3))+
     scale_color_discrete(name  ="Number of Clusters",
                           labels = c(1,2,3))+
    theme(legend.position = "bottom")
}

g_legend = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

ptmp = ggplot(data=pd, aes(y=x, x=mean, xmin=lower, xmax=upper)) +
  geom_point(data=pd, aes(shape = as.factor(shape), color = as.factor(col)), size = 1.7)+
  geom_errorbarh(height=0.05, aes(color = as.factor(col)) )+
  scale_y_continuous(breaks=1:nrow(pd), labels=pd$Yindex) +
  labs(x = ptitle, y = 'Outcome Index') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_hc() + scale_shape_discrete(name  ="Cluster Index",
                                    labels = c(1,2,3))+
  scale_color_discrete(name  ="Number of Clusters",
                       labels = c(1,2,3))+
  theme(legend.position = "bottom")
lg = g_legend(ptmp)

pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p5.pdf", 
    width = 6, height=5.8)
grid.arrange(arrangeGrob(p[[1]] + theme(legend.position = "none"), 
                         p[[2]] + theme(legend.position = "none"),  ncol=2),
             lg, nrow = 2, heights = c(5.6,1))
dev.off()


### plot of rho and alpha
# 1st col truth 2nd col estimate
idxl = list(matrix(c(1,2,1,1), ncol=2),
            matrix(c(1,2,1,2), ncol=2),
            matrix(c(rep(1:2,each = 3), rep(1:3,2)), ncol=2))
popidxl = list(matrix(c(2,2,NA,NA),ncol=2,byrow=T),
               matrix(c(2,2,2,3,NA,NA),ncol=2,byrow=T))
dat = c()
for(modid in 1:3){
  add = c()
  infotab = summtab[[modid]]$rhopar
  add = rbind(add, cbind(infotab, modid, idxl[[modid]],idxl[[modid]][,2]+14, "rho"))
  if(modid==1){
    infotab = summtab[[1]]$poppar
    add = rbind(add, cbind(infotab, modid, 1,1,1+14, "alphat"))
  }  

  if(modid>1){
    infotab = summtab[[modid]]$poppar
    add = rbind(add, cbind(infotab, modid, popidxl[[modid-1]],popidxl[[modid-1]][,2]+14, c(rep("alpha0",modid-1),"alphat")))
  }
  add = as.data.frame(add)
  dat = rbind(dat, add)
}
colnames(dat) = c( "mean", "lower", "upper", "CP", "modid", "true", "est", "shape", "var")
for(j in 1:4) dat[,j] = as.numeric(as.character(dat[,j]))
library(data.table)
dat = as.data.table(dat)
dat$tp = 1; dat[(var == "rho" & modid == 3 & est == 1 & true == 2)|(var == "rho" & modid == 3 & est %in% c(2,3)  & true == 1), tp := 0.3] 

#create forest plot
p = vector("list", 4)
titlevec = c("rho_1","rho_2","alpha_0","alpha_t")
for(j in 1:2){
  pd = dat[dat$true == j & dat$var == "rho",]
  num = nrow(unique(cbind(pd$modid,pd$est)))
  pd$x =c(0.95,0.85,0.78,0.755,0.73)
  pd$col = pd$modid
  
  ptitle = TeX(paste0("Bias in $\\",titlevec[j]))
  p[[j]] = ggplot(data=pd, aes(y=x, x=mean, xmin=lower, xmax=upper)) +
    geom_point(data=pd, aes(shape = as.factor(shape), color = col, alpha=tp), size = 1.5)+
    #geom_point(shape = pd$modid + 14) + 
    geom_errorbarh(height=0.02, aes(color = col) , alpha=pd$tp)+
    #scale_x_discrete(breaks=1:9, labels=1:9)
    #scale_y_continuous(breaks=1:nrow(pd), labels=pd$mean) +
    labs(x = ptitle, y = '') +
    geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
    theme_hc() + scale_shape_discrete(name  ="Cluster Index",
                                           labels = c(1,2,3))+
    scale_color_manual(name  ="Number of Clusters",
                       labels = c(1,2,3),values = c("1" = "#F8766D","2" = "#00BA38", "3" = "#619CFF"))+
    theme(legend.position = "bottom",
          axis.ticks.y=element_blank(),axis.text.y=element_blank())+ylim(c(0.72,1))
}
j=3
pd = dat[dat$var == "alpha0",]
pd$x =  c( 0.85 , 0.77, 0.73)
pd$col = pd$modid
ptitle = TeX(paste0("Bias in $\\",titlevec[j]))
p[[j]] = ggplot(data=pd, aes(y=x, x=mean, xmin=lower, xmax=upper)) +
  geom_point(data=pd, aes(shape = as.factor(shape), color = col), size = 1.2)+
  #geom_point(shape = pd$modid + 14) + 
  geom_errorbarh(height=0.02, aes(color = col))+
  #scale_x_discrete(breaks=1:9, labels=1:9)
  #scale_y_continuous(breaks=1:nrow(pd), labels=pd$Yindex) +
  labs(x = ptitle, y = '') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_hc() + scale_shape_manual(name  ="Cluster Index",
                                         labels = c(2,3), values = c("16"=17,"17"=15))+
  scale_color_manual(name  ="Number of Clusters",
                       labels = c(2,3),values = c("2" = "#00BA38", "3" = "#619CFF"))+
  theme(legend.position = "bottom",
        axis.ticks.y=element_blank(),axis.text.y=element_blank())+ylim(c(0.72,1))

j=4
pd = dat[dat$var == "alphat",]
pd$x =  c( 0.95, 0.85, 0.75)
pd$col = pd$modid
ptitle = TeX(paste0("Bias in $\\",titlevec[j]))
p[[j]] = ggplot(data=pd, aes(y=x, x=mean, xmin=lower, xmax=upper)) +
  geom_point(data=pd, aes(color = col), size = 1.2)+
  #geom_point(shape = pd$modid + 14) + 
  geom_errorbarh(height=0.02, aes(color = col))+
  #scale_x_discrete(breaks=1:9, labels=1:9)
  #scale_y_continuous(breaks=1:nrow(pd), labels=pd$Yindex) +
  labs(x = ptitle, y = '') +
  geom_vline(xintercept=0, color='black', linetype='dashed', alpha=.5) +
  theme_hc()  + 
  scale_color_manual(name  ="Number of Clusters",
                     labels = c(1,2,3),values = c("1" = "#F8766D","2" = "#00BA38", "3" = "#619CFF"))+
  theme(legend.position = "bottom",
        axis.ticks.y=element_blank(),axis.text.y=element_blank()) +ylim(c(0.72,1))

pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p6.pdf", 
    width = 6, height=5.5)
grid.arrange(arrangeGrob(p[[1]] + theme(legend.position = "none"), 
                         p[[2]] + theme(legend.position = "none"),
                         p[[3]] + theme(legend.position = "none"),
                         p[[4]] + theme(legend.position = "none"),ncol=2),
             lg, nrow = 2, heights = c(6.6,1))
dev.off()

### plot of posteriors
dat = c()
for(j in 1:3){
  dat = rbind(dat, cbind(lpall[[j]], j) )
}

dat = as.data.frame(dat); colnames(dat) = c("value", "type")
dat$type = as.factor(dat$type)

pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p8.pdf", 
    width = 6, height=3)
ggplot( data = dat,aes(x=value, fill=type)) +
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity',bins=15) +
  scale_fill_manual(values=c("#69b3a2", "#404080","#E69F00")) +theme_hc()+
  labs(fill="") + xlab("Posterior")
dev.off()

#################################################################
### plot of CP
dat = c()

vnvec = c("gamma2", "gamma3", "betaX0", "betaX1", "betaX2", "betat" )

idxl = list(cbind(1,1),
            rbind(c(2,2),c(1,1)),
            rbind(c(2,2),c(2,3),c(1,1)))
  
dat = c()

for(modid in 1:3){
  add = c()
  infotab = summtab[[modid]]$Ypar
  add = cbind(c(infotab[,4*(1:K)]), rownames(infotab), rep(1:K, each = nrow(infotab)), modid, 1,1)
  dat = rbind(dat, add)
  
  num = ifelse(modid <= 2, 1, 3)
  infotab = summtab[[modid]]$sigpar
  add = cbind(c(infotab[,4*(1:K)]), rownames(infotab), rep(1:K, each = nrow(infotab)), modid, rep(1:2,each=num), rep(1:num, 2))
  dat = rbind(dat, add)
  
  infotab = summtab[[modid]]$rhopar
  add = cbind(infotab[,4], rownames(infotab), 1, modid, rep(1:2,each=num), rep(1:num, 2))
  dat = rbind(dat, add)
  
  infotab = summtab[[modid]]$poppar
  add = cbind(infotab[,4], rownames(infotab), 1, modid, idxl[[modid]][,1], idxl[[modid]][,2])
  dat = rbind(dat, add)
}

colnames(dat) = c( "CP", "var", "Yindex", "modid", "truth", "est")
dat = as.data.frame(dat)
dat[dat$modid==2,]$est = dat[dat$modid==2,]$truth

dat[,1] = as.numeric(as.character(dat[,1]))
dat$est = as.numeric(as.character(dat$est))

### CP of betas gammas

p = vector("list", 6)
vnvec = c("gamma2", "gamma3", "betaX0", "betaX1", "betaX2", "betat" )
titlevec = c("gamma_{2}","gamma_{3}","beta_{0}","beta_{X_1}","beta_{X_2}","beta_{t}")
for(j in 1:length(vnvec)){
  pd = dat[dat$var == vnvec[j],]
  ptitle = TeX(paste0("Coverage of $\\",titlevec[j]))
  p[[j]] = ggplot(data=pd, aes(x=CP, y=Yindex)) + 
    geom_point(data=pd, aes( color = modid), size = 1.8)+
    labs(x = ptitle, y = 'Outcome Index') + xlim(c(0,1))+
    theme(legend.position = "bottom")+theme_hc()+
    scale_color_manual(name  ="Number of Clusters",
                       labels = c(1,2,3),values = c("1" = "#F8766D","2" = "#00BA38", "3" = "#619CFF"))
}


g_legend = function(a.gplot){
  tmp = ggplot_gtable(ggplot_build(a.gplot))
  leg = which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend = tmp$grobs[[leg]]
  return(legend)
}

lg = g_legend(p[[6]])

pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p11.pdf", 
    width = 6, height=7)
grid.arrange(arrangeGrob(p[[1]] + theme(legend.position = "none"), 
                         p[[2]] + theme(legend.position = "none"), 
                         p[[3]] + theme(legend.position = "none"), 
                         p[[4]] + theme(legend.position = "none"), 
                         p[[5]] + theme(legend.position = "none"), 
                         p[[6]] + theme(legend.position = "none"),ncol=2),
             lg, nrow = 2, heights = c(10,1))
dev.off()

### CP of sigma rho and alpha

p  = vector("list",3)
titlevec = c("sigma^2_1","sigma^2_2")
for(j in 1:2){
  pd = dat[which(substr(dat$var, 1, 4)==paste0("Sig",j)),]
  pd$shape = pd$est + 14
  ptitle = TeX(paste0("Coverage of $\\",titlevec[j]))
  pd = as.data.table(pd)
  pd$tp = 1; pd[(modid == 3 & truth == 1 & est > 1) |(modid == 3 & truth == 2 & est == 1), tp := 0.3]
  pd$CP[pd$CP==0 & pd$tp!=1] = pd$CP[pd$CP==0 & pd$tp!=1] - 0.005
  p[[j]] =  ggplot(data=pd, aes(x=CP, y=Yindex)) + 
    geom_point(data=pd, aes(shape = as.factor(shape), color = modid), size = 2, alpha = pd$tp)+
    labs(x = ptitle, y = 'Outcome Index') + xlim(c(-0.03,1))+
    theme(legend.position = "bottom")+theme_hc()+
    scale_color_manual(name  ="Number of Clusters",
                       labels = c(1,2,3),values = c("1" = "#F8766D","2" = "#00BA38", "3" = "#619CFF"))+
    scale_shape_manual(name  ="Cluster Index",
                       labels = c(1,2,3),values = c("15"= 16,"16"= 17,"17"= 15))
  
}

pd = dat[which(substr(dat$var, 1, 3) %in% c("rho","alp")),]
pd$shape = pd$est + 14
pd$type = substr(pd$var, 1, 3); pd$type[pd$type == "alp"] = "alpha"
pd$suf = ""; pd$suf[pd$type == "alpha"] =  substr(pd$var[pd$type == "alpha"] , 6, 6)
pd$vname = paste0("$\\", pd$type, "_", pd$suf)
pd$vname[pd$type == "rho"] = paste0("$\\", pd$type[pd$type == "rho"],"_", pd$truth[pd$type == "rho"])
pd = as.data.table(pd)
pd$tp = 1; pd[(type == "rho" & modid == 3 & truth == 1 & est > 1) |(type == "rho" & modid == 3 & truth == 2 & est == 1), tp := 0.3]
pd$CP[pd$CP==0 & pd$tp!=1] = pd$CP[pd$CP==0 & pd$tp!=1] - 0.005

p[[3]] = ggplot(data=pd, aes(x=CP, y=vname)) + 
  geom_point(data=pd, aes(shape = as.factor(shape), color = modid), size = 2,alpha = pd$tp)+
  labs(x = "Coverage", y = 'Parameters') + xlim(c(-0.03,1))+
  scale_y_discrete(labels = TeX) +
  theme(legend.position = "bottom")+
  scale_color_manual(name  ="Number of Clusters",
                     labels = c(1,2,3),values = c("1" = "#F8766D","2" = "#00BA38", "3" = "#619CFF"))+
  scale_shape_manual(name  ="Cluster Index",
                     labels = c(1,2,3),values = c("15"= 16,"16"= 17,"17"= 15))+
theme(legend.position = "bottom",
      axis.ticks.y=element_blank()) + theme_hc()

lg = g_legend(p[[1]])
pdf("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/draft/p12.pdf", 
    width = 6, height=4.8)
grid.arrange(arrangeGrob(p[[1]] + theme(legend.position = "none"), 
                         p[[2]] + theme(legend.position = "none"), 
                         p[[3]] + theme(legend.position = "none"), nrow=2),
             lg, nrow = 2, heights = c(8,1))
dev.off()