
# define allidx and eachct

eachct = d[, .N, by=Study.ID]
idlist = eachct[N > 1, Study.ID]
#idlist = unique(d$Study.ID)
allidx = which(d[,Study.ID] %in% idlist)
eachct =  d[Study.ID %in% idlist,.N, by=Study.ID]$N

ptm_mclapply <- proc.time()  
getYs = getall_yClst(allidx, eachct,  400,
                        datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                        betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho)
proc.time() - ptm_mclapply # 


save(getYs, file = "BM2_Feb2023_Projection_Y.RData")


if(0){
  # debug
  
  Ni = length(idx) -1
  npred_theta_posterior_c(rep(0, length(idx)), 2 ,
                          datlist$X[idx[1:Ni],], datlist$Z[idx[1:Ni],], 
                          datlist$time[idx], datlist$y[idx[1:Ni],], datlist$y_observed[idx[1:Ni],], datlist$y1on0,
                          betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                          gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  
  npred_theta_gradient_c(rep(0, length(idx)), 1,
                         datlist$X[idx[1:Ni],], datlist$Z[idx[1:Ni],], 
                         datlist$time[idx], datlist$y[idx[1:Ni],], datlist$y_observed[idx[1:Ni],], datlist$y1on0,
                         betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                         gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  
  pred_theta_V_c(rep(0, length(idx)),  1 ,
                 datlist$X[idx[1:Ni],], datlist$Z[idx[1:Ni],], 
                 datlist$time[idx], datlist$y[idx[1:Ni],], datlist$y_observed[idx[1:Ni],], datlist$y1on0,
                 betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                 gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  
  
  Ni = length(idx) -1
  opt = pred_optsol( 1 ,
                     datlist$X[idx[1:Ni],], datlist$Z[idx[1:Ni],], 
                     datlist$time[idx], datlist$y[idx[1:Ni],],  datlist$y[idx[1:Ni],], datlist$y1on0,
                     betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                     gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  
  mvrnormArma(1, opt[[1]], opt[[2]])
  
  ptm_mclapply <- proc.time()  
  tmp = pred_y_persondraw(idx, length(idx) - 1,
                          datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                          betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                          gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  proc.time() - ptm_mclapply # 
  
  s = 3; idx = allidx[1:eachct[1]]
  s = 1; idx = allidx[(eachct[1]+1):(eachct[1]+eachct[2])]
  tmp = pred_y_persondraw(idx-1, length(idx) - 1,
                          datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                          betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                          gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  
}# 0

#########################################################################
# summarize cluster-specific trajectories and make plots for randomly selected individuals
#load("/Users/yizhenxu/BM2_Feb2023_Projection.RData")
load("BM2_Feb2023_Projection_Y.RData")
load("BM2_Feb2023_Projection_P.RData")
summarray = function(x){
  x = x[!is.na(x)]
  sx = sort(x)
  n = length(x)
  r =  c(mean(x, na.rm=T), round(sx[ceiling(n*0.025)],2), round(sx[floor(n*0.975)],2) )
  if(n <= 3) r = c(NA,NA,NA)
  return(r)
}

yc1summ = apply(getYs[[1]], c(2,3), summarray)
yc2summ = apply(getYs[[2]], c(2,3), summarray)

idlist1 = unique(d$Study.ID)
allidx1 = which(d[,Study.ID] %in% idlist1)
eachct1 =  d[Study.ID %in% idlist1,.N, by=Study.ID]$N

tmp = d[, .N, by=Study.ID]
idlist = tmp[N > 1, Study.ID]
#idlist = unique(d$Study.ID)
allidx = which(d[,Study.ID] %in% idlist)
eachct =  d[Study.ID %in% idlist,.N, by=Study.ID]$N

# restrict getPs to those with at least two observations
getPs1 = getPs[which(allidx1 %in% allidx),] 

Psumm = apply(getPs1,1, summarray)


# conversion time: last diagnosis ascertained by last available diagnosis
d1 = d[!is.na(dx),] # length(unique(d1$RID)) == length(unique(d$RID))
fstS = d1[,dx[1],by=Study.ID]$V1
lstS =  d1[,dx[.N],by=Study.ID]$V1
fstage = d1[, ageori[1], by=Study.ID]$V1

firsttime = function(d1, tag){
  dd = copy(d1)
  dd[, tmp := cumsum(!is.na(dx) & dx == tag), by=Study.ID]
  dd[, target := min(1000, ageori[tmp==1]),by=Study.ID]
  rt = dd[, target[1], by=Study.ID]
  rt$V1[rt$V1==1000] = NA
  return(rt)
}
MCIage = firsttime(d1, "MCI")
ADage = firsttime(d1, "AD")

idshow_earlyMCI = MCIage[V1 <=60, Study.ID]
idshow_lateMCI = MCIage[V1 >=80, Study.ID]
tmp = cbind(MCIage, ADage$V1)
idshow_CN = tmp[is.na(V1) & is.na(V2), Study.ID]

library(ggplot2)
set.seed(1)
# pick one of four
#idshow = sort(sample(1:length(idlist), 10))
idshow = which(idlist %in% idshow_earlyMCI)
#idshow = which(idlist %in% sample(idshow_lateMCI,15))
#idshow = which(idlist %in% sample(idshow_CN,15))

for(i in idshow){
   # ith id in idlist, 1 to 289
  id = idlist[i]
  if(i==1){
    idx =1:sum(eachct[1:i])
  } else {
    idx =( sum(eachct[1:(i-1)])+1):sum(eachct[1:i])
  }
  # idx is the index for getPs1 and getYs (3282 rows)
  # idxindat is the corresponding index for the data (3290 rows)
  idxindat = allidx[idx]
  Ni = length(idx)
  if(Ni >= 5){
    # conversion times
    MCIi = MCIage[Study.ID==id,V1]
    ADi = ADage[Study.ID==id,V1]
    
    #time = datlist$time[idxindat]
    agei = d$ageori[idxindat]
    yi = matrix(y[idxindat,],nrow = Ni)
    
    gplist = vector("list", datlist$K+1)
    
    for(j in 1:datlist$K){
      gpd = as.data.frame(cbind( rbind(t(yc1summ[,idx,j]), t(yc2summ[,idx,j]) ), rep(agei,2), c(rep(1,Ni),rep(2,Ni)) ))
      colnames(gpd) = c("mean", "lower", "upper", "age", "cluster")
      rg = range(c(yi[,j], gpd[,1:3]),na.rm=T)
      obsdi = as.data.frame(cbind(agei, yi[,j])); colnames(obsdi) = c("Age","yy")
      gplist[[j+1]] = ggplot() + geom_line(data = gpd, 
                                           aes(x = age, y = mean, group = cluster, color = as.factor(cluster) ))+
        geom_ribbon(data = gpd, 
                    aes(x = age, y = mean, group = cluster,
                        fill = as.factor(cluster), ymin = lower,
                        ymax = upper), alpha = 0.3, show.legend = F) + 
        geom_line(data = obsdi, aes(x = Age, y = yy)) +
        geom_point(data = obsdi, aes(x = Age, y = yy)) + ylim(rg)+ theme(legend.position = "none")+ labs(y = all[j], x = "Age")
      
    }
    gpd = as.data.frame(cbind( t(Psumm[, idx]),agei , 1) )
    colnames(gpd) = c("mean", "lower", "upper", "age", "cluster")
    
    gplist[[1]] = ggplot() + geom_line(data = gpd, 
                                       aes(x = age, y = mean, group = cluster, color = as.factor(cluster) ))+
      geom_ribbon(data = gpd, 
                  aes(x = age, y = mean, group = cluster,
                      fill = as.factor(cluster), ymin = lower,
                      ymax = upper), alpha = 0.3, show.legend = F) +theme(legend.position = "none")+ labs(y = "P(Clsuter 1)", x = "Age")+ylim(c(0,1))
    if(!is.na(MCIi)){
      gplist[[1]]  = gplist[[1]]+ geom_vline(aes(xintercept=MCIi), colour="#BB0000", linetype="dashed")
    } 
    if(!is.na(ADi)){
      gplist[[1]]  = gplist[[1]]+ geom_vline(aes(xintercept=ADi), colour="#BB0000")
    } 
    
    #ggsave(paste0("i",i,".pdf"), gridExtra::marrangeGrob(grobs = gplist, nrow=3, ncol=2))
    ggsave(paste0("i",i,"_eM.pdf"), gridExtra::marrangeGrob(grobs = gplist, nrow=3, ncol=2))
    #ggsave(paste0("i",i,"_lM.pdf"), gridExtra::marrangeGrob(grobs = gplist, nrow=3, ncol=2))
    #ggsave(paste0("i",i,"_CN.pdf"), gridExtra::marrangeGrob(grobs = gplist, nrow=3, ncol=2))
  }

}

#########################################################################
# Calculate accuracy of the last left out observation

yacc = d[allidx, c(all,"Study.ID"), with=F]
yacc$idx = 1:nrow(yacc)

yacclast = yacc[,.SD[.N], by=Study.ID]
ridx = yacclast$idx

getPs_acc = getPs[allidx,]
getPslast = getPs_acc[ridx,]

getYslast = vector("list",2)
for(cl in 1:2){
  getYslast[[cl]] = getYs[[cl]][, ridx,]
}

Ysamp = array(NA, dim = c(length(ridx), npost, nY))
for(i in 1:length(ridx)){
  for(j in 1:npost){
    prob = getPslast[i,j]
    if(!is.na(prob)){
      memb = rbinom(1,1,prob)
      if(memb){
        Ysamp[i,j,] = getYslast[[1]][j,i,]
      } else {
        Ysamp[i,j,] = getYslast[[2]][j,i,]
      }
    } 
  }
}
Ysumm = apply(Ysamp, c(1,3), summarray)

yacclast = yacclast[,all,with=F]

mae = apply(abs(Ysumm[1,,] - yacclast),2, function(x) mean(x,na.rm=T))
bcov = apply(1*(Ysumm[2,,]<=yacclast)*(Ysumm[3,,]>=yacclast),2, function(x) mean(x,na.rm=T))

round(rbind(mae,bcov),2)

if(0){
  # accuracy of the last obs for ppl whose last status was not AD
  
  #eachct = d[, .N, by=Study.ID]
  #idlist = eachct[N > 1, Study.ID]
  d1 = d[!is.na(dx),] # length(unique(d1$RID)) == length(unique(d$RID))
  lstS_D =  d1[,dx[.N],by=Study.ID]
  lstS_ids = lstS_D[V1 %in% c("MCI", "NORMAL"), Study.ID]
  idlist_CM = idlist[idlist %in% lstS_ids]
  allidx_CM = which(d[,Study.ID] %in% idlist_CM)
  
  yacc = d[allidx_CM, c(all,"Study.ID"), with=F]
  yacc$idx = 1:nrow(yacc)
  
  yacclast = yacc[,.SD[.N], by=Study.ID]
  ridx = yacclast$idx
  
  getPs_acc = getPs[allidx_CM,]
  getPslast = getPs_acc[ridx,]
  
  getYslast = vector("list",2)
  for(cl in 1:2){
    getYslast[[cl]] = getYs[[cl]][, ridx,]
  }
  
  Ysamp = array(NA, dim = c(length(ridx), npost, nY))
  for(i in 1:length(ridx)){
    for(j in 1:npost){
      prob = getPslast[i,j]
      if(!is.na(prob)){
        memb = rbinom(1,1,prob)
        if(memb){
          Ysamp[i,j,] = getYslast[[1]][j,i,]
        } else {
          Ysamp[i,j,] = getYslast[[2]][j,i,]
        }
      } 
    }
  }
  Ysumm = apply(Ysamp, c(1,3), summarray)
  
  yacclast = yacclast[,all,with=F]
  
  mae = apply(abs(Ysumm[1,,] - yacclast),2, function(x) mean(x,na.rm=T))
  bcov = apply(1*(Ysumm[2,,]<=yacclast)*(Ysumm[3,,]>=yacclast),2, function(x) mean(x,na.rm=T))
  
  round(rbind(mae,bcov),2)
}

#########################################################################
#########################################################################
#########################################################################


if(0){
  require(gridExtra)
  #png(height = 700, width = 500, file = "tmp.png")
  pdf(height = 4, width = 5, file = "tmp.pdf")
  do.call("grid.arrange",c(gplist, ncol = 2, nrow = 3))
  #grid.arrange(gplist[[1]],gplist[[2]],gplist[[3]],gplist[[4]], ncol=2)
  #grid.arrange(gplist[[5]],gplist[[6]],gplist[[7]],gplist[[8]], ncol=2)
  #grid.arrange(gplist[[9]],gplist[[10]],gplist[[11]],gplist[[12]], ncol=2)
  dev.off()
  
  png(file = "tmp.png")
  print(gplist[[1]])
  #ggpubr::ggarrange(gplist[[1]],gplist[[j]],gplist[[j]], common.legend = T)
  dev.off()
  
  #########################################################
  i = 187
  id = idlist[i]
  load("BM2_Feb2023_Projection.RData")
  getPs1 = getPs[which(allidx1 %in% allidx),] 
  if(i==1){
    idx =1:sum(eachct[1:i])
  } else {
    idx =( sum(eachct[1:(i-1)])+1):sum(eachct[1:i])
  }
  tmp = getPs1[idx,]; apply(tmp,1,summary)
  
  #i1 = which(idlist1 == id)
  #idxindat = which(d$Study.ID %in% id)
  idxindat = allidx[idx]
  
  s = 1
  tmp1 = calc_prob_laplace_persondraw(idxindat-1, s-1,
                                     datlist$X, datlist$Z, datlist$time, datlist$y,  datlist$y_observed, datlist$y1on0,
                                     betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho,lambda)
  
  tmp = getall_probClst(idxindat-1, eachct[i],  400,
                  datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                  betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho, lambda)
  
  cluster = 1; hidx = idxindat[1]; cc = length(hidx)
  c1 = optsol( cluster, matrix(datlist$X[hidx,],nrow=cc), 
          matrix(datlist$Z[hidx,],nrow=cc), datlist$time[hidx], 
          matrix(datlist$y[hidx,],nrow=cc),matrix(datlist$y_observed[hidx,],nrow=cc), datlist$y1on0,
          betaX[s,,], betat[s,], sigma[s,,], gamma1, 
          gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  n1 = - ntheta_posterior_c(c1[[1]],  cluster, 
                         matrix(datlist$X[hidx,],nrow=cc), 
                         matrix(datlist$Z[hidx,],nrow=cc), datlist$time[hidx], 
                         matrix(datlist$y[hidx,],nrow=cc),matrix(datlist$y_observed[hidx,],nrow=cc), datlist$y1on0,
                         betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                         gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  V1 = symmat(c1[[2]])
  l1 = log(lambda[s,cluster]) + n1 + log(det(V1))/2
  
  cluster = 2
  c2 = optsol( cluster, matrix(datlist$X[hidx,],nrow=cc), 
               matrix(datlist$Z[hidx,],nrow=cc), datlist$time[hidx], 
               matrix(datlist$y[hidx,],nrow=cc),matrix(datlist$y_observed[hidx,],nrow=cc), datlist$y1on0,
               betaX[s,,], betat[s,], sigma[s,,], gamma1, 
               gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  n2 = - ntheta_posterior_c(c2[[1]],  cluster, 
                           matrix(datlist$X[hidx,],nrow=cc), 
                           matrix(datlist$Z[hidx,],nrow=cc), datlist$time[hidx], 
                           matrix(datlist$y[hidx,],nrow=cc),matrix(datlist$y_observed[hidx,],nrow=cc), datlist$y1on0,
                           betaX[s,,], betat[s,], sigma[s,,], gamma1, 
                           gamma2[s,], gamma3[s,], as.double(alpha0[s]), alpha[s,], as.vector(rho[s,]) )
  V2 = symmat(c2[[2]])
  l2 = log(lambda[s,cluster]) + n2 + log(det(V2))/2
  exp(l1 - log(exp(l1)+exp(l2)))
  
  #hl 1 cluster 1 
  #numer -3.33512 logdet -1.12566 lc -5.02227 
  #hl 1 cluster 2 
  #numer -6.68131 logdet -3.8972 lc -9.02277 
  
  di = - ifelse(cluster == 1, 0, alpha0[s])  + Z %*% alpha[s,] + theta
  #mu[person time, outcome] 
  mui = X %*% t(betaX[s,,]) + 
    matrix(time ,ncol=1)%*% matrix(betat[s,], nrow=1) + 
    mapply(function(j) gamma1[j]/(1+exp(-gamma2[s,j]*(di-gamma3[s,j]))), 1:K)
  #gamma1/(1 + exp(-gamma2[s,]*(di-gamma3[s,]) ) );
  
  yi = matrix(y,nrow = Ni)
  obs_ptj = matrix(mapply(function(j) -(1/(2*sigma[s, cluster,j]))*(yi[,j] - mui[,j])^2, 1:K), nrow = Ni)
  
  
  Dmat = exp(-(outer(time, time, "-")^2)/(2*rho[s,cluster]^2)) + diag(1e-10,Ni)
  logpj = sum(apply(obs_ptj, 1, function(x) sum(x,na.rm=T))) - t(theta) %*% solve(Dmat) %*% theta/2 # check
  
  }

#########################################################################
# Calculate accuracy of the last left out observation by age category
#load("/Users/yizhenxu/BM2_Feb2023_Projection.RData")
load("BM2_Feb2023_Projection_Y.RData")
load("BM2_Feb2023_Projection_P.RData")

yacc = d[allidx, c(all,"Study.ID","ageori"), with=F]
yacc$idx = 1:nrow(yacc)

yacclast = yacc[,.SD[.N], by=Study.ID]
ridx = yacclast$idx

time = yacclast[, ageori]

getPs_acc = getPs[allidx,]
getPslast = getPs_acc[ridx,]

getYslast = vector("list",2)
for(cl in 1:2){
  getYslast[[cl]] = getYs[[cl]][, ridx,]
}

Ysamp = array(NA, dim = c(length(ridx), npost, nY))
for(i in 1:length(ridx)){
  for(j in 1:npost){
    prob = getPslast[i,j]
    if(!is.na(prob)){
      memb = rbinom(1,1,prob)
      if(memb){
        Ysamp[i,j,] = getYslast[[1]][j,i,]
      } else {
        Ysamp[i,j,] = getYslast[[2]][j,i,]
      }
    } 
  }
}
Ysumm = apply(Ysamp, c(1,3), summarray)

yacclast = yacclast[,all,with=F]

tcat = floor((time-30)/5); tcat[tcat < 0 | tcat > 11 ] = NA # only between age 30 to 89 


mae = apply(abs(Ysumm[1,,] - yacclast),2, function(x) mean(x,na.rm=T))
bcov = apply(1*(Ysumm[2,,]<=yacclast)*(Ysumm[3,,]>=yacclast),2, function(x) mean(x,na.rm=T))


mae = bcov = matrix(NA, nrow = length(0:11), ncol = nY)
rownames(mae) = rownames(bcov) = (0:11)*5+30 #lower bound of the age categories
for(tcatval in 0:11){
  mae[tcatval+1,] = apply(abs(Ysumm[1,,] - yacclast)[tcat == tcatval,],2, function(x) mean(x,na.rm=T))
  bcov[tcatval+1,] = apply((1*(Ysumm[2,,]<=yacclast)*(Ysumm[3,,]>=yacclast))[tcat == tcatval,],2, function(x) mean(x,na.rm=T))
}

mae = mae[5:11,]# 50 55 60 65 70 75 80
bcov = bcov[5:11,]


cog <- c("MMSCORE","logmem", "DSBACK")
mri <- c("biec.thik", "Hippo_dadjust",
         "Ent_dadjust","MTL1","SPARE_AD") 
csf <- c("ttau", "ptau181", "AB42AB40")
all <- c(cog,mri,csf)

pd = data.frame(MAE = c(mae), Coverage = c(bcov), Biomarker = rep(all, each = nrow(mae)), 
                Age = rep(c("[50,55)","[55,60)", "[60,65)","[65,70)", "[70,75)", "[75,80)", "[80,85)"), length(all)))
library(ggplot2)
plot1 = ggplot(pd, aes(x = Age, y = MAE, color = Biomarker, group = Biomarker)) + 
  geom_point() + geom_line() 
plot2 = ggplot(pd, aes(x = Age, y = Coverage, color = Biomarker, group = Biomarker)) + 
  geom_point() + geom_line() 
require(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plot1)  

ggsave("BM2_opt_Feb2023_ProjectForward_new_plot.pdf", 
       grid.arrange(plot1+ theme(legend.position="none"), 
                    plot2+ theme(legend.position="none"), 
                    mylegend,ncol=3,widths = c(1/3,1/3,1/6)),
       width = 10, height = 4, dpi = 300, units = "in")


