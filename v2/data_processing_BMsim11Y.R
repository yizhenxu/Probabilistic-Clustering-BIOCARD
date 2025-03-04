#5010 (20G)
setwd("./AD")

library(MASS) # mvrnorm
library(Matrix) # bdiag

source("BMsim11Y_setting.R")
###############################################################
### load data for obtaining t and Ni
run = 1
if(run == 1){
  load("BIOCARD_mergedtodx_11152021_lag730.rda")
} else {
  load("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/data/BIOCARD_mergedtodx_11152021_lag730.rda")
}

library(dplyr)
library(data.table)

# BIOCARD with replicate

data <- dat[which(!dat$ID %in% exid), ]

## rename variables
colnames(data)[which(names(data) == "ID")] <- "Study.ID"
colnames(data)[which(names(data) == "C1209A")] <- "logmem"
colnames(data)[which(names(data) == "C1208A")] <- "DSST"
colnames(data)[which(names(data) == "biecthick")] <- "biec.thik"
colnames(data)[which(names(data) == "EDUC")] <- "education"
colnames(data)[which(names(data) == "SEX")] <- "SEX"
colnames(data)[which(names(data) == "C1201D")] <- "MMSCORE"
colnames(data)[which(names(data) == "DIAGNOSIS")] <- "DIAG"
data$age <- (data$DIAGDATE - data$DOB.x)/365.25


## fill in educ sex and apoe
dat <- split(data,data$Study.ID)
for (i in 1:length(dat)) {
  temp <- unique(dat[[i]]$apoe)
  if (length(temp)==2)
    temp <- temp[!is.na(unique(temp))]
  dat[[i]]$apoe <- temp
  
  temp <- unique(dat[[i]]$SEX)
  if (length(temp)==2)
    temp <- temp[!is.na(unique(temp))]
  dat[[i]]$SEX <- temp
  
  temp <- unique(dat[[i]]$education)
  if (length(temp)==2)
    temp <- temp[!is.na(unique(temp))]
  dat[[i]]$education <- temp
}

cog <- c("MMSCORE","logmem", "DSST")
mri <- c("biec.thik", "Hippo_dadjust",
         "Ent_dadjust","MTL1","SPARE_AD")
csf <- c("ttau", "ptau181", "AB42AB40")
all <- c(cog,mri,csf)

data$diagyear = format(data$DIAGDATE, format = "%Y")

data <- data[,c(all,  "age","apoe", "SEX", "DIAG", "education","Study.ID","diagyear")]
data <- data[complete.cases(data[, c("age", "apoe",  "SEX", "education")]), ]
data[,c(cog,mri[-5],"AB42AB40")] <- - data[,c(cog,mri[-5],"AB42AB40")]
data$SEX <- data$SEX - 1
data <- unique(data)

data <- data[rowSums(is.na(data[, all])) != ncol(data[ , all]), ]

data <- data[order( data[,"Study.ID"], data[,"age"] ),]
data$ageori <- as.numeric(data$age)
data$age <- (data$age - 65)/10

data = as.data.table(data)

data = data[order(Study.ID, ageori),]
bslage = data[, ageori[1],by=Study.ID]$V1
###############################################################
### bootstrap N*ND baseline age

set.seed(1)

simd  = data.frame(fID = rep(rep(1:N, each = Ni), ND),  dataID = rep(1:ND, each = Ni*N), visit = rep(1:Ni, N*ND),
                   ageori = rep(sample(bslage, N*ND, replace = T),each = Ni)+ 
                     rep( seq(0, gap * (Ni-1), by = gap), N*ND) ) 

simd$age = (simd$ageori - 65)/10

###############################################################
### simulate covariates

Nall = N*ND
simd = as.data.table(simd)
simd[, tmpID := paste0(dataID,"_",fID)]
simd$X1 = rep( rbinom(Nall,1,0.6) , simd[,.N,by=tmpID]$N)
simd$X2 = rep(rnorm(Nall), simd[,.N,by=tmpID]$N)

###############################################################
### simulate gaussian process random effect

theta = matrix(NA, nrow = nrow(simd), ncol = L)

kfun = function(ti,rho){
  timat = expand.grid(x=ti,y=ti)
  timat$k = exp(-(timat$x-timat$y)^2 / (2*rho^2))
  return(matrix(timat$k, nrow=length(ti)))
}

library(MASS)
v = unique(simd$tmpID)
ti = simd$age[simd$tmpID==v[1]]

for(l in 1:L){
  theta[,l] = c(t(mvrnorm(length(v), rep(0,Ni), kfun(ti, rho[l]))))
}


#ggplot(data = simd, aes(x = age, y = theta1, group = fID))+geom_line()
#ggplot(data = simd, aes(x = age, y = theta2, group = fID))+geom_line()

###############################################################
### generate latent scores and outcome

d = matrix(NA, nrow = nrow(simd), ncol = L)
d[,1] = alphat*simd$age + theta[,1]
d[,2] = -alpha0 + alphat*simd$age + theta[,2]

Y = array(NA, dim = c(nrow(simd), K, L))
f = array(NA, dim = c(nrow(simd), K, L))
Xb = matrix(NA, nrow = nrow(simd), ncol =K)
for(k in 1:K){
  for(l in 1:L){
    f[, k, l] = gamma1[k]/(1 + exp(-gamma2[k]*(d[,l] - gamma3[k])))
    Xb[,k] =  beta0[k] + betat[k]*simd$age + betaX1[k]*simd$X1 + betaX2[k]*simd$X2
    Y[, k, l] = rnorm(nrow(simd), Xb[,k] + f[,k,l], sqsig[l,k]) 
  }
}

simd$memb1 = rep( rbinom(Nall,1,lambda1) , simd[,.N,by=tmpID]$N)

for( k in 1:K) 
  simd[, (paste0("Y",k)) := (simd$memb1==1)*Y[,k,1] + (simd$memb1==0)*Y[,k,2]]

simd$group = (simd$memb1==1)*1 + (simd$memb1==0)*2


###############################################################

save(simd, d, file = "simdata11.RData")


if(0){
  
  fp = array(NA, dim = c(nrow(simd), K, L))
  for(k in 1:K){
    for(l in 1:L){
      fp[, k, l] = 1/(1 + exp(-gamma2[k]*(d[,l] - gamma3[k])))
    }
  }
  
  idx = sample(1:nrow(simd),1000)
  dp = d[idx,]; fp = fp[idx,,]  
  ###
  
  plot(dp[,1],fp[,1,1], xlim = c(-8,8), ylim = c(0,1.05), cex = 0.2, col="red"); points(dp[,2],fp[,1,2], cex = 0.4, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3), pch = 19); 
 
  for(k in 2:K){
    points(dp[,1],fp[,k,1], cex = 0.2, col="red"); points(dp[,2],fp[,k,2], cex = 0.4, col = rgb(red = 0, green = 0, blue = 1, alpha = 0.3), pch = 19);  
  }
    
  range(dp[,1]);  range(dp[,2])
  for(l in 1:L){
    rdsd <- density(dp[,l])
    lines(rdsd$x, rdsd$y, lty = l, col = rgb(red = 0, green = 1, blue =0 ))
  }
  
  
  
}



