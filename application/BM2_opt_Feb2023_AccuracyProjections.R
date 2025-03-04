

if(0){
  #fidx = 5
  #source("./AD/BM2_opt_Feb2023_AccuracyProjections.R")
  rm(list = ls())
  setwd("./AD")
  for(fidx in 1:5){
    source("./BM2_opt_Feb2023_AccuracyProjections.R")
  }
  # 1 for 24, 2-4 for 28, 5 for 17
}


library(MASS) # mvrnom

type = 'B'; agesopt = 2 # 


#route = "/users/yxu2/AD/"
route = ""
source("PlotFunctions.R")

fn = paste0("BM2_opt2_Feb2023_noWMH_Accuracy_fold",fidx,".RData")


library(rstan)

type = "B"; ageopt = 2
source(paste0(route, "data_processing_Feb2023_noWMH.R"))

if(file.exists(fn)) load(fn)

library("Rcpp")
Rcpp::sourceCpp("cfun_AD.cpp")
# to compile, (shift command .) to show hidden file, 
# download gfortran from https://mac.r-project.org/tools/
# locate under /usr/local/gfortran/lib the dylib's reported missing in the errors
# create new folder /usr/local/lib/ and move those dylib's into this folder


### trim data to be the test set for the fidx_th fold

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

d = d[Study.ID %in% testid,]

# data list preparation

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


# define allidx and eachct

allidx1 = 1:nrow(d)

eachct1 =  d[,.N, by=Study.ID]$N

ptm_mclapply <- proc.time()  
getPs = getall_probClst(allidx1, eachct1,  400,
                        datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                        betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho, lambda)
proc.time() - ptm_mclapply # 

#IDshow = rep(d$Study.ID, eachct1)
#show = cbind(IDshow, apply(getPs, 1, function(x) mean(x,na.rm=T)),round(getPs[,1:400],2))

save(getPs, file = paste0("BM2_Feb2023_Projection_P",fidx,".RData"))

#########################################################################
#########################################################################

allidx1 = 1:nrow(d)

eachct1 =  d[,.N, by=Study.ID]$N

ptm_mclapply <- proc.time()  
getYs = getall_yClst(allidx1, eachct1,  400,
                     datlist$X, datlist$Z, datlist$time, datlist$y, datlist$y_observed, datlist$y1on0,
                     betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho)
proc.time() - ptm_mclapply # 


save(getYs, file = paste0("BM2_Feb2023_Projection_Y",fidx,".RData"))

#########################################################################
#########################################################################
# Calculate accuracy of the test set over time

time = d$ageori
tcat = floor((time-30)/5); tcat[tcat < 0 | tcat > 11 ] = NA # only between age 30 to 89 

load(paste0("BM2_Feb2023_Projection_Y",fidx,".RData"))
load(paste0("BM2_Feb2023_Projection_P",fidx,".RData"))

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

allidx1 = 1:nrow(d)

eachct1 =  d[,.N, by=Study.ID]$N

Ysamp = array(NA, dim = c(nrow(d), npost, nY))
for(i in 1:nrow(d)){
  for(j in 1:npost){
    prob = getPs[i,j]
    if(!is.na(prob)){
      memb = rbinom(1,1,prob)
      if(memb){
        Ysamp[i,j,] = getYs[[1]][j,i,]
      } else {
        Ysamp[i,j,] = getYs[[2]][j,i,]
      }
    } 
  }
}

Ysumm = apply(Ysamp, c(1,3), summarray) # nrow(d) x nY

yacc = d[allidx1, all, with=F]

mae = bcov = matrix(NA, nrow = length(0:11), ncol = nY)
rownames(mae) = rownames(bcov) = (0:11)*5+30 #lower bound of the age categories
for(tcatval in 0:11){
  mae[tcatval+1,] = apply(abs(Ysumm[1,,] - yacc)[tcat == tcatval,],2, function(x) mean(x,na.rm=T))
  bcov[tcatval+1,] = apply((1*(Ysumm[2,,]<=yacc)*(Ysumm[3,,]>=yacc))[tcat == tcatval,],2, function(x) mean(x,na.rm=T))
}

rl = rbind(mae,bcov)
save(rl, file = paste0("Accuracy",fidx,".RData"))

rm(list = ls())

#########################################################################
#########################################################################
#########################################################################

if(0){
  nfold = 5
  load( paste0("/Users/yizhenxu/Documents/Accuracy",1,".RData"))
  rtab = rl
  for(fidx in 2:nfold){
    load( paste0("/Users/yizhenxu/Documents/Accuracy",fidx,".RData"))
    rtab = rtab + rl
  }
  rtab = rtab / nfold
  rtab;

mae = rtab[1:12,] 
bcov = rtab[13:24,]

#(0:11)*5+30 = 30 35 40 45 50 55 60 65 70 75 80 85

#mae = mae[4:10,]#45 50 55 60 65 70 75
#bcov = bcov[4:10,]
mae = mae[6:9,]#55 60 65 70 
bcov = bcov[6:9,]


cog <- c("MMSCORE","logmem", "DSBACK")
mri <- c("biec.thik", "Hippo_dadjust",
         "Ent_dadjust","MTL1","SPARE_AD") 
csf <- c("ttau", "ptau181", "AB42AB40")
all <- c(cog,mri,csf)
require(latex2exp)
allname = c("MMSE", "Log Mem", "DSBACK", 
           "EC Thick", "Hippo Vol", "EC Vol", "MTL", "SPARE-AD",
           "t-tau", TeX('p-tau$_{181p}$'), TeX("A$\\beta_{42}$/A$\\beta_{40}$"))

pd = data.frame(MAE = c(mae), Coverage = c(bcov), Biomarker = as.factor(rep(all, each = nrow(mae)))  )
#pd$Age = rep(c("[45,50)","[50,55)","[55,60)", "[60,65)","[65,70)", "[70,75)", "[75,80)"), length(all))
pd$Age = rep(c("[55,60)", "[60,65)","[65,70)", "[70,75)"), length(all))
library(ggplot2)
plot1 = ggplot(pd, aes(x = Age, y = MAE, color = Biomarker, group = Biomarker)) + 
  geom_point() + geom_line() +
  scale_color_hue(labels = unname(
    c(TeX("A$\\beta_{42}$/A$\\beta_{40}$"),"EC Thick", "DSBACK","EC Vol",
      "Hippo Vol","Log Mem", "MMSE",  "MTL", TeX('p-tau$_{181p}$'),
      "SPARE-AD","t-tau" )
  )) #levels(pd$Biomarker)
plot2 = ggplot(pd, aes(x = Age, y = Coverage, color = Biomarker, group = Biomarker)) + 
  geom_point() + geom_line() +
  scale_color_hue(labels = unname(
    c(TeX("A$\\beta_{42}$/A$\\beta_{40}$"),"EC Thick", "DSBACK","EC Vol",
      "Hippo Vol","Log Mem", "MMSE",  "MTL", TeX('p-tau$_{181p}$'),
      "SPARE-AD","t-tau" )
  )) #levels(pd$Biomarker)
require(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plot1)  

ggsave("/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/BM2_opt_Feb2023_AccuracyProjection_plot.pdf", 
       grid.arrange(plot1+ theme(legend.position="none"), 
                    plot2+ theme(legend.position="none"), 
                    mylegend,ncol=3,widths = c(1/3,1/3,1/6)),
       width = 7.5, height = 3, dpi = 300, units = "in")

}
