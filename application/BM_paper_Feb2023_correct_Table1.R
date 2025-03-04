
type = 'B'; ageopt = 2 # 

route = ""
source(paste0(route, "PlotFunctions.R"))

rdnum = 2
summ = function(x){
  x = x[!is.na(x)]
  sx = sort(x)
  n = length(x)
  r = paste0(round(mean(x, na.rm=T),rdnum), " (",round(sx[ceiling(n*0.025)],rdnum),", ",round(sx[floor(n*0.975)],rdnum),")")
  if(n <= 3) r = ""
  return(r)
}

load("BIOCARD_align2dx_wrep_022023.rda")
################################################################
#source(paste0(route, "data_processing_Feb2023_noWMH.R"))
# if(run == 1){
#   load("BIOCARD_align2dx_wrep_022023.rda")
# } else {
#   route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/"
#   load(paste0(route, "BIOCARD_align2dx_wrep_022023.rda"))
#   #load("/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/BIOCARD_align2dx_wrep_022023.rda")
# }

library(dplyr)
library(data.table)

## BIOCARD with replicate
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
data <- do.call(rbind, dat)
data <- unique(data)

#data$Whole_brain_WMH = dat$Whole_brain_WMH/dat$ICV

cog <- c("MMSCORE","logmem", "DSBACK")
mri <- c("biec.thik", "Hippo_dadjust",
         "Ent_dadjust","MTL1","SPARE_AD") 
csf <- c("ttau", "ptau181", "AB42AB40")
all <- c(cog,mri,csf)


data$diagyear = format(data$DIAGDATE, format = "%Y")

data <- data[,c(all,  "age","apoe", "SEX", "DIAG", "education","Study.ID","diagyear")]
data <- data[complete.cases(data[, c("age", "apoe",  "SEX", "education")]), ]
#data[,c(cog,mri[-5],"AB42AB40")] <- - data[,c(cog,mri[-5],"AB42AB40")]
data$SEX <- data$SEX - 1
data <- unique(data)

data <- data[rowSums(is.na(data[, all])) != ncol(data[ , all]), ]

data <- data[order( data[,"Study.ID"], data[,"age"] ),]
data$ageori <- data$age

#data$age <- scale(data$age)
if(ageopt == 1) data$age = data$age - 65
if(ageopt == 2) data$age = (data$age - 65)/10

data$eduori = data$education
data$education <- scale(data$education)

## scale biomarkers
# for (i in all) {
#   data[, i] <- scale(data[, i])
# }

data$dx <- NA
data$dx[data$DIAG=="DEMENTIA"] <- "AD" 
data$dx[data$DIAG=="MCI"] <- "MCI" 
data$dx[data$DIAG=="NORMAL"] <- "NORMAL" 
data$dx[data$DIAG=="IMPAIRED NOT MCI"] <- "NORMAL"

data <- data[, c(cog, mri, csf, "age", "apoe", "SEX", "education","dx","Study.ID","ageori","diagyear","DIAG","eduori")]

gamma1 = gamma1_BIOC

library(data.table); library(rstan)
#d = as.data.table(data)
d = copy(data)
dim(d); length(unique(d$Study.ID)) # 299ppl 3312 rows

d1 = d[!is.na(d$dx),] # length(unique(d1$Study.ID)) == length(unique(d$Study.ID))
d1 = as.data.table(d1)
fstS = d1[,dx[1],by=Study.ID]
incID = fstS[V1 == "NORMAL", Study.ID]

d = d[d$Study.ID %in% incID,] # 297ppl 3289 rows

d = as.data.table(d)
################################################################

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
