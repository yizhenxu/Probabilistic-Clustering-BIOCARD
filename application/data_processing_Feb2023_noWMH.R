# ageopt = 1: age - 65
# ageopt = 2: (age - 65)/10

run = 0 # 0 for laptop, 1 for JHPCE
if(run == 1){
  load("gamma1_Feb2023_noWMH.RData") # running gamma1_calculation_Feb2023.R
} else {
  route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/"
  load(paste0(route, "gamma1_Feb2023_noWMH.RData")) # running gamma1_calculation_Feb2023.R
}

if(type == "A"){
  if(run == 1){
    load("ADNImerge_artificialtime_annual_122021_wmtl.rda")
  } else {
    route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/data/"
    load(paste0(route, "ADNImerge_artificialtime_annual_122021_wmtl.rda"))
    #load("/Volumes/GoogleDrive/My Drive/Desktop/2020 spring/Zheyu/data/ADNImerge_artificialtime_annual_122021_wmtl.rda")
  }
  
  #dat_ADNI$WM_HYPOINTENSITIES_VOLUME = dat_ADNI$WM_HYPOINTENSITIES_VOLUME/dat_ADNI$ICV
  
  cog <- c("MMSCORE","logmem", "DSPANBAC") #"DSST"
  mri <- c("biec.thik", "bihippo_dadj","biec.vol_dadj","MTL1")
  csf <- c("TAU", "PTAU", "ABETA")
  all <- c(cog,mri,csf)
  
  dat_ADNI <- dat_ADNI[order(dat_ADNI[, "RID"], dat_ADNI[, "age"] ), ]
  dat_ADNI$PTGENDER <- dat_ADNI$PTGENDER - 1
  dat_ADNI$DSST[which(dat_ADNI$DSST<0)] <- NA
  dat_ADNI$logmem[which(dat_ADNI$logmem<0)] <- NA
  dat_ADNI$AVLTTOTL[which(dat_ADNI$AVLTTOTL<0)] <- NA
  dat_ADNI$education[which(dat_ADNI$education<0)] <- NA
  
  
  ## fill in educ sex and apoe
  dat_ADNI <- dat_ADNI[,c(cog,mri,csf, "age", "apoe", "PTGENDER", "dx", "education","RID","EXAMDATE.x", "Phase")]
  dat <- split(dat_ADNI,dat_ADNI$RID)
  for (i in 1:length(dat)) {
    temp <- unique(dat[[i]]$apoe)
    if (length(temp)==2)
      temp <- temp[!is.na(unique(temp))]
    dat[[i]]$apoe <- temp
    
    temp <- unique(dat[[i]]$PTGENDER)
    if (length(temp)==2)
      temp <- temp[!is.na(unique(temp))]
    dat[[i]]$PTGENDER <- temp
    
    temp <- unique(dat[[i]]$education)
    if (length(temp)==2)
      temp <- temp[!is.na(unique(temp))]
    dat[[i]]$education <- temp
  }
  dat_ADNI <- do.call(rbind, dat)
  dat_ADNI <- unique(dat_ADNI)
  
  ## all subjects (allA == TRUE) or normal subjects (allA == FALSE)
  library(dplyr)
  dat <- dat_ADNI[order(dat_ADNI[, "RID"], dat_ADNI[, "age"]), c("RID", "dx")]
  dat <- dat[complete.cases(dat$dx),]
  if(allA == TRUE){
    output = unique(dat$RID)
  } else {
    dat <- split(dat,dat$RID)
    output <- NULL
    for (i in 1:length(dat)) {
      if(dat[[i]]$dx[1] == 1)
        output <- rbind(output, dat[[i]]$RID[1])
    }
  }
  
  dat_ADNI <- dat_ADNI[dat_ADNI$RID %in% output,]
  
  dat_ADNI <- dat_ADNI[,c(cog,mri,csf, "age", "apoe", "PTGENDER", "dx", "education","RID","EXAMDATE.x")]
  dat_ADNI <- dat_ADNI[complete.cases(dat_ADNI[, c("age",  "PTGENDER", "education","apoe")]), ]
  dat_ADNI[,c(cog,mri,"ABETA")] <- - dat_ADNI[,c(cog,mri,"ABETA")]
  dat_ADNI <- unique(dat_ADNI)
  dat_ADNI[sapply(dat_ADNI, is.nan)] <- NA
  dat_ADNI <- dat_ADNI[rowSums(is.na(dat_ADNI[, all])) != ncol(dat_ADNI[ , all]), ]
  
  dat_ADNI$ageori <- dat_ADNI$age
  
  #dat_ADNI$age <- as.numeric(scale(dat_ADNI$age))
  if(ageopt == 1) dat_ADNI$age = dat_ADNI$age - 65
  if(ageopt == 2) dat_ADNI$age = (dat_ADNI$age - 65)/10
  
  dat_ADNI$eduori = dat_ADNI$education
  dat_ADNI$education <- as.numeric(scale(dat_ADNI$education))
  
  
  ## scale biomarkers
  for (i in all) {
    dat_ADNI[, i] <- scale(dat_ADNI[, i])
  }
  
  dat_ADNI$dx1 <- NA
  dat_ADNI$dx1[dat_ADNI$dx==3] <- "AD" 
  dat_ADNI$dx1[dat_ADNI$dx==2] <- "MCI" 
  dat_ADNI$dx1[dat_ADNI$dx==1] <- "NORMAL" 
  dat_ADNI$dx <- dat_ADNI$dx1
  
  
  dat_ADNI$diagyear = format(dat_ADNI$EXAMDATE.x, format = "%Y")
  
  output <- dat_ADNI[, c(cog, mri, csf, "age", "apoe", "PTGENDER", "education","dx","RID","ageori","diagyear","eduori")]
  
  gamma1 = gamma1_ADNI
  
  library(data.table); library(rstan)
  #d = as.data.table(data)
  d = copy(output)
  dim(d); length(unique(d$RID)) # 785ppl 4539 rows # 2115 ppl 11089 rows
  
  #d[, ct := .N, by=Study.ID]
  #d = d[ct >=2, ] #291 ppl 3304 rows
  y <- as.matrix(d[, which(colnames(d) %in% all)])
  xname = c("apoe", "PTGENDER", "education")
  x <- cbind(1, as.matrix(d[, which(colnames(d) %in% xname)]))
  z<- as.matrix( d[, which(colnames(d) == "apoe")])
  
  yfill0 = y
  yfill0[is.na(yfill0)] = 0
  
  table(z); dim(y)
  
  id = "RID"
}

if(type == "B"){
  if(run == 1){
    load("BIOCARD_align2dx_wrep_022023.rda")
  } else {
    route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/"
    load(paste0(route, "BIOCARD_align2dx_wrep_022023.rda"))
    #load("/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/BIOCARD_align2dx_wrep_022023.rda")
  }
  
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
  data[,c(cog,mri[-5],"AB42AB40")] <- - data[,c(cog,mri[-5],"AB42AB40")]
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
  for (i in all) {
    data[, i] <- scale(data[, i])
  }
  
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
  
  y <- as.matrix(d[, which(colnames(d) %in% all)])
  xname = c("apoe", "SEX", "education")
  x <- cbind(1, as.matrix(d[, which(colnames(d) %in% xname)]))
  z<- as.matrix( d[, which(colnames(d) == "apoe")])
  
  yfill0 = y
  yfill0[is.na(yfill0)] = 0
  
  table(z); dim(y)
  
  id = "Study.ID"
}

d = as.data.table(d)
subjvec = rep(1:length(unique(d[,get(id)])),d[,.N,by=get(id)]$N)

# For BIOCARD
#d[, ct := .N, by=Study.ID]
#d = d[ct >=2, ] #289 ppl 3281 rows


