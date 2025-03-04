#qrsh -l mem_free=30G,h_vmem=30G,h_fsize=30G
#module load conda_R/4.0

setwd("./AD")

route = "/Users/yizhenxu/Library/CloudStorage/GoogleDrive-yizhen_xu@alumni.brown.edu/My Drive/Desktop/2020 spring/Zheyu/code/"
setwd(route)

a = 3;b=1
#x = seq(0,10,0.01)
#y = dgamma(x,shape=a, rate =b)
#plot(x,y,"l")
1-pgamma(8,shape=a, rate =b)#1.4%
pgamma(0.5,shape=a, rate =b)#1.4%

# time correlation use time_orig
scode1 = "
functions {
  matrix Sig_theta_chol(real sig2, real rho, vector t) {
  int r = 1;
  int Q = num_elements(t);
  matrix[Q, Q] G;

  for (i in 1:(Q-1)) {
    G[i, i] = sig2;
    for (j in (i+1):Q) {
      G[i, j] = sig2 * pow(rho, fabs(t[i]-t[j]));
      G[j, i] = G[i, j];
    }
  }
  G[Q, Q] = sig2;

  G = G + diag_matrix(rep_vector(1e-10, Q));
  return cholesky_decompose(G);
  }
  
  matrix submat(matrix x, int[] idx) { 
    int count = size(idx); 
   
    matrix[count, count] result; 
    for(i in 1:count){
      for(j in 1:count){
        result[i, j] = x[idx[i], idx[j]];
      }
    }
     
     return result; 
  }
  
  vector subvec(row_vector x, int[] idx){
    int count = size(idx);
    vector[count] result;
    
    for(i in 1:count){
      result[i] = x[idx[i]];
    }
    return result;
  }
  
  vector subvec1(vector x, int[] idx){
    int count = size(idx);
    vector[count] result;
    
    for(i in 1:count){
      result[i] = x[idx[i]];
    }
    return result;
  }
  
  real[] subvec2(vector x, int pos, int len){
    array[len] real result;

    for(i in 1:len){
      result[i] = x[pos + i - 1];
    }
    return result;
  }
  
  matrix block_diag(matrix qq, int n) {
    int r = 1;
    int sz = rows(qq);
    int Q = sz * n;
    matrix[Q, Q] G;
    for(i in 1:Q){
      for(j in 1:Q){
        G[i,j] = 0;
      }
    }
    while (r < Q) {
      G[r:(r + sz - 1), r:(r + sz - 1)] = qq;
      r += sz;
    }

    return G;
  }
  
  matrix my_block_diag(matrix qq1,matrix qq2,matrix qq3) {
    
    matrix[9,9] G;
    for(i in 1:9){
      for(j in 1:9){
        G[i,j] = 0;
      }
    }
    G[1:3,1:3] = qq1; 
    G[4:6,4:6] = qq2; 
    G[7:9,7:9] = qq3; 

    return G;
  }
}
data {
  int<lower=1> Nsubs; //number of patients
  int<lower=1> Nobs; //total number of records
  int<lower=1> K; //number of outcomes, 9 here
  int<lower=1> Npreds; //number of predictors in X
  int<lower=1> NZ; //number of predictors in Z
  
  matrix[Nobs,Npreds] X;
  matrix[Nobs,NZ] Z;
  
  vector[Nobs] time; //time points for all records
  array[Nsubs] int lsub; //element i being total number of time points for person i
  array[Nsubs] int endidx; // cumulative_sum(lsub);
   
  array[Nobs] int subjvec; // the lth entry is the subject idx to which time[l] belongs to

  matrix[Nobs, K] y_observed; //indicator of non-missing in y
  
  matrix[Nobs, K] y; //fill missingness with zeros would be fine -- IMPORTANT
  
  vector[K] gamma1;
  real rho_lower;
  real rho_upper;// obtained from prior tuning of within person age differences

  int L; //number of mixture components
  
  array[Nobs] int y1on0; // indicator of MMScore on the lower boundary
  
  real<lower=0> a; 
  real<lower=0> b;
  
  array[3] int catlen;  //3,4,3
  array[3] int catsum;  //0,3,7 = cumsum(c(0,catlen))[1:3]
}
parameters {
  
  array[K] real betat; //real betat[K]; 
  array[K] vector[Npreds] betaX; //vector[Npreds] betaX[K];
  
  array[L] vector<lower=0>[K] Sig; //vector<lower=0>[K] Sig[L];
  
  vector[Nobs] z1;
  real<lower=0> rho_theta1;
  vector[Nobs] z2;
  real<lower=0> rho_theta2;
  
  vector[NZ] alpha;
  real alpha0;
  
  simplex[L] lambda;          // mixing proportions
  
  vector<lower=0>[K] gamma2;
  
  vector[catlen[1]] gamma3_cog;
  vector[catlen[2]] gamma3_mri;
  vector[catlen[3]] gamma3_csf;
  
  vector[3] gamma3mu;
}
transformed parameters {
  vector[Nobs] theta1;
  vector[Nobs] theta2;

  array[Nobs,K,L] real mu; //real mu[Nobs,K,L]; 
  vector[K] gamma3;
  
  array[Nobs, L] real d; // real d[Nobs, L];
  
{
    int pos=1;
    for(i in 1:Nsubs){
      array[lsub[i]] real timei; //real timei[lsub[i]];
      matrix[lsub[i], lsub[i]] cov1;
      matrix[lsub[i], lsub[i]] cov2;
      matrix[lsub[i], lsub[i]] L_cov1;
      matrix[lsub[i], lsub[i]] L_cov2;
      timei = subvec2(time, pos, lsub[i]);
      cov1 = gp_exp_quad_cov(timei, 1.0, rho_theta1)+ diag_matrix(rep_vector(1e-10, lsub[i]));
      L_cov1 = cholesky_decompose(cov1);
      cov2 = gp_exp_quad_cov(timei, 1.0, rho_theta2)+ diag_matrix(rep_vector(1e-10, lsub[i]));
      L_cov2 = cholesky_decompose(cov2);
      theta1[pos:(pos+lsub[i]-1)] = L_cov1 * segment(z1, pos, lsub[i]); // theta_i ~ N(0, L_cov_i * L_cov_i')
      theta2[pos:(pos+lsub[i]-1)] = L_cov2 * segment(z2, pos, lsub[i]); 
      pos = pos + lsub[i];
    } 
}

  for(n in 1:Nobs){
    d[n,1] = Z[n] * alpha + theta1[n];
    d[n,2] = -alpha0 + Z[n] * alpha + theta2[n];
  }
 
  
 for(n in 1:catlen[1]){
      gamma3[catsum[1]+n] = gamma3_cog[n] ;
  }
  for(n in 1:catlen[2]){
      gamma3[catsum[2]+n] = gamma3_mri[n] ;
  }
  for(n in 1:catlen[3]){
      gamma3[catsum[3]+n] = gamma3_csf[n] ;
  }
  
  mu = rep_array(0, Nobs, K, L);
  for(k in 1:K){
    for(n in 1:Nobs){
      if(y_observed[n, k]){
        for(l in 1:L){
          mu[n, k, l] = X[n] * betaX[k] + time[n]*betat[k] + gamma1[k]/(1 + exp(-gamma2[k]*(d[n,l]-gamma3[k]) ) );
        }
      }
    }
  }
}
model {

  betat ~ normal(0,10);
  for(k in 1:K){
    betaX[k] ~ normal(0,10); // mean, sd
  }
  
  gamma2 ~ gamma(a,b);
  gamma3_cog ~ normal(gamma3mu[1], 1);
  gamma3_mri ~ normal(gamma3mu[2], 1);
  gamma3_csf ~ normal(gamma3mu[3], 1);
  gamma3mu ~ normal(0,2);
  
  z1 ~ normal(0,1);
  rho_theta1 ~ inv_gamma(rho_lower, rho_upper); // P[rho < l] \approx 0.01, P[rho > u] \approx 0.01
  z2 ~ normal(0,1);
  rho_theta2 ~ inv_gamma(rho_lower, rho_upper); // P[rho < l] \approx 0.01, P[rho > u] \approx 0.01
  
  alpha0 ~ gamma(2,1.5);
  alpha ~ normal(0,1);
  
  for(l in 1:L){
    Sig[l] ~ inv_gamma(0.01, 0.01);
  }
  
  for(i in 1:Nsubs){
    vector[L] lps = log(lambda);
    for(l in 1:L){
      for(n in (endidx[i]-lsub[i]+1):endidx[i] ){
        if(y_observed[n, 1]){
          if(y1on0[n]==1){
            lps[l] += normal_lcdf(y[n,1] | mu[n,1,l], sqrt(Sig[l,1]));
          } else {
            lps[l] += normal_lpdf(y[n,1] | mu[n,1,l], sqrt(Sig[l,1]));
          }
        }
        for(k in 2:K){
          if(y_observed[n, k]){
            lps[l] += normal_lpdf(y[n,k] | mu[n,k,l], sqrt(Sig[l,k]));
          }  
        }//K
      }//Ni
    }//L
    target += log_sum_exp(lps);
  }//i

}
generated quantities {
  matrix[Nsubs, L] logpyz;
  real log_lik = 0;
  
  log_lik += normal_lpdf(z1 | 0,1) ;
  log_lik += normal_lpdf(z2 | 0,1) ;
  
  for(i in 1:Nsubs){
    vector[L] lps = log(lambda);
    for(l in 1:L){
      for(n in (endidx[i]-lsub[i]+1):endidx[i] ){
        if(y_observed[n, 1]){
          if(y1on0[n]==1){
            lps[l] += normal_lcdf(y[n,1] | mu[n,1,l], sqrt(Sig[l,1]));
          } else {
            lps[l] += normal_lpdf(y[n,1] | mu[n,1,l], sqrt(Sig[l,1]));
          }
        }
        for(k in 2:K){
          if(y_observed[n, k]){
            lps[l] += normal_lpdf(y[n,k] | mu[n,k,l], sqrt(Sig[l,k]));
          }  
        }//K
      }//Ni
      logpyz[i,l] = lps[l];
    }//L
    log_lik += log_sum_exp(lps);
  }//i
  
}
"

fn <- "mod1.stan"
if (file.exists(fn))  file.remove(fn)

if(0){
  path = '/users/yxu2/R/4.3.x'
  remove.packages(c("rstan", "StanHeaders"))
  install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")), lib = path)
  library("rstan", lib.loc=path)
  stan_version() # >= 2.26.1
}


#library(rstan)
# Write a new one
source("scode.R")
writeLines(scode0, fn)
rstan:::rstudio_stanc(fn)


###

type = 'B'; ageopt = 2 # 

source("data_processing_Feb2023_noWMH.R")

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

eachct = d[, .N, by=Study.ID]
idlist = eachct[N > 1, Study.ID]
d = d[Study.ID %in% idlist,]
d = d[,.SD[1:(.N-1)], by=Study.ID]
y <- d[, all,with=F]
xname = c("apoe", "SEX", "education")
x <- cbind(1, as.matrix(d[, xname,with=F]))

zname = c("apoe", "age")
z<- matrix( unlist(d[, zname, with=F]), ncol=length(zname))
z = cbind(z, z[,1]*z[,2]) # interaction term


yfill0 = y
yfill0[is.na(yfill0)] = 0
subjvec = rep(1:length(unique(d[,get(id)])),d[,.N,by=get(id)]$N)

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
                 y1on0 = as.numeric(1*(!is.na(y[,1]) & y[,1] == min(y[,1],na.rm=T)) ), 
                 a = a,b=b,
                 catlen = catlen,
                 catsum = cumsum(c(0,catlen))[1:3])


nit = 1000; nwup = 600; nc = 10; seednum = 8


#nit = 200; nwup = 100; nc = 3; seednum = 8

options(mc.cores = nc) #parallel::detectCores()
fit = stan(file = "mod1.stan", data = staninput,
           pars = c("gamma2","gamma3",
                    "betat","betaX","Sig",
                    "alpha0","alpha","rho_theta1","rho_theta2",
                    "lp__","log_lik", "logpyz", "lambda","d","gamma3mu") , include = T,
           chains = nc, iter = nit, warmup = nwup, thin = 1, seed = seednum)

save(fit, file = paste0("BM2_opt2_Feb2023_noWMH_Accuracy.RData"))  #here
#
###############################################################

fn <- "mod2.stan"
if (file.exists(fn))  file.remove(fn)
source("scode.R")
writeLines(scode1, fn)
rstan:::rstudio_stanc(fn)

options(mc.cores = nc) #parallel::detectCores()
fit = stan(file = fn, data = staninput,
           pars = c("gamma2","gamma3",
                    "betat","betaX","Sig",
                    "alpha0","alpha","rho_theta1","rho_theta2",
                    "lp__","log_lik", "logpyz", "lambda","d","gamma3mu") , include = T,
           chains = nc, iter = nit, warmup = nwup, thin = 1, seed = seednum)

save(fit, file = paste0("BM2_opt2_Feb2023_noWMH_Accuracy_newSyntax.RData"))  #here
#

