setwd("./AD")
library(rstan)

scode = "
functions {
  
  real[] subvec2(vector x, int pos, int len){
    real result[len];
    for(i in 1:len){
      result[i] = x[pos + i - 1];
    }
    return result;
  }

}
data {
  int<lower=1> Nsubs; //number of patients
  int<lower=1> Nobs; //total number of records
  int<lower=1> K; //number of outcomes, 3 here
  int<lower=1> Npreds; //number of predictors in X
  int<lower=1> NZ; //number of predictors in Z
  
  matrix[Nobs,Npreds] X;
  matrix[Nobs,NZ] Z;
  
  vector[Nobs] time; //time points for all records
  
  int lsub[Nsubs]; //element i being total number of time points for person i
   
  matrix[Nobs, K] y; //fill missingness with zeros would be fine -- IMPORTANT
  
  vector[K] gamma1;
 
  int catlen[3]; //2,4,3
  int catsum[3]; //0,2,6 = cumsum(c(0,catlen))[1:3] 
}
parameters {
  
  real betat[K]; 
  vector[Npreds] betaX[K];
  
  real<lower=0> Sig[K];
  
  vector[NZ] alpha;
  
  vector[Nobs] z;
  real<lower=0, upper = 1.5> rho_theta;
  
  vector<lower=0>[K] gamma2;
  
  vector[catlen[1]] gamma3_cog;
  vector[catlen[2]] gamma3_mri;
  vector[catlen[3]] gamma3_csf;
  
  vector[3] gamma3mu;
  
}
transformed parameters {

  vector[Nobs] theta;
  
  matrix[Nobs,K] mu = rep_matrix(0, Nobs, K);
  vector[K] gamma3;

  vector[Nobs] d;
  
{
    int pos=1;
    for(i in 1:Nsubs){
      real timei[lsub[i]];
      matrix[lsub[i], lsub[i]] cov;
      matrix[lsub[i], lsub[i]] L_cov;
      timei = subvec2(time, pos, lsub[i]);
      cov = cov_exp_quad(timei, 1.0, rho_theta)+ diag_matrix(rep_vector(1e-10, lsub[i]));
      L_cov = cholesky_decompose(cov);
      theta[pos:(pos+lsub[i]-1)] = L_cov * segment(z, pos, lsub[i]); // theta_i ~ N(0, L_cov_i * L_cov_i')
      pos = pos + lsub[i];
    } 
}
  
  for(l in 1:Nobs){
    d[l] = Z[l] * alpha + theta[l];
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
  
  for (l in 1:Nobs) {
    for(k in 1:K){
        mu[l, k] = X[l] * betaX[k] + time[l]*betat[k]  + gamma1[k]/(1 + exp(-gamma2[k]*(d[l]-gamma3[k]) ) );
    }
  }
  
}
model {
  int pos;
  
  betat ~ normal(0,10);
  for(k in 1:K){
    betaX[k] ~ normal(0,10); // mean, sd
  }
  
  gamma2 ~ gamma(3,1);
  gamma3_cog ~ normal(gamma3mu[1], 1);
  gamma3_mri ~ normal(gamma3mu[2], 1);
  gamma3_csf ~ normal(gamma3mu[3], 1);
  gamma3mu ~ normal(0,2);
  
  alpha ~ normal(0,1);
  
  Sig ~ inv_gamma(0.01, 0.01);
  
  z ~ normal(0,1);
  
  for(n in 1:Nobs){
        for(k in 1:K){
          target += normal_lpdf(y[n,k] | mu[n,k], sqrt(Sig[k])); 
        }//K
  }//Nobs

}
"

fn <- "mod1.stan"
if (file.exists(fn))  file.remove(fn)

# Write a new one
writeLines(scode, fn)
rstan:::rstudio_stanc(fn)


# time correlation use time_orig
scode = "
functions {
 
  real[] subvec2(vector x, int pos, int len){
    real result[len];
    for(i in 1:len){
      result[i] = x[pos + i - 1];
    }
    return result;
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
  int lsub[Nsubs]; //element i being total number of time points for person i
  int endidx[Nsubs]; // cumulative_sum(lsub);
  
  matrix[Nobs, K] y; //fill missingness with zeros would be fine -- IMPORTANT
  
  vector[K] gamma1;

  int L; //number of mixture components
  
  int catlen[3]; //2,4,3
  int catsum[3]; //0,2,6 = cumsum(c(0,catlen))[1:3]
}
parameters {
  
  real betat[K]; 
  vector[Npreds] betaX[K];
  
  vector<lower=0>[K] Sig[L];
  
  vector[Nobs] z[L];
  vector<lower=0, upper = 1.5>[L] rho_theta;
  
  vector[NZ] alpha;
  ordered[L-1] alpha0; 
  
  simplex[L] lambda;          // mixing proportions
  
  vector<lower=0>[K] gamma2;
  
  vector[catlen[1]] gamma3_cog;
  vector[catlen[2]] gamma3_mri;
  vector[catlen[3]] gamma3_csf;
  
  vector[3] gamma3mu;
}
transformed parameters {

  matrix[Nobs,L] theta = rep_matrix(0,Nobs,L);
  
  real mu[Nobs,K,L]; 
  vector[K] gamma3;
  
  real d[Nobs, L];
  
{
    int pos=1;
    for(i in 1:Nsubs){
      real timei[lsub[i]];
      timei = subvec2(time, pos, lsub[i]);
      for(l in 1:L){
        matrix[lsub[i], lsub[i]] covl;
        matrix[lsub[i], lsub[i]] L_covl;
        covl = cov_exp_quad(timei, 1.0, rho_theta[l])+ diag_matrix(rep_vector(1e-10, lsub[i]));
        L_covl = cholesky_decompose(covl);
        theta[pos:(pos+lsub[i]-1), l] = L_covl * segment(z[l], pos, lsub[i]); // theta_i ~ N(0, L_cov_i * L_cov_i')
        //sub_col(theta, pos, l,lsub[i])
      }
      pos = pos + lsub[i];
    } 
}

  for(n in 1:Nobs){
    d[n,1] = Z[n] * alpha + theta[n,1];
    for(l in 2:L){
      d[n,l] = -alpha0[l-1] + Z[n] * alpha + theta[n,l];
    }
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
      for(l in 1:L){
          mu[n, k, l] = X[n] * betaX[k] + time[n]*betat[k] + gamma1[k]/(1 + exp(-gamma2[k]*(d[n,l]-gamma3[k]) ) );
      }
    }
  }
}
model {

  betat ~ normal(0,10);
  for(k in 1:K){
    betaX[k] ~ normal(0,10); // mean, sd
  }
  
  gamma2 ~ gamma(3,1);
  gamma3_cog ~ normal(gamma3mu[1], 0.5);
  gamma3_mri ~ normal(gamma3mu[2], 0.5);
  gamma3_csf ~ normal(gamma3mu[3], 0.5);
  gamma3mu ~ normal(0,2);
  
  alpha0 ~ gamma(2,1.5);
  alpha ~ normal(0,1);
  
  for(l in 1:L){
    z[l] ~ normal(0,1);
    Sig[l] ~ inv_gamma(0.01, 0.01);
  }
  
  for(i in 1:Nsubs){
    vector[L] lps = log(lambda);
    for(l in 1:L){
      for(n in (endidx[i]-lsub[i]+1):endidx[i] ){
        for(k in 1:K){
          lps[l] += normal_lpdf(y[n,k] | mu[n,k,l], sqrt(Sig[l,k])); 
        }//K
      }//Ni
    }//L
    target += log_sum_exp(lps);
  }//i

}
"

fn <- "mod2.stan"
if (file.exists(fn))  file.remove(fn)

# Write a new one
writeLines(scode, fn)
rstan:::rstudio_stanc(fn)