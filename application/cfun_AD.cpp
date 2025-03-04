#include <RcppArmadillo.h> // might cause problem, installing gfortran-6.1.pkg from the R website solves the problem: https://stackoverflow.com/questions/35999874/mac-os-x-r-error-ld-warning-directory-not-found-for-option
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]

static double const log2pi = std::log(2.0 * M_PI);

// [[Rcpp::export]]
arma::mat symmat(const arma::mat& sigma) {
  //lower triangular copied to upper
  int n = sigma.n_cols;
  arma::mat A = sigma;
  for(int i = 0; i < n; i ++){
    for(int j = i+1; j < n; j ++){
      A(i,j) = A(j,i);
    }
  }
  
  return A;
}

// [[Rcpp::export]]
double ntheta_posterior_c(const arma::vec& theta,
                          int cluster,
                          const arma::mat& X,
                          const arma::mat& Z,
                          const arma::vec& time,
                          const arma::mat& y,
                          const arma::mat& yobs,
                          const arma::vec& y1on0,
                          const arma::mat& betaX, 
                          const arma::vec& betat, 
                          const arma::mat& sigma,
                          const arma::vec& gamma1,
                          const arma::vec& gamma2,
                          const arma::vec& gamma3,
                          double alpha0,
                          const arma::vec& alpha,
                          const arma::vec& rho){

  double logpj;
  
  int Ni = X.n_rows;
  int K = y.n_cols;
  
  arma::vec di = zeros<vec>(Ni);
  arma::mat mui(Ni, K);
  arma::mat obs_ptj = zeros(Ni, K);
  arma::mat Dmat(Ni, Ni);
  arma::mat err = arma::eye( Ni,Ni );
  
  if(cluster == 1){
    di = Z * alpha + theta;
  } else {
    di = - alpha0 + Z * alpha + theta;
  }
  
  mui = X * betaX.t() + time * betat.t() ;
  
  for(int j = 0; j < K; j ++){
    mui.col(j) += gamma1(j)/(1+exp(-gamma2(j)*(di-gamma3(j))));
    if(j == 0){
      for(int i = 0; i < Ni; i++){
        if(y1on0[i]==1){
          obs_ptj(i,j) = yobs(i,j) * R::pnorm( y(i,j),  mui(i,j) , sqrt(sigma(cluster-1,j)), 1,1); //logged
        } else {
          obs_ptj(i,j) = yobs(i,j) * ( -(1/(2*sigma(cluster-1,j)))*pow(y(i,j) - mui(i,j),2) );
        }
      }
    } else {
      obs_ptj.col(j) = yobs.col(j) % ( -(1/(2*sigma(cluster-1,j)))*pow(y.col(j) - mui.col(j),2) );
    }
  }
  
  
  for(int i = 0; i < Ni; i++) {
    Dmat.col(i) = pow(time[i] - time,2);
  }
  Dmat = exp(-Dmat/(2*pow(rho(cluster-1),2))) + (1e-10)*err;
  logpj = accu(obs_ptj) - as_scalar(theta.t() * inv(Dmat) * theta )/2;
  return -logpj;
}

// [[Rcpp::export]]
arma::vec ntheta_gradient_c(const arma::vec& theta, 
                            int cluster,
                            const arma::mat& X,
                            const arma::mat& Z,
                            const arma::vec& time,
                            const arma::mat& y,
                            const arma::mat& yobs,
                            const arma::vec& y1on0,
                            const arma::mat& betaX, 
                            const arma::vec& betat, 
                            const arma::mat& sigma,
                            const arma::vec& gamma1,
                            const arma::vec& gamma2,
                            const arma::vec& gamma3,
                            double alpha0,
                            const arma::vec& alpha,
                            const arma::vec& rho){
  
  int Ni = X.n_rows;
  int K = y.n_cols;
  arma::vec di = zeros<vec>(Ni);
  arma::mat mui(Ni, K);
  
  arma::mat E(Ni, K);
  arma::mat A(Ni, K);
  arma::mat B(Ni, K);
  arma::vec C(Ni);

  arma::mat r(Ni, K);
  arma::mat Dmat(Ni, Ni);
  arma::mat err = arma::eye( Ni,Ni );
  arma::vec drvt = zeros<vec>(Ni);
  if(cluster == 1){
    di = Z * alpha + theta;
  } else {
    di = - alpha0 + Z * alpha + theta;
  }
  
  mui = X * betaX.t() + time * betat.t() ;
  
  for(int j = 0; j < K; j ++){
    E.col(j) = exp(-gamma2(j)*(di-gamma3(j)));
    mui.col(j) += gamma1(j)/(1+E.col(j));
    A.col(j) = y.col(j) - mui.col(j);
    B.col(j) = E.col(j)/pow(1+E.col(j),2);
    
    if(j == 0){
      C = exp(-pow(A.col(j),2) / (2*sigma(cluster-1,j)) )/ sqrt(2*M_PI*sigma(cluster-1,j));
      for(int i = 0; i < Ni; i++){
        if(y1on0[i]==1){
          r(i,j) =  - gamma1(j)*gamma2(j)* yobs(i,j) * B(i,j) * C(i) / R::pnorm( y(i,j),  mui(i,j) , sqrt(sigma(cluster-1,j)), 1,0);
        } else {
          r(i,j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs(i,j) * A(i,j) * B(i,j);
        }
      }
    } else {
      r.col(j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs.col(j) % A.col(j) % B.col(j);
    }
    
    drvt += r.col(j);
  }
  
  for(int i = 0; i < Ni; i++) {
    Dmat.col(i) = pow(time[i] - time,2);
  }
  Dmat = exp(-Dmat/(2*pow(rho(cluster-1),2))) + (1e-10)*err;
  drvt -= inv(Dmat) * theta;
  return -drvt;
}

// [[Rcpp::export]]
arma::mat theta_V_c(const arma::vec& theta, 
                    int cluster,
                    const arma::mat& X,
                    const arma::mat& Z,
                    const arma::vec& time,
                    const arma::mat& y,
                    const arma::mat& yobs,
                    const arma::vec& y1on0,
                    const arma::mat& betaX, 
                    const arma::vec& betat, 
                    const arma::mat& sigma,
                    const arma::vec& gamma1,
                    const arma::vec& gamma2,
                    const arma::vec& gamma3,
                    double alpha0,
                    const arma::vec& alpha,
                    const arma::vec& rho){
  int Ni = X.n_rows;
  int K = y.n_cols;
  arma::vec di = zeros<vec>(Ni);
  arma::mat mui(Ni, K);
  
  arma::mat E(Ni, K);
  arma::mat A(Ni, K);
  arma::mat B(Ni, K);
  arma::vec C(Ni);
  double Dij;
  
  arma::mat aderv(Ni, K);
  
  arma::mat Ad(Ni, K);
  arma::mat Bd(Ni, K);
  arma::vec Cd(Ni);
  double Ddij;
  
  arma::mat tmp = zeros(Ni, Ni);
  arma::mat invtmp = zeros(Ni, Ni);
  
  arma::mat Dmat(Ni, Ni);
  arma::mat err = arma::eye( Ni,Ni );
  
  
  if(cluster == 1){
    di = Z * alpha + theta;
  } else {
    di = - alpha0 + Z * alpha + theta;
  }
  
  mui = X * betaX.t() + time * betat.t() ;
  
  for(int j = 0; j < K; j ++){
    E.col(j) = exp(-gamma2(j)*(di-gamma3(j)));
    mui.col(j) += gamma1(j)/(1+E.col(j));
    A.col(j) = y.col(j) - mui.col(j);
    B.col(j) = E.col(j)/pow(1+E.col(j),2);
    Ad.col(j) = -gamma1(j)*gamma2(j)*B.col(j);
    Bd.col(j) = gamma2(j) * B.col(j) % (1 - 2/(1+E.col(j)));
    
    if(j == 0){
      C = exp(-pow(A.col(j),2) / (2*sigma(cluster-1,j)) )/ sqrt(2*M_PI*sigma(cluster-1,j));
      Cd = gamma1(j)*gamma2(j)* A.col(j) % B.col(j) % C / sigma(cluster-1,j);
      for(int i = 0; i < Ni; i++){
        if(y1on0[i]==1){
          Dij = R::pnorm( y(i,j),  mui(i,j) , sqrt(sigma(cluster-1,j)), 1,0);
          Ddij = - gamma1(j)*gamma2(j) * B(i,j) * C(i);
          aderv(i,j) =  gamma1(j)*gamma2(j)* yobs(i,j) * 
            ( B(i,j) * C(i) * Ddij / pow(Dij,2) - (Bd(i,j) * C(i) + B(i,j) * Cd(i))/Dij);
        } else {
          aderv(i,j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs(i,j) * 
                        (Ad(i,j) * B(i,j) + A(i,j) * Bd(i,j));
        }
      }
    } else {
      aderv.col(j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs.col(j) %
                      (Ad.col(j) % B.col(j) + A.col(j) % Bd.col(j));
    }
    
  }
  
  for(int i = 0; i < Ni; i++){
    tmp(i,i) = -accu(aderv.row(i));
    invtmp(i,i) = 1/tmp(i,i);
  }
  
  for(int i = 0; i < Ni; i++) {
    Dmat.col(i) = pow(time[i] - time,2);
  }
  Dmat = exp(-Dmat/(2*pow(rho(cluster-1),2))) + (1e-10)*err;
  
  //return arma::inv(tmp + arma::inv_sympd(Dmat) );
  // inv(B)*inv(inv(A)+inv(B))*inv(A) = inv(A+B)
  return Dmat*arma::inv_sympd(invtmp + Dmat)*invtmp;
}

// [[Rcpp::export]]
Rcpp::List optsol(int cluster,
                  const arma::mat& X,
                  const arma::mat& Z,
                  const arma::vec& time,
                  const arma::mat& y,
                  const arma::mat& yobs,
                  const arma::vec& y1on0,
                  const arma::mat& betaX, 
                  const arma::vec& betat, 
                  const arma::mat& sigma,
                  const arma::vec& gamma1,
                  const arma::vec& gamma2,
                  const arma::vec& gamma3,
                  double alpha0,
                  const arma::vec& alpha,
                  const arma::vec& rho){
  
  Rcpp::List res(2);
  
  // Extract R's optim function
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  int Ni = X.n_rows;

  arma::vec init_val= zeros<vec>(Ni);
  
  // Call the optim function from R in C++ 
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_val,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&ntheta_posterior_c),
                                 Rcpp::_["gr"]     = Rcpp::InternalFunction(&ntheta_gradient_c),
                                 Rcpp::_["method"] = "BFGS",
                                 // Pass in the other parameters as everything
                                 // is scoped environmentally
                                 Rcpp::_["cluster"] = cluster,
                                 Rcpp::_["X"] = X,
                                 Rcpp::_["Z"] = Z,
                                 Rcpp::_["time"] = time,
                                 Rcpp::_["y"] = y,
                                 Rcpp::_["yobs"] = yobs,
                                 Rcpp::_["y1on0"] = y1on0,
                                 Rcpp::_["betaX"] = betaX,
                                 Rcpp::_["betat"] = betat,
                                 Rcpp::_["sigma"] = sigma,
                                 Rcpp::_["gamma1"] = gamma1,
                                 Rcpp::_["gamma2"] = gamma2,
                                 Rcpp::_["gamma3"] = gamma3,
                                 Rcpp::_["alpha0"] = alpha0,
                                 Rcpp::_["alpha"] = alpha,
                                 Rcpp::_["rho"] = rho);
  
  // Extract out the estimated parameter values
  arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);
  
  arma::mat I = theta_V_c(out, cluster,
                          X,Z,time,y,yobs,y1on0,
                          betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho);
  res[0] = out;
  res[1] = I;
  return res;
}

// [[Rcpp::export]]
arma::vec calc_prob_laplace_persondraw( const arma::uvec& idx,
                                        int s,
                                        const arma::mat& X,
                                        const arma::mat& Z,
                                        const arma::vec& time,
                                        const arma::mat& y,
                                        const arma::mat& yobs,
                                        const arma::vec& y1on0,
                                        const arma::cube& betaX, 
                                        const arma::mat& betat, 
                                        const arma::cube& sigma,
                                        const arma::vec& gamma1,
                                        const arma::mat& gamma2,
                                        const arma::mat& gamma3,
                                        const arma::vec& alpha0,
                                        const arma::mat& alpha,
                                        const arma::mat& rho,
                                        const arma::mat& lambda){
  int Ni = idx.n_elem;
  
  double numer;
  double logdet;
  //double tmp;
  arma::mat logtmp(Ni, 2);
  double c0, c1;
  //Rcpp::List output(2);
  //arma::mat probClst(Ni, 2);
  arma::vec probClst(Ni);
  
  for(int hl = 1; hl <= Ni; hl ++){
    for(int cluster = 1; cluster <=2; cluster++){
      Rcpp::List opt = optsol(cluster, X.rows(idx(0), idx(hl-1)), 
                              Z.rows(idx(0), idx(hl-1)),
                              time.subvec(idx(0), idx(hl-1)),
                              y.rows(idx(0), idx(hl-1)), 
                              yobs.rows(idx(0), idx(hl-1)), 
                              y1on0.subvec(idx(0), idx(hl-1)),
                              betaX.row(s), betat.row(s).t(), sigma.row(s), 
                              gamma1, gamma2.row(s).t(), gamma3.row(s).t(), 
                              alpha0(s), alpha.row(s).t(), rho.row(s).t() );
      //double denom = dmvnorm_arma(opt[0], opt[0], opt[1], true)(0);
      Rcout  << "hl "<< hl<< " cluster "<< cluster << " \n";
      numer = - ntheta_posterior_c(opt[0], cluster, 
                                   X.rows(idx(0), idx(hl-1)), 
                                   Z.rows(idx(0), idx(hl-1)),
                                   time.subvec(idx(0), idx(hl-1)),
                                   y.rows(idx(0), idx(hl-1)), 
                                   yobs.rows(idx(0), idx(hl-1)), 
                                   y1on0.subvec(idx(0), idx(hl-1)),
                                   betaX.row(s), betat.row(s).t(), sigma.row(s), 
                                   gamma1, gamma2.row(s).t(), gamma3.row(s).t(), 
                                   alpha0(s), alpha.row(s).t(), rho.row(s).t()  ) ;
      //probClst(hl-1, cluster-1) = lambda(s,cluster-1) * mean(numer/denom);
      arma::mat V = symmat(opt[1]);
      logdet = sum(arma::log(arma::eig_sym(V)));
      logtmp(hl - 1, cluster - 1) = log(lambda(s,cluster-1)) + numer + logdet/2;
      //Rcout  << "numer "<< numer<< " logdet "<< logdet << " lc "<< logtmp(hl - 1, cluster - 1)<<" \n";
    }//cluster
    
    //tmp = accu(probClst.row(hl-1));
    //probClst.row(hl-1) /= tmp;
    c0 =  logtmp(hl - 1, 0) ;
    c1 =  logtmp(hl - 1, 1) ;
    if((c0 < 0) & (c1 < 0) & (exp(c0-c1) == 0)){
      //probClst(hl-1, 0) = 0;
      probClst(hl-1) = 0;
    } else {
      if((c0 < 0) & (c1 < 0) & (exp(c1-c0) == 0)){
        //probClst(hl-1, 0) = 1;
        probClst(hl-1) = 1;
      } else {
        //probClst(hl-1, 0) = exp(logtmp(hl-1, 0) -
        //  log(accu(exp(logtmp(hl-1, 0)) + exp(logtmp(hl-1, 1)) )));
        probClst(hl-1) = exp( c0 - log(exp(c0) + exp(c1)) );
        
      }
    }
    //probClst(hl-1, 1) = 1- probClst(hl-1, 0);
    
  }//hl
  //output[0] = logtmp;
  //output[1] = probClst;
  //return output;
  return probClst;
}

// [[Rcpp::export]]
arma::mat getall_probClst( const arma::uvec& allidx,
                           const arma::vec& eachct,
                           int smax,
                           const arma::mat& X,
                           const arma::mat& Z,
                           const arma::vec& time,
                           const arma::mat& y,
                           const arma::mat& yobs,
                           const arma::vec& y1on0,
                           const arma::cube& betaX, 
                           const arma::mat& betat, 
                           const arma::cube& sigma,
                           const arma::vec& gamma1,
                           const arma::mat& gamma2,
                           const arma::mat& gamma3,
                           const arma::vec& alpha0,
                           const arma::mat& alpha,
                           const arma::mat& rho,
                           const arma::mat& lambda){
  
  arma::mat output(allidx.n_elem, smax);
  int ct = 0;
  int ctl;
  for(uword i = 0; i < eachct.n_elem; i++){
    ctl = ct + eachct(i);
    arma::uvec idx(eachct(i)) ;
    idx = allidx.subvec(ct, ctl-1 );
    
    for(int s = 0; s < smax; s ++){
      output.submat(ct, s, ctl-1, s) =  calc_prob_laplace_persondraw( idx-1, s,
                    X,Z,time,y,yobs, y1on0,
                    betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho, lambda);
      
    }
    ct = ctl;
  }
  return output;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, const arma::vec& mu, const arma::mat& sigma) {
  int q = sigma.n_cols;
  arma::mat Y = arma::randn(q, n); // ncols x n matrix
  // chol.t() * chol = sigma
  arma::mat R(sigma.n_rows,sigma.n_rows);
  arma::mat err = arma::eye(sigma.n_rows,sigma.n_rows);
  arma::mat sigerr = sigma;
  bool success = false;
  int ct = 0;
  while(success == false)
  {
    success = arma::chol(R, sigerr);
    
    if(success == false)
    {
      sigerr += err * 1e-6;
      ct += 1;
      if(ct == 10) return arma::repmat(mu, 1, n); // check?
    }
  }
  return arma::repmat(mu, 1, n) + R.t() * Y;
  //return arma::repmat(mu, 1, n) + arma::chol(sigma).t() * Y;
}

arma::vec Mahalanobis(arma::mat const &x, 
                      arma::vec const &center, 
                      arma::mat const &cov) {
  //https://gallery.rcpp.org/articles/dmvnorm_arma/
  arma::mat x_cen = x.t();
  x_cen.each_col() -= center;
  arma::solve(x_cen, arma::trimatl(chol(cov).t()), x_cen);
  x_cen.for_each( [](arma::mat::elem_type& val) { val = val * val; } );
  return arma::sum(x_cen, 0).t();    
}

// [[Rcpp::export]]
arma::vec dmvnorm_arma(arma::mat const &x, 
                       arma::vec const &mean, 
                       arma::mat const &sigma, 
                       bool const logd = false) { 
  arma::vec const distval = Mahalanobis(x,  mean, sigma);
  double const logdet = sum(arma::log(arma::eig_sym(sigma)));
  arma::vec const logretval = 
    -( (x.n_cols * log2pi + logdet + distval)/2  ) ;
    
    if (logd)
      return logretval;
    return exp(logretval);
}

// [[Rcpp::export]]
double npred_theta_posterior_c(const arma::vec& theta, 
                          int cluster,
                          const arma::mat& X,
                          const arma::mat& Z,
                          const arma::vec& time,
                          const arma::mat& y,
                          const arma::mat& yobs,
                          const arma::vec& y1on0,
                          const arma::mat& betaX, 
                          const arma::vec& betat, 
                          const arma::mat& sigma,
                          const arma::vec& gamma1,
                          const arma::vec& gamma2,
                          const arma::vec& gamma3,
                          double alpha0,
                          const arma::vec& alpha,
                          const arma::vec& rho){
  
  int pN = time.n_elem; //theta and time should be of length pN, others of Ni rows
  double logpj;
  int Ni = X.n_rows;
  int K = y.n_cols;
  arma::vec di = zeros<vec>(Ni);
  arma::mat mui(Ni, K);
  arma::mat obs_ptj = zeros(Ni, K);
  arma::mat Dmat(pN, pN);
  arma::mat err = arma::eye( pN, pN );
  
  if(cluster == 1){
    di = Z * alpha + theta.subvec(0,Ni-1);
  } else {
    di = - alpha0 + Z * alpha + theta.subvec(0,Ni-1);
  }
  
  mui = X * betaX.t() + time.subvec(0,Ni-1) * betat.t() ;
  
  for(int j = 0; j < K; j ++){
    mui.col(j) += gamma1(j)/(1+exp(-gamma2(j)*(di-gamma3(j))));
    if(j == 0){
      for(int i = 0; i < Ni; i++){
        if(y1on0[i]==1){
          obs_ptj(i,j) = yobs(i,j) * R::pnorm( y(i,j),  mui(i,j) , sqrt(sigma(cluster-1,j)), 1,1); //logged
        } else {
          obs_ptj(i,j) = yobs(i,j) * ( -(1/(2*sigma(cluster-1,j)))*pow(y(i,j) - mui(i,j),2) );
        }
      }
    } else {
      obs_ptj.col(j) = yobs.col(j) % ( -(1/(2*sigma(cluster-1,j)))*pow(y.col(j) - mui.col(j),2) );
    }
  }
  
  for(int i = 0; i < pN; i++) {
    Dmat.col(i) = pow(time[i] - time,2);
  }
  Dmat = exp(-Dmat/(2*pow(rho(cluster-1),2))) + (1e-10)*err;
  logpj = accu(obs_ptj) - as_scalar(theta.t() * inv(Dmat) * theta )/2;
  return -logpj;
}

// [[Rcpp::export]]
arma::vec npred_theta_gradient_c(const arma::vec& theta,  
                            int cluster,
                            const arma::mat& X,
                            const arma::mat& Z,
                            const arma::vec& time,
                            const arma::mat& y,
                            const arma::mat& yobs,
                            const arma::vec& y1on0,
                            const arma::mat& betaX, 
                            const arma::vec& betat, 
                            const arma::mat& sigma,
                            const arma::vec& gamma1,
                            const arma::vec& gamma2,
                            const arma::vec& gamma3,
                            double alpha0,
                            const arma::vec& alpha,
                            const arma::vec& rho){
  
  int pN = time.n_elem;
  int Ni = X.n_rows;
  int K = y.n_cols;
  arma::vec di = zeros<vec>(Ni);
  arma::mat mui(Ni, K);
  
  arma::mat E(Ni, K);
  arma::mat A(Ni, K);
  arma::mat B(Ni, K);
  arma::vec C(Ni);
  
  arma::mat r(Ni, K);
  arma::mat Dmat(pN, pN);
  arma::mat err = arma::eye( pN, pN);
  arma::vec drvt = zeros<vec>(pN);
  
  if(cluster == 1){
    di = Z * alpha + theta.subvec(0,Ni-1);
  } else {
    di = - alpha0 + Z * alpha + theta.subvec(0,Ni-1);
  }
  
  mui = X * betaX.t() + time.subvec(0,Ni-1) * betat.t() ;
  
  for(int j = 0; j < K; j ++){
    E.col(j) = exp(-gamma2(j)*(di-gamma3(j)));
    mui.col(j) += gamma1(j)/(1+E.col(j));
    A.col(j) = y.col(j) - mui.col(j);
    B.col(j) = E.col(j)/pow(1+E.col(j),2);
    
    if(j == 0){
      C = exp(-pow(A.col(j),2) / (2*sigma(cluster-1,j)) )/ sqrt(2*M_PI*sigma(cluster-1,j));
      for(int i = 0; i < Ni; i++){
        if(y1on0[i]==1){
          r(i,j) =  - gamma1(j)*gamma2(j)* yobs(i,j) * B(i,j) * C(i) / R::pnorm( y(i,j),  mui(i,j) , sqrt(sigma(cluster-1,j)), 1,0);
        } else {
          r(i,j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs(i,j) * A(i,j) * B(i,j);
        }
      }
    } else {
      r.col(j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs.col(j) % A.col(j) % B.col(j);
    }
    
    drvt.subvec(0,Ni-1) += r.col(j);
  }
  
  for(int i = 0; i < pN; i++) {
    Dmat.col(i) = pow(time[i] - time,2);
  }
  Dmat = exp(-Dmat/(2*pow(rho(cluster-1),2))) + (1e-10)*err;
  drvt -= inv(Dmat) * theta;
  return -drvt;
}

// [[Rcpp::export]]
arma::mat pred_theta_V_c(const arma::vec& theta,  
                    int cluster,
                    const arma::mat& X,
                    const arma::mat& Z,
                    const arma::vec& time,
                    const arma::mat& y,
                    const arma::mat& yobs,
                    const arma::vec& y1on0,
                    const arma::mat& betaX, 
                    const arma::vec& betat, 
                    const arma::mat& sigma,
                    const arma::vec& gamma1,
                    const arma::vec& gamma2,
                    const arma::vec& gamma3,
                    double alpha0,
                    const arma::vec& alpha,
                    const arma::vec& rho){
  int pN = time.n_elem;
  int Ni = X.n_rows;
  int K = y.n_cols;
  
  arma::vec di = zeros<vec>(Ni);
  arma::mat mui(Ni, K);
  
  arma::mat E(Ni, K);
  arma::mat A(Ni, K);
  arma::mat B(Ni, K);
  arma::vec C(Ni);
  double Dij;
  
  arma::mat aderv(Ni, K);
  
  arma::mat Ad(Ni, K);
  arma::mat Bd(Ni, K);
  arma::vec Cd(Ni);
  double Ddij;

  arma::mat tmp = zeros(pN, pN);
  //arma::mat invtmp = zeros(pN, pN);
  
  arma::mat Dmat(pN, pN);
  arma::mat err = arma::eye( pN, pN);
  
  
  if(cluster == 1){
    di = Z * alpha + theta.subvec(0,Ni-1);
  } else {
    di = - alpha0 + Z * alpha + theta.subvec(0,Ni-1);
  }
  mui = X * betaX.t() + time.subvec(0,Ni-1) * betat.t() ;
  
  for(int j = 0; j < K; j ++){
    E.col(j) = exp(-gamma2(j)*(di-gamma3(j)));
    mui.col(j) += gamma1(j)/(1+E.col(j));
    A.col(j) = y.col(j) - mui.col(j);
    B.col(j) = E.col(j)/pow(1+E.col(j),2);
    Ad.col(j) = -gamma1(j)*gamma2(j)*B.col(j);
    Bd.col(j) = gamma2(j) * B.col(j) % (1 - 2/(1+E.col(j)));
    
    if(j == 0){
      C = exp(-pow(A.col(j),2) / (2*sigma(cluster-1,j)) )/ sqrt(2*M_PI*sigma(cluster-1,j));
      Cd = gamma1(j)*gamma2(j)* A.col(j) % B.col(j) % C / sigma(cluster-1,j);
      for(int i = 0; i < Ni; i++){
        if(y1on0[i]==1){
          Dij = R::pnorm( y(i,j),  mui(i,j) , sqrt(sigma(cluster-1,j)), 1,0);
          Ddij = - gamma1(j)*gamma2(j) * B(i,j) * C(i);
          aderv(i,j) =  gamma1(j)*gamma2(j)* yobs(i,j) * 
            ( B(i,j) * C(i) * Ddij / pow(Dij,2) - (Bd(i,j) * C(i) + B(i,j) * Cd(i))/Dij);
        } else {
          aderv(i,j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs(i,j) * 
            (Ad(i,j) * B(i,j) + A(i,j) * Bd(i,j));
        }
      }
    } else {
      aderv.col(j) =  (1/sigma(cluster-1,j))*gamma1(j)*gamma2(j)* yobs.col(j) %
        (Ad.col(j) % B.col(j) + A.col(j) % Bd.col(j));
    }
    
  }
  
  for(int i = 0; i < Ni; i++){
    tmp(i,i) = -accu(aderv.row(i));
    //invtmp(i,i) = 1/tmp(i,i); // this is problembatic, for Ni+1 to pN terms
  }
  
  for(int i = 0; i < pN; i++) {
    Dmat.col(i) = pow(time[i] - time,2);
  }
  Dmat = exp(-Dmat/(2*pow(rho(cluster-1),2))) + (1e-10)*err;
  
  return arma::inv(tmp + arma::inv_sympd(Dmat) );
  // inv(B)*inv(inv(A)+inv(B))*inv(A) = inv(A+B)
  //return Dmat*arma::inv_sympd(invtmp + Dmat)*invtmp;
}

// [[Rcpp::export]]
Rcpp::List pred_optsol(int cluster,
                  const arma::mat& X,
                  const arma::mat& Z,
                  const arma::vec& time,
                  const arma::mat& y,
                  const arma::mat& yobs,
                  const arma::vec& y1on0,
                  const arma::mat& betaX, 
                  const arma::vec& betat, 
                  const arma::mat& sigma,
                  const arma::vec& gamma1,
                  const arma::vec& gamma2,
                  const arma::vec& gamma3,
                  double alpha0,
                  const arma::vec& alpha,
                  const arma::vec& rho){
  
  Rcpp::List res(2);
  
  // Extract R's optim function
  Rcpp::Environment stats("package:stats"); 
  Rcpp::Function optim = stats["optim"];
  
  arma::vec init_val= zeros<vec>(time.n_elem);
  
  // Call the optim function from R in C++ 
  Rcpp::List opt_results = optim(Rcpp::_["par"]    = init_val,
                                 // Make sure this function is not exported!
                                 Rcpp::_["fn"]     = Rcpp::InternalFunction(&npred_theta_posterior_c),
                                 Rcpp::_["gr"]     = Rcpp::InternalFunction(&npred_theta_gradient_c),
                                 Rcpp::_["method"] = "BFGS",
                                 // Pass in the other parameters as everything
                                 // is scoped environmentally
                                 Rcpp::_["cluster"] = cluster,
                                 Rcpp::_["X"] = X,
                                 Rcpp::_["Z"] = Z,
                                 Rcpp::_["time"] = time,
                                 Rcpp::_["y"] = y,
                                 Rcpp::_["yobs"] = yobs,
                                 Rcpp::_["y1on0"] = y1on0,
                                 Rcpp::_["betaX"] = betaX,
                                 Rcpp::_["betat"] = betat,
                                 Rcpp::_["sigma"] = sigma,
                                 Rcpp::_["gamma1"] = gamma1,
                                 Rcpp::_["gamma2"] = gamma2,
                                 Rcpp::_["gamma3"] = gamma3,
                                 Rcpp::_["alpha0"] = alpha0,
                                 Rcpp::_["alpha"] = alpha,
                                 Rcpp::_["rho"] = rho);
  
  // Extract out the estimated parameter values
  arma::vec out = Rcpp::as<arma::vec>(opt_results[0]);
  
  arma::mat I = pred_theta_V_c(out,  cluster,
                          X,Z,time,y,yobs, y1on0,
                          betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho);
  res[0] = out;
  res[1] = I;
  return res;
}

// [[Rcpp::export]]
arma::cube pred_y_persondraw(   const arma::uvec& idx,
                                int Ni,
                                const arma::mat& X,
                                const arma::mat& Z,
                                const arma::vec& time,
                                const arma::mat& y,
                                const arma::mat& yobs,
                                const arma::vec& y1on0,
                                const arma::mat& betaX, 
                                const arma::vec& betat, 
                                const arma::mat& sigma,
                                const arma::vec& gamma1,
                                const arma::vec& gamma2,
                                const arma::vec& gamma3,
                                double alpha0,
                                const arma::vec& alpha,
                                const arma::vec& rho){
 // conditional on the first Ni observations (index being idx.subvec(0,Ni-1))
 // predict on the times at index idx.subvec(Ni, idx.n_elem-1)
  int pN = idx.n_elem; //time.n_elem

  Rcpp::List opt(2);
  int K = y.n_cols;
  
  arma::vec bm = zeros(K);
  
  arma::vec di = zeros<vec>(pN);
  arma::mat mui(pN, K);
  arma::cube ypost(2,pN,K);
  
  for(int cluster = 1; cluster <=2; cluster++){
    opt = pred_optsol(cluster, X.rows(idx(0), idx(Ni-1)), 
                 Z.rows(idx(0), idx(Ni-1)),
                 time.subvec(idx(0), idx(pN-1)),
                 y.rows(idx(0), idx(Ni-1)), 
                 yobs.rows(idx(0), idx(Ni-1)), 
                 y1on0.subvec(idx(0), idx(pN-1)),
                 betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho);
    arma::mat V = symmat(opt[1]);
    //arma::mat thetasamp = mvrnormArma(1, opt[0], V); // Ni x 1
    arma::mat thetasamp = zeros(pN, 1);
    if(Ni == (pN-1)){
      Rcpp::List optsub = optsol(cluster, X.rows(idx(0), idx(Ni-1)), 
                              Z.rows(idx(0), idx(Ni-1)),
                              time.subvec(idx(0), idx(Ni-1)),
                              y.rows(idx(0), idx(Ni-1)), 
                              yobs.rows(idx(0), idx(Ni-1)),
                              y1on0.subvec(idx(0), idx(Ni-1)),
                              betaX, betat, sigma, gamma1, gamma2, gamma3, alpha0, alpha, rho);
      arma::mat V1 = symmat(optsub[1]);
      arma::mat X1 = mvrnormArma(1, optsub[0], V1);

      int l1 = Ni; 
      arma::mat S11 = V.submat(0, 0, l1-1, l1-1);
      arma::mat S12 = V.submat(0, l1, l1-1, V.n_cols-1);
      arma::mat S21 = V.submat(l1, 0, V.n_cols-1, l1-1);
      arma::vec mu = opt[0];

      arma::vec d1 = X1.col(0) - mu.subvec(0,Ni-1);
      arma::vec mu2_1 = mu(Ni, 0) + S21 * inv(S11) * d1;
      arma::mat S2_1  = V(Ni, Ni) - S21 * inv(S11) * S12;
      arma::mat X2 = mvrnormArma(1, mu2_1, S2_1);

      thetasamp.submat(0,0,Ni-1, 0) = X1;
      thetasamp(Ni, 0) = X2(0,0);
    } else {
      thetasamp = mvrnormArma(1, opt[0], V);
      
    }
    
    //Rcout  << "cluster "<< cluster << " step 2 \n";
    //Rcpp::print(opt[0]);  
    //Rcout  << "\n";
    
  if(cluster == 1){
    di = Z.rows(idx(0), idx(pN-1)) * alpha + thetasamp.col(0);
  } else {
    di = - alpha0 + Z.rows(idx(0), idx(pN-1)) * alpha + thetasamp.col(0);
  }
  mui =  X.rows(idx(0), idx(pN-1)) * betaX.t() + time.subvec(idx(0), idx(pN-1)) * betat.t() ;
 
  for(int j = 0; j < K; j ++){
    mui.col(j) += gamma1(j)/(1+exp(-gamma2(j)*(di-gamma3(j))));
  }
  
  ypost( span(cluster-1), span(0, pN-1), span(0, K-1) ) = mui + mvrnormArma(pN, bm, arma::diagmat(sigma.row(cluster-1))).t(); // ypost.row(cluster-1)  = t(K x pN)
  
  }//cluster
  
  return ypost;
}

// [[Rcpp::export]]
Rcpp::List getall_yClst(   const arma::uvec& allidx,
                           const arma::vec& eachct,
                           int smax,
                           const arma::mat& X,
                           const arma::mat& Z,
                           const arma::vec& time,
                           const arma::mat& y,
                           const arma::mat& yobs,
                           const arma::vec& y1on0,
                           const arma::cube& betaX, 
                           const arma::mat& betat, 
                           const arma::cube& sigma,
                           const arma::vec& gamma1,
                           const arma::mat& gamma2,
                           const arma::mat& gamma3,
                           const arma::vec& alpha0,
                           const arma::mat& alpha,
                           const arma::mat& rho){
  Rcpp::List output(2); // two clusters
  int K = y.n_cols;
  arma::cube y1(smax, allidx.n_elem, K);
  arma::cube y2(smax, allidx.n_elem, K);
  int ct = 0;
  int ctl;
  for(uword i = 0; i < eachct.n_elem; i++){
    ctl = ct + eachct(i);
    arma::uvec idx(eachct(i)) ;
    idx = allidx.subvec(ct, ctl-1 );
    
    arma::cube ypost(2,eachct(i),K);
    int hl = eachct(i)-1;//only predict the last observation
    
    for(int s = 0; s < smax; s ++){
      Rcout  << "i " << i << " s "<< s<<"\n";
      Rcout << "idx " << idx.front() << " & " << idx.back() << "\n";
      ypost = pred_y_persondraw(idx-1, hl, X, 
                                Z,
                                time,
                                y, 
                                yobs, 
                                y1on0,
                                betaX.row(s), betat.row(s).t(), sigma.row(s), 
                                gamma1, gamma2.row(s).t(), gamma3.row(s).t(), 
                                alpha0(s), alpha.row(s).t(), rho.row(s).t()  );
      /*ypost = pred_y_persondraw(idx, hl, X.rows(idx(0), idx(hl-1)), 
                     Z.rows(idx(0), idx(hl-1)),
                     time.subvec(idx.front(), idx.back()),
                     y.rows(idx(0), idx(hl-1)), 
                     yobs.rows(idx(0), idx(hl-1)), 
                     K,
                     betaX.row(s), betat.row(s).t(), sigma.row(s), 
                     gamma1, gamma2.row(s).t(), gamma3.row(s).t(), 
                     alpha0(s), alpha.row(s).t(), rho.row(s).t()  );*/
    
      y1( span(s), span(ct, ctl-1), span(0, K-1) ) = ypost.row(0);
      y2( span(s), span(ct, ctl-1), span(0, K-1) ) = ypost.row(1);
    }
    ct = ctl;
  }
  output[0] = y1;
  output[1] = y2;
  return output;
}

