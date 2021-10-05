#include <stdio.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//[[Rcpp::export]]
arma::mat disM(arma::mat xgrid, arma::vec ind){

  int p = xgrid.n_rows;
  int d = xgrid.n_cols;
  int I = ind.n_elem;
  int ii = 0;
  double tmp = 0;
  arma::mat out(I,p);

  for(int i=0; i < I; i++){
    for(int j=0; j < p; j++){
      out(i,j) = 0;
    }
  }

  for(int i=0; i < I; i++){
    ii = ind[i];
    for(int j=0; j < p; j++){
      for(int k=0; k < d; k++){
        tmp = xgrid(ii,k) - xgrid(j,k);
        out(i,j) += tmp * tmp ;
      }
    }
  }

  return out;
}


//[[Rcpp::export]]
arma::mat disM_full(arma::mat xgrid){

  int p = xgrid.n_rows;
  int d = xgrid.n_cols;
  double tmp = 0;
  arma::mat out(p,p);
  out.zeros(p,p);

  for(int i=0; i < p; i++){
    for(int j=i; j < p; j++){
      for(int k=0; k < d; k++){
        tmp = xgrid(j,k) - xgrid(i,k);
        out(i,j) += tmp * tmp ;
      }
      out(j,i) = out(i,j);
    }
  }

  return out;
}

//[[Rcpp::export]]
arma::vec samp_cov(arma::mat X, arma::mat xgrid, arma::vec ind, int n){
  int n0 = ind.n_elem;

  int p = xgrid.n_rows;
  int m = X.n_rows;

  arma::vec out;
  arma::vec len0;
  arma::mat meanX1;
  arma::mat meanX2;
  int tag = 0;

  out.zeros(n);
  len0.zeros(n);

  tag = 0;
  for(int j = 0; j < p; j++){
    for(int i = 0; i <= j; i++){
      for(int k = 0; k < m; k++){
        out(ind[tag]) +=  X(k,i)  *  X(k,j) ;
      }
      len0[ind[tag]]++;
      tag++;
    }
  }

  for(int i = 0; i < n; i++){
    out[i] = out[i]/len0[i]/m;
  }

  return out;
}


//[[Rcpp::export]]
arma::vec samp_cov0(arma::mat X, arma::mat xgrid, arma::vec ind, int n){
  int n0 = ind.n_elem;

  int p = xgrid.n_rows;
  int m = X.n_rows;

  arma::vec out;
  arma::vec len0;
  arma::mat meanX1;
  arma::mat meanX2;

  out.zeros(n);
  len0.zeros(n);
  meanX1.zeros(m,n);
  meanX2.zeros(m,n);


  int tag = 0;
  for(int j = 0; j < p; j++){
    for(int i = 0; i <= j; i++){
      meanX1.col(ind[tag]) += X.col(i);
      meanX2.col(ind[tag]) += X.col(j);
      len0[ind[tag]]++;
      tag++;
    }
  }


  for(int i = 0; i < n; i++){
    meanX1.col(i) = meanX1.col(i)/len0[i];
    meanX2.col(i) = meanX2.col(i)/len0[i];
  }

  tag = 0;
  for(int j = 0; j < p; j++){
    for(int i = 0; i < j; i++){
      for(int k = 0; k < m; k++){
        out(ind[tag]) += ( X(k,i) - meanX1(k,ind[tag]) ) * ( X(k,j) - meanX2(k,ind[tag]));
      }
      tag++;
    }
  }

  for(int i = 0; i < n; i++){
    out[i] = out[i]/len0[i]/m;
  }

  return out;

}


//[[Rcpp::export]]
arma::mat disM_full0(arma::mat xgrid,double rho){

  int p = xgrid.n_rows;
  int d = xgrid.n_cols;
  double tmp = 0;
  arma::mat out(p,p);

  for(int i=0; i < p; i++){
    for(int j=i; j < p; j++){
      out(i,j) = 0;
      for(int k=0; k < d; k++){
        tmp = xgrid(i,k) - xgrid(j,k);
        out(i,j) += tmp * tmp ;
      }
      out(j,i) = out(i,j);
    }
  }

  return sqrt(out)/rho;
}


//[[Rcpp::export]]
arma::mat cal_sumb(arma::mat &b, arma::mat &psi){
  return  (b * psi);
}


//[[Rcpp::export]]
arma::mat cal_S(arma::mat &sumb, double zeta){

  int p = sumb.n_cols;
  int q = sumb.n_rows;

  arma::mat out = sumb;

  for(int j = 0; j < q; j++){
    for(int v = 0; v < p; v++){
      if( (sumb(j,v) < zeta) && (sumb(j,v) > -zeta) ){
        out(j,v) = 0;
      }
    }
  }

  return out;
}


//[[Rcpp::export]]
arma::mat cal_core(arma::mat &X, arma::mat &A, arma::mat &S){
  return (X - A*S);
}



//[[Rcpp::export]]
arma::mat dL_b_sub(arma::mat b, arma::mat X,arma::mat A,
                       arma::vec lambda, arma::mat psi, double epsilon,
                       double zeta, arma::vec sigma, int sizep, int sizen){

  int n = X.n_rows;
  int p = X.n_cols;
  int q = A.n_cols;
  int L = lambda.n_elem;
  arma::mat par(q,sizep);
  IntegerVector subp = Rcpp::sample(p,sizep,false) - 1;
  IntegerVector subn = Rcpp::sample(n,sizen,false) - 1;

  double tmp = 0;
  double tmp2 = 0;
  double tmp3 = 0;
  double ep2 = 0;

  tmp2 = epsilon / 3.1415926535897932;
  ep2 = epsilon * epsilon;

  arma::mat out2(q,L);
  arma::mat out3(q,L);

  arma::mat sumb_sub(q,sizep);
  arma::mat X_core_sub(sizen,sizep);
  arma::mat A_sub(sizen,q);

  arma::mat tmp_sub_a(q,sizep);
  arma::mat tmp_sub_b(q,sizep);
  arma::mat tmp_sub_c(q,sizep);

  arma::mat psi_sub(L,sizep);
  arma::vec sigma_sub(sizep);

  int loc = 0;
  int loc1 = 0;

  for(int v = 0; v < sizep; v++){
    loc = subp[v];
    for(int l = 0; l < L; l++){
      psi_sub(l,v) = psi(l,loc);
    }
    sigma_sub[v] = sigma[loc];
    for(int i = 0; i < sizen; i++){
      loc1 = subn[i];
      X_core_sub(i,v) = X(loc1,loc);
    }
  }

  for(int i = 0; i < sizen; i++){
    loc1 = subn[i];
    for(int j = 0; j < q; j++){
      A_sub(i,j) = A(loc1,j);
    }
  }

  sumb_sub = b * psi_sub;

  tmp_sub_a = ( sumb_sub - zeta ) % ( sumb_sub - zeta ) + ep2;
  tmp_sub_b = ( sumb_sub + zeta ) % ( sumb_sub + zeta ) + ep2;

  tmp_sub_c = tmp2 * (1 / tmp_sub_a - 1 / tmp_sub_b + 1e-20) % sumb_sub;

  int count1 = 0;
  int count2 = 0;
  for(int v = 0; v < sizep; v++){
    for(int j = 0; j < q; j++){
      if( (sumb_sub(j,v) > zeta)||(sumb_sub(j,v) < -zeta) ){
        for(int i = 0; i < sizen; i++){
          X_core_sub(i,v) -=  A_sub(i,j) * sumb_sub(j,v);
        }
        par(j,v) = ( tmp_sub_c(j,v)  + 1  ) / (sigma_sub[v]+1e-20);
      }
      else{
        par(j,v) = ( tmp_sub_c(j,v) ) /  (sigma_sub[v]+1e-20);
      }
    }
  }

  out2 =( (A_sub.t() * X_core_sub) % par ) * psi_sub.t();

  tmp = 1.0 * p*n / sizep/sizen;
  for(int j = 0; j < q; j++){
    for(int l = 0; l < L; l++){
      out2(j,l) = tmp * out2(j,l) - b(j,l) / lambda[l];
    }
  }

  return out2;
}



//[[Rcpp::export]]
void GP_update_b_SGHMC(arma::mat &b, arma::mat &X, arma::mat &A, arma::mat &S, arma::vec &sigma,
                       arma::vec &lambda,  arma::mat &psi, double &epsilon,
                       arma::mat &sumb, double &zeta, // arma::vec &zeta,
                       arma::mat &X_core, double eta, double alpha,
                       int sizep, int sizen, int m, int itr, arma::vec &nu){

  int q = A.n_cols;
  int L = lambda.n_elem;
  int loc = 0;

  double tmp2 = 0;

  arma::mat dl(q,L);

  if( (itr-1)%50==0) {
    nu = rnorm(q*L) * sqrt(eta);
  }

  arma::vec oldnu = nu;

  NumericVector rnorm0(q*L,0.0);

  while(m>0){
    m--;

    dl = dL_b_sub(b, X, A, lambda,  psi, epsilon, zeta, sigma, sizep, sizen);

    rnorm0 = rnorm(q*L) * sqrt(2.0* alpha * eta) ;

    for(int j = 0; j < q; j++){
      for(int l = 0; l < L; l++){
        loc = j+l*q;
        b(j,l) += nu[loc];
        nu[loc] += eta * dl(j,l) - alpha * nu[loc] + rnorm0[loc];
      }
    }

  }

  sumb = cal_sumb(b, psi);
  S = cal_S(sumb, zeta);
  X_core = cal_core(X, A, S);

}



//[[Rcpp::export]]
void GP_update_A(arma::mat &A, arma::vec prior, arma::mat &X_core, arma::mat X,
                 arma::mat S, arma::vec sigma){

  Rcpp::Environment movMFpack("package:movMF");
  Rcpp::Function rrmovMF = movMFpack["rmovMF"];

  int p = X_core.n_cols;
  int n = X_core.n_rows;
  int q = A.n_cols;

  arma::mat out(n,q);
  arma::vec newA_j(n);
  arma::mat par1(n,1);
  arma::mat sumXS(n,q);
  arma::mat sumSS(q,q);
  arma::mat sumAS(n,q);
  arma::mat Sdsigma(q,p);
  arma::mat tmp(n,q);
  NumericVector par(n,0.0);

  double kappa = prior[0];
  arma::vec rho(n);
  double sqrtn = sqrt(1.0*n);

  for(int i = 0 ; i < n; i++){
    rho[i] = prior[i+1] * kappa;
  }

  for(int j = 0; j < q; j++){
    Sdsigma.row(j) = S.row(j) / (sigma.t() + 1e-10);
  }


  sumXS = X  * Sdsigma.t();
  sumSS = S  * Sdsigma.t();

  for(int j = 0; j < q; j++){

    sumAS = A * sumSS;

    tmp = A.col(j) * sumSS.row(j);
    par1 = ( sumXS.col(j) - sumAS.col(j) + tmp.col(j) ) * sqrtn + rho + 1e-50;

    for( int i = 0; i < n; i++){
      par[i] = par1(i,0);
    }

    newA_j = as<arma::vec>( rrmovMF(1,par,1) ) * sqrtn;
    A.col(j) = newA_j;

  }

  X_core = cal_core(X,A,S);
}


//[[Rcpp::export]]
void GP_update_sigma(arma::vec &sigma, arma::mat &X_core, arma::vec prior_sigma){

  int p = X_core.n_cols;
  int n = X_core.n_rows;
  double shape = 0;
  double scale = 0;

  shape = 1.0*n/2 + prior_sigma[0];
  arma::mat X_core2 = X_core % X_core;

  for(int v = 0; v < p; v++){
    scale = sum(X_core2.col(v))/2 + prior_sigma[1];
    sigma[v] = 1 / R::rgamma(shape,1/scale);
  }

}


//[[Rcpp::export]]
double log_p_zeta_Gaussian(double zeta, arma::mat &X_core, arma::vec &sigma){

  double out = 0;

  out = -accu(  (X_core % X_core ) * (1/sigma) )/2;

  return out;
}


//[[Rcpp::export]]
void GP_update_zeta(double &zeta, arma::mat &sumb, arma::mat &X_core, arma::vec &sigma,
                    arma::mat &X, arma::mat &A, arma::mat &S, double stepsize, arma::vec prior,
                    double &count){

  arma::vec logpz(2);
  int utest = 0;
  double u = 0;
  double sd =  stepsize;

  int p = X_core.n_cols;
  int n = X_core.n_rows;
  int q = A.n_cols;

  arma::mat newS(q,p);
  arma::mat newX_core(n, p);
  arma::mat abssumb = abs(sumb);

  std::sort(abssumb.begin(),abssumb.end());

  double max0 = abssumb(q*p*prior[1]-1);

  NumericVector newzeta(1);

  if( (zeta < prior[0]) || (zeta > max0) ){
    newzeta = runif(1,prior[0],max0);
  }else{
    newzeta = zeta + rnorm(1) * sd;
  }

  if( (newzeta[0] < prior[0]) || (newzeta[0] > max0) ){
    utest = 0;
  }
  else{
    newS = cal_S(sumb, newzeta[0]);
    newX_core = cal_core(X,A,newS);

    logpz[0] = log_p_zeta_Gaussian(zeta,X_core,sigma);
    logpz[1] = log_p_zeta_Gaussian(newzeta[0],newX_core,sigma);

    u = R::runif(0,1);

    if(log(u) < (logpz[1] - logpz[0]) ){
      utest = 1;
    }
    else{
      utest = 0;
    }

  }

  if(utest ==1){
    zeta = newzeta[0];
    X_core = newX_core;
    S = newS;
    count++;
  }

}



//[[Rcpp::export]]
double loglk(arma::mat &X, arma::mat &A, arma::mat &S, arma::vec &sigma){
  double out = 0;
  int p = X.n_cols;
  int n = X.n_rows;

  arma::mat X_core(n, p);
  X_core = cal_core(X, A, S);

  for(int i = 0; i < n; i++){
    for(int v = 0; v < p; v++){
      out -= X_core(i,v) * X_core(i,v) / sigma[v];
    }
  }

  out = out / 2;

  for(int v = 0; v < p; v++){
    out -= n/2 * log(2*3.1415926*sigma[v]);
  }

  return out;
}


//[[Rcpp::export]]
List mcmc_bspbss_c(arma::mat &X, arma::mat &A, arma::mat &b,
                  arma::vec &sigma, double zeta, double stepsize_zeta, double subsample_n, double subsample_p,
                  List prior, arma::mat &psi, arma::vec &lambda, double epsilon, double lr, double decay,
                  int MClength, int burn_in, int thin, int show_step){


  int p = X.n_cols;
  int n = X.n_rows;
  int q = A.n_cols;
  int L = lambda.n_elem;
  int sizep = subsample_p * p;
  int sizen = subsample_n * n;

  auto t1 = std::chrono::system_clock::now();
  std::time_t t2 = std::chrono::system_clock::to_time_t(t1);

  double count_zeta = 0;

  int nchain = ( MClength )/thin;

  Rcpp::Dimension A_trace_dim(n,q,nchain);
  Rcpp::Dimension b_trace_dim(q,L,nchain);

  NumericVector A_trace(A_trace_dim);
  NumericVector b_trace(b_trace_dim);

  NumericMatrix sigma_trace(p,nchain);
  NumericVector tau_trace(nchain,0.0);
  NumericVector zeta_trace(nchain,0.0);
  NumericVector loglik_trace(nchain,0.0);

  arma::mat sumb(q,p);
  arma::mat S(q,p);
  arma::mat X_core(n,p);
  arma::vec nu(L*q);

  arma::vec prior_A = as<arma::vec>(  prior["A"] );
  arma::vec prior_sigma = as<arma::vec>( prior["sigma"]);
  arma::vec prior_zeta = as<arma::vec>( prior["zeta"]);

  int tag =  0;

  sumb = cal_sumb(b,psi);
  S = cal_S(sumb,zeta);
  X_core = cal_core(X, A, S);


  for(int itr = 1; itr <= MClength; itr++){

    GP_update_A(A, prior_A, X_core, X, S, sigma);

    GP_update_sigma(sigma, X_core, prior_sigma);

    GP_update_b_SGHMC(b, X, A,S, sigma, lambda, psi, epsilon, sumb, zeta, X_core,lr,decay,sizep,sizen,10,itr,nu);

    GP_update_zeta(zeta, sumb, X_core, sigma, X, A, S, stepsize_zeta, prior_zeta, count_zeta);


    if( ( (itr-burn_in)%thin==0 ) ){
      zeta_trace[tag] = zeta;

      for(int v = 0; v < p; v++){
        sigma_trace(v,tag) = sigma[v];
      }

      for(int i = 0; i < n; i++){
        for(int j = 0; j < q; j++){
          A_trace(i+n*j+tag*q*n) = A(i,j);
        }
      }

      for(int j = 0; j < q; j++){
        for(int l = 0; l < L; l++){
          b_trace(j+q*l+tag*q*L) = b(j,l);
        }
      }

      tag++;
    }

    if(itr%show_step==0){
      t1 = std::chrono::system_clock::now();
      t2 = std::chrono::system_clock::to_time_t(t1);
      std::cout << "iter " << itr << " " << std::ctime(&t2) << std::endl;
      std::cout << " stepsize_zeta " << stepsize_zeta << " count_zeta " << count_zeta << std::endl;
    }

    if(itr < burn_in){
      if(itr%100==0){
        if(count_zeta < 20){
          stepsize_zeta = stepsize_zeta *0.9;
        }
        if(count_zeta > 30){
          stepsize_zeta = stepsize_zeta *1.1;
        }
        count_zeta = 0;
      }
    }

  }

  return Rcpp::List::create(Named("A")=A_trace,
                            Named("b")=b_trace,
                            Named("sigma")=sigma_trace,
                            Named("zeta")=zeta_trace,
                            Named("zeta_stepsize")=stepsize_zeta);

}



//[[Rcpp::export]]
List mcmc_bspbss_pickup_c(arma::mat &X, arma::mat &A, arma::mat &b,
                   arma::vec &sigma, double zeta, double stepsize_zeta, double subsample_n, double subsample_p,
                   List prior, arma::mat &psi, arma::vec &lambda, double epsilon, double lr, double decay,
                   int MClength, int thin, int show_step){


  int p = X.n_cols;
  int n = X.n_rows;
  int q = A.n_cols;
  int L = lambda.n_elem;
  int sizep = subsample_p * p;
  int sizen = subsample_n * n;

  auto t1 = std::chrono::system_clock::now();
  std::time_t t2 = std::chrono::system_clock::to_time_t(t1);

  double count_zeta = 0;

  int nchain =  MClength/thin;

  Rcpp::Dimension A_trace_dim(n,q,nchain);
  Rcpp::Dimension b_trace_dim(q,L,nchain);

  NumericVector A_trace(A_trace_dim);
  NumericVector b_trace(b_trace_dim);

  NumericMatrix sigma_trace(p,nchain);
  NumericVector tau_trace(nchain,0.0);
  NumericVector zeta_trace(nchain,0.0);
  NumericVector loglik_trace(nchain,0.0);

  arma::mat sumb(q,p);
  arma::mat S(q,p);
  arma::mat X_core(n,p);
  arma::vec nu(L*q);

  arma::vec prior_A = as<arma::vec>(  prior["A"] );
  arma::vec prior_sigma = as<arma::vec>( prior["sigma"]);
  arma::vec prior_zeta = as<arma::vec>( prior["zeta"]);

  int tag =  0;

  sumb = cal_sumb(b,psi);
  S = cal_S(sumb,zeta);
  X_core = cal_core(X, A, S);


  for(int itr = 1; itr <= MClength; itr++){

    GP_update_A(A, prior_A, X_core, X, S, sigma);

    GP_update_sigma(sigma, X_core, prior_sigma);

    GP_update_b_SGHMC(b, X, A,S, sigma, lambda, psi, epsilon, sumb, zeta, X_core,lr,decay,sizep,sizen,10,itr,nu);

    GP_update_zeta(zeta, sumb, X_core, sigma, X, A, S, stepsize_zeta, prior_zeta, count_zeta);


    if( itr%thin==0 ){
      zeta_trace[tag] = zeta;

      for(int v = 0; v < p; v++){
        sigma_trace(v,tag) = sigma[v];
      }

      for(int i = 0; i < n; i++){
        for(int j = 0; j < q; j++){
          A_trace(i+n*j+tag*q*n) = A(i,j);
        }
      }

      for(int j = 0; j < q; j++){
        for(int l = 0; l < L; l++){
          b_trace(j+q*l+tag*q*L) = b(j,l);
        }
      }

      tag++;
    }

    if(itr%show_step==0){
      t1 = std::chrono::system_clock::now();
      t2 = std::chrono::system_clock::to_time_t(t1);
      std::cout << "iter " << itr << " " << std::ctime(&t2) << std::endl;
      std::cout << " stepsize_zeta " << stepsize_zeta << " count_zeta " << count_zeta << std::endl;
    }
  }

  return Rcpp::List::create(Named("A")=A_trace,
                            Named("b")=b_trace,
                            Named("sigma")=sigma_trace,
                            Named("zeta")=zeta_trace,
                            Named("zeta_stepsize")=stepsize_zeta);

}


//[[Rcpp::export]]
NumericMatrix smoos(NumericMatrix S, IntegerMatrix xgrid, double smooth){
  Rcpp::Dimension S_dim = S.attr("dim");
  int q = S_dim[0];
  int p = S_dim[1];
  NumericMatrix out(q,p);
  double count = 0;
  double tmp = 0;

  for(int j = 0; j < q; j++){
    for(int v = 0; v < p; v++){
      out(j,v) = 0;
    }
  }

  for(int j = 0; j < q; j++){
    for(int v0 = 0; v0 < p; v0++){
      count = 0;
      tmp = 0;
      for(int v = 0; v < p; v++){
        if( (xgrid(v0,1)-xgrid(v,1))*(xgrid(v0,1)-xgrid(v,1))
              +(xgrid(v0,0)-xgrid(v,0))*(xgrid(v0,0)-xgrid(v,0)) < smooth*smooth ) {
          count = count + 1;
          tmp = tmp + S(j,v);
        }
      }
      out(j,v0) = tmp / count;


    }
  }

  return out;
}
