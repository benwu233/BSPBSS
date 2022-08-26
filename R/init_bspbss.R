#' @title Initial values
#' @description Generate initial values, set up priors and perform kernel decomposition
#' for the MCMC algorithm.
#'
#' @param X Data matrix with n rows (sample) and p columns (voxel).
#' @param coords Cordinate matrix with p rows (voxel) and d columns (dimension).
#' @param standardize If TRUE, standarize each row of X.
#' @param q Number of latent sources.
#' @param dens The initial density level (between 0 and 1) of the latent sources.
#' @param ker_par 2-dimensional vector (a,b) with a>0, b>0, specifing the parameters in the modified exponetial squared kernel.
#' @param num_eigen Number of eigen functions.
#' @param noise Gaussian noise added to the initial latent sources, with mean 0 and standard deviation being noise * sd(S0),
#' where sd(S0) is the standard deviation of the initial latent sources.
#'
#' @return List includes the initial values, priors and eigen functions/eigen values of the kernel.
#' @export
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom stats sd
#' @importFrom glmnet glmnet
#' @importFrom svd propack.svd
#' @import movMF
#' @import gridExtra
#' @import gtools
#' @importFrom stats quantile
#' @importFrom BayesGPfit GP.std.grids
#' @importFrom BayesGPfit GP.eigen.value
#' @importFrom BayesGPfit GP.eigen.funcs.fast
#' @importFrom RandomFieldsUtils matern
#' @importFrom stats kmeans
#' @importFrom stats optim
#'
#'
#' @useDynLib BSPBSS
#'
init_bspbss= function(X, coords, standardize = TRUE, q = 2, dens = 0.5, ker_par = c(0.05, 20), num_eigen = 500, noise = 0.0 ){

  dim = ncol(coords)
  n = nrow(X)
  p = ncol(X)

  if(standardize){
    for(i in 1:n){
      sdx = sd(X[i,])
      if(sdx!=0){
        X[i,] = (X[i,]  )/sdx
      }
    }
  }

  coords0 = GP.std.grids(coords,center=apply(coords,2,mean),scale=NULL,max_range=1)

  lambda0 = GP.eigen.value(1000,ker_par[1],ker_par[2],dim)

  if(num_eigen < length(lambda0)){
    lambda_tmp = lambda0[1:num_eigen]
  }

  L = length(lambda_tmp)

  k = 0
  tag = 0

  while(k<=L){
    tag = tag + 1
    k = choose(tag+dim,dim)
  }

  Psi0 = t( GP.eigen.funcs.fast(coords0,tag,ker_par[1],ker_par[2] ) )
  Psi_tmp = Psi0[1:L,]


  ica_tmp = ICA_imax(X,q)
  A0 = ica_tmp$A * sqrt(n)
  S00 = ica_tmp$S * sqrt(1/n)
  sdS00 = sd(S00)
  S0 = S00 + matrix( rnorm( q*p, mean = 0, sd = noise*sdS00 ) , nrow = q, ncol =p )

  b0 = matrix(0, nrow = q, ncol = num_eigen)

  for(j in 1:q){
    tmp = glmnet(x = t(Psi_tmp), y = S0[j,], family = "gaussian", alpha = 0.5, nlambda = 2,intercept = FALSE)
    b0[j,] = tmp$beta[,length(tmp$lambda)]
  }

  sumb0 = cal_sumb(b0,Psi_tmp)

  zeta0 = quantile(abs(sumb0),1-dens)
  S1 = cal_S(sumb0,zeta0 - (1e-5) )

  init = list()
  init$A = A0
  init$ICA_S = S0
  init$S = S1
  init$zeta = zeta0+1e-10
  init$stepsize_zeta = (zeta0+1e-10) * 0.1
  init$b = b0
  init$sigma = apply(X - A0%*%S1,2,var)

  kernel = list()
  kernel$psi = Psi_tmp
  kernel$lambda = lambda_tmp

  prior = list()
  prior$sigma = c(1,1e-10)
  prior$A = c(1e-10, rep(1/sqrt(n),n) )
  prior$zeta = c(0,1)

  out = list()
  out$init = init
  out$prior = prior
  out$kernel = kernel
  out$coords = coords0
  out$X = X

  return(out)
}




