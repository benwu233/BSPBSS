#' @title Initial values
#' @description Generate initial values, set up priors and perform kernel decomposition
#' for runing MCMC algorithm.
#'
#' @param X Data matrix with n rows (sample) and p columns (voxel).
#' @param xgrid Cordinate matrix with p0 rows (voxel) and d columns (dimension). p0 >= p.
#' @param mask Logical vector of length p0, with p elements being equal to 1 and p0-p elements being equal to 0.
#' @param q Number of latent sources.
#' @param dens The initial density level (between 0 and 1) of the latent sources.
#' @param ker_par 2-dimensional vector (a,b) with a>0, b>0, specifing the parameters in the modified exponetial squared kernel.
#' @param num_eigen Number of the eigen functions.
#' @param noise Gaussian noise added to the initial latent sources, with mean 0 and standard deviation being noise * sd(S0),
#' where sd(S0) is the standard deviation of the initial latent sources.
#'
#' @return List includes the initial values, priors and eigen functions/eigen values of the kernel.
#' @export
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats sd
#' @importFrom glmnet glmnet
#' @import movMF
#' @importFrom stats quantile
#' @importFrom BayesGPfit GP.std.grids
#' @importFrom BayesGPfit GP.eigen.value
#' @importFrom BayesGPfit GP.eigen.funcs.fast
#' @useDynLib BSPBSS
#'
#' @examples
#'
init_bspbss= function(X, xgrid,  mask = rep(1,nrow(xgrid)), q = 2, dens = 0.5, ker_par = c(0.05, 20), num_eigen = 500, noise = 0.0 ){

  dim = ncol(xgrid)
  n = nrow(X)
  p = ncol(X)

  xgrid0 = GP.std.grids(xgrid,center=apply(xgrid,2,mean),scale=NULL,max_range=1)

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

  Psi0 = t( GP.eigen.funcs.fast(xgrid0[mask==1,],tag,ker_par[1],ker_par[2] ) )
  Psi_tmp = Psi0[1:L,]

  ica_tmp = ICA_imax(X,q)
  A0 = ica_tmp$A * sqrt(n)
  S00 = ica_tmp$S * sqrt(1/n)
  sdS00 = sd(S00)
  S0 = S00 + matrix( rnorm( q*p, mean = 0, sd = noise*sdS00 ) , nrow = q, ncol =p )

  b0 = matrix(0, nrow = q, ncol = L)
  S1 = matrix(0, nrow = q, ncol = p)


  sdS0 = rep(0,q)
  for(j in 1:q){

    sdS0[j] =sd(S0[j,])
    S1[j,] = S0[j,]/sdS0[j]

    tmp = glmnet(x = t(Psi_tmp), y = S1[j,], family = "gaussian", alpha = 0.5, nlambda = 2,intercept = FALSE)

    b0[j,] = tmp$beta[,length(tmp$lambda)]

  }


  tau0 = mean(sdS0)
  b0 = b0

  sumb0 = cal_sumb(b0,Psi_tmp)
  zeta0 = quantile(abs(sumb0),1-dens)


  S2 = cal_S(A0,sumb0,zeta0,tau0)
  init = list()

  init$A = A0
  init$ICA_S = S0
  init$S = S2
  init$zeta = zeta0
  init$b = b0
  init$tau = tau0

  kernel = list()
  kernel$psi = Psi_tmp
  kernel$lambda = lambda_tmp

  prior = list()
  prior$sigma = c(1,1e-10)

  prior$A = rep(0,n+1)

  prior$zeta = c(zeta0,zeta0^2)
  prior$tau = c(tau0,tau0^2)

  prior$zeta = c(0,1e20)
  prior$tau = c(1,1e20)

  out = list()
  out$init = init
  out$prior = prior
  out$kernel = kernel


  out$grid = xgrid0

  return(out)
}
