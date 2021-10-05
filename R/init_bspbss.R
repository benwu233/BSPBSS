#' @title Initial values
#' @description Generate initial values, set up priors and perform kernel decomposition
#' for the MCMC algorithm.
#'
#' @param X Data matrix with n rows (sample) and p columns (voxel).
#' @param standardize If TRUE, standarize each row of X.
#' @param xgrid Cordinate matrix with p0 rows (voxel) and d columns (dimension). p0 >= p.
#' @param mask Logical vector of length p0, with p elements being equal to 1 and p0-p elements being equal to 0.
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
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom stats sd
#' @importFrom glmnet glmnet
#' @importFrom svd propack.svd
#' @import movMF
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
#' @examples
#'
init_bspbss= function(X, standardize = TRUE, xgrid,  mask = rep(1,nrow(xgrid)), q = 2, dens = 0.5, kernel="gaussian",ker_par = c(0.05, 20), num_eigen = 500, noise = 0.0 ){

  dim = ncol(xgrid)
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

  xgrid0 = GP.std.grids(xgrid,center=apply(xgrid,2,mean),scale=NULL,max_range=1)

  if(kernel=="gaussian"){

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
  }
  else if(kernel=="matern"){

    if(num_eigen>p){
      num_eigen = p
    }

    xgrid1 = xgrid0[mask==1,]

    n_cluster = num_eigen
    if(n_cluster > p){
      n_cluster = p
    }

    km = kmeans(xgrid1, centers = n_cluster,iter.max = 100)
    ind0 = rep(0,n_cluster)
    for(l in 1:n_cluster){
      ind1 = which(km$cluster==l)
      tmp = abs( xgrid1[km$cluster==l,] - km$centers[l,] )
      ind0[l] = ind1[which.min( apply(tmp,1,sum) )]
    }

    ind0 = floor(seq(1,p,length.out = n_cluster) )

    distMat0 = disM(xgrid1,ind0-1)
    #distMat0 = disM(xgrid1,0:(p-1))
    distMat0 = sqrt(distMat0) #L*p

    rho0 = ker_par[1]
    kappa0 = ker_par[2] #3/2,5/2

    #distMat0 = disM_full(xgrid1,rho0)

    distMat = apply(distMat0/rho0, 2, matern, kappa0)


    distMat_sub = distMat[,ind0]
    eig0 = eigen(distMat_sub)

    Psi0 = sqrt(n_cluster) * t(t(eig0$vectors) %*% distMat / eig0$values)

    Psi_tmp = t(Psi0[,1:num_eigen])
    lambda_tmp = eig0$values[1:num_eigen]/n_cluster

  }

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
  #init$zeta_mu = mean(zeta0)
  init$b = b0

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
  out$grid = xgrid0
  out$X = X

  return(out)
}


fn3 = function(para,sampcov,dist0){
  rho = para[1]
  nu = para[2]
  sigma = para[3]
  n = length(sampcov)
  out = 0
  for(i in 1:n){
    out = out + (sampcov[i] -  matern(dist0[i]/rho,nu))^2
  }
  return(out)
}



