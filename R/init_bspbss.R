#' @title Initial values
#' @description Generate initial values, set up priors and perform kernel decomposition
#' for runing MCMC algorithm.
#'
#' @param X Data matrix with n rows (sample) and p columns (voxel).
#' @param xgrid Cordinate matrix with p0 rows (voxel) and d columns (dimension). p0 >= p.
#' @param mask Logical vector of length p0, with p elements being equal to 1 and p0-p elements being equal to 0.
#' @param q Number of latent sources.
#' @param dens The initial density level (between 0 and 1) of the latent sources.
#' @param kernel "gaussian" or "matern".
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
#' @importFrom svd propack.svd
#' @import movMF
#' @importFrom stats quantile
#' @importFrom BayesGPfit GP.std.grids
#' @importFrom BayesGPfit GP.eigen.value
#' @importFrom BayesGPfit GP.eigen.funcs.fast
#' @importFrom RandomFieldsUtils matern
#' @importFrom stats kmeans
#' @importFrom stats optim
#' @useDynLib BSPBSS
#'
#' @examples
#'
init_bspbss= function(X, xgrid,  mask = rep(1,nrow(xgrid)), q = 2, dens = 0.5, kernel="gaussian",ker_par = c(0.05, 20), num_eigen = 500, noise = 0.0 ){

  dim = ncol(xgrid)
  n = nrow(X)
  p = ncol(X)

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

    L = num_eigen
    if(L>p){
      L = p
    }

    xgrid1 = xgrid0[mask==1,]

    km = kmeans(xgrid1, centers = L,iter.max = 100)
    ind0 = rep(0,L)
    for(l in 1:L){
      ind1 = which(km$cluster==l)
      tmp = abs( xgrid1[km$cluster==l,] - km$centers[l,] )
      ind0[l] = ind1[which.min( apply(tmp,1,sum) )]
    }


    ind0 = floor(seq(1,p,length.out = L) )

    distMat0 = disM(xgrid1,ind0-1)
    #distMat0 = disM(xgrid1,0:(p-1))
    distMat0 = sqrt(distMat0) #L*p

    rho0 = ker_par[1]
    kappa0 = ker_par[2] #3/2,5/2

    #distMat0 = disM_full(xgrid1,rho0)


    # distMat = matrix(0,L,p)
    # for(l in 1:L){
    #   distMat[l,] = matern( distMat0[l,]/rho0, kappa0 )
    # }

    # distMat = matrix(0,p,p)
    # for(l in 1:p){
    #   distMat[l,] = matern( distMat0[l,]/rho0, kappa0 )
    # }

    distMat = apply(distMat0/rho0, 2, matern, kappa0)

    # distMat = matern( distMat0/rho0, kappa0 )
    # dim(distMat) = c(p,p)
    # print("yes")

    #distMat = distMat/distMat[1,1]

    distMat_sub = distMat[,ind0]
    eig0 = eigen(distMat_sub)
    #svd0 = svd(distMat_sub)


    # Psi0 = sqrt(L) * t(t(svd0$u) %*% distMat / svd0$d)
    # Psi_tmp = t(Psi0)
    # lambda_tmp = svd0$d/L

    Psi0 = sqrt(L) * t(t(eig0$vectors) %*% distMat / eig0$values)

    Psi_tmp = t(Psi0)
    lambda_tmp = eig0$values/L



    #svd0 = propack.svd(distMat,neig = L)
    #print("yes")

    #Psi_tmp = t(svd0$u)*sqrt(p)
    #lambda_tmp = svd0$d/p

    #lambda_tmp = lambda_tmp[1:L]
    #Psi_tmp = Psi_tmp[1:L,]
  }

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
  #init$distMat = distMat
  #init$ind = ind0

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



#' Title
#'
#' @param X
#' @param xgrid
#' @param mask
#' @param nu
#' @param thres
#'
#' @return
#' @export
#'
#' @examples
init_matern = function(X,xgrid,mask = rep(1,nrow(xgrid)),q=2,thres=0.75,
                       lower = c(0.01,0.01,0.01), upper = c(1e3,1e3,1e10),
                       initial_values = c(1.0, 1.0, 1.0)){

  xgrid0 = GP.std.grids(xgrid,center=apply(xgrid,2,mean),scale=NULL,max_range=1)
  xgrid1 = xgrid0[mask==1,]

  stdX = X - apply(X,1,mean)

  for(i in 1:nrow(X)){
    stdX[i,] = stdX[i,] / sd(stdX[i,])
  }

  maxX = apply(abs(stdX),2,max)
  th = quantile(maxX,thres[1])

  ind2 = (maxX>=th)

  para = find_Matern_optim_paras_alldat(t(X[,ind2]),xgrid1[ind2,], lower = lower,upper=upper,initial_values = initial_values)

  return( para)

}


init_matern0 = function(X,xgrid,mask = rep(1,nrow(xgrid)),q=2,thres=0.75){

  xgrid0 = GP.std.grids(xgrid,center=apply(xgrid,2,mean),scale=NULL,max_range=1)
  xgrid1 = xgrid0[mask==1,]

  maxX = apply(abs(X),2,max)
  th = quantile(maxX,thres[1])
  th2 = quantile(maxX,thres[2])


  ind2 = (maxX>=th)
  ind3 = (maxX<th2)
  stdX = X[,ind2] - apply(X[,ind2],1,mean)
  var0 = mean( apply(stdX,1,var) ) - mean( apply(X[,ind3],1,var) )

  #noise0 = var(as.numeric(X[,maxX<th2]) )/nrow(X)
  #var0 = (var(as.numeric(X) ) - var(as.numeric(X[,maxX<th2])) )/nrow(X)

  distMat1 = round( disM_full(xgrid1[ind2,]),10)
  distMat_tri = distMat1[upper.tri(distMat1,diag=TRUE)]
  uniqueM = unique(as.numeric(distMat_tri) )

  uniqueM_sorted = sort(uniqueM)

  ind1 = match(distMat_tri,uniqueM_sorted)

  sampcov = samp_cov(stdX,xgrid1[ind2,], ind1-1, length(uniqueM_sorted))

  #sampcov = sampcov/sampcov[uniqueM==0]
  #sampcov = sampcov/var0
  #sampcov = sampcov /3
  #sampcov = sampcov * var0

  thr = uniqueM_sorted[2]*25

  ind3 = which( (uniqueM_sorted<=thr)&(uniqueM_sorted>0) )

  para_est = optim(c(1,1,1),fn3,gr = NULL, sampcov[ind3],uniqueM_sorted[ind3],method = "L-BFGS-B", lower = c(1e-3,1e-3,1e-10),upper = c(1e2,1e2,1e10))

  return(para_est$par)

}


