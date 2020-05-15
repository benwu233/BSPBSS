#' @title Simulate image data using ICA
#' @description The function simulates image data using a probabilistic ICA model
#' whose latent components have specific spatial patterns.
#'
#' @param length The length of the image.
#' @param n sample size.
#' @param sigma variance of the noise.
#' @param smooth smoothness of the latent components.
#'
#' @details
#' The observations are generated using probabilistic ICA:
#' \deqn{  X_i(v) = \sum_{j=1}^q A_{i,j} S_j(v) + \epsilon_i(v) , }
#' where \eqn{S_j, j=1,...,q} are the latent components, \eqn{A_{i,j}} is
#' the mixing coeffecient and \eqn{\epsilon_i} is the noise term.
#' Specifically, the number of components in this function is \eqn{q = 3},
#' with each of them being a specific geometric shape. The mixing coefficient matrix
#' is generated with a von Mises-Fisher distribution with the concentration parameter
#' being zero, which means it is uniformly distributed on the sphere. \eqn{\epsilon_i}
#' is a i.i.d. Gaussian noise term with 0 mean and user-specified variance.
#'
#' @return List that contains the following terms:
#' \describe{
#'   \item{X}{Data matrix with n rows (sample) and p columns (pixel).}
#'   \item{xgrid}{Cordinate matrix with p rows (pixel) and d columns (dimension)}
#'   \item{S}{Latent components.}
#'   \item{A}{Mixing coefficent matrix.}
#'   \item{snr}{Signal-to-noise ratio.}
#' }
#'
#' @export
#'
#' @importFrom stats cov rnorm var
#' @importFrom rstiefel rmf.vector
#'
#' @examples
#'
sim_2Dimage_ICA = function(length = 20,  n = 50, sigma = 2e-3, smooth = 6){


  grid = 1:length
  xgrid = as.matrix(expand.grid(grid,grid))
  voxel = xgrid / max(xgrid)
  p = nrow(xgrid)

  q = 3

  S0 = matrix(0,nrow = q, ncol = p)
  S0[1,] = create.square.in.2D.image(voxel, 0.1,0.1,0.5,0.5)
  S0[2,] = create.circle.in.2D.image(voxel, center = c(0.7,0.7))
  S0[3,] = create.triangle.in.2D.image(voxel, 0.2,0.8,0.2,0.2,0.8,0.5)

  if(smooth>0){
    S = smoos(S0,xgrid,smooth) /sqrt(p)
  }

  else{
    S = S0/sqrt(p)
  }

  A = matrix(0,nrow = n, ncol = q)

  for(j in 1:q){
    A[,j] = rmf.vector(rep(0,n)) * sqrt(n)
  }

  sigma2 = rep(sigma,p)
  ep = matrix(0,nrow = n, ncol = p)

  for(i in 1:p){
    ep[,i] = rnorm(n)
    ep[,i] = ep[,i] * sqrt(sigma2[i])
  }
  AS = A%*%S
  X = AS + ep

  data0 = array(0,dim = c( grid[length(grid)], grid[length(grid)], nrow(X)) )

  tag = 1
  for(i in 1:nrow(data0)){
    for(j in 1:ncol(data0)){
      data0[i,j,] = X[,tag]
      tag = tag + 1
    }
  }

  tmpp = 0
  for(i in 1:n){
    tmpp = tmpp + var(AS[i,])/( var(X[i,]) - var(AS[i,]) )
  }

  snr = tmpp/n

  out = list()

  out$X = X
  out$xgrid = xgrid
  out$A = A
  out$S = S
  out$snr = snr

  return(out)
}
