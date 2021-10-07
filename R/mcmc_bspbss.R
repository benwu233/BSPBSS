#' @title MCMC algorithm for Bayesian spatial blind source separation
#' with the thresholded Gaussian Process prior.
#' @description Performan MCMC algorithm to draw samples from a Bayesian spatial blind source separation
#' model.
#'
#' @param X Data matrix with n rows (sample) and p columns (voxel).
#' @param init List of initial values, see \code{init_mcmc}.
#' @param prior List of priors, see \code{init_mcmc}.
#' @param kernel List including eigenvalues and eigenfunctions of the kernel, see \code{init_TGP}.
#' @param n.iter Total iterations in MCMC.
#' @param n.burn_in Number of burn-in.
#' @param thin Thining interval.
#' @param show_step Frequency for printing the current number of iterations.
#' @param ep Approximation parameter.
#' @param lr Per-batch learning rate.
#' @param subsample_n Size of a mini-batch for samples.
#' @param subsample_p Size of a mini-batch for voxels.
#'
#' @return List that contains MCMC samples for: A, b, sigma, zeta and tau.
#' @export
#'
#' @exampless
mcmc_bspbss = function(X,init,prior,kernel,n.iter,n.burn_in,thin=1,show_step,
                       ep = 0.01,lr = 0.01, decay = 0.01,
                       subsample_n = 0.5,subsample_p=0.5){

  sigma0 = apply(X - init$A%*%init$S,2,var)

  A = init$A * 1
  b = init$b * 1
  zeta =  init$zeta * 1
  stepsize_zeta = init$stepsize_zeta

  if(prior$zeta[1] < 0){
    prior$zeta[1] = 0
  }

  if(prior$zeta[2] > 1){
    prior$zeta[2] = 1
  }


  lr0 = lr / nrow(X) / ncol(X)

  out = mcmc_bspbss_c(X, A, b, sigma0, zeta, stepsize_zeta, subsample_n, subsample_p, prior, kernel$psi, kernel$lambda,
                    ep,lr0, decay, n.iter,n.burn_in,thin, show_step)

  return(out)
}



#' @title MCMC algorithm for Bayesian spatial blind source separation
#' with the thresholded Gaussian Process prior.
#' @description Performan MCMC algorithm to draw samples from a Bayesian spatial blind source separation
#' model.
#'
#' @param X Data matrix with n rows (sample) and p columns (voxel).
#' @param init List of initial values, see \code{init_mcmc}.
#' @param prior List of priors, see \code{init_mcmc}.
#' @param kernel List including eigenvalues and eigenfunctions of the kernel, see \code{init_TGP}.
#' @param ep Approximation parameter.
#' @param lr Per-batch learning rate.
#' @param subsample_n Size of a mini-batch for samples.
#' @param subsample_p Size of a mini-batch for voxels.
#' @param n.iter Total iterations in MCMC.
#' @param n.burn_in Number of burn-in.
#' @param thin Thining interval.
#' @param show_step Frequency for printing the current number of iterations.
#'
#' @return List that contains MCMC samples for: A, b, sigma, zeta and tau.
#' @export
#'
#' @exampless
mcmc_bspbss_pickup = function(X, res, init,prior,kernel,ep = 0.5,lr = 0.1, decay = 0.1, subsample_n,subsample_p,n.iter,thin=1,show_step){

  n0 = ncol(res$sigma)

  sigma0 = res$sigma[,n0] * 1

  A = res$A[,,n0] * 1
  b = res$b[,,n0] * 1
  zeta =  res$zeta[n0] * 1
  stepsize = res$zeta_stepsize

  lr0 = lr / nrow(X) / ncol(X)

  out = mcmc_bspbss_pickup_c(X, A, b, sigma0, zeta, stepsize, subsample_n, subsample_p, prior, kernel$psi, kernel$lambda,
                      ep,lr0, decay, n.iter,thin, show_step)

  return(out)
}

