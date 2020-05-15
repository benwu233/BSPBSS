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
mcmc_bspbss = function(X, init,prior,kernel,ep = 0.5,lr = 0.1,subsample_n,subsample_p,n.iter,n.burn_in,thin=1,show_step){

  sigma0 = apply(X - init$A%*%init$S,2,var)

  A = init$A
  b = init$b
  zeta = init$zeta
  tau = init$tau


  lr0 = lr / nrow(X) / ncol(X)

  out = mcmc_bspbss_c(X, A, b, tau, sigma0, zeta, subsample_n, subsample_p, prior, kernel$psi, kernel$lambda,
                          ep,lr0, n.iter,n.burn_in,thin, show_step)

  return(out)
}
