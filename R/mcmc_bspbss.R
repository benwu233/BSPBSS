#' @title MCMC algorithm for Bayesian spatial blind source separation
#' with the thresholded Gaussian Process prior.
#' @description Performan MCMC algorithm to draw samples from a Bayesian spatial blind source separation
#' model.
#'
#' @param X Data matrix with n rows (sample) and p columns (voxel).
#' @param init List of initial values, see \code{init_bspbss}.
#' @param prior List of priors, see \code{init_bspbss}.
#' @param kernel List including eigenvalues and eigenfunctions of the kernel, see \code{init_bspbss}.
#' @param n.iter Total iterations in MCMC.
#' @param n.burn_in Number of burn-in.
#' @param thin Thining interval.
#' @param show_step Frequency for printing the current number of iterations.
#' @param ep Approximation parameter.
#' @param lr Per-batch learning rate in SGHMC.
#' @param decay Decay parameter in SGHMC.
#' @param subsample_n Mini-batch size of samples.
#' @param subsample_p Mini-batch size of voxels.
#'
#' @return List that contains MCMC samples of: A, b, sigma, and zeta.
#' @export
#'
#' @exampless
mcmc_bspbss = function(X,init,prior,kernel,n.iter,n.burn_in,thin=1,show_step,
                       ep = 0.01,lr = 0.01, decay = 0.01,
                       subsample_n = 0.5,subsample_p=0.5){

  sigma = init$sigma * 1
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

  out = mcmc_bspbss_c(X, A, b, sigma, zeta, stepsize_zeta, subsample_n, subsample_p, prior, kernel$psi, kernel$lambda,
                    ep,lr0, decay, n.iter,n.burn_in,thin, show_step)

  return(out)
}

