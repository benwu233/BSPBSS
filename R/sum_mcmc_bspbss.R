#' @title Summarization of the MCMC result.
#' @description The function summarizes the MCMC results obtained from \code{mcmc_bspbss}.
#'
#' @param res List including MCMC samples, which can be obtained from function \code{mcmc_bspbss}
#' @param X Original data matrix.
#' @param kernel List including eigenvalues and eigenfunctions of the kernel, see \code{init_bspbss}.
#' @param start Start point of the iterations being summarized.
#' @param end End point of the iterations being summarized.
#' @param select_prob Lower bound of the posterior inclusion probability required when summarizing
#' the samples of latent sources.
#'
#' @return List that contains the following terms:
#' \describe{
#'   \item{S}{Estimated latent sources.}
#'   \item{pip}{Voxel-wise posterior inclusion probability for the latent sources.}
#'   \item{A}{Estimated mixing coefficent matrix.}
#'   \item{zeta}{Estimated zeta.}
#'   \item{sigma}{Estimated sigma.}
#'   \item{logLik}{Trace of log-likelihood.}
#'   \item{Slist}{MCMC samples of S.}
#' }
#' @export
#'
#' @examples
#'
#' sim = sim_2Dimage(length = 30, sigma = 5e-4, n = 30, smooth = 6)
#' ini = init_bspbss(sim$X, sim$coords, q = 3, ker_par = c(0.1,50), num_eigen = 50)
#' res = mcmc_bspbss(ini$X,ini$init,ini$prior,ini$kernel,n.iter=200,n.burn_in=100,thin=10,show_step=50)
#' res_sum = sum_mcmc_bspbss(res, ini$X, ini$kernel, start = 11, end = 20, select_p = 0.5)
#'
sum_mcmc_bspbss = function(res, X, kernel, start = 1, end = 100, select_prob = 0.8){

  out = list()

  tmp = mcmc_sum_avgS(res, X, kernel, start, end, select_prob)

  out$S = tmp$S
  out$pip = tmp$pip
  out$A = apply(res$A,c(1,2),mean)
  out$zeta = mean(res$zeta)
  out$sigma = apply(res$sigma,1,mean)
  out$loglik = tmp$loglik
  out$Slist = tmp$Slist

  return(out)
}



mcmc_sum_avgS = function(res, data, kernel, start = 1, end = 1, select_prob = 0.8){

  n = end - start + 1

  q = dim(res$A)[2]
  p = ncol(kernel$psi)

  loglik = rep(0,n)

  Sl = matrix(0,nrow=q,ncol=p)
  spMat = matrix(0,nrow=q,ncol=p)

  Slist=list()

  for(i in start:end){

    sumb = cal_sumb(res$b[,,i],kernel$psi)
    S = cal_S(sumb,res$zeta[i])
    Slist[[i-start + 1]] = S

    loglik[i-start + 1] = loglk(data,res$A[,,i],S,res$sigma[,i])

    Sl = Sl + S
    spMat = spMat + (S!=0)

  }

  Sl = Sl * ( (spMat / n) > select_prob)

  Sl[spMat > 0] = Sl[spMat > 0]/spMat[spMat > 0]

  out = list()
  out$S = Sl
  out$pip = spMat/n
  out$loglik = loglik
  out$Slist =Slist

  return( out )
}



