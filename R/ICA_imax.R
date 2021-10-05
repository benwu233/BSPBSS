#' @importFrom ica icaimax

ICA_imax = function(X,q){
  ica0 = icaimax(t(X),nc=q,center=FALSE)

  IC_initial = ica0$S
  A_mt =  ica0$M

  D = diag(apply(A_mt,2,'norm_vec'))

  A_0 = A_mt%*%solve(D)
  S_0 = t(IC_initial%*%t(D))

  out = list()
  out$S = S_0
  out$A = A_0

  return(out)
}

norm_vec = function(v){
  return(sqrt(sum(v^2)))
}
