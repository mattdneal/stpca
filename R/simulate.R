stpca.simulate.Wknown <- function(stpcaObj, n=1) {
  return(simulate.from.params(n, stpcaObj$W, stpcaObj$mu, stpcaObj$sigSq))
}

stpca.simulate.Wunknown <- function(stpcaObj, n=1) {
  W = t(rmvnorm(stpcaObj$k, sigma=as.matrix(stpcaObj$K)))
  return(simulate.from.params(n, W, stpcaObj$mu, stpcaObj$sigSq))
}

simulate.from.params <- function(n, W, mu, sigSq) {
  V = rmvnorm(n, sigma=diag(ncol(W)))
  X = sweep(tcrossprod(V, W), 2, mu, "+")
  return(list(V=V, X=X))
}
