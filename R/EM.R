#' Expectation step
EM.E <- function(Xc, W, sigSq) {
  if (all(crossprod(W)==0)) {stop(paste("StPCA has failed due to numerical",
    "instability. Try dividing X by it's largest singular value."))}
  M = Matrix(crossprod(W) + sigSq*diag(ncol(W)))
  V = unname(as.matrix(t(solve(M, t(Xc %*% W)))))
  Vvar1 = as.matrix(sigSq*solve(M))
  Vvar = lapply(1:nrow(Xc), function(i_) {
    unname(Vvar1 + tcrossprod(V[i_,]))
  })
  return(list(V=V, Vvar=Vvar))
}

#' Maximization step for sigma^2
EM.M.sigSq <- function(Xc, W, V, Vvar) {
  sigSqNew = (
    norm(Xc, 'F')^2 -
    2*sum(vapply(1:nrow(Xc), function(n_) V[n_,] %*% crossprod(W, Xc[n_,]), numeric(1))) +
    sum(vapply(1:ncol(Xc), function(d_) W[d_,] %*% Reduce('+', Vvar) %*% W[d_,], numeric(1)))
  )/(nrow(Xc)*ncol(Xc))
  return(max(0, sigSqNew))
}

#' Maximization step for W
EM.M.W <- function(Xc, sigSq, V, Vvar, K) {
  xvsum = Matrix(Reduce('+', lapply(1:nrow(Xc), function(i_) tcrossprod(Xc[i_,], V[i_,]))))
  Vvarsum = Matrix(Reduce('+', Vvar))

  A = K/sigSq
  B = solve(Vvarsum)
  C = A %*% t(solve(Vvarsum, t(xvsum)))

  return(sylSolve(A, B, C))
}
