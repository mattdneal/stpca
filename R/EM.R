#' Update theta to be the maximum-a-posteriori value using Expectation
#' Maximisation.
theta_EM <- function(X, W, mu, sigSq, K, maxit=500, bftol=1e-5) {
  stopifnot(bftol > 0)
  stopifnot(maxit > 0)

  k  = ncol(W)
  mu = colMeans(X)
  Xc = sweep(X, 2, mu)

  unNormedLPs = numeric(maxit) # un-normalised log posteriors
  for (iter in seq_len(maxit)) { # TODO: Iterate until zero gradient
    # Expectation step
    E = EM.E(Xc, W, sigSq)

    # Maximization step for W
    W = EM.M.W(Xc, sigSq, E$Vmean, E$Vvar, K)

    # Align W with PCs & fix sign
    W = W %*% svd(W)$v
    W = W %*% diag(sign(W[1,]), nrow=k, ncol=k)

    # Second expectation step using updated W
    E = EM.E(Xc, W, sigSq)

    # Maximization step for sigma^2
    sigSq = EM.M.sigSq(Xc, W, E$Vmean, E$Vvar)

    # ln p(X | \theta) + ln p(\theta | \beta)
    unNormedLPs[iter] = log_likelihood(Xc, W, 0, sigSq) + log_prior(K, W)

    # If the log Bayes' Factor ratio is less than bftol then stop early:
    # log(p(theta^n | X) / p(theta^{n-1} | X)) < bftol
    if (iter>1 && (unNormedLPs[iter]-unNormedLPs[iter-1])<bftol) break
  }

  # Truncate unassigned LPs
  if (iter < maxit) unNormedLPs = unNormedLPs[seq_len(iter)]

  # Assuming convergence, \sigma = \hat\sigma, so the evidence can be approximated
  H = compute_H(X, W, mu, sigSq, K)
  logEvidence = log_evidence(Xc, K, W, 0, sigSq, H)

  # We now know the evidence to normalize the above log posterior values
  LPs = unNormedLPs - logEvidence

  return(list(
    WHat = W,
    muHat = mu,
    sigSqHat = sigSq,
    Vmean = E$Vmean,
    Vvar = E$Vvar,
    logEvidence = logEvidence,
    logPosteriors = LPs,
    H = H
  ))
}

#' Expectation step
EM.E <- function(Xc, W, sigSq) {
  if (all(crossprod(W)==0)) {stop(paste("StPCA has failed due to numerical",
    "instability. Try dividing X by it's largest singular value."))}
  M = Matrix(crossprod(W) + sigSq*diag(ncol(W)))
  Vmean = unname(as.matrix(t(solve(M, t(Xc %*% W)))))
  Vvar1 = as.matrix(sigSq*solve(M))
  Vvar = lapply(1:nrow(Xc), function(i_) {
    unname(Vvar1 + tcrossprod(Vmean[i_,]))
  })
  return(list(Vmean=Vmean, Vvar=Vvar))
}

#' Maximization step for sigma^2
EM.M.sigSq <- function(Xc, W, Vmean, Vvar) {
  sigSqNew = (
    norm(Xc, 'F')^2 -
    2*sum(vapply(1:nrow(Xc), function(n_) Vmean[n_,] %*% crossprod(W, Xc[n_,]), numeric(1))) +
    sum(vapply(1:ncol(Xc), function(d_) W[d_,] %*% Reduce('+', Vvar) %*% W[d_,], numeric(1)))
  )/(nrow(Xc)*ncol(Xc))
  return(max(0, sigSqNew))
}

#' Maximization step for W
EM.M.W <- function(Xc, sigSq, Vmean, Vvar, K) {
  xvsum = Matrix(Reduce('+', lapply(1:nrow(Xc), function(i_) tcrossprod(Xc[i_,], Vmean[i_,]))))
  Vvarsum = as(Matrix(Reduce('+', Vvar)), "dppMatrix")

  A = K/sigSq
  B = solve(Vvarsum)
  C = A %*% t(solve(Vvarsum, t(xvsum)))

  W = sylSolve(A, B, C)

  return(W)
}
