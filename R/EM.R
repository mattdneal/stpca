#' Update theta to be the maximum-a-posteriori value using Expectation
#' Maximisation.
theta_EM <- function(X, W, mu, sigSq, K, maxit=500, bftol=1e-5, sparse=FALSE, b=NULL) {
  stopifnot(bftol > 0)
  stopifnot(maxit > 0)

  k  = ncol(W)
  Xc = sweep(X, 2, mu)

  unNormedLPs = numeric(maxit) # un-normalised log posteriors
  for (iter in seq_len(maxit)) { # TODO: Iterate until zero gradient
    # Expectation step
    E = EM.E(Xc, W, sigSq, sparse=sparse)

    # Maximization step for W
    if (sparse) {
      W = EM.M.W.sparse(Xc, sigSq, E$Vmean, E$Vvar, E$colVmag, E$RtV, K, b)
    } else {
      W = EM.M.W(Xc, sigSq, E$Vmean, E$Vvar, K)
    }

    # Align W with PCs & fix sign
    if (!sparse) {
      W <- W %*% svd(W)$v
      W <- W %*% diag(sign(W[1,]), nrow=k, ncol=k)
    } else {
      magnitudes <- apply(W, 2, function(w) norm(w, "2"))
      W <- W[,order(magnitudes, decreasing=TRUE)]
    }

    # Second expectation step using updated W
    E = EM.E(Xc, W, sigSq, sparse=sparse)

    # Maximization step for sigma^2
    sigSq = EM.M.sigSq(Xc, W, E$Vmean, E$Vvar)

    # ln p(X | \theta) + ln p(\theta | \beta)
    unNormedLPs[iter] = log_likelihood(Xc, W, 0, sigSq) + log_prior(K, W, sigSq)

    # If the log Bayes' Factor ratio is less than bftol then stop early:
    # log(p(theta^n | X) / p(theta^{n-1} | X)) < bftol
    if (iter>1 && (unNormedLPs[iter]-unNormedLPs[iter-1])<bftol) break
  }

  # Truncate unassigned LPs
  if (iter < maxit) unNormedLPs = unNormedLPs[seq_len(iter)]

  EMout <- list(
    WHat = W,
    muHat = mu,
    sigSqHat = sigSq,
    Vmean = E$Vmean,
    Vvar = E$Vvar,
    logPosteriors = unNormedLPs
  )

  if (!sparse) {
    # Assuming convergence, \sigma = \hat\sigma, so the evidence can be approximated
    H = compute_H(X, W, mu, sigSq, K)
    logEvidence = log_evidence(Xc, K, W, 0, sigSq, H)

    # We now know the evidence to normalize the above log posterior values
    EMout$logEvidence = logEvidence
    EMout$logPosteriors = unNormedLPs - logEvidence
    EMout$H = H
  }
  return(EMout)
}

#' Expectation step
EM.E <- function(Xc, W, sigSq, sparse=FALSE) {
  if (all(crossprod(W)==0)) {stop(paste("StPCA has failed due to numerical",
    "instability. Try dividing X by it's largest singular value."))}
  M = Matrix(crossprod(W) + sigSq*diag(ncol(W)))
  Vmean = unname(as.matrix(t(solve(M, t(Xc %*% W)))))
  Vvar1 = as.matrix(sigSq*solve(M))
  Vvar = lapply(1:nrow(Xc), function(i_) {
    unname(Vvar1 + tcrossprod(Vmean[i_,]))
  })
  if (sparse) {
    Minv <- solve(M)
    colVmag <- vapply(1:ncol(W), function(l) {
      sum(diag(sigSq*Minv[l,l]*Diagonal(nrow(Xc)) + tcrossprod(Vmean[,l])))
    }, numeric(1))

    RtV <- Matrix(NA, nrow=nrow(W), ncol=ncol(W))
    for (l in 1:ncol(W)) {
      for (m in 1:nrow(W)) {
        RtV[m,l] <- crossprod(Xc[,m] -
          tcrossprod(Vmean[,-l,drop=FALSE], W[m,-l,drop=FALSE]),
        Vmean[,l])
      }
    }
    return(list(Vmean=Vmean, Vvar=Vvar, colVmag=colVmag, RtV=RtV))
  } else {
    return(list(Vmean=Vmean, Vvar=Vvar))
  }
}

#' Maximization step for sigma^2
EM.M.sigSq <- function(Xc, W, Vmean, Vvar) {
  sigSqNew = (
    norm(Xc, 'F')^2 -
    2*sum(vapply(1:nrow(Xc), function(n_) Vmean[n_,] %*% crossprod(W, Xc[n_,]), numeric(1))) +
    sum(vapply(1:ncol(Xc), function(d_) W[d_,] %*% Reduce('+', Vvar) %*% W[d_,], numeric(1)))
  )/(nrow(Xc)*ncol(Xc) + 2)
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

# TODO: Write unit tests
EM.M.W.sparse <- function(Xc, sigSq, Vmean, Vvar, colVmag, RtV, K, b, reldiff=1e-4) {
  # Initialize from old code.
  W <- EM.M.W(Xc, sigSq, Vmean, Vvar, K)

  # TODO: Tets out thresholding the initialised W: might be better?

  # Useful variables
  n <- nrow(Xc)
  d <- ncol(Xc)
  k <- ncol(Vmean)
  Kinv <- solve(K)

  converged <- FALSE
  lp <- log_likelihood(Xc, W, 0, sigSq) + log_sparse_prior(K, W, sigSq, b)
  while (!converged) {
    for (m in seq_len(d)) {
      for (l in seq_len(k)) {
        W[m,l] <- w_coord_desc(Xc, l, m, W, sigSq, colVmag, RtV, K, b, Kinv)
      }
    }
    lpNew <- log_likelihood(Xc, W, 0, sigSq) + log_sparse_prior(K, W, sigSq, b)
    converged <- (lpNew - lp < reldiff)
    lp <- lpNew
  }

  return(W)
}
