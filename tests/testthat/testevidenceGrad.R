context("Evidence gradient")

test_that("The gradient of the log prior wrt beta matches numerical approximation", {
  grad.analytic <- unname(log_prior_d(stpca$WHat, stpca$beta,
                                      stpca$K, stpca$KD))

  grad.numeric = grad(function(beta_) {
    K = stpca$covFn(locs, beta=beta_)
    log_prior(K, stpca$WHat)
  }, x=stpca$beta)

  expect_equal(grad.analytic, grad.numeric)
})

test_that("The gradient of log(det(H)) matches numerical approximation", {
  grad.analytic = log_det_H_d(stpca$K, stpca$KD, stpca$H)

  grad.numeric = grad(function(beta_) {
    K = stpca$covFn(locs, beta=beta_)
    H = compute_H(X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K)
    sum(vapply(H, function(Hi) { # Log det H
      Matrix::determinant(Hi, logarithm=TRUE)$modulus
    }, numeric(1)))
  }, x=stpca$beta)

  expect_equal(grad.analytic, grad.numeric)
})

test_that("The gradient of the evidence matches numerical approxiation (backend)", {
  beta <- stpca$beta + rnorm(length(stpca$beta))*1e-2

  maxFn <- function(beta_) {
    K_ <- as(stpca$covFn(locs, beta=beta_), "dppMatrix")
    H_ <- compute_H(X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
    log_evidence(X, K_, stpca$WHat, stpca$muHat, stpca$sigSqHat, H_)
  }

  maxFnD <- function(beta_) {
    K_  <- as(stpca$covFn(locs, beta=beta_), "dppMatrix")
    KD_ <- stpca$covFnD(locs, beta=beta_)
    H_  <- compute_H(X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
    log_evidence_d(X, K_, stpca$WHat, stpca$muHat, stpca$sigSqHat,
                   beta_, KD_, H_)
  }

  DD = maxLik::compareDerivatives(maxFn, maxFnD, t0=beta, print=FALSE)
  gradDiff = unname(DD$compareGrad$rel.diff[1,])

  expect_equal(gradDiff, rep(0, length(beta)))
})

test_that("The gradient of the evidence matches numerical approxiation (frontend)", {
  beta <- stpca$beta + rnorm(length(stpca$beta))*1e-2
  stpcaCpy <- stpcaUp$copy()

  maxFn  <- function(beta_) stpcaCpy$set_beta(beta_)$logEvidence
  maxFnD <- function(beta_) stpcaCpy$set_beta(beta_)$logEvidenceD
  DD = maxLik::compareDerivatives(maxFn, maxFnD, t0=beta, print=FALSE)
  gradDiff = unname(DD$compareGrad$rel.diff[1,])

  expect_equal(gradDiff, rep(0, length(beta)))
})
