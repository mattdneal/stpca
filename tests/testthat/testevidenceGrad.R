for (covfun in 1:2) {
  if (covfun==1) {
    covFn  = cov.independent
    covFnD = cov.independent.d
    beta0  = log(1.1)
    context("Evidence gradient with independent covar fun")
  } else {
    covFn  = cov.noisy.SE
    covFnD = cov.noisy.SE.d
    beta0  = c(cov.SE.beta0(X, locs, k), log(1e-2))
    context("Evidence gradient with noisy SE covar fun")
  }

  stpca2 = StpcaModel$new(X, k, beta0, locs, covFn, covFnD)
  stpca2$compute_gradient()

  test_that("The gradient of the log prior wrt beta matches numerical approximation", {
    grad.analytic = unname(log_prior_d(stpca2$WHat, stpca2$beta,
                                       stpca2$K, stpca2$KD))

    grad.numeric = grad(function(beta_) {
      K_ = stpca2$covFn(locs, beta=beta_)
      return(log_prior(K_, stpca2$WHat))
    }, x=stpca2$beta)

    expect_equal(grad.analytic, grad.numeric, tol=1e-6)
  })

  test_that("The gradient of log(det(H)) matches numerical approximation", {
    grad.numeric = grad(function(beta_) {
      K_ = stpca2$covFn(locs, beta=beta_)
      HW_ = compute_H_W(X, stpca2$WHat, stpca2$muHat, stpca2$sigSqHat, K_)
      sum(vapply(HW_, function(Hwi) {
        Matrix::determinant(Hwi, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpca2$beta)

    HW = compute_H_W(stpca2$X, stpca2$WHat, stpca2$muHat, stpca2$sigSqHat, stpca2$K)
    grad.analytic = log_det_H_d(stpca2$K, stpca2$KD, HW)

    expect_equal(grad.analytic, grad.numeric, scale=mean(grad.numeric))
  })

  test_that("Deriv of log(det(H)) wrt beta should be the same as the deriv of log(det(H_W))", {
    grad.numeric = grad(function(beta_) {
      K_ = stpca2$covFn(locs, beta=beta_)
      HW_ = compute_H_W(stpca2$X, stpca2$WHat, stpca2$muHat, stpca2$sigSqHat, K_)
      sum(vapply(HW_, function(Hwi) {
        Matrix::determinant(Hwi, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpca2$beta)

    grad.numeric2 = grad(function(beta_) {
      K_ = stpca2$covFn(locs, beta=beta_)
      H_ = compute_H_W(stpca2$X, stpca2$WHat, stpca2$muHat, stpca2$sigSqHat, K_)
      sum(vapply(H_, function(Hblock) {
        Matrix::determinant(Hblock, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpca2$beta)

    expect_equivalent(grad.numeric, grad.numeric2)
  })

  test_that("The gradient of the evidence matches a numerical approximation", {
    grad.numeric = grad(function(beta_) {
      K_ = stpca2$covFn(locs, beta=beta_)
      H_ = compute_H_W(stpca2$X, stpca2$WHat, stpca2$muHat, stpca2$sigSqHat, K_)
      log_evidence(stpca2$X, K_, stpca2$WHat, stpca2$muHat, stpca2$sigSqHat, H_)
    }, x=stpca2$beta)

    KD = stpca2$covFnD(locs, beta=stpca2$beta)
    grad.analytic = unname(log_evidence_d(
      stpca2$X,
      stpca2$K,
      stpca2$WHat,
      stpca2$muHat,
      stpca2$sigSqHat,
      stpca2$beta,
      KD,
      stpca2$H
    ))

    expect_equal(grad.analytic, grad.numeric, tol=1e-6)
  })
}
