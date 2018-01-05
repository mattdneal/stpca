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

  stpca = StpcaModel$new(X, k, beta0, locs, covFn, covFnD)
  stpca$compute_gradient()

  test_that("The gradient of the log prior wrt beta matches numerical approximation", {
    grad.analytic = log_prior_d(stpca$WHat, stpca$beta,
                                stpca$K, stpca$KD)

    grad.numeric = grad(function(beta_) {
      K_ = stpca$covFn(locs, beta=beta_)
      return(log_prior(K_, stpca$WHat))
    }, x=stpca$beta)

    expect_equivalent(grad.analytic, grad.numeric)
  })

  test_that("The gradient of log(det(H)) matches numerical approximation", {
    grad.numeric = grad(function(beta_) {
      K_ = stpca$covFn(locs, beta=beta_)
      HW_ = compute_H_W(X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
      sum(vapply(HW_, function(Hwi) {
        Matrix::determinant(Hwi, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpca$beta)

    HW = compute_H_W(stpca$X, stpca$WHat, stpca$muHat, stpca$sigSqHat, stpca$K)
    grad.analytic = log_det_H_d(stpca$K, stpca$KD, HW)

    expect_equal(grad.analytic, grad.numeric, scale=mean(grad.numeric))
  })

  test_that("Deriv of log(det(H)) wrt beta should be the same as the deriv of log(det(H_W))", {
    grad.numeric = grad(function(beta_) {
      K_ = stpca$covFn(locs, beta=beta_)
      HW_ = compute_H_W(stpca$X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
      sum(vapply(HW_, function(Hwi) {
        Matrix::determinant(Hwi, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpca$beta)

    grad.numeric2 = grad(function(beta_) {
      K_ = stpca$covFn(locs, beta=beta_)
      H_ = compute_H_W(stpca$X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
      sum(vapply(H_, function(Hblock) {
        Matrix::determinant(Hblock, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpca$beta)

    expect_equivalent(grad.numeric, grad.numeric2)
  })

  test_that("The gradient of the evidence matches a numerical approximation", {
    grad.numeric = grad(function(beta_) {
      KD_ = stpca$covFnD(locs, beta=beta_)
      K_  = covFn(locs, beta=beta_)
      return(log_evidence(stpca$X, K_, stpca$WHat, stpca$muHat, stpca$sigSqHat))
    }, x=stpca$beta)

    KD = stpca$covFnD(locs, beta=stpca$beta)
    grad.analytic = unname(log_evidence_d(
      stpca$X,
      stpca$K,
      stpca$WHat,
      stpca$muHat,
      stpca$sigSqHat,
      stpca$beta,
      KD
    ))

    expect_equivalent(grad.analytic, grad.numeric)
  })
}
