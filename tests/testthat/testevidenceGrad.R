for (covfun in 1:2) {
  if (covfun==1) {
    covar.fun   = cov.independent
    covar.fun.d = cov.independent.d
    beta0       = log(1.1)
    context("Evidence gradient with independent covar fun")
  } else {
    covar.fun   = cov.noisy.SE
    covar.fun.d = cov.noisy.SE.d
    beta0       = c(cov.SE.beta0(X, locations, k), log(1e-2))
    context("Evidence gradient with SE covar fun")
  }

  stpcaObj = stpca.init(X, k, locations, covar.fun, covar.fun.d, beta0, trace=0)

  test_that("The gradient of the log prior wrt beta matches numerical approximation", {
    dK = stpcaObj$covar.fn.d(locations, beta=stpcaObj$beta)
    grad.analytic = stpca.log_prior_d(stpcaObj$W, stpcaObj$beta,
                                      stpcaObj$K, dK)

    grad.numeric = grad(function(beta_) {
      K = stpcaObj$covar.fn(locations, beta=beta_)
      return(stpca.log_prior(K, stpcaObj$W))
    }, x=stpcaObj$beta)

    expect_equivalent(grad.analytic, grad.numeric)
  })

  test_that("The gradient of log(det(H)) matches numerical approximation", {
    grad.numeric = grad(function(beta_) {
      K_ = stpcaObj$covar.fn(locations, beta=beta_)
      HW_ = with(stpcaObj, stpca.H.W(X, W, mu, sigSq, K_))
      sum(vapply(HW_, function(Hwi) {
        Matrix::determinant(Hwi, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpcaObj$beta)

    dK = stpcaObj$covar.fn.d(locations, beta=stpcaObj$beta)
    K  = stpcaObj$covar.fn(locations, beta=stpcaObj$beta)
    HW = with(stpcaObj, stpca.H.W(X, W, mu, sigSq, K))
    grad.analytic = stpca.log_det_H_d(stpcaObj$K, dK, HW)

    expect_equal(grad.analytic, grad.numeric, scale=mean(grad.numeric))
  })

  test_that("Deriv of log(det(H)) wrt beta should be the same as the deriv of log(det(H_W))", {
    grad.numeric = grad(function(beta_) {
      K_ = stpcaObj$covar.fn(locations, beta=beta_)
      HW_ = with(stpcaObj, stpca.H.W(X, W, mu, sigSq, K_))
      sum(vapply(HW_, function(Hwi) {
        Matrix::determinant(Hwi, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpcaObj$beta)

    grad.numeric2 = grad(function(beta_) {
      K_ = stpcaObj$covar.fn(locations, beta=beta_)
      H_ = with(stpcaObj, stpca.H.W(X, W, mu, sigSq, K_))
      sum(vapply(H_, function(Hblock) {
        Matrix::determinant(Hblock, logarithm=TRUE)$modulus
      }, numeric(1)))
    }, x=stpcaObj$beta)

    expect_equivalent(grad.numeric, grad.numeric2)
  })

  test_that("The gradient of the evidence matches a numerical approximation", {
    grad.numeric = grad(function(beta_) { with(stpcaObj, {
      dK_ = covar.fn.d(locations, beta=beta_)
      K_  = covar.fn(locations, beta=beta_)
      return(stpca.log_evidence(Xc, K_, W, rep(0,ncol(Xc)), sigSq))
    })}, x=stpcaObj$beta)

    dK = stpcaObj$covar.fn.d(locations, beta=stpcaObj$beta)
    grad.analytic = unname(stpca.log_evidence_d(
      stpcaObj$Xc,
      stpcaObj$K,
      stpcaObj$W,
      rep(0,d),
      stpcaObj$sigSq,
      stpcaObj$beta,
      dK
    ))

    expect_equivalent(grad.analytic, grad.numeric)
  })
}
