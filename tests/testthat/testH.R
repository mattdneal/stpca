context("Numerical vs analytic blocks of H")

test_that("Analytic H_{\\sigma^2} is equal to numeric H_{\\sigma^2}", {
  HsigSq.analytic = unname(as.matrix(stpca$H$sigSq))
  HsigSq.numeric = numDeriv::hessian(function(sigSq_) {
    -log_likelihood(X, stpca$WHat, stpca$muHat, sigSq_)
  }, x=stpca$sigSqHat)
  expect_equal(HsigSq.analytic, HsigSq.numeric)
})

test_that("Analytic H_{\\mu} is equal to numeric H_{\\mu}", {
  Hmu.analytic = unname(as.matrix(stpca$H$mu))
  Hmu.numeric = numDeriv::hessian(function(mu_) {
    -(log_likelihood(stpca$X, stpca$WHat, mu_, stpca$sigSqHat) +
      log_prior(stpca$K, stpca$WHat, stpca$sigSqHat))
  }, x=stpca$muHat)
  expect_equal(Hmu.analytic, Hmu.numeric, tolerance=1e-6)
})

test_that("Analytic H_{w_i} are all equal to numeric H_{w_i}", {
  for (i in 1:k) {
    Hwi.analytic = unname(as.matrix(stpca$H[[paste("w", i, sep='')]]))
    Hwi.numeric = numDeriv::hessian(function(wi) {
      W_ = stpca$WHat
      W_[,i] = wi
      -(log_likelihood(stpca$X, W_, stpca$muHat, stpca$sigSqHat) +
        log_prior(stpca$K, W_, stpca$sigSqHat))
    }, x=stpca$WHat[,i])

    expect_equal(as.matrix(Hwi.analytic), Hwi.numeric, tolerance=1e-5,
                 check.attributes=FALSE)
  }
})

test_that("Building H triggers an informative error if K cannot be inverted", {
  expect_error({
    K_ = cov.SE(locs, beta=log(c(1, 100000)))
    compute_H(X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
  }, "Could not invert K.")

  # More informative if beta attr is present
  expect_error({
    beta <- log(c(1, 100000))
    K_ <- cov.SE(locs, beta=beta)
    attr(K_, "beta") <- beta
    compute_H(X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
  }, "Could not invert K. beta=0,11.513")
})
