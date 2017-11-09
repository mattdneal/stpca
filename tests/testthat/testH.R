context("Numerical vs analytic blocks of H")

test_that("Analytic H_{\\sigma^2} is equal to numeric H_{\\sigma^2}", {
  HsigSq.analytic = H$sigSq
  HsigSq.numeric = Matrix(numDeriv::hessian(function(sigSq_) {
    -stpca.log_posterior(X, stpcaObj$K, stpcaObj$W, stpcaObj$mu, sigSq_)
  }, x=stpcaObj$sigSq))
  expect_equal(HsigSq.analytic, HsigSq.numeric)
})

test_that("Analytic H_{\\mu} is equal to numeric H_{\\mu}", {
  Hmu.analytic = H$mu
  Hmu.numeric = Matrix(numDeriv::hessian(function(mu_) {
    -stpca.log_posterior(X, stpcaObj$K, stpcaObj$W, mu_, stpcaObj$sigSq)
  }, x=stpcaObj$mu))
  expect_equal(Hmu.analytic, Hmu.numeric, tolerance=1e-6)
})

test_that("Analytic H_{w_i} are all equal to numeric H_{w_i}", {
  for (i in 1:k) {
    Hwi.analytic = H[[paste("w", i, sep='')]]
    Hwi.numeric = Matrix(numDeriv::hessian(function(wi) {
      W = stpcaObj$W
      W[,i] = wi
      -stpca.log_posterior(X, stpcaObj$K, W, stpcaObj$mu, stpcaObj$sigSq)
    }, x=stpcaObj$W[,i]))

    expect_equal(Hwi.analytic, Hwi.numeric)
  }
})
