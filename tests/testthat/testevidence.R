context("Evidence calculation")

test_that("Log evidence returns -Inf if K cannot be inverted", {
  logEv = with(stpcaObj.it.t, {
    K_ = cov.SE(locations, beta=log(c(1, 10000))) # Cannot invert this!
    stpca.log_evidence(X, K_, W, mu, sigSq)
  })
  expect_equal(logEv, -Inf)
})
