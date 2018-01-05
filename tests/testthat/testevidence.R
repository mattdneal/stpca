context("Evidence calculation")

test_that("Log evidence returns -Inf if K cannot be inverted", {
  K_ = cov.SE(locs, beta=log(c(1, 10000))) # Cannot invert this!
  logEv = log_evidence(X, K_, stpca$WHat, stpca$muHat, stpca$sigSqHat)
  expect_equal(logEv, -Inf)
})
