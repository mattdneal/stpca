context("Evidence calculation")

test_that("Log evidence returns -Inf if K cannot be inverted", {
  logEv = with(stpcaObj.it.t, {
    K_ = cov.SE(locations, beta=log(c(1, 10000))) # Cannot invert this!
    stpca.log_evidence(X, K_, W, mu, sigSq)
  })
  expect_equal(logEv, -Inf)
})

#test_that("Real HPs are within 95%CI of inferred HPs", {
#
#  set.seed(1)
#
#  require(mvtnorm)
#
#  # Generate some data
#  dim   = c(9, 9)
#  d     = prod(dim)
#  beta  = log(c(2))
#  cov.fn = functional::Curry(cov.independent, beta=beta)
#  sigSq = 1e-4
#
#  dataset = synthesize_data_kern(1000, d, dim, kern=cov.fn, noisesd=sqrt(sigSq))
#
#  locations = dataset$grid
#  X = as.matrix(dataset$X)
#  W = as.matrix(dataset$W)
#
#  # Function mapping HP values to log evidence
#  evF = function(beta_) {
#    K_ = cov.independent(locations, beta=beta_)
#    stpca.log_evidence(X, K_, W, 0, sigSq)
#  }
#
#  # Find maxima of log evidence wrt HP
#  OO = optimise(f=evF, interval=log(c(1.5, 2.5)), maximum=TRUE)
#  betaHat = OO$maximum
#
#  sd = sqrt(-1/as.numeric(numDeriv::hessian(evF, x=betaHat)))
#  ci95 = betaHat + 2.576*c(sd, -sd)
#
#  x = log(seq(1.5, 2.5, length=10))
#  y = vapply(x, evF, numeric(1))
#  plot(exp(x), y)
#})
