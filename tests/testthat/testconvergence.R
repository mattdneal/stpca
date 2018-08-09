context("Global convergence")

test_that("At convergence, the value of theta maximises the posterior", {
  thetaHatList = list(
    W     = stpcaUp$WHat,
    mu    = stpcaUp$muHat,
    sigSq = stpcaUp$sigSqHat
  )

  theta_hat_condition = function(theta) {
    thetaList <- relist(theta, thetaHatList)
    W      <- thetaList$W
    mu     <- thetaList$mu
    sigSq  <- thetaList$sigSq
    (log_likelihood(stpcaUp$X, W, mu, sigSq) +
     log_prior(stpcaUp$K, W, sigSq))
  }

  thetaGrad = grad(theta_hat_condition, unlist(thetaHatList))

  expect_equal(thetaGrad, rep(0, length(thetaGrad)), tol=1e-4)
})

test_that("At convergence, the value of beta maximises the approximate evidence", {
  beta_hat_condition = function(beta_) {
    K_ <- stpcaUp$covFn(stpcaUp$locs, beta=beta_)
    H_ <- compute_H(stpcaUp$X, stpcaUp$WHat, stpcaUp$muHat,
                    stpcaUp$sigSqHat, K_)
    log_evidence(stpcaUp$X, K_, stpcaUp$WHat, stpcaUp$muHat,
                 stpcaUp$sigSqHat, H_)
  }

  betaGrad = grad(beta_hat_condition, stpcaUp$beta)

  expect_equal(betaGrad, rep(0, length(betaGrad)), tol=5e-3)
})

