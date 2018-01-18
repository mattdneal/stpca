context("Log evidence")

test_that("Evidence calculation is consistent between frontend and backend", {
  evBack <- function(beta) {
    K_   <- stpca$covFn(stpca$locs, beta=beta)
    vals <- theta_EM(stpca$X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_, maxit=100)
    vals$logEvidence
  }

  evFrnt <- function(beta) {
    unname(stpca$copy()$set_beta(beta)$update_theta(100)$logEvidence)
  }

  beta = stpca$beta + rnorm(length(stpca$beta))*1e-1
  expect_equal(evFrnt(beta), evBack(beta))
})
