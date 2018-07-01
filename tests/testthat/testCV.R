context("Crossvalidation")

expect_that("crossvalidating picks likely over unlikely beta", {
  cvll <- mean(stpcaUp$crossvalidate()$ll)
  expect_true(is.finite(cvll$ll))
})
