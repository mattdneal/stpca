context("Crossvalidation")

expect_that("crossvalidating gives sensible output", {
  cvll <- mean(stpcaUp$crossvalidate()$ll)
  expect_true(is.finite(cvll))
})

expect_that("crossvalidating picks most likely beta", {
  stpcaCpy     <- stpcaUp$copy()$set_beta(beta+10)
  cvllLikely   <- mean(stpcaUp$crossvalidate()$ll)
  cvllUnlikely <- mean(stpcaCpy$crossvalidate()$ll)
  expect_gt(cvllLikely, cvllUnlikely)
})

expect_that("crossvalidating does not change stpca object", {
  stpcaCpy <- stpcaUp$copy()
  stpcaCpy$crossvalidate(3, 3)
  expect_equal(stpcaUp, stpcaCpy)
})
