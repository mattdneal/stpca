#test_that("tune_beta() method optimises beta", {
#  stpca2 = stpca$copy()
#  stpca2$tune_beta()
#})

test_that("set_beta() plays nice with extreme hyperparameter values", {
  stpca2 = stpca$copy()

  # No error with small values
  expect_silent(stpca2$set_beta(stpca$beta + 10000))
  expect_equal(stpca2$logEvidence, -Inf)

  # No error with big values
  expect_silent(stpca2$set_beta(stpca$beta - 10000))
  expect_equal(stpca2$logEvidence, -Inf)

  # Returning to acceptable values works fine
  expect_silent(stpca2$set_beta(stpca$beta))
  expect_that(is.finite(stpca2$logEvidence), is_true())
})
