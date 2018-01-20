context("High level tests of StpcaModel")

test_that("update() changes fields as expected", {
  stpcaCpy = stpca$copy()$update(1, 100, constraints=constr)

  expect_false(isTRUE(all.equal(stpca$beta, stpcaCpy$beta)))
  expect_false(isTRUE(all.equal(stpca$K,    stpcaCpy$K)))
  expect_false(isTRUE(all.equal(stpca$KD,   stpcaCpy$KD)))
  expect_false(isTRUE(all.equal(stpca$WHat, stpcaCpy$WHat)))
  expect_false(isTRUE(all.equal(stpca$Vmean,stpcaCpy$Vmean)))
  expect_false(isTRUE(all.equal(stpca$Vvar ,stpcaCpy$Vvar)))
  expect_false(isTRUE(all.equal(stpca$convergence,   stpcaCpy$convergence)))
  expect_false(isTRUE(all.equal(stpca$logEvidence,   stpcaCpy$logEvidence)))
  expect_false(isTRUE(all.equal(stpca$logPosteriors, stpcaCpy$logPosteriors)))
  expect_false(isTRUE(all.equal(stpca$H,     stpcaCpy$H)))
  expect_false(isTRUE(all.equal(stpca$maxim, stpcaCpy$maxim)))
})

test_that("set_beta() gives informative errors for bad hyperparameter values", {
  stpca2 = stpca$copy()

  expect_error(stpca2$set_beta(stpca$beta + 10000), "non-finite")
  expect_error(stpca2$set_beta(stpca$beta - 10000), "non-finite")

  # Contents of stpca2 are unchanged by bad beta setting
  expect_equal(stpca2$beta, stpca$beta)
})
