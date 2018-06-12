context("Constraints")

test_that("Can start beta at boundary of feasible region", {
  ss <- stpca$
          copy()$
          set_beta(-constr$ineqB[1:3])$
          update_beta(constraints=constr)
  expect_gt(ss$logEvidence, stpca$logEvidence)
})
