context("sylvester equation solving")

test_that("sylSolve matches analytic solutions", {
  set.seed(1)
  # B = identity
  A = crossprod(matrix(rnorm(d*d), ncol=d))
  C = matrix(rnorm(d*k), nrow=d, ncol=k)
  expect_equal(sylSolve(A, diag(k), C), solve(A + diag(d), C))

  # A = identity
  B = crossprod(matrix(rnorm(k*k), ncol=k))
  C = matrix(rnorm(d*k), nrow=d, ncol=k)
  expect_equal(sylSolve(diag(d), B, C), t(solve(B + diag(k), t(C))))
})

test_that("sylSolve solution for W satisfies AW + WB = C", {
  set.seed(1)
  for (i in 1:10) {
    A = Matrix(crossprod(matrix(rnorm(d*d), ncol=d)))
    B = Matrix(crossprod(matrix(rnorm(k*k), ncol=k)))
    C = Matrix(matrix(rnorm(d*k), nrow=d, ncol=k))
    W = sylSolve(A, B, C)
    expect_equal(A%*%W + W%*%B, C)
  }
})
