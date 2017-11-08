context("expectation")

n = 10
d = 20
k = 3
Xc = scale(matrix(rnorm(n*d), nrow=n, ncol=d), scale=FALSE)
sigSq = 1
W     = diag(rep(1, k), nrow=d, ncol=k)

test_that("The expectation function behaves as expected", {
  expectations = EM.E(Xc, W, sigSq)

  expect_equal(expectations$V, 0.5*Xc[,1:k])

  for (i in 1:n) {
    expect_is(expectations$Vvar[[i]], 'matrix')
    expect_equal(dim(expectations$Vvar[[i]]), c(k, k))
    expect_equal(expectations$Vvar[[i]], 0.5*diag(k) + 0.25*tcrossprod(Xc[i,1:k]))
  }
})

## Generate the data
library(functional)

n     = 15
k     = 4
dim   = c(11, 11)
d     = prod(dim)
beta  = log(c(2, 0.2))
k_se  = Curry(cov.SE, beta=beta)
sigSq = 1.8

dataset = synthesize_data_kern(n, k, dim, kern=k_se, noisesd=sqrt(sigSq))
locations = dataset$grid
X = dataset$X

maxit.inner=3
maxit.outer=1

beta0 = cov.SE.beta0(Xc, locations, k)

stpcaObj = stpca(X, k, locations, cov.SE, cov.SE.d, beta0, trace=0, maxit.inner=10, maxit.outer=0)
stpcaObj = stpca.iterate.beta(stpcaObj) # Update beta, now the thetas need updating

context("maximisation in sigma^2")

test_that("Iterating Expectation and sigma^w maximisation finds zero-gradient point", {
  Xc = stpcaObj$Xc
  K = stpcaObj$K
  W = stpcaObj$W
  mu = rep(0, d)
  sigSq = stpcaObj$sigSq
  for (i in 1:200) {
    E     = EM.E(Xc, stpcaObj$W, sigSq)
    sigSq = EM.M.sigSq(Xc, stpcaObj$W, E$V, E$Vvar)
  }

  lp = stpca.log_posterior(Xc, K, W, mu, sigSq)

  delta = 1e-8
  lpUp = stpca.log_posterior(Xc, K, W, mu, sigSq+delta)
  lpDn = stpca.log_posterior(Xc, K, W, mu, sigSq-delta)
  gradient = (lpUp-lpDn)/(2*delta)
  expect_lt(gradient, 1e-3)
})

test_that("sigma^2 maximisation stage increases complete log posterior", {
  clp1 = with(stpcaObj, stpca.complete_log_posterior(X, V, W, mu, sigSq, K))
  sigSqNew = with(stpcaObj, EM.M.sigSq(Xc, W, V, Vvar))
  clp2 = with(stpcaObj, stpca.complete_log_posterior(X, V, W, mu, sigSqNew, K))
  expect_gt(clp2, clp1)
})

context("maximisation in W")

#Wrong! It increases the COMPLETE Log posterior!
#test_that("W maximisation stage increases log posterior", {
#  Wmax  = with(stpcaObj, EM.M.W(Xc, sigSq, V, Vvar, K))
#  lpOld = stpcaObj$lp
#  lpMax = with(stpcaObj, stpca.log_posterior(Xc, K, Wmax, rep(0,ncol(Xc)), sigSq))
#  expect_gt(lpMax, lpOld)
#})

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
  A = Matrix(crossprod(matrix(rnorm(d*d), ncol=d)))
  B = Matrix(crossprod(matrix(rnorm(k*k), ncol=k)))
  C = Matrix(matrix(rnorm(d*k), nrow=d, ncol=k))
  W = sylSolve(A, B, C)
  expect_equal(A%*%W + W%*%B, C)
})
