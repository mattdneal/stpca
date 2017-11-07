context("stpca variable types")

library(fields)
data(ozone2)

# Missing data: Replace missing values by column means
X = ozone2$y
for (col in 1:ncol(X)) {
  ind = is.na(X[,col])
  X[ind,col] = mean(X[,col], na.rm=TRUE)
}
X = X/sd(X) # Scale for numerical reasons

locations = ozone2$lon.lat
locations = apply(locations, 2, function(col) (col-min(col))/(max(col)-min(col)))


d = ncol(X)
n = nrow(X)
k = 3

beta0 = cov.SE.beta0(X, locations, k)

stpca.init = stpca.init(X, k, locations, cov.SE, cov.SE.d, beta0, trace=0)
stpca.iter = stpca.iterate(stpca.init, maxit.inner=10, maxit.outer=1)

test_that("stpca.init returns a valid stpca object", {
  expect_is(stpca.init$X,     "matrix")
  expect_is(stpca.init$W,     "matrix")
  expect_is(stpca.init$sigSq, "numeric")
  expect_is(stpca.init$mu,    "numeric")
  expect_is(stpca.init$V,     "matrix")
  expect_is(stpca.init$ll,    "numeric")
  expect_is(stpca.init$lp,    "numeric")
  expect_is(stpca.init$lps,   "numeric")
  expect_is(stpca.init$bic,   "numeric")
  expect_true(length(stpca.init$beta0)==0 | is(stpca.init$beta0, "numeric"))
  expect_is(stpca.init$D,     "Matrix")
  expect_is(stpca.init$K,     "Matrix")
  expect_is(stpca.init$covar.fn,   "function")
  expect_is(stpca.init$covar.fn.d, "function")
  expect_is(stpca.init$locations,  "matrix")

  expect_equal(dim(stpca.init$X), c(n, d))
  expect_equal(dim(stpca.init$W), c(d, k))
  expect_length(stpca.init$sigSq, 1)
  expect_length(stpca.init$mu, d)
  expect_equal(dim(stpca.init$V), c(n, k))
  expect_length(stpca.init$ll, 1)
  expect_length(stpca.init$lp, 1)
  expect_length(stpca.init$lps, 1)
  expect_length(stpca.init$bic, 1)
  expect_equal(dim(stpca.init$D), c(d, d))
  expect_equal(dim(stpca.init$K), c(d, d))
  expect_equal(dim(locations)[1], d)
})

test_that("stpca.iterate returns a valid stpca object", {
  expect_is(stpca.iter$X,     "matrix")
  expect_is(stpca.iter$W,     "matrix")
  expect_is(stpca.iter$sigSq, "numeric")
  expect_is(stpca.iter$mu,    "numeric")
  expect_is(stpca.iter$V,     "matrix")
  expect_is(stpca.iter$ll,    "numeric")
  expect_is(stpca.iter$lp,    "numeric")
  expect_is(stpca.iter$lps,   "numeric")
  expect_is(stpca.iter$bic,   "numeric")
  expect_true(length(stpca.iter$beta0)==0 | is(stpca.iter$beta0, "numeric"))
  expect_is(stpca.iter$D,     "Matrix")
  expect_is(stpca.iter$K,     "Matrix")
  expect_is(stpca.iter$covar.fn,   "function")
  expect_is(stpca.iter$covar.fn.d, "function")
  expect_is(stpca.iter$locations,  "matrix")

  expect_equal(dim(stpca.iter$X), c(n, d))
  expect_equal(dim(stpca.iter$W), c(d, k))
  expect_length(stpca.iter$sigSq, 1)
  expect_length(stpca.iter$mu, d)
  expect_equal(dim(stpca.iter$V), c(n, k))
  expect_length(stpca.iter$ll, 1)
  expect_length(stpca.iter$lp, 1)
  expect_length(stpca.iter$lps, 21) # 1 + maxit.inner*(maxit.outer + 1)
  expect_length(stpca.iter$bic, 1)
  expect_equal(dim(stpca.iter$D), c(d, d))
  expect_equal(dim(stpca.iter$K), c(d, d))
  expect_equal(dim(locations)[1], d)
})
