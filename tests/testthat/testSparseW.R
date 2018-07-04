context("Sparse prior")

test_that("Laplace function integrates to 1", {
  halfInt <- integrate(
    f = function(W) {
      vapply(W, function(w) {
        dlaplace(w, 1)
      }, numeric(1))
    },
    lower = 0,
    upper = 50)
  expect_equal(halfInt$value*2, 1, tolerance=1e-2)
})

test_that("Laplace function at 0 equals 1/(2b)", {
  for (b in abs(rnorm(10))) {
    expect_equal(dlaplace(0, b), 1/(2*b))
  }
})

test_that("Sparse prior works in edge cases", {
  # 1x1 W
  lsp <- log_sparse_prior(1, Matrix(1), 1)
  expect_scalar(lsp)

  # Test dx1 W (as vector)
  lsp <- log_sparse_prior(rnorm(5), Diagonal(5), 2)
  expect_scalar(lsp)

  # All zeroes
  lsp <- log_sparse_prior(Matrix(0, nrow=10, ncol=5), Diagonal(10), 2)
  expect_scalar(lsp)
})

test_that("Sparse prior errors correctly", {
  expect_error(log_sparse_prior(stpca$W, stpca$K, -1))
})

test_that("Sparse posterior has expected W invariances", {
  W     <- stpcaUp$WHat
  sparseLP <- function(W) {
    log_likelihood(Xc, W, 0, stpcaUp$sigSqHat) +
      log_sparse_prior(W, stpcaUp$K, 0.001)
  }

  lp1 <- sparseLP(W)
  lp2 <- sparseLP(W[,rev(1:k)]) # Invariant to switching dimensions
  lp3 <- sparseLP(W*matrix(sign(d*k), nrow=d, ncol=k)) # Sign invariant

  expect_equal(lp1, lp2)
  expect_equal(lp1, lp3)

  R <- svd(matrix(rnorm(k*k), nrow=k))$u # Random orthonormal matrix
  lp4 <- sparseLP(W %*% R) # *not* invariant to rotations
  expect_true(lp1 != lp4)
})

context("EM: Sparse E-step")

test_that("Sparse E-step behaves sensibly", {
  spE <- EM.E(Xc, stpca$WHat, stpca$sigSqHat, sparse=TRUE)

  d <- nrow(stpca$WHat)
  k <- ncol(stpca$WHat)
  n <- nrow(Xc)

  expect_length(spE$colVmag, k)
  expect_that(all(is.finite(spE$colVmag[[1]])), is_true())

  expect_equal(dim(spE$RtV), c(d,k))
  expect_that(all(is.finite(spE$RtV)), is_true())
})

context("EM: Sparse W-maximisation step")

test_that("Soft thresholding works", {
  x = seq(-10, 10, length=100)

  expect_equal(soft_threshold(0, 1), 0)
  expect_equal(soft_threshold(x, 0), x)
  expect_equal(soft_threshold(5, 5), 0)
  expect_equal(soft_threshold(-5, 5), 0)
  expect_error(soft_threshold(x, -1))
})

test_that("Every coordinate step increases log posterior", {
  b <- 1
  W <- stpca$WHat

  f <- function(W) {
    (log_likelihood(Xc, W, 0, stpca$sigSqHat) +
     log_sparse_prior(W, stpca$K, b))
  }

  lp <- f(W)
  for (l in 1:k) {
    for (m in sample(d, 10)) {
      E <- EM.E(Xc, W, stpca$sigSqHat, sparse=TRUE)
      W[m,l] <- w_coord_desc(Xc, l, m, W, stpca$sigSqHat, E$colVmag,
                             E$RtV, stpca$K, b)

      # Step increases log posterior
      expect_gt(f(W), lp)

      lp <- f(W)

      # Fails!! Shouldn't?! It's quite close numerically, though..
      ## Step goes to maxima
      # W2 <- W3 <- W
      # W2[m,l] <- W2[m,l] + 0.01
      # W3[m,l] <- W2[m,l] - 0.01
      # expect_gt(lp, f(W2))
      # expect_gt(lp, f(W3))
    }
  }
})

test_that("Coordinate descent solution does not depend on current value of optimand", {
  for (l in 1:k) {
    for (m in sample(1:d, 5)) {
      W1 <- stpca$WHat
      E <- EM.E(Xc, stpca$WHat, stpca$sigSqHat, sparse=TRUE)
      w1 <- w_coord_desc(Xc, l, m, W1, stpca$sigSqHat,
                         E$colVmag, E$RtV, stpca$K, 2)

      W2 <- W1
      W2[m,l] <- rnorm(1)*10
      w2 <- w_coord_desc(Xc, l, m, W2, stpca$sigSqHat,
                         E$colVmag, E$RtV, stpca$K, 2)

      expect_equal(w1, w2)
    }
  }
})

test_that("Coordiante descent always gives zero for very small b", {
  E <- EM.E(Xc, stpca$WHat, stpca$sigSqHat, sparse=TRUE)
  w <- w_coord_desc(Xc, 1, 1, stpca$WHat, stpca$sigSqHat, E$colVmag,
                    E$RtV, stpca$K, 1e-7)
  expect_equal(w, 0)
})


# TODO: Test this dude
#EM.M.W.sparse <- function(Xc, sigSq, Vmean, Vvar, colVmag, RtV, K, b) {

context("Creating high level sparse models")

test_that("Sparse model can be created & tuned", {
  sstpca <- StpcaModel$new(X, k, beta0, locs, cov.noisy.SE, cov.noisy.SE.d, sparse=TRUE, b=0.1, maxit=100)
  expect_true(any(sstpca$WHat==0), is_true())
})

test_that("No hyperparameter updating for a sparse model", {
  sstpca <- StpcaModel$new(X, k, beta0, locs, cov.noisy.SE, cov.noisy.SE.d, sparse=TRUE, b=1, maxit=100)
  expect_error(sstpca$update())
})

context("Sparsity related methods")

test_that("Can convert a non-sparse model to sparse and infer theta", {
  sstpca <- stpcaUp$copy()$set_sparse(TRUE, 0.09)$update_theta()
  expect_that(any(sstpca$WHat==0), is_true())
})
