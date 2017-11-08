#' @import Matrix
#' @import fields
distanceMatrix <- function(X, X2=NA, max.dist=Inf, max.points=NA) {
  stopifnot(!(all(is.na(nrow(X)))))

  symmetric=FALSE
  if (missing(X2) || all(is.na(X2))) {
    symmetric=TRUE
    X2 = X
    d = nrow(X)
  } else {
    d = max(nrow(X), nrow(X2))
  }

  if (!is.infinite(max.dist)) {
    # Take a guess at the maximum number of points needed. Calculate the volume
    # of data-space neighboring a single point. The maximum number of ponts is
    #    #points x (#points x neighboring.volume)
    #  =
    #    #points x (expected #neighbours)
    if (is.na(max.points)) {
      d_s        = ncol(X) # Spatial dimensions
      #nBallVol   = pi^(0.5*d_s) * max.dist^d_s / gamma(0.5*d_s+1)
      #max.points = ceiling(nBallVol*d*d)
      nCubeVol   = (2*max.dist)^d_s
      max.points = ceiling(nCubeVol*d)*d*2
      #max.points = max(ceiling(nCubeVol*d)*d*2, nrow(X) * nrow(X2)) # TODO: Is this better?
    }

    # Use fields.rdist.near to find all distances within max.dist, stored as triples
    Dtriples = fields.rdist.near(X, X2, delta=max.dist, max.points = max.points)

    if (length(Dtriples$ra) == nrow(X) * nrow(X2)) {
      # Covariance matrix is dense: all points lie within maximum distance
      if (symmetric) {
        return(Matrix(rdist(X)))
      } else {
        return(Matrix(rdist(X, X2)))
      }
    }

    if (symmetric) {
      upInd = Dtriples$ind[,1] <= Dtriples$ind[,2]
      D = sparseMatrix(i=Dtriples$ind[upInd,1], j=Dtriples$ind[upInd,2],
                       dims=c(d,d), x=Dtriples$ra[upInd], symmetric=TRUE,
                       index1=TRUE)
    } else {
      D = sparseMatrix(i=Dtriples$ind[,1], j=Dtriples$ind[,2],
                       dims=c(nrow(X), nrow(X2)), x=Dtriples$ra, index1=TRUE)
    }
  } else {
    if (symmetric) {
      D = Matrix(rdist(X))
    } else {
      D = Matrix(rdist(X, X2))
    }
  }
  return(D)
}

#' Solve the sylvester equation AW + WB = C for W.
#'
#' @param A d x d positive definite matrix
#' @param B k x k positive definite matrix
#' @param C d x k matrix
#' @return The solution W
#' @import Matrix
sylSolve <- function(A, B, C) {
  Beig = eigen(B, symmetric=TRUE)
  Ctil = C %*% Beig$vectors
  Wtil = vapply(1:ncol(C), function(i_) {
    lhs = A
    Matrix::diag(lhs) = Matrix::diag(lhs) + Beig$values[i_]
    return(as.vector(Matrix::solve(lhs, Ctil[,i_])))
  }, numeric(nrow(C)))
  W = Wtil %*% t(Beig$vectors)
  return(W)
}
