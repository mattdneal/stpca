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

cholSolve <- function(R, X) {
  # Where R is the pivoted cholesky decomposition of the matrix A, solve
  # the system A^{-1}X
  # TODO: Don't transpose R (use flags in forward/backsolve)
  pivot = attr(R, "pivot")
  unpivot = order(pivot)
  if (is.vector(X) | ncol(X)==1) {
    working = matrix(forwardsolve(t(R), X[unpivot,]), ncol=1)
    return(matrix(backsolve(R, working), ncol=1)[pivot,,drop=FALSE])
  } else {
    return(backsolve(R, forwardsolve(t(R), X[pivot,]))[unpivot,])
  }
}
