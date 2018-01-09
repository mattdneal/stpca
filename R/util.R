#' @import Matrix
#' @import fields
distanceMatrix <- function(X, X2=NA, max.dist=Inf) {
  stopifnot(!(all(is.na(nrow(X)))))

  stopifnot(all(is.finite(X)))

  symmetric=FALSE
  if (missing(X2) || all(is.na(X2))) {
    symmetric=TRUE
    X2 = X
    d = nrow(X)
  } else {
    d = max(nrow(X), nrow(X2))
  }

  if (max.dist == Inf) {
    if (symmetric) {
      D = Matrix(rdist(X), sparse=FALSE)
    } else {
      D = Matrix(rdist(X, X2), sparse=FALSE)
    }
  } else {
    # Take a guess at the maximum number of points needed. Calculate the volume
    # of data-space neighboring a single point. The maximum number of ponts is
    # nPoints x (nPoints x neighboring.volume) = nPoints x (expected neighbours)
    d_s        = ncol(X) # Spatial dimensions
    nCubeVol   = (2*max.dist)^d_s
    max.points = ceiling(nCubeVol*d)*d*2
    #max.points = max(ceiling(nCubeVol*d)*d*2, nrow(X) * nrow(X2)) # TODO: Is this better?

    # Use fields.rdist.near to find all distances < max.dist, stored as triples
    continue = FALSE
    while (!continue) {
      Dtriples = try(fields.rdist.near(X, X2, delta=max.dist, max.points = max.points), silent=TRUE)
      if (inherits(Dtriples, "try-error")) {
        if (grepl("Ran out of space, increase max.points", gettext(Dtriples), ignore.case=TRUE)) {
          max.points = max.points * 2
        } else {
          stop(Dtriples)
        }
      } else {
         continue=TRUE
      }
    }

    if (length(Dtriples$ra) == nrow(X) * nrow(X2)) {
      # Covariance matrix is dense: all points lie within maximum distance
      if (symmetric) {
        return(Matrix(rdist(X), sparse=FALSE))
      } else {
        return(Matrix(rdist(X, X2), sparse=FALSE))
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

getij <- function(A) {
  stopifnot(is(A, "Matrix"))

  if (is(A, "denseMatrix")) {
    grid = expand.grid(1:nrow(A), 1:ncol(A))
    return(list(i=grid$Var1, j=grid$Var2))
  }

  if (.hasSlot(A, "i")) {
    i = A@i+1
  } else {
    i = findInterval(seq(A@x)-1,A@p[-1])+1
  }

  if (.hasSlot(A, "j")) {
    j = A@j+1
  } else {
    j = findInterval(seq(A@x)-1,A@p[-1])+1
  }

  return(list(i=i, j=j))
}
