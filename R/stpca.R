#' Structured PCA Model
#'
#' @import Matrix
#' @import foreach
#' @import maxLik
#' @export
StpcaModel <- setRefClass("StpcaModel",
  fields = list(
    X        = "matrix",
    n        = "integer",
    k        = "integer",
    d        = "integer",
    beta     = "numeric",
    K        = "Matrix",
    KD       = "list",
    locs     = "matrix",
    covFn    = "function",
    covFnD   = "function",
    muHat    = "numeric",
    WHat     = "matrix",
    sigSqHat = "numeric",
    Vmean    = "matrix",
    Vvar     = "list",
    logEvidence   = "numeric",
    logPosteriors = "numeric",
    logEvidenceD  = "numeric",
    H        = "list",
    maxim    = "maxim",
    thetaConv = "logical",
    betaHist  = "matrix",
    convergence = "data.frame",
    sparse   = "logical",
    b        = "numeric"
  ),
  methods = list(
    initialize = function(X=matrix(nrow=0, ncol=0), k=1, beta0=numeric(0),
                          locs=matrix(), covFn=function() NULL,
                          covFnD=function() NULL, maxit=500, sparse=FALSE,
                          b=-Inf, ...) {
      X      <<- X
      n      <<- nrow(X)
      k      <<- as.integer(k)
      d      <<- ncol(X)
      locs   <<- locs
      covFn  <<- covFn
      covFnD <<- covFnD
      sparse <<- sparse
      b      <<- b
      betaHist <<- matrix(nrow=0, ncol=length(beta0))
      convergence <<- data.frame(
        logEvidence=numeric(0),
        logPosterior=numeric(0),
        updateStep=factor(c(), levels="beta", "theta")
      )

      matrix(nrow=0, ncol=3,
        dimnames=list(NULL, c('logEvidence', 'logPosterior', 'updateStep')))

      emptyMaxim <- list()
      class(emptyMaxim) <- "maxim"
      maxim  <<- emptyMaxim

      if (nrow(X)>0) {
        thetaInit <- initialize_from_ppca(X, k)
        muHat    <<- thetaInit$mu
        sigSqHat <<- thetaInit$sigSq
        WHat     <<- thetaInit$W
        set_beta(beta0)
        update_theta(maxit)
      }

      callSuper(...)
    },
    update_theta = function(maxit=500, bftol=1e-5) {
      tryCatch({
        vals <- theta_EM(X, WHat, muHat, sigSqHat, K, maxit=maxit,
                         bftol=bftol, sparse=sparse, b=b)
      }, error = function(err) {
        err$message <- paste0("Error in updating theta:\n", err$message)
        stop(err)
      })

      # Theta-related variables to update
      WHat          <<- vals$WHat
      muHat         <<- vals$muHat
      sigSqHat      <<- vals$sigSqHat
      Vmean         <<- vals$Vmean
      Vvar          <<- vals$Vvar
      logPosteriors <<- vals$logPosteriors
      thetaConv     <<- TRUE
      if (!sparse) {
        logEvidence   <<- vals$logEvidence
        H             <<- vals$H
      }

      invisible(.self)
    },
    update_beta = function(...) {
      "Method docs go here"

      if (sparse) stop(paste("Cannot tune beta under a sparse model since",
        "the laplace approximation does not apply. This will probably",
        "never be implemented. Sorry!"))

      # If beta begins on the boundary of the feasible region, shift it
      # a small amount within the region.
      if (!is.null(list(...)$constr)) {
        A <- list(...)$constr$ineqA
        B <- list(...)$constr$ineqB
        for (i in 1:nrow(A)) {
          if (A[i,]%*%beta + B[i] == 0) {
            beta <<- beta + A[i,]/crossprod(A[i,])*1e-4
          }
        }
      }

      maxFn    <- function(beta_) set_beta(beta_)$logEvidence
      maxFnD   <- function(beta_) set_beta(beta_)$logEvidenceD
      betaMax  <- maxLik(maxFn, grad=maxFnD, start=beta, method="bfgs", ...)
      set_beta(coef(betaMax))
      maxim   <<- betaMax #TODO: Is the class of maxim "maxlik" or "maxim"?

      invisible(.self)
    },
    set_beta = function(betaNew) {
      # Only recalculate if neccesary (beta has changed, or theta has changed
      # since last time beta was set).
      if (!identical(beta, betaNew) || thetaConv) {

        # Sanity checks; error if bad
        K_ <- NULL
        try(K_ <- covFn(locs, beta=betaNew), silent=TRUE) # As dppMatrix?
        if (any(!is.finite(betaNew) | !is.numeric(betaNew))) {
            stop(paste0("Bad beta proposed:", paste(betaNew, collapse=',')))
        }
        if (is.null(K_)) stop("Could not construct K with the given beta & locs")
        if (any(!is.finite(K_@x))) stop("K contains non-finite values")
        if (nnzero(K_)==0) stop("K is entirely 0's!")
        if (is(K_, "denseMatrix")) {
          out <- try(K_ <- as(K_, "dppMatrix"))
          if (is(out, "try-error")) stop(paste("Could not cast K to a",
            "dppMatrix as it is semidefinite"))
            # TODO: #"dppMatrix as it is semidefinite. Beta=",paste(round(beta,2),
            #collapse=','))
        }

        # Attach beta to K as an attribute
        attr(K_, "beta") <- betaNew

        # Checks are passed; start assigning variables
        beta <<- betaNew
        K    <<- K_
        thetaConv <<- FALSE

        KD <<- covFnD(locs, beta=beta)
        H  <<- compute_H(X, WHat, muHat, sigSqHat, K)
        logEvidence  <<- log_evidence(X, K, WHat, muHat, sigSqHat, H)
        logEvidenceD <<- log_evidence_d(X, K, WHat, muHat,
                                        sigSqHat, beta, KD, H)
        logPosteriors <<- log_likelihood(X, WHat, muHat, sigSqHat) +
                          log_prior(K, WHat) - logEvidence
      }
      invisible(.self)
    },
    set_b = function(bNew) {
      if (!sparse) stop("This is not a sparse model; 'b' does not apply")
      stopifnot(bNew>0)
      b <<- bNew
      invisible(.self)
    },
    set_sparse = function(spNew, b) {
      stopifnot(is.logical(spNew))
      sparse <<- spNew
      if (spNew) {
        H <<- list()
        logEvidence <<- numeric(0)
        set_b(b)
      }
      invisible(.self)
    },
    update = function(tune.maxit=10, tune.tol=1e-5, EM.maxit=500, EM.bftol=1e-5, ...) {
      for (iter in seq_len(tune.maxit)) {
        # Beta-update
        update_beta(...)
        betaHist <<- rbind(betaHist, beta)
        newConvRow = data.frame(
          'logEvidence'=logEvidence,
          'logPosterior'=tail(logPosteriors, 1),
          'updateStep'='beta')
        convergence <<- rbind(convergence, newConvRow)

        # Theta-update
        update_theta(maxit=EM.maxit, bftol=EM.bftol)
        newConvRow = data.frame(
          'logEvidence'=logEvidence,
          'logPosterior'=tail(logPosteriors, 1),
          'updateStep'='theta')
        convergence <<- rbind(convergence, newConvRow)

        if (iter>=1) {
          lpDiff = abs(Reduce('-', tail(convergence$logPosterior, 2)))
          leDiff = abs(Reduce('-', tail(convergence$logEvidence, 2)))
          if (lpDiff < tune.tol & leDiff < tune.tol) break
        }
      }
      invisible(.self)
    },
    simulate = function(n=1, Wknown=TRUE) {
      Vnew = matrix(rnorm(n*k), nrow=n, ncol=k)
      if (Wknown) { # Simulate from likelihood
        Xnew = sweep(tcrossprod(Vnew, WHat), 2, muHat, "+")
      } else { # Simulate from exact marginal likelihood
        R = chol(K)
        Xnew = vapply(seq_len(n), function(i) {
          W = R %*% matrix(rnorm(d*k), nrow=d, ncol=k)
          as.numeric(W%*%Vnew[i,]) + muHat
        }, numeric(d))
      }
      Xnew = Xnew + rnorm(n*d, sd=sqrt(sigSqHat))
      return(list(X=Xnew, V=Vnew))
    },
    encode = function(Xnew) {
      stopifnot(ncol(Xnew)==d)
      Xnew <- sweep(Xnew, 2, muHat) # Center
      Vnew <- EM.E(Xnew, WHat, sigSqHat)
      return(Vnew) # Stochastic encoder: contains Vmean and Vvar
    },
    decode = function(Vnew) {
      stopifnot(ncol(Vnew)==k)
      Xrec <- tcrossprod(Vnew, WHat) + t(replicate(nrow(Vnew), muHat))
      return(Xrec) # Deterministic decoder, only contains mean because cov mat is big.
    },
    set_X = function(Xnew) {
      stopifnot(ncol(Xnew) == d)
      X <<- Xnew
      n <<- nrow(Xnew)
      muHat <<- colMeans(X)
      invisible(.self)
    },
    crossvalidate = function(nFolds=3, nThreads=1) {
      stopifnot(nFolds>=2)
      stopifnot(nThreads>0)

      # Set up threading if neccesary
      `%doOp%` <- `%do%`
      if (nThreads>1) {
        registerDoMC(nThreads)
        `%doOp%` <- `%dopar%`
      }

      # Logical indices for training set partitions
      trInd <- lapply(1:nFolds, function(fold) {
        ((((1:nrow(X))-1) %% nFolds)+1) != fold
      })

      cvout <- (
        foreach(fold=1:nFolds, .combine=rbind)
      ) %doOp% {
        # Split data into train/test
        Xtr <- X[ trInd[[fold]],]
        Xte <- X[!trInd[[fold]],]

        # Fit model to training data; init from current model
        stpcaTr <- .self$copy()$set_X(Xtr)$update_theta()

        # Compute likelihood of test set given trained model
        llTe <- log_likelihood(Xte, stpcaTr$WHat,
                               stpcaTr$muHat, stpcaTr$sigSqHat)

        foldout <- data.frame(fold=fold, ll=llTe)
        return(foldout)
      }
      return(cvout)
    }# \methods
  )
) # \class def
