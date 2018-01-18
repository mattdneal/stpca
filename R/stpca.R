#' Structured PCA Model
#'
#' @import Matrix
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
    logEvidences  = "numeric",
    logPosteriors = "numeric",
    logEvidenceD  = "numeric",
    H        = "list",
    maxim    = "maxim",
    thetaConv = "logical"
  ),
  methods = list(
    initialize = function(X=matrix(nrow=0, ncol=0), k=1, beta0=numeric(0),
                          locs=matrix(), covFn=function() NULL,
                          covFnD=function() NULL, maxit=50, ...) {
      X      <<- X
      n      <<- nrow(X)
      k      <<- as.integer(k)
      d      <<- ncol(X)
      locs   <<- locs
      covFn  <<- covFn
      covFnD <<- covFnD

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
    update_theta = function(maxit=50, bftol=1e-5) {
      # TODO: Better convergence criteria (so it can change between outer iters)
      tryCatch({
        vals <- theta_EM(X, WHat, muHat, sigSqHat, K, maxit=maxit, bftol=bftol)
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
      logEvidence   <<- vals$logEvidence
      logEvidences  <<- c(logEvidences, vals$logEvidence)
      logPosteriors <<- vals$logPosteriors
      H             <<- vals$H
      thetaConv     <<- TRUE

      invisible(.self)
    },
    update_beta = function(...) {
      "Method docs go here"

      # TODO: Don't compute K, H twice. Use set_beta.
      #maxFn    <- function(beta_) {
      #  K_ <- as(covFn(locs, beta=beta_), "dppMatrix")
      #  H_ <- compute_H(X, WHat, muHat, sigSqHat, K_)
      #  log_evidence(X, K_, WHat, muHat, sigSqHat, H_)
      #}
      #maxFnD   <- function(beta_) {
      #  K_  <- as(covFn(locs, beta=beta_), "dppMatrix")
      #  KD_ <- covFnD(locs, beta=beta_)
      #  H_  <- compute_H(X, WHat, muHat, sigSqHat, K_)
      #  log_evidence_d(X, K_, WHat, muHat, sigSqHat, beta_, KD_, H_)
      #}
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
        }

        # Checks are passed; start assigning variables
        beta <<- betaNew
        K    <<- K_
        thetaConv <<- FALSE

        KD <<- covFnD(locs, beta=beta)
        H  <<- compute_H(X, WHat, muHat, sigSqHat, K)
        logEvidence  <<- log_evidence(X, K, WHat, muHat, sigSqHat, H)
        logEvidenceD <<- log_evidence_d(X, K, WHat, muHat,
                                        sigSqHat, beta, KD, H)
      }
      invisible(.self)
    },
    update = function(nIterOuter, EM.maxit=50, EM.bftol=1e-5, ...) {
      for (iter in seq_len(nIterOuter)) {
        update_beta(...)
        update_theta(maxit=EM.maxit, bftol=EM.bftol)
      }
      invisible(.self)
    }
  )
)
