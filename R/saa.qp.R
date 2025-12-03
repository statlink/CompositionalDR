saa.qp <- function(x, k, maxiter = 1000, tol = 1e-6, ridge = 1e-8) {
  n <- dim(x)[1]  ;  D <- dim(x)[2]
  # Initialize W and H randomly and normalize to simplex
  W <- matrix(Rfast2::Runif(n * k), nrow = n, ncol = k)
  W <- W / Rfast::rowsums(W)
  H <- matrix(Rfast2::Runif(k * D), nrow = k, ncol = D)
  H <- H / Rfast::rowsums(H)
  prev.obj <- NULL
  # Setup for H update (each column on simplex)
  Aeq_H   <- rep(1, k)
  Aineq_H <- diag(k)
  Amat_H  <- cbind(Aeq_H, Aineq_H)
  bvec_H  <- c(1, rep(0, k))
  meq     <- 1
  ridgeK <- ridge * Aineq_H
  # Setup for W update (each row on simplex)
  Aeq_W   <- rep(1, k)
  Aineq_W <- diag(k)
  Amat_W  <- cbind(Aeq_W, Aineq_W)
  bvec_W  <- c(1, rep(0, k))

  for ( it in 1:maxiter ) {
    ## Update H: each COLUMN on simplex
    WtW <- crossprod(W)           # k x k
    WtX <- crossprod(W, x)        # k x D
    DmatH <- 2 * WtW + ridgeK
    for ( j in 1:D ) {
      dvec <- 2 * WtX[, j]        # gradient for column d
      sol <- try( quadprog::solve.QP( Dmat = DmatH, dvec = dvec, Amat = Amat_H,
                                      bvec = bvec_H, meq = meq)$solution, silent = TRUE )
      if ( inherits(sol, "try-error") || any( !is.finite(sol) ) )  sol <- rep(1/k, k) # fallback to barycenter
      H[, j] <- sol               # update column d
    }
    ## Update W: each row on simplex
    HHt <- tcrossprod(H)
    DmatW <- 2 * HHt + ridgeK
    for ( i in 1:n ) {
      dvec <- 2 * drop(H %*% x[i, ])
      sol <- try( quadprog::solve.QP( Dmat = DmatW, dvec = dvec, Amat = Amat_W,
                                      bvec = bvec_W, meq = meq)$solution, silent = TRUE )
      if ( inherits(sol, "try-error") || any( !is.finite(sol) ) )  sol <- rep(1/k, k)
      W[i, ] <- sol
    }
    ## Check convergence
    obj <- sum( (x - W %*% H)^2 )
    if ( is.null(prev.obj) ) {
      relchg <- Inf
    } else {
      relchg <- abs(prev.obj - obj) / (1 + prev.obj)
      if ( !is.finite(relchg) )  relchg <- Inf
    }

    if ( relchg < tol )  break
    prev.obj <- obj
  }

  list(W = W, H = H, Z = W %*% H, obj = obj, iters = it)
}
