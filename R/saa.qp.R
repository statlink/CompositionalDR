saa.qp <- function(x, k, maxiter = 10000, tol = 1e-6, ridge = 1e-6) {
  n <- dim(x)[1]  ;  D <- dim(x)[2]
  # Initialize W,H on simplex
  W <- matrix( Rfast2::Runif(n * k), nrow = n, ncol = k)
  W <- W / Rfast::rowsums(W)
  H <- matrix( Rfast2::Runif(k * D), nrow = k, ncol = D)
  H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
  error <- numeric(maxiter)

  for ( it in 1:maxiter ) {
    # ----------------- W update -----------------
    G_W <- 2 * tcrossprod(H) + diag(ridge, k)
    A_W <- rbind(rep(1, k), diag(k))
    b_W <- c(1, rep(0, k))

    for (i in 1:n) {
      g_W <- 2 * (H %*% X[i, ])
      sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W, Amat = t(A_W), bvec = b_W, meq = 1)
      W[i, ] <- sol$solution
    }
    # ----------------- H update -----------------
    WH <- W %*% H
    A_H <- rbind(rep(1, D), diag(D))
    b_H <- c(1, rep(0, D))

    for ( j in 1:k ) {
      # Residual without row j's contribution
      WH_minus_j <- WH - W[, j] %*% H[j, , drop = FALSE]
      Rj <- X - WH_minus_j
      # Weighted least squares with weights W[,j]^2
      w3 <- W[, j]^3    ;  w4 <- W[, j]^4
      # Linear term: 2 * sum_i W[i,j]^3 * Rj[i,d] for each d
      g_H <- 2 * Rfast::colsums(w3 * Rj)
      # Quadratic term: 2 * (sum_i W[i,j]^4) * I_D
      G_H <- 2 * sum(w4) * diag(D) + diag(ridge, D)
      sol <- quadprog::solve.QP(Dmat = G_H, dvec = g_H, Amat = t(A_H), bvec = b_H, meq = 1)
      H[j, ] <- sol$solution
    }
    # ----------------- Error -----------------
    err <- sum( (X - W %*% H)^2 )
    error[it] <- err
    if ( it > 1 && abs(error[it - 1] - err) < tol ) {
      break
    }
  }

  list(W = W, H = H, Z = W %*% H, obj = errors[it], iters = it)
}
