saa <- function(x, k, lr_w = 0.1, lr_h = 0.1, maxiter = 10000, tol = 1e-5) {
  # x: n x D compositional data (rows sum to 1)
  # k: number of archetypes
  # lr_w, lr_h: EG learning rates for W and H
  # returns W (n x k), H (k x D), Z = W %*% H, objective trace info

  n <- dim(x)[1]  ;   D <- dim(x)[2]
  # initialize W and H randomly and normalize rows
  W <- matrix( Rfast2::Runif(n * k), nrow = n, ncol = k)
  W <- W / Rfast::rowsums(W)
  H <- matrix( Rfast2::Runif(k * D), nrow = k, ncol = D)
  H <- H / Rfast::rowsums(H)

  prev_obj <- Inf
  obj <- Inf

  for ( iter in 1:maxiter ) {
    Z <- W %*% H           # n x D
    E <- Z - x             # n x D
    grad_w <- tcrossprod(E, H)
    W <- W * exp(- lr_w * grad_w)
    W <- W / Rfast::rowsums(W) # +eps?
    Z <- W %*% H
    E <- Z - x
    grad_h <- crossprod(W, E)
    H <- H * exp(- lr_h * grad_h)
    H <- H / Rfast::rowsums(H) ## +eps?
    Z <- W %*% H
    recon_err <- sum( (x - Z)^2 )
    obj <- recon_err

    # convergence check
    if ( abs(prev_obj - obj) < tol ) {
      cat("Converged at iter", iter, "obj =", obj, "recon_err =", recon_err, "\n")
      break
    }
    prev_obj <- obj
  }

  list(W = W, H = H, Z = Z, recon_err = recon_err, iters = iter)
}


