saakl.eg2 <- function(x, k, W = NULL, H = NULL, lr_w = 0.1, maxiter = 10000, tol = 1e-6, clip_exp = 50) {

  runtime <- proc.time()
  n <- dim(x)[1]  ;  D <- dim(x)[2]
  # Initialize W, H randomly on simplex
  if ( is.null(W) ) {
    W <- matrix( Rfast2::Runif(n * k), nrow = n, ncol = k )
    W <- W / Rfast::rowsums(W)
  }
  if ( is.null(H) ) {
    H <- matrix( Rfast2::Runif(k * D), nrow = k, ncol = D )
    H <- H / Rfast::rowsums(H)  
  } 
  prev_obj <- Inf
  Z <- W %*% H              # model

  for ( iter in 1:maxiter ) {
    R <- x / Z                # elementwise
    # ---- Gradient for KL divergence ----
    grad_w <- tcrossprod(- R, H)
    expo_w <- -lr_w * grad_w
    expo_w <- pmax( pmin(expo_w, clip_exp), -clip_exp )
    W <- W * exp(expo_w)
    W <- W / Rfast::rowsums(W)
	mod <- try( Compositional::tflr.irls(x, W), silent = TRUE)
	if ( identical( class(mod), "try-error" ) ) {
      H <- codalm::codalm(x, W)
      obj <- sum( x * log(x / Z), na.rm = TRUE )
	} else {
      H <- mod$be
      obj <- mod$kl
	}  
    Z <- W %*% H              # model
    if ( abs(prev_obj - obj ) < tol) {
      break
    }
    prev_obj <- obj
  }

  runtime <- proc.time() - runtime

  list(W = W, H = H, Z = Z, obj = obj, iters = iter, runtime = runtime)
}