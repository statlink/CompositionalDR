# Dirichlet Component Analysis (DCA)
# Based on Wang et al. (2008)
dca <- function(x, k, maxiter = 1000, pop_size = 50) {
  M <- dim(x)[1]  ;  N <- dim(x)[2]
  # Initialize population of balanced rearrangement matrices
  population <- list()
  for ( i in 1:pop_size )  population[[ i ]] <- .random_balanced_rearrangement(N, k)

  best_R <- NULL
  best_alpha <- Inf
  # Genetic algorithm
  for ( iter in 1:maxiter ) {
    # Evaluate fitness for each individual
    alphas <- numeric(pop_size)
    for ( i in 1:pop_size ) {
      R <- population[[ i ]]
      # Apply transformation
      xtrans <- tcrossprod(x, R)
      # Apply regularization
      xreg <- .regularize_data(xtrans)
      # Estimate Dirichlet correlation
      alphas[i] <- .diri_a0(xreg)
    }

    # Track best solution
    min_idx <- which.min(alphas)
    if ( alphas[min_idx] < best_alpha ) {
      best_alpha <- alphas[min_idx]
      best_R <- population[[ min_idx ]]
    }
    # Check convergence
    if ( iter > 10 && sd(alphas) < 1e-6 ) {
      break
    }

    # Selection and crossover
    fitness <- .compute_fitness(alphas)
    fitness[fitness == 0] <- 1e-10  # Avoid zero probabilities
    probs <- fitness / sum(fitness)
    # Create next generation
    new_population <- list()
    new_population[[ 1 ]] <- best_R
    # Reduce population size
    current_size <- max( 5, pop_size - floor(iter / 10) )

    for ( i in 2:current_size ) {
      # Sample two parents
      parents <- sample( 1:pop_size, 2, prob = probs, replace = TRUE )
      R1 <- population[[ parents[1] ]]
      R2 <- population[[ parents[2] ]]
      # Crossover: weighted average
      weight <- fitness[ parents[1] ] / ( fitness[ parents[1] ] + fitness[ parents[2] ] )
      R_new <- weight * R1 + (1 - weight) * R2
      new_population[[ i ]] <- R_new
    }

    population <- new_population
    pop_size <- length(population)
  }

  # Return results
  R <- t( best_R )
  Z <- x %*% R
  Z <- .regularize_data(Z)
  list( R = R, alpha = best_alpha, iters = iter, Z = Z )
}



# Function to estimate Dirichlet precision parameter (alpha_0) using Newton-Raphson
.diri_a0 <- function(x, tol = 1e-6) {
  n <- dim(x)[1]   ;   D <- dim(x)[2]
  a <- 0   ;   ea <- 1
  slx <- sum( Rfast::Log(x) )
  lik1 <- n * ( lgamma(D * ea) - D * lgamma(ea) ) + (ea - 1) * slx
  grad <- n * D * ea * ( digamma(D * ea) - digamma(ea) ) + ea * slx
  hess <- n * trigamma(D * ea) * (D * ea)^2 - n * D * trigamma(ea) * ea^2 + grad
  a <- a - grad/hess
  ea <- exp(a)
  lik2 <- n * ( lgamma(D * ea) - D * lgamma(ea) ) + (ea - 1) * slx

  while ( lik2 - lik1 > tol ) {
    lik1 <- lik2
    grad <- n * D * ea * ( digamma(D * ea) - digamma(ea) ) + ea * slx
    hess <- n * trigamma(D * ea) * (D * ea)^2 - n * D * trigamma(ea) * ea^2 + grad
    a <- a - grad/hess
    ea <- exp(a)
    lik2 <- n * ( lgamma(D * ea) - D * lgamma(ea) ) + (ea - 1) * slx
  }
  ea
}

# Function to apply regularization operator
.regularize_data <- function(X) {
  # X: M x N matrix of compositional data
  # Returns regularized data
  delta <- min(X)
  if ( delta >= 0 )  return(X)  # No regularization needed
  X_reg <- X - delta
  # Radial projection back to simplex
  X_reg <- X_reg / Rfast::rowsums(X_reg)
  return(X_reg)
}

# Function to check if matrix is a balanced rearrangement
.is_balanced_rearrangement <- function(R, N, K, tol = 1e-6) {
  # Check non-negativity
  if ( any(R < -tol) )  return(FALSE)
  # Check column sums = 1
  col_sums <- Rfast::colsums(R)
  if ( any( abs(col_sums - 1) > tol) )  return(FALSE)
  # Check row sums = N/K
  row_sums <- Rfast::rowsums(R)
  if ( any( abs(row_sums - N/K) > tol ) )  return(FALSE)
  return(TRUE)
}

# Function to generate random balanced rearrangement matrix
.random_balanced_rearrangement <- function(N, K) {
  # Initialize with uniform probabilities
  R <- matrix( Rfast2::Runif(K * N), K, N )
  # Make it satisfy constraints using iterative proportional fitting
  for ( iter in 1:1000 ) {
    # Normalize columns to sum to 1
    R <- Rfast::eachrow(R, Rfast::colsums(R), oper="/")
    # Normalize rows to sum to N/K
    R <- R / Rfast::rowsums(R) * (N / K)
    if ( .is_balanced_rearrangement(R, N, K) )  break
  }
  R
}

# Fitness function for genetic algorithm
.compute_fitness <- function(alpha_values) {
  median_alpha <- median(alpha_values)
  -log( pmin(alpha_values / median_alpha, 1) )
}





