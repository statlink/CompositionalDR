index_min <- function(y, mod) {
  ina <- mod$cluster
  centers <- mod$centers
  ni <- tabulate(ina)
  k <- length(ni)
  ind <- numeric(9)
  n <- dim(y)[1]

  ## Banfeld_Raftery
  m <- centers * ni
  m2 <- rowsum(y^2, ina)
  s <- Rfast::rowsums( (m2 - m^2/ni) )
  ind[1] <- sum( ni * log( s / ni ) )

  ## Davies_Bouldin
  dk <- numeric( k )
  for ( i in 1:k )  dk[i] <- mean( Rfast::dista(y[ina == i, ], centers[i, , drop = FALSE] ) )
  Dkk <- Rfast::Dist(mod$centers)
  dkk <- matrix(0, nrow = k, ncol = k)
  for ( i in 1:(k - 1) ) {
    for ( j in (i + 1):k ) {
      dkk[i, j] <- dkk[j, i] <- dk[i] + dk[j]
    }
  }
  a <- dkk / Dkk
  diag(a) <- 0
  ind[2] <- sum( Rfast::rowMaxs(a, TRUE) ) / k

  suppressWarnings({
  ## Det_Ratio
  mat1 <- mat2 <- 0
  for ( i in 1:k ) {
    ssk <- Rfast::cova(y[ina == i, ]) * ( ni[i] - 1)
    mat1 <- mat1 + ssk
    mat2 <- mat2 + ni[i] * determinant( ssk / ni[i], logarithm = TRUE )$modulus
  }
  ss <- Rfast::cova(y)
  ind[3] <- det( ss ) / det( mat1/(n - 1) )

  ## Log_Det_ratio
  ind[4] <- n * log( ind[3] )
  })

  ## "Log_SS_Ratio"
  sw <- sum(s)
  st <- sum( sum( diag(ss) ) * ( n - 1) )
  ind[5] <- log( (st - sw)/sw )

  ## "Ray-Turi"
  dmin <- min( Rfast::Dist(mod$centers, result = "vector") )
  ind[7] <- sw / (n * dmin^2)

  ## "Scott_Symons" not quite the same
  ind[8] <- as.numeric(mat2)

  ## "Xie_Beni"
  dmin <- 1e+6
  for ( i in 1:(k - 1) ) {
    yi <- y[ina == i, ]
    for ( j in (i + 1):k ) {
      dmin <- min( dmin, min( Rfast::dista(yi, y[ina == j, ], square = TRUE) ) )
    }
  }
  ind[9] <- sw / n / dmin

  ## "McClain_Rao" not always the same quantity
  sw <- 0
  for ( i in 1:k )  sw <- sw + Rfast::Dist(y[ina == i, ], result = "sum")
  sb <- 0
  for ( i in 1:(k - 1) ) {
    yi <- y[ina == i, ]
    for ( j in (i + 1):k ) {
      sb <- sb + Rfast::dista(yi, y[ina == j, ], result = "sum")
    }
  }
  nw <- 0.5 * ( sum(ni^2) - n )
  nb <- 0.5 * n * (n - 1) - nw
  ind[6] <- nb/nw * sw/sb

  names(ind) <- c( "Banfeld_Raftery", "Davies_Bouldin", "Det_Ratio", "Log_Det_Ratio",
                   "Log_SS_Ratio", "McClain_Rao", "Ray_Turi", "Scott_Symons", "Xie_Beni" )
  ind
}



