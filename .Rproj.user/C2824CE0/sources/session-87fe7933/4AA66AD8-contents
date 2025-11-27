index_max <- function(y, mod) {
  ina <- mod$cluster
  centers <- mod$centers
  ni <- tabulate(ina)
  k <- length(ni)
  ind <- numeric(24)
  n <- dim(y)[1]

  ## "Calinski_Harabasz"
  m <- centers * ni
  m2 <- rowsum(y^2, ina)
  com <- m2 - m^2/ni
  sw <- sum( com )
  s <- Rfast::colVars(y) * ( n - 1)
  st <- sum(s)
  ind[2] <- (n - k ) / (k - 1) * (st - sw)/sw

  ## "Ratkowsky_Lance"
  sw <- Rfast::colsums( com )
  ind[22] <- sqrt( mean( 1 - sw / s ) / k )

  suppressWarnings({
  ## "Ksq_DetW"
  Wg <- 0
  for ( i in 1:k )  Wg <- Wg + Rfast::cova(y[ina == i, ]) * ( ni[i] - 1)
  ind[19] <- k^2 * det(Wg)
  })
  ## "Trace_WiB"
  suppressWarnings({
  Bg <- Rfast::cova(y) * (n - 1) - Wg
  ind[24] <- sum( diag( solve(Wg, Bg) ) )
  })

  dbh <- dpbm <- D1k <- D2k <- D3k <- numeric(k)
  for ( i in 1:k )  {
    yi <- y[ina == i, ]
    dic <- Rfast::dista(yi, centers[i, , drop = FALSE] )
    dbh[i] <- mean(dic^2)
    dpbm[i] <- sum(dic)
    di <- Rfast::Dist(yi, result = "vector")
    D1k[i] <- max(di)
    D2k[i] <- mean(di) * 0.5
    D3k[i] <- 2 * dpbm[i] / ni[i]
  }

  ## Ball_Hall
  ind[1] <- sum(dbh) / k

  ## "PBM" not always the same
  ew <- sum(dpbm)
  et <- sum( Rfast::dista(y, matrix( Rfast::colsums(centers * ni)/n, nrow = 1) ) )
  db <- max( Rfast::Dist(centers) )
  ind[20] <- ( et/ew * db / k)^2

  d4 <- Rfast::Dist(centers)
  d4[lower.tri(d4, diag = TRUE)] <- NA
  xi <- Rfast::Dist(y[ina == 1, ], result = "vector")
  dmax <- max( yi )
  sw <- sum(yi)
  for ( i in 2:k ) {
    yi <- Rfast::Dist(y[ina == i, ], result = "vector")
    dmax <- max(dmax, yi)
    sw <- sw + sum(yi)
  }

  d1 <- d2 <- d3 <- d5 <- matrix(nrow = k, ncol = k)
  sb <- 0
  dmin <- dmax
  for ( i in 1:(k - 1) ) {
    yi <- y[ina == i, ]
    for ( j in (i + 1):k ) {
      dij <- Rfast::dista(yi, y[ina == j, ])
      d1[i, j] <- min(dij)
      d2[i, j] <- max(dij)
      d3[i, j] <- mean(dij)
      d5[i, j] <- (D3k[i] * 0.5 * ni[i] + D3k[j] * 0.5 * ni[j] ) / (ni[i] + ni[j])
      sb <- sb + sum(dij)
      dmin <- min( dmin, min( dij ) )
    }
  }

  ## "Point_biserial" not always the same
  nw <- 0.5 * ( sum(ni^2) - n )
  nb <- 0.5 * n * (n - 1) - nw
  ind[21] <- (sw/nw - sb/nb) * sqrt(nw * nb) / (nw + nb)

  ## "Dunn" not always the same
  ind[3] <- dmin / dmax

  ## GDI
  ind[4] <- min(d1, na.rm = TRUE) / max(D1k)
  ind[5] <- min(d1, na.rm = TRUE) / max(D2k)
  ind[6] <- min(d1, na.rm = TRUE) / max(D3k)
  ind[7] <- min(d2, na.rm = TRUE) / max(D1k)
  ind[8] <- min(d2, na.rm = TRUE) / max(D2k)
  ind[9] <- min(d2, na.rm = TRUE) / max(D3k)
  ind[10] <- min(d3, na.rm = TRUE) / max(D1k)
  ind[11] <- min(d3, na.rm = TRUE) / max(D2k)
  ind[12] <- min(d3, na.rm = TRUE) / max(D3k)
  ind[13] <- min(d4, na.rm = TRUE) / max(D1k)
  ind[14] <- min(d4, na.rm = TRUE) / max(D2k)
  ind[15] <- min(d4, na.rm = TRUE) / max(D3k)
  ind[16] <- min(d5, na.rm = TRUE) / max(D1k)
  ind[17] <- min(d5, na.rm = TRUE) / max(D2k)
  ind[18] <- min(d5, na.rm = TRUE) / max(D3k)

  ## "Silhouette"
  ind[23] <- mean( Rfast2::silhouette(y, ina)$stats[,4] )

  names(ind) <- c( "Ball_Hall", "Calinski_Harabasz", "Dunn", "GDI11", "GDI12", "GDI13", "GDI21",
                   "GDI22", "GDI23", "GDI31", "GDI32", "GDI33", "GDI41", "GDI42", "GDI43", "GDI51",
                   "GDI52", "GDI53", "Ksq_DetW", "PBM", "Point_biserial", "Ratkowsky_Lance",
                   "Silhouette", "Trace_WiB" )
  ind
}



















