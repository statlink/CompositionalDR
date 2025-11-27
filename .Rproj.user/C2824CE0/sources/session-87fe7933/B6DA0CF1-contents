alfa.cikmeans <- function(x, ncl = 10, trim = 0, a = seq(-1, 1, by = 0.1), max.iters = 50, nstart = 10) {

  if ( min(x) == 0 )  a <- a[a > 0]
  la <- length(a)
  best_min <- min_crit <- matrix(nrow = 9, ncol = la)
  best_max <- max_crit <- matrix(nrow = 24, ncol = la)
  colnames(best_min) <- colnames(min_crit) <- paste("alpha=", a, sep = "")
  colnames(best_max) <- colnames(max_crit) <- paste("alpha=", a, sep = "")

  for ( i in 1:la ) {
    y <- Compositional::alfa(x, a[i])$aff
    y <- Rfast::standardise(y)
    mod <- cikmeans(y, ncl = ncl, trim = trim, max.iters = max.iters, nstart = nstart)
    best_min[, i] <- mod$best_min
    min_crit[, i] <- Rfast::colMins(mod$min_crit, TRUE)
    best_max[, i] <- mod$best_max
    max_crit[, i] <- Rfast::colMaxs(mod$max_crit, TRUE)
  }
  rownames(best_min) <- rownames(min_crit) <- names(mod$best_min)
  rownames(best_max) <- rownames(max_crit) <- names(mod$best_max)

  id <- Rfast::rowMins(min_crit)
  idmin <- numeric( dim(min_crit)[1] )
  for ( i in 1:dim(min_crit)[1] )  idmin[i] <- best_min[i, id[i]]
  minim <- cbind( a[id], idmin, Rfast::rowMins(min_crit, TRUE) )
  colnames(minim) <- c("alpha", "clusters", "value")
  rownames(minim) <- rownames(best_min)

  id <- Rfast::rowMaxs(max_crit)
  idmax <- numeric( dim(max_crit)[1] )
  for ( i in 1:dim(max_crit)[1] )  idmax[i] <- best_max[i, id[i]]
  maxim <- cbind( a[id], idmax, Rfast::rowMaxs(max_crit, TRUE) )
  colnames(maxim) <- c("alpha", "clusters", "value")
  rownames(maxim) <- rownames(best_max)

  list(minim = minim, maxim = maxim)
}


