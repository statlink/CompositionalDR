cikmeans <- function(y, ncl = 10, trim = 0, max.iters = 50, nstart = 10, all = FALSE) {

  minc <- c( "Banfeld_Raftery", "Davies_Bouldin", "Det_Ratio", "Log_Det_Ratio",
             "Log_SS_Ratio", "McClain_Rao", "Ray_Turi", "Scott_Symons", "Xie_Beni")

  maxc <- c( "Ball_Hall", "Calinski_Harabasz", "Dunn", "GDI11", "GDI12", "GDI13", "GDI21",
             "GDI22", "GDI23", "GDI31", "GDI32", "GDI33", "GDI41", "GDI42", "GDI43", "GDI51",
             "GDI52", "GDI53", "Ksq_DetW", "PBM", "Point_biserial", "Ratkowsky_Lance",
             "Silhouette", "Trace_WiB")

  min_mat <- matrix(nrow = ncl - 1, ncol = 9)
  max_mat <- matrix(nrow = ncl - 1, ncol = 24)

  colnames(min_mat) <- minc
  rownames(min_mat) <- 2:ncl
  colnames(max_mat) <- maxc
  rownames(max_mat) <- 2:ncl

  p <- dim(y)[2]

  if ( all ) {
    cl_mat <- matrix(nrow = dim(y)[1], ncol = ncl - 1)
    for ( i in 2:ncl ) {
      if ( trim == 0 ) {
        mod <- kmeans(y, i, iter.max = max.iters, nstart = nstart )
      } else {
        m <- lowmemtkmeans::tkmeans(y, i, alpha = trim, iter = max.iters, nstart = nstart )
        cluster <- as.vector( lowmemtkmeans::nearest_cluster(y, m) )
	      mod <- list()
	      mod$centers <- m
        mod$cluster <- cluster
      }
      if ( min( tabulate(mod$cluster) ) > p ) {
        min_mat[i - 1, ] <- index_min(y, mod)
        max_mat[i - 1, ] <- index_max(y, mod)
      }
      cl_mat[, i - 1] <- mod$cluster
    }
    best_min <- Rfast::colMins(min_mat) + 1
    names(best_min) <- minc
    best_max <- Rfast::colMaxs(max_mat) + 1
    names(best_max) <- maxc
    result <- list( "min_crit" = min_mat, "best_min" = best_min,
                    "max_crit" = max_mat, "best_max" = best_max, "cluster" = cl_mat )

  } else {
    for ( i in 2:ncl ) {
      if ( trim == 0 ) {
        mod <- kmeans(y, i, iter.max = max.iters, nstart = nstart )
      } else {
        m <- lowmemtkmeans::tkmeans(y, i, alpha = trim, iter = max.iters, nstart = nstart )
        cluster <- as.vector( lowmemtkmeans::nearest_cluster(y, m) )
        mod <- list()
	      mod$centers <- m
        mod$cluster <- mod$cluster
      }
      if ( min( tabulate(mod$cluster) ) > p ) {
        min_mat[i - 1, ] <- index_min(y, mod)
        max_mat[i - 1, ] <- index_max(y, mod)
      }
    }
    best_min <- Rfast::colMins(min_mat) + 1
    names(best_min) <- minc
    best_max <- Rfast::colMaxs(max_mat) + 1
    names(best_max) <- maxc
    result <- list( "min_crit" = min_mat, "best_min" = best_min,
                    "max_crit" = max_mat, "best_max" = best_max )
  }

  result
}
