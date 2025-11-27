clust.plot <- function(mod, x) {
  x <- Rfast::standardise(x)
  factoextra::fviz_cluster( mod, data = x, palette = "jco", geom = "point",
                            ellipse.type = "convex")
}
