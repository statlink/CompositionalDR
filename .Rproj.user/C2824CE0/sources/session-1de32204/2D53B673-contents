clust.plot <- function(mod, x) {
  x <- Rfast::Standardize(x)
  factoextra::fviz_cluster( mod, data = x, palette = "jco", geom = "point",
                            ellipse.type = "convex", ggtheme = theme_minimal() )
}
