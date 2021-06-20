### ------------------
### Find attribute from grid
### ------------------

## g is a data frame of grid of size n x 3
## p is a data frame of points of size n x 2
grid_to_point <- function(g,p) {
  pz <- rep(0,nrow(p))
  for (i in 1:nrow(p)) {
    idx_min <- which.min(abs(p[i,1] - g[,1]) + abs(p[i,2] - g[,2]))
    pz[i] <- g[idx_min,3]
  }
  p[,3] <- pz
  p
}

