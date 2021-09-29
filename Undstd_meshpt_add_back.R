X <- subset(Voronoi$dirsgs, (ind1 == 2 | ind2 == 2) )
X
X <- matrix(c(X$x1, X$x2, X$y1, X$y2, X$bp1, X$bp2), ncol = 3)
X <- unique(X)


if (sum(X[, 3]) > 0) {
  # because if the pt being triangulated is on the boundary of the convex hull,
  # the region will be unbounded, 
  # them won't be able to calculate the area
  
  # The region Vi is unbounded iff pi belongs to the boundary of the convex hull of P.
  X <- rbind(X, c(mesh_locs[2, ], 0))
}

plot(X[, 1:2])



X <- subset(Voronoi$dirsgs, (ind1 == 1 | ind2 == 1) )
X

X <- matrix(c(X$x1, X$x2, X$y1, X$y2, X$bp1, X$bp2), ncol = 3)
X

X <- unique(X)
plot(X)
c(mesh_locs[1, ], 0) # [1] -129.58324   36.82639    0.00000
# just a little bit higher than bps, to avoid the unbounded situation 
# happened when the triangulation p is on the boundary. 



