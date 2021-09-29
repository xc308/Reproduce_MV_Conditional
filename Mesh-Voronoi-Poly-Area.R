##======================#
# Mesh-Voronoi-Poly-Area
##======================#

install.packages("deldir")
install.packages("gpclib")

library(deldir)
library(gpclib)


#---------------------#
# From mesh to voronoi
#---------------------#

Voronoi <- deldir(x = mesh_locs[, 1],
       y = mesh_locs[, 2], 
       plot = T, 
       rw = c(min(mesh_locs[, 1]) - 0.0001,
              max(mesh_locs[, 1]) + 0.0001,
              min(mesh_locs[, 2]) - 0.0001,
              max(mesh_locs[, 2]) + 0.0001))


Voronoi$summary$dir.area
# what we want for integration


str(Voronoi)
#  $ dirsgs  :'data.frame':	6052 obs. of  10 variables:
#..$ x1     : num [1:6052] -134 -128 -128 -129 -128 ... end pts of edge of tiles
#..$ y1     : num [1:6052] 36.8 37.1 38.3 37.9 37.6 ...
#..$ x2     : num [1:6052] -132 -128 -129 -129 -128 ...
#..$ y2     : num [1:6052] 39.1 37.1 38.3 37.7 37.7 ...
#..$ ind1   : int [1:6052] 1375 803 840 840 840 932 932 932 971 995 ...# indices of pt1 being triangualted
#..$ ind2   : int [1:6052] 522 580 440 581 803 440 650 781 440 512 ...
#..$ bp1    : logi [1:6052] TRUE FALSE FALSE FALSE FALSE FALSE ...indicate if pt1 is bdpt, on the rw
#..$ bp2    : logi [1:6052] FALSE FALSE FALSE FALSE FALSE FALSE ...
#..$ thirdv1: num [1:6052] -1 1208 1092 1404 1208 ... indices of the 3 vertex 
#..$ thirdv2: num [1:6052] 1275 1357 1404 1208 582 ...
    # if -1, lies outside, lies outside of the rectangular window rw. 


#$ summary :'data.frame':	2071 obs. of  9 variables:
#..$ x       : num [1:2071] -130 -114 -111 -111 -114 ...
#..$ y       : num [1:2071] 36.8 36.8 40.1 52.6 55.7 ...
#..$ n.tri   : num [1:2071] 2 2 2 2 2 3 2 2 6 6 ...the number of Delaunay triangles emanating from this point (x, y)
#..$ del.area: num [1:2071] 0.0767 0.0876 0.1005 0.0853 0.0684 ...
#..$ del.wts : num [1:2071] 0.000182 0.000208 0.000239 0.000203 
#..$ n.tside : num [1:2071] 3 3 3 3 3 4 3 3 6 6 ...
#..$ nbpt    : num [1:2071] 2 2 2 2 2 2 2 2 0 0 ...
#..$ dir.area: num [1:2071] 0.112 0.137 0.146 0.128 0.111 ...
#..$ dir.wts : num [1:2071] 0.000246 0.000302 0.00032 0.000282 

# (x, y) -> delaunany triangle -> n.tri -> 1/3 of total area of del tri
# note: dir.area: the area of Dirichelt tile surrounding the point
# 1/3 is because the different oreder of 3 pts create the same tri, so one set of 3 pts, 3 same tris

# ..$ dir.area: num [1:2071] 0.112 0.137 0.146 0.128 0.111 ...
# is the same as what we self calculated area.tess

#plot(Voronoi$dirsgs$ind1)

#which((Voronoi$dirsgs$ind1 = 1:nrow(mesh_locs)) == TRUE)




#------------------------#
# from Voronoi to Polygon
#------------------------#
Polys <- vector("list", nrow(mesh_locs))
for (i in 1:nrow(mesh_locs)) {
  X <- subset(Voronoi$dirsgs, (ind1 == i | ind2 == i) )
  X <- matrix(c(X$x1, X$x2, X$y1, X$y2, X$bp1, X$bp2), ncol = 3)
  X <- unique(X)
  
  if (sum(X[, 3]) > 0) {
    # becuase if the end pts are on the boudary, 
    # the triangulated pts will be on the bd as well, 
    # then the convex hull will unbounded,
    # then there will be no area
    # so need to add back the pt used to tri
    X <- rbind(X, c(mesh_locs[i, ], 0))
  }
  
  edges <- X[, 1:2]
  edges <- edges[chull(edges), ]
  
  Polys[[i]] <- as(edges, "gpc.poly")
}

str(Polys)



#------------------------#
# from Polygon to area
#------------------------#

area.tess <- rep(0, nrow(mesh_locs))
for (i in nrow(mesh_locs)) {
  area.tess[i] <- area.poly(Polys[[i]])
}

str(area.tess)








