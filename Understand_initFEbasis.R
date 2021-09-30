##==============#
## initFEbasis()
##==============#
install.packages("deldir")
library(deldir)

install.packages("gpclib")
library(gpclib)

# fns ref:https://github.com/andrewzm/bicon/blob/master/R/Finite_elements.R
# voronoi tesselation ref: https://cran.r-project.org/web/packages/deldir/deldir.pdf
# general polygon clip gpc ref: https://cran.r-project.org/web/packages/gpclib/gpclib.pdf



#---------------------------------------#
# play: construct voronoi tess from mesh
#---------------------------------------#

Voro <- deldir(x = mesh_locs[, 1],
       y = mesh_locs[, 2], 
       plotit = T, 
       sort = F, 
       rw = c(min(mesh_locs[, 1]) - 0.00001,
              max(mesh_locs[, 1]) + 0.00001,
              min(mesh_locs[, 2]) - 0.00001, 
              max(mesh_locs[, 2]) + 0.00001))

str(Voro)

# $ dirsgs  :'data.frame':	6052 obs. of  10 variables:
#..$ x1     : num [1:6052] -125 -128 -126 -130 -124 ... coords of end pts of tile
#..$ y1     : num [1:6052] 49.9 49.7 48.8 50.9 49.5 ...
#..$ x2     : num [1:6052] -125 -128 -126 -130 -123 ...
#..$ y2     : num [1:6052] 49.9 49.7 48.8 50.9 49.3 ...
#..$ ind1   : int [1:6052] 15 17 21 23 27 27 30 32 33 34 ...indices of pts being triagulated
#..$ ind2   : int [1:6052] 14 16 20 22 18 26 18 19 30 26 ...
#..$ bp1    : logi [1:6052] FALSE FALSE FALSE FALSE FALSE FALSE ...
#..$ bp2    : logi [1:6052] FALSE FALSE FALSE FALSE FALSE FALSE

#The entry in column bp1 indicates whether the first endpoint of
#the corresponding edge of a tile is a boundary point 
#(a point on the boundary of the rectangular window). 


sum(Voro$dirsgs$bp1 + Voro$dirsgs$bp2)  # [1] 158

# $ summary :'data.frame':	2071 obs. of  9 variables:
#..$ x       : num [1:2071] -130 -114 -111 -111 -114 ...
#..$ y       : num [1:2071] 36.8 36.8 40.1 52.6 55.7 ...
#..$ n.tri   : num [1:2071] 2 2 2 2 2 3 2 2 6 6 ... number of triangles from this pt
#..$ del.area: num [1:2071] 0.0767 0.0876 0.1005 0.0853 0.068  
    # (1/3 of the total area of all the Delaunay triangles emanating from the points


plot(Voro, wlines = "tess")
plot(Voro, wlines = "triang")



?chull 
# Computes the subset of points which lie on the convex hull of the set of points specified.
# Value:
#An integer vector giving the indices of the unique points lying on the convex hull, in clockwise order



#---------------------------------------#
# play: construct polygon from voronoi
#---------------------------------------#

for (i in 1:nrow(mesh_locs)) {
  X <- subset(Voro$dirsgs,(ind1 == i | ind2 == i))
}

head(X)
tail(X)
dim(X)  # 5 10
str(X)

X <- matrix(c(X$x1,X$x2,X$y1,X$y2,X$bp1,X$bp2),ncol=3)
head(X) # vecx, vecy, vecbp

X <- unique(X)
if(sum(X[,3])>0) {
  X <- rbind(X,c(p[i,],0))
}


edges <- X[,(1:2)]
head(edges)
#             x1       y1
# 2042 -129.3418 37.01165
# 2169 -129.6275 37.12063
# 2739 -129.9218 36.82638

chull(edges) # 3 2 1
# An integer vector giving the indices of the unique points 
# lying on the convex hull, in clockwise order.

  
edges <- edges[chull(edges), ]
as(edges,"gpc.poly")
#GPC Polygon
#Num. Contours:  1 
#Num. Vertices:  3 
#BBox (X):  -129.9218 --> -129.3418 
#BBox (Y):  36.82638 --> 37.12063 


#-----------------
# from Poly to area
#-----------------
area.poly(as(edges,"gpc.poly"))  # [1] 0.05807568



#======================================

n <- dim(mesh_locs)[1]
Poly <- vector("list", n)
for (i in 1: n) {
  X <- subset(Voro$dirsgs, (ind1 == i | ind2 ==i))
  X <- matrix(c(X$x1, X$x2, X$y1, X$y2, X$bp1, X$bp2), ncol = 3)
  X <- unique(X)
  
  if (sum(X[, 3]) > 0) {
    # find endpts is bp pt, set ind of bp to 0
   X <- rbind(X, c(mesh_locs[i, ], 0))
  }
  
  edges <- X[, 1:2]
  edges <- edges[chull(edges), ]  # the smallest convex set that contains it. 
  Poly[[i]] <- as(edges, "gpc.poly") # gpc.poly A class for representing polygons composed of multiple contours, some of which may be holes.

}

str(Poly)
# List of 2071
#$ :Formal class 'gpc.poly' [package "gpclib"] with 1 slot
#.. ..@ pts:List of 1
#.. .. ..$ :List of 3
#.. .. .. ..$ x   : num [1:4] -129 -130 -130 -129
#.. .. .. ..$ y   : num [1:4] 36.8 36.8 37.1 37
#.. .. .. ..$ hole: logi FALSE



#---------------------------------------#
# play: from Polygon to area
#---------------------------------------#
area.tess <- rep(0, nrow(mesh_locs))
for (i in 1:nrow(mesh_locs)) {
  area.tess[i] <- area.poly(Poly[[i]])  # : Compute and return the sum of the areas of all contours in a "gpc.poly" object.
}

str(area.tess)  # num [1:2071] 0.112 0.137 0.146 0.128 0.111 ...

head(Voro$summary$dir.area, 10)
head(area.tess, 10)
head(round(area.tess, 6), 10) 

A.tess <- round(area.tess, 6)

A.tess == Voro$summary$dir.area # all true



#----------------------
# play: collect into df
#----------------------

FEbasis_df <- data.frame(x = mesh_locs[, 1], 
           y = mesh_locs[, 2], 
           n = 1:nrow(mesh_locs),
           area.tess = area.tess)

head(FEbasis_df)






initFEbasis = function(p,t,K) {
  fn <- pars <- list()
  pars$p <- p
  pars$t <- t
  pars$K <- K
  df <- data.frame(x = pars$p[,1],
                   y = pars$p[,2],
                   n = 1:nrow(p))
  pars$vars <- df
  # Do tessellation
  Voronoi <- deldir(pars$p[,1],
                    pars$p[,2],
                    plotit='F',
                    sort=F,
                    rw=c(min(pars$p[,1])-0.00001,
                         max(pars$p[,1])+0.00001,
                         min(pars$p[,2])-0.00001,
                         max(pars$p[,2])+.00001))
  pars$pol <- PolygonfromVoronoi(Voronoi,pars$p)
  
  pars$vars$area_tess = rep(0,nrow(p))
  for (i in 1:nrow(p)) {
    pars$vars$area_tess[i] <- area.poly(pars$pol[[i]])
  }
  this_basis <- new("FEBasis", pars=pars, n=nrow(p), fn=fn)
  return(this_basis)
}



PolygonfromVoronoi <- function(Voronoi,p) {
  n = dim(p)[1]
  polygons <- vector("list",n)
  for (i in 1:n)  {
    X <- subset(Voronoi$dirsgs,(ind1 ==i | ind2 == i))
    X <- matrix(c(X$x1,X$x2,X$y1,X$y2,X$bp1,X$bp2),ncol=3)
    X <- unique(X)
    if(sum(X[,3])>0) {
      X <- rbind(X,c(p[i,],0))
    }
    edges <- X[,(1:2)]
    
    edges <- edges[chull(edges), ]
    polygons[[i]] <- as(edges,"gpc.poly")
  }
  return(polygons)
  
}


#dirsgs A data frame with 10 columns. 
  #The first 4 entries of each row are the coordinates of the endpoints of 
  #one the edges of a Dirichlet tile, in the order (x1,y1,x2,y2).
  #The 5th 6th entries, in the columns named ind1 and ind2, are the indices
  #of the two points, in the set being triangulated, which are separated by that edge. 

