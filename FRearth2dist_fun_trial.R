#***********#
# Distances
#***********#

# -compare distances calculated using different functions


#============
# RFearth2dis
#============

D1 <- RFearth2dist(coord = as.matrix(mesh_locs))
str(D1)  # 'dist' num [1:2143485] 

D1_mat <- as.matrix(D1)
str(D1_mat) #num [1:2071, 1:2071] 



#===========================
# RFearth2cartesian + dist()
#===========================

mesh_cart <- RFearth2cartesian(coord = as.matrix(mesh_locs))
D2 <- dist(mesh_cart)
str(D2)  # 'dist' num [1:2143485] 
D2_mat <- as.matrix(D2)
str(D2_mat) # num [1:2071, 1:2071] 

d <- D1 - D2 
all(d == 0)  # [1] TRUE


#===================
# distances package 
#===================

install.packages("distances")
library(distances)

D3 <- distances(as.matrix(mesh_locs))
str(D3)
# 'distances' num [1:2, 1:2071] 0 15.453 18.98 24.436 24.455 15.453 0 4.592 16.062 18.844 18.98 4.592 0 12.483 15.904 24.436 16.062 12.483 0 4.40| __truncated__ ...
#- attr(*, "normalization")= num [1:2, 1:2] 1 0 0 1
#- attr(*, "weights")= num [1:2, 1:2] 1 0 0 1

D3_dist <- distance_matrix(D3)
str(D3_dist)  # 'dist' num [1:2143485] 
D3_matrix <- as.matrix(D3_dist)
str(D3_matrix)
# num [1:2071, 1:2071]
