#===================#
# bisquare function
#===================#

#str(h)  # num [1:4289041, 1:2]

bisq_2D <- function(h, delta = c(0, 0), r = 1, A = 1) {
  # note: r = 1 is inappropriate
  # r: the dist at which weight is set to 0
  y <- t(t(h) - delta)  # h: num [1:4289041, 1:2]
  y <- sqrt(y[, 1]^2 + y[, 2]^2)   # num [1:4289041] 
  A * (1 - (y / r)^2)^2 * (y < r)  # num [1:4289041] 
}




bisq_B <- function(h, delta = c(0, 0), r = 1, A = 1, 
                   area = rep(1, length(h[,1])), n1 = 10L, n2 = 10L) {
  BA <- area * bisq_2D(h = h, delta = delta, r = r, A = A) 
  matrix(BA, n2, n1, byrow = T)
}



# TRY out

#str(h) # num [1:4289041, 1:2]

#which(bisq_2D(hvec = h, r = 160) != 0)

#y <- h
#y <- sqrt(y[, 1]^2 + y[, 2]^2)
#y
#str(y)   # num [1:4289041] 
#range(y) # 147.7454 191.8233


#r <- 160
#length(which((1 - (y / r)^2)^2 * (y < r) != 0)) # [1] 609276


#bisq <- bisq_2D(h = h, r = 160)
#str(bisq)  # num [1:4289041]
#area <- rep(1, length(h[,1]))
#str(area)
#BA <- area*bisq
#BA2 <- bisq * area
#str(BA)
#str(BA2)

#b <- bisq_B(h = h, r = 160, n1 = n1, n2 = n2)
#str(b)  # num [1:2071, 1:2071]

#bisq_B <- bisq * Area
#str(bisq_B)  # num [1:4289041]

#B <- matrix(bisq_B, n2, n1, byrow = T)
#str(B)  #  num [1:2071, 1:2071] 0 0.00113 0.00201 0.00201 0.00116 ...