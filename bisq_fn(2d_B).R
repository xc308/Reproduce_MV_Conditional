#===================#
# bisquare function
#===================#

str(h)  # num [1:4289041, 1:2]

bisq_2D <- function(h, delta = c(0, 0), r = 1, A = 1) {
  # note: r = 1 is inappropriate
  # r: the dist at which weight is set to 0
  y <- t(t(h) - delta)
  y <- sqrt(y[, 1]^2 + y[, 2]^2)   # num [1:4289041] 
  A * (1 - (y / r)^2)^2 * (y < r)  # num [1:4289041] 
}


bisq_B <- function(h, delta = c(0, 0), r = 1, A = 1, area = 1, n1 = 10L, n2 = 10L) {
  BA <- bisq_2D(h, delta, r, A) * area
  matrix(BA, n2, n1, byrow = T)
}



