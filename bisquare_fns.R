#===================#
# Bisquare function
#===================#

install.packages("roxygen2")
library(roxygen2)
devtools::load_all()

source("Covariance_Matrix_Construction.R")

#' @title Bisquare function
#' @name bisquare_fns
#' @aliases bisquare_1
#' @aliases bisquare_1d
#' @aliases bisquare_2d
#' @aliases bisquare_B
#'
#' @description The bisquare function
#' @param h displacement (1d)
#' @param h1 first component of displacement vector (2d)
#' @param h2 second component of displacement vector (2d)
#' @param delta shift parameter(s)
#' @param r aperture parameter
#' @param A gain parameter (amplitude)
#' @param areas area associated with each column in B
#' @param n1 number of rows
#' @param n2 number of columns
#' 
#' @details The bisquare function (shifted by \eqn{\Delta}) is given by \deqn{b(s,v) \equiv \left\{\begin{array}{ll} A\{1 - (|v- s - \Delta|/r)^2\}^2, &| v -s  - \Delta| \le r \\ 0, & \textrm{otherwise}. \end{array} \right.}{b(s,v) =  A{1 - (|v- s - d|/r)^2}^2  if | v -s  - d| <= r,  0 otherwise}
#' The function \code{bisquare_1d} accepts \code{h} in any shape or size, while for \code{bisquare_2d}, \code{h1} and \code{h2} need to be vectors of the same size. The parameter \code{delta} needs to be equal to one or two in length, depending on the spatial dimension.
#' 
#' The function \code{bisquare_B} is used to construct the matrix \eqn{B} given the bisquare parameters. It is meant to be used in problems of 2 dimensions.
#' @export
#' @examples
#' h <- seq(-10,10,0.1)
#' y <- bisquare_1d(h=h,delta=3,r=4)
#' y <- bisquare_1d(h=h,delta=3,r=4)
#' plot(h,y)
#' 
#' 
bisquare_1d <- function(h, delta = 0, r = 1, A = 1) {
  y <- abs(h - delta)
  A * (1 - (y / r)^2)^2 * (y < r)
}


bisquare_2d <- function(h1, h2, delta = c(0, 0), r = 1, A = 1) {
  bisquare_call(h1, h2, delta, r, A)
}


## NOT via C using R (slow)
bisquare_notC <- function(h, delta, r = 1, A = 1) {
  y <- t(t(h) - delta)
  y <- sqrt(y[1]^2 + y[2]^2)
  A * (1 - (y/r)^2)^2 *(y < r)
}



#' @export
#' @rdname bisquare_fns
bisquare_B <- function(h1, h2, delta = c(0, 0), r = 1, A = 1,
                       areas = 1, n1 = 10L, n2 = 10L) {
  this_z <- bisquare_2d(h1, h2, delta, r, A) * areas
  matrix(this_z, n2, n1, byrow = T) # faster than transpose
  # 1st decide the dimension: n2 * n1
  # 2nd decide the content: byrow = T
}





