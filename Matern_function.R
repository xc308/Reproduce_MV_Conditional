#=====================#
# Mat{\'e}rn functions
#=====================#



Matern32 <- function(r,var,kappa) { # r is d or h
  X <- covMat3_call(r,var,kappa)
  matrix(X,nrow=sqrt(length(r)))  # cov matrix is n*n
}

Matern21 <- function(r, var, kappa) {
  X <- covMat1_call(r, var, kappa)
  matrix(X, nrow = sqrt(length(r))) 
}



"Matern" <- function(d, scale = 1, phi = scale, smoothness = 0.5, nu = smoothness, 
         range = 1, alpha = 1/range) {
  # Matern's covariance function transcribed from Stein's book P31
  
  if (any(d) < 0) 
    stop("distance argument must be nonnegative")
  
  d <- d * alpha
  # avoid sending exact zeros to BasselK
  d[d == 0] <- 1e-10
  
  # the hairy constant
  con <- 2^(nu - 1) * gamma(nu)
  con <- 1/con
  
  # call to Bessel function from R base package
  return(phi * con * (d)^nu * besselK(d, nu))
  # return vector
}




# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"Matern" <- function(d, scale = 1, range = 1, alpha = 1/range,
                     smoothness = 0.5, nu = smoothness, phi = scale) {
  # Matern covariance function transcribed from Stein's book page 31                                 # nu==smoothness, alpha ==  1/range                                                                # GeoR parameters map to kappa==smoothness and phi == range                                        # check for negative distances                                                                
  if (any(d < 0))
    stop("distance argument must be nonnegative")
  d <- d * alpha
  # avoid sending exact zeroes to besselK
  d[d == 0] <- 1e-10
  # the hairy constant ...                                                                           # this is different from Stein to make this a correlation function when
  con <- (2^(nu - 1)) * gamma(nu)
  con <- 1/con
  # call to  Bessel function from R base package
  return(phi * con * (d^nu) * besselK(d, nu))
}





