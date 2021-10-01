#================#
# Matern Function
#================#

# - Matern
# - Matern 1/2
# - Matern 3/2

#----------#
# Matern
#----------#

Matern_fn <- function(phi = 1, Alpha = 1, nu = 0.05, d) {
  if (any(d) < 0) 
    stop("distance d must be non-neg.")
  
  d <- Alpha * d
  d[d == 0] <- 1e-6   # avoid sending exact 0 to basselK
  
  const <- 2^(nu - 1) * gamma(mu)
  const <- 1/const
  
  X <- phi * const * (d)^nu * besselK(d, nu)  # a vector
  matrix(X, nrow = sqrt(length(X)))
  
}



#------------#
# Matern 1/2
#------------#

Matern_fn_12 <- function(phi = 1, Alpha = 1, d) {
  X <- phi * exp(- Alpha * d)
  matrix(X, nrow = sqrt(length(d)))  # 
}



#------------#
# Matern 2/3
#------------#

Matern_fn_32 <- function(phi = 1, Alpha = 1, d) {
  X <- phi * (1 + Alpha * d) * exp(- Alpha * d)
  matrix(X, nrow = sqrt(length(X)))
}































