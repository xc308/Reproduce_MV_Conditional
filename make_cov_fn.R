#=======================#
# Make cov-var function
#=======================#

# - makeS
# -


#------#
# makeS
#------#

source("Matern_fn_2132.R")

makeS <- function(d, Var = phi, Kappa = Alpha, nu) {
  if (nu == 1/2) {
    S <- Matern_fn_12 (d, phi, Alpha)
  } else (nu == 3/2) {
    S <- Matern_fn_32(d, phi, Alpha )
  } else {
    S <- Matern_fn(d, phi, Alpha, nu)
  }
  S
}







