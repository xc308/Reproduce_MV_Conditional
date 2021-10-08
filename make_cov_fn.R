#=======================#
# Make cov-var function
#=======================#

# - makeC



#------#
# makeC
#------#

source("Matern_fn_2132.R")


makeC <- function(d, Var, Kappa, nu) {
  if (nu == 3/2) {
    C <- Matern_fn_32(d, phi = Var, Alpha = Kappa)
  } else if (nu == 1/2) {
    C <- Matern_fn_12(d, phi = Var, Alpha = Kappa)
  } else {
    C <- Matern_fn(d, phi = Var, Alpha = Kappa, nu)
  }
  C
}

























