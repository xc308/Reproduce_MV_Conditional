#=======================#
# Make cov-var function
#=======================#

# - makeC



#------#
# makeC
#------#

source("Matern_fn_2132.R")


makeC <- function(Var, Kappa, nu, D) {
  if (nu == 3/2) {
    C <- Matern_fn_32(phi = Var, Alpha = Kappa, d = D)
  } else if (nu == 1/2) {
    C <- Matern_fn_12(phi = Var, Alpha = Kappa, d = D)
  } else {
    C <- Matern_fn(phi = Var, Alpha = Kappa, nu = nu, d = D)
  }
  C
}

























