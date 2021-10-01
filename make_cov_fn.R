#=======================#
# Make cov-var function
#=======================#

# - makeC
# -


#------#
# makeC
#------#

source("Matern_fn_2132.R")


makeC <- function(d, Var, Kappa, nu) {
  if (nu == 3/2) {
    C <- Matern_fn_32(phi = Var, Alpha = Kappa, nu)
  } else if (nu == 1/2){
    C <- Matern_fn_12 (d, phi = Var, Alpha = Kappa)
  } else {
    C <- Matern_fn(d, phi = Var, Alpha = Kappa, nu)
  }
  C
}








makeS <- function(r,var,kappa,nu) { # r is distance d, differ from aperture parameter
  if(nu == 3/2)  {
    S <- Matern32(r,var,kappa)
  } else if (nu == 1/2) {
    S <- Matern12(r,var,kappa)
  } else {
    svec <- Matern(d=r,scale=var,alpha=kappa,smoothness = nu)
    S <- matrix(svec, nrow = sqrt(length(r)))
  }
  S
}











