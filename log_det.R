#=========#
# log_det
#=========#

log_det <- function(U) {
  # U is a chol factor
  Uii <- diag(U)
  return(2 * sum(log(Uii)))
}

