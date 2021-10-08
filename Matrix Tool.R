#=============
# Matrix Tool 
#============

# -Chol2inv
# -logdet

#--------------------#
# Example of chol2inv
#--------------------#
ma <- cbind(1, 1:3, c(1,3,7))
chma <- chol(ma)
t(chma) %*% chma  # ma
ma %*% chol2inv(chma)   # I


#---------#
# logdet
#---------#

logdet <- function(U) {
  # U is a chol factor
  # 1st extract the diag of U, get Uii
  diagU <- diag(U)

  return(2 * sum(log(diag(diagU))))
  
}











