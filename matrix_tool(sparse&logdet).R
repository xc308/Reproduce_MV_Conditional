###--------------------
### Matrix tools
###--------------------

#' @title Create a sparse identity matrix
#'
#' @description Creates a sparse identity matrix of size n x n
#' @param n size of matrix
#' @export
#' @examples
#' require(Matrix)
#' Q <- Imat(4)
#' 
Imat <- function(n) {
  require(Matrix)
  
  sparseMatrix(i = 1:n, j = 1:n, x = 1)
}



#' @title Find the log determinant
#'
#' @description Find the log determinant of a matrix Q from its Cholesky factor L.
#' @param L the Cholesky factor of Q
#' @export
#' @examples
#' require(Matrix)
#' Q <- sparseMatrix(i=c(1,1,2,2),j=c(1,2,1,2),x=c(0.1,0.2,0.2,1))
#' logdet(chol(Q))


log_det <- function(U) {
  # Find the log-det(Sigma) from its Cholesky U
  # U is the cholesky factor chol(Sigma)
  
  digaU <- diag(U) 
  return(2 * sum(log(diagU)))
}

#Sig <- sparseMatrix(i = c(1, 1, 2, 2), j = c(1, 2, 1, 2), x = c(1, 2, 2, 1))
#Sig
#log_det(chol(Sig))




