#==============================#
# Covariance Matrix Construction
#==============================#

# under conditional approach

source("Matern_function.R")


#' @title Construct full (bivariate) covariance/precision matrix
#' @name makeCov
#' @aliases makeSY
#' @aliases makeQY
#'
#' @description Construct the covariance or precision matrix for the bivariate model constructed using the conditional approach.
#' @param r vector of distances
#' @param var variance of C
#' @param var1 variance of C_{11}
#' @param var2 variance of C_{2|1}
#' @param kappa scale of C
#' @param kappa1 scale of C_{11}
#' @param kappa2 scale of C_{2|1}
#' @param B interaction matrix
#' @return Covariance (or precision) matrix
#' @details  Both C_{11} and C_{2|1} are Matern covariance functions with smoothness parameter equal to 3/2. Covariance matrices are computed from Matern covariance functions using the vector of distances \code{r}, so that Sigma[1,1] = cov(Y1(s),Y1(s+r[1])), Sigma[1 + n,1] = cov(Y2(s),Y1(s+r[1])) and so on. Currently the grids on which Y1 and Y2 are evaluated need to be identical.
#'
#' The matrix \eqn{B} is the interaction matrix. The full covariance matrix returned is \deqn{\Sigma =   \left( \begin{array}{cc}\Sigma_{11} & \Sigma_{11}B' \\ B \Sigma_{11} & \Sigma_{2\mid 1}+B\Sigma_{11}B'  \end{array}\right).}{Sigma = [Sigma_{11} & Sigma_{11}B' ; B Sigma_{11} & Sigma_{2|1} + B Sigma_{11}B'].}
#' @export
#' @examples
#' s <- 0 : 99
#' D <- as.matrix(dist(s)) [1:100, 1:100]
#' r <- as.vector(D) [1:10000]
#'


# Example 
B <- diag(100)
Sigma <- makeSY(r=r,var1=1,var2=1,kappa1=0.5,kappa2=0.1,B=B)
image(Sigma)


makeSY <- function(r, var1, var2, kappa1, kappa2, B, nu1 = 3/2, nu2 = 3/2) {
  S11  <- makeS(r, var = var1, kappa = kappa1)  # Sigma_11
  S2_1 <- makeS(r, var = var2, kappa = kappa2)  # Sigma_2|1
  BS_21 <- B %*% S11                            # BSigma_11
  
  # Sigma
  rbind(cbind(S11, t(BS_21)), cbind(BS_21), BS_21 %*% t(B) + S2_1)
  #rbind(cbind(S11, t(BS_21)), cbind(BS_21), tcrossprod(BS_21, B) + S2_1)
}



#' @rdname makeCov
#' @export
makeQY <- function(r, var1, var2, kappa1, kapp2, nu1 = 3/2, nu2 = 3/2) {
  Q11  <- makeQ(r, var = var1, kappa = kapp1, nu = nu1)   # Sigma_11^(-1)
  Q2_1 <- makeQ(r, var = var2, kappa = kappa2, nu = nu2) # Sigam_2|1^(-1)
  BQ2_1 <- crossprod(B, Q2_1)                             # B^T Sigam_2|1^(-1)
  
  # Sigma^(-1)
  rbind(cbind(crossprod(chol(Q2_1) %*% B) + Q11, -BQ2_1), cbind(-t(BQ2_1), Q2_1))
}


#' @rdname makeCov
#' @export
makeQ <- function(r, var, kappa, nu) {
  S <- makeS(r = r, var = var, kappa = kappa, nu = nu)
  chol2inv(chol(S))   # S^(-1) via cholesky decomposition
}


# Example of chol2inv
# 1st chol
# 2nd chol2inv
cma <- chol(ma  <- cbind(1, 1:3, c(1,3,7)))
t(cma) %*% cma #=ma
ma %*% chol2inv(cma)  # = I so chol2inv gives sigma^(-1) via cholesky decomp



#' @rdname makeCov
#' @export
makeS <- function(r, var, kappa, nu) {
  if (nu == 3/2) {
    S <- Matern32(r, var, kappa)
  } else if (nu == 1/2) {
    S <- Mater21(r, var, kappa)
  } else {
    svec <- Matern(d = r, phi = var, alpha = kappa, smoothness = nu)
    S <- matrix(svec, nrow = sqrt(length(r)))
  }
  S
}



#' @rdname makeCov
#' @export
makeS <- function(r,var,kappa,nu) {
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