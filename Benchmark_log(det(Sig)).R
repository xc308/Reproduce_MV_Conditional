#======================================
# bench mark different methods for log(det(Sigma))
#======================================
install.packages("microbenchmark")
library(microbenchmark)


# Method1: det(Sigma) in base, log(det(Sigma))
# Method2: determinant(Sigma, log = T)$modulus in base
# Method3: 
        # Sigma is symmetric, must have real eigen values, and must be diagonalizable
          # to a diag matrix with eigen values on the main diag
          # so, det(Sigma) = prod(eigen(Sigma)$val)
          # log(prod(eigen(Sigma)$val)) = sum(log(eigen(Sigma)$val))

# Method4:
        # |Sigma| = |L||U|, L, U are chol lower and upper factor
                # = |t(U)||U|
                # = (prod(Uii))^2
        # log(|Sigma|) = 2 * sum(log(diag(U)))


log_eigen <- function(Sigma) {
  
  sum(log(eigen(Sigma)$val))
}

log_chol <- function(Sigma) {
  
  chol_L <- chol(Sigma)
  2 * sum(log(diag(chol_L)))
}



microbenchmark(
  log(det(Sigma)),
  determinant(Sigma, logarithm = T)$modulus, 
  log_eigen(Sigma),
  log_chol(Sigma)
)

#Unit: milliseconds
#       expr       min        lq      mean
#log(det(Sigma))  8.580209  8.609167  8.784936
#determinant(..)  8.575709  8.595751  8.652189
#log_eigen(Sigma) 70.037042 70.239001 72.802494
#log_chol(Sigma)  6.324376  6.349105  6.452497

# median        uq        max neval
#8.639314  8.767314  13.347292   100
#8.617021  8.642397   9.070751   100
#70.565230 72.934000 195.356543   100
#6.378251  6.500438   8.423334   100


#-----------
# Conclusion
#-----------
# the fastest way to calculate log(det(Sigma)) is 
  # 2 * sum(log(diag(chol(Sigma))))

# the 2nd fastest is use determinant(Sigma, log = T)$modulus

# the slowest: use sum(log(eigen(Sigma)$val))





