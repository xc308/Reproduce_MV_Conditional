#========================================
# Algorithm for Stepwise propagated Sigma 
#========================================

a <- 1
b <- 5
c <- 3
d <- 6

M <- rbind(cbind(a, b), cbind(c, d))

M <- cbind(M, 0)
M <- rbind(M, 0)
M[2, 3] <- a
M[3, 2] <- b
M[3, 3] <- c
M


Sgm <- matrix(NA, nrow = 5, ncol = 5)
Sgm[1, 1] <- 2
for(l in seq(2, 5)) {
  for(k in seq(1, 4)) {
    if (abs(k - l) >= 2) {
      Sgm_kl = Sgm_lk = 0
      
      Sgm[k, l] <- Sgm_kl
      Sgm[l, k] <- Sgm_kl
    } else {
      Sgm_kl = Sgm_lk = 1
      Sgm_ll = 2
      
      Sgm[k, l] <- Sgm_kl
      Sgm[l, k] <- Sgm_lk
      Sgm[l, l] <- Sgm_ll
    }
  }
}
Sgm


#--------------------------------------#
# with Cross-cov and var update formula
#--------------------------------------#

Sgm <- matrix(NA, 5, 5)
Sgm[1, 1] <- 2
b_lk <- 1
phi_l <- 1
for(l in seq(2, 5)) {
  for(k in seq(1, 4)) {
    if (abs(k - l) >= 2) { # k !in pa(l)
      Sgm_kl = Sgm_lk = 0
      Sgm[k, l] = Sgm_kl
      Sgm[l, k] = Sgm_lk
      
    } else{
      i = k - (l - 1)  # see notes below
      Sgm_kk <- Sgm[k-i, k-i]
      Sgm_kl = Sgm_lk = Sgm_kk * b_lk
      Sgm_ll = phi_l + b_lk * Sgm_kl
      
      Sgm[k, l] <- Sgm_kl
      Sgm[l, k] <- Sgm_lk
      Sgm[l, l] <- Sgm_ll
    }
  }
}
Sgm

# Notes
# set a count i = k - (l - 1)
# is to ensure within each l, no matter how k changes, 
# the corresponding Sigma_kk will remain the same, 
# and won't progagate as k changes
# but with different l, Sigma_kk updates 1 lag forward

# so when l = 2, k-i === 1
# when l = 3, k-i === 2, so on


## try sparse structure
sparse_Sgm <- as(Sgm, "sparseMatrix")
sparse_Sgm

sparse_Sgm <- Matrix(Sgm, sparse = T)

print(object.size(Sgm))
# 416 bytes

print(object.size(sparse_Sgm))
summary(sparse_Sgm)


#--------------------
# Block matrix Sigma
#--------------------

n <- 2
Sgm <- matrix(NA, 5*n, 5*n)
Sgm[1:(1*n), 1:(1*n)] <- matrix(2, nrow = n, ncol = n)
#B_lk <- I_mat(n)
#Phi_l <- I_mat(n)
for(l in seq(2, 5)) {
  for(k in seq(1, 4)) {
    if (abs(k - l) >= 2) { # k !in pa(l)
      #Sgm_kl = Sgm_lk = O_mat(n)  
      Sgm_kl = Sgm_lk = matrix(0, 2, 2)
      Sgm[((k-1)*n + 1):(k*n), ((l-1)*n + 1):(l*n)] = Sgm_kl
      Sgm[((l-1)*n + 1):(l*n), ((k-1)*n +1):(k*n)] = Sgm_lk
      
    } else{
      i = k - (l - 1)  # see notes below
      Sgm_kk <- Sgm[((k-i)*n -1):((k-i)*n), ((k-i)*n - 1):((k-i)*n)]
      
      B_lk <- cbind(c(1, 0), c(0, 1))
      Phi_l <- cbind(c(1, 0), c(0, 1))
      Sgm_kl = Sgm_lk = Sgm_kk %*% B_lk
      Sgm_ll = Phi_l + B_lk %*% Sgm_kl
      
      Sgm[(n*k - 1):(n*k), (n*l - 1):(n*l)] <- Sgm_kl
      Sgm[(n*l - 1):(n*l), (n*k - 1):(n*k)] <- Sgm_lk
      Sgm[(n*l - 1):(n*l), (n*l - 1):(n*l)] <- Sgm_ll
      
    }
  }
}
Sgm


#-----------
# Questions
#-----------
# Why got error like length missmatch for matrix replacement
# when adopt sparseMatrix?


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sparse matrix element replace
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(Matrix)

# Example 1
n <- 3
K <- Matrix(0, n*n, n*(n-1)/2, sparse = T )
str(K)

for (i in 1:(n - 1)) {
  K[((i - 1) * (n + 1) + 2) : (i * n), 
    (1 + (i - 1) * (n - i/2)) : (i * (n - i) * (i + 1))/2] <- diag(n - i)
  
}
K









#----------------------------
# Sparse Unitilites functions
#----------------------------

I_mat <- function(n) {
  sparseMatrix(i = 1:n, j = 1:n, x = 1)
}
I_mat(2)  
 
B_lk <- cbind(c(1, 0), c(0, 1))
Try <- I_mat(2) %*% B_lk
str(Try)



O_mat <- function(n) {
  sparseMatrix(i = 1:n, j = 1:n, x = 0)
}
O_mat(2)

