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
      Sgm_kl = Sgm_lk = matrix(0, n, n)
      Sgm[((k-1)*n + 1):(k*n), ((l-1)*n + 1):(l*n)] = Sgm_kl
      Sgm[((l-1)*n + 1):(l*n), ((k-1)*n +1):(k*n)] = Sgm_lk
      
    } else{
      i = k - (l - 1)  # see notes below
      Sgm_kk <- Sgm[((k-i)*n - (n-1)):((k-i)*n), ((k-i)*n - (n-1)):((k-i)*n)]
      
      B_lk <- diag(n)  # n * n identity
      Phi_l <- diag(n)
      Sgm_kl = Sgm_lk = Sgm_kk %*% B_lk
      Sgm_ll = Phi_l + B_lk %*% Sgm_kl
      
      Sgm[(n*k - (n-1)):(n*k), (n*l - (n-1)):(n*l)] <- Sgm_kl
      Sgm[(n*l - (n-1)):(n*l), (n*k - (n-1)):(n*k)] <- Sgm_lk
      Sgm[(n*l - (n-1)):(n*l), (n*l - (n-1)):(n*l)] <- Sgm_ll
      
    }
  }
}
Sgm


#-----------
# Questions
#-----------
# Why got error like length missmatch for matrix replacement
# when adopt sparseMatrix?
  # as they are different in structure  
  # sparse matrix can join normal calculation, e.g. addition and multiplication
  # but is cannot be assigned to a regular matrix unless they align their format first
  # either both in sparse or both in regular



#--------------------------------------#
# with Cross-cov and var update formula 
# Sparse matrix element replace
#--------------------------------------#

Sgm_sp <- Matrix(0, 5, 5, sparse = T)
Sgm_sp[1, 1] <- 2
b_lk <- 1
phi_l <- 1
for(l in seq(2, 5)) {
  for(k in seq(1, 4)) {
    if (abs(k - l) >= 2) { # k !in pa(l)
      Sgm_kl = Sgm_lk = 0
      Sgm_sp[k, l] = Sgm_kl
      Sgm_sp[l, k] = Sgm_lk
      
    } else{
      i = k - (l - 1)  # see notes below
      Sgm_kk <- Sgm_sp[k-i, k-i]
      Sgm_kl = Sgm_lk = Sgm_kk * b_lk
      Sgm_ll = phi_l + b_lk * Sgm_kl
      
      Sgm_sp[k, l] <- Sgm_kl
      Sgm_sp[l, k] <- Sgm_lk
      Sgm_sp[l, l] <- Sgm_ll
    }
  }
}
Sgm_sp



#--------------------
# Block matrix Sigma
# sparse matrix
#--------------------

n <- 2
n <- 3
n <- 5
#n <- 200
Sgm_sp <- Matrix(0, 5*n, 5*n, sparse = T) #str(Sgm_sp) # Formal class 'ddiMatrix'
Sgm_sp[1:(1*n), 1:(1*n)] <- Matrix(2, nrow = n, ncol = n, sparse = T)
B_lk <- I_mat(n)   # "dgCMatrix" can be product with "ddiMatrix"
Phi_l <- I_mat(n)  # ok
#B_lk <- Matrix(diag(n), sparse = T)  # ok
#Phi_l <- Matrix(diag(n), sparse = T) # ok
for(l in seq(2, 5)) {
  for(k in seq(1, 4)) {
    if (abs(k - l) >= 2) { # k !in pa(l)
      Sgm_kl = Sgm_lk = Matrix(0, n, n, sparse = T)  
      #Sgm_kl = Sgm_lk = matrix(0, n, n)
      Sgm_sp[((k-1)*n + 1):(k*n), ((l-1)*n + 1):(l*n)] = Sgm_kl 
      Sgm_sp[((l-1)*n + 1):(l*n), ((k-1)*n +1):(k*n)] = Sgm_lk  
      
    } else{
      i = k - (l - 1)  # see notes below
      Sgm_kk <- Sgm_sp[((k-i)*n - (n-1)):((k-i)*n), ((k-i)*n - (n-1)):((k-i)*n)]
      
      #B_lk <- diag(n)  # n * n identity
      #Phi_l <- diag(n)
      Sgm_kl = Sgm_kk %*% t(B_lk)
      Sgm_lk = B_lk %*% Sgm_kk
      Sgm_ll = Phi_l + B_lk %*% Sgm_kl
      
      Sgm_sp[(n*k - (n-1)):(n*k), (n*l - (n-1)):(n*l)] <- Sgm_kl
      Sgm_sp[(n*l - (n-1)):(n*l), (n*k - (n-1)):(n*k)] <- Sgm_lk
      Sgm_sp[(n*l - (n-1)):(n*l), (n*l - (n-1)):(n*l)] <- Sgm_ll
      
    }
  }
}

Sgm_sp
print(object.size(Sgm_sp), units = "Kb")
# 2.9 Kb

image(Sgm_sp)
rm(Sgm_sp)


#----------------------------
# Sparse Unitilites functions
#----------------------------

I_mat <- function(n) {
  sparseMatrix(i = 1:n, j = 1:n, x = 1)
}
I_mat(2)  
 

B_lk <- I_mat(2)
Try <- I_mat(2) %*% B_lk
str(Try)

M <- Matrix(1, 10, 10, sparse = T)
str(M)

M[1:3, 1:3]

n = 2
Matrix(diag(n), sparse = T)



O_mat <- function(n) {
  sparseMatrix(i = 1:n, j = 1:n, x = 0)
}
O_mat(2)

