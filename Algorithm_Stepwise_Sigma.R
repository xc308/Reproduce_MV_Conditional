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
n <- 3
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
# 1. Why got error like length missmatch for matrix replacement?
  # Because due to different format of sparse matrix
# 2. when adopt sparseMatrix?
  # as they are different in structure  
  # sparse matrix can join normal calculation, e.g. addition and multiplication
  # but is cannot be assigned to a regular matrix unless they align their format first
  # either both in sparse or both in regular

# 3. When will the additional storage of spare matrix (e.g. index) is compensated?
  # only when the # of non-zero < (m (n − 1) − 1) / 2 if M is m by n


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



#-----------------------------------------
# Block matrix Sigma
# sparse matrix
# location s itself in 1-d format
# multi-dimension for a variable pair (k,l)
#-----------------------------------------


##~~~~~~~~~~~~~~~~~~~~~~~~
# set up simulation domain
##~~~~~~~~~~~~~~~~~~~~~~~~

# D <- [-1, 1]
# spacing eta_i = 0.1, i = 1, ... 20
# n_l: # of grid cells for process vec(Y_k)
# n_k: # of grid cells for process vec(Y_l)
# k = 1, ..., 4; l = 2, ..., 5

# in this study, n_l = n_k = 20
# n = Sum_{i = 1}^{i = 5} (n_i) = 100


##~~~~~~~~~~~~~~~~~~~~~~~~
# construct process grid
##~~~~~~~~~~~~~~~~~~~~~~~~

ds = 0.1
s <- seq(-1 + ds/2, 1 - ds/2, by = ds)
str(s) # num [1:20] -0.95 -0.85 

df <- data.frame(s = s, area = ds)
n_l <- n_k <- nrow(data.frame(df))


##~~~~~~~~~~~~~~~~~~~~~~~~
# True process parameters
##~~~~~~~~~~~~~~~~~~~~~~~~

sigma2_11 <- 1       # marginal variance for initial C11 (or Sgm_11)
sigma2_21 <- 0.2     # marginal variance for the conditional variance C2|1 (Sgm_2|1)

kappa_11 <- 25      # large scale parameter for C11
kappa_21 <- 75      # large scale parameter for C2|1

nu_11 <- nu_21 <- 1.5 # special Matern form


A <- 1              # amplitude for b(.;,)
delta <- -0.3       # fix shift parameter for all variable pair
r <- 0.3            # aperture parameter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct b and B(include grid area for integration)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

h <- outer(df$s, df$s, FUN = "-")
str(h) # num [1:20, 1:20] 

H <- t(outer(df$s, df$s, FUN = "-"))

b_kl <- function(A = 1, delta, r, d) {
  y <- abs(d - delta)  # deliberately exaggerate asymmetry
  A * (1- (y^2))^2 * (y < r)
}

B_kl <- A * b_kl(delta = delta, r = r, d = H) * ds
str(B_kl) # num [1:20, 1:20]

B_kl_sp <- Matrix(B_kl, sparse = T)
# 20 x 20 sparse Matrix of class "dtCMatrix"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct required known matrices
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C11 (Sgm_11)
# C2|1(Sgm_21)

D <- abs(H)     # distance must be positive
D_vec <- as.double(c(D))   # vectorize distance

source('Matern_32.R')
C11 <- Matern_32(Var = sigma2_11, Kappa = kappa_11, d_vec = D_vec)
C2_1 <- Matern_32(Var = sigma2_21, Kappa = kappa_21, d_vec = D_vec)



##~~~~~~~~~~~~~~~~~~~~~~~
## Test symmetry and pd
##~~~~~~~~~~~~~~~~~~~~~~~

source("fun_test_sym_pd.R")

Test_sym_pd(C11)
Test_sym_pd(C2_1)
# [1] "Symmetric: Yes"
# [1] "p.d.: Yes"


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Construct full covariance matrix
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

n <- n_l
Sgm_sp <- Matrix(0, 5*n, 5*n, sparse = T) #str(Sgm_sp) # Formal class 'ddiMatrix'
Sgm_sp[1:(1*n), 1:(1*n)] <- C11
Phi_l <- C2_1
B_lk <- t(B_kl_sp)

for(l in seq(2, 5)) {
  for(k in seq(1, 4)) {
    if (abs(k - l) >= 2) { # k !in pa(l)
      Sgm_kl = Sgm_lk = Matrix(0, n, n, sparse = T)  
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
# 60.9 Kb

image(Sgm_sp)
## 1. obvious variance matrix Sgm_kk and cross cov matrix
  ## Sgm_kl propagation, and the later the order of variables
  ## the more information it borrowed from its predecessors

## Q: does this mean the last variable have the smallest bias?
  ## need to calibrate with real data
  ## esp. for SS, OM

## 2. within each cross-cov matrix, the asymmetry 
  ## is displayed as either lower or upper diagonal
## 3. although Sgm_lk = t(Sgm_kl) as a whole in order 
  ## to satisfy the condition for a valid multivaraite
  ## joint covariance matrix.



#----------------------------
# Sparse Unitilites functions (depretched)
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


n = 3
O_mat <- function(n) {
  sparseMatrix(i = 1:n, j = 1:n, x = 0)
}
O_mat(3)


Matrix(0, n, n, sparse = T)  
identical(O_mat(3), Matrix(0, n, n, sparse = T))

