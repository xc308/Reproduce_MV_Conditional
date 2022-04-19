#**********************
# Algo_Stepwise_Sigma
#**********************
# 1-d simulation
# b(.,.) Matern
# b(.,.) remains the same for all pairs of variables



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
# Construct b and B (include grid area for integration)
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

p <- 5
n <- n_l
Sgm_sp <- Matrix(0, p*n, p*n, sparse = T) #str(Sgm_sp) # Formal class 'ddiMatrix'
Sgm_sp[1:(1*n), 1:(1*n)] <- C11
Phi_l <- C2_1
B_lk <- t(B_kl_sp)

for(l in seq(2, p)) {
  for(k in seq(1, (p-1))) {
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
    








