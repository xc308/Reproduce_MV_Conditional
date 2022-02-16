##************************
## Reproduce 3.1 (Cressie)
##************************

##==========
## Packages
##==========
# for numeric operation
install.packages("Matrix")
library(Matrix) # for numeric operation

# for data manipulation
library(dplyr)  # for data manipulation
library(tidyr)  

# for plotting and figure arrangement into panels for publication
library(ggplot2)
library(gridExtra)
library(gird)
library(extrafont)
loadfonts()


##========
## Setup
##========
img_path <- "./Results/Fig"  # where to save the fig
show_figs <- 1               # show fig in documents
print_figs <- 0              # print fig to file


##--------------------------
## Set up simulation domain 
##--------------------------

# D = [-1, 1]
# spacing eta_i = 0.01, i = 1,..., 200
# collect the grid info into df, to which extra colums will be added.
# define: n1 the number of grid cells for Yvec_1; 
#         n2 the number of grid cells for Yvec_2;
# in this study, n1 = n2 = 200
# so total n = n1 + n2 = 400.


##----------------
## Construct grid
##----------------

ds <- 0.01  # space difference
s <- seq(-1 + ds/2, 1 - ds/2, by = ds) # result in a set of centroids
df <- data.frame(s = s, area = ds)

n1 <- n2 <- nrow(df)
n <- n1 + n2


##----------------------------------------
## True process and observation parameters
##----------------------------------------
sigma2_11 <- 1      ## Marginal variance of C11(.)
sigma2_21 <- 0.2    ## Marginal variance of C2|1(.)

kappa1 <- 25        ## land scale of C11(.)
kappa2 <- 75        ## land scale of C2|1(.)

nu11 <- nu2_1 <- 1.5  ## fixed smooth parameter for Matern

A <- 5              ## Amplitude of b(.)
r <- 0.3            ## Aperture of b(.)
delta <- -0.3       ## shift of b(.)

sigmav <- 0.5       ## obs error std



##===================================
## Matrix construction and simulation
##===================================

# after setting the required parameters, 
# now construct the full covariance matrix with block matrices
# C11, C12, C21, C22
# but first need to construct matrix B: 
  # a b(,) evaluated at grid cells multiplied by grid spacing (to approximate integration)
  # i.e.B(j, k) = eta_k * b(s_j, v_k)


#-----------------------------
# Construct b and B(include grid area for integration)
#-----------------------------
#h <- outer(df$s, df$s, FUN = "-")
#h[1:10, 1:10]

H <- t(outer(df$s, df$s, FUN = "-"))    # displacement
#H[1:10, 1:10] # each row is the dist btw si with the rest, i = 1..n


bisq_1d <- function(A = 1, d, delta = 0, r) {
  y <- abs(d - delta)
  A * (1 - (y/r)^2)^2 * (y < r)
}

B <- A * bisq_1d(d = H, delta = delta, r = r) * ds  # find B

#B_his <- A * bisquare_1d(h = H, delta = delta, r = r) * ds
#all(B == B_his) # [1] TRUE


#-----------------------------
# Construct required matrices
#-----------------------------

D <- abs(H)             # distance must be postive
Dvec <- as.double(c(D)) # vectorize distance

C11 <- Matern_32(Var = sigma2_11, Kappa = kappa1, d = Dvec)
C2_1 <- Matern_32(Var = sigma2_21, Kappa = kappa2, d = Dvec)

C12 <- C11 %*% t(B)  # t(C21)
C21 <- B %*% C11

C22 <- C2_1 + B %*% C11 %*% t(B)


##*************************
## Test if symmetry and pd
source("fun_test_sym_pd.R")

Test_sym_pd(C11)
Test_sym_pd(C2_1)
Test_sym_pd(C22)
##***********************

Sigma <- rbind(cbind(C11, C12), cbind(C21, C22))
Test_sym_pd(Sigma) # YES


##*****************
## To extract individual marginal and cross-covariance functions at 
## the midpoint location of D, i.e. s = 0

# to find the mid-point of D which is the 0 btw [-1,1]
D[n1/2, n1/2]  # [1] 0

all(diag(D) == 0)  # [1] TRUE

# the individual marginal cov function and cross-cov function
Cov11 <- Sigma[n1/2, 1:n1] # cov(Y1(s_mid), Y1(si)), i = 1,...n
Cov12 <- Sigma[n1/2, (n1 + 1):n]
Cov21 <- Sigma[n1 + n2/2, 1:n1]
Cov22 <- Sigma[n1 + n2/2, (n1 + 1):n]



















































