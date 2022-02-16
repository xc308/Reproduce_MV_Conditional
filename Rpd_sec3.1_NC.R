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
str(Cov11)  # num [1:200]



#=============================
# Plot the covariance function
#=============================
Cov_df <- expand.grid(df$s, proc1 = c("Y1", "Y2"), 
                      proc2 = c("Y1", "Y2"))

Cov_df$cov <- c(Cov11, Cov21, Cov12, Cov22) # df 800 rows


g_cov <- LinePlotTheme() + 
  geom_line(data = Cov_df, aes(s, cov)) +
  facet_grid(proc1 ~ proc2)


if(print_figs) ggsave(g_cov,
                      filename = file.path(img_path, "cov_functions.png"),
                      width = 12, height = 10, family = "Arial")
# print_figs <- 0 so don't print, just save

# save first, then decide show or not
if(show_figs) print(g_cov, width = 12, height = 10)



#====================
# Generate noisy data
#====================

str(df) # data.frame':	200 obs. of  2 variables: s and area

# Given the full Joint Sigma, can simulate from the bivariate
# field jointly. 
# Obs are obtained by adding Gaussina error to the generated fields
# the simulations are added to the df


#---------------
# Generate data
#---------------
set.seed(16-02-22)

# since Y ~ MVN (mu, Sigma)
# Sigma^(-1/2)(Y - mu) = X ~ MVN(0, I) 
# so Y = Sigma^(1/2)X + mu
# only need to sample from X~ MVN(0, I)

sample_Y <- t(chol(Sigma)) %*% rnorm(n) # jointly sample Y1&Y2 1:400
# value of chol() The upper triangular factor of the Choleski 
# n = n1 + n2; n1 = n2 = 200

df <- df %>%
  mutate(sample_Y1 = sample_Y[1:n1],
         sample_Y2 = sample_Y[-(1:n1)],
         Z1 = sample_Y1 + sigmav * rnorm(n1),
         Z2 = sample_Y2 + sigmav * rnorm(n2))

Z <- matrix(c(df$Z1, df$Z2)) # concatenate/vectorize obs in to col matrix Z
# 400 by 1


#=======================================
# Demostrate the benefits of cokriging
#=======================================
# only keep half of Z1 on positive part
# inference on Y1 in the negative domain will be 
# facilitate via obs Z2

keep_z1 <- 101:200
keep_z2 <- 1:200
















































































