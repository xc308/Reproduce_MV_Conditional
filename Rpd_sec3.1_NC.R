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
library(grid)
library(extrafont)
loadfonts()


##========
## Setup
##========
img_path <- "./Results/Fig"  # where to save the fig
show_figs <- 1               # show fig in documents
#print_figs <- 0              # print fig to file
print_figs <- 1               # save to the Results folder


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


##-----------------------
## Construct process grid
##-----------------------

ds <- 0.01  # space difference
s <- seq(-1 + ds/2, 1 - ds/2, by = ds) # result in a set of centroids
df <- data.frame(s = s, area = ds)

n1 <- n2 <- nrow(df)  # process 
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

source("Matern_32.R")
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


##----------------------
# use Fun_MakeSigmaTst.R
##----------------------
# to generate the full Sigma and test its sym and pd in one go

source("Fun_MakeSigmaTst.R")
Sigma_lst <- Make_Sigma_Tst(Var1 = sigma2_11, Var2 = sigma2_21, 
               Kappa1 = kappa1, Kappa2 = kappa2, 
               d = Dvec, B = B, Test = T)

str(Sigma_lst)
# List of 5
#$ Full_Sigma: num [1:400, 1:400] 1 0.974 0.91 0.827 0.736 ...
#$ C11       : num [1:200, 1:200] 1 0.974 0.91 0.827 0.736 ...
#$ C12       : num [1:200, 1:200] 0 0 0 0 0 0 0 0 0 0 ...
#$ C21       : num [1:200, 1:200] 0 0.000215 0.00104 0.002809 0.005788 ...
#$ C22       : num [1:200, 1:200] 0.2 0.1653 0.1116 0.0685 0.0398 ...



##*****************
## To extract individual marginal and cross-covariance functions at 
## the midpoint location of D, i.e. s = 0

# to find the mid-point of D which is the 0 btw [-1,1]
D[n1/2, n1/2]  # [1] 0  not D_vec

all(diag(D) == 0)  # [1] TRUE

## get full Sigma from above list
str(Sigma_lst[[1]])  # num [1:400, 1:400]
Sigma <- Sigma_lst[[1]]


# the individual marginal cov function and cross-cov function
Cov11 <- Sigma[n1/2, 1:n1] # cov(Y1(s_mid), Y1(si)), i = 1,...n
Cov12 <- Sigma[n1/2, (n1 + 1):n]
Cov21 <- Sigma[n1 + n2/2, 1:n1]
Cov22 <- Sigma[n1 + n2/2, (n1 + 1):n]
str(Cov11)  # num [1:200]


#=============================
# Plot the covariance function
#=============================
Cov_df <- expand.grid(s = df$s, proc1 = c("Y1", "Y2"), 
                      proc2 = c("Y1", "Y2"))

Cov_df$cov <- c(Cov11, Cov21, Cov12, Cov22) # df 800 rows
# the order of Cov21 and Cov12 matters
head(Cov_df)


g_cov <- LinePlotTheme() + 
  geom_line(data = Cov_df, aes(s, cov)) +
  facet_grid(proc1 ~ proc2)


if(print_figs) ggsave(g_cov,
                      filename = file.path(img_path, "cov_functions.png"),
                      width = 12, height = 10, family = "Arial")
# print_figs <- 0 so don't print, just save

# save first, then decide show or not
if(show_figs) print(g_cov, width = 12, height = 10)
print(g_cov)



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


#----------
# Cokrig
#---------

# both proc Y1 and Y2 are simulated from MVN(0, Sigma),
  # the cokrig of Y1 at s0, and Y2 at s0 for s0 in D is
  # via simple cokrig

head(df) 
# s, area, Y1, Y2, Z1, Z2

co_krig <- function(df, A, delta, r, obs_ind, name = NULL) {
  
  ## Form B matrix 
  B <- A * bisq_1d(d = H, delta = delta, r = r) *ds
  
  
  ## Construct full Sigma
  source("Fun_MakeSigmaTst.R")
  Sigma_lst <- Make_Sigma_Tst(Var1 = sigma2_11, Var2 = sigma2_21,
                 Kappa1 = kappa1, Kappa2 = kappa2,
                 d = Dvec, B = B)
  Sigma <- Sigma_lst$Full_Sigma
  
  
  ## Compute precision, only include where obs is available
  source("SparseMat_Imat.R")
  Q <- solve(Sigma[obs_ind, obs_ind] + sigmav^2 * Imat(length(obs_ind)))
  # 300
  
  ## cokrig formula
  mu <- Sigma[, obs_ind] %*% Q %*% Z[obs_ind, ] # 400*300 %*% 300*300 %*% 300*1 = 400*1
  sd <- diag(Sigma - Sigma[, obs_ind] %*% Q %*% t(Sigma[, obs_ind]))  # diag(400 * 400) = 400
  
  
  ## save results
  df[paste0(name, "_mu1")] <- mu[1:n1]
  df[paste0(name, "_mu2")] <- mu[-(1:n1)]
  df[paste0(name, "_sd1")] <- sd[1:n1]
  df[paste0(name, "_sd2")] <- sd[-(1:n1)]
  
  
  ## return 
  df
}


#--------------------
# Specify obs to keep
#--------------------
keep_z1 <- 101:200
keep_z2 <- 1:200

df$keep_Z1 <- 1:nrow(df) %in% keep_z1 
# $ keep_Z1  : logi  FALSE FALSE FALSE
df$keep_Z2 <- 1:nrow(df) %in% keep_z2

obs_ind <- c(keep_z1, keep_z2 + n1) # lower z, 101-200-201-400
# int [1:300] 101 102 103

str(obs_ind)


#-----------------------
# 3 senarios for cokrig
#-----------------------

# 1. Kriging predicator using only Z1 (Y1*): A set to 0; Y1 and Y2 indepandent
# 2. Cokrig predictor using both Z1 and Z2 under misspecified model (Y1@):
    # A and r are found using maximum likelihood with Delta = 0 (no asymmetry)
    # although we know true process has asymetry with delta = -.3
# 3. Cokrig using both Z1 and Z2 under true model (Y1_true): A and r are fixed to their true values

# For 1, equivalent to simple univariate kriging
# For 2, use log-likelihood to find the parameters by optimization 1sr and then cokrig
  # the optimized parameters are stored in "non_symm_par"


df2 <- subset(df, s > 0)
str(df2) # 'data.frame':	100 obs. of  8 variables:
  # $ keep_Z1  : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
  # $ keep_Z2  : logi  TRUE TRUE TRUE TRUE TRUE TRUE

H2 <- t(outer(df2$s, df2$s, FUN = "-")) 
D2 <- abs(H2)
Dvec2 <- as.double(c(D2))

Z_2 <- matrix(c(df2$Z1, df2$Z2))
str(Z_2) # num [1:200, 1] 



#---------------------
# define loglik_model
#---------------------
loglik_Model <- function(theta, model_num, i = NULL) {
  # theta[1] <- A
  # theta[2] <- r
  ## theta[3] <- delta, below cokrig set to fix 0, misspecified
  # fitted using only postive data Z1, Z2
  
  ## subset data only on postive domain
  df2 <- subset(df, s > 0)
  
  ## find displacement, each row is si with the rest of that row
  H2 <- t(outer(df2$s, df2$s, FUN = "-")) 
  
  ## find distance, distance can only be positive
  D2 <- abs(H2) 
  
  ## vectorize the distance for further use
  Dvec2 <- as.double(c(D2))  
  
  ## concatenate new obs data to get the parameters
  Z_2 <- matrix(c(df2$Z1, df2$Z2)) #str(Z_2) num [1:200, 1] 
  
  
  if(theta[2] < 0.005) { # do not allow aperture too small
    return(Inf)
  } else {
    ## construct B for further use
    ## 1. bisq(.) feed into displacement H, both +/- 
    ## differ from cov matrix make fun feed into distance D, only +
    ## 2. misspecified model with delta = 0
    B <- theta[1] * bisq_1d(A = 1, d = H2, r = theta[2], delta = theta[3]) * ds
    
    ## Construct full Sigma = Sigma_Y + Sigma_0 
    source("Fun_MakeSigmaTst.R")
    Sigma_Y <- Make_Sigma_Tst(Var1 = sigma2_11, Var2 = sigma2_21,
                              Kappa1 = kappa1, Kappa2 = kappa2, 
                              d = Dvec2, B = B)$Full_Sigma
    
    
    source("SparseMat_Imat.R")
    Sigma <- Sigma_Y +  sigmav^2 * Imat(nrow(df2) * 2)  #nrow(df2) = 100
    
    
    ## Cholesky the full Sigma
    cholZ <- chol(Sigma)
    
    ## neg logliklihood
    source("log_det.R") # fast compute log(det(Sigma))
    
    neg_loglik <- -(- 0.5 * nrow(Z_2) * log(2 * pi) - 
                      0.5 * log_det(cholZ) - 
                      0.5 * t(Z_2) %*% chol2inv(cholZ) %*% Z_2) %>% as.numeric()
    
    return(neg_loglik)
    
  }
}


#----------------
## optimization
#----------------
# to get the parameters
non_symm_par <- optim(par = c(1, 1, 1), # initial guess
      fn = loglik_Model, 
      hessian = F, 
      control = list(trace = 6, 
                     pgtol = 0, 
                     maxit = 3000))$par
str(non_symm_par)
# num [1:2] 4.266 0.802

# num [1:3] 5.184 0.61 -0.153


##----------------------------------
## cokrig with estimated parameters
##----------------------------------

## since we've got the parameters for A and r, 
  ## can now carry out (co)kriging

## can pipe original df through co_krig(.) using different
  ## values of A, delta, r

df <- df %>%
  co_krig(A = 0, delta = 0, r = r, obs_ind = obs_ind, name = "indp_model") %>%
  co_krig(A = non_symm_par[1], delta = 0, r = non_symm_par[2], obs_ind = obs_ind, name = "sym_model") %>%
  co_krig(A = A, delta = delta, r = r, obs_ind = obs_ind, name = "true_model") %>%
  co_krig(A = non_symm_par[1], delta = non_symm_par[3], r = non_symm_par[2], obs_ind = obs_ind, name = "asym_model")


str(df)
head(df)


##=======
## Plot
##=======

t <- df %>%
  dplyr::select(s, Z1, Z2, keep_Z1, keep_Z2) %>%
  gather(obs, z, Z1:Z2)
#         s keep_Z1 keep_Z2 obs          z
#1  -0.995   FALSE    TRUE  Z1  1.9798787
#2  -0.985   FALSE    TRUE  Z1  2.2279586
#3  -0.975   FALSE    TRUE  Z1  1.8493380
#4  -0.965   FALSE    TRUE  Z1  0.5935049
#5  -0.955   FALSE    TRUE  Z1  0.6083048
#6  -0.945   FALSE    TRUE  Z1  0.1505113
#7  -0.935   FALSE    TRUE  Z1  1.1363691
#8  -0.925   FALSE    TRUE  Z1  0.7891095

df_obs <- df %>% 
  dplyr::select(s, Z1, Z2, keep_Z1, keep_Z2) %>%
  gather(obs, z, Z1:Z2) %>%
  filter((keep_Z2 == T & obs == "Z2") | (keep_Z1 == T & obs == "Z1"))

str(df_obs)
# 'data.frame':	300 obs. of  5 variables:
#$ s      : num  0.005 0.015 0.025 0.035 0.045 ...
#$ keep_Z1: logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#$ keep_Z2: logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#$ obs    : chr  "Z1" "Z1" "Z1" "Z1" ...
#$ z      : num  0.9264 0.0143 2.0957 0.1707 0.5508 ...


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
t <- df %>%
  dplyr::select(s, sample_Y1, indp_model_mu1, sym_model_mu1, true_model_mu1) %>%
  gather(process, z, sample_Y1, indp_model_mu1, sym_model_mu1, true_model_mu1)
# sample_Y1 is the ture hidden process
# the rest 3 are results from co-krig function
# altogether 800

t <- df %>% 
  dplyr::select(s, sample_Y1, indp_model_mu1, sym_model_mu1, true_model_mu1) %>%
  gather(process, z, sample_Y1, indp_model_mu1, sym_model_mu1, true_model_mu1) %>%
  mutate(group = substr(process, 1, 4))
  # subtr(x, start, stop) extract/replace  substrings in a char vector

head(t, 30)
tail(t, 30)
str(t)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

df_estY1 <- df %>% 
  dplyr::select(s, sample_Y1, indp_model_mu1, sym_model_mu1, true_model_mu1, asym_model_mu1) %>%
  gather(process, z, sample_Y1, indp_model_mu1, sym_model_mu1, true_model_mu1, asym_model_mu1) %>%
  mutate(group = substr(process, 1, 4))


df_estY2 <- df %>%
  dplyr::select(s, sample_Y2, indp_model_mu2, sym_model_mu2, true_model_mu2, asym_model_mu2) %>%
  gather(process, z, sample_Y2, indp_model_mu2, sym_model_mu2, true_model_mu2, asym_model_mu2) 
  #mutate(group = subtr(process, 1, 4))



#----------------
# ggplot for obs
#----------------
## for obs ##
head(df_obs)
#      s keep_Z1 keep_Z2 obs          z
#1 0.005    TRUE    TRUE  Z1 0.92638116
#2 0.015    TRUE    TRUE  Z1 0.01431249

obs_plt <- LinePlotTheme() +
  geom_point(data = df_obs, aes(x = s, y = z, shape = obs),
             size = 3, alpha = 1) +
  theme(legend.title = element_blank(),
        plot.margin = grid::unit(c(3, 0, 0, 0), units = "mm")) + 
  scale_shape_manual(values = c(1, 20), # different shape of pts
                     labels = c(expression(Z[1]), expression(Z[2]))) +
  ylab("Z")

plot(obs_plt)


##----------------------------------------------------------------
# ggplot for 3 senarios predictor for Y1 + sampling process of Y1
##----------------------------------------------------------------
## for 3 senarios for Y1 + true Y1 ##
str(df_estY1)
df_estY1$process <- as.factor(df_estY1$process)
str(df_estY1)
#'data.frame':	800 obs. of  4 variables:
#$ s      : num  -0.995 -0.985 -0.975 -0.965 -0.955 -0.945 -0.935 -0.925 -0.915 -0.905 ...
#$ process: Factor w/ 4 levels "indp_model_mu1",..: 2 2 2 2 2 2 2 2 2 2 ...

df_estY1$process <- relevel(df_estY1$process, ref = 3) # reorder the 2nd level 
# $ process: Factor w/ 4 levels "sample_Y1","indp_model_mu1",..: 1 1 1 1 1 1 1 1 1 1 ...


est_plotY1_no_CIs <- LinePlotTheme() + 
  geom_line(data = df_estY1, 
            aes(s, z, colour = process,linetype = process,size = process)) +
  theme(legend.title = element_blank(),
        plot.margin = grid::unit(c(3, 0, 0, 0), unit = "mm")) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "twodash"),
                        labels = c(expression(Y[1]),
                                   expression(tilde(Y)[1]),
                                   expression(Y[1]^"degga"),
                                   expression(hat(Y)[1]),
                                   expression(Y[1]^"mle"))) + 
  scale_size_manual(values = c(0.4, 0.8, 0.8, 0.8, 0.5), guide = 'none') +
  scale_color_manual(values = c("black", "black", "black", "black", "grey"), guide = "none", name = "") +
  ylab("Y")

plot(est_plotY1_no_CIs)
## mle estimated lline (grey twodashed) approximates the true sample

est_plotY1 <- est_plotY1_no_CIs + 
  geom_ribbon(data = df, aes(s, ymax = indp_model_mu1 + indp_model_sd1,
                             ymin = indp_model_mu1 - indp_model_sd1),
              alpha = .2, fill = "red") +
  geom_ribbon(data = df, aes(s, ymax = true_model_mu1 + true_model_sd1,
                             ymin = true_model_mu1 - true_model_sd1),
              alpha = .2, fill = "black") +
  geom_ribbon(data = df, aes(s, ymax = sym_model_mu1 + sym_model_sd1,
                             ymin = sym_model_mu1 - sym_model_sd2),
              alpha = .2, fill = "green") +
  geom_ribbon(data =df, aes(s, ymax = asym_model_mu1 + asym_model_sd1,
                            ymin = asym_model_mu1 - asym_model_sd1),
              alpha = .2, fill = "blue")



plot(est_plotY1)


## df_estY2

est_plotY2 <- LinePlotTheme() + 
  geom_line(data = df_estY2, 
            aes(s, z, color = process, linetype = process, size = process)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "twodash"),
                        guide = "none") + 
  scale_size_manual(values = c(1, 1.3, 1.3, 1.3, 0.5), guide = "none") + 
  scale_color_manual(values = c("black", "red", "green", "blue", "grey"),
                     labels = c("Y2", "IM", "TM"), 
                     name="") +
  ylab("")
  
plot(est_plotY2)


##------------
# save & print
##------------
# cokriging using spatial covariance defined by conditional approach
# Top panel: The simulated obs Z1 (open circle), and Z2 (dots);
# Bottom panel: The hidden true process value Y1 (solid line) and 
  # the kriging predictor tilde(Y1) (dashed)
  # the misspecified (delta = 0) cokriging dagger(Y1) (dotted)
  # the true cokrig predictor with parameters fixed to true value hat(Y1) (dotted dash))
  # the mle delta = -.15 cokringing Y^(mle) (twodashed)

# Prediction error intervals are shown using different shadings.
if(print_figs) ggsave(obs_plt, 
                      filename = file.path(img_path, "sim_obs.eps"),
                      width = 7, height = 4)

if(print_figs) ggsave(est_plotY1_no_CIs,
                      filename = file.path(img_path, "sim_est_no_CIs.eps"),
                      width = 7, height = 4)


if(print_figs) ggsave(est_plotY1,
                      filename = file.path(img_path, "sim_est_Y1.eps"),
                      width = 7, height = 4)

if(print_figs) ggsave(est_plotY2,
                      filename = file.path(img_path, "sim_est_Y2.eps"),
                      width = 7, height = 4)

if(show_figs) grid.arrange(obs_plt, est_plotY1, ncol = 1)



#============================================
# Plot the covariance and cross - covariance 
#============================================
# obtained from makeSY or make_Sigma_Tst

t <- expand.grid(s1 = df$s, comp1 = c("Y1", "Y2"), 
            s2 = df$s, comp2 = c("Y1", "Y2"))
# produce df containing all combinations of vars
  # with the 1st col runs fastest, 2nd col 2nd fast,...

Sigma_df <- expand.grid(s1 = df$s, comp1 = c("Y1", "Y2"),
            s2 = df$s, comp2 = c("Y1", "Y2")) %>%
  mutate(cov = c(Sigma))

str(Sigma_df) #'data.frame':	160000 obs. of  5 variables:

Sigma_plot <- LinePlotTheme() + 
  geom_tile(data = Sigma_df, aes(s2, s1, fill = cov)) + 
  facet_grid(comp1 ~ comp2) + 
  scale_fill_gradient(low = "white", high = "black") +
  coord_fixed() +
  ylab("s") + xlab("u") + scale_y_reverse() + 
  theme(panel.spacing = grid::unit(1, "lines"))


if (print_figs) ggsave(Sigma_plot,
                       filename = file.path(img_path, "Sigma.eps"),
                       width = 8, height = 7)

if (show_figs) print(Sigma_plot, width = 16, height = 7, family = "Arial")


if (print_figs) {
  g_all <- grid.arrange(Sigma_plot, 
               arrangeGrob(obs_plt, est_plotY1_no_CIs, ncol = 1),
               ncol = 2, widths = c(1, 1))
  
  cairo_ps(filename = file.path(img_path, "Fig1.eps"),
           width = 16, height = 7)
  
  grid.draw(g_all)  # draw a grid grob
  #dev.off()
  
}



??grid.draw()





































