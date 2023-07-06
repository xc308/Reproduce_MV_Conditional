#============================
# Effects of different kappa
#============================
ds = 0.01

A <- 5
delta <- -0.3
r <- 0.3

sigma2_11 <- 1       # marginal variance for C11
sigma2_21 <- 0.25    # marginal variance for C2|1

nu11 = nu2_1 = 3/2


H <- t(outer(df$s, df$s, "-"))
Dvec <- as.double(c(abs(H)))


#---------------------------------
# Make Sigma for different Kappas
#---------------------------------

make_Sigma <- function(Kappa11, Kappa21) {
  
  # construct b(;) and B = b(;) * ds
  B = A * b_1d(d = H, A = 1, delta = delta, r = r) * ds
  
  # construct C11, C12, C21, C22
  source("Matern_32.R")
  C11 <- Matern_32(Var = sigma2_11, Kappa = Kappa11, d_vec = Dvec)
  C2_1 <- Matern_32(Var = sigma2_21, Kappa = Kappa21, d_vec = Dvec)
  
  C12 <- C11 %*% t(B)
  C21 <- t(C12)
  C22 <- C21 %*% t(B) + C2_1
  
  Sig <- rbind(cbind(C11, C12), cbind(C21, C22))
  
  return(Sig)
}


#-------------------------
# Plot each block of sigma 
#-------------------------

#----------
# Option 1
#----------
plt_Sig <- function(Sigma) {
  
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma)
}

plt_Sig(Sigma)


#----------
# Option 2
#----------
plt_Sigma <- function(Sigma) {
  Sigma_df <- expand.grid(s1 = df$s, comp1 = c("Y1", "Y2"), 
                          s2 = df$s, comp2 = c("Y1", "Y2")) %>%
    mutate(cov = c(Sigma))
  
  Sigma_plt <- ggplot() + 
    geom_tile(data = Sigma_df, aes(x = s1, s2, fill = cov)) + 
    facet_grid(comp1 ~ comp2) +
    scale_fill_gradient(low = "white", high = "black") +
    ylab("s") + xlab("u") + 
    scale_y_reverse()
  
  print(Sigma_plt)
}

print(plt_Sigma(Sigma))



#==========================
# Try on different kappas
#==========================

kappa_11 <- 25
kappa_21 <- 75

Sigma_1 <- make_Sigma(Kappa11 = 25, Kappa21 = 75)
plt_Sig(Sigma_1)  # correct
title(main = paste0("Kappa11= ", kappa_11, ", Kappa21= ", kappa_21),
      cex.main = 0.8)

print(plt_Sigma(Sigma_1))


#===========================================
# Joint function reflect Kappa in main title
#===========================================

make_plt_Sigma <- function(Kappa11, Kappa21) {
  
  ## build Sigma  
  B = A * b_1d(d = H, A = 1, delta = delta, r = r) * ds
  
  # construct C11, C12, C21, C22
  source("Matern_32.R")
  C11 <- Matern_32(Var = sigma2_11, Kappa = Kappa11, d_vec = Dvec)
  C2_1 <- Matern_32(Var = sigma2_21, Kappa = Kappa21, d_vec = Dvec)
  
  C12 <- C11 %*% t(B)
  C21 <- t(C12)
  C22 <- C21 %*% t(B) + C2_1
  
  Sigma <- rbind(cbind(C11, C12), cbind(C21, C22))
  
  
  ## Plot Sigma
  # to reverse order of cols in Sigma
  rev_Sigma <- Sigma[, ncol(Sigma):1]
  
  # Plot the matrix with reversed y-axis scale
  image(1:nrow(Sigma), 1:ncol(Sigma), rev_Sigma, 
        xlab = "s", ylab = "u")
  title(main = paste0("Kappa11= ", Kappa11, ", Kappa21= ", Kappa21),
        cex.main = 0.8)
  
  # rename the axis of image
  #x_lab <- "Locations s"
  #y_lab <- "Locations u"
  #title(xlab = x_lab)
  #title(ylab = y_lab)
  
}


#=====================
# Try on real numbers
#=====================

Kappa11 <- 50
Kappa21 <- 50


Kappa11 <- 5
Kappa21 <- 5

Kappa11 <- 2.5
Kappa21 <- 2.5


Kappa11 <- 75  # very rough
Kappa21 <- 2.5 # very smooth

Kappa11 <- 2.5 # rho is large, corre large, very smooth
Kappa21 <- 75  # rho is very small, corr small, very rough


Kappa11 <- 75  # very rough
Kappa21 <- 1.5 # very smooth

Kappa11 <- 1.5 # rho is large, corre large, very smooth
Kappa21 <- 75  # rho is very small, corr small, very rough


#-----------
# save plots
#-----------

image_path <- "./Results/Fig/"

png(paste0(image_path, "different_kappas.png"), width = 8, height = 7, units = "in", res = 300)
par(mfrow = c(1, 2))

Kappa11 <- 75  # very rough
Kappa21 <- 1.5 # very smooth
make_plt_Sigma(Kappa11 = Kappa11, Kappa21 = Kappa21)

Kappa11 <- 1.5 # rho is large, corre large, very smooth
Kappa21 <- 75  # rho is very small, corr small, very rough
make_plt_Sigma(Kappa11 = Kappa11, Kappa21 = Kappa21)

dev.off()








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

Sigma_plot <- ggplot() + 
  geom_tile(data = Sigma_df, aes(s2, s1, fill = cov)) + 
  facet_grid(comp1 ~ comp2) + 
  ylab("s") + xlab("u") + scale_y_reverse()  
  
print(Sigma_plot)








make_plt_Sigma <- function(Kappa11, Kappa21) {
  
  # construct b(;) and B = b(;) * ds
  B = A * b_1d(d = H, A = 1, delta = delta, r = r) * ds
  
  # construct C11, C12, C21, C22
  C11 <- Matern_32(Dvec, sigma2 = sigma2_11, kappa = Kappa11, nu = nu11)
  C2_1 <- Matern_32(Dvec, sigma2 = sigma2_21, kappa = Kappa21, nu = nu2_1)
  
  C12 <- C11 %*% t(B)
  C21 <- t(C12)
  C22 <- C21 %*% t(B) + C2_1
  
  Sig <- rbind(cbind(C11, C12), cbind(C21, C22))
  
  image(Sig, main = paste0("kappa1: ", Kappa11, ", kappa2: ", Kappa21))
}



kappa_11 <- 25
kappa_21 <- 75

make_plt_Sigma(Kappa11 = kappa_11, Kappa21 = kappa_21)

kappa_11 <- 1/25     # 0.04
kappa_21 <- 1/75     # 0.01333333


