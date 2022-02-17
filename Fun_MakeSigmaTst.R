##******************************
## Function Make Sigma and Test
##******************************

# the function makes the full sigma for bivariate process Y1 and Y2

# the univariate covariance matrix on the diag block is 
  # modelled as stationary Matern with nu = 3/2
  # these block matrices are required to by symmetric and pd

# the off-diag covariance matrix is obtained based on cressie's bincon method
  # required to be transpose to each other

# the full joint var-cov matrix is required to be symmtry and pd

# if Test set to TRUE, then the symmetry and pd properties will be
  # tested for diag block matrices and the full Sigma


Make_Sigma_Tst <- function(Var1, Var2, Kappa1, Kappa2, d, B, Test = FALSE) {
  source("Matern_32.R")
  
  C11 <- Matern_32(Var = Var1, Kappa = Kappa1, d)
  C2_1 <- Matern_32(Var = Var2, Kappa = Kappa2, d)
  
  C12 <- C11 %*% t(B)
  C21 <- B %*% C11
  
  C22 <- C2_1 + B %*% C11 %*% t(B)
  
  Sigma <- rbind(cbind(C11, C12), cbind(C21, C22))
  
  
  if (Test == TRUE) {
    source("fun_test_sym_pd.R")
    
    print("Test the sysmmetry and positive definite of C11: ")
    Test_sym_pd(C11)
    
    print("Test the sysmmetry and positive definite of C22: ")
    Test_sym_pd(C22)
    
    print("Test the sysmmetry and positive definite of Sigma: ")
    Test_sym_pd(Sigma)
  }
  
  list(Full_Sigma = Sigma, 
             C11 = C11, C12 = C12, 
             C21 = C21, C22 = C22)
}