#------------
# Test symmetric and pd of a matrix with option print or not
#------------
# Aim:
    # This function tests the symmetricity and p.d of a matrix

# Funtion and Results:
    # when Print set to default T
      # print out the test result
    # When Print is set to F
      # collect all the test results into a data frame with two variables
      # symmetric and pd

# Arguments:
    # test.mat: matrix to be tested
    # Print: or not

Test_sym_pd_V2 <- function(test.mat, Print = T) {
  
  Sym <- NULL
  PD <- NULL
  
  if (Print) {
    if (isSymmetric(test.mat) == TRUE) {
      print("Symmetric: Yes")
    } else print("Symmetric: No")
    
    if (all(eigen(test.mat)$values > 0)) {
      print("p.d.: Yes")
    } else print("p.d.: No")
  } else {
    if (isSymmetric(test.mat) == TRUE) {
      Sym <- rbind(Sym, 1)
    } else {Sym <- rbind(Sym, 0)}
    
    if (all(eigen(test.mat)$values > 0)) {
      PD <- rbind(PD, 1)
    } else{ PD <- rbind(PD, 0)}
    
    data.frame(Symmetric = Sym, PD = PD)
  }
}



#Res <- Test_sym_pd_V2(Sigma, Print = F)

#str(Res)
# 'data.frame':	1 obs. of  2 variables:
# $ Symmetric: num 1
# $ PD       : num 1

#RES <- Test_sym_pd_V2(Sigma, Print = T)
RES <- Test_sym_pd_V2(Sigma, Print = F)
rm(RES)
