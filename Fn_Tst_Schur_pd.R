#========================================
# Test the p.d. of any schur complements
#========================================
# Aim: This function test the p.d. of Schur complements of a given block matrix

# Resons: 
    # The sufficient conditions for a block matrix to be p.d. are:
        # 1. all its block diagonal matrices must be p.d.
        # 2. all of its Schur complents must be p.d.

# Arguments: 
# idex_choice: the row and col indices to form the sub-block matrics and 
# corresponding Schur complements of MatrixA_tst;
# MatrixA_tst: a block matrix whose Schur complements to be tested.

Tst_Schur_pd <- function(indx_choice, matrixA_tst) {
  
  n <- nrow(matrixA_tst)
  if (any(indx_choice >= n) || any(indx_choice <= 1)) {
    print("Each indx_choic element must be in (1, n), n = nrow(matrixA_tst).")
  } else {
    C <- length(indx_choice)
    R <- data.frame()
    
    for (c in 1:C) {
      # row and col ind must be the same to form sub block diag matrices
      row_ind <- indx_choice
      col_ind <- indx_choice
      
      source("Fn_Generate_Schur.R")
      Schur <- Generate_Schur(Sigma_M, row_ind = 1:row_ind[c], col_ind = 1:col_ind[c])
      
      #print(paste0("row, col index ", "1:", indx_choice[c]))
      
      # Test p.d. of Schur
      #source("fun_test_sym_pd.R")
      #Test_sym_pd(Schur)
      cat(c, "\n")
      
      source("Fn_test_sym_pd_v2.R")
      r <- Test_sym_pd_V2(Schur, Print = F)  # a df sym, pd
      res <- cbind(indx = paste0("1:", indx_choice[c]), r)
      
      R <- rbind(R, res)
    }
    
    return(R)
  }
}


#--------
# Examples
#--------
#saveRDS(Sigma_M, "Sigma_M.rds")
#readRDS("Sigma_M.rds")

#Sigma_M : 400 * 400

## EX1:
#choice <- c(2, 20, 100, 200, 250, 300, 350, 398)
#Res_choice1 <- Tst_Schur_pd(indx_choice = choice, matrixA_tst = Sigma_M)
#Res_choice1
#str(Res_choice1)
  # 'data.frame':	8 obs. of  3 variables:

#   indx Symmetric PD
#1   1:2         1  1
#2  1:20         0  1
#3 1:100         0  1
#4 1:200         0  1
#5 1:250         0  1
#6 1:300         0  1
#7 1:350         0  1
#8 1:398         0  1


## EX2:
#choice2 <- seq(2, 399, by = 1)
#Res_choice2 <- Tst_Schur_pd(indx_choice = choice2, matrixA_tst = Sigma_M)
#str(Res_choice2) 
# data.frame':	398 obs. of  3 variables:
#$ indx     : chr  "1:2" "1:3" "1:4" "1:5" ...
#$ Symmetric: num  1 1 0 0 0 0 0 0 0 0 ...
#$ PD       : num  1 1 1 1 1 1 1 1 1 1 ...

#all(Res_choice2$PD == 1)
# [1] TRUE




