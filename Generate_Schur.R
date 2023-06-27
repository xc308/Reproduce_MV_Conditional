#===========================
# Generate Schur complement 
#===========================
# Aim: 
  # This function generates Schur complement for arbitray rows and cols of a matrix

# Args:
  # matrixA : the matrix to be partitioned
  # row_ind, col_ind: indices from 2:(n-1), n = nrow(matrixA)
    # such that the partition at least be A[1:2, 1:2]
      # at most be A[1:(n-1), 1:(n-1)]


Generate_Schur <- function(matrixA, row_ind, col_ind) {
  
  A11 <- matrixA[row_ind, col_ind]
  A12 <- matrixA[row_ind, -(col_ind)]
  A21 <- matrixA[-(row_ind), col_ind]
  A22 <- matrixA[-(row_ind), -(col_ind)]
  
  A11_inv <- solve(A11)
  Schur <- A22 - A21 %*% A11_inv %*% A12
  
  return(Schur)
  #source("fun_test_sym_pd.R")
  #Test_sym_pd(Schur)
  
}
