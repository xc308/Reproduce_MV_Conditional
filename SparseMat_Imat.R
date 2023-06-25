#=======================#
# Sparse Identity Matrix
#=======================#

# - Imat
# - Sparse matrices are used to efficiently represent matrices 
    # that have a large number of zero values.


library(Matrix)
Imat <- function(n) {
  sparseMatrix(i = 1:n, j = 1:n, x = 1)
}











