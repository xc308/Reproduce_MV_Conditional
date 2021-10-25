#=======================#
# Sparse Identity Matrix
#=======================#

# - Imat

Imat <- function(n) {
  sparseMatrix(i = 1:n, j = 1:n, x = 1)
}











