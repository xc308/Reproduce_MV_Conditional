#=====================================
# Subset cols of matrix according to t
#=====================================

# when lenght(t) > 1
# want to subset B and C matrices in R

A <- matrix(c(1, 2, 3, 4), nrow = 2)
B <- matrix(c(5, 6, 7, 8), nrow = 2)
C <- matrix(c(9, 10, 11, 12), nrow = 2)

R <- cbind(A, B, C)

t <- c(2, 3)
n <- nrow(A)

t <- 2

Subset_cols <- function(t) {
  start_col <- (t - 1) * n + 1
  end_col <- t * n
  
  result <- R[, start_col:end_col]
  
}

result_lst <- lapply(t, FUN = Subset_cols)
do.call(cbind, result_lst)





