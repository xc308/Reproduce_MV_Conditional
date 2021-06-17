#==================================================
# dgTMatrix-class: Sparse matrices in triplet form
#==================================================

# The "dgTMatrix" class is the class of sparse matrices stored as (possibly redundant) triplets.

# Objects can be created by calls of the form new("dgTMatrix", ...), 
# but more typically via as(*, "dgTMatrix"), spMatrix(), or sparseMatrix(*, repr = "T").


# slots:
  #@ i: row indices of non-zero entries in 0-base for each column, i.e., must be in 0:(nrow(.)-1).
  #@ j: column indices of non-zero entries. Must be the same length as slot i and 0-based as well, i.e., in 0:(ncol(.)-1).
  #@ x: numeric vector - the (non-zero) entry at position (i,j). Must be the same length as slot i
  #@ Dim: the dimensions of the matrix.


## Hence the triplet representations:
# the simplest representation of a sparse matrix is as a triplet 
# of an integer vector i giving the row numbers,
# an integer vector j giving the column numbers, 
# and a numeric vector x giving the non-zero values in the matrix



# Examples:
m <- Matrix(0+1:28, nrow = 4)
m

m[-3,c(2,4:5,7)] <- m[ 3, 1:4] <- m[1:3, 6] <- 0
m
#4 x 7 Matrix of class "dgeMatrix"
#[,1] [,2] [,3] [,4] [,5] [,6] [,7]
#[1,]    1    0    9    0    0    0    0
#[2,]    2    0   10    0    0    0    0
#[3,]    0    0    0    0   19    0   27
#[4,]    4    0   12    0    0   24    0


mT <- as(m, "dgTMatrix")
mT 

#4 x 7 sparse Matrix of class "dgTMatrix"

#[1,] 1 .  9 .  .  .  .
#[2,] 2 . 10 .  .  .  .
#[3,] . .  . . 19  . 27
#[4,] 4 . 12 .  . 24  .

str(mT)
mT[4, drop = FALSE]


identical(mT[lower.tri(mT)], m[lower.tri(m)]) # [1] TRUE

mT[lower.tri(mT, diag = T)]
mT[lower.tri(mT, diag = T)] <- 0
mT


## Triplet representation with repeated (i,j) entries

T2 <- new("dgTMatrix",
          i = as.integer(c(1,1,0,3,3)),
          j = as.integer(c(2,2,4,0,0)), 
          x = 10 * (1:5), 
          Dim = 4:5)

T2

nnzero(T2) # number of non-zero values
# 3



T3 <- new("dgTMatrix",
    i = as.integer(c(0, 1, 3, 5)),
    j = as.integer(c(2, 1, 4, 1)),
    x = as.double(2:5), 
    Dim = 6:5)

T3



