#*************************
# Topics: Sparse Matrix
#*************************
# 1. Basics
# 2. Values in sparse matrix replacement


library(Matrix)


#====================================
# 1. Basics: 
#    Sparse Matrix vs regular matrix
#====================================
# Ref: https://thecoatlessprofessor.com/programming/r/benefits-of-using-a-sparse-matrix/

# use binorm dist (p = .06) to generate entries for matrix
standard_matrix = matrix(
  rbinom(100*100, size = 1, prob = .06),
  nrow = 100
)

print(object.size(standard_matrix), units = "Kb")

standard_matrix = matrix(
  rbinom(100*100, size = 1, prob = .06), 
  nrow = 100
)
# 39.3 Kb

standard_matrix[1:5, 1:5]


## convert regular to sparse matrix
sparse_matrix <- as(standard_matrix, "sparseMatrix")
sparse_matrix[1:5, 1:5]

#5 x 5 sparse Matrix of class "dgCMatrix"

#[1,] . . . . .
#[2,] . . . . .
#[3,] . . . . .
#[4,] . . . . 1
#[5,] . . . . .

# As a result of discarding zero entries, 
# the size of the matrix decreases greatly.

print(object.size(sparse_matrix), units = "Kb")
# 8.9 Kb
# vs # 39.3 Kb standard_mat


# because sparseMatrix only save coords for non-zero entries
# i is row index;
# j is column index;
# x is the value (optional if values other than 1 preferred); and
# dim is the dimension of the matrix.


# with scheme at hand, one can construct sparseMatrix directly
# using sparseMatrix() in Matrix package
sample_sparse_matrix = sparseMatrix(
  i = c(1, 4, 6),
  j = c(2, 3, 5), 
  x = c(rbinom(3, 10, p = .06)),
  dims = list(10, 10))

sample_sparse_matrix


#---------------------------------------
# Alternative way to create sparseMatrix
#---------------------------------------
# Example 2
m <- Matrix(toeplitz(c(10, 0, 1, 0, 3)), sparse = TRUE)
m
toeplitz(c(10, 1, 0, 0, 0)) # tri-diagnoal
View(toeplitz)


#====================================
# 2. Sparse Matrix value replacement
#====================================

## Aim: row standardize the matrix 
# ref: https://stat.ethz.ch/pipermail/r-help/2010-December/262365.html

set.seed(14-04-2022)

M <- sparseMatrix( 
  # sample 1000 # from 1:5000 with replacement
  i = sample((1:5000), size = 1000, replace = T),
  j = sample((1:5000), size = 1000, replace = T),
  x = rnorm(1000),
  dims = c(5000, 5000)
)

str(M)
range(M@i) # [1]   15 4997

str(rs <- rowSums(M, na.rm = T))
# num [1:5000] 0 0 0 0 0 0 0 0 0 0 ...

range(rs)
head(rs[M@i + 1])

standardize_mat <- sparseMatrix(
  i = M@i, p = M@p, 
  x = M@x / rs[M@i + 1L], 
  index1 = F
  )

str(standardize_mat)

table(rowSums(standardize_mat))
#    0    1 
#  4104  894
# 894 rows that is 1 after standardzation

proc.time()


#-----------------------------------------
# Assign value to a colum of sparseMatrix
#-----------------------------------------
# ref: https://stackoverflow.com/questions/40907960/assign-value-to-r-sparse-matrix

str(M)
range(M@p)

# assign a specific value to a certain col
system.time(M[, 1] <- 0)
#   user  system elapsed 
#  0.004   0.001   0.008 


summary(M)
#5000 x 5000 sparse Matrix of class "dgCMatrix", with 1000 entries 
#       i    j            x
#1    350    2 -0.169810533
#2   1210    6  0.284641801
#3   3576   14 -0.200787354
#4   3529   15 -1.223419501
#5   3753   26 -0.478074405


# a bit slow, try to nullify (invalid) col from the dense storage encoding of the matrix 
nullify_col <- function(M, i) {
  M.dense <- summary(M)
  filter <- M.dense$j != i
  return(sparseMatrix(i = M.dense$i[filter],
               j = M.dense$j[filter],
               x = M.dense$x[filter]))
  
}


build_rd_sparse_mat <- function(n, p, q) {
  i <- sample(1:n, size = q, replace = T)
  j <- sample(1:p, size = q, replace = T)
  s <- rnorm(q)^2
  M <- sparseMatrix(i = i, j = j, x = s)
  return(M)
}


t0 <- Sys.time()  

n <- 1000
p <- 500
sparse.ratio <- 0.0001
q <- n*p*sparse.ratio

t1 <- Sys.time()
A <- build_rd_sparse_mat(n = n, p = p, q = q)
str(A)
B <- build_rd_sparse_mat(n = n, p = 1, q = q*2) # with more cols to be nullified for non-zero values
str(B)

see <- sample(1:1, q*2, replace = T)
head(see)
str(see)

M <- cbind(B, A) ##???
str(M)
t2 <- Sys.time()

time_diff <- round(as.numeric(difftime(t2, t1, units = "secs")), 2)
print(paste0("Build sparse matrix took ", time_diff, " seconds", sep = ""))

str(nullify_col(M, 1))


#---------------------------------------
# Speed up Sparse Matrix mulitplication
#---------------------------------------

# ref:https://stackoverflow.com/questions/53307953/speed-up-sparse-matrix-multiplication-in-r

install.packages("tictoc")
library(tictoc)

set.seed(14-04-2022)

S <- sample(1e4)
str(S)
# int [1:10000] 5260 8825 8616 8783 185 1272 2620 1860 6569 9953 ...


m <- Matrix(
  sample(c(0, 1), 
         size = length(S)^2, 
         replace = T, 
         prob = c(0.99, 0.01)),
  nrow = length(S), 
  ncol = length(S), 
  sparse = F)

m_sparse <- Matrix(m, sparse = T)

tic("dense")
x <- m %*% S
toc()
# dense: 24.694 sec elapsed
str(x)

tic("sparse")
y <- m_sparse %*% S
toc()
# sparse: 21.121 sec elapsed
str(y)


#---------------------------
# Sparse Matrix manipulation
#---------------------------
# ref: https://www.geeksforgeeks.org/working-with-sparse-matrices-in-r-programming/


sp_mat_1 <- Matrix(0, nrow = 1000, ncol = 1000, sparse = T)
str(sp_mat_1)
summary(sp_mat_1)

# setting the value at 1st row & 1st col to be 1
sp_mat_1[1][1] <- 5
summary(sp_mat_1)

print(paste0("size of sparse mat1: ", object.size(sp_mat_1), sep = ""))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Constructing Sparse Matrices From Dense
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Syntax: as(dense_matrix, type = )
  # type : Default evaluates to dgCMatrix, 
  # compressed sparse column( CSC ) format. 
  # The other type available is the dgRMatrix, 
  # which converts the dense matrix in sparse row format.


# create a dense matrix: 
  # sample 0 with prob 0.8
  # sample 6 with prob 0.1
  # sample 8 with prob 0.1

set.seed(0)
rows <- 4L
cols <- 6L
vals <- sample(
  c(0, 6, 8), 
  prob = c(0.8, 0.1, 0.1), 
  size = rows * cols, 
  replace = T)
  
dense_mat <- matrix(vals, nrow = rows)
print("dense matrix")
print(dense_mat)
print(paste0("dense matrix size: ", object.size(dense_mat), " bytes"))


## convert to sparse
sparse_mat <- as(dense_mat, "sparseMatrix")
print("sparse matrix")
print(sparse_mat)
print(paste0("sparse matrix size: ", object.size(sparse_mat), " bytes"))



#~~~~~~~~~~~
# operation
#~~~~~~~~~~~

# Addition and subtraction by Scalar Value
  # scalar values are added or subtracted to all the elements of the sparse matrix
  # resultant matrix is a dense matrix


rows <- 4L
cols <- 6L
vals <- sample(
  x = c(0, 9), 
  prob = c(0.85, 0.15),
  replace = T, 
  size = rows * cols
)

dense_mat2 <- matrix(vals, nrow = rows)

# convert to sparse
sparse_mat2 <- as(dense_mat2, "sparseMatrix")
print("sparseMatrix:")
print(sparse_mat2)

# addition
print(sparse_mat2 + 5)

# subtraction
print(sparse_mat2 - 1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multiplication or Division by Scalar
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# performed on all the non-zero elements of the matrix
  # resultant matrix is a sparse matrix

rows <- 4L
cols <- 6L
vals <- sample(
  x = c(0, 6),
  prob = c(0.85, 0.15),
  replace = T,
  size = rows * cols
)

dense_mat3 <- matrix(vals, nrow = rows)

# convert to sparse matrix
sparse_mat3 <- as(dense_mat3, "sparseMatrix")

# multiplication
sparse_mat3 * 10
object.size(sparse_mat3* 10)

# division
sparse_mat3 / 10


## ``````
## Matrices can be multiplied with each other, 
# irrespective of sparse or dense.
# columns of the first matrix should be equal to rows of the second
## ``````

mul_mat <- sparse_mat3 %*% t(sparse_mat3)
print(mul_mat)

Matrix(sparse_mat3, sparse = F)


##````````
## multiplication by vector
##````````
# the first row multiplied by the 1st element of vec
# the 2nd row by the 2nd element of vec
# till the end of length of vec
# then recycle

set.seed(0)
rows <- 4L
cols <- 6L
vals <- sample(
  x = c(0, 10), 
  prob = c(0.85, 0.15), 
  size = rows * cols, 
  replace = TRUE
)

dense_mat4 <- matrix(vals, nrow = rows)
sparse_mat4 <- as(dense_mat4, "sparseMatrix")

vec <- c(3, 2)
vec2 <- c(3, 2, 4)  # ??
vec3 <- c(3, 2, 4, 6)
vec4 <- c(3, 2, 4, 6, 7, 8) #?? rules??
sparse_mat4 * vec
sparse_mat4 * vec2
sparse_mat4 * vec3
sparse_mat4 * vec4


#`````````````````````````
# Combination of Matrices
#```````````````````````````

# matrices can be combined with vectors
# using column bind cbind( ) or row bind rbind( ) operations
# The resultant matrices rows are the summation of the rows of the input matrices in rbind() function 
# columns are the summation of the columns of the input matrices in cbind()

rows <- 4L
cols <- 6L
vals <- sample(
  x = c(0, 10),
  prob = c(0.8, 0.2),
  replace = T, 
  size = rows * cols
)

dense_mat5 <- matrix(vals, nrow = rows)
sparse_mat5 <- as(dense_mat5, "sparseMatrix")

print(sparse_mat5)
rbind(sparse_mat5, sparse_mat5)
# 8 x 6 sparse Matrix of class "dgCMatrix"

cbind(sparse_mat5, sparse_mat5)
# 4 x 12 sparse Matrix of class "dgCMatrix"

rbind(sparse_mat5, dense_mat5)
# 8 x 6 sparse Matrix of class "dgCMatrix"


##```````````````````````````````
# Properties of Sparse Matrices`
#````````````````````````````````

# NA values are not considered equivalent to sparsity 
# and therefore are treated as non-zero values
# yet they don’t participate in any sparse matrix operations.


mat_NA <- matrix(c(5.5, 0, NA,
         0, 0, NA), nrow = 3)

print(mat_NA)

sparse_matNA <- as(mat_NA, "sparseMatrix")
# 3 x 2 sparse Matrix of class "dgCMatrix"

#[1,] 5.5  .
#[2,] .    .
#[3,]  NA NA


sparse_matNA[3,] <- 0
sparse_matNA
# 3 x 2 sparse Matrix of class "dgCMatrix"

#[1,] 5.5 .
#[2,] .   .
#[3,] .   .








