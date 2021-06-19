#==============#
# cbind & cBind
#==============#

library(Matrix)

(a <- matrix(c(2:1,1:2), 2,2))
#cbind(0, rBind(a, 7)) # remains traditional matrix
cbind(0, rbind(a, 7))

D <- Diagonal(2)
#cBind(4, a, D, -1, D, 0) # a sparse Matrix
cbind(4, a, D, -1, D, 0) # a sparse Matrix