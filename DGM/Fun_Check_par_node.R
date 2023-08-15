#=====================
# Checking parent node  (15 Aug 2023)
#=====================

#rm(list = ls())


# Aim: 
  # learn how to select parent nodes for my automation algo

# Steps:
  # 1. Construct the Hierarchical Structure:
      # a data frame with columns like "node_id", "par_id".
  
  # 2. Define a Function:
      # input: desired query input node_id; data with hierarchical strucuture
      # output: the parent node id corresponds to the input one


#------------------------
# hierarchy data strutcture
#------------------------

hierarchy_data <- data.frame(
  node_id = c(1, 2, 3, 4),
  par_id = c(NA, 1, 2, 1)
  #ch_id = list(list(2, 4), list(3), list(), list())
)

head(hierarchy_data)



#----------------------------------------------------------
# function to check if one node is a parent node of another
#----------------------------------------------------------

# input: query node_id; data with hierarchy structure
# output: parent node id for the query node


Check_par_node <- function(node_id, data = hierarchy_data) {
  
  if (!all(colnames(data) == c("node_id", "par_id"))) {
    stop("Data must have two columns with the first col name 'node_id' and the second name 'par_id'")
  }
  
  # input node_id, return its parent id
  par_index <- data[node_id, ]$par_id
  
  if (is.na(par_index)) {
    return("The input node has no parent.")
  } else {
    return(par_index)
  }
  
}


#-----
# Test
#------
Check_par_node(node_id = 3, data = hierarchy_data) # [1] 2
Check_par_node(node_id = 4, data = hierarchy_data) # [1] 1



