#=============#
# tidyr:gather
#=============#

# ref:https://tidyr.tidyverse.org/reference/gather.html



stocks <- tibble(time = as.Date("2021-10-15") + 0:9,
                 X = rnorm(10, 0, 1),
                 Y = rnorm(10, 0, 2),
                 Z = rnorm(10, 0, 4))

head(stocks)


gather(stocks, key = "stock", value = "price", -time)
gather(stocks, key = "stock", value = "price", -time, convert = T)



head(iris)
# get the 1st obs for each specise
mini_iris <- iris[c(1, 51, 101), ]
head(mini_iris)


gather_a <- gather(mini_iris, key = "attr", value = "measurement", 
       -Species)


gather_b <- gather(mini_iris, key = "attr", value = "measure",
                   Sepal.Length, Sepal.Width, Petal.Length, Petal.Width)

gather_b

# two are the same 



# using pipe operator
mini_iris2 <- iris %>% group_by(Species) %>% slice(1)
mini_iris2
head(mini_iris2) # the same as mini_iris


mini_iris2 %>% gather(key = "attr", value = "measure",
                      -Species)




















