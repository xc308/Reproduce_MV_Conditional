## Do some variogram
####################
library(sp)
library(gstat)

data(meuse)
class(meuse)  # "data.frame"

names(meuse)
coordinates(meuse) = ~x+y
class(meuse)  # "SpatialPointsDataFrame"
# attr(,"package")
# [1] "sp"

summary(meuse)

coordinates(meuse)[1:5]
# promotes the data.frame meuse into a SpatialPoints- DataFrame, 
# which knows about its spatial coordinates;


################################3
# Spatial data on a regular grid
################################3

data("meuse.grid")
summary(meuse.grid)
class(meuse.grid)  # "data.frame"
coordinates(meuse.grid) = ~x+y
class(meuse.grid)  # "SpatialPointsDataFrame"
                    #attr(,"package")  [1] "sp"
gridded(meuse.grid) = TRUE
summary(meuse.grid)

class(meuse.grid)  # "SpatialPixelsDataFrame"

image(meuse.grid["dist"])
title("distance to river (red = 0)")



#############
## Variograms
#############
# first argument: log(zinc)~1 means that we assume a constant trend for the variable log(zinc)
lzn.vgm = variogram(log(zinc)~1, meuse)
lzn.vgm

lzn.fit = fit.variogram(lzn.vgm, model = vgm(1, "Sph", 900, 1))
lzn.fit

lzn.fit2 = fit.variogram(lzn.vgm, model = vgm(1, "Exp", 1, 1))
lzn.fit2
plot(lzn.vgm, lzn.fit)



#######################
## try on weather data
#######################

head(weather)
names(weather)
coordinates(weather) = ~lon+lat

class(weather)  # [1] "SpatialPointsDataFrame"
p.vgm <- variogram(pressure ~ temperature, weather)
p.vgm
p.fit <- fit.variogram(p.vgm, model = vgm(psill = 1, "Exp", 1, 1))
#(still not work)
p.fit <- autofitVariogram(pressure ~ temperature,  weather, model = "Exp")
plot(p.vgm, p.fit)


##==========
# Error
##===========
#No convergence after 200 iterations: try different initial values?
# need to regress on some covariates, not just 1
#p.fit <- fit.variogram(p.vgm, model = vgm(1, "Exp", 800, 1))


install.packages("automap")
library(automap)
fit.vgm = autofitVariogram(pressure ~1, weather, model = "Exp")
plot(p.vgm, fit.vgm$var_model)


library(bicon)
data(temps)

head(temps)
plot(temps$lon,temps$lat)


Y2 <- left_join(temps,X) %>%
  group_by(DATE) %>%
  mutate(TMAX_anom = TMAX/10 - mean(TMAX/10),
         TMIN_anom = TMIN/10 - mean(TMIN/10)) %>% 
  data.frame()



Y2 <- temps %>%
  mutate(TMAX_anom = maxT/10 - mean(maxT/10),
         TMIN_anom = minT/10 - mean(minT/10)) %>%
  data.frame()


coordinates(Y2)<- ~lon+lat
class(Y2)
head(Y2)

minT.vgm <- variogram(TMIN_anom ~ 1, Y2)
minT.fit <- fit.variogram(minT.vgm, model = vgm(psill = 1, model = "Exp", range = 1, nugget = 1))



coordinates(Ysp) <- ~lon + lat
maxT.vgm = variogram(maxT~ELEVATION, Ysp)
maxT.fit = fit.variogram(maxT.vgm, model = vgm(1, "Exp", 1, 1))
plot(maxT.vgm, maxT.fit)





