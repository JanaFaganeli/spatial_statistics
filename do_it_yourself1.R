library(gstat)
library(rgdal)

library(lattice)

library(sp)
library(rspatial)


##air quality California data

p <- readRDS('./california_ozone.rds')



names(p)[names(p) == 'LATITUDE'] <- 'Y'
names(p)[names(p) == 'LONGITUDE'] <- 'X'

coordinates(p) <- ~X+Y
proj4string(p) <- CRS('+proj=longlat +datum=NAD83')

#find the best estimates at the points
#LATITUDE 2.8364 34.1072 37.3483 37.6825 40.7769
#LONGITUDE -117.1289 -116.7238 -121.8947 -121.4405 -124.1775


######grid from my example can be used


