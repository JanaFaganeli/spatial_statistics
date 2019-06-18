library(spatstat)
library(spatstat.data)
library(maptools)




bf<-readRDS('./big_foot.rds')


bf <- bf[,c(1,2)]

pbf<- ppp(bf[,1],bf[,2],poly)
poly <- owin(xrange=range(bf$lon), 
             yrange=range(bf$lat))


# Do Big foots live in clusters?

#where is big boot near the Great lakes?

# plot US map
data(wrld_simpl)
plot(wrld_simpl, add=TRUE)