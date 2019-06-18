library(gstat)
library(rgdal)

library(lattice)

library(spgwr)
library(sp)




counties <- readRDS('./counties.rds')
p<-readRDS('./precipitation.rds')

names(p)[names(p) == 'LAT'] <- 'Y'
names(p)[names(p) == 'LONG'] <- 'X'

coordinates(p) <- ~X+Y


proj4string(p) <- CRS('+proj=longlat +datum=NAD83')

countiesTA <- spTransform(counties, CRS('+proj=longlat +datum=NAD83'))

#grid for kriging

ca <- readOGR("./CA_State/CA_State_TIGER2016.shp")
ca <- spTransform(ca, CRS('+proj=longlat +datum=NAD83'))
grd <- makegrid(ca, n = 5000)
colnames(grd) <- c('X','Y')
grd_pts <- SpatialPoints(coords = grd, proj4string=CRS(proj4string(ca)))
grd_pts_in <- grd_pts[ca, ]

aq <- p

hist(aq$NOV)
spplot(aq, "NOV",  colorkey = TRUE, main = "precipitation")
#log because of skewed data
prec <- log(aq$NOV)


#model variogram for ordinry kriging
v <- variogram(prec~1, locations= coordinates(aq), data=aq)


v.fit1 <-  fit.variogram(v, vgm("Sph"),fit.kappa =TRUE)
v.fit2 <-  fit.variogram(v, vgm("Mat"),fit.kappa =TRUE)
v.fit3 <-  fit.variogram(v, vgm("Exp"),fit.kappa =TRUE)

print(attr(v.fit1,"singular"))
print(attr(v.fit1,"SSErr"))
print(attr(v.fit3,"singular"))
print(attr(v.fit2,"SSErr"))
print(attr(v.fit3,"singular"))
print(attr(v.fit3,"SSErr"))

#with trend
v2 <- variogram(prec~ X + Y, locations= coordinates(aq), data=aq)

v.fit <-  fit.variogram(v, vgm(c("Sph","Mat","Exp")),fit.kappa =TRUE)

v.fit_trend <- fit.variogram(v2, vgm(c("Sph","Mat","Exp")))

print(attr(v.fit_trend,"singular"))
print(attr(v.fit_trend,"SSErr"))


#use estimated for kriging
k <- krige(prec~1, aq, grd_pts_in, v.fit2)

k2 <- krige(prec~ (X +Y), aq, grd_pts_in, model = v.fit_trend) 

#plot results
spplot(k["var1.pred"], main = "Ordinary kriging predictions",cuts=15)
spplot(k["var1.var"], main = "Ordinary kriging variance",cuts=15)

spplot(k2["var1.pred"], main = "Universal kriging predictions", cuts=15)
spplot(k2["var1.var"], main = "Universal kriging variance", cuts=15)

#backtransforming log values

YY = exp(k$var1.pred + 0.5 * (k$var1.var))
mu_bt <- mean(YY)
mu_original <- mean(exp(prec))

print(mu_bt)
print(mu_original)

YY_corr <- YY* (mu_original/mu_bt)

##visualization of back-transformed data
coords<-data.frame(k$X,k$Y)
data   <- as.data.frame(YY_corr)
crs <- CRS('+proj=longlat +datum=NAD83')
no_log<- SpatialPointsDataFrame(coords = coords,
                               data = data, 
                               proj4string = crs)
spplot(no_log["YY_corr"], main = "Ordinary kriging non-log", cuts=15)


##evaluate kriging fit with variogramwith cross validation

set.seed(123)


a <-c(rep(1, length(prec)/5),rep(2, length(prec)/5),rep(3, length(prec)/5),rep(4, length(prec)/5),rep(5, (length(prec)/5)+1))

a<- sample(a)

RMSE <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm=TRUE))
}

MAE <- function(observed, predicted) {
  mean(abs(predicted - observed))
}

krigmae1<- krigmae2<-krigrmse1 <- krigrmse2 <- rep(NA, 5)
for (i in 1:5) {
  train <- aq[a!=i,]
  test <- aq[a==i,]
  

  m <- gstat(formula= log(train$NOV)~1, locations=train, model=v.fit2)
 
  p2 <- predict(m, newdata=test, debug.level=0)$var1.pred
  krigrmse1[i] <-  RMSE(log(test$NOV), p2)
  krigmae1[i] <-  MAE(log(test$NOV), p2)
  
  m <- gstat(formula= log(train$NOV)~ X+Y, locations=train, model=v.fit_trend)
  p2 <- predict(m, newdata=test, debug.level=0)$var1.pred
  krigrmse2[i] <-  RMSE(log(test$NOV), p2)
  krigmae2[i] <-  MAE(log(test$NOV), p2)
  
}

#with implemented function
cv1 <- krige.cv(log(NOV)~1, aq, v.fit2, nfold=5)

# mean error, ideally 0:
mean(cv1$residual)
# MSPE, ideally small
mean(cv1$residual^2)
# Mean square normalized error, ideally close to 1
mean(cv1$zscore^2)
# correlation observed and predicted, ideally 1
cor(cv1$observed, cv1$observed - cv1$residual)
# correlation predicted and residual, ideally 0
cor(cv1$observed - cv1$residual, cv1$residual)

cv2 <- krige.cv(log(NOV)~X+Y, aq, v.fit3, nfold=5)

# mean error, ideally 0:
mean(cv2$residual)
# MSPE, ideally small
mean(cv2$residual^2)
# Mean square normalized error, ideally close to 1
mean(cv2$zscore^2)
# correlation observed and predicted, ideally 1
cor(cv2$observed, cv2$observed - cv2$residual)
# correlation predicted and residual, ideally 0
cor(cv2$observed - cv2$residual, cv2$residual)

#####cokriging



v.lpb <- variogram(log(JAN)~1, aq)
v.fit_pb <- fit.variogram(v.lpb, vgm(c("Sph","Mat","Exp")))


prec.g <- gstat(id = "NOV", formula = log(NOV) ~ 1, data = aq)
prec.g <- gstat(prec.g, "DEC", log(DEC) ~ 1, aq)
prec.g <- gstat(prec.g, model = v.fit_pb, fill.all = TRUE)

x <- variogram(prec.g, cutoff = 1000)
prec.fit <- fit.lmc(x, prec.g)

out <- gstat.cv(prec.fit, nmax = 40, nfold = c(rep(1,100), rep(2,55)))
summary(out)
mean(out$residual)
mean(out$residual^2)
mean(out$zscore^2)
cor(out$observed, out$observed - out$residual)


k.c2 <- predict(prec.fit, grd_pts_in, debug.level=0)

spplot(k.c2["NOV.pred"], main = "Cokriging", cuts=15)

#geographycally weighted regression

bw <- gwr.sel(NOV ~ ALT, data=aq)
gwr.model <- gwr(NOV ~ ALT, data=aq, bandwidth=bw, hatmatrix=T,se.fit=T)

spplot(gwr.model$SDF, "ALT", cuts=quantile(gwr.model$SDF$ALT))

t <- gwr.model$SDF$ALT/gwr.model$SDF$ALT_se
sig.map = SpatialPointsDataFrame(aq, data.frame(t))

spplot(sig.map)




