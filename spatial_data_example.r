library(raster)
library(sp)


#read image into raster
r0 <- raster('/home/jana/Ispra_delo/2013/slovenija/izvoz_ploskovne_Slovenija/emiss2_dis_pm25.tif')

ll <- projectRaster(r0,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#ll <- r0

#slovenija
b <- as(extent(13.2, 16.7, 45.3, 47), 'SpatialPolygons')
crs(b) <- crs(r0)

rbn <- crop(ll, b)


###################points

to <- read.csv('/home/jana/Ispra_delo/2013/slovenija/emiss_2013_tockovne_tone.txt',header=TRUE)
ta <- read.csv('/home/jana/Ispra_delo/2013/slovenija/emiss_2013_tockovne_tone.txt',header=TRUE)
coordinates(ta) <- ~lon + lat 

#to make a Spatial Points object, used for prediction points in gstat
ngrd <- data.frame(to$lon,to$lat)
new_pts <- SpatialPoints(coords = ngrd, proj4string=CRS('+proj=longlat +datum=NAD83'))

#to make a raster
pts<- cbind(to$lon,to$lat)
rt0 <- raster()

extent(rt0)<-extent(r0)
crs(rt0) <- crs(r0)

############## joining different data in a single raster

rb0 <- raster()


a<- c(1,2,3,4,5,6,7,8,9,10)
par <-data.frame()

for (i in a){
  if(i %in% c(1,3,6,9)){
    
    
    rast2 <- rasterize(pts, rbn, to[to$snap==i,]$pm25, fun=sum) 
    rb0 <- stack(rb0,rast2,bands=i)
    
    par <- rbind(par,data.frame(sec=i,part=sum(rast2[],na.rm=TRUE)))
  }
  
  else{
    r <- raster(paste('/home/jana/Ispra_delo/2013/slovenija/izvoz_ploskovne_Slovenija/emiss',i,'_dis_pm25.tif',sep=''))
    
    ll <- projectRaster(r,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",method="ngb")

    rb <- crop(ll, b)
 
    rb0 <- stack(rb0,rb,bands=i)
    par <- rbind(par,data.frame(sec=i,part=sum(rb[],na.rm=TRUE)))
  }
}

names(rb0) <- c('Sector 1','Sector 2','Sector 3','Sector 4','Sector 5','Sector 6','Sector 7','Sector 8','Sector 9','Sector 10')

#read shape file

#rsh <-readOGR('/home/jana/delavnica/Rskripte/Mreza_5km_x_5km/Mreza_5km.shp')


