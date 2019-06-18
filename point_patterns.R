library(spatstat)
library(spatstat.data)
library(maptools)



cc<- cells
summary(cells)
plot(density(cells))


M <- quadrat.test(cc, nx=3, ny=3)
plot(cc,main='Cell data-quadrat count')
plot(M, add=TRUE )


Y <- redwood
plot(Y)
data(Y)
summary(Y)
plot(density(Y))



DD <- lansing
data(lansing)
plot(DD)
plot(split(DD))
yy <- lansing$y
xx <- lansing$x


plot(split(DD))

dd <- split(DD)$whiteoak


### G frunction
par(mfrow=c(3,2))
plot(Fest(split(DD)$blackoak))
plot(Fest(split(DD)$hickory))
plot(Fest(split(DD)$maple))
plot(Fest(split(DD)$misc))
plot(Fest(split(DD)$redoak))
plot(Fest(split(DD)$whiteoak))



par(mfrow=c(3,2))
plot(Fest(split(DD)$blackoak),. ~theo)
plot(Fest(split(DD)$hickory),. ~theo)
plot(Fest(split(DD)$maple),. ~theo)
plot(Fest(split(DD)$misc), .~theo)
plot(Fest(split(DD)$redoak),. ~theo)
plot(Fest(split(DD)$whiteoak),. ~theo)

par(mfrow=c(1,1))
plot(Fest(cc))


############G function
par(mfrow=c(3,2))
plot(Gest(split(DD)$blackoak))
plot(Gest(split(DD)$hickory))
plot(Gest(split(DD)$maple))
plot(Gest(split(DD)$misc))
plot(Gest(split(DD)$redoak))
plot(Gest(split(DD)$whiteoak))

par(mfrow=c(3,2))
plot(Gest(split(DD)$blackoak),. ~theo)
plot(Gest(split(DD)$hickory),. ~theo)
plot(Gest(split(DD)$maple),. ~theo)
plot(Gest(split(DD)$misc),. ~theo)
plot(Gest(split(DD)$redoak),. ~theo)
plot(Gest(split(DD)$whiteoak),. ~theo)

par(mfrow=c(1,1))
plot(Gest(cc))

###K function

par(mfrow=c(3,2))
plot(Kest(split(DD)$blackoak))
plot(Kest(split(DD)$hickory))
plot(Kest(split(DD)$maple))
plot(Kest(split(DD)$misc))
plot(Kest(split(DD)$redoak))
plot(Kest(split(DD)$whiteoak))



par(mfrow=c(3,2))
plot(Lest(split(DD)$blackoak))
plot(Lest(split(DD)$hickory))
plot(Lest(split(DD)$maple))
plot(Lest(split(DD)$misc))
plot(Lest(split(DD)$redoak))
plot(Lest(split(DD)$whiteoak))

par(mfrow=c(1,1))
plot(pcf(split(DD)$blackoak))

#######envelopes
par(mfrow=c(3,2))
plot(envelope(split(DD)$blackoak, Kest, nsim=100, rank=1))
plot(envelope(split(DD)$hickory, Kest, nsim=100, rank=1))
plot(envelope(split(DD)$maple, Kest, nsim=100, rank=1))
plot(envelope(split(DD)$misc, Kest, nsim=100, rank=1))
plot(envelope(split(DD)$redoak, Kest, nsim=100, rank=1))
plot(envelope(split(DD)$whiteoak, Kest, nsim=100, rank=1))

par(mfrow=c(1,1))
plot(envelope(split(DD)$whiteoak, Kest, nsim=100, rank=1))


par(mfrow=c(1,1))
plot(envelope(split(DD)$whiteoak, Lest, nsim=100, rank=1, global=TRUE))

par(mfrow=c(1,1))
plot(envelope(cc, Kest, nsim=100, rank=1))

par(mfrow=c(1,1))
plot(envelope(cc, Lest, nsim=100, rank=1, global=TRUE))


############kernel density estimation and nearest neighbor distance

#cells
plot(density(cells, sigma=bw.diggle(cells), edge= TRUE, diggle=TRUE))
plot(nndensity(cells))
#trees
plot(density(split(DD)$maple, sigma=bw.diggle(split(DD)$maple), edge= TRUE, diggle=TRUE, SE=TRUE))
plot(nndensity(split(DD)$maple))



##likelihood maximization of an inhomogenious Poisson process
ppm_es <- ppm(cells, ~1)
ppm_es1
plot(ppm_es1)


ppm_es1 <- ppm(split(DD)$maple, trend= ~x+y)
ppm_es1
plot(ppm_es1)

ppm_es2 <- ppm(split(DD)$maple, trend= ~polynom(x,y,2))
ppm_es2
plot(ppm_es2)

#with additional covariates
data(bei)
el <- bei.extra$elev
ppm_es3 <- ppm(bei, ~el ,covariates= list(el))

#checking the Poisson model
diagnose.ppm(ppm_es)
diagnose.ppm(ppm_es1)
diagnose.ppm(ppm_es2)
diagnose.ppm(ppm_es3)


qqplot.ppm(ppm_es4,nsim=50)




