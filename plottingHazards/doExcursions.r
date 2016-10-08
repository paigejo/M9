setwd("~/git/M9")
source("loadTestData.r")
source("fitExpVarioParams.RData")
library(fields)

#we now have s, n, r for Cov model:
#C(d | s,n,r) = (s - n)*exp(-d/r), d>0

covFun = function(d) {
  (s - n)*exp(-d/r)
}

R = 3959 #radius of the Earth in miles
lonExtent= c(235.79781, 235.82087)
latExtent = c(41.739671,41.762726)
CCLon = mean(lonExtent)
CCLat = mean(latExtent)
lonDist = cos(CCLat)*2*pi/360*R
latDist = 2*pi*R/360

X = lon*lonDist + min(lon)
Y = lat*latDist + min(lat)
coords = cbind(X, Y)

#get mean flood field
mu.x = apply(allHMax, 1, mean)

#clear memory
rm(allHMax, X, Y, VG, coVG)

# get precision matrix
mat = rdist(coords, coords) #get distance mat
mat = covFun(mat) #get covariance mat
mat = solve(mat) #get precision matrix

#do excursions for "wet tennis shoes" flood level
out.0 = excursions(alpha=.1, u=0.01, mu=mu.x, type="=", Q=mat, verbose=1)
par(mfrow=c(2,1))
nX = dim(allHMax)[2]
nY = dim(allHMax)[3]
quilt.plot(X, Y, out.0$F, nx=nX, ny=nY, main="Excursion Function")
quilt.plot(X, Y, mu.x, nx=nX, ny=nY, main="Mean Field")

