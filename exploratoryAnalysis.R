# This file loads all the data and functions for the analysis and sets the working directory

library(fields)
setwd("~/git/M9/")
source('plotSubfault.R')
source('loadTestData.r')
source('fitModel.R')

##### plot slip rates inferred from GPS data
slipDat = read.table("Cascadia-lockingrate-hdr.txt", header=TRUE)
names(slipDat) = c("lon", "lat", "slip", "Depth", "slipErr")

attach(slipDat)

quilt.plot(lon, lat, slip)
map("world", "Canada", add=TRUE)
US(add=TRUE)
quilt.plot(lon, lat, Depth)
map("world", "Canada", add=TRUE)
world(add=TRUE)
quilt.plot(lon, lat, slipErr)
map("world", "Canada", add=TRUE)
world(add=TRUE)

faultGeom = read.csv("CSZe01.csv")
faultGeom$longitude = faultGeom$longitude - 360
kmCols = c(4, 6, 7)
faultGeom[,kmCols] = faultGeom[,kmCols] * 10^3

plot(lon, lat, type="n")
quilt.plot(lon, lat, rep(NA, length=length(Depth)), zlim=c(range(Depth)))
#plotSubfault(faultGeom[1,], varRange=range(Depth))
#apply(faultGeom, 1, plotSubfault, varRange=range(Depth))
#plotSubfaultGeom(faultGeom[1,], varRange=range(Depth))
plotFault(faultGeom, new=FALSE, plotData=FALSE)
world(add=TRUE)

quilt.plot(lon, lat, Depth, nx=80, ny=80)
world(add=TRUE, lwd=2)
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)

quilt.plot(lon, lat, slip, nx=120, ny=120, main="Locking Rate (mm/yr)", ylab="Latitude", xlab="Longitude")
world(add=TRUE, lwd=2)
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)
quilt.plot(lon, lat, Depth, nx=120, ny=120, main="Depth (km)", ylab="Latitude", xlab="Longitude")
world(add=TRUE, lwd=2)
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)
quilt.plot(lon, lat, slipErr, nx=120, ny=120, main="Locking Rate SE (mm/yr)", ylab="Latitude", xlab="Longitude")
world(add=TRUE, lwd=2)
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)

##### plot subsidence inferred from turbidite data
load("DR1.RData")
dr1$bpSE = dr1$bpMOE/qnorm(.975)
dr1$rcybpSE = dr1$rcybpMOE/qnorm(.975)

attach(dr1)

plot(dr1$Lon, dr1$Lat, col="red", pch="+")
map("world", "Canada", add=TRUE)
US(add=TRUE)

# plot histogram of different estimated ages for events (everything's so muddled...)
hist(c(dr1$bpCntr, dr1$rcybpCntr) + 66, breaks=200)

events = as.character(event)
uniqueEvents = unique(events)
sortI = rev(c(1, 11, 10, 13, 2, 14, 3, 4, 17, 5, 12, 6, 18, 7, 19, 9, 20, 8, 15, 16)) # T1 is most recent, so reverse
uniqueEvents = uniqueEvents[sortI]

latRange = range(Lat)
subRange = range(range(subsidence - Uncertainty), range(subsidence + Uncertainty))
par(mar=c(4, 4.5, 1.093333, 0.50000))
for(i in 1:length(uniqueEvents)) {
  # get event data for the given event from table
  e = uniqueEvents[i]
  eventDat = dr1[events == e,]
  print(paste0("Amount of data for event ", e, ": ", nrow(eventDat)))
  
  # plot the data (reverse order for consistency with Leonard et al. 2010 paper) with error bars
  plot(eventDat$Lat, eventDat$subsidence, main=e, xlab="Latitude", ylab="Subsidence (m)", 
       xlim=rev(latRange), ylim=rev(subRange), type = "n")
  arrows(eventDat$Lat, eventDat$subsidence, eventDat$Lat, 
         eventDat$subsidence + eventDat$Uncertainty, 
         angle=90, col="grey")
  arrows(eventDat$Lat, eventDat$subsidence, eventDat$Lat, 
         eventDat$subsidence - eventDat$Uncertainty, 
         angle=90, col="grey")
  points(eventDat$Lat, eventDat$subsidence, pch=19, cex=.7)
}
plot(rev(dr1$Lat), rev(dr1$subsidence), main="All events", 
     xlab="Latitude", ylab="Subsidence (m)", xlim=rev(latRange), ylim=rev(subRange), type="n")
arrows(dr1$Lat, dr1$subsidence, dr1$Lat, 
       dr1$subsidence + dr1$Uncertainty, 
       angle=90, col="grey")
arrows(dr1$Lat, dr1$subsidence, dr1$Lat, 
       dr1$subsidence - dr1$Uncertainty, 
       angle=90, col="grey")
points(rev(dr1$Lat), rev(dr1$subsidence), pch=19, cex=.7)


##### plot subsidence inferred from turbidite data with slip rates from GPS data

# first show where all the data is
plot(lon, lat, pch=19, cex=.1, col="blue")
world(add=TRUE)
points(dr1$Lon, dr1$Lat, pch=3, cex=.5, col="red")

# now compute closest slip data point for each subsidence data point
distMat = rdist(dr1[,c(4, 3)], slipDat[,1:2])
slipI = apply(distMat, 1, which.min)
slipDatPaired = slipDat[slipI,]

# check to make sure it makes sense
points(slipDatPaired[,1:2])
lonRange = range(slipDatPaired$lon)
latRange = range(slipDatPaired$lat)
slipRange = range(slipDatPaired$slip)
quilt.plot(lon, lat, slip, xlim=lonRange, ylim=latRange, zlim=slipRange)
world(add=TRUE)
quilt.plot(slipDatPaired[,1:2], slipDatPaired$slip, xlim=lonRange, ylim=latRange, zlim=slipRange)
world(add=TRUE)

# now plot subsidence versus slip rate
# NOTE: This looks terrible
plot(Lat, subsidence, pch=19, cex=.5)
points(Lat, slipDatPaired$slip/15, col="blue", cex=.7)

plot(Lat, subsidence[event == "T1"], pch=19, cex=.5)
points(Lat, slipDatPaired$slip/15, col="blue", cex=.7)

plot(Lat[event != "T1"], subsidence[event != "T1"], pch=19, cex=.5)
points(Lat, slipDatPaired$slip/15, col="blue", cex=.7)

# plot slip rate with subsidence as a second plot on the side
par(mai=c(1.360000, 1.093333, 1.093333, 0.560000))
par(mar=c(4, 2.5, 1.093333, 0.50000))

# layout(matrix(c(1, 2), 1, 2), widths=c(3, 1), heights=2)
par(mfrow=c(1,2), mar=c(4, 3.5, 2, 0.50000))
quilt.plot(lon, lat, slip, xlim=lonRange, ylim=latRange, nx = 120, ny=120, 
           main="Slip rate (mm/yr)", xlab="Longitude", ylab="Latitude")
world(add=TRUE)
points(dr1$Lon, dr1$Lat, pch=3, cex=.7)
plot(dr1$subsidence, dr1$Lat, pch=19, cex=.3, ylim=latRange, 
     main="All Subsidence", ylab="", xlab="Subsidence (m)")

quilt.plot(lon, lat, slip, xlim=lonRange, ylim=latRange, nx = 120, ny=120, 
           main="Slip rate (mm/yr)", xlab="Longitude", ylab="Latitude")
world(add=TRUE)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], pch=3, cex=.7)
plot(dr1$subsidence[event == "T1"], dr1$Lat[event == "T1"], pch=19, cex=.3, ylim=latRange, 
     main="T1 Subsidence", ylab="", xlab="Subsidence (m)")

# final plots including fault geometry
lonRange = range(lon)
latRange = range(lat)
slipRange = range(slip)

pdf(file="lockRateVsSub.pdf", width=4, height=5)
quilt.plot(lon, lat, slip, xlim=lonRange, ylim=latRange, nx = 120, ny=120, 
           main="Locking rate (mm/yr)", xlab="Longitude", ylab="Latitude")
world(add=TRUE, lwd=2)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)
plot(dr1$subsidence[event == "T1"], dr1$Lat[event == "T1"], pch=19, cex=.3, ylim=latRange, 
     main="T1 Subsidence", ylab="", xlab="Subsidence (m)")
dev.off()

par(mfrow=c(1,3))
pdf(file="lockRateVsSubFull.pdf", width=4, height=5)
quilt.plot(lon, lat, slip, xlim=lonRange, ylim=latRange, nx = 120, ny=120, 
           main="Locking rate (mm/yr)", xlab="Longitude", ylab="Latitude")
world(add=TRUE, lwd=2)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)
plot(dr1$subsidence[event == "T1"], dr1$Lat[event == "T1"], pch=19, cex=.3, ylim=latRange, 
     main="T1 Subsidence", ylab="", xlab="Subsidence (m)")
plot(slipDatPaired$slip/15, slipDatPaired$lat, pch=19, cex=.3, ylim=latRange, 
     main="GPS lock rate / 15", ylab="", xlab="Subsidence (m)")
dev.off()

pdf(file="depthVsSub.pdf", width=4, height=5)
quilt.plot(lon, lat, Depth, xlim=lonRange, ylim=latRange, nx = 120, ny=120, 
           main="Depth (km)", xlab="Longitude", ylab="Latitude")
world(add=TRUE, lwd=2)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)
plot(dr1$subsidence[event == "T1"], dr1$Lat[event == "T1"], pch=19, cex=.3, ylim=latRange, 
     main="T1 Subsidence", ylab="", xlab="Subsidence (m)")
dev.off()

pdf(file="slipErrVsSub.pdf", width=4, height=5)
quilt.plot(lon, lat, slipErr, xlim=lonRange, ylim=latRange, nx = 120, ny=120, 
           main="Locking rate SE (mm/yr)", xlab="Longitude", ylab="Latitude")
world(add=TRUE, lwd=2)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)
plot(dr1$subsidence[event == "T1"], dr1$Lat[event == "T1"], pch=19, cex=.3, ylim=latRange, 
     main="T1 Subsidence", ylab="", xlab="Subsidence (m)")
dev.off()

par(mfrow=c(1,1))

##### test Okada
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
tmp = okadaSubfault(faultGeom[1,], lonGrid, latGrid, slip=10)

quilt.plot(lon, lat, rep(NA, length(lon)), xlim=lonRange, ylim=latRange, 
           main="Seafloor deformation (m)", xlab="Longitude", ylab="Latitude")
plotSubfault(faultGeom[1,], plotData=FALSE)
quilt.plot()
testRow = faultGeom[1,]
testRow[1,4] = testRow[1,4] * 10^3

# set up the test case from okada.ipynb
nx = 200
ny=  200
lonGrid = seq(-.5, 1.5, l=nx)
latGrid = seq(-1, 1, l=ny)
testRow = data.frame(rbind(c(Fault=30.1, longitude=0, latitude=0, depth=5e3, 
                             strike=0, length=100e3, width=50e3, dip=10)))
dZ = okadaSubfault(testRow, lonGrid, latGrid, slip=1)
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(dZ)))
# It works!!!

# now test across entire CSZ geometry
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)

testFault = faultGeom

dZ = okada(testFault, lonGrid, latGrid, slips=10)
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(dZ)), nx=120, ny=120)
world(add=TRUE, lwd=2)

##### now plot versus subsidence data
# calculate the simulated subsidence at the data locations
# round the data locations to the lon lat grid
roundToRange = function(coordVec, coordRange) {
  inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
  inds = round(inds*(length(coordRange)-1)) + 1
  roundedCoordVec = coordRange[inds]
  return(roundedCoordVec)
}
roundLon = roundToRange(dr1$Lon, lonGrid)
roundLat = roundToRange(dr1$Lat, latGrid)
roundCoords = cbind(roundLon, roundLat)

# find indices of grid coords corresponding to rounded data coords
findIndex = function(rCoords, gCoords) {
  return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
}
coordGrid = make.surface.grid(list(lonGrid, latGrid))
inds = findIndex(roundCoords, coordGrid)

simDef = c(t(dZ))[inds]

par(mfrow=c(1,3))

pdf(file="slip10Subsidence.pdf", width=4, height=5)
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(dZ)), nx=120, ny=120)
world(add=TRUE, lwd=2)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
apply(faultGeom, 1, plotData=FALSE, new=FALSE, lwd=2)
plot(dr1$subsidence[event == "T1"], dr1$Lat[event == "T1"], pch=19, cex=.3, ylim=latRange, 
     main="T1 Subsidence", ylab="", xlab="Subsidence (m)", xlim=range(-simDef))
plot(-simDef, coordGrid[inds,2], pch=19, cex=.3, ylim=latRange, 
     main="Simulated Subsidence", ylab="", xlab="Subsidence (m)")
dev.off()

par(mfrow=c(1,1))

##### test for subdivided faults
plotFault(testFault, xlim=c(-128, -123.75), ylim=c(40, 50))
world(add=TRUE)
test = divideSubfault(testRow)
plotFault(test, xlim=c(0, .45), ylim=c(-.45, .45))
plotSubfault(testRow, xlim=c(0, .45), ylim=c(-.45, .45), new=TRUE)
points(test$longitude, test$latitude)

test = divideFault(testFault)
plotFault(testFault, xlim=c(-128, -123.75), ylim=c(40, 50))
plotFault(test, xlim=c(-128, -123.75), ylim=c(40, 50))

# make sure the subdivided depths are correct
quilt.plot(test$longitude, test$latitude, test$depth, xlim=c(-128, -123.75), ylim=c(40, 50), zlim=c(0, 20000))
quilt.plot(testFault$longitude, testFault$latitude, testFault$depth, xlim=c(-128, -123.75), ylim=c(40, 50), zlim=c(0, 20000))
#It works!

# test okada on subfaults
csz = divideFault(testFault)
test = okada(csz, lonGrid, latGrid, slips=10)
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(test)), nx=120, ny=120)
world(add=TRUE, lwd=2)

##### Assume the slip is proportional to the GPS data.  Plot the Okada results
distMat = rdist(cbind(csz$longitude, csz$latitude), slipDat[,1:2])
slipI = apply(distMat, 1, which.min)
slipDatCSZ = slipDat[slipI,]
lockRates = slipDatCSZ$slip

# plot them to make sure they're correct (NOTE: the large number of bottom red locations 
# are a consequence of the nearest neighbor prediction. Hopefully that won't cause problems)
quilt.plot(csz$longitude, csz$latitude, lockRates, xlim=c(-128, -123.75), ylim=c(40, 50))
world(add=TRUE, lwd=2)
quilt.plot(lon, lat, slip, nx=120, ny=120, main="Locking Rate (mm/yr)", ylab="Latitude", xlab="Longitude")
world(add=TRUE, lwd=2)
plotFault(faultGeom, plotData=FALSE, new=FALSE, lwd=2)

gpsTest = okada(csz, lonGrid, latGrid, slips=lockRates)
simDef = c(t(gpsTest))[inds]

par(mfrow=c(2,2))
pdf(file="lockRateOkadaFull.pdf", width=4, height=5)
# quilt.plot(lon, lat, slip, xlim=lonRange, ylim=latRange, nx = 120, ny=120, 
#            main="Locking rate (mm/yr)", xlab="Longitude", ylab="Latitude")
# world(add=TRUE, lwd=2)
# plotFault(testFault, cols=c(NA, NA), new=FALSE, lwd=2)

# slip using GPS locking rates
propFac = 1/2.25
quilt.plot(csz$longitude, csz$latitude, lockRates*propFac, xlim=lonRange, ylim=latRange, nx = 30, ny=30, 
           main="Slip (m)", xlab="Longitude", ylab="Latitude")
map("world", "Canada", add=TRUE)
US(add=TRUE)
plotFault(testFault, new=FALSE, lwd=2, plotData=FALSE)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
# seadef from Okada model generated from the slips
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(gpsTest)), nx=120, ny=120, main="Seafloor Deformation (m)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
plotFault(testFault, new=FALSE, lwd=2, plotData=FALSE)
# T1 subsidence data
plot(dr1$subsidence[event == "T1"], dr1$Lat[event == "T1"], pch=19, cex=.3, ylim=latRange, 
     xlim=range(-simDef*propFac), main="T1 Subsidence", ylab="", xlab="Subsidence (m)")
# simulated subsidence data from Okada model using GPS locking rates
plot(-simDef*propFac, coordGrid[inds,2], pch=19, cex=.3, ylim=latRange, 
     main="Simulated Subsidence", ylab="", xlab="Subsidence (m)")
dev.off()

par(mfrow=c(1,1))

##### profile okada functions to figure out how to speed them up.  This will be helpful
##### when recomputing the okada output for different Poisson ratios and/or rakes
library(profvis)
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
profvis({dZ1 <- okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slips=1)})
profvis({dZ2 <- okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slips=1, rakes=100)})
system.time(dZ1 <- okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slips=1))
system.time(dZ2 <- okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slips=1, rakes=100))

# also make sure the results are right
# here we see small changes in rake have little affect on subsidence so we can probably 
# ignore rake.  Biggest difference is at near 0 subsidence in Cali.
totRange = range(c(colSums(dZ1), colSums(dZ2)))
quilt.plot(Lon, Lat, colSums(dZ1), zlim=totRange)
quilt.plot(Lon, Lat, colSums(dZ2), zlim=totRange)

# 128 seconds -> 10 seconds

##### test affect of poisson ratio on results
dZ1 <- okada(csz, lonGrid, latGrid, slips=1, poisson=.25)
dZ2 <- okada(csz, lonGrid, latGrid, slips=1, poisson=.1)
dZ3 <- okada(csz, lonGrid, latGrid, slips=1, poisson=1)
totRange = range(dZ1, dZ2, dZ3)
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(dZ1)), nx=120, ny=120, 
         main="Seafloor Deformation (m)", zlim=totRange)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(dZ2)), nx=120, ny=120, 
         main="Seafloor Deformation (m)", zlim=totRange)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(dZ3)), nx=120, ny=120, 
         main="Seafloor Deformation (m)", zlim=totRange)
points(dr1$Lon[event == "T1"], dr1$Lat[event == "T1"], cex=1, pch=3, col="red")
# higher poisson ratio tends to have more uplift.  Effect seems to be very nonlinear, though

##### Test differences between CSZ fault geometry that Randy uses to Pollitz GPS geometry

# subset GPS data so we only have data within boundaries of CSZ geometry
getGPSSubfaultDepths = function() {
  
  # helper function for determining if GPS data is within a specific subfault geometry
  getSubfaultGPSDat = function(i) {
    row = csz[i,]
    geom = calcGeom(row)
    corners = geom$corners[,1:2]
    in.poly(cbind(slipDat$lon, slipDat$lat), corners)
  }
  
  # construct logical matrix, where each column is the result of getSubfaultGPSDat(j)
  # Hence, row represents data index, column represents subfault index.
  inSubfaults = sapply(1:nrow(csz), getSubfaultGPSDat)
  
  # get the depth of each subfault by taking mean
  subfaultDepth = function(inds) {
    subfaultDat = slipDat[inds,]
    mean(subfaultDat$Depth)
  }
  apply(inSubfaults, 1, subfaultDepth)
}

par(mfrow=c(1, 2))
origDepths = getFaultCenters(csz)[,3]
depthsGPS = getGPSSubfaultDepths()
slipDatCSZ = getFaultGPSDat()
depthRange = range(c(depthsGPS, origDepths, slipDatCSZ$Depth), na.rm=TRUE)
plotFault(csz, plotVar=origDepths, main="Original Geometry", varRange=depthRange)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(csz, plotVar=depthsGPS, main="GPS Geometry", varRange=depthRange)

xRange = range(csz$longitude)
yRange = range(csz$latitude)
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$Depth, 
           main="GPS Data", xlab="Langitude", ylab="Latitude", 
           zlim=depthRange, nx=80, ny=160)
plotFault(csz, plotData = FALSE, new = FALSE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(csz, plotVar=depthsGPS, main="GPS Geometry", 
          xlim=xRange, ylim=yRange, varRange=depthRange)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)



##### Functions testing spline basis
library(splines)

# generate spline basis functions in subsidence space
latRange = c(40, 50)
lats = seq(40, 50, l=50)
Xi = bs(lats, knots=5, intercept=TRUE, Boundary.knots=latRange)

# test Xi
par(mfrow=c(1,1))
matplot(lats, Xi, main="Spline Basis", xlab="Latitude")

latRange = c(40, 50)
nKnots=2
knots = seq(40, 50, l=nKnots+2)[2:(nKnots+1)]
# Xi = ns(dr1$Lat, knots=knots, intercept=TRUE, Boundary.knots=latRange)
Xi = ns(lats, knots=knots, intercept=FALSE, Boundary.knots=latRange)

# test Xi
par(mfrow=c(1,1))
matplot(lats, Xi, main="Spline Basis", xlab="Latitude")

latRange = c(40, 50)
nKnots=3
knots = seq(40, 50, l=nKnots+2)[2:(nKnots+1)]
Xi = ns(dr1$Lat, knots=knots, intercept=FALSE, Boundary.knots=latRange)
Xi = cbind(rep(1, nrow(Xi)), Xi)
mod = lm(dr1$subsidence ~ Xi - 1)
Xi = ns(lats, knots=knots, intercept=FALSE, Boundary.knots=latRange)
Xi = cbind(rep(1, nrow(Xi)), Xi)
preds = Xi %*% coef(mod)
plot(dr1$Lat, dr1$subsidence)
lines(lats, preds)

knots = seq(40, 50, l=nKnots+2)[2:(nKnots+1)]
Xi = ns(dr1$Lat, knots=knots, intercept=TRUE, Boundary.knots=latRange)
mod = lm(dr1$subsidence ~ Xi - 1)
Xi = ns(lats, knots=knots, intercept=TRUE, Boundary.knots=latRange)
preds = Xi %*% coef(mod)
plot(dr1$Lat, dr1$subsidence)
lines(lats, preds)

df=4
Xi = bs(dr1$Lat, df=df, intercept=FALSE, Boundary.knots=latRange)
mod = lm(dr1$subsidence ~ Xi)
Xi = bs(lats, df=df, intercept=FALSE, Boundary.knots=latRange)
preds = coef(mod)[1] + Xi %*% coef(mod)[2:length(coef(mod))]
plot(dr1$Lat, dr1$subsidence)
lines(lats, preds)
matplot(Xi)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(csz$depth, lambda=fixedFitMVN$lambdaMLE)

# embed spline basis on CSZ fault plane
GQR = qr(G)
Q1 = qr.Q(GQR)
R1 = qr.R(GQR)
R1Inv = backsolve(r = R1, x = diag(ncol(R1)))
embedMat = diag(tvec^(-1)) %*% R1Inv %*% t(Q1)
XiPrime = embedMat %*% Xi
hist(log10(abs(XiPrime)), breaks=100)

par(mfrow=c(3,2))
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="Observed Locking Rate (mm/yr)", 
           xlab="Longitude", ylab="Latitude", xlim=lonRange, ylim=latRange)
plotFault(csz, new=FALSE, plotData=FALSE)
plotFault(csz, plotVar=XiPrime[,5])
for(i in 1:nrow(Xi)) {
  plotFault(csz, plotVar=XiPrime[,i], main=paste0("Spline ", i), xlab="Longitude", ylab="Latitude", 
            xlim=lonRange, ylim=latRange)
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}

# reconstruct Xi
XiReconstruct = G %*% diag(tvec) %*% XiPrime

par(mfrow=c(2,1))
matplot(dr1$Lat, Xi, main="Splines", xlab="Latitude")
matplot(dr1$Lat, XiReconstruct, main="Reconstructed Splines (with QR)", xlab="Latitude")

test = diag(tvec^(-1)) %*% backsolve(r=R1, x=Xi)
XiPrime[1,]
test[1,]
# since that looks so terrible, maybe it's because of conditioning?
# check the condition number
test = svd(G)
sqrt(test$d[1]/test$d[240])
test = svd(R1)
sqrt(test$d[1]/test$d[240])

# now try using svd instead of QR
test = svd(G %*% diag(tvec))
par(mfrow=c(1,1))
hist(test$d)
plot(1:240, cumsum(test$d)/sum(test$d), ylim=c(0,1), type="l", ylab="Condition Number", xlab="Rank of G %*% T", 
     main="99% Variance Explained")
abline(h=.99, col="red")
plot(1:240, sqrt(test$d[1]/test$d), type="l", log="y", ylab="Condition Number", xlab="Rank of G %*% T", 
     main="Condition Number is 1000")
abline(h=1000, col="red")
rank = min(which(cumsum(test$d)/sum(test$d) > .99999)) - 1
pseudoInv = test$v %*% diag(c(test$d[1:rank]^(-1), rep(0, 240-rank))) %*% t(test$u)

# generate OLS model of splines fit to subsidence data
mod = lm(dr1$subsidence ~ Xi-1)
summary(mod)

# make plots
par(mfrow=c(3,2))
for(i in 1:nrow(Xi)) {
  plotFault(csz, plotVar=pseudoInv %*% Xi[,i], main=paste0("Spline ", i), xlab="Longitude", ylab="Latitude", 
            xlim=lonRange, ylim=latRange)
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}
plotFault(csz, plotVar=pseudoInv %*% Xi %*% cbind(coef(mod)), main="Fitted Splines", 
          xlim=lonRange, ylim=latRange)

##### Do the same analysis but with low res fault

# generate spline basis functions in subsidence space
library(splines)
Xi = bs(dr1$Lat, df=5, intercept=TRUE)
par(mfrow=c(1,1))
matplot(dr1$Lat, Xi, main="Spline Basis", xlab="Latitude")

# make low res fault
fault = divideFault(faultGeom, nDown=1, nStrike=2)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

# embed spline basis on fault fault plane
GQR = qr(G)
Q1 = qr.Q(GQR)
R1 = qr.R(GQR)
R1Inv = backsolve(r = R1, x = diag(ncol(R1)))
embedMat = diag(tvec^(-1)) %*% R1Inv %*% t(Q1)
XiPrime = embedMat %*% Xi
hist(log10(abs(XiPrime)), breaks=100)

par(mfrow=c(3,2))
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="Observed Locking Rate (mm/yr)", 
           xlab="Longitude", ylab="Latitude", xlim=lonRange, ylim=latRange)
plotFault(fault, new=FALSE, plotData=FALSE)
for(i in 1:ncol(Xi)) {
  plotFault(fault, plotVar=XiPrime[,i], main=paste0("Spline ", i), xlab="Longitude", ylab="Latitude", 
            xlim=lonRange, ylim=latRange)
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}

# reconstruct Xi
XiReconstruct = G %*% diag(tvec) %*% XiPrime

# this works for:
# n=80 (2,2)
# doesn't work for:
# n=120 (2,3), n=120 (3,2)
par(mfrow=c(2,1))
matplot(dr1$Lat, Xi, main="Spline Basis", xlab="Latitude")
matplot(dr1$Lat, XiReconstruct, main="Reconstructed Spline Basis (with QR)", xlab="Latitude")

test = diag(tvec^(-1)) %*% backsolve(r=R1, x=Xi)
XiPrime[1,]
test[1,]
# since that looks so terrible, maybe it's because of conditioning?
# check the condition number of G %*% diag(tvec)
# n = 80 (2,2,8760120)
# check the condition number of G
# n = 80 (2,2,1144600)
test = svd(G)
sqrt(test$d[1]/test$d[n])
test = svd(G %*% diag(tvec))
sqrt(test$d[1]/test$d[n])
test = svd(R1)
sqrt(test$d[1]/test$d[n])

# now try using svd instead of QR
test = svd(G %*% diag(tvec))
par(mfrow=c(1,1))
hist(test$d)
plot(1:n, cumsum(test$d)/sum(test$d), ylim=c(0,1), type="l", ylab="Condition Number", xlab="Rank of G %*% T", 
     main="99% Variance Explained")
abline(h=.99, col="red")
plot(1:n, sqrt(test$d[1]/test$d), type="l", log="y", ylab="Condition Number", xlab="Rank of G %*% T", 
     main="Condition Number is 1000")
abline(h=1000, col="red")
# rank = min(which(cumsum(test$d)/sum(test$d) > .9999)) - 1
rank = n
sqrt(test$d[1]/test$d[rank]) # condition number
pseudoInv = test$v %*% diag(c(test$d[1:rank]^(-1), rep(0, n-rank))) %*% t(test$u)

# generate OLS model of splines fit to subsidence data
mod = lm(dr1$subsidence ~ Xi-1)
summary(mod)

# plot back-projected splines
par(mfrow=c(3,2))
for(i in 1:ncol(Xi)) {
  plotFault(fault, plotVar=pseudoInv %*% Xi[,i], main=paste0("Spline ", i), xlab="Longitude", ylab="Latitude", 
            xlim=lonRange, ylim=latRange)
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}
plotFault(fault, plotVar=pseudoInv %*% Xi %*% cbind(coef(mod)), main=paste0("Fitted Splines, rank ", rank), 
          xlim=lonRange, ylim=latRange)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# plot reconstructed spline basis
par(mfrow=c(2,1))
matplot(dr1$Lat, Xi, main="Spline Basis", xlab="Latitude")
matplot(dr1$Lat, G %*% diag(tvec) %*% pseudoInv %*% Xi, main="Reconstructed Spline Basis (with SVD)", xlab="Latitude")

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=pseudoInv %*% (-dr1$subsidence), xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE, na.rm = TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
newSigma = diag(pseudoInv %*% diag(dr1$Uncertainty) %*% t(pseudoInv))
plotFault(fault, plotVar=newSigma, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          logScale=TRUE, na.rm = TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Get minimum standard error
# ~2800 (2,2, rank=47), 237 (1,2, rank=40), 7.9 (1,1, rank=20)!
min(newSigma)

##### Now try penalized regression
library(lars)

# generate spline basis functions in subsidence space
library(splines)
Xi = bs(dr1$Lat, df=5, intercept=TRUE)
par(mfrow=c(1,1))
matplot(dr1$Lat, Xi, main="Spline Basis", xlab="Latitude")

# make low res fault
fault = divideFault(faultGeom, nDown=2, nStrike=2)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

##### do LASSO regression for regularized back-projection
GT = G %*% diag(tvec)
XiPrime = matrix(nrow=n, ncol=ncol(Xi))
for(i in 1:ncol(Xi)) {
  lassoMod = lars(GT, as.matrix(Xi[,i], ncol=1), type="lasso", intercept=FALSE, trace=TRUE, max.steps=5000)
  XiPrime[,i] = coef(lassoMod)[nrow(coef(lassoMod)),]
}

# plot back-projected splines
par(mfrow=c(3,2))
for(i in 1:ncol(XiPrime)) {
  plotFault(fault, plotVar=XiPrime[,i], main=paste0("Spline ", i), xlab="Longitude", ylab="Latitude", 
            xlim=lonRange, ylim=latRange)
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}
plotFault(fault, plotVar=XiPrime %*% cbind(coef(mod)), main=paste0("Fitted Splines, rank ", rank), 
          xlim=lonRange, ylim=latRange)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# plot reconstructed spline basis
par(mfrow=c(2,1))
matplot(dr1$Lat, Xi, main="Spline Basis", xlab="Latitude")
matplot(dr1$Lat, G %*% diag(tvec) %*% XiPrime, main="Reconstructed Spline Basis (with LASSO)", xlab="Latitude")

##### do ridge regression for regularized back-projection
library(glmnet)

# matrix for storing back-projected splines
XiPrime = matrix(nrow=n, ncol=ncol(Xi))

# use leave one out cross validation to get regularization param for ridge regression
# getRidgeLambda = function(lambdaGrid) {
getRidgeLambda = function(lambdaGrid) {
  allRMSE = matrix(nrow=m, ncol=0)
  for(i in 1:ncol(Xi)) {
    print(paste0("i: ", i))
    
    # leave one out cross-validation.  Returns rmse for each value of lambda
    LOOCV = function(j, lambdaGrid) {
      print(paste0("j: ", j))
      #       ridgeMod = glmnet(GT[-j,], as.matrix(Xi[-j,i], ncol=1), family="gaussian", alpha=0, 
      #                         intercept=FALSE, lambda=lambdaGrid)
      
      getSqErr = function(k) {
        ridgeMod = glmnet(GT[-j,], as.matrix(Xi[-j,i], ncol=1), family="gaussian", alpha=0, 
                          intercept=FALSE, lambda=lambdaGrid[k])
        pred = predict(ridgeMod, rbind(GT[j,]), type="link")
        (pred - Xi[j,i])^2
      }
      sapply(1:m, getSqErr)
    }
    
    # get rmse versus lambda matrix
    CVMat = sapply(1:nrow(Xi), LOOCV, lambdaGrid=lambdaGrid)
    allRMSE = cbind(allRMSE, rowMeans(CVMat))
  }
  # choose the best lambda
  rmse = rowMeans(allRMSE)
  minI = which.min(rmse)
  lambda = lambdaGrid[minI]
  list(lambda=lambda, rmse=rmse)
}

# range for lambda regularization parameter based on some previous results
lambdaMax = 100
lambdaMin = .01
m=100
lambdaGrid = rev(exp(seq(log(lambdaMin), log(lambdaMax), length=m)))

# get best value for lambda
# long story short: optimal lambda is near 0, so no point in ridge regression
out = getRidgeLambda(lambdaGrid)
plot(lambdaGrid, out$rmse, log="x", main="RMSE", xlab="lambda")

##### Test backprojecting subsidence data from just one event
event = dr1[dr1$event=="T5",]

# make low res fault
fault = divideFault(faultGeom, nDown=1, nStrike=1)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(event$Lon, event$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

# try using lm function (WLS)
wts = (1/event$Uncertainty^2)/sum(1/event$Uncertainty^2)
mod = lm(-event$subsidence ~ G %*% diag(tvec) - 1, weights=wts)
summary(mod)
ests = coef(summary(mod))[,1]
SEs = coef(summary(mod))[,2]

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(event$Lon, event$Lat, -event$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(event$Lon, event$Lat, event$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=ests, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=SEs, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

min(SEs)

##### now try using svd instead of QR
test = svd(G %*% diag(tvec))
# par(mfrow=c(1,1))
# hist(test$d)
# plot(1:n, cumsum(test$d)/sum(test$d), ylim=c(0,1), type="l", ylab="Condition Number", xlab="Rank of G %*% T", 
#      main="99% Variance Explained")
# abline(h=.99, col="red")
# plot(1:n, sqrt(test$d[1]/test$d), type="l", log="y", ylab="Condition Number", xlab="Rank of G %*% T", 
#      main="Condition Number is 1000")
# abline(h=1000, col="red")
# rank = min(which(cumsum(test$d)/sum(test$d) > .999)) - 1
rank = min(n, nrow(event))
sqrt(test$d[1]/test$d[rank]) # condition number
pseudoInv = test$v %*% diag(c(test$d[1:rank]^(-1), rep(0, n-rank))) %*% t(test$u)

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(event$Lon, event$Lat, -event$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(event$Lon, event$Lat, event$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=pseudoInv %*% (-event$subsidence), xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
newSigma = diag(pseudoInv %*% diag(event$Uncertainty) %*% t(pseudoInv))
plotFault(fault, plotVar=newSigma, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Get minimum standard error
# 122.9 (1,1,rank=full, T1), 6.6e34 (1,1,rank=full,T2), 1306 (1,1,full,T5)
min(newSigma)
# 5 (1,1, T1), 23/FAIL (1,1,T2), 15 (1,1,T5)
min(SEs)

##### try using generalized least squares
library(GLSME)

# get fit MLEs
lambda = fixedFitMVN$MLEs[1]
muZeta = fixedFitMVN$MLEs[2]
sigmaZeta = fixedFitMVN$MLEs[3]
lambda0 = fixedFitMVN$MLEs[4]
muXi = fixedFitMVN$MLEs[5]

# set other relevant parameters
nuZeta = 3/2 # Matern smoothness
phiZeta = 232.5722 # fit from fitGPSCovariance()

# get CSZ covariance matrix
xp = cbind(fault$longitude, fault$latitude)
Vt = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                    theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
Vt = G %*% diag(tvec) %*% Vt %*% diag(tvec) %*% t(G) # why negative eignvalues?
Ve = diag(event$Uncertainty^2)
glsMod = GLSME(-event$subsidence, G %*% diag(tvec), Vt, Ve, 0, 0, c(FALSE, FALSE), 
               CenterPredictor = FALSE, Vttype="Matrix")
estsGLS = glsMod$GLSestimate
SEsGLS = glsMod$errorGLSestim

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(event$Lon, event$Lat, -event$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(event$Lon, event$Lat, event$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=estsGLS, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=SEsGLS, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Make same plots except don't use log scale, only plot reasonable values
quilt.plot(event$Lon, event$Lat, -event$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(event$Lon, event$Lat, event$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=estsGLS, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", varRange=c(0, 100))
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=SEsGLS, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

min(SEsGLS)
estsGLS

##### Test backprojecting subsidence data from all events

# make low res fault
fault = divideFault(faultGeom, nDown=1, nStrike=1)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

# try using lm function (WLS)
wts = (1/dr1$Uncertainty^2)/sum(1/dr1$Uncertainty^2)
mod = lm(-dr1$subsidence ~ G %*% diag(tvec) - 1, weights=wts)
summary(mod)
ests = coef(summary(mod))[,1]
SEs = coef(summary(mod))[,2]

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=ests*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=SEs*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

min(SEs)

##### now try using svd instead of QR
test = svd(G %*% diag(tvec))
# par(mfrow=c(1,1))
# hist(test$d)
# plot(1:n, cumsum(test$d)/sum(test$d), ylim=c(0,1), type="l", ylab="Condition Number", xlab="Rank of G %*% T", 
#      main="99% Variance Explained")
# abline(h=.99, col="red")
# plot(1:n, sqrt(test$d[1]/test$d), type="l", log="y", ylab="Condition Number", xlab="Rank of G %*% T", 
#      main="Condition Number is 1000")
# abline(h=1000, col="red")
# rank = min(which(cumsum(test$d)/sum(test$d) > .999)) - 1
rank = min(n, nrow(dr1))
sqrt(test$d[1]/test$d[rank]) # condition number
pseudoInv = test$v %*% diag(c(test$d[1:rank]^(-1), rep(0, n-rank))) %*% t(test$u)

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=(pseudoInv %*% (-dr1$subsidence))*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
newSigma = diag(pseudoInv %*% diag(dr1$Uncertainty) %*% t(pseudoInv))
plotFault(fault, plotVar=newSigma*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Get minimum standard error
min(newSigma)
min(SEs)

##### try using generalized least squares
library(GLSME)

# get fit MLEs
lambda = fixedFitMVN$MLEs[1]
muZeta = fixedFitMVN$MLEs[2]
sigmaZeta = fixedFitMVN$MLEs[3]
lambda0 = fixedFitMVN$MLEs[4]
muXi = fixedFitMVN$MLEs[5]

# set other relevant parameters
nuZeta = 3/2 # Matern smoothness
phiZeta = 232.5722 # fit from fitGPSCovariance()

# get CSZ EQ slip covariance matrix
xp = cbind(fault$longitude, fault$latitude)
SigmaZeta = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                           theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
slipParams = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G, tvec)
Vt = slipParams$Sigma

# residuals variance matrix
Ve = diag(dr1$Uncertainty^2)

# do GLS
glsMod = GLSME(-dr1$subsidence, G %*% diag(tvec), Vt, Ve, 0, 0, c(FALSE, FALSE), 
               CenterPredictor = FALSE, Vttype="Matrix")
estsGLS = glsMod$GLSestimate
SEsGLS = glsMod$errorGLSestim
# save(glsMod, file="fixedMVNGLSMod.RData")

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=estsGLS*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=SEsGLS*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Make same plots except don't use log scale, only plot reasonable values
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=estsGLS*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=SEsGLS*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

min(SEsGLS)
estsGLS

##### do same as above but with nonnegative weighted least squares (nnwls)
library(nnls)

# make low res fault
fault = divideFault(faultGeom, nDown=1, nStrike=1)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

# scale A and y by sqrt of weights (1/sigma)
# W^(-1/2) y = W^(-1/2) A beta + W^(-1/2) epsilon
scalings = dr1$Uncertainty
A = diag(scalings) %*% G %*% diag(tvec)
y = diag(scalings) %*% (-dr1$subsidence)
nnMod = nnls(A, y)

nnEsts = nnMod$x
nnSEs = nnlsBootstrapSE(A, y, nnMod$fitted, nSamples=10000)

min(nnSEs) # .49

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
tmp = nnEsts*tvec
tmp[nnEsts < .1] = .1
plotFault(fault, plotVar=tmp, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=nnSEs*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          varRange=c(0.1, max(nnSEs)), logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Make same plots except don't use log scale, only plot reasonable values
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=nnEsts*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=nnSEs*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

##### do nonnegative GLS
library(nnls)

# make low res fault
fault = divideFault(faultGeom, nDown=1, nStrike=1)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

# compute subsidence variance matrix and eigenDecomp
subVarMat = getSubsidenceVarianceMat(fixedFitMVN$MLEs, fault=fault, G=G)
V = subVarMat$decomp$vectors
Lambda = subVarMat$decomp$values
Sigma = subVarMat$Sigma

# compute GLS scaling matrix
# NOTE: t(V) = V^(-1) since Sigma symmetric and real
# scalingMat = V %*% diag(Lambda^(-0.5)) %*% t(V)
# scalingMat = Lambda^(-0.5)
# scalingMat[abs(Lambda) < 10^(-7)] = 0 # pseudoInverse so set certain parts to 0
# scalingMat = V %*% diag(scalingMat) %*% t(V)
# tmp = sqrt(pseudoinverse(diag(Lambda)))
# tmp = V %*% tmp %*% t(V)
# library(expm)
# scalingMat = pseudoinverse(Sigma, 10^(-7)) %^% 1/2
tmp = eigen(Sigma)
d = tmp$values
d[d < 10^(-7)] = 0
d12 = d^(-1/2)
d12[d < 10^-7] = 0
scalingMat = tmp$vectors %*% diag(d12) %*% t(tmp$vectors)

# scale A and y by inverse sqrt of Variance matrix
# Sigma^(-1/2) y = Sigma^(-1/2) A beta + Sigma^(-1/2) epsilon
A = scalingMat %*% G %*% diag(tvec)
y = scalingMat %*% (-dr1$subsidence)
nnMod = nnls(A, y)

nnEsts = nnMod$x
nnSEs = nnlsBootstrapSE(A, y, nnMod$fitted, nSamples=10000)

min(nnSEs) # 523 for expm code, 36 for new code

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
tmp = nnEsts*tvec
tmp[nnEsts < .1] = .1
plotFault(fault, plotVar=tmp, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)", logScale=TRUE)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=nnSEs*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10, 
          varRange=c(0.1, max(nnSEs)))
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Make same plots except don't use log scale, only plot reasonable values
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=nnEsts*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift Data (m)")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=nnSEs*tvec, xlim=lonRange, ylim=latRange, 
          main="Back-Projected Coastal Uplift SE (m)", legend.mar = 10)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

##### Nonnegative LASSO regression
library(penalized)
penalProfile = profL1(y, A, positive=TRUE)

##### Using INLA, setting lognormal prior and exponential priors for coefs
library(INLA)

# make low res fault
fault = divideFault(faultGeom, nDown=1, nStrike=1)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

# compute subsidence variance matrix and eigenDecomp
subVarMat = getSubsidenceVarianceMat(fixedFitMVN$MLEs, fault=fault, G=G)
V = subVarMat$decomp$vectors
Lambda = subVarMat$decomp$values
Sigma = subVarMat$Sigma

# compute GLS scaling matrix
# NOTE: t(V) = V^(-1) since Sigma symmetric and real
# scalingMat = V %*% diag(Lambda^(-0.5)) %*% t(V)
# scalingMat = Lambda^(-0.5)
# scalingMat[abs(Lambda) < 10^(-7)] = 0 # pseudoInverse so set certain parts to 0
# scalingMat = V %*% diag(scalingMat) %*% t(V)
# tmp = sqrt(pseudoinverse(diag(Lambda)))
# tmp = V %*% tmp %*% t(V)
# library(expm)
# scalingMat = pseudoinverse(Sigma, 10^(-7)) %^% 1/2
tmp = eigen(Sigma)
d = tmp$values
d[d < 10^(-7)] = 0
d12 = d^(-1/2)
d12[d < 10^-7] = 0
scalingMat = tmp$vectors %*% diag(d12) %*% t(tmp$vectors)

# scale A and y by inverse sqrt of Variance matrix
# Sigma^(-1/2) y = Sigma^(-1/2) A beta + Sigma^(-1/2) epsilon
A = scalingMat %*% G %*% diag(tvec)
y = scalingMat %*% (-dr1$subsidence)

# define exponential prior log-density
exponential = "expression:
scale = 1;
return(-x/scale - log(scale));"

# define lognormal prior log-density
lognormal = "expression:
mu=0; s=1;
logFrac = -log(sqrt(2*pi)*x*s);
logPropDens = -(log(x)-mu)^2/(2*s^2);
return(logFrac + logPropDens);"

# get scale parameter for fault grid cells
# NOTE: tau corresponds to before tapering, since zeta itself is not tapered
Ai = faultGeom$length * faultGeom$width
area = sum(Ai)
mu = 4e10 # rigidty (value as assumed by Randy)
S = mean(M0, na.rm=TRUE)/(mu * sum(Ai*tvec)) # pre-taper avg-slip
Si = S*tvec # post-taper avg-slip (not used)
tau = S

# get mean and variance parameters for lognormal prior
# pre-taper:
mu = mean(log(M0), na.rm=TRUE) - log(mu * sum(Ai*tvec))
s = sd(log(M0), na.rm=TRUE)
# post-taper (not used):
mui = mu + log(tvec)

# set prior
exp.prior = list(prec = list(params=tau))
ln.prior = list(prec = list(params=c(mu,s)))

# do INLA
ids = 1:ncol(A)
expFormula = y ~ -1 + f(ids, model="iid", hyper = exp.prior)
expMod = inla(expFormula, data = list(y=y, ids=ids), control.predictor = list(A=A)) 
summary(expMod)

# this is the prior for the precision of beta
param.beta = list(prec = list(param = c(1.0e-3, 1.0e-3)))

AFrame = data.frame(A)
AFrame$beta1 = rep(1, nrow(Aframe))
AFrame$beta2 = rep(2, nrow(Aframe))
AFrame$beta3 = rep(3, nrow(Aframe))
AFrame$beta4 = rep(4, nrow(Aframe))
AFrame$beta5 = rep(5, nrow(Aframe))
AFrame$beta6 = rep(6, nrow(Aframe))
AFrame$beta7 = rep(7, nrow(Aframe))
AFrame$beta8 = rep(8, nrow(Aframe))
AFrame$beta9 = rep(9, nrow(Aframe))
AFrame$beta10 = rep(10, nrow(Aframe))
AFrame$beta11 = rep(11, nrow(Aframe))
AFrame$beta12 = rep(12, nrow(Aframe))
AFrame$beta13 = rep(13, nrow(Aframe))
AFrame$beta14 = rep(14, nrow(Aframe))
AFrame$beta15 = rep(15, nrow(Aframe))
AFrame$beta16 = rep(16, nrow(Aframe))
AFrame$beta17 = rep(17, nrow(Aframe))
AFrame$beta18 = rep(18, nrow(Aframe))
AFrame$beta19 = rep(19, nrow(Aframe))
AFrame$beta20 = rep(20, nrow(Aframe))
expFormula = y ~ -1 + f(beta1, X1, model="iid", values = 1:20,
                        hyper = param.beta) +
  f(beta2, X2, copy="beta1", fixed=T) +
  f(beta3, X3, copy="beta1", fixed=T) +
  f(beta4, X4, copy="beta1", fixed=T) +
  f(beta5, X5, copy="beta1", fixed=T) +
  f(beta6, X6, copy="beta1", fixed=T) +
  f(beta7, X7, copy="beta1", fixed=T) +
  f(beta8, X8, copy="beta1", fixed=T) +
  f(beta9, X9, copy="beta1", fixed=T) +
  f(beta10, X10, copy="beta1", fixed=T) +
  f(beta11, X11, copy="beta1", fixed=T) +
  f(beta12, X12, copy="beta1", fixed=T) +
  f(beta13, X13, copy="beta1", fixed=T) +
  f(beta14, X14, copy="beta1", fixed=T) +
  f(beta15, X15, copy="beta1", fixed=T) +
  f(beta16, X16, copy="beta1", fixed=T) +
  f(beta17, X17, copy="beta1", fixed=T) +
  f(beta18, X18, copy="beta1", fixed=T) +
  f(beta19, X19, copy="beta1", fixed=T) +
  f(beta20, X20, copy="beta1", fixed=T)

##### Using Stan to estimate backprojection model

## first transform data so it is uncorrelated standard normal
# make low res fault
fault = divideFault(faultGeom, nDown=1, nStrike=1)
n = nrow(fault)

# get Okada linear transformation matrix
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# get taper values
load("fixedFit_MVN.RData")
tvec = taper(fault$depth, lambda=fixedFitMVN$lambdaMLE)

# compute subsidence variance matrix and eigenDecomp
subVarMat = getSubsidenceVarianceMat(fixedFitMVN$MLEs, fault=fault, G=G)
V = subVarMat$decomp$vectors
Lambda = subVarMat$decomp$values
Sigma = subVarMat$Sigma

# compute GLS scaling matrix
# NOTE: t(V) = V^(-1) since Sigma symmetric and real
# scalingMat = V %*% diag(Lambda^(-0.5)) %*% t(V)
# scalingMat = Lambda^(-0.5)
# scalingMat[abs(Lambda) < 10^(-7)] = 0 # pseudoInverse so set certain parts to 0
# scalingMat = V %*% diag(scalingMat) %*% t(V)
# tmp = sqrt(pseudoinverse(diag(Lambda)))
# tmp = V %*% tmp %*% t(V)
# library(expm)
# scalingMat = pseudoinverse(Sigma, 10^(-7)) %^% 1/2
tmp = eigen(Sigma)
d = tmp$values
d[d < 10^(-7)] = 0
d12 = d^(-1/2)
d12[d < 10^-7] = 0
scalingMat = tmp$vectors %*% diag(d12) %*% t(tmp$vectors)

# scale A and y by inverse sqrt of Variance matrix
# Sigma^(-1/2) y = Sigma^(-1/2) A beta + Sigma^(-1/2) epsilon
A = scalingMat %*% G %*% diag(tvec)
y = scalingMat %*% (-dr1$subsidence)

## now calculate prior for Stan with Goldfinger 2012 data
# load rstan and precompile parallel MC code for faster fitting in parallel
library(rstan)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

# first try exponential prior:
goldfinger = read.csv("goldfinger2012Table8.csv", header=TRUE)
dr1Events = c(1:10, 12:17, 19:20, 22, 24, 28:29)
dr1Events = 1:31 # only keep events after and including T13
goldfinger = goldfinger[dr1Events,]
M0 = goldfinger$seismicMoment

# convert from dyne cm to N m
M0 = M0 * 10^(-7)

# get taper values
# 2 is a reasonable value for lambda (taper parameter) based on all the fits
lambda = fixedFitMVN$lambdaMLE
tvec = taper(faultGeom$depth, lambda=lambda)

# get rate parameter for fault grid cells
# NOTE: this is before tapering, since zeta itself is not tapered
Ai = faultGeom$length * faultGeom$width
area = sum(Ai)
mu = 4e10 # rigidty (value as assumed by Randy)
S = mean(M0, na.rm=TRUE)/(mu * sum(Ai*tvec)) # pre-taper avg-slip
Si = S*tvec # post-taper avg-slip
tau = 1/S

# set up required data variables in the R environment for Stan
n = length(y)
p = nrow(faultGeom)
X = A
y = as.vector(y)

# fit the regression model with Stan (exponential model)
expFit <- stan(file = "backprojectExp.stan", iter=100000)

# show results
print(expFit, digits = 1)

# fit the regression model with Stan (lognormal model)
priorMu = mean(log(M0/(mu * sum(Ai*tvec))), na.rm=TRUE)
priorSigma = sd(log(M0/(mu * sum(Ai*tvec))), na.rm=TRUE)
lnFit <- stan(file = "backprojectLN.stan", iter=100000)

# show results
print(lnFit, digits = 1)
# NOTE: results seem highly dependent on prior.  Exponential prior 
#       seems to get more stable and reasonable results.  Lognormal 
#       results aren't terribly unreasonable, though, since the 
#       uncertainty in the coefficients is very high and reasonable 
#       results can't be ruled out.

# save results
# save(expFit, file="tmpExpFit.RData")
# save(lnFit, file="tmpLNFit.RData")
save(expFit, file="expFit.RData")
save(lnFit, file="LNFit.RData")

# load results back, if necessary
# load("tmpExpFit.RData")
# load("tmpLNFit.RData")
load("expFit.RData")
load("LNFit.RData")

# extract the samples of the relevant parameters
expTab <- extract(expFit, permuted = TRUE) # return a list of arrays 
expBetaTab = expTab$beta
expLogBetaTab = expTab$logBeta

expBetaEsts = colMeans(expBetaTab)
expBetaSD = apply(expBetaTab, 2, sd)
expBeta025 = apply(expBetaTab, 2, quantile, probs=0.025)
expBeta975 = apply(expBetaTab, 2, quantile, probs=0.975)

# plot results
par(mfrow=c(1,1))
# plot biggest parameters first
hist(expBetaTab[,6], freq=FALSE)
hist(expBetaTab[,10], freq=FALSE)
hist(expBetaTab[,15], freq=FALSE)
# now the smallest
hist(expBetaTab[,1], freq=FALSE)
hist(expBetaTab[,7], freq=FALSE)
hist(expBetaTab[,8], freq=FALSE)
# plot biggest parameters first
hist(expLogBetaTab[,6], freq=FALSE)
hist(expLogBetaTab[,10], freq=FALSE)
hist(expLogBetaTab[,15], freq=FALSE)
# now the smallest
hist(expLogBetaTab[,1], freq=FALSE)
hist(expLogBetaTab[,7], freq=FALSE)
hist(expLogBetaTab[,8], freq=FALSE)

# plot back-projected subsidence data with uncertainty
par(mfrow=c(2,2))
# GPS dat
# quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="GPS Observations (mm/yr)", 
#            xlim=lonRange, ylim=latRange)
# map("world", "Canada", add=TRUE, lwd=1.5)
# US(add=TRUE, lwd=1.5)
# plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence dat
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=expBetaEsts, xlim=lonRange, ylim=latRange, 
          main="Estimated Coastal Uplift Data (m)")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=expBetaSD, xlim=lonRange, ylim=latRange, 
          main="Estimated Coastal Uplift SD (m)", legend.mar = 10)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Now plot lower 2.5th quantile and upper 97.5th quantile
quilt.plot(dr1$Lon, dr1$Lat, -dr1$subsidence, main="Observed Uplift (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# subsidence unvertainty
quilt.plot(dr1$Lon, dr1$Lat, dr1$Uncertainty, main="Subsidence Uncertainty (m)", 
           xlim=lonRange, ylim=latRange, nx=40, ny=40)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected subsidence data
plotFault(fault, plotVar=expBeta025, xlim=lonRange, ylim=latRange, 
          main="Coastal Uplift 2.5th Percentile (m)")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(fault, plotData=FALSE, new=FALSE)
# backprojected uncertainty
plotFault(fault, plotVar=expBeta975, xlim=lonRange, ylim=latRange, 
          main="Coastal Uplift 97.5th Percentile (m)", legend.mar = 10)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

##### Let's take a step back and plot the G matrices.  Each row corresponds to a fault grid 
##### cell and each column represents a subsidence observation.

MLEs = fixedFitMVN$MLEs

nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(faultGeom, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=MLEs[4])
Gcsz = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=MLEs[4])

# The images of both G and Gcsz show a ridiculous amount of sparsity.
image(G)
image(cszG)

# Plot infinity norm of each row vector
GMax = apply(abs(G), 2, max)
cszGMax = apply(abs(cszG), 2, max)

# plot the max value (infinity norm) of each row
plot(GMax, main="Inifinty norm for rows of G", ylab="Infinity norm", xlab="Row index", 
     cex=.4, pch=19)
abline(h=.005, col="blue")
plot(cszGMax, main="Inifinty norm for rows of G", ylab="Infinity norm", xlab="Row index", 
     cex=.4, pch=19)
abline(h=.005, col="blue")

# plot the ecdfs of the infinity norms of G
FG = ecdf(GMax)
FcszG = ecdf(cszGMax)
xs = seq(0, max(c(GMax, cszGMax)), l=500)
plot(xs, FG(xs), main="ECDF for inifinty norm for rows of G", xlab="Infinity norm", 
     ylab="Probability", pch=19, cex=.4)
abline(h=.005, col="blue")
plot(xs, FcszG(xs), main="ECDF for inifinty norm for rows of G", xlab="Infinity norm", 
     ylab="Probability", pch=19, cex=.4)
abline(h=.005, col="blue")

hist(GMax, breaks=100)
abline(v=.005, col="blue")
mean(GMax < .005)
sum(GMax < .005)

hist(cszGMax, breaks=100)
abline(v=.005, col="blue")
mean(cszGMax < .005)
sum(cszGMax < .005)

##### Compare full conditional predictions to the GPS conditional predictions
load("fixedFit_MVN.RData")
MLEs = fixedFitMVN$MLEs

# compute taper
tvec = taper(slipDatCSZ$Depth, lambda=MLEs[1])
tvecX = tvec
tvecB = taper(csz$depth, lambda=MLEs[1])

# get GPS data and CSZ prediction coordinates
xd = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
xp = cbind(csz$longitude, csz$latitude)

# generate the predictions
gpsPreds = predsGivenGPS(MLEs, 10)
fullPreds = genFullPredsMVN(MLEs, 1000)
fullPreds10 = genFullPredsMVN(MLEs, 1000, 100) # add additional uncertainty to gps dat SD

# compute the conditional means (areal and pointwise)
gpsMucX = exp(gpsPreds$mucGPS+diag(gpsPreds$SigmaGPS)/2) * tvec
gpsMucB = gpsPreds$meanSlip
fullMucX = fullPreds$meanSlipGPS
fullMucB = fullPreds$meanSlip
fullMucX10 = fullPreds10$meanSlipGPS
fullMucB10 = fullPreds10$meanSlip

# compute the conditional SDs for the normals
gpsSDX = sqrt(diag(gpsPreds$SigmaGPS))
gpsSDB = sqrt(diag(gpsPreds$Sigma))
fullSDX = sqrt(diag(fullPreds$SigmaGPS))
fullSDB = sqrt(diag(fullPreds$Sigma))
fullSDX10 = sqrt(diag(fullPreds10$SigmaGPS))
fullSDB10 = sqrt(diag(fullPreds10$Sigma))

# compute the conditional SDs for the lognormal tapered slips (m is LN mean not N mean)
lnSD = function(m, s) { sqrt((exp(s^2)-1)*m^2) }
gpsSDX = lnSD(gpsMucX, gpsSDX) * tvecX
gpsSDB = lnSD(gpsMucB, gpsSDB) * tvecB
fullSDX = lnSD(fullMucX, fullSDX) * tvecX
fullSDB = lnSD(fullMucB, fullSDB) * tvecB
fullSDX10 = lnSD(fullMucX10, fullSDX10) * tvecX
fullSDB10 = lnSD(fullMucB10, fullSDB10) * tvecB

## plot results
# GPS coordinate plots
par(mfrow=c(2,2))
quilt.plot(xd, gpsMucX, main="GPS Mean")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(faultGeom, plotData=FALSE, new=FALSE)

quilt.plot(xd, gpsSDX, main="GPS SD")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(faultGeom, plotData=FALSE, new=FALSE)

quilt.plot(xd, fullMucX10, main="Full Mean")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(faultGeom, plotData=FALSE, new=FALSE)

quilt.plot(xd, fullSDX10, main="Full SD")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(faultGeom, plotData=FALSE, new=FALSE)

# fault cell plots
par(mfrow=c(2,2))
plotFault(csz, gpsMucB, main="GPS Mean")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

plotFault(csz, gpsSDB, main="GPS SD")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

plotFault(csz, fullMucB10, main="Full Mean")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

plotFault(csz, fullSDB10, main="Full SD")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# compute moment magnitudes
mags = apply(fullPreds$slipSims, 2, getMomentFromSlip)
mags10 = apply(fullPreds10$slipSims, 2, getMomentFromSlip)

# plot slip simulations
par(mfrow=c(2,2))
plotFault(csz, fullPreds$slipSims[,1], main=paste0("Simulation 1, Mo=", mags[1]))
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

plotFault(csz, fullPreds$slipSims[,2], main=paste0("Simulation 1, Mo=", mags[2]))
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

plotFault(csz, fullPreds$slipSims[,3], main=paste0("Simulation 1, Mo=", mags[3]))
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

plotFault(csz, fullPreds$slipSims[,4], main=paste0("Simulation 1, Mo=", mags[4]))
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# Generate subsidence predictions
gpsSubs = predsToSubsidence(MLEs, gpsPreds)
fullSubs = predsToSubsidence(MLEs, fullPreds, useMVNApprox = FALSE)
fullSubs10 = predsToSubsidence(MLEs, fullPreds10, useMVNApprox = FALSE)
comparePredsToSubs(MLEs, slipPreds=fullPreds, subPreds=fullSubs, nsim=1000, 
                   plotNameRoot="full", savePlots=FALSE)
comparePredsToSubs(MLEs, slipPreds=fullPreds10, subPreds=fullSubs10, nsim=1000, 
                   plotNameRoot="full", savePlots=FALSE)

##### test updateMu function (6665 seconds for 5000 iter)
system.time(test <- updateMu(MLEs, fault=csz, niter=10))
system.time(newMu <- updateMu(MLEs, fault=csz, niter=5000))
save(newMu, file="newMu.RData")
load("newMu.RData")

clim1=range(c(newMu$muMat[,19], newMu$newMu))
plotFault(csz, newMu$muMat[,19], varRange = clim1)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(csz, newMu$newMu, varRange = clim1)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

clim2=range(c(sqrt(diag(newMu$Sigmas[[19]])), sqrt(diag(newMu$varMat))))
plotFault(csz, sqrt(diag(newMu$Sigmas[[19]])), varRange = clim2)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
plotFault(csz, sqrt(diag(newMu$varMat)), varRange = clim2)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

#
sdMat = newMu$muMat
for(i in 1:ncol(newMu$muMat)) {
  sdMat[,i] = sqrt(diag(newMu$Sigmas[[i]]))
}
VarMat = sdMat^2
W = 1/VarMat

## calculate new mean estimate based only on subsidence and its variance matrix
sampleSize=19
W = W[,1:(sampleSize-1)]
W = sweep(W, 1, rowSums(W), "/")
newMuSub = rowSums(W*newMu$muMat[,1:(sampleSize-1)])
varMatSub = matrix(0, nrow=nrow(csz), ncol=nrow(csz))
for(k in 1:(sampleSize-1)) {
  thisVar = newMu$Sigmas[[k]]
  thisWeights = W[,k]
  thisVar = sweep(sweep(thisVar, 1, thisWeights, "*"), 2, thisWeights, "*")
  varMatSub = varMatSub + thisVar
}

# GPS mean
clim1=range(c(newMu$muMat[,19], newMu$newMu))
plotFault(csz, newMu$muMat[,19], varRange = clim1)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
# overall mean
plotFault(csz, newMu$newMu, varRange = clim1)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
# subsidence mean
plotFault(csz, newMuSub)
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
# subsidence SD
plotFault(csz, sqrt(diag(varMatSub)))
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)

# compare predictions versus subsidence data
par(mfrow=c(1,2))
T1 = dr1[events == "T1",]
plotFault(csz, exp(newMu$muMat[,18]), main="Median T1 Predicted Slip (m)")
map("world", "Canada", add=TRUE, lwd=1.5)
US(add=TRUE, lwd=1.5)
points(T1$Lon, T1$Lat, pch="+", col="red")

quilt.plot(T1$Lon, T1$Lat, T1$subsidence)
map("world", "Canada", add=TRUE, lwd=1.5, main="Observed T1 Subsidence (m)", 
    xlab="Longitude")
US(add=TRUE, lwd=1.5)

# get seismic moments for imputed earthquakes
mags = c()
usedEvents = c(uniqueEvents[-c(1,2,5)], "GPS", "Avg")
for(i in 1:ncol(newMu$muMat)) {
  mags = c(mags, getMomentFromSlip(exp(newMu$muMat[,i])))
}
mags = c(mags, print(getMomentFromSlip(exp(newMu$newMu))))
print(cbind(usedEvents, mags))

# do traceplots of the stan results:
traceplot(newMu$allStanResults[[1]])
traceplot(newMu$allStanResults[[2]], inc_warmup=TRUE, 
          pars=c("beta[181]", "beta[182]", "beta[183]", 
                 "beta[184]", "beta[185]", "beta[186]", 
                 "beta[187]", "beta[188]", "beta[189]"))

system.time(test <- fitModelIterative(maxIter=2, niterMCMC=250)) #params fit maxIter times, mean fit maxIter-1 times
save(testInputs, file="testInputs.RData")
load("testInputs.RData")
parFit = doFixedFit(testInputs[[1]], testInputs[[2]], testInputs[[3]], testInputs[[4]], testInputs[[5]], testInputs[[6]])

# check results of iterative model fit:
out = load("fullIterFit2.RData")
muMatCSZ = state$muMatCSZ
muMatGPS = state$muMatGPS
parMat = state$parMat

# plot both pointwise mean slips and block averages over CSZ geometry
for(i in 1:ncol(muMatGPS)) {
  meanVec = muMatGPS[,i]
  sigmaZeta = parMat[3,i]
  lambda = parMat[1,i]
  untaperedMeanSlip = exp(meanVec + sigmaZeta^2/2)
  taperedMeanSlip = taper(slipDatCSZ$Depth, lambda)*untaperedMeanSlip
  quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, taperedMeanSlip, main=paste0("Mean slip at iteration ", i), 
             zlim=c(0,30))
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
  plotFault(csz, plotData = FALSE, new=FALSE)
}
for(i in 1:ncol(muMatCSZ)) {
  meanVec = muMatCSZ[,i]
  sigmaZeta = parMat[3,i]
  lambda = parMat[1,i]
  untaperedMeanSlip = exp(meanVec + sigmaZeta^2/2)
  taperedMeanSlip = taper(csz$depth, lambda)*untaperedMeanSlip
  plotFault(csz, taperedMeanSlip, main=paste0("Mean slip at iteration ", i), varRange=c(0, 30))
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}

# plot 95th percentiles of slip over both pointwise and block CSZ geometry
for(i in 1:ncol(muMatGPS)) {
  meanVec = muMatGPS[,i]
  sigmaZeta = parMat[3,i]
  lambda = parMat[1,i]
  untaperedHiSlip = exp(meanVec + qnorm(.95)*sigmaZeta)
  taperedHiSlip = taper(slipDatCSZ$Depth, lambda)*untaperedHiSlip
  quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, taperedHiSlip, main=paste0("Slip 95th percentile at iteration ", i), 
             zlim=c(0,100))
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
  plotFault(csz, plotData = FALSE, new=FALSE)
}
for(i in 1:ncol(muMatCSZ)) {
  meanVec = muMatCSZ[,i]
  sigmaZeta = parMat[3,i]
  lambda = parMat[1,i]
  untaperedHiSlip = exp(meanVec + qnorm(.95)*sigmaZeta)
  taperedHiSlip = taper(csz$depth, lambda)*untaperedHiSlip
  plotFault(csz, taperedHiSlip, main=paste0("Slip 95th percentile at iteration ", i), 
            varRange=c(0,100))
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}

# get magnitudes of the mean and 95th percentile slips for each iteration
meanMags = c()
hiMags = c()
for(i in 1:ncol(muMatCSZ)) {
  meanVec = muMatCSZ[,i]
  sigmaZeta = parMat[3,i]
  lambda = parMat[1,i]
  untaperedMeanSlip = exp(meanVec + sigmaZeta^2/2)
  taperedMeanSlip = taper(csz$depth, lambda)*untaperedMeanSlip
  untaperedHiSlip = exp(meanVec + qnorm(.95)*sigmaZeta)
  taperedHiSlip = taper(csz$depth, lambda)*untaperedHiSlip
  meanMags[i] = getMomentFromSlip(taperedMeanSlip)
  hiMags[i] = getMomentFromSlip(taperedHiSlip)
}
print(cbind(meanMags, hiMags))

# explore likelihood of the fit model through the iterations:
print(tail(state$optimTables[[1]][,4:6]))
for(i in 2:length(state$optimTables)) {
  print(tail(state$optimTables[[i]][,3:5]))
}
state$logLiks
# NOTE: the log likelihoods seem to reach a minimum for some reason
parMat = state$parMat
rownames(parMat) = c("lambda", "muZeta", "sigmaZeta", "lambda0", "muXi")
parMat
# parameters seem to converge (with the exception of muXi, 
# since the mean vector is changing).  sigmaZeta decreases monotonically (?), 
# as expected

##### Check to see if imputed data makes sense:
load("fullIterFitPrior.RData")

nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

tvec = taper(csz$depth, lambda=0.25)

test = updateMuGivenStan(subFunIns$params, state$lastStanResults, subFunIns$muVec, G=G)
# Only use events that have at least 5 observations (leaves out T12, T11, and T9a)
minObs = 5
numObs = table(subDat$event)
threshUniqueEvents = uniqueEvents[numObs >= minObs]

for(i in 1:length(threshUniqueEvents)) {
  thisEvent = threshUniqueEvents[i]
  subDatInds = events == thisEvent
  subDat = dr1[subDatInds,]
  thisG = G[subDatInds,]
  thisMean = exp(test$muMat[,i] + diag(test$Sigmas[[i]])/2)
  thisSubMean = thisG %*% (tvec * thisMean)
  
  thisRange = range(c(-subDat$subsidence, thisSubMean))
  plot(subDat$Lat, -subDat$subsidence, main=thisEvent, ylim=thisRange)
  points(subDat$Lat, thisSubMean, col="blue")
}

goldfinger = read.csv("goldfinger2012Table8.csv", header=TRUE)
dr1Events = c(1:10, 12:17, 19:20, 22, 24, 28:29)
dr1Events = 1:31 # only keep events after and including T13
goldfinger = goldfinger[dr1Events,]
goldCompareI = c(22, 19, 17, 16, 15, 14, 13:12, 10:4, 2:1)
threshCompareI = 2:length(threshCompareI)
goldfingerComp = goldfinger[goldCompareI,c(1, 4, 6, 7)]
#####
# test taper types:
dStar = 26000
test = getInitialSplineEsts(MLEs[2], MLEs[3], MLEs[1], G)
testNoNorm = getInitialSplineEsts(MLEs[2], MLEs[3], MLEs[1], G, normalizeTaper = FALSE)
testStar = getInitialSplineEsts(MLEs[2], MLEs[3], MLEs[1], G, dStar=dStar)
plotFault(csz, taper(csz$depth, lambda=test$lambdaVecCSZ), legend.mar = 10)
plotFault(csz, taper(csz$depth, lambda=testNoNorm$lambdaVecCSZ, normalize=FALSE), legend.mar = 10)
plotFault(csz, taper(csz$depth, lambda=testStar$lambdaVecCSZ, dStar=dStar), legend.mar=10)

# see what points can either have subsidence or uplift
hasUp = apply(subSims>0, 1, any)
hasSub = apply(subSims<0, 1, any)
sum(hasUp & hasSub)
which(hasUp & hasSub) # 328:336
range(dr1$Lat[328:336])

# examine distribution of predicted subsidences, see if it's exponential
subDat=dr1
subDat$Uncertainty = subDat$Uncertainty/2
reals = preds(params, nsim=10000, fault=csz, tvec=tvec)
subs = predsToSubsidence(params, reals, useMVNApprox=FALSE, G=G, subDat=subDat)
subSims = subs$noiseSims
i=50
hist(subSims[i,], breaks=250, freq=F)
logSims = log(abs(subSims[i,]))
cleanLogSims = logSims[is.finite(logSims)]
xs = seq(min(subSims[i,]), max(subSims[i,]), l=500)
# lines(xs, dlnorm(abs(xs), mean(cleanLogSims), sd(cleanLogSims)), col="green")
# lines(xs, dexp(abs(xs), 1/mean(abs(subSims[i,]))), col="red")
lines(xs, dnorm(xs, mean(subSims[i,]), sd(subSims[i,])), col="blue")

inRange = function(x, low, hi) {
  tmp = x[x > low]
  tmp[tmp < hi]
}
inRangeI = function(x, low, hi) {
  tmpI = x > low
  tmpI & (x < hi)
}

# test asymmetric laplace distribution fit
test = mleALD(subSims[i,])
parALD = test$par
lines(xs, dALD(xs, parALD[1], parALD[2], parALD[3]))

# how long does full MLE fit take? (roughly 200s for full dataset)
system.time(MLEs <- apply(subSims, 1, mleALD))

# try using MMEs
parASL = approxASL(subSims, breaks=1000)
mi = parASL$mHats[i]
lami = parASL$lambdaHats[i]
kappai = parASL$kappaHats[i]
lines(xs, dASL(xs, mi, lami, kappai), col="cyan")

skewVal = mean((subSims[i,] - mean(subSims[i,]))^3)
meanDiff = mean(subSims[i,])-mi
V = var(subSims[i,])
lamiTest = sqrt(2/(V - meanDiff^2))
y = lamiTest*meanDiff
kappaiTest = (-y+sqrt(y^2 + 4))/2
lines(xs, dASL(xs, mi, lamiTest, kappaiTest), col="green")

# Correlation scale of subsidence data is shorter than slip data.  This makes sense 
# because it is heavily influenced by how close the fault geometry is to the shore.
meanLon = mean(dr1$Lon)
minLonI = which.min(dr1$Lon)
phiZeta = 232.5722 # MLE based on fitGPSCovariance result
nuZeta = 3/2 # we assume this is the Matern smoothness parameter
coords = cbind(dr1$Lon, dr1$Lat)
dists = rdist.earth.vec(matrix(coords[minLonI,], ncol=2), coords)
plot(coords[,2], dists)
plot(coords[,2], Matern(dists/phiZeta, smoothness = nuZeta))


# test fit of splines to subsidence data in latitude space
# compute spline basis coefficients for mean of subsidence data
latRange = c(40, 50)
latSeq = seq(latRange[1], latRange[2], l=100)
nKnotsMean=4
nKnotsVar=4
# Xi = bs(dr1$Lat, df=nKnots, intercept=FALSE, Boundary.knots=latRange)
# Xi = cbind(rep(1, nrow(Xi)), Xi)
XiMean = bs(latSeq, df=nKnotsMean, intercept=FALSE, Boundary.knots=latRange)
XiMean = cbind(rep(1, nrow(XiMean)), XiMean)
XiVar = bs(latSeq, df=nKnotsVar, intercept=FALSE, Boundary.knots=latRange)
XiVar = cbind(rep(1, nrow(XiVar)), XiVar)

test = fitSplinesToSubDat(nKnotsMean=nKnotsMean, nKnotsVar=nKnotsVar)
meanSplinePar = test$meanSplinePar
varSplinePar = test$varSplinePar
meanVec = XiMean %*% meanSplinePar
sdVec = sqrt(XiVar %*% varSplinePar)
hi = meanVec + qnorm(.975)*sdVec
low = meanVec + qnorm(.025)*sdVec
plot(dr1$Lat, dr1$subsidence, pch=19, cex=.5, col="blue", 
     xlab="Latitude", ylab="Subsidence (m)", main="Spline Fit")
lines(latSeq, meanVec)
lines(latSeq, low, lty=2)
lines(latSeq, hi, lty=2)

# now make same plot but with observation noise added to confidence intervals
XiMean = bs(dr1$Lat, df=nKnotsMean, intercept=FALSE, Boundary.knots=latRange)
XiMean = cbind(rep(1, nrow(XiMean)), XiMean)
XiVar = bs(dr1$Lat, df=nKnotsVar, intercept=FALSE, Boundary.knots=latRange)
XiVar = cbind(rep(1, nrow(XiVar)), XiVar)

test = fitSplinesToSubDat(nKnotsMean=nKnotsMean, nKnotsVar=nKnotsVar)
meanSplinePar = test$meanSplinePar
varSplinePar = test$varSplinePar
meanVec = XiMean %*% meanSplinePar
sdVec = sqrt(XiVar %*% varSplinePar + dr1$Uncertainty^2)
hi = meanVec + qnorm(.975)*sdVec
low = meanVec + qnorm(.025)*sdVec
plot(dr1$Lat, dr1$subsidence, pch=19, cex=.5, col="blue", 
     xlab="Latitude", ylab="Subsidence (m)", main="Spline Fit")
sortI = sort(dr1$Lat, index.return=TRUE)$ix
lines(dr1$Lat[sortI], meanVec[sortI])
lines(dr1$Lat[sortI], low[sortI], lty=2)
lines(dr1$Lat[sortI], hi[sortI], lty=2)

# make sure prediction bands contain ~95% of data
mean(dr1$subsidence > hi) + mean(dr1$subsidence < low)

# now try moving average variance estimator (1 degree in width, .5 in each direction)
Xi = bs(dr1$Lat, df=nKnotsMean, intercept=FALSE, Boundary.knots=latRange)
Xi = cbind(rep(1, nrow(Xi)), Xi)
datMeans = Xi %*% meanSplinePar
varStats = (dr1$subsidence - datMeans)^2 - dr1$Uncertainty^2
# for any observation's index, estimate the true field's variance at that latitude
getVar = function(i) {
  thisLat = dr1$Lat[i]
  nearbyLats = rdist(thisLat, dr1$Lat) < .5
  varEst = mean(varStats[nearbyLats])
  
  # make sure varEst > 0
  if(varEst < 0)
    varEst = .001
  
  return(varEst)
}
# estimate variances for each data latitude here:
newVarVec = sapply(1:length(dr1$subsidence), getVar)

# plot new fit:
newSdVec = sqrt(newVarVec)
hi = datMeans + qnorm(.975)*newSdVec
low = datMeans + qnorm(.025)*newSdVec
plot(dr1$Lat, dr1$subsidence, pch=19, cex=.5, col="blue", 
     xlab="Latitude", ylab="Subsidence (m)", main="Moving Average + Spline Fit")
sortI = sort(dr1$Lat, index.return=TRUE)$ix
lines(dr1$Lat[sortI], datMeans[sortI])
lines(dr1$Lat[sortI], low[sortI], lty=2)
lines(dr1$Lat[sortI], hi[sortI], lty=2)

# now plot with additional variance from data noise:
newSdVec = sqrt(newVarVec + dr1$Uncertainty^2)
hi = datMeans + qnorm(.975)*newSdVec
low = datMeans + qnorm(.025)*newSdVec
plot(dr1$Lat, dr1$subsidence, pch=19, cex=.5, col="blue", 
     xlab="Latitude", ylab="Subsidence (m)", main="Moving Average + Spline Fit")
sortI = sort(dr1$Lat, index.return=TRUE)$ix
lines(dr1$Lat[sortI], datMeans[sortI])
lines(dr1$Lat[sortI], low[sortI], lty=2)
lines(dr1$Lat[sortI], hi[sortI], lty=2)

# make sure 95% confidence intervals contains ~95% of data
mean(dr1$subsidence > hi) + mean(dr1$subsidence < low)

## now do same analysis but with moving average mean and variance
# for any observation's index, estimate the true field's mean at that latitude
getMean = function(i) {
  thisLat = dr1$Lat[i]
  nearbyLats = rdist(thisLat, dr1$Lat) < .5
  mean(dr1$subsidence[nearbyLats])
}
datMeans = sapply(1:length(dr1$subsidence), getMean)
varStats = (dr1$subsidence - datMeans)^2 - dr1$Uncertainty^2
# for any observation's index, estimate the true field's mean/variance at that latitude
getVar = function(i) {
  thisLat = dr1$Lat[i]
  nearbyLats = rdist(thisLat, dr1$Lat) < .5
  varEst = mean(varStats[nearbyLats])
  
  # make sure varEst > 0
  if(varEst < 0)
    varEst = .001
  
  return(varEst)
}
# estimate variances for each data latitude here:
newVarVec = sapply(1:length(dr1$subsidence), getVar)

# plot new fit:
newSdVec = sqrt(newVarVec)
hi = datMeans + qnorm(.975)*newSdVec
low = datMeans + qnorm(.025)*newSdVec
plot(dr1$Lat, dr1$subsidence, pch=19, cex=.5, col="blue", 
     xlab="Latitude", ylab="Subsidence (m)", main="Moving Average Fit")
sortI = sort(dr1$Lat, index.return=TRUE)$ix
lines(dr1$Lat[sortI], datMeans[sortI])
lines(dr1$Lat[sortI], low[sortI], lty=2)
lines(dr1$Lat[sortI], hi[sortI], lty=2)

# now plot with additional variance from data noise:
newSdVec = sqrt(newVarVec + dr1$Uncertainty^2)
hi = datMeans + qnorm(.975)*newSdVec
low = datMeans + qnorm(.025)*newSdVec
plot(dr1$Lat, dr1$subsidence, pch=19, cex=.5, col="blue", 
     xlab="Latitude", ylab="Subsidence (m)", main="Moving Average Fit")
sortI = sort(dr1$Lat, index.return=TRUE)$ix
lines(dr1$Lat[sortI], datMeans[sortI])
lines(dr1$Lat[sortI], low[sortI], lty=2)
lines(dr1$Lat[sortI], hi[sortI], lty=2)

# make sure 95% confidence intervals contains ~95% of data
mean(dr1$subsidence > hi) + mean(dr1$subsidence < low)

## test linear interpolated coastline subsidence levels (too jagged?  
# within bounds of fit splines in subsidence space?)

# get subsidences from MLEs
# precompute GTest
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
GTest = okadaAll(csz, lonGrid, latGrid, allCoastCoords, slip=1, poisson=0.25)

# get taper
tvec = taper(csz$depth, lambda=MLEs[1])

# get zeta expectation
expectZeta = rep(exp(MLEs[2] + (MLEs[3]^2)/2), nrow(csz))*tvec
subOrig = G %*% expectZeta
subInterp = GTest %*% expectZeta

# plot original and interpolated subsidence along coast
par(mfrow=c(1,2))
latRange = c(40, 50)
subRange = range(c(subOrig, subInterp))
plot(subOrig, dr1$Lat, pch=19, col="blue", cex=.6, ylim=latRange, 
     xlim=subRange, main="Original Uplift Predictions", xlab="Uplift (m)", 
     ylab="Latitude")
plot(subInterp, allCoastLat, pch=19, col="red", cex=.6, ylim=latRange, 
     xlim=subRange, main="Interpolated Uplift Predictions", xlab="Uplift (m)", 
     ylab="Latitude")

# now plot with splines:
plot(subOrig, dr1$Lat, pch=19, col="blue", cex=.6, ylim=latRange, 
     xlim=subRange, main="Original Uplift Predictions", xlab="Uplift (m)", 
     ylab="Latitude")
lines(-meanVec, latSeq)
lines(-low, latSeq, lty=2)
lines(-hi, latSeq, lty=2)
plot(subInterp, allCoastLat, pch=19, col="red", cex=.6, ylim=latRange, 
     xlim=subRange, main="Interpolated Uplift Predictions", xlab="Uplift (m)", 
     ylab="Latitude")
lines(-meanVec, latSeq)
lines(-low, latSeq, lty=2)
lines(-hi, latSeq, lty=2)


##### optimization results
opt = list()
opt2 = list()
opt3 = list()
opt4 = list()
opt5 = list()
opt6 = list()
opt7 = list()
opt8 = list()
opt9 = list()
opt10 = list()
opt11 = list()
opt$par = c(2.338752,  2.04888,      12.01887, 6.779287, 17.78633, -22.77749, -14.14914)
opt2$par = c(2.188923,  2.088136,    15.09685, -2.153181, 23.38081, -30.15017, -14.09528)
opt3$par = c(2.675348,  1.968168,    17.89174, -9.98194, 30.99193, -37.8903, -14.33926)
opt4$par = c(2.575059,  1.989285,    22.5119, -25.45767, 45.53586, -51.1616, -15.96006)
opt5$par = c(3.254444,  1.715885,    27.61135, -47.44486, 45.46336, -58.6099, -18.85374)
opt6$par = c(3.074955,  1.747611, 29.56489, -52.63478, 46.45042, -61.34435, -20.66533)
opt7$par = c(3.10461,  1.745752, 32.04086, -58.023, 49.90471, -66.43788, -21.88596)
opt8$par = c(3.058119,  1.757687, 32.45259, -58.95977, 50.50164, -67.3081, -22.08851)
opt9$par = c(3.085913,  1.750588, 32.66448, -59.48048, 50.81436, -67.71813, -22.18403)
opt10$par = c(3.087351,  1.751152, 32.89431, -59.97645, 51.1278, -68.19173, -22.29191)
opt11$par = c(3.102128,  1.748229, 33.10879, -60.44, 51.4165, -68.62971, -22.39055)

testOpt = opt11
splinePar = testOpt$par[3:7]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)
points(dr1$Lon, dr1$Lat, pch="+", col="red", cex=.8)
sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))

# estimate muXi MLE with inverse variance weighting
logX = log(slipDatCSZ$slip)
ci = 1/sigmaXi^2
ci = ci/sum(ci)
muXi = sum((logX-testOpt$par[1])*ci)
params=c(NA, testOpt$par[1:2], 0.25, muXi, testOpt$par[3:7])

comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "testOpt11") ###### MAKE SURE TO CHANGE PLOT TITLES
#plotSlipDistribution(params, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=100000, fileNameRoot = "testOpt11")
#plotSlipDistribution(params, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=100000, logScale=TRUE, fileNameRoot = "testOpt11Log")


##### try interpolating GPS depths onto grid

# first make the grid
latRange=c(40, 50)
lonRange=c(-128, -122.5)
nx = 200
ny = 600
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
lonLatGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))

phiZeta = 232.5722
out = fastTps(cbind(slipDat$lon, slipDat$lat), slipDat$Depth, m=3, theta=phiZeta, lon.lat=TRUE, Dist.args=list(miles=FALSE, method="greatcircle"))

# plot the interpolated versus true observations
par(mfrow=c(1,2))
surface(out, zlim=range(slipDat$Depth), main="Interpolated CSZ Depth (m)", levels=c(10000, 25000, 50000))
quilt.plot(slipDat$lon, slipDat$lat, slipDat$Depth, main="Observed CSZ Depth (m)")

# make predictions on a grid
depths = predict(out, lonLatGrid)
negDepth = depths < 0
depths[negDepth] = 0
quilt.plot(lonLatGrid, depths)
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$Depth)

##### test pointwise predictions
## first get the parameters
testOpt = opt
splinePar = testOpt$par[3:7]
tvec = getTaperSpline(splinePar, nKnots=5, dStar=26000)
sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))

# estimate muXi MLE with inverse variance weighting
logX = log(slipDatCSZ$slip)
ci = 1/sigmaXi^2
ci = ci/sum(ci)
muXi = sum((logX-testOpt$par[1])*ci)
params=c(NA, testOpt$par[1:2], 0.25, muXi, testOpt$par[3:7])

# now get the predictions
outAreal = preds(params, nsim=2, tvec=tvec)
outPoint = predsPoint(params, nsim=2, dStar=26000)

## make plots
par(mfrow=c(1,2))
# areal
plotFault(csz, outAreal$meanSlip)
map("world", "Canada", add=TRUE)
US(add=TRUE)
# pointwise
quilt.plot(outPoint$lonLatGrid, outPoint$meanSlip)
map("world", "Canada", add=TRUE)
US(add=TRUE)

##### plot different spline basis
# plot original spline basis
# latRange=c(40,50)
latRange = c(42.5, 47.5)
lats = seq(latRange[1], latRange[2], l=500)
nKnots=5
splineMat = getSplineBasis(lats=lats, nKnots=nKnots, latRange=latRange)
matplot(lats, splineMat, type="l", lwd=2, xlim=latRange)
abline(v=dr1$Lat, col="red")

# custom spline basis
cutoffs=c(43.25, 46.75)
latRange=c(40,50)
lats = seq(latRange[1], latRange[2], l=100)
nKnots=4
splineMat = getCustomSplineBasis(lats=lats, cutoffs=cutoffs)
matplot(lats, splineMat, xlim=range(dr1$Lat))
abline(v=dr1$Lat, col="red")

# plot normal basis
# latRange=c(40,50)
latRange=range(csz$latitude)
lats = seq(latRange[1], latRange[2], l=500)
nKnots=3
sds=1.5
# splineMat = getNormalBasis(lats=lats, nKnots=nKnots, sds=sds, latRange=latRange, intercept=FALSE)
splineMat = getNormalBasis(lats=lats, nKnots=nKnots, latRange=latRange, intercept=FALSE)
matplot(lats, splineMat, type="l", lwd=2, xlim=latRange)
abline(v=dr1$Lat, col="red")

##### Infer events
# number of observations per event:
#  T1   T10 T10R1   T11   T12    T2    T3   T3a    T4   T4a    T5   T5a   T5b    T6   T6a    T7   T7a    T8   T8a    T9   T9a 
# 196     7     6     3     2    24    25    13    42    13    51    27     6    30    15    27     9    11     7     8     1 

opt = list()
opt11 = list()
opt$par = c(2.338752, 2.04888, 12.01887, 6.779287, 17.78633, -22.77749, -14.14914)
opt11$par = c(3.102128, 1.748229, 33.10879, -60.44, 51.4165, -68.62971, -22.39055)
dStar=26000
nKnots=5

# get the full set of parameters in the required order and the taper vector
testOpt = opt11
splinePar = testOpt$par[3:7]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))

# estimate muXi MLE with inverse variance weighting
logX = log(slipDatCSZ$slip)
ci = 1/sigmaXi^2
ci = ci/sum(ci)
muXi = sum((logX-testOpt$par[1])*ci)
params=c(NA, testOpt$par[1:2], 0.25, muXi, testOpt$par[3:7])

# T1
isT1 = events=="T1"
T1Dat = dr1[isT1,]
GT1 = G[isT1,]

eventPreds = predsGivenSubsidence(params, subDat=T1Dat, niter=2000, G=GT1, prior=FALSE, tvec=tvec)

# areal values of zeta
muAreal = eventPreds$zetaEsts * tvec
sdAreal = eventPreds$zetaSD * tvec
medAreal = eventPreds$zetaMed * tvec
l95Areal = eventPreds$zeta025 * tvec
u95Areal = eventPreds$zeta975 * tvec

# get simulations
tab <- extract(eventPreds$predResults, permuted = TRUE)
zetaSims = t(tab$zeta)
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# # recompute tvec for GPS data
# xs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
# gpsDat = data.frame(list(latitude=slipDatCSZ$lat, depth=slipDatCSZ$Depth))
# tvec = getTaperSpline(splinePar, fault=gpsDat, dStar=dStar)
# 
# # pointwise values of zeta
# muPoint = eventPreds$zetaPointEsts * tvec
# sdPoint = eventPreds$zetaPointSD * tvec
# medPoint = eventPreds$zetaPointMed * tvec
# l95Point = eventPreds$zetaPoint025 * tvec
# u95Point = eventPreds$zetaPoint975 * tvec

# plot results:
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1")
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=TRUE, 
              fileNameRoot="T1Log")
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", fileNameRoot="T1", subDat=T1Dat)
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", fileNameRoot="T1Log", 
                   subDat=T1Dat, logScale=TRUE)
# Note the marginal SD of zeta (pre-tapaer) is approximately 680: 
# 
# mean(sqrt((exp(diag(Sigma)) - 1)*exp(2*muZeta + diag(Sigma))))

# look at 4 spline basis knots optimum:
load("splineFit26kGrad4knots.RData")
## first get the parameters
testOpt = out
params=testOpt$MLEs
splinePar = params[6:length(params)]
tvec=getTaperSpline(splinePar, nKnots=4, dStar=26000)

comparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="4 Knots ", fileNameRoot="4knots", subDat=T1Dat)
comparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="4 Knots", fileNameRoot="4knotsLog", 
                   subDat=dr1, logScale=TRUE)

# generate T1 predictions

# plot results:
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1")
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=TRUE, 
              fileNameRoot="T1Log")
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", fileNameRoot="T1", subDat=T1Dat)
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", fileNameRoot="T1Log", 
                   subDat=T1Dat, logScale=TRUE)

##### plot redone constant lambda fit
load("fixedFit_MVN.RData")
testOpt = fixedFitMVN
params=testOpt$MLEs
splinePar = params[6:length(params)]
tvec=getTaperSpline(splinePar, nKnots=1, dStar=26000)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)

comparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="Consant Fit", fileNameRoot="const", subDat=dr1)
comparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="Constant Fit", fileNameRoot="constLog", 
                   subDat=dr1, logScale=TRUE)

##### plot additional points

library(maps)
library(mapproj)
library(mapdata)
library(rgeos)
library(maptools)
library(shapefiles)
library(ggmap)

westCoast <- readShapeLines('~/Desktop/Western/Western',
                            repair = T, delete_null_obj=TRUE)
map("world",
    xlim = c(-128.5, -120), 
    ylim = c(40, 50), 
    fill = T,
    col = 'grey',
    resolution = 0,
    bg = 'white',
    mar = c(6,3,2,1),
    type = 'n'
)

plot(westCoast, add = T)
box()

# most shapes are tiny islands we don't care about.  remove those
# > .025: 1512 polygons
# > .05: 733 polygons
# > .075: 442
# > .1: 287 polygons
goodShapes = which(westCoast$Shape_Leng > .1)
testCoast = westCoast
newData = list()
for(i in 1:14) {
  # testCoast[[i]] = testCoast[[i]][goodShapes]
  tmp = testCoast@data[[i]]
  newData[[i]] = tmp[goodShapes]
}
newData = data.frame(newData)
names(newData) = names(westCoast)
testCoast@data = newData
names(testCoast) = names(westCoast)
testCoast@lines = testCoast@lines[goodShapes]
for(i in 1:length(goodShapes)) {
  testCoast@lines[[i]]@ID = as.character(i)
}

# now try plotting
map("world",
    xlim = c(-128.5, -120), 
    ylim = c(40, 50), 
    fill = T,
    col = 'grey',
    resolution = 0,
    bg = 'white',
    mar = c(6,3,2,1),
    type = 'n'
)
plot(testCoast, add = T)
box()

# now lower resolution
# 100000: 110.5 Mb
# 1: same
# .00001: 127.5 Mb
lowRes = gSimplify(westCoast, tol=1)
westCoastShape = read.shapefile("~/Desktop/Western/Western")
newRes = dp(westCoastShape, tol=1)
write.shapefile(newRes, "lowRes")

map("world",
    xlim = c(-123.963, -123.923), 
    ylim = c(45.20, 45.238), 
    fill = T,
    col = 'grey',
    resolution = 0,
    bg = 'white',
    mar = c(6,3,2,1),
    type = 'n'
)
plot(lowRes, add = T)
box()

sizes <- sapply(ls(), function(n) object.size(get(n)), simplify = FALSE)
print(sapply(sizes[order(as.integer(sizes))], function(s) format(s, unit = 'auto')))

plot(dr1$Lon, dr1$Lat, pch="+", col="red", ylim=c(40, 50), xlim=c(-128, -123.5))
map("world", "Canada", add=TRUE, resolution=0)
US(add=TRUE)
plotFault(csz, new=FALSE, plotData=FALSE)

fauxObs = matrix(c(
  -126.998236, 49.818024, 
  -126.658771, 49.580435, 
  -126.395718, 49.411433, 
  -126.144274, 49.273292, 
  -125.541669, 48.920783, 
  -125.116728, 48.731546, 
  -124.698432, 48.590345, 
  -124.726965, 48.386018, 
  -124.717353, 48.143782, 
  -124.637417, 47.909926, 
  -124.444246, 47.754555, 
  -124.362284, 47.564145, 
  -124.413923, 42.510306, 
  -124.363353, 42.185168, 
  -124.211893, 41.873892, 
  -124.080846, 41.547055, 
  -124.098906, 41.242677, 
  -124.363634, 40.275861, 
  -124.102275, 40.093722
), byrow = TRUE, ncol=2) 


points(fauxObs[,1], fauxObs[,2], pch="+", col="blue")

# testMap = get_openstreetmap(bbox=c(left=-128.5, bottom=40, right=-122, left=-128.5), 
# format="pdf", color="bw")

library(maps)
library(maptools)  ## For map2SpatialPolygons()
library(mapdata)

## Convert data from a "maps" object to a "SpatialPolygonsDataFrame" object
mp <- map("state", fill = TRUE, plot = FALSE)
SP <- map2SpatialPolygons(mp, IDs = mp$names, 
                          proj4string = CRS("+proj=longlat +datum=WGS84"))
DATA <- data.frame(seq_len(length(SP)), row.names = names(SP))
SPDF <- SpatialPolygonsDataFrame(SP, data = DATA)

## Plot it
spplot(SPDF, col.regions = "transparent", colorkey = FALSE,
       par.settings = list(axis.line = list(col = "transparent")))


mp <- map("world", fill = TRUE, plot = FALSE)
SP <- map2SpatialPolygons(mp, IDs = mp$names, 
                          proj4string = CRS("+proj=longlat +datum=WGS84"))
DATA <- data.frame(seq_len(length(SP)), row.names = names(SP))
SPDF <- SpatialPolygonsDataFrame(SP, data = DATA)
mp2 <- map("Canada", fill = TRUE, plot = FALSE)
SP2 <- map2SpatialPolygons(mp, IDs = mp$names, 
                           proj4string = CRS("+proj=longlat +datum=WGS84"))
DATA2 <- data.frame(seq_len(length(SP)), row.names = names(SP))
SPDF2 <- SpatialPolygonsDataFrame(SP, data = DATA)

## Plot it
spplot(SPDF, col.regions = "transparent", colorkey = FALSE,
       par.settings = list(axis.line = list(col = "transparent")), 
       xlim=c(-128.5, -122), ylim=c(40, 50))

library(rworldmap)
library(rworldxtra)
newmap <- getMap(resolution = "high")
plot(newmap, xlim=c(-128.5, -122), ylim=c(40, 50), asp = 1)

westCoast <- readShapeLines('~/Desktop/Western/Western',
                            repair = T)
map("world",
    xlim=c(-128.5, -122), ylim=c(40, 50), 
    fill = T,
    col = 'grey',
    resolution = 0,
    bg = 'white',
    mar = c(6,3,2,1),
    type = 'n'
)

plot(westCoast, add = T)
box()

library(rgeos)
lowRes = gSimplify(westCoast, tol=10000)
map("world",
    xlim=c(-128.5, -122), ylim=c(40, 50), 
    fill = T,
    col = 'grey',
    resolution = 0,
    bg = 'white',
    mar = c(6,3,2,1),
    type = 'n'
)

plot(lowRes, add = T)
box()

# try using google maps (ggmap does not support different projections)
library(ggmap)
library(ggplot2)

### Set a range
lon = c(-128.5, -121.5)
lat = c(40, 50)

### Get a map
# https://mapstyle.withgoogle.com/
# map <- get_map(location = c(lon[1], lat[1], lon[2], lat[2]), zoom = 6,
#                maptype = "terrain", source = "google")
style1=c(feature="administrative", element="labels", visibility="off")
style2=c("&style=", feature="road", element="geometry", visibility="off")
style3=c("&style=", feature="poi", element="labels", visibility="off")
style4=c("&style=", feature="landscape", element="labels", visibility="off")
style5=c("&style=", feature="administrative", element="geometry.stroke", color="black")
style6=c("&style=", feature="administrative", element="geometry.stroke", weight=.75)
map <- get_googlemap(center=c(lon = mean(lon), lat = mean(lat)), zoom=5,
                     style=c(style1, style2, style3, style4, style5, style6))
# map <- get_googlemap(center=c(lon = mean(lon), lat = mean(lat)), zoom=6,
# style='feature:administrative|element:labels|visibility:off&style=feature:road|element:labels|visibility:off')

### When you draw a figure, you limit lon and lat.
fauxObs = data.frame(list(Lon=fauxObs[,1], Lat=fauxObs[,2]))
mapPoints <- ggmap(map) +
  geom_point(aes(x = Lon, y = Lat), data=dr1, size=3, shape="+", color="red")
foo <- mapPoints + 
  scale_x_continuous(limits = lon, expand = c(0, 0)) +
  scale_y_continuous(limits = lat, expand = c(0, 0)) + 
  geom_point(aes(x=Lon, y=Lat), data=fauxObs, size=3, shape="+", color="blue") + 
  ggtitle("Real and `Faux' Observations") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Longitude", y="Latitude")

foo

##### try to fit model with prior
nKnots=5
dStar=21000
constrLambda=FALSE
useSubPrior=TRUE
useSlipPrior=TRUE
splineFit21k5 = doFitSpline(initParams=c(MLEs[2], MLEs[3], MLEs[6], rep(0, nKnots-1)), dStar=dStar, useMVNApprox=TRUE, 
                            useGrad=TRUE, nKnots=nKnots, maxit=100, useSubPrior=useSubPrior, useSlipPrior=useSlipPrior, 
                            G=G, constrLambda=constrLambda)

# keep optimization in case we're not done yet
init = c(splineFit21k15$MLEs[2], splineFit21k15$MLEs[3], splineFit21k15$splineParMLE)
init = c(opt2$MLEs[2], opt2$MLEs[3], opt2$splineParMLE)
opt3= doFitSpline(initParams=init, dStar=dStar, useMVNApprox=TRUE, 
                  useGrad=TRUE, nKnots=nKnots, maxit=100, useSubPrior=TRUE, useSlipPrior=TRUE, G=G)

load("splineFit21k25.RData")
splineFit21k15 = opt2
params = splineFit21k5$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "21k5FrechetNoConstr")

splineFit21k5 = splineFit21k
save(splineFit21k25, file="splineFit21k25.RData")

# now generate predictions
load("splineFit21k25.RData")
params = splineFit21k15$MLEs
splinePar = params[6:length(params)]
par(mfrow=c(1,1))
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="15 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "21k15Prior")

##### try out different error inflation levels (multiply SEs by factor)
inflateDR1 = dr1
# inflates = rev(c(.5, .75, 1, 1.25, 1.5, 1.75, 2))
# inflates = rev(c(1, 1.25, 1.5, 1.75, 2))
inflates=2
startQS = c()
startQY = c()
endQS = c()
endQY = c()
LLs = c()
for(i in 1:length(inflates)) {
  inflate = inflates[i]
  inflateDR1 = dr1
  inflateDR1$Uncertainty = dr1$Uncertainty * inflate
  
  nKnots=5
  dStar=21000
  initPar=c(MLEs[2], MLEs[3], MLEs[6], rep(0, nKnots-1))
  startQS = c(startQS, getQS(initPar))
  statQY = c(startQY, getQY(initPar, getSplineBasis(nKnots=nKnots), G, fauxG, subDat=inflateDR1))
  splineFit21k5 = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                              useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=TRUE, useSlipPrior=TRUE, G=G, 
                              constrLambda=FALSE, subDat=inflateDR1)
  endPar = splineFit21k5$MLEs
  endPar = c(endPar[2], endPar[3], endPar[6:length(endPar)])
  endQS = c(endQS, getQS(endPar))
  endQY = c(endQY, getQY(endPar, getSplineBasis(nKnots=nKnots), G, fauxG, subDat=inflateDR1))
  LLs = c(LLs, splineFit21k5$logLikMLE)
  
  # save results
  save(splineFit21k5, file=paste0("21k5NoPriorUnconstrInflateBy", inflate, ".RData"))
  
  ##### plot marginal distribution
  params = splineFit21k5$MLEs
  splinePar = params[6:length(params)]
  tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
  comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, 
                     fileNameRoot = paste0("21k5PriorInflateBy", inflate), subDat=inflateDR1)
  
  ##### plot predictive distribution
  
  # T1
  isT1 = events=="T1"
  T1Dat = inflateDR1[isT1,]
  GT1 = G[isT1,]
  
  eventPreds = predsGivenSubsidence(params, subDat=T1Dat, niter=500, G=GT1, prior=FALSE, tvec=tvec)
  save(eventPreds, file=paste0("21k5PriorUnconstrInflateBy", inflate, "EventPreds.RData"))
  
  # areal values of zeta
  muAreal = eventPreds$zetaEsts * tvec
  sdAreal = eventPreds$zetaSD * tvec
  medAreal = eventPreds$zetaMed * tvec
  l95Areal = eventPreds$zeta025 * tvec
  u95Areal = eventPreds$zeta975 * tvec
  
  # get simulations
  tab <- extract(eventPreds$predResults, permuted = TRUE)
  zetaSims = t(tab$zeta)
  slipSims = sweep(zetaSims, 1, tvec, "*")
  slipPreds = list(meanSlip=muAreal, slipSims=slipSims)
  
  # plot results:
  plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=TRUE, 
                fileNameRoot=paste0("21k5PriorInflateBy", inflate))
  comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                     subDat=T1Dat, logScale=TRUE, fileNameRoot=paste0("21k5PriorInflateBy", inflate, "PredT1"))
}
save(inflates=inflates, startQS=startQS, startQY=startQY, endQS=endQS, endQY=endQY, LLs=LLs, 
     file="inflateResults.RData")


####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
####################################################################################################
# test including prior in the prediction but not the original parameter estimation
nKnots=5
dStar=21000
initPar=c(MLEs[2], MLEs[3], MLEs[6], rep(0, nKnots-1))
niter=100
splineFit21k5 = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                            useGrad=TRUE, nKnots=nKnots, maxit=niter, useSubPrior=FALSE, useSlipPrior=FALSE, G=G, 
                            constrLambda=FALSE, subDat=dr1)
endPar = splineFit21k5$MLEs
endPar = c(endPar[2], endPar[3], endPar[6:length(endPar)])
endQS = c(endQS, getQS(endPar))
endQY = c(endQY, getQY(endPar, getSplineBasis(nKnots=nKnots), G, fauxG, subDat=dr1))

# save results
save(splineFit21k5, file=paste0("21k5NoPrior", niter, "iter.RData"))

##### plot marginal distribution
params = splineFit21k5$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
# comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, 
#                    fileNameRoot = paste0("21k5PriorInflateBy", inflate), subDat=dr1)

##### plot predictive distribution
# T1
# eventPredsGumbel = predsGivenSubsidence2(params, subDat=dr1, niter=30, G=G, fauxG=fauxG, tvec=tvec, ev="T1")
# eventPredsFrechet = predsGivenSubsidence2(params, subDat=dr1, niter=50, G=G, fauxG=fauxG, tvec=tvec, ev="T1", FrechetFGumbelT=FALSE)
# eventPredsThresh = predsGivenSubsidence3(params, subDat=dr1, niter=400, G=G, fauxG=fauxG, tvec=tvec, ev="T1", 
#                                          slipThresh=80, adaptDelta=.8)
eventPredsThresh = predsGivenSubsidence3(params, subDat=dr1, niter=100, G=G, fauxG=fauxG, tvec=tvec, ev="T1", 
                                         slipThresh=80, alg="NUTS", adaptDelta=.999999)
eventPreds=eventPredsThresh
pairs(eventPreds$predResults, pars = c("logZetaAreal[50]", "logZetaAreal[100]", "lp__"), las = 1)
save(eventPreds, file="21k5PredPriorT1Preds.RData")

# areal values of zeta
muAreal = eventPreds$zetaEsts * tvec
sdAreal = eventPreds$zetaSD * tvec
medAreal = eventPreds$zetaMed * tvec
l95Areal = eventPreds$zeta025 * tvec
u95Areal = eventPreds$zeta975 * tvec

# get simulations
tab <- extract(eventPreds$predResults, permuted = TRUE)
zetaSims = t(tab$zeta)
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)
LLs = tab$lp__
maxSlips = tab$maxSlip
hist(maxSlips)

# plot results:
isT1 = events=="T1"
GT1 = G[isT1,]
T1Dat = dr1[isT1,]
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=TRUE, 
              fileNameRoot=paste0("21k5PriorInflateBy", inflate))
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                   subDat=T1Dat, logScale=TRUE, fileNameRoot="21k5PredThreshSlipT1.RData")

library(rstan)
funnel <- stan_demo("funnel", seed = 12345)   # has 5 divergent transitions
pairs(funnel, pars = c("y", "x[1]", "lp__"), las = 1) # below the diagonal
funnel_reparam <- stan_demo("funnel_reparam") # has no divergent transitions

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
# get subset of observations along coast for computing max subsidence.  Use a greedy algorithm:
# namely, at each iteration consider a pair (x,y) s.t. d(x,y) is minimized.  Remove either x or 
# y, and repeat until min_{x,y} d(x,y) > tol
fauxDat = getFauxDat()
tmp = data.frame(list(Lon=dr1$Lon, Lat=dr1$Lat, Uncertainty=dr1$Uncertainty, event=dr1$event, 
                      subsidence=dr1$subsidence))
fullDat = rbind(tmp, fauxDat)
text(fullDat$Lon, fullDat$Lat, labels=1:nrow(fullDat))

greedyThin = function(allDat, tol=.1) {
  coords = cbind(allDat$Lon, allDat$Lat)
  distMat = rdist(coords, coords) + diag(2*tol, nrow=nrow(coords))
  rowMat = row(distMat)
  nRM = 0
  
  # get minimum distance and index
  minI = which.min(distMat)
  minDist = distMat[minI]
  
  # remove observations
  while(minDist < tol) {
    print(paste0("minDist: ", minDist))
    print(paste0("num removed: ", nRM))
    # get minimum distance and index
    minI = which.min(distMat)
    minDist = distMat[minI]
    
    # find what observation it is and remove from set
    rowI = rowMat[minI]
    rowMat = rowMat[-rowI, -rowI]
    distMat = distMat[-rowI, -rowI]
    nRM = nRM+1
    
    # get new minimum distance and index
    minI = which.min(distMat)
    minDist = distMat[minI]
  }
  
  # find which points we're left with and plot them
  includes = rowMat[1:nrow(rowMat),]
  includeDat = allDat[includes,]
  par(mfrow=c(1,2))
  plot(allDat$Lon, allDat$Lat, pch="+")
  plot(includeDat$Lon, includeDat$Lat, pch="+")
  
  return(list(includes=includes, includeDat=includeDat))
}

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##### test what predictions look like if I cut off gaps in data:
hiLatDat = 46.75
lowLatDat = 43.25
inRangeDat = (dr1$Lat < hiLatDat) & (dr1$Lat > lowLatDat)
fauxDat = getFauxDat()
inRangeFauxDat = (fauxDat$Lat < hiLatDat) & (fauxDat$Lat > lowLatDat)
rangeDat = dr1[inRangeDat,]
rangeFauxDat = fauxDat[inRangeFauxDat,]
fauxObs = getFauxObs()
rangeFauxObs = fauxObs[inRangeFauxDat,]

# plot data in range
plot(dr1$Lon[!inRangeDat], dr1$Lat[!inRangeDat], col="red")
points(rangeDat$Lon, rangeDat$Lat)
map("world", "Canada", add=TRUE)
US(add=TRUE)

# subset fault to be within fault
hiLatFault = 47.5
lowLatFault = 42.5
faultLatRange=c(lowLatFault, hiLatFault)
inRangeFault = (csz$latitude < hiLatFault) & (csz$latitude > lowLatFault)
faultRange = csz[inRangeFault,]
GRange = G[inRangeDat, inRangeFault]
fauxGRange = fauxG[inRangeFauxDat, inRangeFault]

# plot fault in range
plotFault(faultRange, plotDat=FALSE, new=FALSE)

# inflate data if necessary
rangeDat$Uncertainty = rangeDat$Uncertainty*2

# now fit model and generate predictions
nKnots=3
dStar=21000
initPar=c(MLEs[2], MLEs[3], MLEs[6], rep(0, nKnots-1))
splineFit21k5Inflate = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                                   useGrad=TRUE, nKnots=nKnots, maxit=200, useSubPrior=FALSE, useSlipPrior=FALSE, G=GRange, 
                                   fauxG=fauxGRange, constrLambda=FALSE, subDat=rangeDat, fault=faultRange, latRange=faultLatRange)
endPar = splineFit21k5Inflate$MLEs
endPar = c(endPar[2], endPar[3], endPar[6:length(endPar)])
endQS = getQS(endPar, nKnots=nKnots, fault=faultRange)
endQY = getQY(endPar, getSplineBasis(nKnots=nKnots, fault=faultRange, latRange=latRangeFault), 
              GRange, fault=faultRange, fauxGRange, subDat=rangeDat, fauxObs=rangeFauxObs)

# save results
save(splineFit21k5Inflate, file=paste0("21k5NoPriorRangedInflate.RData"))

##### plot marginal distribution
params = splineFit21k5Inflate$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=faultRange, latRange=faultLatRange)
comparePredsToSubs(params, G=GRange, plotNameRoot="3 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, 
                   fileNameRoot=paste0("21k5PriorRangeInflateby2"), subDat=rangeDat, fault=faultRange, latRange=faultLatRange)

##### plot predictive distribution

# T1
isT1 = events=="T1"
T1DatRange = dr1[isT1 & inRangeDat,]
GT1Range = G[isT1 & inRangeDat, inRangeFault]

eventPreds = predsGivenSubsidence(params, fault=faultRange, subDat=T1DatRange, niter=500, G=GT1Range, prior=FALSE, tvec=tvec)
save(eventPreds, file=paste0("21k5PriorRangeEventPreds.RData"))

# areal values of zeta
muAreal = eventPreds$zetaEsts * tvec
sdAreal = eventPreds$zetaSD * tvec
medAreal = eventPreds$zetaMed * tvec
l95Areal = eventPreds$zeta025 * tvec
u95Areal = eventPreds$zeta975 * tvec

# get simulations
tab <- extract(eventPreds$predResults, permuted = TRUE)
zetaSims = t(tab$zeta)
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=TRUE, 
              fileNameRoot=paste0("21k5PriorRange"))
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1Range, tvec=tvec, plotNameRoot="T1", 
                   subDat=T1DatRange, logScale=TRUE, fileNameRoot=paste0("21k5PriorRangePredT1"), 
                   fault=faultRange, latRange=faultLatRange)

# 4/13 meeting
par(mfrow=c(1,2))
hiLatDat = 46.75
lowLatDat = 43.25
inRangeDat = (dr1$Lat < hiLatDat) & (dr1$Lat > lowLatDat)
fauxDat = getFauxDat()
inRangeFauxDat = (fauxDat$Lat < hiLatDat) & (fauxDat$Lat > lowLatDat)
rangeDat = dr1[inRangeDat,]
rangeFauxDat = fauxDat[inRangeFauxDat,]
fauxObs = getFauxObs()
rangeFauxObs = fauxObs[inRangeFauxDat,]

# plot data in range
plot(dr1$Lon[!inRangeDat], dr1$Lat[!inRangeDat], col="red", xlab="Longitude", ylab="Latitude", pch="+", main="Observations Included in Test")
points(rangeDat$Lon, rangeDat$Lat, pch="+")
map("world", "Canada", add=TRUE)
US(add=TRUE)

# subset fault to be within fault
hiLatFault = 47.5
lowLatFault = 42.5
faultLatRange=c(lowLatFault, hiLatFault)
inRangeFault = (csz$latitude < hiLatFault) & (csz$latitude > lowLatFault)
faultRange = csz[inRangeFault,]
GRange = G[inRangeDat, inRangeFault]
fauxGRange = fauxG[inRangeFauxDat, inRangeFault]

# plot fault in range
plotFault(faultRange, plotDat=FALSE, new=FALSE)

# plot spline basis
nKnots=3
splineMat = getNormalBasis(lats=lats, nKnots=nKnots, latRange=faultLatRange, intercept=FALSE)
matplot(lats, splineMat, type="l", lwd=2, xlim=latRange, main="Gaussian Spline Basis", xlab="Latitude")




##### test estimators of gps data scaling factor under normal model

getGPSLik = function(muZeta, sigmaZeta, muXi, corPar=NULL, normalModel=FALSE) {
  # get correlation parameters if necessary
  if(is.null(corPar)) {
    corPar = getCorPar(normalModel = normalModel)
    phiZeta = corPar$phiZeta
    nuZeta = corPar$nuZeta
  }
  
  # compute correlation matrix for GPS data (or log GPS data if lognormal model)
  coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                             onlyUpper=TRUE, smoothness=nuZeta, 
                             Distance="rdist.earth", Dist.args=list(miles=FALSE))
  SigmaZeta = sigmaZeta^2 * corMatGPS
  
  # get measurement error and observations for the model
  if(!normalModel) {
    sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
    x = log(gpsDat$slip)
    xCntr = x - muXi - muZeta
  }
  else {
    sigmaXi = gpsDat$slipErr
    x = gpsDat$slip
    xCntr = x - muXi*muZeta
    
    # we also need to modify variance accounting for scaling factor muXi
    SigmaZeta = muXi^2 * SigmaZeta
  }
  
  # get likelihood
  Sigma = SigmaZeta + diag(sigmaXi^2)
  return(logLikGP(xCntr, chol(Sigma)))
}


# first test the previous estimator transformed to the current scale (This one seems to work the best)
gpsDat = slipDatCSZ
sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
muZeta= MLEs[2]
x = log(gpsDat$slip)
ci = 1/sigmaXi^2
ci = ci/sum(ci)
muXi = sum((x-muZeta)*ci)
# muZeta = exp(MLEs[2] + MLEs[3]^2/2) # MLEs[3] is just too big to be reasonable here
muZeta = exp(MLEs[2])
# sigmaZeta = sqrt((exp(MLEs[3]^2) - 1)*exp(2*MLEs[2] + MLEs[3]^2))
# sigmaZeta = MLEs[3]
sigmaZeta = sd(gpsDat$slip)
# sigmaZeta = sqrt(var(gpsDat$slip) - mean(gpsDat$slipErr^2))
scaleFac = exp(muXi)
scaledDat = gpsDat$slip * scaleFac
hist(scaledDat)
getGPSLik(muZeta, sigmaZeta, scaleFac, normalModel=TRUE)

# now test a mean-based estimator
gpsDat = slipDatCSZ
sigmaXi = gpsDat$slipErr
# muZeta = exp(MLEs[2] + MLEs[3]^2/2)
muZeta = exp(MLEs[2])
# sigmaZeta = sqrt((exp(MLEs[3]^2) - 1)*exp(2*MLEs[2] + MLEs[3]^2))
sigmaZeta = sd(gpsDat$slip)
# sigmaZeta = sqrt(var(gpsDat$slip) - mean(gpsDat$slipErr^2))
scaleFac = mean(gpsDat$slip)/muZeta
scaledDat = gpsDat$slip * scaleFac
hist(scaledDat)
getGPSLik(muZeta, sigmaZeta, scaleFac, normalModel=TRUE)

# now test a WLS mean-based estimator
gpsDat = slipDatCSZ
sigmaXi = gpsDat$slipErr
# muZeta = exp(MLEs[2] + MLEs[3]^2/2)
muZeta = exp(MLEs[2])
# sigmaZeta = sqrt((exp(MLEs[3]^2) - 1)*exp(2*MLEs[2] + MLEs[3]^2))
sigmaZeta = sd(gpsDat$slip)
# sigmaZeta = sqrt(var(gpsDat$slip) - mean(gpsDat$slipErr^2))
ci = 1/gpsDat$slipErr^2
ci = ci/sum(ci)
scaleFac = sum(gpsDat$slip*ci)/muZeta
scaledDat = gpsDat$slip * scaleFac
hist(scaledDat)
getGPSLik(muZeta, sigmaZeta, scaleFac, normalModel=TRUE)

# now test a GLS mean-based estimator (that ignores noise in X)
gpsDat = slipDatCSZ
sigmaXi = gpsDat$slipErr
corPar = getCorPar(normalModel=TRUE)
phiZeta = corPar$phiZeta
nuZeta = corPar$nuZeta
coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                           onlyUpper=TRUE, smoothness=nuZeta, 
                           Distance="rdist.earth", Dist.args=list(miles=FALSE))
# muZeta = exp(MLEs[2] + MLEs[3]^2/2)
muZeta = exp(MLEs[2])
# sigmaZeta = sqrt((exp(MLEs[3]^2) - 1)*exp(2*MLEs[2] + MLEs[3]^2))
sigmaZeta = sd(gpsDat$slip)
# sigmaZeta = sqrt(var(gpsDat$slip) - mean(gpsDat$slipErr^2))
xmat = cbind(rep(1, nrow(coords)))
# meanEst = (t(xmat) %*% solve(corMatGPS) %*% xmat)^(-1) * t(xmat) %*% solve(corMatGPS) %*% gpsDat$slip
meanEst = sum(solve(corMatGPS, xmat))^(-1) * sum(solve(corMatGPS, gpsDat$slip))
scaleFac = meanEst/muZeta
scaledDat = gpsDat$slip * scaleFac
hist(scaledDat)
getGPSLik(muZeta, sigmaZeta, scaleFac, normalModel=TRUE)

# now test a GLS mean-based estimator setting marginal variances to measurement noise
gpsDat = slipDatCSZ
sigmaXi = gpsDat$slipErr
corPar = getCorPar(normalModel=TRUE)
phiZeta = corPar$phiZeta
nuZeta = corPar$nuZeta
coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                           onlyUpper=TRUE, smoothness=nuZeta, 
                           Distance="rdist.earth", Dist.args=list(miles=FALSE))
# muZeta = exp(MLEs[2] + MLEs[3]^2/2)
muZeta = exp(MLEs[2])
# sigmaZeta = sqrt((exp(MLEs[3]^2) - 1)*exp(2*MLEs[2] + MLEs[3]^2))
sigmaZeta = sd(gpsDat$slip)
# sigmaZeta = sqrt(var(gpsDat$slip) - mean(gpsDat$slipErr^2))
covMatGPS = sweep(sweep(corMatGPS, 1, sigmaXi, "*"), 2, sigmaXi, "*")
xmat = cbind(rep(1, nrow(coords)))
meanEst = sum(solve(covMatGPS, xmat))^(-1) * sum(solve(covMatGPS, gpsDat$slip))
scaleFac = meanEst/muZeta
scaledDat = gpsDat$slip * scaleFac
hist(scaledDat)
getGPSLik(muZeta, sigmaZeta, scaleFac, normalModel=TRUE)

## test the closed form estimator assuming X = C*(zeta + xi)
muZeta = exp(MLEs[2])
# sigmaZeta = sd(gpsDat$slip)
sigmaZeta = sqrt(var(gpsDat$slip) - mean(gpsDat$slipErr^2))
# first compute V, the variance matrix of zeta plus the measurement error
coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
corPar = getCorPar(normalModel = TRUE)
phiZeta = corPar$phiZeta
nuZeta = corPar$nuZeta
corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                           onlyUpper=TRUE, smoothness=nuZeta, 
                           Distance="rdist.earth", Dist.args=list(miles=FALSE))
SigmaZeta = sigmaZeta^2 * corMatGPS
V = SigmaZeta + diag(slipDatCSZ$slipErr^2)

# now compute estimator
VInv = solve(V)
oneVInv = apply(VInv, 1, sum)
p = ncol(V)
muzp = muZeta
x = slipDatCSZ$slip
muXi = c((muzp * oneVInv %*% x + sqrt((muzp^2 * oneVInv %*% x)^2 + 4*t(x) %*% VInv %*% x))/2)
muXi = c(-oneVInv %*% x + sqrt((muZeta * oneVInv %*% x)^2 + 4* t(x) %*% VInv %*% x))/(2)
getGPSLik(muZeta, sigmaZeta, muXi, normalModel=TRUE)


# now get areally averaged correlation matrix
# arealCSZCorNormal = arealZetaCov(normalModel=TRUE)
# > arealCSZCorNormal2 = arealZetaCov(normalModel=TRUE, nDown1=4, nStrike1=5)
# > mean((arealCSZCorNormal - arealCSZCorNormal2)^2)
# [1] 6.106394e-10
# > mean(diag((arealCSZCorNormal - arealCSZCorNormal2))^2)
# [1] 7.687583e-09
# > arealCSZCorNormal3 = arealZetaCov(normalModel=TRUE, nDown1=4, nStrike1=6)
# > mean((arealCSZCorNormal - arealCSZCorNormal3)^2)
# [1] 9.499436e-10
# > mean(diag((arealCSZCorNormal - arealCSZCorNormal3))^2)
# [1] 1.24706e-08
# > max((arealCSZCorNormal - arealCSZCorNormal3)^2)
# [1] 2.500993e-08
# > max((arealCSZCorNormal - arealCSZCorNormal2)^2)
# [1] 1.532497e-08
# > arealCSZCorNormal4 = arealZetaCov(normalModel=TRUE, nDown1=6, nStrike1=8)
# > max((arealCSZCorNormal - arealCSZCorNormal4)^2)
# [1] 5.538091e-08
# > max((arealCSZCorNormal3 - arealCSZCorNormal4)^2)
# [1] 5.95771e-09
arealCSZCorNormal5 = arealZetaCov(normalModel=TRUE, nDown1=9, nStrike1=12)
# > max(abs(arealCSZCorNormal - arealCSZCorNormal5))
# [1] 0.0002787949
# > max(abs(arealCSZCorNormal4 - arealCSZCorNormal5))
# [1] 4.346344e-05
# > max(abs(diag(arealCSZCorNormal4 - arealCSZCorNormal5)))
# [1] 4.346344e-05
# > max(abs(arealCSZCorNormal - arealCSZCorNormal5))
# [1] 0.0002787949
# > max(abs(diag(arealCSZCorNormal - arealCSZCorNormal5)))
# [1] 0.0002787949
arealCSZCor = arealCSZCorNormal5
save(arealCSZCor, file="arealCSZCorNormal.RData")

##### now try fitting normal model (without gradient):
##### test what predictions look like if I cut off gaps in data:
hiLatDat = 46.75
lowLatDat = 43.25
inRangeDat = (dr1$Lat < hiLatDat) & (dr1$Lat > lowLatDat)
fauxDat = getFauxDat()
inRangeFauxDat = (fauxDat$Lat < hiLatDat) & (fauxDat$Lat > lowLatDat)
rangeDat = dr1[inRangeDat,]
rangeFauxDat = fauxDat[inRangeFauxDat,]
fauxObs = getFauxObs()
rangeFauxObs = fauxObs[inRangeFauxDat,]

# plot data in range
plot(dr1$Lon[!inRangeDat], dr1$Lat[!inRangeDat], col="red")
points(rangeDat$Lon, rangeDat$Lat)
map("world", "Canada", add=TRUE)
US(add=TRUE)

# subset fault to be within fault
hiLatFault = 47.5
lowLatFault = 42.5
faultLatRange=c(lowLatFault, hiLatFault)
inRangeFault = (csz$latitude < hiLatFault) & (csz$latitude > lowLatFault)
faultRange = csz[inRangeFault,]
GRange = G[inRangeDat, inRangeFault]
fauxGRange = fauxG[inRangeFauxDat, inRangeFault]

# plot fault in range
plotFault(faultRange, plotDat=FALSE, new=FALSE)
rangeDatCopy = rangeDat

# inflate data if necessary
inflate=1.75
rangeDat$Uncertainty = rangeDatCopy$Uncertainty * inflate 
# inflate=2/3: -1994.4
# inflate=1:   -1794.9, -1791.6
# inflate=1.5: -1696.8
# inflate=7/4: -1696.3
# inflate=2:   -1706.6
# inflate=2.5: -1740.3
# inflate=3:   -1778.8


# now fit model and generate predictions
# NOOOOOOOOOOOTE: make sure to modify spline basis if necessary
nKnots=5
dStar=21000
splineFit21k5Ranged = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                                  useGrad=TRUE, nKnots=nKnots, maxit=200, useSubPrior=FALSE, useSlipPrior=FALSE, G=GRange, 
                                  fauxG=fauxGRange, constrLambda=FALSE, subDat=rangeDat, fault=faultRange, latRange=faultLatRange, 
                                  normalModel=TRUE)
endPar = splineFit21k5Ranged$MLEs

# test to see if we're at the optimum
newInit = endPar[c(2,3,(length(endPar)-nKnots+1):length(endPar))]
splineFit21k5RangedNew = doFitSpline(initParams=newInit, dStar=dStar, useMVNApprox=TRUE, 
                                     useGrad=TRUE, nKnots=nKnots, maxit=200, useSubPrior=FALSE, useSlipPrior=FALSE, G=GRange, 
                                     fauxG=fauxGRange, constrLambda=FALSE, subDat=rangeDat, fault=faultRange, latRange=faultLatRange, 
                                     normalModel=TRUE)

# generate subsidence predictions
params = splineFit21k5RangedNew$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=faultRange, latRange=faultLatRange)
plotFault(faultRange, tvec)
comparePredsToSubs(params, G=GRange, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2, 
                   fileNameRoot=paste0("21k5PriorRangeNormalGradInflate2.5"), subDat=rangeDat, fault=faultRange, 
                   latRange=faultLatRange, normalModel=TRUE)

# check to find probability all slips are positive (negative slips are all negative)
library(mvtnorm)
muZeta = params[2]
sigmaZeta = params[3]
arealCSZCor = getArealCorMat(faultRange, normalModel=TRUE)
arealCSZCov = arealCSZCor * sigmaZeta^2
pmvnorm(upper=rep(0, nrow(faultRange)), mean=rep(-muZeta, nrow(faultRange)), sigma=arealCSZCov)
# could also simulate distribution of minimum slip


# now fit model and generate predictions
nKnots=5
dStar=21000

LLsLN = c()
LLsSubLN = c()
LLsN = c()
LLsSubN = c()
pPos = c()
# inflates = c(1, seq(1.5,2.5, l=7))
inflates = seq(1, 3, l=13)
arealCSZCor = getArealCorMat(faultRange, normalModel=TRUE)
# NOOOOOOOOOOOTE: make sure to modify spline basis if necessary
for(i in 1:length(inflates)) {
  # inflate data if necessary
  inflate=inflates[i]
  rangeDat$Uncertainty = rangeDatCopy$Uncertainty * inflate
  
  # fit the normal model
  initPar = c(20,20, rep(1,nKnots))
  splineFit21k5Normal = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                                    useGrad=TRUE, nKnots=nKnots, maxit=200, useSubPrior=FALSE, useSlipPrior=FALSE, G=GRange, 
                                    fauxG=fauxGRange, constrLambda=FALSE, subDat=rangeDat, fault=faultRange, latRange=faultLatRange, 
                                    normalModel=TRUE)
  LLsN = c(LLsN, splineFit21k5Normal$logLikMLE)
  subLL = splineFit21k5Normal$optimTable[nrow(splineFit21k5Normal$optimTable), ncol(splineFit21k5Normal$optimTable)-3]
  LLsSubN = c(LLsSubN, subLL)
  
  # get the probability all slips are positive
  library(mvtnorm)
  params = splineFit21k5Normal$MLEs
  sigmaZeta = params[3]
  arealCSZCov = arealCSZCor * sigmaZeta^2
  posProb = pmvnorm(upper=rep(0, nrow(faultRange)), mean=rep(-muZeta, nrow(faultRange)), sigma=arealCSZCov)
  pPos = c(pPos, posProb)
  
  # fit the normal model
  initPar = c(1,1, rep(1,nKnots))
  splineFit21k5LN = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                                useGrad=TRUE, nKnots=nKnots, maxit=200, useSubPrior=FALSE, useSlipPrior=FALSE, G=GRange, 
                                fauxG=fauxGRange, constrLambda=FALSE, subDat=rangeDat, fault=faultRange, latRange=faultLatRange, 
                                normalModel=FALSE)
  LLsLN = c(LLsLN, splineFit21k5LN$logLikMLE)
  subLL = splineFit21k5LN$optimTable[nrow(splineFit21k5LN$optimTable), ncol(splineFit21k5LN$optimTable)-3]
  LLsSubLN = c(LLsSubLN, subLL)
}
# save and plot results
# inflateDat = list(inflates=inflates, LLsN=LLsN, pPos=pPos, LLsLN=LLsLN)
inflateDat = list(inflates=inflates, LLsN=LLsN, LLsSubN=LLsSubN, pPos=pPos, LLsLN=LLsLN, LLsSubLN=LLsSubLN)
save(inflateDat, file="inflateDat.RData")
par(mfrow=c(1,2))
plot(inflates, LLsN, type="l", main="Inflation LL (Normal)", xlab="SD Inflation", ylab="LL")
abline(v=inflates[which.max(LLsN)], col="red", lty=2)
plot(inflates, LLsLN, type="l", main="Inflation LL (Log-Normal)", xlab="SD Inflation", ylab="LL", col="blue")
abline(v=inflates[which.max(LLsLN)], col="red", lty=2)
par(mfrow=c(1,1))
subLLRange = range(c(LLsSubLN, LLsSubN))
plot(inflates, LLsSubN, type="l", main="Subsidence Log-Likelihood", xlab="SD Inflation", ylab="LL", ylim=subLLRange)
abline(v=inflates[which.max(LLsSubLN)], col="red", lty=2)
abline(v=inflates[which.max(LLsSubN)], col="red", lty=2)
lines(inflates, LLsSubLN, col="blue")
plot(inflates, pPos, type="l", main="All Positive Slip Probability (Normal)", xlab="SD Inflation", ylab="LL", col="blue")
abline(v=inflates[which.max(LLsSubN)], col="red", lty=2)

##### plot covariogram fit
par(mfrow=c(1,1))
locs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
n <- nrow(locs)
is = rep(1:n, n)
js = rep(1:n, rep(n, n))
ind <- is > js
id <- cbind(is, js)[ind, ]
distVec <- rdist.earth.vec(locs[id[, 1], ], locs[id[, 2],], miles=FALSE)
vg = vgram(locs, slipDatCSZ$slip, lon.lat=TRUE, d=distVec, type="correlogram", breaks=seq(0, max(distVec), l=60))
pdf(file="normalCorrelogram.pdf", width=5, height=5)
plot(vg, N=75)
corPar = getCorPar(TRUE)
nu = corPar$nuZeta
phi = corPar$phiZeta
lambda = corPar$lambda # nugget to sill ratio
xs = seq(0, max(distVec), l=100)
cors = Matern(xs, range=phi, smoothness = 3/2) * (1 - lambda)
lines(xs, cors, col="blue")
dev.off()


#####
# take a look at Rob McCaffrey's GPS data, compare with Pollitz/Evans data
# https://dl.dropboxusercontent.com/u/106218507/site/defnode/tdefnode_io.html#nod
# above link has structure for pn1d.nod
test = read.table("pn2d.nod")
colnames(test) = c("FaultName", "FaultNum", "Node_X_Index", "Node_Z_Index", "Hanging_Wall_Block_Name", 
                   "Foot_Wall_Block_Name", "Lon", "Lat", "Depth", "Phi", 
                   "Phi_SD", "SlipERate", "SlipNRate", "SlipERateSD", "SlipNRateSD", 
                   "SlipRateNECorr", "SlipEDefRate", "SlipNDefRate", "Slip_Azimuth", "Along_Strike_Dist", 
                   "Across_Strike_Dist", "Downdip_Dist", "Strike", "Dip", "Moment_Rate")
test$Lon = test$Lon - 360
Cascadia = test[test$FaultName == "Cascadia",]

# compute total slip deficit rate and SD:
library(shotGroups)
# slipDefRate = sqrt(Cascadia$SlipERate^2 + Cascadia$SlipNRate^2)
getSlipDefRate = function(i) {
  qmvnEll(.5, sigma=diag(c(Cascadia$SlipERateSD[i]^2, Cascadia$SlipNRateSD[i]^2)), mu=c(Cascadia$SlipERate[i], Cascadia$SlipNRate[i]))
}
slipDefRate1 = sapply(1:nrow(Cascadia), getSlipDefRate)
slipDefRate2 = sqrt(Cascadia$SlipERate^2 + Cascadia$SlipNRate^2)
slipDefRateSD = 4

# plot the fault depth
par(mfrow=c(1,1))
quilt.plot(Cascadia$Lon, Cascadia$Lat, Cascadia$Depth)
map("world", "Canada", add=TRUE)
US(add=TRUE)

# plot the fault slip rate estimates
par(mfrow=c(1,2))
quilt.plot(Cascadia$Lon, Cascadia$Lat, Cascadia$SlipEDefRate)
map("world", "Canada", add=TRUE)
US(add=TRUE)
quilt.plot(Cascadia$Lon, Cascadia$Lat, Cascadia$SlipNDefRate)
map("world", "Canada", add=TRUE)
US(add=TRUE)

# plot the fault depth
quilt.plot(Cascadia$Lon, Cascadia$Lat, Cascadia$Depth)
map("world", "Canada", add=TRUE)
US(add=TRUE)


##############################################################################################################
# plot results from the best inflated normal and lognormal model
hiLatDat = 46.75
lowLatDat = 43.25
inRangeDat = (dr1$Lat < hiLatDat) & (dr1$Lat > lowLatDat)
fauxDat = getFauxDat()
inRangeFauxDat = (fauxDat$Lat < hiLatDat) & (fauxDat$Lat > lowLatDat)
rangeDat = dr1[inRangeDat,]
rangeFauxDat = fauxDat[inRangeFauxDat,]
fauxObs = getFauxObs()
rangeFauxObs = fauxObs[inRangeFauxDat,]

# plot data in range
plot(dr1$Lon[!inRangeDat], dr1$Lat[!inRangeDat], col="red")
points(rangeDat$Lon, rangeDat$Lat)
map("world", "Canada", add=TRUE)
US(add=TRUE)

# subset fault to be within fault
hiLatFault = 47.5
lowLatFault = 42.5
faultLatRange=c(lowLatFault, hiLatFault)
inRangeFault = (csz$latitude < hiLatFault) & (csz$latitude > lowLatFault)
faultRange = csz[inRangeFault,]
GRange = G[inRangeDat, inRangeFault]
fauxGRange = fauxG[inRangeFauxDat, inRangeFault]

# plot fault in range
plotFault(faultRange, plotDat=FALSE, new=FALSE)
rangeDatCopy = rangeDat

# inflate data if necessary
inflate=1.75
rangeDat$Uncertainty = rangeDatCopy$Uncertainty * inflate 
# inflate=2/3: -1994.4
# inflate=1:   -1794.9, -1791.6
# inflate=1.5: -1696.8
# inflate=7/4: -1696.3
# inflate=2:   -1706.6
# inflate=2.5: -1740.3
# inflate=3:   -1778.8


# now fit model and generate predictions for the ranged model
# NOOOOOOOOOOOTE: make sure to modify spline basis if necessary
nKnots=5
dStar=21000
initPar=c(20,15, rep(1, nKnots))
splineFit21k5RangedN = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                                   useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=GRange, 
                                   fauxG=fauxGRange, constrLambda=FALSE, subDat=rangeDat, fault=faultRange, latRange=faultLatRange, 
                                   normalModel=TRUE)
endParN = splineFit21k5RangedN$MLEs

initPar=c(2, 1.5, rep(1, nKnots))
splineFit21k5RangedLN = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                                    useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=GRange, 
                                    fauxG=fauxGRange, constrLambda=FALSE, subDat=rangeDat, fault=faultRange, latRange=faultLatRange, 
                                    normalModel=FALSE)
endParLN = splineFit21k5RangedLN$MLEs

# generate subsidence predictions
params = splineFit21k5RangedN$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=faultRange, latRange=faultLatRange)
plotFault(faultRange, tvec)
comparePredsToSubs(params, G=GRange, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5RangeNormalInflate1.75"), subDat=rangeDat, fault=faultRange, 
                   latRange=faultLatRange, normalModel=TRUE)
comparePredsToSubs(params, G=GRange, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5RangePosNormalInflate1.75"), subDat=rangeDat, fault=faultRange, 
                   latRange=faultLatRange, normalModel=TRUE, posNormalModel=TRUE)
params = splineFit21k5RangedLN$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=faultRange, latRange=faultLatRange)
plotFault(faultRange, tvec)
comparePredsToSubs(params, G=GRange, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5RangeLogNormalInflate1.75"), subDat=rangeDat, fault=faultRange, 
                   latRange=faultLatRange, normalModel=FALSE)

# check to find probability all slips are positive (negative slips are all negative)
library(mvtnorm)
params = splineFit21k5N$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=faultRange, latRange=faultLatRange)
arealCSZCor = getArealCorMat(faultRange, normalModel=TRUE)
arealCSZCov = arealCSZCor * sigmaZeta^2
arealSlipCov = sweep(sweep(arealCSZCov, 1, tvec, "*"), 2, tvec, "*")
meanSlips = tvec * muZeta
par(mfrow=c(1,1))
plotMinMVN(meanSlips, arealSlipCov, fault=faultRange)

##### do the same thing for the full fault model
# NOOOOOOOOOOOTE: make sure to modify spline basis if necessary
inflate=1.75
inflateDr1 = dr1
inflateDr1$Uncertainty = dr1$Uncertainty*inflate
nKnots=5
dStar=21000
initPar=c(20,15, 1, rep(0, nKnots-1))
splineFit21k5N = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                             useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=G, 
                             fauxG=fauxG, constrLambda=FALSE, subDat=inflateDr1, fault=csz, 
                             normalModel=TRUE)
endParN = splineFit21k5N$MLEs

initPar=c(2, 1.5, 1, rep(0, nKnots-1))
splineFit21k5LN = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                              useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=G, 
                              fauxG=fauxG, constrLambda=FALSE, subDat=inflateDr1, fault=csz, 
                              normalModel=FALSE)
endParLN = splineFit21k5LN$MLEs

# generate subsidence predictions
params = splineFit21k5N$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5NormalInflate1.75"), subDat=inflateDr1, fault=csz, 
                   normalModel=TRUE)

params = splineFit21k5N$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5PosNormalInflate1.75"), subDat=inflateDr1, fault=csz, 
                   normalModel=TRUE, posNormalModel=TRUE)

params = splineFit21k5LN$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5LogNormalInflate1.75"), subDat=inflateDr1, fault=csz, 
                   normalModel=FALSE)

# check to find probability all slips are positive (negative slips are all negative)
library(mvtnorm)
params = splineFit21k5N$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
arealCSZCor = getArealCorMat(csz, normalModel=TRUE)
arealCSZCov = arealCSZCor * sigmaZeta^2
arealSlipCov = sweep(sweep(arealCSZCov, 1, tvec, "*"), 2, tvec, "*")
meanSlips = tvec * muZeta
par(mfrow=c(1,1))
plotMinMVN(meanSlips, arealSlipCov)

# compute truncated normal mean and covariance
library(tmvtnorm)
params = splineFit21k5N$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
arealCSZCor = getArealCorMat(csz, normalModel=TRUE)
arealCSZCov = arealCSZCor * sigmaZeta^2
arealSlipCov = sweep(sweep(arealCSZCov, 1, tvec, "*"), 2, tvec, "*")
n=nrow(csz)
n=20 # 3 seconds for n=10, 48 seconds for n=20 (==> O(n^4))
meanTest = rep(muZeta, n)
sigmaTest = arealCSZCov[1:n,1:n]
system.time(out <- mtmvnorm(mean=meanTest, sigma=sigmaTest, lower=rep(0, n), doComputeVariance = TRUE))

# now just compute the mean
n=100 # .3 seconds for n=10, 3 seconds for n=20, 6 for n=30, 11 for n=40, 18 for n=50, 27 for n=60, 51 for n=80, 86 for n=100
# ==> O(n^2) to O(n^3)
meanTest = rep(muZeta, n)
sigmaTest = arealCSZCov[1:n,1:n]
system.time(out <- mtmvnorm(mean=meanTest, sigma=sigmaTest, lower=rep(0, n), doComputeVariance = FALSE))

#####
# test predictions given subsidence for normal and pos normal models

# T1
isT1 = events=="T1"
# T1DatRange = dr1[isT1 & inRangeDat,]
# GT1Range = G[isT1 & inRangeDat, inRangeFault]
# T1Dat = dr1[isT1,]
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=1000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE)
posNormalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=1000, G=GT1, prior=FALSE, tvec=tvec, 
                                      normalModel=TRUE, posNormalModel=TRUE)

# areal values of zeta
muAreal = normalPreds$zetaEsts * tvec
sdAreal = normalPreds$zetaSD * tvec
medAreal = normalPreds$zetaMed * tvec
l95Areal = normalPreds$zeta025 * tvec
u95Areal = normalPreds$zeta975 * tvec

# get simulations
tab <- normalPreds$predResults
zetaSims = tab$zeta
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=FALSE, 
              fileNameRoot=paste0("21k5NormalT1"))
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("21k5NormalT1"), 
                   fault=csz, normalModel=TRUE, useMVNApprox=FALSE)

plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=FALSE, 
              fileNameRoot=paste0("21k5PosNormalT1"))
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("21k5PosNormalT1"), 
                   fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE)



#####
## cross-validation

# first fit the models
inflate=1.75
inflateDr1 = dr1
inflateDr1$Uncertainty = dr1$Uncertainty*inflate
nKnots=5
dStar=21000
initPar=c(20,15, 1, rep(0, nKnots-1))
splineFit21k5N = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                             useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=G, 
                             fauxG=fauxG, constrLambda=FALSE, subDat=inflateDr1, fault=csz, 
                             normalModel=TRUE)
endParN = splineFit21k5N$MLEs

initPar=c(2, 1.5, 1, rep(0, nKnots-1))
splineFit21k5LN = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                              useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=G, 
                              fauxG=fauxG, constrLambda=FALSE, subDat=inflateDr1, fault=csz, 
                              normalModel=FALSE)
endParLN = splineFit21k5LN$MLEs

### now do the cv
nFold=20
bySite=FALSE
if(bySite) { nFold = 21 } # since 21 different sites in T1 quake

# normal model
params = splineFit21k5N$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
outN = doCVSub(params, testEvent="T1", G=G, niter=500, tvec=tvec, normalModel=TRUE, nFold=nFold, seed=123, 
               subDat=inflateDr1, bySite=bySite)
plot(outN$MSEs)
outN$MSE

# positive normal model
params = splineFit21k5N$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
outPN = doCVSub(params, testEvent="T1", G=G, niter=500, tvec=tvec, normalModel=TRUE, posNormalModel=TRUE, 
                nFold=nFold, seed=123, subDat=inflateDr1, bySite=bySite)
plot(outPN$MSEs)
outPN$MSE

# lognormal model
params = splineFit21k5LN$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
outLN = doCVSub(params, testEvent="T1", G=G, niter=500, tvec=tvec, normalModel=FALSE, 
                nFold=nFold, seed=123, subDat=inflateDr1, bySite=bySite)
plot(outLN$MSEs)
outLN$MSE

# give all results
NCol = c(outN$MSE, outN$bias, outN$variance)
PNCol = c(outPN$MSE, outPN$bias, outPN$variance)
LNCol = c(outLN$MSE, outLN$bias, outLN$variance)
resTabSummary = cbind(NCol, PNCol, LNCol)
resTab = rbind(resTab, c(outN$MSE, outPN$MSE, outLN$MSE))
colnames(resTabSummary) = c("normal", "positive normal", "lognormal")
rownames(resTabSummary) = c("MSE", "bias", "variance")
colnames(resTab) = c("normal", "positive normal", "lognormal")
rownames(resTab) = c(paste0("split", 1:nFold), "Avg")
resTab
resTabSummary

# plot results
matplot(resTab[1:(nFold-1),], type="l", xlab="Split", ylab="Subsidence prediction MSE", lty=1, main="Cross-validation")
legend("top", c("Normal", "Positive Normal", "Lognormal"), col=c("black", "red", "green"), lty=1)

# save results
if(!bySite) { save(resTab, resTabSummary, file=paste0(nFold, "foldCV.RData")) }
if(bySite) { save(resTab, resTabSummary, file=paste0(nFold, "foldCVBySite.RData")) }

# plot points with different color by site
sites = unique(dr1$Site)
siteLats = aggregate(dr1$Lat, list(dr1$Site), mean)[,2]
sortI = sort(siteLats, index.return=TRUE)$ix
sites = sites[sortI]
# cols = tim.colors(length(sites))
cols = rainbow(length(sites))
lonRange = range(dr1$Lon)
latRange = range(dr1$Lat)
for(i in 1:length(sites)) {
  thisSite = sites[i]
  siteDat = dr1[dr1$Site == thisSite,]
  
  if(i == 1) {
    plot(siteDat$Lon, siteDat$Lat, pch="+", col=cols[i], xlab="Longitude", ylab="Latitude", main="Subsidence data", 
         xlim=lonRange, ylim=latRange)
  }
  else {
    points(siteDat$Lon, siteDat$Lat, pch="+", col=cols[i])
  }
}
map("world", "Canada", add=TRUE)
US(add=TRUE)

##### check traceplots of lognormal mcmc
# 1 magnitude
# 4 random areal slips
# 4 random pointwise slips

params = splineFit21k5LN$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)

##### plot predictive distribution

# T1
isT1 = events=="T1"
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1,]

# get T1 predictions
eventPreds = predsGivenSubsidence(params, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec)
predResults = eventPreds$predResults

# extract results from all chains
zetaAreal = extract(predResults, pars="beta", permuted=FALSE)
zetaPoint = extract(predResults, pars="zeta", permuted=FALSE)
mag = extract(predResults, pars="Mw", permuted=FALSE)

# extract specific par to plot
arealI = sample(1:dim(zetaAreal)[3], 4)
pointI = sample(1:dim(zetaPoint)[3], 4)
zetaAreaSelect = zetaAreal[,,arealI]
zetaPointSelect = zetaPoint[,,pointI]
magSelect = mag

### make plot
pdf("traceplotsT1n5000.pdf", width=8, height=8)
par(mfrow=c(3,3))
cols = rainbow(4)
# add magnitude traces to plot
plot(magSelect[,1,1], main="Magnitude", xlab="Iteration", ylab="Magnitude", 
     type="l", col=cols[1])
for(chain in 1:4) {
  lines(magSelect[,chain,1], col=cols[chain])
}
#legend("topleft", paste0("chain ", 1:4), col=cols, lty=1)
# add areal traces to plot
for(i in 1:4) {
  plot(zetaAreaSelect[,1,i], main=paste0("Zeta (Areal) ", arealI[i]), xlab="Iteration", ylab="Zeta", 
       type="l", col=cols[1])
  for(chain in 1:4) {
    lines(zetaAreaSelect[,chain,i], col=cols[chain])
  }
  #legend("topleft", paste0("chain ", 1:4), col=cols, lty=1)
}
# add pointwise traces to plot
for(i in 1:4) {
  plot(zetaPointSelect[,1,i], main=paste0("Zeta (Pointwise) ", pointI[i]), xlab="Iteration", ylab="Zeta", 
       type="l", col=cols[1])
  for(chain in 1:4) {
    lines(zetaPointSelect[,chain,i], col=cols[chain])
  }
  #legend("topleft", paste0("chain ", 1:4), col=cols, lty=1)
}
dev.off()



### now do the cv (but for marginal fits instead of predictive fits)

# normal model
niter=1000
outN = getEventMSE(csz, inflateDr1, "T1", niter, G, fauxG, TRUE, FALSE)
outN$MSE
outN$bias
outN$variance

# positive normal model
outPN = getEventMSE(csz, inflateDr1, "T1", niter, G, fauxG, TRUE, TRUE)
outPN$MSE
outPN$bias
outPN$variance

# lognormal model
outLN = getEventMSE(csz, inflateDr1, "T1", niter, G, fauxG, FALSE, FALSE)
outLN$MSE
outLN$bias
outLN$variance

# give all results
resTab = cbind(outN, outPN, outLN)
colnames(resTab) = c("normal", "positive normal", "lognormal")
rownames(resTab) = c("MSE", "bias", "variance")
resTab

# save results
save(resTab, file="marginalCV.RData")

#####
# display quality of data

# first plot histogram of quality levels (1 is best, 3 is worst)
# can look up quality levels in LeonardDatMod.xlsx:
# “Estimate quality: 
# 1, high-quality estimate from statistical microfossil analysis; 
# 2, medium-quality estimate - relative organic content, with macrofossil data on both sides of the contact and/or relative fresh/brackish diaton concentration; 
# 3, low-quality estimate - e.g., relative organic content, with macrofossil data for one/no sides of contact.

hist(as.numeric(dr1$quality), main="Subsidence data quality (1 is best)", freq=F, breaks=(1:4) - .5)

# now display qualities spatially
lonRange = range(dr1$Lon)
latRange = range(dr1$Lat)
qual1 = as.numeric(dr1$quality) == 1
qual2 = as.numeric(dr1$quality) == 2
qual3 = as.numeric(dr1$quality) == 3
plot(dr1$Lon[qual3], dr1$Lat[qual3], col="red", pch="+", xlim=lonRange, ylim=latRange, 
     main="Subsidence data quality", xlab="Longitude", ylab="Latitude")
points(dr1$Lon[qual2], dr1$Lat[qual2], col="blue", pch="+")
points(dr1$Lon[qual1], dr1$Lat[qual1], col="green", pch="+")
map("world", "Canada", add=TRUE)
US(add=TRUE)
legend("left", c("High quality", "Medium quality", "Low quality"), pch="+", col=c("green", "blue", "red"))

# get histograms of uncertainties for each type of quality
par(mfrow=c(3,1))
sdRange = c(0, max(dr1$Uncertainty))
aggregate(dr1$Uncertainty, list(dr1$quality), mean)
breakSeq = seq(0, max(dr1$Uncertainty), l=20)
hist(dr1$Uncertainty[qual1], freq=FALSE, main="High Quality Reported MOEs", 
     xlab="MOE", xlim=sdRange, breaks=breakSeq)
hist(dr1$Uncertainty[qual2], freq=FALSE, main="Medium Quality Reported MOEs", 
     xlab="MOE", xlim=sdRange, breaks=breakSeq)
hist(dr1$Uncertainty[qual2], freq=FALSE, main="Low Quality Reported MOEs", 
     xlab="MOE", xlim=sdRange, breaks=breakSeq)

# Why are the medium and low quality histograms so similar?  Let's eyeball the data
which(qual2)
dr1[28:33,]
dr1[34:38,]
dr1[37:40,]
dr1[54:60,]

##### try out using different error inflation levels for different types of data (multiply SEs by factor)
highQual = as.numeric(dr1$quality) == 1
lowQual = as.numeric(dr1$quality) != 1

inflateDR1 = dr1
inflates = rev(c(.5, .75, 1, 1.25, 1.5, 1.75, 2))
# inflates = rev(c(1, 1.25, 1.5, 1.75, 2))
# inflates=2
LLs = matrix(nrow=length(inflates), ncol=length(inflates))
LLsSub = matrix(nrow=length(inflates), ncol=length(inflates))
for(i in 1:length(inflates)) {
  highInflate = inflates[i]
  
  for(j in 1:length(inflates)) {
    lowInflate = inflates[j]
    
    inflateDR1 = dr1
    inflateDR1$Uncertainty[highQual] = dr1$Uncertainty[highQual] * highInflate
    inflateDR1$Uncertainty[lowQual] = dr1$Uncertainty[lowQual] * lowInflate
    
    nKnots=5
    dStar=21000
    initPar=c(20,15, 1, rep(0, nKnots-1))
    splineFit21k5N = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                                 useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=G, 
                                 fauxG=fauxG, constrLambda=FALSE, subDat=inflateDR1, fault=csz, 
                                 normalModel=TRUE)
    endPar = splineFit21k5$MLEs
    endPar = c(endPar[2], endPar[3], endPar[6:length(endPar)])
    
    # get full log likelihood and also subsidence log-likelihood
    LLs[i,j] = splineFit21k5N$logLikMLE
    LLsSub[i,j] = splineFit21k5N$optimTable[nrow(splineFit21k5N$optimTable), ncol(splineFit21k5N$optimTable)-3]
    
    # print results:
    print(paste0("Low inflate: ", lowInflate, ".  High inflate: ", highInflate, ".  LnLik: ", LLs[i,j], ".  LnLikSub: ", LLsSub[i,j]))
  }
}
save(inflates=inflates, LLs=LLs, LLsSub, 
     file="inflateResultsByQualFullFault.RData")

# now plot results.  Rows of LLs correspond to highInflates
par(mfrow=c(1,1))
highInflatesMat = matrix(rep(rev(inflates), length(inflates)), nrow=length(inflates))
lowInflatesMat = t(highInflatesMat)
LLs2 = t(matrix(rev(LLs), nrow=nrow(LLs)))
# image(lowInflatesMat, highInflatesMat, LLs2, col=tim.colors())
image.plot(rev(inflates), rev(inflates), LLs2, col=tim.colors(), 
           main="Normal Model Log-Likelihood", xlab="Low Quality Inflation", 
           ylab="High Quality Inflation")
max(LLs)
maxI = which.max(LLs)
highI = row(LLs)[maxI]
lowI = col(LLs)[maxI]
highInflate=inflates[highI]
lowInflate=inflates[lowI]
points(lowInflate, highInflate, col="green", pch="x")

# NOTE: the LnLik's given below are subsidence LnLiks.
#[1] "Low inflate: 2.  High inflate: 2.  LnLik: -1838.00087599582.  LnLikSub: -337.009699598503"
#[1] "Low inflate: 1.75.  High inflate: 2.  LnLik: -1826.90116210535.  LnLikSub: -326.125678681082"
#[1] "Low inflate: 1.5.  High inflate: 2.  LnLik: -1830.44536232665.  LnLikSub: -329.845110970792"
#[1] "Low inflate: 1.  High inflate: 1.  LnLik: -1932.17075038411.  LnLikSub: -429.495187972112"
#[1] "Low inflate: 2.  High inflate: 0.75.  LnLik: -1832.81652581612.  LnLikSub: -332.071334287872"
#[1] "Low inflate: 1.25.  High inflate: 0.75.  LnLik: -1868.04156777317.  LnLikSub: -368.13506057415"
#[1] "Low inflate: 0.75.  High inflate: 0.75.  LnLik: -1949.35637416385.  LnLikSub: -442.38016909445"

# examine whether subLL or GPSLL affect likelihood more
# muZeta sigmaZeta  beta1 beta2   beta3   beta4   beta5      LL   subLL   GPSLL priorLL LLSE
# newRow 160.88    133.08 4.6551 -7.1729 -8.8015 -1.4103 -8.2592 -2418.8 -919.02 -1499.8       0    0
# newRow 259.27    374.73 25.017 -31.59 -33.516 -15.082 -45.698 -2222.2 -713.01 -1509.2       0    0


## look at the locking rates divided by the fit taper function:

# first fit the model
inflate=1.75
inflateDr1 = dr1
inflateDr1$Uncertainty = dr1$Uncertainty*inflate
nKnots=5
dStar=21000
initPar=c(20,15, 1, rep(0, nKnots-1))
splineFit21k5N = doFitSpline(initParams=initPar, dStar=dStar, useMVNApprox=TRUE, 
                             useGrad=TRUE, nKnots=nKnots, maxit=500, useSubPrior=FALSE, useSlipPrior=FALSE, G=G, 
                             fauxG=fauxG, constrLambda=FALSE, subDat=inflateDr1, fault=csz, 
                             normalModel=TRUE)

params = splineFit21k5N$MLEs

# now get taper function values at locking rate estimate locations
splinePar = params[6:length(params)]
tempFault = list(depth=slipDatCSZ$Depth, latitude=slipDatCSZ$lat)
tvec = getTaperSpline(splinePar, tempFault, nKnots=nKnots, dStar=dStar)

# plot taper function
par(mfrow=c(1,1))
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, tvec, xlab="Longitude", ylab="Latitude", main="Taper function")
map("world", "Canada", add=TRUE)
US(add=TRUE)

# compute and plot renormalized locking rates (hopefully it looks relatively stationary)
normalizedVals = slipDatCSZ$slip/tvec
finite = is.finite(normalizedVals)
normalizedVals = normalizedVals[finite]

quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, normalizedVals, xlab="Longitude", ylab="Latitude", main="Renormalized Locking Rate")
map("world", "Canada", add=TRUE)
US(add=TRUE)

##### test difference between areal correlation matrix and non-areally average correlation matrix
#
# first get the areal corMat for lognormal and normal models
arealMatLN = getArealCorMat(csz)
corParLN = getCorPar()
corMatLN = stationary.cov(cbind(csz$longitude, csz$latitude), Covariance="Matern", 
                          theta=corParLN$phiZeta, smoothness=corParLN$nuZeta, 
                          Distance="rdist.earth", Dist.args=list(miles=FALSE))

arealMatN = getArealCorMat(csz, normalModel = TRUE)
corParN = getCorPar(normalModel = TRUE)
corMatN = stationary.cov(cbind(csz$longitude, csz$latitude), Covariance="Matern", 
                         theta=corParN$phiZeta, smoothness=corParN$nuZeta, 
                         Distance="rdist.earth", Dist.args=list(miles=FALSE))

# now compute difference between areal and non-areal corMats
max(abs(corMatLN - arealMatLN))
max(abs(diag(corMatLN) - diag(arealMatLN)))
mean(abs(corMatLN - arealMatLN))
mean(abs(diag(corMatLN) - diag(arealMatLN)))
max(abs(corMatN - arealMatN))
max(abs(diag(corMatN) - diag(arealMatN)))
mean(abs(corMatN - arealMatN))
mean(abs(diag(corMatN) - diag(arealMatN)))
# > max(abs(corMatLN - arealMatLN))
# [1] 0.01051342
# > max(abs(diag(corMatLN) - diag(arealMatLN)))
# [1] 0.002419194
# > mean(abs(corMatLN - arealMatLN))
# [1] 0.002289694
# > mean(abs(diag(corMatLN) - diag(arealMatLN)))
# [1] 0.001599329
# > max(abs(corMatN - arealMatN))
# [1] 0.01417535
# > max(abs(diag(corMatN) - diag(arealMatN)))
# [1] 0.004218268
# > mean(abs(corMatN - arealMatN))
# [1] 0.002383222
# > mean(abs(diag(corMatN) - diag(arealMatN)))
# [1] 0.002797084

# compute percentage difference
max(abs((corMatLN - arealMatLN)/arealMatLN))
max(abs((diag(corMatLN) - diag(arealMatLN))/diag(arealMatLN)))
mean(abs((corMatLN - arealMatLN)/arealMatLN))
mean(abs((diag(corMatLN) - diag(arealMatLN))/diag(arealMatLN)))
max(abs((corMatN - arealMatN)/arealMatN))
max(abs(diag(corMatN) - diag(arealMatN)))
mean(abs(corMatN - arealMatN))
mean(abs(diag(corMatN) - diag(arealMatN)))

# test what percent quadratic forms involving the non-areal correlation matrix will be from the areal
outDiffLN = eigen(corMatLN - arealMatLN)
outArealLN = eigen(arealMatLN)
outNonArealLN = eigen(corMatLN)

# NOTE: technically, should compute the Gaussian approximation to the LN covariance.
max(abs(outDiffLN$values))/max(abs(outArealLN$values))
max(abs(outNonArealLN$values))/max(abs(outArealLN$values))
# > max(abs(outDiffLN$values))/max(abs(outArealLN$values))
# [1] 0.004630316
# > max(abs(outNonArealLN$values))/max(abs(outArealLN$values))
# [1] 1.003401

outDiffN = eigen(solve(corMatN) - solve(arealMatN))
outArealN = eigen(arealMatN)
outNonArealN = eigen(corMatN)

max(abs(outDiffN$values))/max(abs(outArealN$values))
max(abs(outNonArealN$values))/max(abs(outArealN$values))
# > max(abs(outDiffN$values))/max(abs(outArealN$values))
# [1] 878.6425
# > max(abs(outNonArealN$values))/max(abs(outArealN$values))
# [1] 1.004137

# now test to see if we can do better in approximating the covariance:
areas = csz$length*csz$width
ds = sqrt(diag(arealMatN))
cs = (1/ds - 1)^2/sqrt(areas)
hist(cs)
C = mean(cs)

approxSDs = 1/(1 + sqrt(C*areas))
approxMat = sweep(sweep(corMatN, 1, approxSDs, "*"), 2, approxSDs, "*")
max(abs(approxMat - arealMatN))
max(abs(diag(approxMat) - diag(arealMatN)))
mean(abs(approxMat - arealMatN))
mean(abs(diag(approxMat) - diag(arealMatN)))

# take 2, just divide by the ds
approxMat = sweep(sweep(corMatN, 1, ds, "*"), 2, ds, "*")
max(abs(approxMat - arealMatN))
max(abs(diag(approxMat) - diag(arealMatN)))
mean(abs(approxMat - arealMatN))
mean(abs(diag(approxMat) - diag(arealMatN)))

# > mean(abs(approxMat - arealMatN))
# [1] 0.002356769
# > max(abs(approxMat - arealMatN))
# [1] 0.01214008
# > max(abs(diag(approxMat) - diag(arealMatN)))
# [1] 2.962075e-13
# > mean(abs(approxMat - arealMatN))
# [1] 0.002356769
# > mean(abs(diag(approxMat) - diag(arealMatN)))
# [1] 4.03344e-14

##### test how long it will take to compute the pointwise correlation matrix for fault vs. lock rates
corParN = getCorPar(normalModel = TRUE)
system.time(corMatN <- stationary.cov(cbind(csz$longitude, csz$latitude), Covariance="Matern", 
                                      theta=corParN$phiZeta, smoothness=corParN$nuZeta, 
                                      Distance="rdist.earth", Dist.args=list(miles=FALSE)))

system.time(corMatN <- stationary.cov(cbind(slipDatCSZ$lon, slipDatCSZ$lat), Covariance="Matern", 
                                      theta=corParN$phiZeta, smoothness=corParN$nuZeta, 
                                      Distance="rdist.earth", Dist.args=list(miles=FALSE)))


# plot slipDat subset to different max depths
maxDepth=30000
goodDepth = slipDat$Depth < maxDepth
quilt.plot(lon[goodDepth], lat[goodDepth], slip[goodDepth], xlab="Longitude", ylab="Latitude", 
           main="Locking Rate (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
quilt.plot(slipDat$lon[goodDepth], slipDat$lat[goodDepth], slipDat$Depth[goodDepth], 
           xlab="Longitude", ylab="Latitude", main="Depth (m)")
map("world", "Canada", add=TRUE)
world(add=TRUE)

# do the same for the data subset over the fault geometry
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, xlab="Longitude", ylab="Latitude", 
           main="Locking Rate (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$Depth, 
           xlab="Longitude", ylab="Latitude", main="Depth (m)")
map("world", "Canada", add=TRUE)
world(add=TRUE)

sum(goodDepth) # 964 for 30000

system.time(corMatN <- stationary.cov(cbind(slipDat$lon[goodDepth], slipDat$lat[goodDepth]), Covariance="Matern", 
                                      theta=corParN$phiZeta, smoothness=corParN$nuZeta, 
                                      Distance="rdist.earth", Dist.args=list(miles=FALSE)))

#####
# all mapping packages and setup:
library(ggmap)
library(ggplot2)
library(mapdata)
library(maptools)
library(maps)
library(RColorBrewer)
library(sp)
library(gstat)

lon = c(-128, -122)
lat = c(40, 50)

### Get a map
# https://mapstyle.withgoogle.com/
# map <- get_map(location = c(lon[1], lat[1], lon[2], lat[2]), zoom = 6,
#                maptype = "terrain", source = "google")
style1=c(feature="administrative", element="labels", visibility="off")
style2=c("&style=", feature="road", element="geometry", visibility="off")
style3=c("&style=", feature="poi", element="labels", visibility="off")
style4=c("&style=", feature="landscape", element="labels", visibility="off")
style5=c("&style=", feature="administrative", element="geometry.stroke", color="black")
style6=c("&style=", feature="administrative", element="geometry.stroke", weight=.75)
map <- get_googlemap(center=c(lon = mean(lon), lat = mean(lat)), zoom=5,
                     style=c(style1, style2, style3, style4, style5, style6))

states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
canada <- map_data("world", "Canada")
coastMap = rbind(canada, west_coast)

#####
# try using google maps to plot data (note: ggmap does not support different projections)
library(ggmap)
library(ggplot2)

### Set a range
# lon = c(-128.5, -121.5)
lon = c(-128, -122)
lat = c(40, 50)

### Get a map
# https://mapstyle.withgoogle.com/
# map <- get_map(location = c(lon[1], lat[1], lon[2], lat[2]), zoom = 6,
#                maptype = "terrain", source = "google")
style1=c(feature="administrative", element="labels", visibility="off")
style2=c("&style=", feature="road", element="geometry", visibility="off")
style3=c("&style=", feature="poi", element="labels", visibility="off")
style4=c("&style=", feature="landscape", element="labels", visibility="off")
style5=c("&style=", feature="administrative", element="geometry.stroke", color="black")
style6=c("&style=", feature="administrative", element="geometry.stroke", weight=.75)
map <- get_googlemap(center=c(lon = mean(lon), lat = mean(lat)), zoom=5,
                     style=c(style1, style2, style3, style4, style5, style6))
# map <- get_googlemap(center=c(lon = mean(lon), lat = mean(lat)), zoom=6,
# style='feature:administrative|element:labels|visibility:off&style=feature:road|element:labels|visibility:off')

### When you draw a figure, you limit lon and lat, scramble site labels by latitude so nearby sites 
### have very different colors
fauxObs = data.frame(list(Lon=dr1$Lon, Lat=dr1$Lat))
goodScramble = c(19, 2, 7, 5, 3, 20, 12, 9, 4, 6, 16, 11, 17, 15, 21, 14, 10, 18, 13, 1, 8)
sites = unique(dr1$Site)
siteLats = aggregate(dr1$Lat, list(dr1$Site), mean)[,2]
sortI = sort(siteLats, index.return=TRUE)$ix
sites = sites[sortI]
sites = sites[goodScramble]

mapPoints <- ggmap(map) + 
  geom_point(aes(x = Lon, y = Lat, color=factor(Site, levels=sites)), data=dr1, size=3, shape="+")
foo <- mapPoints + 
  scale_x_continuous(limits = lon, expand = c(0, 0)) +
  scale_y_continuous(limits = lat, expand = c(0, 0)) + 
  ggtitle("Subsidence Estimates") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Longitude", y="Latitude")

foo + guides(color=FALSE)

# now test with faultGeom polygons

foo + ggplotFault(faultGeom) + guides(color=FALSE)

bg = ggmap(map) + 
  scale_x_continuous(limits = lon, expand = c(0, 0)) +
  scale_y_continuous(limits = lat, expand = c(0, 0)) + 
  ggtitle("Subsidence Estimates") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Longitude", y="Latitude")

mapPoints = geom_point(aes(x = Lon, y = Lat, color=factor(Site, levels=sites)), data=dr1, size=3, shape="+")

faultPoly = ggplotFault(faultGeom, color=rgb(.2,.2,.2))

bg + faultPoly + mapPoints + guides(color=FALSE)

##### now plot the locking rate data using grey background map

library(mapdata)
library(maptools)
library(maps)
library(RColorBrewer)

# get map data
states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
canada <- map_data("world", "Canada")
coastMap = rbind(canada, west_coast)

# plot it (choose background color with fields.color.picker())
ggplot(data = coastMap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = lon,  ylim = lat, ratio = 1.3) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='lightblue1'))

# now interpolate and overlay locking rate data
# try using inverse distance weighting inerpolation
library(sp)
library(gstat)

bg = ggplot(data = coastMap) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black") + 
  coord_fixed(xlim = lon,  ylim = lat, ratio = 1.3) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='lightblue1'))

# make grid of prediction locations
predGrid = expand.grid(x = seq(from = lon[1], to = lon[2], by = 0.05),
                       y = seq(from = lat[1], to = lat[2], by = 0.05))  # expand points to grid
coordinates(predGrid) <- ~x + y
gridded(predGrid) <- TRUE

# interpolate data using inverse weighted distance
preds <- idw(formula = slip ~ 1, locations = ~lon+lat, 
             newdata = predGrid, data=slipDatCSZ)
preds = as.data.frame(preds)
preds$lon = preds$x
preds$lat = preds$y
preds$slip = preds$var1.pred
preds$x = NULL
preds$y = NULL
preds$var1.pred = NULL
preds$var1.var = NULL

# subset data to only be within fault geometry:
predsCSZ = getFaultGPSDat(preds)

# lockDat = slipDatCSZ
# lockDat$long = lockDat$lon
# bg + geom_raster(aes(fill=slip, x=lon, y=lat), data=lockDat, interpolate=FALSE) + 
#   coord_fixed(xlim=lon,ylim=lat, ratio=1.3)

# overlay interpolated data on plot:
bg + geom_tile(data = predsCSZ, aes(x = lon, y = lat, fill = slip)) + 
  geom_point(aes(x=lon, y=lat), pch=20, col="black", size=.1, data=slipDatCSZ)

# subset slip data to have depth < 30km
slipDat30 = slipDat[slipDat$Depth < 30000,]

##### test new model
inflate=2
inflateDr1 = dr1
inflateDr1$Uncertainty = dr1$Uncertainty*inflate
nKnots=5
dStar=27000 # must be larger than the maximum depth of the GPS locking rate data at the very least
# initPar=c(20,15, 1/20000, rep(0, nKnots-1), 175)
initPar=c(20,15, 1, rep(0, nKnots-1), 175)
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=slipDatCSZ, 
                        useGrad=FALSE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=inflateDr1, fault=csz, 
                        normalModel=TRUE, normalizeTaper=TRUE)
endParN = fitTwo21k5N$MLEs

### inflate=1.75, 1.25
## using gradient
# dStar=25000, depthThresh=25000
# initPar=c(20,15, 1, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta beta1  beta2   beta3   beta4  beta5    phi      LL  subLL   GPSLL priorLL LLSE
# newRow 12.239    14.985 1.533 -3.002 0.20675 -2.1421 -1.307 164.11 -1959.3 -398.8 -1560.5       0    0
# NOTE: taper has very little variation

# dStar=25000, depthThresh=21000
#        muZeta sigmaZeta  beta1   beta2   beta3 beta4   beta5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow 14.684    13.843 1.6286 -3.7164 -2.3272 -1.43 -1.6193 141.47 -1840.9 -385.99 -1454.9       0    0
# NOTE: taper has some variation, niter=1000 (2000 simulations) takes ~ 11 minutes

# dStar=30000, depthThresh=21000
#        muZeta sigmaZeta  beta1   beta2  beta3   beta4   beta5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow 12.555    14.279 2.1703 -4.6894 2.3972 -3.7668 -1.6706 199.86 -1827.9 -401.81 -1426.1       0    0
# NOTE: taper has some variation, but predictions take too long for posNormal model

# dStar=21000, depthThresh=21000
#        muZeta sigmaZeta   beta1   beta2    beta3   beta4  beta5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow 15.727    17.674 0.16658 0.75713 -0.89657 0.15486 1.5713 114.04 -1998.7 -389.98 -1608.7       0    0
# NOTE: basically no variation

# dStar=27000, gpsDat=slipDatCSZ
#        muZeta sigmaZeta  beta1   beta2  beta3   beta4   beta5   phi      LL   subLL   GPSLL priorLL LLSE
# newRow 11.915    14.511 1.8903 -3.9699 1.0144 -2.7666 -1.6573 179.6 -1944.6 -401.16 -1543.4       0    0
# NOTE: taper has minor variations, predictions take a pretty long time (~11?)

# dStar=25000, depthThresh=21000, downsample to nrow(dr1) GPS observations
#        muZeta sigmaZeta  beta1   beta2 beta3   beta4   beta5    phi      LL   subLL GPSLL priorLL LLSE
# newRow 14.637    13.364 1.6243 -3.6334 -2.67 -1.2649 -1.6643 135.34 -1670.1 -382.16 -1288       0    0
# NOTE: taper has significant variation

# dStar=28000
# initPar=c(20,15, 2.5, rep(0, nKnots-1), 200)
#        muZeta sigmaZeta  beta1  beta2  beta3   beta4  beta5    phi   LL      subLL    GPSLL priorLL LLSE
# newRow 11.499    13.641 2.0134 -4.366 1.3107 -3.0586 -1.739 182.26 -1928.5 -404.09 -1524.4       0    0

# dStar=100000
# initPar=c(20,15, 5, rep(0, nKnots-1), 200)
#         muZeta sigmaZeta  beta1   beta2  beta3  beta4   beta5   phi      LL   subLL   GPSLL priorLL LLSE
# newRow 8.6401    10.768 6.9745 -10.802 2.7437 -7.803 -8.6106 221.6 -1894.9 -404.77 -1490.2       0    0

# dStar = 45000
# initPar=c(20,15, 3, rep(0, nKnots-1), 200)
#        muZeta sigmaZeta  beta1   beta2  beta3   beta4  beta5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow 9.7421     11.73 3.1756 -5.0768 1.5366 -3.6636 -3.721 212.98 -1896.4 -405.35 -1491.1       0    0

# normalizeTaper=FALSE
# initPar=c(20,15, 1/21000, rep(0, nKnots-1), 200)
#        muZeta sigmaZeta      beta1      beta2      beta3       beta4       beta5    phi      LL   subLL   GPSLL
# newRow 8.6342    10.794 6.9459e-05 -0.0001056 2.6778e-05 -7.8383e-05 -8.5569e-05 221.72 -1895.1 -404.63 -1490.5

# normalizeTaper=FALSE
# initPar=c(20,15, 1/40000, rep(0, nKnots-1), 200)
#        muZeta sigmaZeta     beta1      beta2      beta3       beta4       beta5    phi      LL   subLL   GPSLL
# newRow 8.6336    10.794 6.946e-05 -0.0001056 2.6779e-05 -7.8393e-05 -8.5571e-05 221.73 -1895.1 -404.63 -1490.5

# normalizeTaper=FALSE
# initPar=c(20,15, 1/40000, rep(0, nKnots-1), 150)
#        muZeta sigmaZeta     beta1      beta2      beta3       beta4       beta5    phi      LL   subLL   GPSLL
# newRow 8.6339    10.794 6.946e-05 -0.0001056 2.6779e-05 -7.8392e-05 -8.5571e-05 221.72 -1895.1 -404.63 -1490.5

# normalizeTaper=FALSE
# initPar=c(20,15, 1/20000, rep(0, nKnots-1), 300)
#        muZeta sigmaZeta     beta1      beta2      beta3       beta4       beta5    phi      LL   subLL   GPSLL
# newRow 8.6337    10.794 6.946e-05 -0.0001056 2.6779e-05 -7.8392e-05 -8.5568e-05 221.73 -1895.1 -404.63 -1490.5

# normalizeTaper=FALSE
# maxit=500, reltol=1e-9
# initPar = c(11.876,    19.821, -6.2747e-05, 4.5077e-05, 0.00013932, 9.3883e-05, 5.8227e-05, 319.48)
#        muZeta sigmaZeta       beta1      beta2     beta3     beta4      beta5    phi      LL   subLL   GPSLL priorLL
# newRow 9.1936    11.793 -6.1316e-05 5.0756e-05 0.0001365 6.845e-05 5.3795e-05 229.13 -1896.1 -405.38 -1490.7       0
# maxit=750, reltol=1e-10:
#        muZeta sigmaZeta       beta1      beta2      beta3      beta4      beta5    phi    LL   subLL   GPSLL priorLL
# newRow 9.2059    11.693 -6.1696e-05 5.1795e-05 0.00013657 7.0527e-05 4.1256e-05 226.33 -1896 -405.15 -1490.9       0
# maxit=750, reltol=1e-11:
#       muZeta sigmaZeta       beta1      beta2      beta3      beta4      beta5    phi    LL   subLL   GPSLL 
# newRow 9.2059    11.693 -6.1696e-05 5.1795e-05 0.00013657 7.0527e-05 4.1256e-05 226.33 -1896 -405.15 -1490.9
# maxit=1000, reltol=1e-11:

## without using gradient
# initPar=c(20,15, 1/20000, rep(0, nKnots-1), 150)
#        muZeta sigmaZeta       beta1      beta2      beta3      beta4      beta5    phi      LL   subLL   GPSLL
# newRow 11.876    19.821 -6.2747e-05 4.5077e-05 0.00013932 9.3883e-05 5.8227e-05 319.48 -1901.1 -411.89 -1489.2
# 

### inflate=.5
#~        muZeta sigmaZeta      beta1       beta2       beta3     beta4       beta5    phi      LL   subLL   GPSLL
#~ newRow 31.538    35.781 0.00010245 -5.9259e-05 -5.2124e-05 -0.000133 -0.00014906 97.665 -2531.2 -915.08 -1616.1

### inflate = 1
#~         muZeta sigmaZeta       beta1      beta2       beta3      beta4      beta5    phi      LL   subLL   GPSLL
#~ newRow  39.92    38.264 -6.7397e-05 4.0907e-05 -3.8734e-05 0.00012672 2.9629e-05 126.54 -1987.7 -462.15 -1525.5

### inflate = 2
#~        muZeta sigmaZeta      beta1      beta2      beta3       beta4       beta5    phi      LL   subLL   GPSLL
#~ newRow 8.0731    10.283 6.9332e-05 -0.0001055 2.2026e-05 -7.1259e-05 -8.3263e-05 228.94 -1908.8 -419.12 -1489.7

initPar=c(2, 1.5, 1, rep(0, nKnots-1), 250)
fitTwo21k5LN = fitModel2(initParams=initPar, dStar=dStar, 
                         useGrad=FALSE, nKnots=nKnots, maxit=500, G=G, 
                         fauxG=fauxG, subDat=inflateDr1, fault=csz, 
                         normalModel=FALSE)
endParLN = fitTwo21k5LN$MLEs

# generate subsidence predictions
params = fitTwo21k5N$MLEs
splinePar = params[6:(length(params)-1)]
dStar=100000
splinePar=c(1/40000, rep(0, 4))
normalizeTaper=FALSE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, normalize=FALSE)
plotFault(csz, tvec, varRange=c(0,1))
comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5TwoNormalInflate1.75"), subDat=inflateDr1, fault=csz, 
                   normalModel=TRUE, taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=normalizeTaper)

params = fitTwo21k5N$MLEs
splinePar = params[6:(length(params)-1)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5TwoPosNormalInflate1.75"), subDat=inflateDr1, fault=csz, 
                   normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=normalizeTaper)

params = fitTwo21k5LN$MLEs
splinePar = params[6:(length(params)-1)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="5 knots ", savePlots=TRUE, tvec=tvec, nsim=2000, 
                   fileNameRoot=paste0("21k5TwoLogNormalInflate1.75"), subDat=inflateDr1, fault=csz, 
                   normalModel=FALSE, taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=normalizeTaper)


##### now get inflation parameters

# determine data quality
highQual = as.numeric(dr1$quality) == 1
lowQual = as.numeric(dr1$quality) != 1

# inflates = rev(c(1, 1.25, 1.5, 1.75, 2))
# inflates=2
inflates = rev(c(.5, .75, 1, 1.25, 1.5, 1.75, 2))
LLs = matrix(nrow=length(inflates), ncol=length(inflates))
LLsSub = matrix(nrow=length(inflates), ncol=length(inflates))
for(i in 1:length(inflates)) {
  highInflate = inflates[i]
  
  for(j in 1:length(inflates)) {
    lowInflate = inflates[j]
    
    inflateDR1 = dr1
    inflateDR1$Uncertainty[highQual] = dr1$Uncertainty[highQual] * highInflate
    inflateDR1$Uncertainty[lowQual] = dr1$Uncertainty[lowQual] * lowInflate
    
    nKnots=5
    dStar=25000 # dStar doesn't matter since we're not normalizing the taper
    initPar=c(20,15, 1, rep(0, nKnots-1), 175)
    fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                            useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                            fauxG=fauxG, subDat=inflateDR1, fault=csz, 
                            normalModel=TRUE, normalizeTaper=TRUE, doHess=FALSE)
    endPar = fitTwo21k5N$MLEs
    endPar = c(endPar[2], endPar[3], endPar[6:length(endPar)])
    
    # get full log likelihood and also subsidence log-likelihood
    LLs[i,j] = fitTwo21k5N$logLikMLE
    LLsSub[i,j] = fitTwo21k5N$optimTable[nrow(fitTwo21k5N$optimTable), ncol(fitTwo21k5N$optimTable)-3]
    
    # print results:
    print(paste0("Low inflate: ", lowInflate, ".  High inflate: ", highInflate, ".  LnLik: ", LLs[i,j], ".  LnLikSub: ", LLsSub[i,j]))
  }
}
save(inflates=inflates, LLs=LLs, LLsSub, 
     file="inflateResultsByQualFullFaultTwo.RData")

# now plot results.  Rows of LLs correspond to highInflates
# 1.75 for low, 1.25 for high quality inflation, no matter whether using the full or subsidence log-likelihood!
par(mfrow=c(1,1))
highInflatesMat = matrix(rep(rev(inflates), length(inflates)), nrow=length(inflates))
lowInflatesMat = t(highInflatesMat)
tmp = LLsSub
LLs2 = t(matrix(rev(tmp), nrow=nrow(tmp)))
# image(lowInflatesMat, highInflatesMat, LLs2, col=tim.colors())
image.plot(rev(inflates), rev(inflates), LLs2, col=tim.colors(), 
           main="Normal Model Log-Likelihood", xlab="Low Quality Inflation", 
           ylab="High Quality Inflation")
max(tmp)
maxI = which.max(tmp)
highI = row(tmp)[maxI]
lowI = col(tmp)[maxI]
highInflate=inflates[highI]
lowInflate=inflates[lowI]
points(lowInflate, highInflate, col="green", pch="x")


## test predictions given subsidence for normal and pos normal models

# subset GPS data to be shallower than some threshold
# depthThresh = max(getFaultCenters(csz,3)[,3]) 
depthThresh=21000
set.seed(123)
threshSlipDat = slipDat[slipDat$Depth<depthThresh,]
# ids = sample(1:nrow(threshSlipDat), round(nrow(inflateDr1)*.333))
ids = sample(1:nrow(threshSlipDat), 189)
threshSlipDat = threshSlipDat[ids,]
quilt.plot(threshSlipDat$lon, threshSlipDat$lat, threshSlipDat$slip)
plotFault(csz, plotData=F, new=F)

# first fit the model
# highInflate=1.00
# lowInflate=1.75
# highQual = as.numeric(dr1$quality) == 1
# lowQual = as.numeric(dr1$quality) != 1
# inflateDr1 = dr1
inflateDr1$Uncertainty[highQual] = dr1$Uncertainty[highQual]*highInflate
inflateDr1$Uncertainty[lowQual] = dr1$Uncertainty[lowQual]*lowInflate
nKnots=5
dStar=25000 # must be larger than the maximum depth of the GPS locking rate data at the very least
initPar=c(20,15, 1/20000, rep(0, nKnots-1), 175)
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=inflateDr1, fault=csz, 
                        normalModel=TRUE, normalizeTaper=FALSE)

# tried (with depth threshold at max depth of fault geometry):
# dStar=24000
# dStar=25000

# tried with depth threshold at 25000

# tried with slipDatCSZ
# dStar=27000 (FAIL: possible, but would take ~11 )
initPar=c(20,15, 1, rep(0, nKnots-1), 175)
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=inflateDr1, fault=csz, 
                        normalModel=TRUE, normalizeTaper=TRUE)

# T1
isT1 = events=="T1"
# T1DatRange = dr1[isT1 & inRangeDat,]
# GT1Range = G[isT1 & inRangeDat, inRangeFault]
# T1Dat = dr1[isT1,]
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

params = fitTwo21k5N$MLEs
splinePar = params[6:(length(params)-1)]
Xi = getSplineBasis(csz, c(40,50), 5)
lambdas = Xi %*% splinePar
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, normalize=TRUE)
plotFault(csz, tvec)

normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=1000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=slipDatCSZ)
posNormalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=1000, G=GT1, prior=FALSE, tvec=tvec, 
                                      normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                      dStar=dStar, gpsDat=slipDatCSZ)

# areal values of zeta
muAreal = normalPreds$zetaEsts * tvec
sdAreal = normalPreds$zetaSD * tvec
medAreal = normalPreds$zetaMed * tvec
l95Areal = normalPreds$zeta025 * tvec
u95Areal = normalPreds$zeta975 * tvec

# get simulations
tab <- normalPreds$predResults
zetaSims = tab$zeta
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=FALSE, 
              fileNameRoot=paste0("27k5TwoNormalT1ThreshSlipDatCSZ"))
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("27k5TwoNormalT1SlipDatCSZ"), 
                   fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                   dStar=dStar, normalizeTaper=FALSE)

# areal values of zeta
muAreal = posNormalPreds$zetaEsts * tvec
sdAreal = posNormalPreds$zetaSD * tvec
medAreal = posNormalPreds$zetaMed * tvec
l95Areal = posNormalPreds$zeta025 * tvec
u95Areal = posNormalPreds$zeta975 * tvec

# get simulations
tab <- posNormalPreds$predResults
zetaSims = tab$zeta
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=FALSE, 
              fileNameRoot=paste0("25k5TwoPosNormalT121k"))
comparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("27k5TwoPosNormalT1SlipDatCSZ"), 
                   fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                   dStar=dStar, normalizeTaper=FALSE)


plotFault(csz, getFaultCenters(csz, 3)[,3])
quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$Depth)
plotFault(csz, plotData = FALSE, new=FALSE)

quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$Depth)
plotFault(csz, plotData = FALSE, new=FALSE)
quilt.plot(threshSlipDat$lon, threshSlipDat$lat, threshSlipDat$Depth)
plotFault(csz, plotData = FALSE, new=FALSE)


plotFault(csz, muc[-(1:nrow(threshSlipDat))])
plotFault(csz, sqrt(diag(Sigmac[-(1:nrow(threshSlipDat)),-(1:nrow(threshSlipDat))])))

# look at variation in taper values along the middle of a fault
hist(csz$depth, breaks=50)
sum(csz$depth < 3000)
sum(csz$depth<5000 & (csz$depth > 3000))
sum(csz$depth<10000 & (csz$depth > 5000))
sum(csz$depth<12000 & (csz$depth > 10000))
sum(csz$depth<15000 & (csz$depth > 12000))
sum(csz$depth<12000)
lowFault = csz$depth < 3000
lowMidFault = csz$depth<5000 & (csz$depth > 3000)
midLowFault = csz$depth<10000 & (csz$depth > 5000)
midFault = csz$depth<12000 & (csz$depth > 10000)
midHighFault = csz$depth<15000 & (csz$depth > 12000)
highFault = csz$depth > 15000
faultFilt = lowMidFault

plot(csz$latitude[faultFilt], tvec[faultFilt], type="l", col="blue", main="Taper Variation", 
     xlab="Latitude", ylab="Taper", ylim=c(0,1))

xs = cbind(csz$latitude[lowFault], csz$latitude[lowMidFault], csz$latitude[midLowFault], 
           csz$latitude[midFault], csz$latitude[midHighFault], csz$latitude[highFault])
ys = cbind(tvec[lowFault], tvec[lowMidFault], tvec[midLowFault], 
           tvec[midFault], tvec[midHighFault], tvec[highFault])
matplot(xs, ys, type="l", main="Taper Variation Vs. Depth", xlab="Latitude", ylab="Taper", ylim=c(0,1), lty=1)


# compute effective degrees of freedom for the gps data

gpsDat = slipDatCSZ
gpsDat = slipDat
gpsDat = slipDat[slipDat$Depth < 21000,]
out = Krig(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, GCV=TRUE, method="GCV", Distance="rdist.earth")
out2 = MLESpatialProcess(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, 
                         cov.args=list(Distance="rdist.earth", Dist.args=list(miles=FALSE)), 
                         theta.start=175, theta.range=c(75, 250))
out$eff.df # too many
test = out$matrices

# get covariance and cross-covariance matrices
MLEs = out$best.model
lambda = MLEs[1] # tausq/rho
tausq = MLEs[2]
rho = MLEs[3] # or sigmasq, partial sill

# fit spatial model to gps data
thetas = seq(100, 500, l=30)
maxLnLik = -Inf
maxPar = NULL
maxTheta = NULL
for(i in 1:length(thetas)) {
  print(paste0("i: ", i))
  out = myKrig(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, theta=thetas[i], GCV=TRUE, method="GCV",
               cov.args=list(Distance="rdist.earth", Dist.args=list(miles=FALSE)), m=2,
               lambda.grid=seq(0, 5e-2, l=200), weights=1/gpsDat$slipErr^2)
  # out = myKrig(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, theta=thetas[i], GCV=TRUE, method="GCV", 
  #              cov.args=list(Distance="rdist.earth", Dist.args=list(miles=FALSE)), m=2, 
  #              lambda.grid=seq(0, 5e-2, l=200))
  bestModI = which.min(out$gcv.grid$GCV)
  lnLik = -out$gcv.grid$`-lnLike Prof`[bestModI]
  print(paste0("lnLik: ", lnLik))
  print(paste0("effDF: ", out$eff.df))
  if(lnLik > maxLnLik) {
    maxPar = out$best.model
    maxTheta = thetas[i]
    maxLnLik = lnLik
    maxBeta = out$d
    bestMod = out
  }
}
maxTheta
lnLik
maxTausq = maxPar[2]
maxRho = maxPar[3]

## estimate effective degrees of freedom based on linear smoothers estimator
# get covariance/cross-covariance
nuggetPart = diag(maxTausq*gpsDat$slipErr^2)
nuggetPart = maxTausq * diag(nrow=nrow(gpsDat))
SigmaY = maxRho*stationary.cov(cbind(gpsDat$lon, gpsDat$lat), Distance="rdist.earth", 
                               Dist.args=list(miles=FALSE), theta=maxTheta) + nuggetPart
SigmaZY = maxRho*stationary.cov(cbind(gpsDat$lon, gpsDat$lat), Distance="rdist.earth", 
                                Dist.args=list(miles=FALSE), theta=maxTheta)
SigmaYInv = solve(SigmaY)
mult = SigmaZY %*% SigmaYInv
# X = cbind(1, gpsDat$lon, gpsDat$lat)
# X = matrix(1, nrow=nrow(gpsDat), ncol=1)
X = do.call(bestMod$null.function.name, c(bestMod$null.args, list(x=bestMod$x, Z=bestMod$Z)))
Smat1 = (diag(nrow(gpsDat)) - mult) %*% X %*% solve(t(X) %*% SigmaYInv %*% X) %*% t(X) %*% SigmaYInv
Smat2 = mult
Smat = Smat1 + Smat2
nrow(gpsDat) - sum(diag(Smat))

out = Tps(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, m=2, lon.lat=TRUE, miles=FALSE, 
          df=200)
out = fastTps(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, m=2, lon.lat=TRUE, 
              cov.args=list(Dist.args=list(miles=FALSE)), lambda=0, theta=175)

# NOTE: fault geometry is ~200km wide at widest point
quilt.plot(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip)
plotFault(csz, new=F, plotData=F)
rdist.earth(cbind(c(-126, -124), c(47, 48)), miles=F)
theta = 175
theta=233
theta=200 # results in 189 EffDF for weighted fit with linear trend term
out = myKrig(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, theta=theta, GCV=TRUE, method="GCV", 
             cov.args=list(Distance="rdist.earth", Dist.args=list(miles=FALSE)), m=2, 
             lambda.grid=seq(0, 5e-2, l=200), weights=1/gpsDat$slipErr^2)
# out = myKrig(cbind(gpsDat$lon, gpsDat$lat), gpsDat$slip, theta=theta, GCV=TRUE, method="GCV", 
#              cov.args=list(Distance="rdist.earth", Dist.args=list(miles=FALSE)), m=2, 
#              lambda.grid=seq(0, 5e-2, l=200))
maxPar = out$best.model
beta = out$d
tausq = maxPar[2]
rho = maxPar[3]
out$eff.df
nuggetPart = diag(maxTausq*gpsDat$slipErr^2)
nuggetPart = maxTausq * diag(nrow=nrow(gpsDat))
SigmaY = maxRho*stationary.cov(cbind(gpsDat$lon, gpsDat$lat), Distance="rdist.earth", 
                               Dist.args=list(miles=FALSE), theta=theta) + nuggetPart
SigmaZY = maxRho*stationary.cov(cbind(gpsDat$lon, gpsDat$lat), Distance="rdist.earth", 
                                Dist.args=list(miles=FALSE), theta=maxTheta)

## other estimates of EffDF based on Bretherton 1992

# Nmm (not possible since requires multiple timesteps)
# E = sum(gpsDat$slip^2)
# 2*mean(E)/var(E)

# Neff
C = SigmaY
sum(diag(C))^2/sum(C^2)


initPar=c(14.637, 13.364, 1.6243, -3.6334, -2.67, -1.2649, -1.6643, 135.34)
initPar=c(14.626, 13.379, 5.0016, -3.9156, -2.9471, -1.3274, -1.6608, 133.51)
initPar=c(14.626, 13.379, 5.0016, -3, -3, -3, -3, 133.51)

# removing GPS data:
# initPar = c(14.637, 13.364, 1.6243, -3.6334, -2.67, -1.2649, -1.6643, 135.34)
#        muZeta sigmaZeta  beta1   beta2   beta3   beta4   beta5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow 14.626    13.379 1.0016 -3.9156 -2.9471 -1.3274 -1.6608 133.51 -364.46 -355.72 -8.7397       0    0
# optGrad = 0

# initPar=c(20,15, 1, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta beta1   beta2   beta3     beta4     beta5    phi      LL   subLL  GPSLL priorLL LLSE
# newRow 19.948    15.127 1.984 0.55992 0.36835 -0.090501 -0.021691 169.75 -411.33 -402.88 -8.446       0    0
# optGrad = 0

# initPar=c(20,15, 1, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta  beta1   beta2   beta3    beta4    beta5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow 19.929    15.097 1.6395 0.53446 0.33206 -0.11605 -0.02371 146.16 -372.25 -363.79 -8.4593       0    0

# with unnormalized taper:
# initPar=c(20,15, 1/20000, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta      beta1      beta2     beta3       beta4       beta5    phi      LL   subLL   GPSLL priorLL  LLSE
# newRow 19.552    15.714 5.1016e-05 2.2037e-05 7.739e-05 -2.2877e-05 -5.9492e-06 161.14 -371.54 -363.02 -8.5164       0     0

# with updated par scales:
# muZeta sigmaZeta      beta1      beta2      beta3       beta4      beta5   phi      LL   subLL   GPSLL priorLL LLSE
# newRow 18.191    17.883 5.1355e-05 2.3333e-05 7.5563e-05 -2.6648e-05 -6.266e-06 160.9 -365.66 -356.87 -8.7887       0    0

# initPar=c(15,10, 1, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta  beta1  beta2   beta3    beta4     beta5    phi      LL   subLL   GPSLL
# newRow  14.75     11.65 1.6138 0.5853 0.26336 -0.31077 -0.036934 164.07 -372.12 -367.67 -4.4532
# optGrad = 0

# with correlated gps data:
# initPar=c(15,10, 1, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta   beta1   beta2   beta3   beta4    beta5    phi      LL   subLL  GPSLL
# newRow 17.296    12.639 0.61984 -2.5301 -1.6332 0.19315 -0.92437 128.04 -757.03 -392.23 -364.8
# fitTwo21k5N$optGrad ~ 1e-4 to 1e-6
# > sqrt(diag(solve(-fitTwo21k5N$hess)))
# [1] 0.541487232 0.232573154 0.011843770 0.190905291 0.128490947 0.005485611 0.026827882 2.729749857

### now throw out most of GPS data:
## using nelder-mead:
# initPar=c(14.637, 13.364, -0.0400, 0.3994, -1.0057, -1.9691, 1.6643, 135.34)
# muZeta sigmaZeta   beta1  beta2  beta3   beta4  beta5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow  18.35    5.0468 -0.1796 2.9964 -3.979 -2.3043 6.7951 28.375 -288.01 -284.93 -3.0827       0    0
# optPar: c(18.3503740, 5.0468239, -0.1796016, 2.9964006, -3.9789832, -2.3042896,  6.7951043, 28.3749764)

#### using contrOptim (makes robust gradients I guess even when constraint doesn't end up mattering?):
### throw out GPS data further south than 41 lat, <20km depth, dStar=25000
# initPar=c(15,10, 1, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta   beta1    beta2  beta3   beta4  beta5    phi      LL  subLL   GPSLL priorLL LLSE
# newRow 17.511    8.5407 0.15173 -0.96288 2.0067 -4.7649 4.4382 107.64 -646.51 -359.1 -287.41       0    0
# > fitTwo21k5N$optGrad
# [1] -0.0019511721  0.0049398012  0.0021155045  0.0038308683  0.0053781487  0.0040767852 -0.0123335494 -0.0006049244
# > fitTwo21k5N$optPar
# [1]  17.5107121   8.5407163   0.1517300  -0.9628848   2.0066847  -4.7649010   4.4381949 107.6359089
# subsidences in south look better, good amount of taper variation.  Still too much uplift in middle

### <20km depth, dStar=25000
# initPar=c(15,10, 1, rep(0, nKnots-1), 175)
# looks just like how it did before using constrOptim (since lambda never got too large in this model)

### throw out GPS data, dStar=25000 (final fit, reltol=1e-13)
# initPar=c(15,10, 1, rep(0, nKnots-1), 175)
#        muZeta sigmaZeta  beta1   beta2  beta3   beta4    beta5        phi      LL   subLL   GPSLL priorLL LLSE
# newRow 18.776    8.1792 4.9302 -9.2719 1.1842 -4.6474 -0.23904 1.7878e-06 -283.92 -280.95 -2.9657       0    0
# looks great, but parameter SEs undefined since hessian is singular

### fitting seperate taper for GPS locking:
# dStar=25000, depthThresh=21000
#        muZeta sigmaZeta  beta1   beta2   beta3   beta4   beta5 beta'1  beta'2 beta'3  beta'4  beta'5    phi      LL   subLL   GPSLL priorLL LLSE
# newRow 16.659     8.503 6.4756 -11.862 0.35207 -6.7466 -1.3266 6.2185 -14.354  6.624 -10.821 0.10078 129.92 -619.47 -296.26 -323.21       0    0
# > fitTwo21k5N$optPar
# [1]  16.6588421   8.5029908   6.4756429 -11.8615043   0.3520656  -6.7465646  -1.3265508   6.2185023 -14.3541668
# [10]   6.6240185 -10.8209327   0.1007761 129.9181050
# solvable hessian:
# > sqrt(diag(solve(-fitTwo21k5N$hess)))
# [1] 0.319614553 0.138211858 0.094384816 0.240452501 0.007534719 0.166993043 0.034014780 0.139482193 0.343755171
# [10] 0.237865019 0.278334616 0.005274246 2.328260466
# > fitTwo21k5N$optGrad
# [,1]         [,2]         [,3]         [,4]          [,5]          [,6]         [,7]          [,8]
# [1,] -1.551975e-05 3.561577e-05 7.314052e-05 -2.17231e-06 -7.951076e-05 -2.685484e-05 0.0001637986 -0.0002109157
# [,9]         [,10]         [,11]         [,12]         [,13]
# [1,] -9.461934e-06 -3.444455e-05 -2.836473e-05 -0.0001713171 -7.972185e-06

# dStar=30000, depthThresh=25000
# muZeta sigmaZeta  beta1   beta2  beta3   beta4   beta5 beta'1  beta'2 beta'3  beta'4   beta'5    phi      LL
# newRow 12.239     10.98 7.3593 -13.051 0.2848 -7.3295 -1.3703 7.3576 -15.882 6.3633 -12.229 0.059539 172.91 -575.68
#   subLL   GPSLL priorLL LLSE
# -305.24 -270.44       0    0
# > fitTwo21k5N$optPar
# [1]  12.23922965  10.97972655   7.35926481 -13.05137878   0.28479619  -7.32950803  -1.37030017   7.35756168 -15.88177493
# [10]   6.36327637 -12.22941627   0.05953882 172.90810594
# > fitTwo21k5N$optGrad
# [,1]         [,2]          [,3]          [,4]          [,5]         [,6]         [,7]         [,8]
# [1,] 4.081068e-05 6.228362e-05 -0.0003448118 -0.0002510627 -0.0002670846 -0.000308497 2.334847e-05 -0.000205005
# [,9]         [,10]         [,11]        [,12]         [,13]
# [1,] -1.25745e-05 -0.0001531681 -0.0002325278 0.0001081064 -1.029363e-06
# > sqrt(diag(solve(-fitTwo21k5N$hess)))
# [1] 0.129772685         NaN 0.038976831         NaN         NaN 0.051574127         NaN 0.086992206 0.177649780
# [10]         NaN         NaN 0.001755541 3.338553733

# unnormalized taper, 30km depthThresh
# muZeta sigmaZeta      beta1       beta2      beta3       beta4       beta5     beta'1      beta'2     beta'3
# newRow 7.1347     13.67 0.00022586 -0.00039608 7.5186e-06 -0.00021186 -3.0821e-05 0.00023675 -0.00032791 -0.0001729
#      beta'4      beta'5    phi      LL   subLL   GPSLL priorLL LLSE
# -8.5986e-05 -8.4484e-05 277.24 -528.37 -317.52 -210.85       0    0
# > fitTwo21k5N$optPar
# [1]  7.134677e+00  1.367034e+01  2.258622e-04 -3.960813e-04  7.518637e-06 -2.118641e-04 -3.082099e-05  2.367532e-04
# [9] -3.279103e-04 -1.729013e-04 -8.598647e-05 -8.448366e-05  2.772408e+02
# > fitTwo21k5N$optGrad
# [,1]         [,2]      [,3]      [,4]      [,5]      [,6]     [,7]     [,8]     [,9]     [,10]     [,11]
# [1,] -0.0003568282 0.0001111812 -14.66225 -5.128251 -1.792358 -1.827856 -3.93076 16.48993 8.559512 0.5567106 0.4588311
# [,12]         [,13]
# [1,] 5.207347 -8.865153e-06
# > sqrt(diag(solve(-fitTwo21k5N$hess)))
# [1]          NaN 1.124518e-01 1.540353e-06          NaN          NaN          NaN 3.478149e-07          NaN 1.865825e-06
# [10]          NaN 1.276251e-06          NaN          NaN


# highInflate=1.00
# lowInflate=1.75
# highQual = as.numeric(dr1$quality) == 1
# lowQual = as.numeric(dr1$quality) != 1
# inflateDr1 = dr1
# inflateDr1$Uncertainty[highQual] = dr1$Uncertainty[highQual]*highInflate
# inflateDr1$Uncertainty[lowQual] = dr1$Uncertainty[lowQual]*lowInflate
depthThresh=15000
dStar=25000
nKnots=5
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
latThresh = 41
threshSlipDat = threshSlipDat[threshSlipDat$lat>latThresh,]
minLat = min(c(csz$latitude, threshSlipDat$lat)) - .001
maxLat = max(c(csz$latitude, threshSlipDat$lat)) + .001
par(mfrow=c(1,1))
quilt.plot(threshSlipDat$lon, threshSlipDat$lat, threshSlipDat$slip)
US(add=TRUE)
# initPar=c(14.637, 13.364, 1.6243, -3.6334, -2.67, -1.2649, -1.6643, 135.34) # this is under old basis
initPar=c(14.637, 13.364, -0.0400, 0.3994, -1.0057, -1.9691, 1.6643, 135.34)
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=inflateDr1, fault=csz, latRange=c(minLat, maxLat), 
                        normalModel=TRUE, normalizeTaper=TRUE, corGPS=TRUE)

# throw out GPS data:
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
threshSlipDat = threshSlipDat[sample(1:nrow(threshSlipDat), 1),]
minLat = min(c(csz$latitude, threshSlipDat$lat)) - .001
maxLat = max(c(csz$latitude, threshSlipDat$lat)) + .001
nKnots=5
initPar=c(14.637, 13.364, -0.0400, 0.3994, -1.0057, -1.9691, 1.6643, 135.34)
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=inflateDr1, fault=csz, latRange=c(minLat, maxLat), 
                        normalModel=TRUE, normalizeTaper=TRUE, corGPS=TRUE, doHess=TRUE)

# throw out subsidence data:
depthThresh=21000
dStar=25000
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
testSubDat = inflateDr1
testSubDat = testSubDat[sample(1:nrow(inflateDr1), 50),]
testG = okadaAll(csz, lonGrid, latGrid, cbind(testSubDat$Lon, testSubDat$Lat), slip=1, poisson=0.25)
initPar=c(15,10, 1, rep(0, nKnots-1), 175)
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=testG, 
                        fauxG=fauxG, subDat=testSubDat, fault=csz, 
                        normalModel=TRUE, normalizeTaper=TRUE, corGPS=TRUE)

# inflate gps errors by a factor of 3:
depthThresh=21000
dStar=25000
nKnots=5
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
threshSlipDat$slipErr = threshSlipDat$slipErr*3
# latThresh = 41
# threshSlipDat = threshSlipDat[threshSlipDat$lat>latThresh,]
testSubDat = inflateDr1
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=testSubDat, fault=csz, latRange=c(minLat, maxLat), 
                        normalModel=TRUE, normalizeTaper=TRUE, corGPS=TRUE, doHess=TRUE)

# fit model allowing GPS data to have extra taper basis, multiplying its uncertainty by 3
depthThresh=25000
dStar=30000
nKnots=5
nKnotsGPS=5
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
threshSlipDat$slipErr = threshSlipDat$slipErr*3
par(mfrow=c(1,1))
quilt.plot(threshSlipDat$lon, threshSlipDat$lat, threshSlipDat$slip)
US(add=TRUE)
# latThresh = 41
# threshSlipDat = threshSlipDat[threshSlipDat$lat>latThresh,]
testSubDat = inflateDr1
minLat = min(c(csz$latitude, threshSlipDat$lat)) - .001
maxLat = max(c(csz$latitude, threshSlipDat$lat)) + .001
initPar=c(15,10, 1, rep(0, nKnots-1), rep(0, nKnotsGPS), 175)
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=testSubDat, fault=csz, latRange=c(minLat, maxLat), 
                        normalModel=TRUE, normalizeTaper=TRUE, corGPS=TRUE, doHess=TRUE, 
                        diffGPSTaper=TRUE, nKnotsGPS=nKnotsGPS)

# try normalize=FALSE and redo the above code:
initPar=c(20,15, 1/20000, rep(0, nKnots-1), rep(0, nKnotsGPS), 175)
depthThresh=30000
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
threshSlipDat$slipErr = threshSlipDat$slipErr*3
minLat = min(c(csz$latitude, threshSlipDat$lat)) - .001
maxLat = max(c(csz$latitude, threshSlipDat$lat)) + .001
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=testSubDat, fault=csz, latRange=c(minLat, maxLat), 
                        normalModel=TRUE, normalizeTaper=FALSE, corGPS=TRUE, doHess=TRUE, 
                        diffGPSTaper=TRUE, nKnotsGPS=nKnotsGPS)

# compare subsidence and GPS tapers:
params = fitTwo21k5N$MLEs
splinePar = params[6:(5+nKnots)]
splineParGPS = params[(5+nKnots+1):(length(params)-1)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, latRange=c(minLat, maxLat))
GPSMat1 = getSplineBasis(csz, nKnots=nKnots, latRange=c(minLat, maxLat))
GPSMat2 = getSplineBasis(csz, nKnots=nKnotsGPS, latRange=c(minLat, maxLat))
lambdaGPS = GPSMat1 %*% splinePar - GPSMat2 %*% splineParGPS
tvecGPS = taper(getFaultCenters(csz)[,3], lambda=lambdaGPS, alpha=2, normalize=TRUE, dStar=dStar)
par(mfrow=c(1,2))
plotFault(csz, tvecGPS, main="GPS taper", varRange=c(0,1))
plotFault(csz, tvec, main="Subsidence taper", varRange=c(0,1))

# look at results:
params = fitTwo21k5N$MLEs
splinePar = params[6:(5+nKnots)]
splineParGPS = params[(5+nKnots):(length(params)-1)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, latRange=c(minLat, maxLat), normalize=FALSE)
GPSMat1 = getSplineBasis(csz, splinePar, nKnots=nKnots, latRange=c(minLat, maxLat))
plotFault(csz, tvec)
ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="Marginal ", nsim = 5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("GPSx3DiffTaper30km25km"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE)

## make ggplot plots:
fit = fitTwo21k5N
isT1 = events=="T1"
# T1DatRange = dr1[isT1 & inRangeDat,]
# GT1Range = G[isT1 & inRangeDat, inRangeFault]
# T1Dat = dr1[isT1,]
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

params = fit$MLEs
splinePar = params[6:(5+nKnots)]
Xi = getSplineBasis(csz, c(minLat,maxLat), nKnots)
lambdas = Xi %*% splinePar
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, normalize=TRUE)
plotFault(csz, tvec)

normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=FALSE, 
                                   dStar=dStar, gpsDat=threshSlipDat)

# areal values of zeta
muAreal = normalPreds$zetaEsts * tvec
sdAreal = normalPreds$zetaSD * tvec
medAreal = normalPreds$zetaMed * tvec
l95Areal = normalPreds$zeta025 * tvec
u95Areal = normalPreds$zeta975 * tvec

# get simulations
tab <- normalPreds$predResults
zetaSims = tab$zeta
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
ggplotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1 ", logScale=FALSE, 
                fileNameRoot=paste0("GPSx3DiffTaper30km25kmT1"))
ggComparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("GPSx3DiffTaper30km25kmT1"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE)

# how to switch between bases for flipped intercept basis fxn
nKnots=5
lats = csz$latitude
latRange=c(40,50)
intercept=FALSE

splineMat1 = cbind(1, bs(lats, df=nKnots-1, intercept=intercept, Boundary.knots=latRange))
splineMat2 = cbind(1, bs(latRange[1] + latRange[2] - lats, df=nKnots-1, intercept=intercept, Boundary.knots=latRange))

par1=c(1.6243, -3.6334, -2.67, -1.2649, -1.6643)
ys1 = splineMat1 %*% par1
mod = lm(ys1 ~ splineMat2-1)
par2 = coef(mod)
ys2 = splineMat2 %*% par2
sortI = sort(lats, index.return=TRUE)$ix

plot(lats[sortI], ys1[sortI], type="l", col="black")
lines(lats[sortI], ys2[sortI], col="red")

initPar=c(14.637, 13.364, 1.6243, -3.6334, -2.67, -1.2649, -1.6643, 135.34)
initPar=c(14.637, 13.364, -0.0400, 0.3994, -1.0057, -1.9691, 1.6643, 135.34)

# how to switch between bases for different numbers of knots
nKnots=5
nNewKnots=7
splineMat1 = getSplineBasis(csz, c(minLat, maxLat), nKnots)
splineMat2 = getSplineBasis(csz, c(minLat, maxLat), nNewKnots)

par1=c(-0.0400, 0.3994, -1.0057, -1.9691, 1.6643)
ys1 = splineMat1 %*% par1
mod = lm(ys1 ~ splineMat2-1)
par2 = coef(mod)
ys2 = splineMat2 %*% par2
sortI = sort(lats, index.return=TRUE)$ix

par(mfrow=c(1,1))
plot(lats[sortI], ys1[sortI], type="l", col="black")
lines(lats[sortI], ys2[sortI], col="red")

par2

# nNewKnots=7:
# initPar = c(14.637, 13.364, -0.04000000, 0.17190035, 0.09304813, -0.96061037, -1.78129373, 0.22882262, 1.66430000, 135.34)

# look at adding extra southern basis element that is thinner to accomodate sharp change in south
nKnots=5
splineMat1 = getSplineBasis(csz, c(minLat, maxLat), nKnots)
splineMat2 = getSplineBasis(csz, c(minLat, maxLat), nKnots, southernFun=TRUE)
par(mfrow=c(1,1))
sortI = sort(csz$latitude, index.return=TRUE)$ix
matplot(csz$latitude[sortI], splineMat2[sortI,], type="l")

# calculate uncertainty of parameters assuming fixed correlation length parameter, phi:
p = nrow(fitTwo21k5N$hess)
Sigma = solve(-fitTwo21k5N$hess[-p,-p])
SigmaUU = Sigma[1:(p-1), 1:(p-1)]
sigmaSqV = Sigma[p,p]
SigmaUV = Sigma[1:(p-1), p]
newSigma = SigmaUU - outer(SigmaUV, SigmaUV, "*")/sigmaSqV
origSEs = sqrt(diag(Sigma))
newSEs = c(sqrt(diag(newSigma)), 0)
ests = fitTwo21k5N$optPar
tab = rbind(ests, origSEs, newSEs)
colnames(tab) = c("mu", "sigma", paste0("beta", 1:nKnots), "phi")
rownames(tab) = c("MLEs", "original SEs", "new SEs")

#####
# test positive normal but with a naive mean adjustment
library(tmvtnorm)
params = fitTwo21k5N$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:(5+nKnots)]
phiZeta = params[length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)
corMatCSZ = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
                           onlyUpper=FALSE, distMat=rdist.earth(cbind(csz$longitude, csz$latitude), miles=FALSE), smoothness=3/2)
covMatCSZ = corMatCSZ * sigmaZeta^2
slipCov = sweep(sweep(covMatCSZ, 1, tvec, "*"), 2, tvec, "*")
n=nrow(csz)
n=10 # 3 seconds for n=10, 48 seconds for n=20 (==> O(n^4))
meanTest = rep(muZeta, n)
sigmaTest = covMatCSZ[1:n,1:n]
system.time(out <- mtmvnorm(mean=meanTest, sigma=sigmaTest, lower=rep(0, n), doComputeVariance = TRUE))

# now just compute the mean
n=20 # .3 seconds for n=10, 3 seconds for n=20, 6 for n=30, 11 for n=40, 18 for n=50, 27 for n=60, 51 for n=80, 86 for n=100, 680 for n=240
# ==> O(n^2) to O(n^3)
meanTest = rep(muZeta, n)
sigmaTest = covMatCSZ[1:n,1:n]
system.time(out <- mtmvnorm(mean=meanTest, sigma=sigmaTest, lower=rep(0, n), doComputeVariance = FALSE))
mean(out$tmean)

getPosNormMuN = function(muZeta, covMatCSZ, n=10, newMuInit=muZeta) {
  sigmaTest = covMatCSZ[1:n,1:n]
  
  newMu = newMuInit
  muDiff = Inf
  while(muDiff > .01) {
    meanTest = rep(newMu, n)
    out <- mtmvnorm(mean=meanTest, sigma=sigmaTest, lower=rep(0, n), doComputeVariance = FALSE)
    adjustedMu = mean(out$tmean)
    muDiff = adjustedMu - muZeta
    newMu = newMu - muDiff
    
    print(paste0("muDiff: ", muDiff, "; newMu: ", newMu))
  }
  newMu
}

getPosNormMu = function(muZeta, covMatCSZ, startN=20) {
  n = startN
  maxN = nrow(covMatCSZ)
  newMu = muZeta
  
  while(n < maxN) {
    newMu = getPosNormMuN(muZeta, covMatCSZ, n, newMuInit=newMu)
    n = n*2
    if(n > maxN)
      n = maxN
  }
  
  # one last fit at maxN:
  newMu = getPosNormMuN(muZeta, covMatCSZ, maxN, newMuInit=newMu)
  
  newMu
}



