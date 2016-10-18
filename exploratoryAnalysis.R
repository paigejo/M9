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
Xi = bs(dr1$Lat, df=5, intercept=TRUE)

# test Xi
par(mfrow=c(1,1))
matplot(dr1$Lat, Xi, main="Spline Basis", xlab="Latitude")

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
rstan_options(auto_write = TRUE)
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
out = load("fullIterFit.RData")
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
             zlim=c(0,35))
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
             zlim=c(0,80))
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
  plotFault(csz, taperedHiSlip, main=paste0("Slip 95th percentile at iteration ", i), varRange=c(0,80))
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
for(i in 1:length(state$optimTables)) {
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


