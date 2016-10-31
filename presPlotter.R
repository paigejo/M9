# this script generates the plots for the results in the presentation

## first plot the data

# plot the GPS data with the fault geometry
par(mfrow=c(2,2))
quilt.plot(lon, lat, slip, nx=150, ny=150, main="Locking Rate (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
plotFault(faultGeom, plotData=FALSE, new=FALSE)
quilt.plot(lon, lat, Depth, nx=150, ny=150, main="CSZ Fault Depth (m)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
plotFault(faultGeom, plotData=FALSE, new=FALSE)
quilt.plot(lon, lat, slipErr, nx=150, ny=150, main="Locking Rate SD (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
plotFault(faultGeom, plotData=FALSE, new=FALSE)

# plot just GPS data locking rate with standard error
par(mfrow=c(1,2))
quilt.plot(lon, lat, slip, nx=150, ny=150, main="Locking Rate (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
quilt.plot(lon, lat, slipErr, nx=150, ny=150, main="Locking Rate SD (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)

# plot the subsidence data with the fault geomtry
lonRange=c(-127, -122.5)
latRange=c(40, 50)
par(mfrow=c(1,2))
plot(dr1$Lon, dr1$Lat, main="Data Locations", pch="+", col="red", 
     xlim=lonRange, ylim=latRange, xlab="Longitude", ylab="Latitude", cex=.7)
map("world", "Canada", add=TRUE)
US(add=TRUE)
# plotFault(faultGeom, plotData=FALSE, new=FALSE)
plot(dr1$subsidence, dr1$Lat, ylim=latRange, main="Estimated Subsidence (m)", 
     cex=.3, pch=19, col="blue", xlab="Subsidence (m)", ylab="Latitude")

# plot the fault geometry
par(mfrow=c(1,1))
lonRange=c(-127, -122.5)
latRange=c(40, 50)
plotFault(csz, plotData=FALSE)
map("world", "Canada", add=TRUE)
US(add=TRUE)

# precompute G
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)

# first load the data
out = load("FullIterFitPrior.RData")
muMatCSZ = state$muMatCSZ
muMatGPS = state$muMatGPS
parMat = state$parMat

## first look at stationary model
params = parMat[,1]
lambda = params[1]
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

comparePredsToSubs(params, G=G, plotNameRoot="Stationary", savePlots=FALSE)
comparePredsToGPS(params)

# generate a few realizations from the marginal distribution:
reals = preds(params, nsim=4, fault=csz)
reals = reals$slipSims
par(mfrow=c(2,2))
for(i in 1:4) {
  mag=getMomentFromSlip(reals[,i])
  plotFault(csz, reals[,i], 
            main=paste0("Sample magnitude ", round(mag, 2), " coseismic slips (m)"))
  map("world", "Canada", add=TRUE)
  US(add=TRUE)
}

# test the nonstationary model:
iter=4
muVecGPS = muMatGPS[,iter-1]
muVecCSZ = muMatCSZ[,iter-1]
muVec=c(muVecGPS, muVecCSZ)
params = parMat[,iter]
comparePredsToGPS(params, muVec=muVec)
comparePredsToSubs(params, G=G, plotNameRoot="Nonstationary", savePlots=FALSE, muVec=muVec)

# plot realizations from the marginal distribution:
reals = preds(params, nsim=4, fault=csz, muVec=muVec)
reals = reals$slipSims
par(mfrow=c(2,2))
for(i in 1:4) {
  mag=getMomentFromSlip(reals[,i])
  plotFault(csz, reals[,i], 
            main=paste0("Sample magnitude ", round(mag, 2), " coseismic slips (m)"))
  map("world", "Canada", add=TRUE)
  US(add=TRUE)
}

# test spline fit model
load("splineFit.RData")
params = splineFit$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]
splinePar = splineFit$splineParMLE
nKnots=4
dStar=26000
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)

comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "spline")       
comparePredsToGPS(params)

# generate a few realizations from the marginal distribution:
reals = preds(params, nsim=4, fault=csz, tvec=tvec)
reals = reals$slipSims
par(mfrow=c(2,2))
varRange=c(.1, max(reals))
for(i in 1:4) {
  mag=getMomentFromSlip(reals[,i])
  plotFault(csz, reals[,i], logScale = FALSE, varRange=varRange, 
            main=paste0("Sample slip (m), Mw=", round(mag, 2)), 
            xlim=c(-127.8, -122.5))
  map("world", "Canada", add=TRUE)
  US(add=TRUE)
}

# get magnitudes
par(mfrow=c(1,1))
reals = preds(params, nsim=25000, fault=csz, tvec=tvec)
reals = reals$slipSims
mags = apply(reals, 2, getMomentFromSlip)
hist(mags, main="Earthquake Magnitude", breaks=50, freq=F, xlab="Magnitude")
mean(mags > 9.3)
mean(mags > 9.2)
mean(mags > 9)

# look at historic quakes
load("historicQuakes.RData")
muMat = historicQuakes$muMat
sdMat = sqrt(historicQuakes$varMat)
minObs = 5
numObs = table(dr1$event)
threshUniqueEvents = uniqueEvents[numObs >= minObs]

# something's wrong with the matrices, so recompute them using Stan results
sdMat = matrix(nrow=nrow(muMat), ncol=ncol(muMat))
muMat = sdMat
for(i in 1:length(threshUniqueEvents)) {
  tab <- extract(historicQuakes$allStanResults[[i]], permuted = TRUE)
  zetaTab = tab$zeta
  sdMat[,i] = apply(zetaTab, 2, sd)
  muMat[,i] = apply(zetaTab, 2, mean)
}
muMatTaper = sweep(muMat, 1, tvec, "*")
sdMatTaper = sweep(sdMat, 1, tvec, "*")

# Only use events that have at least 5 observations (leaves out T12, T11, and T9a)

par(mfrow=c(1,2))
for(i in 1:length(threshUniqueEvents)) {
  thisEvent = threshUniqueEvents[i]
  plotFault(csz, muMat[,i]*tvec, main=paste0("Inferred ", thisEvent, " Mean"), logScale=F)
  plotFault(csz, sdMat[,i]*tvec, main=paste0("Inferred ", thisEvent, " SD"), logScale=F)
}

# consider their magnitudes
par(mfrow=c(1,1))
for(i in 1:length(threshUniqueEvents)) {
  thisEvent = threshUniqueEvents[i]
  tab <- extract(historicQuakes$allStanResults[[i]], permuted = TRUE) # return a list of arrays 
  Mw = tab$Mw
  hist(Mw, main=paste0("Imputed ", thisEvent, " Magnitude"), xlab="Magnitude", freq=F, breaks=30)
}


# plot Matern correlation:
ds = seq(0, 4, l=100)
par(mfrow=c(1,1))
plot(ds, Matern(ds, nu=3/2, range=1), type="l", main="Matern Correlation", 
     xlab="Distance", ylab="Correlation", ylim=c(0,1))
