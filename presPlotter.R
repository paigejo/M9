# this script generates the plots for the results in the presentation

## first plot the data

# plot just GPS locking rate data product with standard error
par(mfrow=c(1,2))
quilt.plot(lon, lat, slip, nx=100, ny=100, main="Locking Rate (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)
quilt.plot(lon, lat, slipErr, nx=100, ny=100, main="Locking Rate SD (mm/yr)")
map("world", "Canada", add=TRUE)
US(add=TRUE)

# plot the subsidence data with the coarse fault geometry
par(mfrow=c(1,2))

sites = unique(dr1$Site)
siteLats = aggregate(dr1$Lat, list(dr1$Site), mean)[,2]
sortI = sort(siteLats, index.return=TRUE)$ix
sites = sites[sortI]
# cols = tim.colors(length(sites))
cols = rainbow(length(sites))
lonRange=c(-127, -122.5)
latRange=c(40, 50)
for(i in 1:length(sites)) {
  thisSite = sites[i]
  siteDat = dr1[dr1$Site == thisSite,]
  
  if(i == 1) {
    plot(siteDat$Lon, siteDat$Lat, pch="+", col=cols[i], xlab="Longitude", ylab="Latitude", main="Subsidence sites", 
         xlim=lonRange, ylim=latRange)
  }
  else {
    points(siteDat$Lon, siteDat$Lat, pch="+", col=cols[i])
  }
}
map("world", "Canada", add=TRUE)
US(add=TRUE)
plotFault(faultGeom, plotData=FALSE, new=FALSE)

for(i in 1:length(sites)) {
  thisSite = sites[i]
  siteDat = dr1[dr1$Site == thisSite,]
  
  if(i == 1) {
    plot(siteDat$subsidence, siteDat$Lat, pch="+", col=cols[i], xlab="Longitude", ylab="Latitude", main="Subsidence (m)", 
         ylim=latRange)
  }
  else {
    points(siteDat$subsidence, siteDat$Lat, pch="+", col=cols[i])
  }
}

########
# plot the fault geometry
par(mfrow=c(1,1))
lonRange=c(-127, -122.5)
latRange=c(40, 50)
plotFault(csz, plotData=FALSE)
map("world", "Canada", add=TRUE)
US(add=TRUE)

# plot the taper function
library(latex2exp)
par(mfrow=c(1,1))
ds = seq(0, 1, l=100)
plot(ds, taper(ds, lambda=1, dStar=1), type="l", main="Taper", 
     xlab=TeX("Depth ($d/d^*$)"), ylab="Taper")
lines(ds, taper(ds, lambda=2, dStar=1), col="red")
lines(ds, taper(ds, lambda=3, dStar=1), col="blue")
lines(ds, taper(ds, lambda=5, dStar=1), col="purple")
legend("topright", c(TeX("$\\lambda = 1$"), TeX("$\\lambda = 2$"), TeX("$\\lambda = 3$"), 
                     TeX("$\\lambda = 5$")), col=c("black", "red", "blue", "purple"), 
       lty=1)

# plot spline basis
nKnots=4
lats = seq(40, 50, l=100)
XiLat = bs(lats, df=nKnots, intercept=FALSE, Boundary.knots=latRange)
XiLat = cbind(rep(1, nrow(XiLat)), XiLat)
matplot(lats, XiLat, main="B-spline basis", xlab="Latitude", type="l", lty=1, 
        pch=19, cex=.4, lwd=4, ylab="B-spline basis", ylim=c(0, 1.1))

# plot why not to use MSE example
lonRange=c(-127, -122.5)
latRange=c(40, 50)
par(mfrow=c(1,2))
plotFault(csz, plotData=FALSE, ylim=latRange, main="Cascadia Fault")
plot(dr1$Lon, dr1$Lat, main="Data Locations", pch="+", col="red", 
     xlim=lonRange, ylim=latRange, xlab="Longitude", ylab="", cex=.7)
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

for(i in 1:length(threshUniqueEvents)) {
  thisEvent = threshUniqueEvents[i]
  pdf(paste0("infered", thisEvent, "vsSub.pdf"), width=9, height=4.5)
  par(mfrow=c(1,3))
  plotFault(csz, muMat[,i]*tvec, main=paste0("Inferred ", thisEvent, " Mean"), logScale=F, 
            ylim=c(40,50), xlim=c(-128, -122.5))
  map("world", "Canada", add=TRUE)
  US(add=TRUE)
  plotFault(csz, sdMat[,i]*tvec, main=paste0("Inferred ", thisEvent, " SD"), logScale=F, 
            ylab="", ylim=c(40,50), xlim=c(-128, -122.5))
  map("world", "Canada", add=TRUE)
  US(add=TRUE)
  plot(dr1[events == thisEvent,"subsidence"], dr1[events == thisEvent,"Lat"], 
       main="Observed Subsidence (m)", ylim=c(40,50), xlab="Subsidence (m)", 
       ylab="")
  dev.off()
}

for(i in 1:length(threshUniqueEvents)) {
  thisEvent = threshUniqueEvents[i]
  pdf(paste0("infered", thisEvent, ".pdf"), width=9, height=8)
  par(mfrow=c(1,2))
  plotFault(csz, muMat[,i]*tvec, main=paste0("Inferred ", thisEvent, " Mean"), logScale=F)
  plotFault(csz, sdMat[,i]*tvec, main=paste0("Inferred ", thisEvent, " SD"), logScale=F)
  dev.off()
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

# make table of coefficient estimates and standard errors:
library(xtable)
finalEsts = c(232.5722, params[c(2, 3, 5:9)])
covMat = hessianToCovMat(splineFit$hess)
SEs = c(76.47129, sqrt(diag(covMat)))
resultsTab = rbind(finalEsts, SEs)
colnames(resultsTab) = c("phi", "mu_zeta", "sigma_zeta", "beta_0", "beta_1", "beta_2", 
                        "beta_3", "beta_4")
rownames(resultsTab) = c("Estimates", "SE")
xtable(resultsTab, digits=2)
corMat = hessianToCorrMat(splineFit$hess)
splineCorMat = corMat[3:7, 3:7]
rownames(splineCorMat) = c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4")
colnames(splineCorMat) = c("beta_0", "beta_1", "beta_2", "beta_3", "beta_4")
print(splineCorMat, digits=2)
xtable(splineCorMat, digits=2)

splineCovMat = covMat[3:7, 3:7]
plotSplineUncertainty(splinePar, splineCovMat)

# Plot an example realization of zeta:
par(mfrow=c(1,1))
reals = preds(params, nsim=2, fault=csz, tvec=rep(1, nrow(csz)))
plotFault(csz, reals$slipSims[,1], main=TeX("Sample Realization of $\\zeta$"))
map("world", "Canada", add=TRUE)
US(add=TRUE)

# ASL fit
load("splineFitASL.RData")
params = splineFitASL$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]
splinePar = splineFitASL$splineParMLE
nKnots=4
dStar=26000
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=FALSE, tvec=tvec, nsim=1000, fileNameRoot = "splineASL") 
