# plotter functions for the BIRS poster

# makes plots comparing the predicted subsidence levels to the subsidence data.
# if prediction inputs are left NULL, they are computed using params
margDistPlotter = function(params, slipPreds=NULL, slipPredsGPS=NULL, subPreds=NULL, 
                           subPredsGPS=NULL, nsim=100, plotNameRoot="full", 
                           savePlots=TRUE, G=NULL, fileNameRoot=plotNameRoot, 
                           muVec=NULL, useGPS=FALSE, tvec=NULL, subDat=dr1, 
                           logScale=FALSE, fault=csz, latRange=c(40, 50), 
                           posNormalModel=FALSE, normalModel=posNormalModel, doGPSPred=FALSE, 
                           useMVNApprox=FALSE) {
  # get parameters
  if(is.null(muVec)) {
    lambdaMLE = params[1]
    muZetaMLE = params[2]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
    muZetaCSZ = rep(muZetaMLE, nrow(fault))
  }
  else {
    lambdaMLE = params[1]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaMLE = muVec
    muZetaGPS = muVec[1:nrow(slipDatCSZ)]
    muZetaCSZ = muVec[(nrow(slipDatCSZ)+1):length(muVec)]
  }
  
  # get taper if necessary
  if(is.null(tvec))
    tvec = taper(fault$depth, lambda=lambdaMLE)
  
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = preds(params, nsim=nsim, fault=fault, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec, 
                      posNormalModel=posNormalModel, normalModel=normalModel)
  if(is.null(slipPredsGPS) && doGPSPred)
    slipPredsGPS = predsGivenGPS(params, nsim=nsim, muVec=c(muZetaGPS, muZetaCSZ), fault=fault, tvec=tvec, 
                                 posNormalModel=posNormalModel, normalModel=normalModel)
  if(is.null(subPreds)) {
    if(is.null(G))
      subPreds = predsToSubsidence(params, slipPreds, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPreds = predsToSubsidence(params, slipPreds, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  if(is.null(subPredsGPS) && doGPSPred) {
    if(is.null(G))
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  meanSlip = slipPreds$meanSlip
  meanSub = subPreds$meanSub
  u95 = subPreds$u95
  l95 = subPreds$l95
  u95Noise = subPreds$u95Noise
  l95Noise = subPreds$l95Noise
  slipSD = apply(slipPreds$slipSims, 1, sd)
  if(doGPSPred) {
    meanSlipGPS = slipPredsGPS$meanSlip
    meanSubGPS = subPredsGPS$meanSub
    u95GPS = subPredsGPS$u95
    l95GPS = subPredsGPS$l95
    u95NoiseGPS = subPredsGPS$u95Noise
    l95NoiseGPS = subPredsGPS$l95Noise
  }
  
  ##### generate mean seaDef field from Okada model
  # set Okada subsidence grid
  lonRange=c(-128, -122.5)
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  coordGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
  meanSeaDef = okada(fault, lonGrid, latGrid, slips=meanSlip)
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubsidencePredictionsPoster.pdf"), width=8, height=10, bg="white")
  par(mfrow=c(2,2))
  if(!logScale) {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  }
  else {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, varRange=c(.1, max(meanSlip)), logScale=TRUE)
  }
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotFault(fault, plotVar = slipSD, legend.mar=6, main=paste0(plotNameRoot, "Slip SD (m)"), 
            xlim=lonRange, ylim=latRange)
  map("world", "Canada", add=TRUE, lwd=1)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  plotFault(csz, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, subDat$subsidence))
  # simulated subsidence data from Okada model using marginal distribution
  plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, "95% Prediction Band"), ylab="Latitude", 
       xlab="Subsidence (m)")
  points(-meanSub, subDat$Lat, pch=19, cex=.3, col="blue")
  ord = order(subDat$Lat)
  lines(-l95Noise[ord], subDat$Lat[ord])
  lines(-u95Noise[ord], subDat$Lat[ord])
  # plot magnitude distribution
  mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault)
  cleanMags = mags[is.finite(mags)]
  hist(cleanMags, breaks=50, main=paste0(plotNameRoot, "95% Magnitude CI"), xlab="Magnitude", 
       freq=F)
  xs = seq(min(cleanMags), max(cleanMags), l=100)
  lines(xs, dnorm(xs, mean(cleanMags), sd(cleanMags)), col="blue")
  abline(v=mean(cleanMags), col="black")
  abline(v=quantile(cleanMags, probs=.975), col="purple", lty=2)
  abline(v=quantile(cleanMags, probs=.025), col="purple", lty=2)
  
  if(savePlots)
    dev.off()
  
  invisible(NULL)
}

##### plot the data

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
    plot(siteDat$subsidence, siteDat$Lat, pch="+", col=cols[i], xlab="Subsidence", ylab="Latitude", main="Subsidence (m)", 
         ylim=latRange)
  }
  else {
    points(siteDat$subsidence, siteDat$Lat, pch="+", col=cols[i])
  }
}

# plot all data at once
par(mfrow=c(1,3))
quilt.plot(lon, lat, slip, nx=100, ny=100, main="Locking Rate (mm/yr)", ylab="Latitude", xlab="Longitude")
map("world", "Canada", add=TRUE)
US(add=TRUE)

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
    plot(siteDat$Lon, siteDat$Lat, pch="+", col=cols[i], xlab="Longitude", ylab="", main="Subsidence sites", 
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
    plot(siteDat$subsidence, siteDat$Lat, pch="+", col=cols[i], xlab="Subsidence", ylab="", main="Subsidence (m)", 
         ylim=latRange)
  }
  else {
    points(siteDat$subsidence, siteDat$Lat, pch="+", col=cols[i])
  }
}

##### plot the taper function
library(latex2exp)
par(mfrow=c(1,1))
ds = seq(0, 1, l=100)
plot(ds, taper(ds, lambda=1, dStar=1), type="l", main="Double-exponential taper", 
     xlab=TeX("Depth ($d/d^*$)"), ylab="")
lines(ds, taper(ds, lambda=2, dStar=1), col="red")
lines(ds, taper(ds, lambda=3, dStar=1), col="blue")
lines(ds, taper(ds, lambda=5, dStar=1), col="purple")
legend("topright", c(TeX("$\\lambda = 1$"), TeX("$\\lambda = 2$"), TeX("$\\lambda = 3$"),
                     TeX("$\\lambda = 5$")), col=c("black", "red", "blue", "purple"), lty=1)

##### plot spline basis
nKnots=5
lats = seq(min(csz$latitude), max(csz$latitude), l=100)
splineMat = getSplineBasis(csz, nKnots=nKnots, lats=lats)
matplot(lats, splineMat, type="l", lwd=2, lty=1, xlim=latRange, main="B-Spline Basis", xlab="Latitude", ylab="")


#####
# test predictions for normal and pos normal models

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

params = splineFit21k5N$MLEs
muZeta = params[2]
sigmaZeta = params[3]
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)

# get marginal predictions
margDistPlotter(params, G=G, tvec=tvec, plotNameRoot="Marginal ", nsim=2000, 
                subDat=dr1, logScale=FALSE, fileNameRoot=paste0("21k5Normal"), 
                fault=csz, normalModel=TRUE, useMVNApprox=FALSE)
margDistPlotter(params, G=G, tvec=tvec, plotNameRoot="Marginal ", nsim=2000, 
                subDat=dr1, logScale=FALSE, fileNameRoot=paste0("21k5PosNormal"), 
                fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE)

# get T1 predictions
isT1 = events=="T1"
# T1DatRange = dr1[isT1 & inRangeDat,]
# GT1Range = G[isT1 & inRangeDat, inRangeFault]
# T1Dat = dr1[isT1,]
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=2000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE)
posNormalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=2000, G=GT1, prior=FALSE, tvec=tvec, 
                                      normalModel=TRUE, posNormalModel=TRUE)

# areal values of zeta
muArealN = normalPreds$zetaEsts * tvec
muArealPN = posNormalPreds$zetaEsts * tvec

# get simulations
tab <- normalPreds$predResults
zetaSims = tab$zeta
slipSimsN = sweep(zetaSims, 1, tvec, "*")
tab <- posNormalPreds$predResults
zetaSims = tab$zeta
slipSimsPN = sweep(zetaSims, 1, tvec, "*")

slipPredsN = list(meanSlip=muArealN, slipSims=slipSimsN)
slipPredsPN = list(meanSlip=muArealPN, slipSims=slipSimsPN)

# plot results:
margDistPlotter(params, slipPreds=slipPredsN, G=GT1, tvec=tvec, plotNameRoot="1700 ", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("21k5NormalT1"), 
                   fault=csz, normalModel=TRUE, useMVNApprox=FALSE)
margDistPlotter(params, slipPreds=slipPredsPN, G=GT1, tvec=tvec, plotNameRoot="1700 ", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("21k5PosNormalT1"), 
                   fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE)

# do the same for lognormal
params = splineFit21k5LN$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz)

margDistPlotter(params, G=G, tvec=tvec, plotNameRoot="Marginal ", nsim=2000, 
                subDat=dr1, logScale=FALSE, fileNameRoot=paste0("21k5LogNormal"), 
                fault=csz, normalModel=FALSE, posNormalModel=FALSE, useMVNApprox=FALSE)

logNormalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=2000, G=GT1, prior=FALSE, tvec=tvec, 
                                      normalModel=FALSE, posNormalModel=FALSE)
muArealLN = logNormalPreds$zetaEsts * tvec
tab <- logNormalPreds$predResults
zetaSims = tab$zeta
slipSimsLN = sweep(zetaSims, 1, tvec, "*")
slipPredsLN = list(meanSlip=muArealLN, slipSims=slipSimsLN)

margDistPlotter(params, slipPreds=slipPredsLN, G=GT1, tvec=tvec, plotNameRoot="1700 ", 
                subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("21k5LogNormalT1"), 
                fault=csz, normalModel=FALSE, posNormalModel=FALSE, useMVNApprox=FALSE)


##### make the CV table
setwd("/Users/johnpaige/git/M9")
load("foldCV.RData")
load("marginalCV.RData")
fullTab = rbind(resTab, resTabSummary)
rownames(fullTab) = c(paste0("Marginal ", rownames(resTab)), 
                      paste0("1700 ", rownames(resTabSummary)))

# round the values
fullTab = matrix(sapply(fullTab, signif, digits=2), ncol=3)
colnames(fullTab) = colnames(resTab)
rownames(fullTab) = c(paste0("Marginal ", rownames(resTab)), 
                      paste0("1700 ", rownames(resTabSummary)))

# convert to latex
library(xtable)
xtable(fullTab)


