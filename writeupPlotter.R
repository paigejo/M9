# plots for writeup

library(ggmap)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(latex2exp)

setwd("~/git/M9/writeup/mathGeosci/")

##### plot the subsidence estimates and the fault geometry
### Set a range
# lon = c(-128.5, -121.5)
lon = c(-128, -122)
lat = c(39.5, 50.5)

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

bg = ggmap(map) + 
  scale_x_continuous(limits = lon, expand = c(0, 0)) +
  scale_y_continuous(limits = lat, expand = c(0, 0)) + 
  ggtitle("Subsidence Data Sites") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Longitude", y="Latitude")

mapPoints = geom_point(aes(x = Lon, y = Lat, color=factor(Site, levels=sites)), data=dr1, size=3, shape="+")

faultPoly = ggplotFault(faultGeom, color=rgb(.2,.2,.2))

bg + faultPoly + mapPoints + guides(color=FALSE)

##### Let's try with greyed out land:
library(maps)
library(mapdata)

# choose color if necessary (lightblue1 or white)
fields.color.picker()

# get relevant map data
states <- map_data("state")
west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
canada = map_data("world", "Canada")

fullMap = rbind(canada, west_coast)

# try plotting on an interpolated grid
library(gstat)
library(sp)
library(maptools)
library(RColorBrewer)

# plot it (choose background color with fields.color.picker())
bg = ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
  coord_fixed(xlim = lon,  ylim = lat, ratio = 1.3, expand=FALSE) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='lightblue1')) + 
  labs(x="Longitude", y="Latitude")

##### do the same thing but for subsetting locking rates to depths < 30km

# construct bounding polygon around data:
#make concave hull around prediction mask to get final prediction points
source("model1/seaDefAnis.r")
locking30km = slipDat
locking30km = locking30km[slipDat$Depth < 30000,]
ashape30km = ahull(locking30km$lon, locking30km$lat, alpha=2)
indx=ashape30km$arcs[,"end1"]  
hullPts <- cbind(locking30km$lon, locking30km$lat)[indx,]                  # extract the boundary points from maskXY

#plot hull and data to make sure it works
plotSubPoly(rbind(hullPts, hullPts[1,]), cbind(locking30km$lon, locking30km$lat))

# now subset prediction data frame to only be predictions within the polygon from alphahull
library(akima)

# interpolate between observations for maximum purdyness
predGrid = make.surface.grid(list(lon=seq(lon[1], lon[2], by=.01), lat=seq(lat[1], lat[2], by=.01)))
preds = interpp(locking30km$lon, locking30km$lat, locking30km$slip, predGrid[,1], predGrid[,2])
maskFinalPreds = in.poly(cbind(preds$x, preds$y), hullPts, convex.hull=FALSE)
preds = data.frame(preds)
finalPreds = preds[maskFinalPreds,]

p1 = bg + geom_tile(data = finalPreds, aes(x = x, y = y, fill = z)) + 
  scale_fill_distiller("", palette = "Spectral", direction=-1) + 
  geom_point(aes(x=lon, y=lat), pch=20, col="black", size=.1, data=locking30km) + 
  ggtitle("Locking Rates (mm/yr)") + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)
p1

##### now add plot of the standard errors
# must first interpolate standard errors in the same way as the locking rates
predsSE <- interpp(locking30km$lon, locking30km$lat, locking30km$slipErr, predGrid[,1], predGrid[,2])
predsSE = data.frame(predsSE)
finalPredsSE = predsSE[maskFinalPreds,]

p2 = bg + geom_tile(data = finalPredsSE, aes(x = x, y = y, fill = z)) + 
  scale_fill_distiller("", palette = "Spectral", direction=-1) + 
  geom_point(aes(x=lon, y=lat), pch=20, col="black", size=.1, data=locking30km) + 
  ggtitle("Locking Rate SEs (mm/yr)") + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)
p2

##### plot both at the same time!
p12 = multiplot(p1, p2, cols=2)

pdf("lockingRates.pdf", width=7)
p12 = multiplot(p1, p2, cols=2)
print(p12)
dev.off()

# try the same thing but without plotting the data points
p1 = bg + geom_tile(data = finalPreds, aes(x = x, y = y, fill = z)) + 
  scale_fill_distiller("", palette = "Spectral", direction=-1) + 
  ggtitle("Locking Rates (mm/yr)") + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)

p2 = bg + geom_tile(data = finalPredsSE, aes(x = x, y = y, fill = z)) + 
  scale_fill_distiller("", palette = "Spectral", direction=-1) + 
  ggtitle("Locking Rate SEs (mm/yr)") + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)

p12 = multiplot(p1, p2, cols=2)
p12

##### compare google maps plotting of sites with grey map plot equivalent

# google maps plot:
style1=c(feature="administrative", element="labels", visibility="off")
style2=c("&style=", feature="road", element="geometry", visibility="off")
style3=c("&style=", feature="poi", element="labels", visibility="off")
style4=c("&style=", feature="landscape", element="labels", visibility="off")
style5=c("&style=", feature="administrative", element="geometry.stroke", color="black")
style6=c("&style=", feature="administrative", element="geometry.stroke", weight=.75)
map <- get_googlemap(center=c(lon = mean(lon), lat = mean(lat)), zoom=5,
                     style=c(style1, style2, style3, style4, style5, style6))

# When you draw a figure, you limit lon and lat, scramble site labels by latitude so nearby sites 
# have very different colors
fauxObs = data.frame(list(Lon=dr1$Lon, Lat=dr1$Lat))
# goodScramble = c(19, 2, 7, 5, 3, 20, 12, 9, 4, 6, 16, 11, 17, 15, 21, 14, 10, 18, 13, 1, 8)
goodScramble = c(6, 8, 22, 5, 1, 3, 10, 12, 17, 18, 14, 13, 19, 16, 11, 21, 4, 2, 7, 23, 20, 15, 9)
sites = unique(dr1$Site)
siteLats = aggregate(dr1$Lat, list(dr1$Site), mean)[,2]
sortI = sort(siteLats, index.return=TRUE)$ix
sites = sites[sortI]
sites = sites[goodScramble]

bg = ggmap(map) + 
  coord_fixed(xlim = lon,  ylim = lat, ratio = 1.3, expand=FALSE) + 
  ggtitle("Subsidence Data Sites") + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(x="Longitude", y="Latitude")

mapPoints = geom_point(aes(x = Lon, y = Lat, color=factor(Site, levels=sites)), data=dr1, size=3, shape="+")

faultPoly = ggPlotFault(faultGeom, color=rgb(.2,.2,.2))

bg + faultPoly + mapPoints + guides(color=FALSE)

# grey maps plot:
bg = ggplot() + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
  coord_fixed(xlim = lon,  ylim = lat, ratio = 1.3, expand=FALSE) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='lightblue1')) + 
  labs(x="Longitude", y="Latitude") + 
  ggtitle("Subsidence Data Sites")

# bg + faultPoly + mapPoints + guides(color=FALSE)

bg + faultPoly + mapPoints + guides(color=FALSE) + scale_colour_hue(c=400, l=70)
# bg + faultPoly + mapPoints + guides(color=FALSE) + scale_fill_manual(values=hsv(seq(0, 1, l=23)))

# now add in subsidence data
ggplot() + geom_point(aes(x=dr1$subsidence, y=dr1$Lat), pch=19, color="blue", size=.5) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white')) + 
  labs(x="Subsidence", y="Latitude") + 
  ggtitle("Subsidence Estimates (m)") + 
  coord_fixed(ylim = lat, expand=FALSE)

# plot them together
p1 = bg + faultPoly + mapPoints + guides(color=FALSE) + scale_colour_hue(c=400, l=70)
p2 = ggplot() + geom_point(aes(x=dr1$subsidence, y=dr1$Lat), pch=19, color="blue", size=.5) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white', color="black"), aspect.ratio=2.5/1) + 
  labs(x="Subsidence", y="") + 
  ggtitle("Subsidence Estimates (m)") + 
  coord_fixed(ylim = lat, expand=FALSE)

multiplot(p1, p2, cols=2)

pdf("subsidenceData.pdf", width=7)
p12 = multiplot(p1, p2, cols=2)
print(p12)
dev.off()



##### taper plots
lambdas = c(.1, 1, 2, 5, 10)
ds = seq(0, 1, l=100)
taperMat = sapply(lambdas, taper, dStar=1, d=ds)

pdf("Taper.pdf", width=6)
matplot(ds, taperMat, xlab=TeX("Depth fraction (d/d^*)"), ylab="Taper", 
        main="", type="l", lty=1, col=rainbow(length(lambdas)))
legend("topright", c(TeX(paste0("$\\lambda = ", lambdas[1], "$")), 
                     TeX(sprintf("$\\lambda = %.0f$", lambdas[2:length(lambdas)]))), 
       lty=1, col=rainbow(length(lambdas)))
dev.off()

##### B-spline basis
lats = seq(40,50, l=100)
Xi = getSplineBasis(csz, c(40,50), nKnots=5, lats=lats)
pdf("BSplineBasis.pdf", width=6)
matplot(lats, Xi, xlab="Latitude", ylab="B-spline basis", 
        main="", type="l", lty=1, col=rainbow(ncol(Xi)), ylim=c(0, 1.2))
dev.off()


##### Make parameter table

# setup: SD inflation and setting depth thresholds
highInflate = 1.25
lowInflate = 1.75
highQual = as.numeric(dr1$quality) == 1
lowQual = as.numeric(dr1$quality) != 1
depthThresh = 21000
dStar = 25000
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]

# comparison model (no latitudal variation in taper)
initPar=c(20,15, 1, 175)
nKnots=1
fitTwo21k1N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=inflateDr1, fault=csz, 
                        normalModel=TRUE, normalizeTaper=TRUE)

# proposed model
initPar=c(20,15, 1, rep(0, nKnots-1), 175)
nKnots=5
fitTwo21k5N = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                        useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                        fauxG=fauxG, subDat=inflateDr1, fault=csz, 
                        normalModel=TRUE, normalizeTaper=TRUE)

params = fitTwo21k5N$MLEs
params = fitTwo21k1N$MLEs
params = fit$MLEs
MLEs = c(params[c(2, 3, 5, 6:(5+nKnots), length(params))], lowInflate, highInflate)
SEs = sqrt(diag(solve(-fitTwo21k5N$hess)))
SEs = sqrt(diag(solve(-fitTwo21k1N$hess)))
SEs = sqrt(diag(solve(-fit$hess)))
SEs = c(SEs[1:2], NA, SEs[3:length(SEs)], NA, NA)
tab = rbind(MLEs, SEs)
colnames(tab) = c("muzeta", "sigmazeta", "gamma", paste0("beta", 1:nKnots), "phi", "psil", "psih")
rownames(tab) = c("MLEs", "SEs")
library(xtable)
xtable(tab, digits=3)

splinePar = params[6:(length(params)-1)]
Xi = getSplineBasis(csz, c(40,50), nKnots)
lambdas = Xi %*% splinePar
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, normalize=TRUE)
plotFault(csz, tvec)

# plot marginal fits
ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="Marginal ", 
                   subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("finalN"), 
                   fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                   dStar=dStar, normalizeTaper=TRUE, nsim=10000)

ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="Marginal ", 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("finalPN"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=TRUE, posNormalModel=TRUE, nsim=10000)

# T1 data and predictions
isT1 = events=="T1"
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat)
posNormalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                      normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                      dStar=dStar, gpsDat=threshSlipDat)

# areal values of zeta
muArealN = normalPreds$zetaEsts * tvec
sdArealN = normalPreds$zetaSD * tvec
medArealN = normalPreds$zetaMed * tvec
l95ArealN = normalPreds$zeta025 * tvec
u95ArealN = normalPreds$zeta975 * tvec

# get simulations
tab <- normalPreds$predResults
zetaSimsN = tab$zeta
slipSimsN = sweep(zetaSimsN, 1, tvec, "*")
slipPredsN = list(meanSlip=muArealN, slipSims=slipSimsN)

# plot results:
ggplotFixedSlip(muArealN, NULL, l95ArealN, u95ArealN, sdArealN, event="T1", plotNameRoot="T1 ", logScale=FALSE, 
              fileNameRoot=paste0("finalNT1"))
ggComparePredsToSubs(params, slipPreds=slipPredsN, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("finalNT1"), 
                   fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                   dStar=dStar, normalizeTaper=FALSE, nsim=5000)

# areal values of zeta
muArealPN = posNormalPreds$zetaEsts * tvec
sdArealPN = posNormalPreds$zetaSD * tvec
medArealPN = posNormalPreds$zetaMed * tvec
l95ArealPN = posNormalPreds$zeta025 * tvec
u95ArealPN = posNormalPreds$zeta975 * tvec

# get simulations
tab <- posNormalPreds$predResults
zetaSimsPN = tab$zeta
slipSimsPN = sweep(zetaSimsPN, 1, tvec, "*")
slipPredsPN = list(meanSlip=muArealPN, slipSims=slipSimsPN)

ggplotFixedSlip(muArealPN, NULL, l95ArealPN, u95ArealPN, sdArealPN, event="T1", plotNameRoot="T1 ", logScale=FALSE, 
                fileNameRoot=paste0("finalPNT1"))
ggComparePredsToSubs(params, slipPreds=slipPredsPN, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                   subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("finalPNT1"), 
                   fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                   dStar=dStar, normalizeTaper=FALSE, posNormalModel=TRUE, nsim=5000)

##### Make a plot comparing subsidence predictions for the different models (marginal and 1700)

## first generate subsidences
# Normal marginal
slipPreds = preds(params, nsim=10000, fault=csz, tvec=tvec, 
                  posNormalModel=FALSE, normalModel=TRUE, phiZeta=params[length(params)])
subPreds1 = predsToSubsidence(params, slipPreds, G=G, useMVNApprox=FALSE, subDat=dr1, 
                             posNormalModel=FALSE, normalModel=TRUE, tvec=tvec)

# Positive Normal marginal
slipPreds = preds(params, nsim=10000, fault=csz, tvec=tvec, 
                  posNormalModel=TRUE, normalModel=TRUE, phiZeta=params[length(params)])
subPreds2 = predsToSubsidence(params, slipPreds, G=G, useMVNApprox=FALSE, subDat=dr1, 
                              posNormalModel=TRUE, normalModel=TRUE, tvec=tvec)

# Normal T1 predictive
subPreds3 = predsToSubsidence(params, slipPredsN, G=GT1, useMVNApprox=FALSE, subDat=T1Dat, 
                              posNormalModel=FALSE, normalModel=TRUE, tvec=tvec)

# Positive Normal T1 predictive
subPreds4 = predsToSubsidence(params, slipPredsPN, G=GT1, useMVNApprox=FALSE, subDat=T1Dat, 
                              posNormalModel=TRUE, normalModel=TRUE, tvec=tvec)

## now compare them with plot:
ggCompareSubs(params, subPreds1, subPreds2, subPreds3, subPreds4, dr1, dr1, T1Dat, T1Dat, 
              tvec, fileNameRoot="finalSubComparison", 
              plotNameRoot1="Normal Marginal ", plotNameRoot2="Pos. Normal Marginal ", 
              plotNameRoot3="Normal T1 (1700) ", plotNameRoot4="Pos. Normal T1 (1700) ", 
              noTitle=FALSE)

###################################
###################################
###################################
###################################
# compare tapers of different models

# first modify GPS model taper
fitGPS = fitSub = fitDiff
splineParGPS = fitGPS$optPar[(3+nKnots):(2+nKnots+nKnotsGPS)]
XiGPS = getSplineBasis(csz, c(minLat, maxLat), nKnotsGPS)
fitGPS$tvec = fitGPS$tvec - XiGPS %*% splineParGPS


plotSplineFit = function(fit, nKnotsGPS=0, diffGPSTaper=FALSE) {
  covMat = solve(-fit$hess)
  splineParI = 3:(2+nKnots+nKnotsGPS)
  splineCovMat = covMat[splineParI, splineParI]
  ggplotSplineUncertainty(fit$optPar[splineParI], splineCovMat, nKnots, 
                        diffGPSTaper=diffGPSTaper, nKnotsGPS=5, latsOnX=FALSE, main="", 
                        uncertaintyBands=TRUE)
}

pdf(file="taperComparison.pdf", width=8, height=10)
p1= plotSplineFit(fitComb, main)
p2 = ggPlotFaultDat(csz, fitComb$tvec, c(0,1), main="")

p3=plotSplineFit(fitSub)
p4=ggPlotFaultDat(csz, fitSub$tvec, c(0,1), main="")

p5=plotSplineFit(fitGPS)
p6=ggPlotFaultDat(csz, fitGPS$tvec, c(0,1), main="")

multiplot(p1, p2, p3, p4, p5, p6, byrow=TRUE, cols=2)
dev.off()


###################################
###################################
###################################
###################################
# compare model slip means, sub preds, mags

fitList = list(fitComb, fitSub, fitGPS)
ggCompareModels(fitList, nsim=5000, G=G, latRange=c(minLat, maxLat))
ggCompareModels(fitList, nsim=5000, G=G, latRange=c(minLat, maxLat), 
                fileNameRoot="posN", posNormalModel=rep(TRUE, 3))

###################################
###################################
###################################
###################################
# give summaries for slip and simulations for for slip and subs.  All for combined and sub models

# first genereate slip and subsidence simulations:
slipPredsComb = preds(fitComb$MLEs, nsim=10000, fault=csz, tvec=fitComb$tvec, 
                      posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitComb$MLEs[length(fitComb$MLEs)])
subPredsComb = predsToSubsidence(fitComb$MLEs, slipPreds, G=G, useMVNApprox=FALSE, subDat=dr1, 
                              posNormalModel=FALSE, normalModel=TRUE, tvec=fitComb$tvec)

slipPredsSub = preds(fitSub$MLEs, nsim=10000, fault=csz, tvec=fitSub$tvec, 
                      posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitSub$MLEs[length(fitSub$MLEs)])
subPredsSub = predsToSubsidence(fitSub$MLEs, slipPreds, G=G, useMVNApprox=FALSE, subDat=dr1, 
                                 posNormalModel=FALSE, normalModel=TRUE, tvec=fitSub$tvec)

# now compute summary statistics for slip for combinaed and subsidence models
lowComb = apply(slipPredsComb$slipSims, 1, quantile, probs=.025)
hiComb = apply(slipPredsComb$slipSims, 1, quantile, probs=.975)
meanComb = rowMeans(slipPredsComb$slipSims)

lowSub = apply(slipPredsSub$slipSims, 1, quantile, probs=.025)
hiSub = apply(slipPredsSub$slipSims, 1, quantile, probs=.975)
meanSub = rowMeans(slipPredsSub$slipSims)

slipMat = cbind(lowComb, meanComb, hiComb, lowSub, meanSub, hiSub)
ggplotSlipGrid(slipMat, nc=3, fileNameRoot="Summary")

# now plot some simulations against the data:

ggplotSlipGrid(slipPredsComb$slipSims[,1:9], nc=3, fileNameRoot="comb")
ggplotSubsidenceGrid(allSubsComb$subSims[,1:9], nc=3, fileNameRoot="comb")

ggplotSlipGrid(slipPredsSub$slipSims[,1:9], nc=3, fileNameRoot="sub")
ggplotSubsidenceGrid(allSubsSub$subSims[,1:9], nc=3, fileNameRoot="sub")

###################################
###################################
###################################
###################################
# Moving on to T1 predictions

# subset data to only be T1 data (and the Okada matrix)
isT1 = events=="T1"
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

## combined model (normal):
tvec = fitComb$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitComb$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
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
ggComparePredsToSubs(fitComb$MLEs, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("combT1"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE)

## combined model (posnormal):
tvec = fitComb$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitComb$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
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
ggComparePredsToSubs(fitComb$MLEs, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("combT1posN"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE)

## subsidence model (normal):
tvec = fitSub$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitSub$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
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
ggComparePredsToSubs(fitSub$MLEs, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("subT1"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE)

## subsidence model (posnormal):
tvec = fitSub$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitSub$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
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
ggComparePredsToSubs(fitSub$MLEs, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("subT1"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE)

### now plot simulations from combined model only:

tvec = fitComb$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitComb$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
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
subPreds = GT1 %*% slipSims

ggplotSlipGrid(slipPreds[,1:9], nc=3, fileNameRoot="combT1")
ggplotSubsidenceGrid(subPreds[,1:9], nc=3, fileNameRoot="combT1")



