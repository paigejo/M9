# this script is for cross-validation

library(fields)
library(rstan)
library(corpcor) # for fast computation of pseudoinverse
setwd("~/git/M9/")
source("taper.R")
source("okada.R")
source('predictions.R')
source('plotSubfault.R')
source('loadTestData.r')
source('fitModel.R')
source("loadFloodDat.R")
source("test.R")
source("exploratoryAnalysisFuns.R") # -418.9, 319
source("splines.R")
source("priors.R")
source("myMLESpatialProcess.R")
source("fitModel2.R")
library(splines)
library(abind)
library(numDeriv)
library(VGAM)
library(ggmap)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(latex2exp)
library(maps)
library(mapdata)
library(gstat)
library(sp)
library(maptools)
library(gridExtra)
library(foreach)
library(doParallel)

# precompute G
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)
fauxG = getFauxG()

# set up other variables
depthThresh=21000
nKnots = 5
nKnotsGPS = 5
dStar=25000
set.seed(123)
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
threshSlipDat$slipErr = threshSlipDat$slipErr*3
minLat = min(c(csz$latitude, threshSlipDat$lat)) - .001
maxLat = max(c(csz$latitude, threshSlipDat$lat)) + .001
highQual = as.numeric(dr1$quality) == 1
lowQual = as.numeric(dr1$quality) != 1
lowInflate=1.75
lowInflateComb=1.75
lowInflateDiff=1.75
highInflate=1.25
highInflateComb=1.25
highInflateDiff=1.25
inflateDr1=dr1
inflateDr1$Uncertainty[lowQual] = inflateDr1$Uncertainty[lowQual]*lowInflateComb
inflateDr1$Uncertainty[highQual] = inflateDr1$Uncertainty[highQual]*highInflateComb

# load model fits and computed parameters
load("finalFitComb.RData")
load("finalFitDiff.RData")

fitGPS = fitSub = fitDiff
splinePar = fitGPS$optPar[3:(2+nKnots)]
splineParGPS = fitGPS$optPar[(3+nKnots):(2+nKnots+nKnotsGPS)]
Xi = getSplineBasis(csz, c(minLat, maxLat), nKnotsGPS)
fitGPS$tvec = c(taper(getFaultCenters(csz)[,3], Xi %*% (splinePar-splineParGPS), dStar=dStar, normalize=TRUE))
XiGPS = getSplineBasis(NULL, c(minLat, maxLat), nKnotsGPS, lats=threshSlipDat$lat)
fitSub$tvecGPS = c(taper(threshSlipDat$Depth, XiGPS %*% splinePar, dStar=dStar, normalize=TRUE))
fitComb$tvec = c(fitComb$tvec)
fitSub$tvec = c(fitSub$tvec)

load("adjustedMuComb.RData")
load("adjustedMuSub.RData")
load("adjustedMuGPS.RData")

### perform CV for positive normal models (unadjusted):
# let's only include full-margin earthquakes with >=20 observations
allEvents = c("T1", "T2", "T3", "T4", "T5", "T6", "T7")

# for combined taper model:
params = fitComb$MLEs
tvec = fitComb$tvec
MSEsCombPN = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesCombPN = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsCombPN = matrix(NA, nrow=1, ncol=length(allEvents))
weightCombPN = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar, 
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat)
  MSEsCombPN[ev] = out$MSE
  biasesCombPN[ev] = out$bias
  nObsCombPN[ev] = out$nObs
  weightCombPN[ev] = out$weight
  
  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsCombPN[length(allEvents)+1] = sum(MSEsCombPN*weightCombPN)/sum(weightCombPN)
biasesCombPN[length(allEvents)+1] = sum(biasesCombPN*weightCombPN)/sum(weightCombPN)
save(MSEsCombPN, biasesCombPN, nObsCombPN, weightCombPN, file="marginalCombPNCV.RData")
load("marginalCombPNCV.RData")

# for subsidence taper model:
params = fitSub$MLEs
tvec = fitSub$tvec
MSEsSubPN = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesSubPN = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsSubPN = matrix(NA, nrow=1, ncol=length(allEvents))
weightSubPN = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar, 
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat)
  MSEsSubPN[ev] = out$MSE
  biasesSubPN[ev] = out$bias
  nObsSubPN[ev] = out$nObs
  weightSubPN[ev] = out$weight
  
  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSubPN[length(allEvents)+1] = sum(MSEsSubPN*weightSubPN)/sum(weightSubPN)
biasesSubPN[length(allEvents)+1] = sum(biasesSubPN*weightSubPN)/sum(weightSubPN)
save(MSEsSubPN, biasesSubPN, nObsSubPN, weightSubPN, file="marginalSubPNCV.RData")
load("marginalSubPNCV.RData")

# for GPS/locking taper model:
params = fitGPS$MLEs
tvec = fitGPS$tvec
MSEsGPSPN = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesGPSPN = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsGPSPN = matrix(NA, nrow=1, ncol=length(allEvents))
weightGPSPN = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar, 
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat)
  MSEsGPSPN[ev] = out$MSE
  biasesGPSPN[ev] = out$bias
  nObsGPSPN[ev] = out$nObs
  weightGPSPN[ev] = out$weight
  
  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsGPSPN[length(allEvents)+1] = sum(MSEsGPSPN*weightGPSPN)/sum(weightGPSPN)
biasesGPSPN[length(allEvents)+1] = sum(biasesGPSPN*weightGPSPN)/sum(weightGPSPN)
save(MSEsGPSPN, biasesGPSPN, nObsGPSPN, weightGPSPN, file="marginalGPSPNCV.RData")
load("marginalGPSPNCV.RData")

### Now do the same thing but for positive normal models (unadjusted):
# for combined taper model:
params = fitComb$MLEs
tvec = fitComb$tvec
params[2] = adjustedMuComb
MSEsCombPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesCombPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsCombPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
weightCombPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar, 
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat)
  MSEsCombPNAdj[ev] = out$MSE
  biasesCombPNAdj[ev] = out$bias
  nObsCombPNAdj[ev] = out$nObs
  weightCombPNAdj[ev] = out$weight
  
  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsCombPNAdj[length(allEvents)+1] = sum(MSEsCombPNAdj*weightCombPNAdj)/sum(weightCombPNAdj)
biasesCombPNAdj[length(allEvents)+1] = sum(biasesCombPNAdj*weightCombPNAdj)/sum(weightCombPNAdj)
save(MSEsCombPNAdj, biasesCombPNAdj, nObsCombPNAdj, weightCombPNAdj, file="marginalCombPNAdjCV.RData")
load("marginalCombPNAdjCV.RData")

# for subsidence taper model:
params = fitGPS$MLEs
params[2] = adjustedMuSub
tvec = fitGPS$tvec
MSEsSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
weightSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar, 
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat)
  MSEsSubPNAdj[ev] = out$MSE
  biasesSubPNAdj[ev] = out$bias
  nObsSubPNAdj[ev] = out$nObs
  weightSubPNAdj[ev] = out$weight
  
  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSubPNAdj[length(allEvents)+1] = sum(MSEsSubPNAdj*weightSubPNAdj)/sum(weightSubPNAdj)
biasesSubPNAdj[length(allEvents)+1] = sum(biasesSubPNAdj*weightSubPNAdj)/sum(weightSubPNAdj)
save(MSEsSubPNAdj, biasesSubPNAdj, nObsSubPNAdj, weightSubPNAdj, file="marginalSubPNAdjCV.RData")
load("marginalSubPNAdjCV.RData")

# for GPS/locking taper model:
params = fitGPS$MLEs
params[2] = adjustedMuGPS
tvec = fitGPS$tvec
MSEsGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
weightGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar, 
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat)
  MSEsGPSPNAdj[ev] = out$MSE
  biasesGPSPNAdj[ev] = out$bias
  nObsGPSPNAdj[ev] = out$nObs
  weightGPSPNAdj[ev] = out$weight
  
  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsGPSPNAdj[length(allEvents)+1] = sum(MSEsGPSPNAdj*weightGPSPNAdj)/sum(weightGPSPNAdj)
biasesGPSPNAdj[length(allEvents)+1] = sum(biasesGPSPNAdj*weightGPSPNAdj)/sum(weightGPSPNAdj)
save(MSEsGPSPNAdj, biasesGPSPNAdj, nObsGPSPNAdj, weightGPSPNAdj, file="marginalGPSPNAdjCV.RData")
load("marginalGPSPNAdjCV.RData")
