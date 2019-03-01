# this script is for cross-validation

library(fields)
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
library(foreach)
library(doParallel)
library(tmvtnorm)

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
load("finalFitCombRevised.RData")
load("finalFitDiffRevised.RData")

fitGPS = fitSub = fitDiff
splinePar = fitGPS$optPar[3:(2+nKnots)]
splineParGPS = fitGPS$optPar[(3+nKnots):(2+nKnots+nKnotsGPS)]
Xi = getSplineBasis(csz, c(minLat, maxLat), nKnotsGPS)
fitGPS$tvec = c(taper(getFaultCenters(csz)[,3], Xi %*% (splinePar-splineParGPS), dStar=dStar, normalize=TRUE))
XiGPS = getSplineBasis(NULL, c(minLat, maxLat), nKnotsGPS, lats=threshSlipDat$lat)
fitSub$tvecGPS = c(taper(threshSlipDat$Depth, XiGPS %*% splinePar, dStar=dStar, normalize=TRUE))
fitComb$tvec = c(fitComb$tvec)
fitSub$tvec = c(fitSub$tvec)

load("adjustedMuCombRevised.RData")
load("adjustedMuSubRevised.RData")
load("adjustedMuGPSRevised.RData")

##### do Marginal dist'n Cross-Validation by site:
##### [(Marginal Normal, Marginal Pos. Normal, Marginal Pos. Normal Adjusted) x (Comb, Sub, GPS)] x (T1, T2, ..., AVG)

### let's only include full-margin earthquakes with >=20 observations
allEvents = c("T1", "T2", "T3", "T4", "T5", "T6", "T7")

# for combined taper model:
params = fitComb$MLEs
tvec = fitComb$tvec
MSEsComb = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesComb = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsComb = matrix(NA, nrow=1, ncol=length(allEvents))
weightComb = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, FALSE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsComb[ev] = out$MSE
  biasesComb[ev] = out$bias
  nObsComb[ev] = out$nObs
  weightComb[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsComb[length(allEvents)+1] = sum(MSEsComb[1:length(allEvents)]*weightComb)/sum(weightComb)
biasesComb[length(allEvents)+1] = sum(biasesComb[1:length(allEvents)]*weightComb)/sum(weightComb)
save(MSEsComb, biasesComb, nObsComb, weightComb, file="marginalCombCVRevised.RData")
load("marginalCombCVRevised.RData")

# for subsidence taper model:
params = fitSub$MLEs
tvec = fitSub$tvec
MSEsSub = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesSub = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsSub = matrix(NA, nrow=1, ncol=length(allEvents))
weightSub = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev)) # 244 sec for parallel, # 457 sec for sequential (~2x speedup)
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, FALSE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsSub[ev] = out$MSE
  biasesSub[ev] = out$bias
  nObsSub[ev] = out$nObs
  weightSub[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSub[length(allEvents)+1] = sum(MSEsSub[1:length(allEvents)]*weightSub)/sum(weightSub)
biasesSub[length(allEvents)+1] = sum(biasesSub[1:length(allEvents)]*weightSub)/sum(weightSub)
save(MSEsSub, biasesSub, nObsSub, weightSub, file="marginalSubCVRevised.RData")
load("marginalSubCVRevised.RData")

# for GPS/locking taper model:
params = fitGPS$MLEs
tvec = fitGPS$tvec
MSEsGPS = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesGPS = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsGPS = matrix(NA, nrow=1, ncol=length(allEvents))
weightGPS = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, FALSE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsGPS[ev] = out$MSE
  biasesGPS[ev] = out$bias
  nObsGPS[ev] = out$nObs
  weightGPS[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsGPS[length(allEvents)+1] = sum(MSEsGPS[1:length(allEvents)]*weightGPS)/sum(weightGPS)
biasesGPS[length(allEvents)+1] = sum(biasesGPS[1:length(allEvents)]*weightGPS)/sum(weightGPS)
save(MSEsGPS, biasesGPS, nObsGPS, weightGPS, file="marginalGPSCVRevised.RData")
load("marginalGPSCVRevised.RData")

### Now do the same thing but for positive normal models (unadjusted):
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
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, TRUE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsCombPN[ev] = out$MSE
  biasesCombPN[ev] = out$bias
  nObsCombPN[ev] = out$nObs
  weightCombPN[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsCombPN[length(allEvents)+1] = sum(MSEsCombPN[1:length(allEvents)]*weightCombPN)/sum(weightCombPN)
biasesCombPN[length(allEvents)+1] = sum(biasesCombPN[1:length(allEvents)]*weightCombPN)/sum(weightCombPN)
save(MSEsCombPN, biasesCombPN, nObsCombPN, weightCombPN, file="marginalCombPNCVRevised.RData")
load("marginalCombPNCVRevised.RData")

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
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, TRUE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsSubPN[ev] = out$MSE
  biasesSubPN[ev] = out$bias
  nObsSubPN[ev] = out$nObs
  weightSubPN[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSubPN[length(allEvents)+1] = sum(MSEsSubPN[1:length(allEvents)]*weightSubPN)/sum(weightSubPN)
biasesSubPN[length(allEvents)+1] = sum(biasesSubPN[1:length(allEvents)]*weightSubPN)/sum(weightSubPN)
save(MSEsSubPN, biasesSubPN, nObsSubPN, weightSubPN, file="marginalSubPNCVRevised.RData")
load("marginalSubPNCVRevised.RData")

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
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, TRUE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsGPSPN[ev] = out$MSE
  biasesGPSPN[ev] = out$bias
  nObsGPSPN[ev] = out$nObs
  weightGPSPN[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsGPSPN[length(allEvents)+1] = sum(MSEsGPSPN[1:length(allEvents)]*weightGPSPN)/sum(weightGPSPN)
biasesGPSPN[length(allEvents)+1] = sum(biasesGPSPN[1:length(allEvents)]*weightGPSPN)/sum(weightGPSPN)
save(MSEsGPSPN, biasesGPSPN, nObsGPSPN, weightGPSPN, file="marginalGPSPNCVRevised.RData")
load("marginalGPSPNCVRevised.RData")

### Now do the same thing but for positive normal models (adjusted):
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
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, TRUE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsCombPNAdj[ev] = out$MSE
  biasesCombPNAdj[ev] = out$bias
  nObsCombPNAdj[ev] = out$nObs
  weightCombPNAdj[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsCombPNAdj[length(allEvents)+1] = sum(MSEsCombPNAdj[1:length(allEvents)]*weightCombPNAdj)/sum(weightCombPNAdj)
biasesCombPNAdj[length(allEvents)+1] = sum(biasesCombPNAdj[1:length(allEvents)]*weightCombPNAdj)/sum(weightCombPNAdj)
save(MSEsCombPNAdj, biasesCombPNAdj, nObsCombPNAdj, weightCombPNAdj, file="marginalCombPNAdjCVRevised.RData")
load("marginalCombPNAdjCVRevised.RData")

# for subsidence taper model:
params = fitSub$MLEs
params[2] = adjustedMuSub
tvec = fitSub$tvec
MSEsSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
weightSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, TRUE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsSubPNAdj[ev] = out$MSE
  biasesSubPNAdj[ev] = out$bias
  nObsSubPNAdj[ev] = out$nObs
  weightSubPNAdj[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSubPNAdj[length(allEvents)+1] = sum(MSEsSubPNAdj[1:length(allEvents)]*weightSubPNAdj)/sum(weightSubPNAdj)
biasesSubPNAdj[length(allEvents)+1] = sum(biasesSubPNAdj[1:length(allEvents)]*weightSubPNAdj)/sum(weightSubPNAdj)
save(MSEsSubPNAdj, biasesSubPNAdj, nObsSubPNAdj, weightSubPNAdj, file="marginalSubPNAdjCVRevised.RData")
load("marginalSubPNAdjCVRevised.RData")

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
  out = getEventMSE(params, csz, inflateDr1, thisEvent, 20000, G, NULL, TRUE, TRUE, 123,
                    TRUE, tvec, dStar, TRUE, threshSlipDat, anisotropic=TRUE)
  MSEsGPSPNAdj[ev] = out$MSE
  biasesGPSPNAdj[ev] = out$bias
  nObsGPSPNAdj[ev] = out$nObs
  weightGPSPNAdj[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsGPSPNAdj[length(allEvents)+1] = sum(MSEsGPSPNAdj[1:length(allEvents)]*weightGPSPNAdj)/sum(weightGPSPNAdj)
biasesGPSPNAdj[length(allEvents)+1] = sum(biasesGPSPNAdj[1:length(allEvents)]*weightGPSPNAdj)/sum(weightGPSPNAdj)
save(MSEsGPSPNAdj, biasesGPSPNAdj, nObsGPSPNAdj, weightGPSPNAdj, file="marginalGPSPNAdjCVRevised.RData")
load("marginalGPSPNAdjCVRevised.RData")

##### do Predictive dist'n Cross-Validation by site:
##### [(Predictive Normal, Predictive Pos. Normal, Predictive Pos. Normal Adjusted) x (Comb, Sub, GPS)] x (T1, T2, ..., AVG)

### let's only include full-margin earthquakes with >=20 observations
allEvents = c("T1", "T2", "T3", "T4", "T5", "T6", "T7")

# for combined taper model:
params = fitComb$MLEs
tvec = fitComb$tvec
MSEsComb = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesComb = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsComb = matrix(NA, nrow=1, ncol=length(allEvents))
weightComb = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar,
                tvec, NULL, TRUE, FALSE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE, 
                anisotropic=TRUE)
  MSEsComb[ev] = out$MSE
  biasesComb[ev] = out$bias
  nObsComb[ev] = out$nObs
  weightComb[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsComb[length(allEvents)+1] = sum(MSEsComb[1:length(allEvents)]*weightComb)/sum(weightComb)
biasesComb[length(allEvents)+1] = sum(biasesComb[1:length(allEvents)]*weightComb)/sum(weightComb)
save(MSEsComb, biasesComb, nObsComb, weightComb, file="predictiveCombCVRevised.RData")
load("predictiveCombCVRevised.RData")

# for subsidence taper model:
params = fitSub$MLEs
tvec = fitSub$tvec
MSEsSub = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesSub = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsSub = matrix(NA, nrow=1, ncol=length(allEvents))
weightSub = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev)) # 244 sec for parallel, # 457 sec for sequential (~2x speedup)
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar,
                tvec, NULL, TRUE, FALSE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE, 
                anisotropic=TRUE)
  MSEsSub[ev] = out$MSE
  biasesSub[ev] = out$bias
  nObsSub[ev] = out$nObs
  weightSub[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSub[length(allEvents)+1] = sum(MSEsSub[1:length(allEvents)]*weightSub)/sum(weightSub)
biasesSub[length(allEvents)+1] = sum(biasesSub[1:length(allEvents)]*weightSub)/sum(weightSub)
save(MSEsSub, biasesSub, nObsSub, weightSub, file="predictiveSubCVRevised.RData")
load("predictiveSubCVRevised.RData")

# for GPS/locking taper model:
params = fitGPS$MLEs
tvec = fitGPS$tvec
MSEsGPS = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesGPS = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsGPS = matrix(NA, nrow=1, ncol=length(allEvents))
weightGPS = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar,
                tvec, NULL, TRUE, FALSE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE, 
                anisotropic=TRUE)
  MSEsGPS[ev] = out$MSE
  biasesGPS[ev] = out$bias
  nObsGPS[ev] = out$nObs
  weightGPS[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsGPS[length(allEvents)+1] = sum(MSEsGPS[1:length(allEvents)]*weightGPS)/sum(weightGPS)
biasesGPS[length(allEvents)+1] = sum(biasesGPS[1:length(allEvents)]*weightGPS)/sum(weightGPS)
save(MSEsGPS, biasesGPS, nObsGPS, weightGPS, file="predictiveGPSCVRevised.RData")
load("predictiveGPSCVRevised.RData")

### Now do the same thing but for positive normal models (unadjusted):
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
  print(system.time(out <- doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar,
                    tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=FALSE, fastPNSim=TRUE, 
                    anisotropic=TRUE)))
  MSEsCombPN[ev] = out$MSE
  biasesCombPN[ev] = out$bias
  nObsCombPN[ev] = out$nObs
  weightCombPN[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsCombPN[length(allEvents)+1] = sum(MSEsCombPN[1:length(allEvents)]*weightCombPN)/sum(weightCombPN)
biasesCombPN[length(allEvents)+1] = sum(biasesCombPN[1:length(allEvents)]*weightCombPN)/sum(weightCombPN)
save(MSEsCombPN, biasesCombPN, nObsCombPN, weightCombPN, file="predictiveCombPNCVRevised.RData")
load("predictiveCombPNCVRevised.RData")

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
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE, 
                anisotropic=TRUE)
  MSEsSubPN[ev] = out$MSE
  biasesSubPN[ev] = out$bias
  nObsSubPN[ev] = out$nObs
  weightSubPN[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSubPN[length(allEvents)+1] = sum(MSEsSubPN[1:length(allEvents)]*weightSubPN)/sum(weightSubPN)
biasesSubPN[length(allEvents)+1] = sum(biasesSubPN[1:length(allEvents)]*weightSubPN)/sum(weightSubPN)
save(MSEsSubPN, biasesSubPN, nObsSubPN, weightSubPN, file="predictiveSubPNCVRevised.RData")
load("predictiveSubPNCVRevised.RData")

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
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE, 
                anisotropic=TRUE)
  MSEsGPSPN[ev] = out$MSE
  biasesGPSPN[ev] = out$bias
  nObsGPSPN[ev] = out$nObs
  weightGPSPN[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsGPSPN[length(allEvents)+1] = sum(MSEsGPSPN[1:length(allEvents)]*weightGPSPN)/sum(weightGPSPN)
biasesGPSPN[length(allEvents)+1] = sum(biasesGPSPN[1:length(allEvents)]*weightGPSPN)/sum(weightGPSPN)
save(MSEsGPSPN, biasesGPSPN, nObsGPSPN, weightGPSPN, file="predictiveGPSPNCVRevised.RData")
load("predictiveGPSPNCVRevised.RData")

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
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE, 
                anisotropic=TRUE)
  MSEsCombPNAdj[ev] = out$MSE
  biasesCombPNAdj[ev] = out$bias
  nObsCombPNAdj[ev] = out$nObs
  weightCombPNAdj[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsCombPNAdj[length(allEvents)+1] = sum(MSEsCombPNAdj[1:length(allEvents)]*weightCombPNAdj)/sum(weightCombPNAdj)
biasesCombPNAdj[length(allEvents)+1] = sum(biasesCombPNAdj[1:length(allEvents)]*weightCombPNAdj)/sum(weightCombPNAdj)
save(MSEsCombPNAdj, biasesCombPNAdj, nObsCombPNAdj, weightCombPNAdj, file="predictiveCombPNAdjCVRevised.RData")
load("predictiveCombPNAdjCVRevised.RData")

# for subsidence taper model:
params = fitSub$MLEs
params[2] = adjustedMuSub
tvec = fitSub$tvec
MSEsSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
biasesSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
nObsSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
weightSubPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
for(ev in 1:length(allEvents)) {
  thisEvent = allEvents[ev]
  print(paste0("Performing CV for event ", ev))
  out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar,
                tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE, 
                anisotropic=TRUE)
  MSEsSubPNAdj[ev] = out$MSE
  biasesSubPNAdj[ev] = out$bias
  nObsSubPNAdj[ev] = out$nObs
  weightSubPNAdj[ev] = out$weight

  print(paste0("Event MSE: ", out$MSE))
  print(paste0("Event bias: ", out$bias))
}
MSEsSubPNAdj[length(allEvents)+1] = sum(MSEsSubPNAdj[1:length(allEvents)]*weightSubPNAdj)/sum(weightSubPNAdj)
biasesSubPNAdj[length(allEvents)+1] = sum(biasesSubPNAdj[1:length(allEvents)]*weightSubPNAdj)/sum(weightSubPNAdj)
save(MSEsSubPNAdj, biasesSubPNAdj, nObsSubPNAdj, weightSubPNAdj, file="predictiveSubPNAdjCVRevised.RData")
load("predictiveSubPNAdjCVRevised.RData")
# 
# # for GPS/locking taper model:
# params = fitGPS$MLEs
# params[2] = adjustedMuGPS
# tvec = fitGPS$tvec
# MSEsGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
# biasesGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents)+1)
# nObsGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
# weightGPSPNAdj = matrix(NA, nrow=1, ncol=length(allEvents))
# for(ev in 1:length(allEvents)) {
#   thisEvent = allEvents[ev]
#   print(paste0("Performing CV for event ", ev))
#   out = doCVSub(params, params[2], csz, inflateDr1, thisEvent, 10000, G, FALSE, TRUE, dStar, 
#                 tvec, NULL, TRUE, TRUE, NULL, 123, TRUE, TRUE, threshSlipDat, inPar=TRUE)
#   MSEsGPSPNAdj[ev] = out$MSE
#   biasesGPSPNAdj[ev] = out$bias
#   nObsGPSPNAdj[ev] = out$nObs
#   weightGPSPNAdj[ev] = out$weight
#   
#   print(paste0("Event MSE: ", out$MSE))
#   print(paste0("Event bias: ", out$bias))
# }
# MSEsGPSPNAdj[length(allEvents)+1] = sum(MSEsGPSPNAdj[1:length(allEvents)]*weightGPSPNAdj)/sum(weightGPSPNAdj)
# biasesGPSPNAdj[length(allEvents)+1] = sum(biasesGPSPNAdj[1:length(allEvents)]*weightGPSPNAdj)/sum(weightGPSPNAdj)
# save(MSEsGPSPNAdj, biasesGPSPNAdj, nObsGPSPNAdj, weightGPSPNAdj, file="predictiveGPSPNAdjCVRevised.RData")
# load("predictiveGPSPNAdjCVRevised.RData")