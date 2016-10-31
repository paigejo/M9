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
library(splines)

load("fixedFit_MVN.RData")
MLEs = fixedFitMVN$MLEs

# fullFit = doFullFit(c(1.557515, 1.165601, 0.7544919, 8.88326), nsim=50000)
# fixedFit = doFixedFit(c(1.557515, 1.165601, 0.7544919), nsim=10000)
# 
# fullFit = doFullFit(nsim=50000)
# fixedFit = doFixedFit(nsim=50000)
# 
# fullFitMVN = doFullFit(useMVNApprox=TRUE)
# fixedFitMVN = doFixedFit(useMVNApprox=TRUE)
# 
# # save(fullFit, file="fullFit_50000sims.RData")
# # save(fixedFit, file="fixedFit_50000sims.RData")
# # save(fullFitMVN, file="fullFit_MVN.RData")
# # save(fixedFitMVN, file="fixedFit_MVN.RData")
# # save(fullFitMVN, file="fullFitAreal_MVN.RData")
# # save(fixedFitMVN, file="fixedFitAreal_MVN.RData")
# 
# # load("fullFit_50000sims.RData")
# # load("fixedFit_50000sims.RData")
# # load("fullFit_MVN.RData")
# # load("fixedFit_MVN.RData")
# # load("fullFitAreal_MVN.RData")
# # load("fixedFitAreal_MVN.RData")
# 
# fullSlips = predsGivenGPS(fullFit$MLEs)
# fullSubs = predsToSubsidence(fullFit$MLEs, fullSlips)
# 
# fixedSlips = predsGivenGPS(fixedFit$MLEs)
# fixedSubs = predsToSubsidence(fixedFit$MLEs, fixedSlips)
# 
# fullSlipsMVN = predsGivenGPS(fullFitMVN$MLEs)
# fullSubsMVN = predsToSubsidence(fullFitMVN$MLEs, fullSlipsMVN)
# 
# fixedSlipsMVN = predsGivenGPS(fixedFitMVN$MLEs)
# fixedSubsMVN = predsToSubsidence(fixedFitMVN$MLEs, fixedSlipsMVN)
# 
# comparePredsToSubs(fullFit$MLEs, fullSlips, fullSubs, plotNameRoot="full")
# comparePredsToSubs(fixedFit$MLEs, fixedSlips, fixedSubs, plotNameRoot="fixed")
# comparePredsToSubs(fullFitMVN$MLEs, fullSlipsMVN, fullSubsMVN, plotNameRoot="fullMVN")
# comparePredsToSubs(fixedFitMVN$MLEs, fixedSlipsMVN, fixedSubsMVN, plotNameRoot="fixedMVN",savePlots = F)

print(system.time(testResults <- fitModelIterative(maxIter=3, niterMCMC=250))) #params fit maxIter times, mean fit maxIter-1 times
save(testResults, file="testResults.RData")

testResults <- fitModelIterative(maxIter=20, niterMCMC=500, saveFile="fullIterFit.RData") #params fit maxIter times, mean fit maxIter-1 times
testResults <- fitModelIterative(maxIter=20, niterMCMC=500, saveFile="fullIterFit2.RData", loadDat="fullIterFit2.RData")

testResults <- fitModelIterative(maxIter=30, niterMCMC=500, saveFile="fullIterFitPrior.RData", usePrior=TRUE) #params fit maxIter times, mean fit maxIter-1 times

testSpline = doFitSpline(useMVNApprox = TRUE) # typical intial params: 1, log(20), 1
testSpline2 = doFitSpline(initParams = c(MLEs[2], MLEs[3], rep(MLEs[1], 5)), useMVNApprox = TRUE)

testSpline3 = doFitSpline(nKnots = 4) # typical intial params: 1, log(20), 1
testSpline4 = doFitSpline(initParams = c(MLEs[2], MLEs[3], rep(MLEs[1], 4)), nKnots=4, useMVNApprox = TRUE)

test = testSpline3
nKnots=length(test$splineParMLE)
lats = seq(40, 50, l=200)
splineMat = bs(lats, df=nKnots, intercept=TRUE, Boundary.knots=latRange)
splinePar = test$splineParMLE

# first plot the lambdas curve
lambdas = splineMat %*% splinePar
plot(lats, lambdas, type="l")

# now plot the taper for the minimum and maximum taper values
minLambda = min(lambdas)
maxLambda = max(lambdas)
depths = seq(1, 21000, l=200)
plot(depths, taper(depths, lambda=minLambda), col="red", type="l")
lines(depths, taper(depths, lambda=maxLambda), col="green")

# now plot the taper values across the fault
lambdasCSZ = getTaperSpline(splinePar, nKnots=nKnots)
# lambdasCSZ = seq(1, 5, l=240)
tvec = taper(csz$depth, lambda=lambdasCSZ)
plotFault(csz, tvec, varRange=c(0, 1))

dStar = 26000
nKnots=4
splineInit = getInitialSplineEsts(MLEs[2], MLEs[3], MLEs[1], G, dStar=dStar)
initParams = c(MLEs[2:3], splineInit$betaHat)
splineFit = doFitSpline(initParams=initParams, dStar=dStar, useMVNApprox = TRUE)
splinePar = splineFit$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
save(splineFit, file="splineFit.RData")

# now impute historic quakes:
historicQuakes = updateMu(params, G=G, tvec=tvec, niter=1000)
save(historicQuakes, file="historicQuakes.RData")



