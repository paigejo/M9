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
splineFit26k = doFitSpline(initParams=initParams, dStar=dStar, useMVNApprox = TRUE)
splinePar = splineFit26k$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
save(splineFit26k, file="splineFit26k.RData")

dStar = 40000
nKnots=4
splineInit = getInitialSplineEsts(MLEs[2], MLEs[3], MLEs[1], G, dStar=dStar)
tvec = getTaperSpline(splineInit$betaHat, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
initParams = c(MLEs[2:3], splineInit$betaHat)
splineFit = doFitSpline(initParams=initParams, dStar=dStar, useMVNApprox = TRUE)
splinePar = splineFit$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
save(splineFit, file="splineFit40k.RData")

# plot taper
latRange = c(40, 50)
lats = seq(40, 50, l=100)
Xi = bs(lats, df=nKnots, intercept=FALSE, Boundary.knots=latRange)
Xi = cbind(rep(1, nrow(Xi)), Xi)
plot(lats, Xi %*% splinePar, type="l", col="blue")
varMat = solve(-splineFit$hess)
plotSplineUncertainty = function(splinePar, covMat, nKnots=4, niter=1000)

# test out the ASL approximation method:
dStar = 26000
nKnots=4
splineInit = getInitialSplineEsts(MLEs[2], MLEs[3], MLEs[1], G, dStar=dStar)
initParams = c(MLEs[2:3], splineInit$betaHat)
initParams[1] = 2.5
initParams[2]= .226
# splineFitASL = doFitSpline(initParams=initParams, dStar=dStar, useASLApprox=TRUE, nsim=20000)
splineFitASL = doFitSpline(initParams=initParams, dStar=dStar, useMVNApprox=TRUE)
splinePar = splineFitASL$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
save(splineFitASL, file="splineFitASL.RData")

# now impute historic quakes:
load("splineFit.RData")
tvec = c(splineFit$tvec)
params=splineFit$MLEs
historicQuakes = updateMu(params, G=G, tvec=tvec, niter=3000)
save(historicQuakes, file="historicQuakes.RData")

# preload zeta correlation matrix
load("arealCSZCor.RData")
corMatCSZ = arealCSZCor

# get initial parameters at which to start optimization
dStar = 26000
nKnots=4
splineInit = getSplineEstsMomentMatch(MLEs[2], MLEs[3], MLEs[1], corMatCSZ, G=G, dStar=dStar)
splineParInit = splineInit$betaHat
initParams = c(MLEs[2:3], splineParInit)
tvec = getTaperSpline(splineParInit, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
comparePar = c(MLEs[1:3], 0.25, MLEs[5], splineParInit)
comparePredsToSubs(comparePar, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "splineInit26k")       

# begin optimization
splineFit26k = doFitSplineCoord(initParams=initParams, dStar=dStar, useMVNApprox=TRUE, G=G)
splinePar = splineFit26k$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
save(splineFit26k, file="splineFit26k.RData")
params = splineFit26k$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

# plot optimization results
comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "splineCoord26k")       
plotSplineUncertainty(splinePar, splineFit26k$hess, niter=10000)

# get initial parameters at which to start optimization using gradient for mu and sigma
dStar = 26000
nKnots=5
splineInit = getSplineEstsMomentMatch(MLEs[2], MLEs[3], MLEs[1], corMatCSZ, G=G, dStar=dStar, 
                                      initPar = c(MLEs[2], MLEs[3], MLEs[1], rep(0, nKnots-1)))
splineParInit = splineInit$betaHat
initParams = c(splineInit$muZeta, splineInit$sigmaZeta, splineParInit)
tvec = getTaperSpline(splineParInit, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
comparePar = c(MLEs[1], splineInit$muZeta, splineInit$sigmaZeta, 0.25, MLEs[5], splineParInit)
comparePredsToSubs(comparePar, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "splineInit26k") 

# now do additional optimization without gradient
splineInit2 = getSplineEstsMomentMatch(MLEs[2], MLEs[3], MLEs[1], corMatCSZ, G=G, dStar=dStar, 
                                       initPar = initParams, useGrad=F)
splineParInit2 = splineInit2$betaHat
initParams2 = c(splineInit$muZeta, splineInit$sigmaZeta, splineParInit2)
tvec2 = getTaperSpline(splineParInit2, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec2, varRange=c(0, 1))
comparePar2 = c(MLEs[1], splineInit2$muZeta, splineInit2$sigmaZeta, 0.25, MLEs[5], splineParInit2)
comparePredsToSubs(comparePar2, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec2, nsim=1000, fileNameRoot = "splineInit26k2")

# OMG THIS WORKS?

# now do final likelihood fit
splineFit26k = doFitSpline(initParams=initParams2, dStar=dStar, useMVNApprox=TRUE)
splinePar = splineFit26k$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
save(splineFit26k, file="splineFit26k.RData")
params = splineFit26k$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

# try likelihood with gradient
nKnots=5
dStar=26000
# splineFit26k = doFitSpline(initParams=initParams2, dStar=dStar, useMVNApprox=TRUE, useGrad=TRUE, nKnots=nKnots)
splineFit26k = doFitSpline(initParams=c(MLEs[2], MLEs[3], MLEs[1], rep(0, nKnots-1)), dStar=dStar, useMVNApprox=TRUE, 
                           useGrad=TRUE, nKnots=nKnots)
splinePar = splineFit26k$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)
save(splineFit26k, file="splineFit26kGrad.RData")
params = splineFit26k$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "splineInit26kGrad")
solve(-splineFit26k$hess)

# try likelihood with gradient, but this time don't include GPS 
# data in likelihood
nKnots=5
dStar=26000
# splineFit26k = doFitSpline(initParams=initParams2, dStar=dStar, useMVNApprox=TRUE, useGrad=TRUE, nKnots=nKnots)
splineFit26kSub = doFitSpline(initParams=c(MLEs[2], MLEs[3], MLEs[1], rep(0, nKnots-1)), dStar=dStar, useMVNApprox=TRUE, 
                           useGrad=TRUE, nKnots=nKnots)
splinePar = splineFit26kSub$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)
save(splineFit26kSub, file="splineFit26kGradSub.RData")
params = splineFit26kSub$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "splineInit26kGradSub")
solve(-splineFit26kSub$hess)

test = doFitSpline(initParams=c(2, 1, c(.5, rep(0,4))), dStar=dStar, useMVNApprox=TRUE, 
                              useGrad=TRUE, nKnots=nKnots)
splinePar = test$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)
points(dr1$Lon, dr1$Lat, pch="+", col="red", cex=.8)
save(test, file="test")
params = test$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "test")

solve(-test$hess)

test = doFitSpline(initParams=opt5$par, dStar=dStar, useMVNApprox=TRUE, 
                   useGrad=TRUE, nKnots=nKnots)

##### try using 4 basis functions in spline

nKnots=4
dStar=26000
# splineFit26k = doFitSpline(initParams=initParams2, dStar=dStar, useMVNApprox=TRUE, useGrad=TRUE, nKnots=nKnots)
splineFit26k = doFitSpline(initParams=c(MLEs[2], MLEs[3], MLEs[1], rep(0, nKnots-1)), dStar=dStar, useMVNApprox=TRUE, 
                           useGrad=TRUE, nKnots=nKnots, maxit=1500)
splinePar = splineFit26k$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)
save(splineFit26k, file="splineFit26kGrad4knots.RData")
params = splineFit26k$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

comparePredsToSubs(params, G=G, plotNameRoot="", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "splineInit26kGrad")
solve(-splineFit26k$hess)

out = out5
muZeta = out$MLEs[2]
sigmaZeta = out$MLEs[3]
splinePar = out$MLEs[6:length(out$MLEs)]

tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)

solve(-out$hess)

out5 = doFitSpline(initParams=c(muZeta, sigmaZeta, splinePar), 
                   dStar=dStar, useMVNApprox=TRUE, useGrad=TRUE, nKnots=nKnots, maxit=300)

# > out2$MLEs
# [1]          NA   1.2601723   2.1999320   0.2500000   0.8373628  38.9823708 -47.9278810 -31.4851532 -40.9806075
# > out3$MLEs
# [1]         NA   1.429537   2.168051   0.250000   0.869244  42.064117 -53.442496 -32.883078 -44.986308
# > out4$MLEs
# [1]          NA   1.4473459   2.1668688   0.2500000   0.8704259  43.1912683 -55.4989727 -33.3519979 -46.4644689
# > out5$MLEs
# [1]          NA   1.4473484   2.1669040   0.2500000   0.8703907  43.1913053 -55.4989778 -33.3519969 -46.4645080

##### test while ensuring lambda > 0
nKnots=4
dStar=26000
# splineFit26k = doFitSpline(initParams=initParams2, dStar=dStar, useMVNApprox=TRUE, useGrad=TRUE, nKnots=nKnots)
splineFit26k = doFitSpline(initParams=c(MLEs[2], MLEs[3], MLEs[1], rep(0, nKnots-1)), dStar=dStar, useMVNApprox=TRUE, 
                           useGrad=TRUE, nKnots=nKnots, maxit=100)
splinePar = splineFit26k$splineParMLE
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec, varRange=c(0, 1))
map("world", "Canada", add=TRUE)
US(add=TRUE)
save(out, file="splineFit26kGrad4knots.RData")
params = splineFit26k$MLEs
muZeta = params[2]
sigmaZeta = params[3]
lambda0 = params[4]
muXi = params[5]

##### retry constant 1 fit
nKnots=1
dStar=26000
splineFitConst = doFitSpline(initParams=c(MLEs[2], MLEs[3], MLEs[1]), dStar=dStar, useMVNApprox=TRUE, 
                             useGrad=TRUE, nKnots=nKnots, maxit=100)
params = splineFitConst$MLEs

##### try with constrained spline basis
nKnots=25
dStar=21000
splineFit21k25 = doFitSpline(initParams=c(MLEs[2], MLEs[3], MLEs[6], rep(0, nKnots-1)), dStar=dStar, useMVNApprox=TRUE, 
                           useGrad=TRUE, nKnots=nKnots, maxit=500)
load("splineFit21k25.RData")
params = splineFit21k25$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec)
comparePredsToSubs(params, G=G, plotNameRoot="15 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "21k15")

splineFit21k5 = splineFit21k
save(splineFit21k25, file="splineFit21k25.RData")

# now generate predictions
load("splineFit21k25.RData")
params = splineFit21k25$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec)

# T1
isT1 = events=="T1"
T1Dat = dr1[isT1,]
GT1 = G[isT1,]

eventPreds21k25 = predsGivenSubsidence(params, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec)
save(eventPreds21k25, file="eventPreds21k25.RData")

# areal values of zeta
load("eventPreds21k15.RData")
eventPreds = eventPreds21k25
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
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1 25 knots", fileNameRoot = "T1_21k25")
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1 25 knots", logScale=TRUE, 
              fileNameRoot="T1_21k25Log")


##### try predictions with prior of 500 on max slip
nKnots=25
dStar=21000
priorMax=300

# now generate predictions
load("splineFit21k25.RData")
params = splineFit21k25$MLEs
splinePar = params[6:length(params)]
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar)
plotFault(csz, tvec)
sqrt((exp(params[3]^2) - 1)*exp(2*params[2] + params[3]^2))

# T1
isT1 = events=="T1"
T1Dat = dr1[isT1,]
GT1 = G[isT1,]

eventPreds21k25 = predsGivenSubsidence(params, subDat=T1Dat, niter=500, G=GT1, prior=FALSE, tvec=tvec, priorMaxSlip=priorMax)
save(eventPreds21k25, file="eventPreds21k25Max300.RData")

# areal values of zeta
load("eventPreds21k15.RData")
eventPreds = eventPreds21k25
muAreal = eventPreds$zetaEsts * tvec
sdAreal = eventPreds$zetaSD * tvec
medAreal = eventPreds$zetaMed * tvec
l95Areal = eventPreds$zeta025 * tvec
u95Areal = eventPreds$zeta975 * tvec
getMomentFromSlip(muAreal)

# get simulations
tab <- extract(eventPreds$predResults, permuted = TRUE)
zetaSims = t(tab$zeta)
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1 25 knots", fileNameRoot = "T1_21k25Max")
plotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1 25 knots", logScale=TRUE, 
              fileNameRoot="T1_21k25LogMax")
comparePredsToSubs(params, slipPreds=slipPreds, G=G, plotNameRoot="15 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, fileNameRoot = "21k25Max300")
comparePredsToSubs(params, slipPreds=slipPreds, G=G, plotNameRoot="15 knots ", savePlots=TRUE, tvec=tvec, nsim=1000, 
                   fileNameRoot = "21k25Max300Log", logScale=TRUE)
