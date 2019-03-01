# final fitting process:

## do the initial setup
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

## set up down-dip slip limit and gps dataset: set threshold, match error with 
## empirical estimates from Pollitz and Evans (2017)
depthThresh=21000
nKnots = 5
nKnotsGPS = 5
dStar=25000
dStarGPS=40000
set.seed(123)
threshSlipDat = slipDatCSZ[slipDatCSZ$Depth<depthThresh,]
threshSlipDat$slipErr = threshSlipDat$slipErr*3
minLat = min(c(csz$latitude, threshSlipDat$lat)) - .001
maxLat = max(c(csz$latitude, threshSlipDat$lat)) + .001
nKnotsGamma=7

## calculate inflations for high and low quality subsidence data (under combined model)
highQual = as.numeric(dr1$quality) == 1
lowQual = as.numeric(dr1$quality) != 1

inflates = rev(c(.5, .75, 1, 1.25, 1.5, 1.75, 2))
LLs = matrix(nrow=length(inflates), ncol=length(inflates))
LLsSub = matrix(nrow=length(inflates), ncol=length(inflates))
for(i in 1:length(inflates)) {
  highInflate = inflates[i]
  
  for(j in 1:length(inflates)) {
    lowInflate = inflates[j]
    
    print(paste0("lowInflate: ", lowInflate, "; highInflate: ", highInflate))
    print(paste0("highInflate: ", i, "/", length(inflates), "; lowInflate: ", j, "/", length(inflates)))
    print(paste0("run ", (i - 1) * length(inflates) + j, "/", length(inflates)^2))
    
    inflateDR1 = dr1
    inflateDR1$Uncertainty[highQual] = dr1$Uncertainty[highQual] * highInflate
    inflateDR1$Uncertainty[lowQual] = dr1$Uncertainty[lowQual] * lowInflate
    
    initPar=c(20,15, 1, rep(0, nKnots-1), 175, 1)
    fit = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                    useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                    fauxG=fauxG, subDat=inflateDR1, fault=csz, 
                    normalModel=TRUE, normalizeTaper=TRUE, doHess=FALSE, 
                    latRange=c(minLat, maxLat), corGPS=TRUE, anisotropic=TRUE, 
                    verbose=FALSE)
    endPar = fit$MLEs
    endPar = c(endPar[2], endPar[3], endPar[6:length(endPar)])
    
    # get full log likelihood and also subsidence log-likelihood
    LLs[i,j] = fit$logLikMLE
    LLsSub[i,j] = fit$optimTable[nrow(fit$optimTable), ncol(fit$optimTable)-3]
    
    # print results:
    print(paste0("Low inflate: ", lowInflate, ".  High inflate: ", highInflate, ".  LnLik: ", LLs[i,j], ".  LnLikSub: ", LLsSub[i,j]))
    print(endPar)
  }
}
save(inflates=inflates, LLs=LLs, LLsSub, 
     file="inflatesFinalCombRevised.RData")
# load("inflatesFinalComb.RData")

# now plot results, calculate MLEs for inflation rates.  Rows of LLs correspond to highInflates
# 1.75 for low, 1 for high quality inflation, no matter whether using the full or subsidence log-likelihood!
par(mfrow=c(1,1))
highInflatesMat = matrix(rep(rev(inflates), length(inflates)), nrow=length(inflates))
lowInflatesMat = t(highInflatesMat)
tmp = LLs
LLs2 = t(matrix(rev(tmp), nrow=nrow(tmp)))
# image(lowInflatesMat, highInflatesMat, LLs2, col=tim.colors())
pdf(file="inflationCombRevised.pdf", width=5, height=5)
image.plot(rev(inflates), rev(inflates), LLs2, col=tim.colors(), 
           main="", xlab="Low Quality Inflation", 
           ylab="High Quality Inflation")
max(tmp)
maxI = which.max(tmp)
highI = row(tmp)[maxI]
lowI = col(tmp)[maxI]
highInflateComb=inflates[highI]
lowInflateComb=inflates[lowI]
points(lowInflateComb, highInflateComb, col="green", pch="x")
dev.off()
# 1.25 for high quality, 1.75 for low quality

## calculate inflations for high and low quality subsidence data (under different tapers model)
LLs = matrix(nrow=length(inflates), ncol=length(inflates))
LLsSub = matrix(nrow=length(inflates), ncol=length(inflates))
for(i in 1:length(inflates)) {
  highInflate = inflates[i]
  
  for(j in 1:length(inflates)) {
    lowInflate = inflates[j]
    
    print(paste0("lowInflate: ", lowInflate, "; highInflate: ", highInflate))
    print(paste0("highInflate: ", i, "/", length(inflates), "; lowInflate: ", j, "/", length(inflates)))
    print(paste0("run ", (i - 1) * length(inflates) + j, "/", length(inflates)^2))
    
    inflateDR1 = dr1
    inflateDR1$Uncertainty[highQual] = dr1$Uncertainty[highQual] * highInflate
    inflateDR1$Uncertainty[lowQual] = dr1$Uncertainty[lowQual] * lowInflate
    
    initPar=c(20,15, 1, rep(0, nKnots-1 + nKnotsGPS), 175, 1)
    fit = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                    useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, 
                    fauxG=fauxG, subDat=inflateDR1, fault=csz, corGPS=TRUE, 
                    normalModel=TRUE, normalizeTaper=TRUE, latRange=c(minLat, maxLat), 
                    doHess=FALSE, diffGPSTaper=TRUE, nKnotsGPS=nKnotsGPS, anisotropic=TRUE, 
                    verbose=FALSE)
    endPar = fit$MLEs
    endPar = c(endPar[2], endPar[3], endPar[6:length(endPar)])
    
    # get full log likelihood and also subsidence log-likelihood
    LLs[i,j] = fit$logLikMLE
    LLsSub[i,j] = fit$optimTable[nrow(fit$optimTable), ncol(fit$optimTable)-3]
    
    # print results:
    print(paste0("Low inflate: ", lowInflate, ".  High inflate: ", highInflate, ".  LnLik: ", LLs[i,j], ".  LnLikSub: ", LLsSub[i,j]))
  }
}
save(inflates=inflates, LLs=LLs, LLsSub, 
     file="inflatesFinalDiffRevised.RData")

# now plot results, calculate MLEs for inflation rates.  Rows of LLs correspond to highInflates
# 1.75 for low, 1 for high quality inflation, no matter whether using the full or subsidence log-likelihood!
par(mfrow=c(1,1))
highInflatesMat = matrix(rep(rev(inflates), length(inflates)), nrow=length(inflates))
lowInflatesMat = t(highInflatesMat)
tmp = LLsSub
LLs2 = t(matrix(rev(tmp), nrow=nrow(tmp)))
# image(lowInflatesMat, highInflatesMat, LLs2, col=tim.colors())
pdf(file="inflationDiff.pdf", width=5, height=5)
image.plot(rev(inflates), rev(inflates), LLs2, col=tim.colors(), 
           main="", xlab="Low Quality Inflation", 
           ylab="High Quality Inflation")
max(tmp)
maxI = which.max(tmp)
highI = row(tmp)[maxI]
lowI = col(tmp)[maxI]
highInflateDiff=inflates[highI] # 1.25
lowInflateDiff=inflates[lowI] # 1.75
points(lowInflateDiff, highInflateDiff, col="green", pch="x")
dev.off()

#### Now we perform the main analysis:

### for combined model:
# inflate the subsidence uncertainty
highQual = as.numeric(dr1$quality) == 1
lowQual = as.numeric(dr1$quality) != 1
# lowInflate=1.75
# lowInflateComb=1.75
# lowInflateDiff=1.75
# highInflate=1.25
# highInflateComb=1.25
# highInflateDiff=1.25
inflateDr1=dr1
inflateDr1$Uncertainty[lowQual] = inflateDr1$Uncertainty[lowQual]*lowInflateComb
inflateDr1$Uncertainty[highQual] = inflateDr1$Uncertainty[highQual]*highInflateComb

## fit the model (twice, with second time on smaller parameter scale to ensure 
## near optimum for hessian calculations)
initPar=c(20,15, 1, rep(0, nKnots-1), 175, 1)
# initPar=c(16.81,13.63, -1.168, -0.238, 3.871, 2.477, 1.005, 153, 1)
fit = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat,
                useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, corGPS=TRUE, 
                fauxG=fauxG, subDat=inflateDr1, fault=csz, latRange=c(minLat, maxLat), 
                normalModel=TRUE, normalizeTaper=TRUE, doHess=FALSE, anisotropic=TRUE, 
                doGammaSpline=FALSE, nKnotsGamma=nKnotsGamma)
fit = fitModel2(initParams=fit$optPar, dStar=dStar, gpsDat=threshSlipDat, 
                useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, corGPS=TRUE, 
                fauxG=fauxG, subDat=inflateDr1, fault=csz, latRange=c(minLat, maxLat), 
                normalModel=TRUE, normalizeTaper=TRUE, doHess=TRUE, 
                finalFit=TRUE, anisotropic=TRUE, doGammaSpline=TRUE, nKnotsGamma=nKnotsGamma, 
                dStarGPS=dStarGPS)
#fitComb = fit
save(fitComb, file="finalFitCombRevised.RData")
load("finalFitCombRevised.RData")

# make parameter table
params = fitComb$MLEs
MLEs = c(params[c(2, 3, 5, 6:(5+nKnots))], rep(NA, nKnotsGPS), params[(length(params) - 1):length(params)], lowInflate, highInflate)
SEs = sqrt(diag(solve(-fitComb$hess)))
SEs = c(SEs[1:2], NA, SEs[3:(length(SEs)-2)], rep(NA, nKnotsGPS), SEs[(length(SEs) - 1):length(SEs)], NA, NA)
tabComb = rbind(MLEs, SEs)
colnames(tabComb) = c("muZ", "sigmaZ", "gamma", paste0("beta", 1:nKnots), paste0("beta'", 1:nKnotsGPS), "phi", "alpha", "psil", "psih")
rownames(tabComb) = c("MLEs", "SEs")
library(xtable)
xtable(tabComb, digits=3)

## make ggplot plots:
isT1 = events=="T1"
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

params = fitComb$MLEs
splinePar = params[6:(length(params)-2)]
Xi = getSplineBasis(csz, c(40,50), 5)
lambdas = Xi %*% splinePar
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, normalize=TRUE)
plotFault(csz, tvec)

normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=1000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat, anisotropic=TRUE)

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
                fileNameRoot=paste0("finalT1CombRevised"))
ggComparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("finalT1CombRevised"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=TRUE, anisotropic=TRUE)

### for different taper model:
# inflate the subsidence uncertainty
# highQual = as.numeric(dr1$quality) == 1
# lowQual = as.numeric(dr1$quality) != 1
# lowInflate=1.75
# highInflate=1
inflateDr1=dr1
inflateDr1$Uncertainty[lowQual] = inflateDr1$Uncertainty[lowQual]*lowInflateDiff
inflateDr1$Uncertainty[highQual] = inflateDr1$Uncertainty[highQual]*highInflateDiff

## fit the model (twice, with second time on smaller parameter scale to ensure 
## near optimum for hessian calculations)
initPar=c(20,15, 1, rep(0, nKnots + nKnotsGPS -1), 175, 1)
fit = fitModel2(initParams=initPar, dStar=dStar, gpsDat=threshSlipDat, 
                useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, corGPS=TRUE, 
                fauxG=fauxG, subDat=inflateDr1, fault=csz, latRange=c(minLat, maxLat), 
                normalModel=TRUE, normalizeTaper=TRUE, doHess=FALSE, 
                diffGPSTaper=TRUE, nKnotsGPS=nKnotsGPS, anisotropic=TRUE)
fit = fitModel2(initParams=fit$optPar, dStar=dStar, gpsDat=threshSlipDat, 
                useGrad=TRUE, nKnots=nKnots, maxit=500, G=G, corGPS=TRUE, 
                fauxG=fauxG, subDat=inflateDr1, fault=csz, latRange=c(minLat, maxLat), 
                normalModel=TRUE, normalizeTaper=TRUE, doHess=TRUE, 
                diffGPSTaper=TRUE, nKnotsGPS=nKnotsGPS, finalFit=TRUE, 
                anisotropic=TRUE, doGammaSpline=TRUE, nKnotsGamma=nKnotsGamma)
fitDiff = fit
save(fitDiff, file="finalFitDiffRevised.RData")
load("finalFitDiffRevised.RData")

## make parameter table
params = fitDiff$MLEs
MLEs = c(params[c(2, 3, 5, 6:length(params))], lowInflate, highInflate)
SEs = sqrt(diag(solve(-fitDiff$hess)))
SEs = c(SEs[1:2], NA, SEs[3:length(SEs)], NA, NA)
tabDiff = rbind(MLEs, SEs)
colnames(tabDiff) = c("muZ", "sigmaZ", "gamma", paste0("beta", 1:nKnots), paste0("beta'", 1:nKnotsGPS), "phi", "alpha", "psil", "psih")
rownames(tabDiff) = c("MLEs", "SEs")
library(xtable)
xtable(tabDiff, digits=3)

# combine both tables:
summaryTab = cbind(t(tabComb), t(tabDiff))
xtable(summaryTab, digits=3)

##### 
###################################
###################################
###################################
###################################
# compare tapers of different models

# first modify GPS model taper
fitGPS = fitSub = fitDiff
splinePar = fitGPS$optPar[3:(2+nKnots)]
splineParGPS = fitGPS$optPar[(3+nKnots):(2+nKnots+nKnotsGPS)]
Xi = getSplineBasis(csz, c(minLat, maxLat), nKnotsGPS)
fitGPS$tvec = c(taper(getFaultCenters(csz)[,3], Xi %*% (splinePar-splineParGPS), dStar=dStar, normalize=TRUE))
XiGPS = getSplineBasis(NULL, c(minLat, maxLat), nKnotsGPS, lats=threshSlipDat$lat)
fitSub$tvecGPS = c(taper(threshSlipDat$Depth, XiGPS %*% splinePar, dStar=dStar, normalize=TRUE))
fitComb$tvec = c(fitComb$tvec)
fitSub$tvec = c(fitSub$tvec)

plotSplineFit = function(fit, nKnotsGPS=0, useDiffGPSTaper=FALSE, uncertaintyBands=TRUE) {
  covMat = solve(-fit$hess)
  splineParI = 3:(2+nKnots+nKnotsGPS)
  splineCovMat = covMat[splineParI, splineParI]
  ggplotSplineUncertainty(fit$optPar[splineParI], splineCovMat, nKnots,
                          diffGPSTaper=useDiffGPSTaper, nKnotsGPS=nKnotsGPS, latsOnX=FALSE, main="",
                          uncertaintyBands=uncertaintyBands)
  # ggplotTaperDepthUncertainty(fit$optPar[splineParI], splineCovMat, nKnots,
  #                             diffGPSTaper=useDiffGPSTaper, nKnotsGPS=nKnotsGPS, latsOnX=FALSE,
  #                             uncertaintyBands=uncertaintyBands, depthFrac=.9, confLevel=.95, dStar=dStar)
}

pdf(file="taperComparisonRevised.pdf", width=8, height=10)
p1= plotSplineFit(fitComb, uncertaintyBands=FALSE) +
  coord_fixed(xlim = c(-1.5,1.75),  ylim = c(40,50), ratio = 1, expand=FALSE) +
  scale_x_continuous(TeX("$\\lambda$"), c(-1, 0, 1), labels=c("-1", "0", "1"), limits=c(-360, 360)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.margin=unit(c(0,-4,0,-1), "cm"))
p2 = ggPlotFaultDat(csz, fitComb$tvec, c(0,1), main="", ylim=c(40,50), clab="Taper", lwd=.5) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='lightblue1'), 
        plot.margin=unit(c(0,-1,0,-3), "cm"))

p3=plotSplineFit(fitSub, uncertaintyBands=FALSE) +
  coord_fixed(xlim = c(-1,6.5),  ylim = c(40,50), ratio = 2.5, expand=FALSE) +
  scale_x_continuous(TeX("$\\lambda$"), seq(0,6, by=2), labels=c("0", "2", "4", "6"), limits=c(-360, 360)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.margin=unit(c(0,-4,0,-1), "cm"))
p4=ggPlotFaultDat(csz, fitSub$tvec, c(0,1), main="",  ylim = c(40,50), clab="Taper", lwd=.5) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='lightblue1'), 
        plot.margin=unit(c(0,-1,0,-3), "cm"))

p5=plotSplineFit(fitGPS, nKnotsGPS, useDiffGPSTaper=TRUE, uncertaintyBands=FALSE) +
  coord_fixed(xlim = c(-1.5,1.5),  ylim = c(40,50), ratio = 1, expand=FALSE) +
  scale_x_continuous(TeX("$\\lambda$"), c(-1, 0, 1), labels=c("-1", "0", "1"), limits=c(-360, 360)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill='white'),
        plot.margin=unit(c(0,-4,0,-1), "cm"))
p6=ggPlotFaultDat(csz, fitGPS$tvec, c(0,1), main="",  ylim = c(40,50), clab="Taper", lwd=.5) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='lightblue1'), 
        plot.margin=unit(c(0,-1,0,-3), "cm"))

multiplot(p1, p2, p3, p4, p5, p6, byrow=TRUE, cols=2)
dev.off()


###################################
###################################
###################################
###################################
# plot normalized residuals

pdf(file="residualsRevised.pdf", width=8, height=10)
pl1 = ggplotLockingResiduals(fitComb, fitComb$tvecGPS, threshSlipDat, c(minLat, maxLat), doGammaSpline=TRUE) + 
  scale_x_continuous("Normalized Residuals", c(-1, 0, 1), labels=c("-1", "0", "1"), limits=c(-360, 360)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white'), 
        plot.margin=unit(c(.5,-6,.5,-1), "cm"))
pl2 = ggplotSubsidenceResiduals(fitComb, fitComb$tvec, inflateDr1, G, c(minLat, maxLat)) + 
  scale_x_continuous("Normalized Residuals", (-1):4, labels=c("-1", "0", "1", "2", "3", "4"), limits=c(-360, 360)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white'), 
        plot.margin=unit(c(.5,-1,.5,-1), "cm"))

pl3 = ggplotLockingResiduals(fitSub, fitSub$tvecGPS, threshSlipDat, c(minLat, maxLat), doGammaSpline=TRUE) + 
  scale_x_continuous("Normalized Residuals", c(-1, 0, 1), labels=c("-1", "0", "1"), limits=c(-360, 360)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white'), 
        plot.margin=unit(c(.5,-6,.5,-1), "cm"))
pl4 = ggplotSubsidenceResiduals(fitSub, fitSub$tvec, inflateDr1, G, c(minLat, maxLat)) + 
  scale_x_continuous("Normalized Residuals", (-1):4, labels=c("-1", "0", "1", "2", "3", "4"), limits=c(-360, 360)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white'), 
        plot.margin=unit(c(.5,-1,.5,-1), "cm"))

pl5 = ggplotLockingResiduals(fitGPS, fitGPS$tvecGPS, threshSlipDat, c(minLat, maxLat), doGammaSpline=TRUE) + 
  scale_x_continuous("Normalized Residuals", c(-1, 0, 1), labels=c("-1", "0", "1"), limits=c(-360, 360)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white'), 
        plot.margin=unit(c(.5,-6,.5,-1), "cm"))
pl6 = ggplotSubsidenceResiduals(fitGPS, fitGPS$tvec, inflateDr1, G, c(minLat, maxLat)) + 
  scale_x_continuous("Normalized Residuals", (-1):4, labels=c("-1", "0", "1", "2", "3", "4"), limits=c(-360, 360)) + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_rect(fill='white'), 
        plot.margin=unit(c(.5,-1,.5,-1), "cm"))
multiplot(pl1, pl2, pl3, pl4, pl5, pl6, byrow=TRUE, cols=2)
dev.off()

# save progress
save.image("finalAnalysis1.RData")
load("finalAnalysis1.RData")

###################################
###################################
###################################
###################################
# compare model slip means, sub preds, mags

fitList = list(fitComb, fitSub, fitGPS)
ggCompareModels(fitList, nsim=20000, G=G, latRange=c(minLat, maxLat), magRange=c(8.5, 9.5), lwd=.2, 
                magTicks=seq(8.6, 9.4, by=.2), plotSubPreds=TRUE, plotNameRoot="subPredsRevised", 
                subPredMeanRange=c(-3,2.3), varRange=c(0, 18), anisotropic=TRUE)
ggCompareModels(fitList, nsim=5000, G=G, latRange=c(minLat, maxLat), magRange=c(8.5, 9.5), lwd=.2, 
                magTicks=seq(8.6, 9.4, by=.2), plotSubPreds=TRUE, subPredMeanRange=c(-4,2.3),
                plotNameRoot="subPredsPNRevised", posNormalModel=rep(TRUE, 3), varRange=c(0, 25), 
                anisotropic=TRUE)

# the positive normal mean adjusted results are calculated on line 885

###################################
###################################
###################################
###################################
# give summaries for slip and simulations for for slip and subs.  All for combined and sub models

# first genereate slip and subsidence simulations:
slipPredsComb = preds(fitComb$MLEs, nsim=10000, fault=csz, tvec=fitComb$tvec, 
                      posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitComb$MLEs[length(fitComb$MLEs)-1], 
                      anisotropic=TRUE, taperedGPSDat=TRUE)
subPredsComb = predsToSubsidence(fitComb$MLEs, slipPredsComb, G=G, useMVNApprox=FALSE, subDat=inflateDr1, 
                                 posNormalModel=FALSE, normalModel=TRUE, tvec=fitComb$tvec, dStar=dStar)

slipPredsSub = preds(fitSub$MLEs, nsim=10000, fault=csz, tvec=fitSub$tvec, 
                     posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitSub$MLEs[length(fitSub$MLEs)-1], 
                     anisotropic=TRUE, taperedGPSDat=TRUE)
subPredsSub = predsToSubsidence(fitSub$MLEs, slipPredsSub, G=G, useMVNApprox=FALSE, subDat=inflateDr1, 
                                posNormalModel=FALSE, normalModel=TRUE, tvec=fitSub$tvec, 
                                dStar=dStar)

# now compute summary statistics for slip for combined and subsidence models
lowComb = apply(slipPredsComb$slipSims, 1, quantile, probs=.2)
hiComb = apply(slipPredsComb$slipSims, 1, quantile, probs=.8)
meanComb = rowMeans(slipPredsComb$slipSims)

lowSub = apply(slipPredsSub$slipSims, 1, quantile, probs=.2)
hiSub = apply(slipPredsSub$slipSims, 1, quantile, probs=.8)
meanSub = rowMeans(slipPredsSub$slipSims)

slipMat = cbind(lowComb, meanComb, hiComb, lowSub, meanSub, hiSub)
allPlots = ggplotSlipGrid(slipMat, nc=3, fileNameRoot="Summary", savePlots = FALSE, lwd=.2)
for(i in 1:length(allPlots)) {
  allPlots[[i]] = allPlots[[i]] + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='lightblue1'), 
          plot.margin=unit(c(0,0,0,0), "cm"))
}
pdf(file=paste0("Summary", "SlipGridRevised.pdf"), width=8, height=10)
multiplot(plotlist=allPlots, cols=3, byrow=TRUE)
dev.off()

# now plot some simulations against the data
# (`subsidences' are uplifts so must take negative):

ggplotSlipGrid(slipPredsComb$slipSims[,1:9], nc=3, fileNameRoot="combRevised")
ggplotSubsidenceGrid(-subPredsComb$subSims[,1:9], nc=3, fileNameRoot="combRevised")

ggplotSlipGrid(slipPredsSub$slipSims[,1:9], nc=3, fileNameRoot="subRevised")
ggplotSubsidenceGrid(-subPredsSub$subSims[,1:9], nc=3, fileNameRoot="subRevised")

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
                                   dStar=dStar, gpsDat=threshSlipDat, anisotropic=TRUE)

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
slipPredsComb = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
ggComparePredsToSubs(fitComb$MLEs, slipPreds=slipPredsComb, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("combT1Revised"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, magLimits=c(8.8, 9.2), 
                     anisotropic=TRUE)

## combined model (posnormal):
tvec = fitComb$tvec

# generate T1 predictive distributions (niter * 2 = 20000 simulations)
posnormalPreds = predsGivenSubsidence(fitComb$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat, anisotropic=TRUE)

# areal values of zeta
muAreal = posnormalPreds$zetaEsts * tvec
sdAreal = posnormalPreds$zetaSD * tvec
medAreal = posnormalPreds$zetaMed * tvec
l95Areal = posnormalPreds$zeta025 * tvec
u95Areal = posnormalPreds$zeta975 * tvec

# get simulations
tab <- posnormalPreds$predResults
zetaSims = tab$zeta
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
ggComparePredsToSubs(fitComb$MLEs, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("combT1posNRevised"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE, anisotropic=TRUE)

## subsidence model (normal):
tvec = fitSub$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitSub$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat, anisotropic=TRUE)

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
slipPredsSub = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
ggComparePredsToSubs(fitSub$MLEs, slipPreds=slipPredsSub, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("subT1Revised"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE, magLimits=c(8.8, 9.2), 
                     anisotropic=TRUE)

## subsidence model (posnormal):
tvec = fitSub$tvec

# generate T1 predictive distributions
posnormalPreds = predsGivenSubsidence(fitSub$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat, anisotropic=TRUE)

# areal values of zeta
muAreal = posnormalPreds$zetaEsts * tvec
sdAreal = posnormalPreds$zetaSD * tvec
medAreal = posnormalPreds$zetaMed * tvec
l95Areal = posnormalPreds$zeta025 * tvec
u95Areal = posnormalPreds$zeta975 * tvec

# get simulations
tab <- posnormalPreds$predResults
zetaSims = tab$zeta
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
ggComparePredsToSubs(fitSub$MLEs, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("subT1posNRevised"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE, anisotropic=TRUE)

## locking taper model (normal): (just as a test to see that it doesn't work)
tvec = fitGPS$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitGPS$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat, anisotropic=TRUE)

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
slipPredsGPS = list(meanSlip=muAreal, slipSims=slipSims)

# plot results:
ggComparePredsToSubs(fitSub$MLEs, slipPreds=slipPredsGPS, G=GT1, tvec=tvec, plotNameRoot="T1 ", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("gpsT1Revised"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE, noTitle=TRUE, magLimits=c(8.8, 9.2), 
                     anisotropic=TRUE)

## make a table summarizing model mean magnitudes, 95% confidence band in magnitude
slipSimsComb = slipPredsComb$slipSims
slipSimsSub = slipPredsSub$slipSims
slipSimsGPS = slipPredsGPS$slipSims
magsComb = apply(slipSimsComb, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsSub = apply(slipSimsSub, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsGPS = apply(slipSimsGPS, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
cleanMagsComb = magsComb[is.finite(magsComb)]
cleanMagsSub = magsSub[is.finite(magsSub)]
cleanMagsGPS = magsGPS[is.finite(magsGPS)]
tab = rbind(
  c(mean(cleanMagsComb), quantile(cleanMagsComb, probs=.025), quantile(cleanMagsComb, probs=.975)), 
  c(mean(cleanMagsSub), quantile(cleanMagsSub, probs=.025), quantile(cleanMagsSub, probs=.975)), 
  c(mean(cleanMagsGPS), quantile(cleanMagsGPS, probs=.025), quantile(cleanMagsGPS, probs=.975))
)
rownames(tab) = c("Combined", "Subsidence", "Locking")
colnames(tab) = c("Mean", "95\\% CI", "")
xtable(tab, digits=2)

### now plot simulations from combined model only:

tvec = fitComb$tvec

# generate T1 predictive distributions
normalPreds = predsGivenSubsidence(fitComb$MLEs, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat, anisotropic=TRUE)

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

subPreds = predsToSubsidence(fitComb$MLEs, slipPreds, csz, useMVNApprox = FALSE, G=GT1, 
                             subDat=T1Dat, normalModel=TRUE, tvec=tvec, normalizeTaper=TRUE, 
                             dStar=dStar)

ggplotSlipGrid(slipPreds$slipSims[,1:9], nc=3, fileNameRoot="combT1Revised", lwd=.2)
ggplotSubsidenceGrid(-subPreds$noiseSims[,1:9], nc=3, fileNameRoot="combT1Revised", subDat=T1Dat)


# save progress
save.image("finalAnalysis2.RData")
load("finalAnalysis2.RData")

##### compute percent chance all slips positive for each model
library(mvtnorm)
library(tmvtnorm)
params = fitComb$MLEs
# params = fitDiff$MLEs
muZeta = params[2]
sigmaZeta = params[3]
alpha = params[length(params)]
phiZeta = params[length(params) - 1]
nuZeta=3 / 2
# arealCSZCor = stationary.cov(cbind(csz$longitude, csz$latitude), Covariance="Matern", 
#                              theta=phiZeta, smoothness=3/2)

### Rather than training the fault, we redefine an axis to be the strike access in Euclidean space
### using a Lambert projection and PCA
out = straightenFaultLambert()
faultGeomStraight = out$fault
scale = out$scale
parameters = out$projPar
transformation = out$transformation

cszStraight = divideFault2(faultGeomStraight)
centers = getFaultCenters(csz)[,1:2]
newCenters = transformation(centers)
cszStraight$centerX = newCenters[,1]
cszStraight$centerY = newCenters[,2]

# calculate along strike and along dip squared distances in kilometers
strikeCoordsCSZ = cbind(0, cszStraight$centerY)
dipCoordsCSZ = cbind(cszStraight$centerX, 0)
squareStrikeDistCsz = rdist(strikeCoordsCSZ)^2
squareDipDistCsz = rdist(dipCoordsCSZ)^2

# compute gps, fault, and cross distance matrices
distMatCSZComb = sqrt(alpha^2 * squareStrikeDistCsz + alpha^(-2) * squareDipDistCsz)
coordsCSZ = cbind(cszStraight$longitude, cszStraight$latitude)
arealCSZCor = stationary.cov(coordsCSZ, Covariance="Matern", theta=phiZeta, 
                           smoothness=nuZeta, distMat = distMatCSZComb)
arealCSZCov = arealCSZCor * sigmaZeta^2
pmvnorm(upper=rep(0, nrow(csz)), mean=rep(-muZeta, nrow(csz)), sigma=arealCSZCov, abseps=.00001) # P(all pos) = (fitComb: 0.555679, fitDiff: 0.8202905)
1- pmvnorm(upper=rep(0, nrow(csz)), mean=rep(-muZeta, nrow(csz)), sigma=arealCSZCov, abseps=.00001) # P(any neg) = (fitComb: 0.4457753, fitDiff: 0.1816371)

##### do a quick chisq test for the extra taper parameters
1-pchisq(2*(fitDiff$logLikMLE - fitComb$logLikMLE), 5) # p-value is numerically 0...

##### Could try naively adjusting mean of positive normal models and looking at those predictions here

### for combined marginal model:
# get the parameters
params = fitComb$MLEs
muZeta = params[2]
sigmaZeta = params[3]
tvec = fitComb$tvec
nuZeta=3 / 2

# compute covariance matrix
# rather than using below code, use the same covariance as above
# arealCSZCor = stationary.cov(cbind(csz$longitude, csz$latitude), Covariance="Matern", 
#                              theta=params[length(params)], smoothness=3/2, 
#                              Distance="rdist.earth", Dist.args=list(miles=FALSE))
# arealCSZCov = arealCSZCor * sigmaZeta^2

# compute adjusted mean parameter
adjustedMuComb = getPosNormMu(muZeta, arealCSZCov, startN=150, initNewMu=8)
tempPar = params
tempPar[2] = adjustedMuComb
save(adjustedMuComb, file="adjustedMuCombRevised.RData")
load("adjustedMuCombRevised.RData")

# compute probability all positive
pmvnorm(upper=rep(0, nrow(csz)), mean=rep(-adjustedMuComb, nrow(csz)), sigma=arealCSZCov, abseps=.00001) # P(all pos) (0.1458336)
1- pmvnorm(upper=rep(0, nrow(csz)), mean=rep(-adjustedMuComb, nrow(csz)), sigma=arealCSZCov, abseps=.00001) # P(any neg) (0.8534005)

# get summary of predictions

ggComparePredsToSubs(tempPar, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("combPNAdjustedRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("combPNUnadjustedRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("combNRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=FALSE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

### for subsidence marginal model:
# get the parameters
params = fitSub$MLEs
muZeta = params[2]
sigmaZeta = params[3]
tvec = fitSub$tvec
phiZeta = params[length(params) - 1]
alpha = params[length(params)]
nuZeta=3 / 2

# compute covariance matrix
# arealCSZCor = stationary.cov(cbind(csz$longitude, csz$latitude), Covariance="Matern", 
#                              theta=params[length(params)], smoothness=3/2, 
#                              Distance="rdist.earth", Dist.args=list(miles=FALSE))
# arealCSZCov = arealCSZCor * sigmaZeta^2
distMatCSZDiff = sqrt(alpha^2 * squareStrikeDistCsz + alpha^(-2) * squareDipDistCsz)
coordsCSZ = cbind(cszStraight$longitude, cszStraight$latitude)
arealCSZCor = stationary.cov(coordsCSZ, Covariance="Matern", theta=phiZeta, 
                             smoothness=nuZeta, distMat = distMatCSZDiff)
arealCSZCov = arealCSZCor * sigmaZeta^2

# compute adjusted mean parameter
adjustedMuSub = getPosNormMu(muZeta, arealCSZCov, startN=nrow(csz), initNewMu=15)
tempPar = params
tempPar[2] = adjustedMuSub
adjustedMuGPS = adjustedMuSub
save(adjustedMuSub, file="adjustedMuSubRevised.RData")
save(adjustedMuGPS, file="adjustedMuGPSRevised.RData")
load("adjustedMuSubRevised.RData")
load("adjustedMuGPSRevised.RData")

# compute probability all positive
adjustedMu = adjustedMuGPS
pmvnorm(upper=rep(0, nrow(csz)), mean=rep(-adjustedMu, nrow(csz)), sigma=arealCSZCov, abseps=.00001) # P(all pos) (0.000941532)
1- pmvnorm(upper=rep(0, nrow(csz)), mean=rep(-adjustedMu, nrow(csz)), sigma=arealCSZCov, abseps=.00001) # P(any neg) (0.9990585)

# get summary of predictions

ggComparePredsToSubs(tempPar, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("subPNAdjustedRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("subPNUnadjustedRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("subNRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=FALSE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

### for GPS marginal model:
# get the parameters
params = fitGPS$MLEs
muZeta = params[2]
sigmaZeta = params[3]
tvec = fitGPS$tvec

# get adjusted mean parameter
load("adjustedMuGPS.RData")
tempPar = params
tempPar[2] = adjustedMuGPS


# get summary of predictions

ggComparePredsToSubs(tempPar, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("gpsPNAdjustedRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("gpsPNUnadjustedRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=TRUE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="", nsim=5000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("gpsNRevised"), 
                     fault=csz, normalModel=TRUE, posNormalModel=FALSE, useMVNApprox=FALSE, 
                     taperedGPSDat=TRUE, dStar=dStar, normalizeTaper=TRUE, noTitle=TRUE, 
                     anisotropic=TRUE)

# finish the magnitude results started on line 439
fitList = list(fitComb, fitSub, fitGPS)
adjustedMeans = c(adjustedMuComb, adjustedMuSub, adjustedMuGPS)
ggCompareModels(fitList, nsim=20000, G=G, latRange=c(minLat, maxLat), magRange=c(8.5, 9.5), lwd=.2, 
                magTicks=seq(8.6, 9.4, by=.2), plotSubPreds=TRUE, subPredMeanRange=c(-3,2.3),
                plotNameRoot="subPredsPNAdjRevised", posNormalModel=rep(TRUE, 3), adjustedMeans=adjustedMeans, 
                varRange=c(0, 18), anisotropic=TRUE)

# save progress
save.image("finalAnalysis3.RData")
load("finalAnalysis3.RData")

##### do Marginal dist'n Cross-Validation by site:
##### [(Marginal Normal, Marginal Pos. Normal, Marginal Pos. Normal Adjusted) x (Comb, Sub, GPS)] x (T1, T2, ..., AVG)

### first see how many observations and sites there are for each earthquake
nUnique = function(xs) {
  length(unique(xs))
}
nSites = tapply(dr1$Site, dr1$event, nUnique)
nObs = table(dr1$event)
eventNs = rbind(nSites, nObs)
theseEvents = c(21, 20, 19, 17, 15, 12, 10)
theseEventNs = eventNs[,theseEvents]
theseEventNs
xtable(theseEventNs)

# now do the CV
source("posNormCV.R")

# collect marginal CV results into table
load("marginalCombCVRevised.RData")
margBiasTab = biasesComb
margMSETab = MSEsComb
load("marginalSubCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesSub)
margMSETab = rbind(margMSETab, MSEsComb)
load("marginalGPSCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesGPS)
margMSETab = rbind(margMSETab, MSEsGPS)
load("marginalCombPNCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesCombPN)
margMSETab = rbind(margMSETab, MSEsCombPN)
load("marginalSubPNCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesSubPN)
margMSETab = rbind(margMSETab, MSEsSubPN)
load("marginalGPSPNCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesGPSPN)
margMSETab = rbind(margMSETab, MSEsGPSPN)
load("marginalCombPNAdjCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesCombPNAdj)
margMSETab = rbind(margMSETab, MSEsCombPNAdj)
load("marginalSubPNAdjCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesSubPNAdj)
margMSETab = rbind(margMSETab, MSEsSubPNAdj)
load("marginalGPSPNAdjCVRevised.RData")
margBiasTab = rbind(margBiasTab, biasesGPSPNAdj)
margMSETab = rbind(margMSETab, MSEsGPSPNAdj)
colnames(margBiasTab) = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "Avg")
rownames(margBiasTab) = c("Comb", "Sub", "GPS", "CombPN", "SubPN", "GPSPN", 
                          "CombPNAdj", "SubPNAdj", "GPSPNAdj")
round(margBiasTab, digits=3)
xtable(t(margBiasTab), digits=3)
colnames(margMSETab) = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "Avg")
rownames(margMSETab) = c("Comb", "Sub", "GPS", "CombPN", "SubPN", "GPSPN", 
                          "CombPNAdj", "SubPNAdj", "GPSPNAdj")
round(margMSETab, digits=3)
xtable(t(margMSETab), digits=3)

# collect predictive CV results into table
load("predictiveCombCV.RData")
predBiasTab = biasesComb
predMSETab = MSEsComb
load("predictiveSubCV.RData")
predBiasTab = rbind(predBiasTab, biasesSub)
predMSETab = rbind(predMSETab, MSEsComb)
load("predictiveGPSCV.RData")
predBiasTab = rbind(predBiasTab, biasesGPS)
predMSETab = rbind(predMSETab, MSEsGPS)
load("predictiveCombPNCV.RData")
predBiasTab = rbind(predBiasTab, biasesCombPN)
predMSETab = rbind(predMSETab, MSEsCombPN)
load("predictiveSubPNCV.RData")
predBiasTab = rbind(predBiasTab, biasesSubPN)
predMSETab = rbind(predMSETab, MSEsSubPN)
load("predictiveGPSPNCV.RData")
predBiasTab = rbind(predBiasTab, biasesGPSPN)
predMSETab = rbind(predMSETab, MSEsGPSPN)
load("predictiveCombPNAdjCV.RData")
predBiasTab = rbind(predBiasTab, biasesCombPNAdj)
predMSETab = rbind(predMSETab, MSEsCombPNAdj)
load("predictiveSubPNAdjCV.RData")
predBiasTab = rbind(predBiasTab, biasesSubPNAdj)
predMSETab = rbind(predMSETab, MSEsSubPNAdj)
load("predictiveGPSPNAdjCV.RData")
predBiasTab = rbind(predBiasTab, biasesGPSPNAdj)
predMSETab = rbind(predMSETab, MSEsGPSPNAdj)
colnames(predBiasTab) = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "Avg")
rownames(predBiasTab) = c("Comb", "Sub", "GPS", "CombPN", "SubPN", "GPSPN",
                          "CombPNAdj", "SubPNAdj", "GPSPNAdj")
round(predBiasTab, digits=3)
xtable(t(predBiasTab), digits=3)
colnames(predMSETab) = c("T1", "T2", "T3", "T4", "T5", "T6", "T7", "Avg")
rownames(predMSETab) = c("Comb", "Sub", "GPS", "CombPN", "SubPN", "GPSPN",
                         "CombPNAdj", "SubPNAdj", "GPSPNAdj")
round(predMSETab, digits=3)
xtable(t(predMSETab), digits=3)

# save progress
save.image("finalAnalysis4.RData")
load("finalAnalysis4.RData")

## marginal magnitude table
# combined taper model
params = fitComb$MLEs
slipPredsComb = preds(params, nsim=10000, fault=csz, tvec=fitComb$tvec, 
                      posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitComb$MLEs[length(fitComb$MLEs)])
slipPredsCombPN = preds(params, nsim=10000, fault=csz, tvec=fitComb$tvec, 
                      posNormalModel=TRUE, normalModel=TRUE, phiZeta=fitComb$MLEs[length(fitComb$MLEs)])
params[2] = adjustedMuComb
slipPredsCombPNAdj = preds(params, nsim=10000, fault=csz, tvec=fitComb$tvec, 
                           posNormalModel=TRUE, normalModel=TRUE, phiZeta=fitComb$MLEs[length(fitComb$MLEs)])
slipSimsComb = slipPredsComb$slipSims
slipSimsCombPN = slipPredsCombPN$slipSims
slipSimsCombPNAdj = slipPredsCombPNAdj$slipSims

magsComb = apply(slipSimsComb, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsCombPN = apply(slipSimsCombPN, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsCombPNAdj = apply(slipSimsCombPNAdj, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)

cleanMagsComb = magsComb[is.finite(magsComb)]
cleanMagsComb = cleanMagsComb[cleanMagsComb != 0]
cleanMagsCombPN = magsCombPN[is.finite(magsCombPN)]
cleanMagsCombPNAdj = magsCombPNAdj[is.finite(magsCombPNAdj)]

# subsidence taper model
params = fitSub$MLEs
slipPredsSub = preds(params, nsim=10000, fault=csz, tvec=fitSub$tvec, 
                      posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitSub$MLEs[length(fitSub$MLEs)])
slipPredsSubPN = preds(params, nsim=10000, fault=csz, tvec=fitSub$tvec, 
                        posNormalModel=TRUE, normalModel=TRUE, phiZeta=fitSub$MLEs[length(fitSub$MLEs)])
params[2] = adjustedMuSub
slipPredsSubPNAdj = preds(params, nsim=10000, fault=csz, tvec=fitSub$tvec, 
                           posNormalModel=TRUE, normalModel=TRUE, phiZeta=fitSub$MLEs[length(fitSub$MLEs)])
slipSimsSub = slipPredsSub$slipSims
slipSimsSubPN = slipPredsSubPN$slipSims
slipSimsSubPNAdj = slipPredsSubPNAdj$slipSims

magsSub = apply(slipSimsSub, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsSubPN = apply(slipSimsSubPN, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsSubPNAdj = apply(slipSimsSubPNAdj, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)

cleanMagsSub = magsSub[is.finite(magsSub)]
cleanMagsSubPN = magsSubPN[is.finite(magsSubPN)]
cleanMagsSubPNAdj = magsSubPNAdj[is.finite(magsSubPNAdj)]

# GPS taper model
params = fitGPS$MLEs
slipPredsGPS = preds(params, nsim=10000, fault=csz, tvec=fitGPS$tvec, 
                      posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitGPS$MLEs[length(fitGPS$MLEs)])
slipPredsGPSPN = preds(params, nsim=10000, fault=csz, tvec=fitGPS$tvec, 
                        posNormalModel=TRUE, normalModel=TRUE, phiZeta=fitGPS$MLEs[length(fitGPS$MLEs)])
params[2] = adjustedMuGPS
slipPredsGPSPNAdj = preds(params, nsim=10000, fault=csz, tvec=fitGPS$tvec, 
                           posNormalModel=TRUE, normalModel=TRUE, phiZeta=fitGPS$MLEs[length(fitGPS$MLEs)])
slipSimsGPS = slipPredsGPS$slipSims
slipSimsGPSPN = slipPredsGPSPN$slipSims
slipSimsGPSPNAdj = slipPredsGPSPNAdj$slipSims

magsGPS = apply(slipSimsGPS, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsGPSPN = apply(slipSimsGPSPN, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)
magsGPSPNAdj = apply(slipSimsGPSPNAdj, 2, getMomentFromSlip, fault=csz, dStar=dStar, normalizeTaper=TRUE)

cleanMagsGPS = magsGPS[is.finite(magsGPS)]
cleanMagsGPSPN = magsGPSPN[is.finite(magsGPSPN)]
cleanMagsGPSPNAdj = magsGPSPNAdj[is.finite(magsGPSPNAdj)]

# now make the table
tab = rbind(
  c(mean(cleanMagsComb), quantile(cleanMagsComb, probs=.025), quantile(cleanMagsComb, probs=.975)), 
  c(mean(cleanMagsCombPN), quantile(cleanMagsCombPN, probs=.025), quantile(cleanMagsCombPN, probs=.975)), 
  c(mean(cleanMagsCombPNAdj), quantile(cleanMagsCombPNAdj, probs=.025), quantile(cleanMagsCombPNAdj, probs=.975)), 
  c(mean(cleanMagsSub), quantile(cleanMagsSub, probs=.025), quantile(cleanMagsSub, probs=.975)), 
  c(mean(cleanMagsSubPN), quantile(cleanMagsSubPN, probs=.025), quantile(cleanMagsSubPN, probs=.975)), 
  c(mean(cleanMagsSubPNAdj), quantile(cleanMagsSubPNAdj, probs=.025), quantile(cleanMagsSubPNAdj, probs=.975)), 
  c(mean(cleanMagsGPS), quantile(cleanMagsGPS, probs=.025), quantile(cleanMagsGPS, probs=.975)), 
  c(mean(cleanMagsGPSPN), quantile(cleanMagsGPSPN, probs=.025), quantile(cleanMagsGPSPN, probs=.975)), 
  c(mean(cleanMagsGPSPNAdj), quantile(cleanMagsGPSPNAdj, probs=.025), quantile(cleanMagsGPSPNAdj, probs=.975))
)
colnames(tab) = c("Mean", "95 CI", "")
rownames(tab) = paste(c("Combined", "Combined", "Combined", "Subsidence", "Subsidence", "Subsidence", 
                         "Locking", "Locking", "Locking"), rep(c("Gaussian", "Pos. Gaussian", "Adj. Pos. Gaussian")))
round(tab, digits=2)
xtable(tab, digits=2)

## now do the same but for the predictive distribution
isT1 = events=="T1"
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

# combined model
params = fitComb$MLEs
tvec = fitComb$tvec
slipPredsComb = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                   normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                   dStar=dStar, gpsDat=threshSlipDat)
slipPredsCombPN = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                     normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                     dStar=dStar, gpsDat=threshSlipDat)
# params[2] = adjustedMuComb
# slipPredsCombPNAdj = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
#                                      normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
#                                      dStar=dStar, gpsDat=threshSlipDat)

cleanMagsComb = slipPredsComb$predResults$Mw
cleanMagsCombPN = slipPredsCombPN$predResults$Mw
# cleanMagsCombPNAdj = slipPredsCombPNAdj$predResults$Mw

# subsidence model
params = fitSub$MLEs
tvec = fitSub$tvec
slipPredsSub = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                     normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                     dStar=dStar, gpsDat=threshSlipDat)
slipPredsSubPN = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                       normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                       dStar=dStar, gpsDat=threshSlipDat)
# params[2] = adjustedMuSub
# slipPredsSubPNAdj = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
#                                      normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
#                                      dStar=dStar, gpsDat=threshSlipDat)

cleanMagsSub = slipPredsSub$predResults$Mw
cleanMagsSubPN = slipPredsSubPN$predResults$Mw
# cleanMagsSubPNAdj = slipPredsSubPNAdj$predResults$Mw

# combined model
params = fitGPS$MLEs
tvec = fitGPS$tvec
slipPredsGPS = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                     normalModel=TRUE, posNormalModel=FALSE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                     dStar=dStar, gpsDat=threshSlipDat)
slipPredsGPSPN = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
                                       normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
                                       dStar=dStar, gpsDat=threshSlipDat)
# params[2] = adjustedMuGPS
# slipPredsGPSPNAdj = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=5000, G=GT1, prior=FALSE, tvec=tvec, 
#                                      normalModel=TRUE, posNormalModel=TRUE, taperedGPSDat=TRUE, normalizeTaper=TRUE, 
#                                      dStar=dStar, gpsDat=threshSlipDat)

cleanMagsGPS = slipPredsGPS$predResults$Mw
cleanMagsGPSPN = slipPredsGPSPN$predResults$Mw
# cleanMagsGPSPNAdj = slipPredsGPSPNAdj$predResults$Mw

# now make the table
# tab = rbind(
#   c(mean(cleanMagsComb), quantile(cleanMagsComb, probs=.025), quantile(cleanMagsComb, probs=.975)), 
#   c(mean(cleanMagsCombPN), quantile(cleanMagsCombPN, probs=.025), quantile(cleanMagsCombPN, probs=.975)), 
#   c(mean(cleanMagsCombPNAdj), quantile(cleanMagsCombPNAdj, probs=.025), quantile(cleanMagsCombPNAdj, probs=.975)), 
#   c(mean(cleanMagsSub), quantile(cleanMagsSub, probs=.025), quantile(cleanMagsSub, probs=.975)), 
#   c(mean(cleanMagsSubPN), quantile(cleanMagsSubPN, probs=.025), quantile(cleanMagsSubPN, probs=.975)), 
#   c(mean(cleanMagsSubPNAdj), quantile(cleanMagsSubPNAdj, probs=.025), quantile(cleanMagsSubPNAdj, probs=.975)), 
#   c(mean(cleanMagsGPS), quantile(cleanMagsGPS, probs=.025), quantile(cleanMagsGPS, probs=.975)), 
#   c(mean(cleanMagsGPSPN), quantile(cleanMagsGPSPN, probs=.025), quantile(cleanMagsGPSPN, probs=.975)), 
#   c(mean(cleanMagsGPSPNAdj), quantile(cleanMagsGPSPNAdj, probs=.025), quantile(cleanMagsGPSPNAdj, probs=.975))
# )
# colnames(tab) = c("Mean", "95 CI", "")
# rownames(tab) = paste(c("Combined", "Combined", "Combined", "Subsidence", "Subsidence", "Subsidence", 
#                         "Locking", "Locking", "Locking"), rep(c("Gaussian", "Pos. Gaussian", "Adj. Pos. Gaussian")))
tab = rbind(
  c(mean(cleanMagsComb), quantile(cleanMagsComb, probs=.025), quantile(cleanMagsComb, probs=.975)), 
  c(mean(cleanMagsCombPN), quantile(cleanMagsCombPN, probs=.025), quantile(cleanMagsCombPN, probs=.975)), 
  c(mean(cleanMagsSub), quantile(cleanMagsSub, probs=.025), quantile(cleanMagsSub, probs=.975)), 
  c(mean(cleanMagsSubPN), quantile(cleanMagsSubPN, probs=.025), quantile(cleanMagsSubPN, probs=.975)), 
  c(mean(cleanMagsGPS), quantile(cleanMagsGPS, probs=.025), quantile(cleanMagsGPS, probs=.975)), 
  c(mean(cleanMagsGPSPN), quantile(cleanMagsGPSPN, probs=.025), quantile(cleanMagsGPSPN, probs=.975))
)
colnames(tab) = c("Mean", "95 CI", "")
rownames(tab) = paste(c("Combined", "Combined", "Subsidence", "Subsidence", 
                        "Locking", "Locking"), rep(c("Gaussian", "Pos. Gaussian")))
round(tab, digits=2)
xtable(tab, digits=2)

# subsidence taper

tab = rbind(
  c(mean(cleanMagsComb), quantile(cleanMagsComb, probs=.025), quantile(cleanMagsComb, probs=.975)), 
  c(mean(cleanMagsSub), quantile(cleanMagsSub, probs=.025), quantile(cleanMagsSub, probs=.975)), 
  c(mean(cleanMagsGPS), quantile(cleanMagsGPS, probs=.025), quantile(cleanMagsGPS, probs=.975))
)
rownames(tab) = c("Combined", "Subsidence", "Locking")
colnames(tab) = c("Mean", "95\\% CI", "")
xtable(tab, digits=2)

# save progress
save.image("finalAnalysis5.RData")
load("finalAnalysis5.RData")

##### Other plotting commands I once thought would be helpful:




## make ggplot plots (once for subsidence-based taper, once for GPS-based taper)

# for subsidence-based taper:
isT1 = events=="T1"
# T1DatRange = dr1[isT1 & inRangeDat,]
# GT1Range = G[isT1 & inRangeDat, inRangeFault]
# T1Dat = dr1[isT1,]
T1Dat = inflateDr1[isT1,]
GT1 = G[isT1, ]

params = fitDiff$MLEs
splinePar = params[6:(5 + nKnots)]
splineParGPS = params[(6 + nKnots):(length(params)-1)]
Xi = getSplineBasis(csz, c(minLat, maxLat), 5)
lambdas = Xi %*% splinePar
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, normalize=TRUE)
plotFault(csz, tvec)

# first plot marginal fit
ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="Marginal ", nsim = 10000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("finalDiffSub"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE)

# now generate predictions for 1700 event
normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
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
ggplotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=FALSE, 
                fileNameRoot=paste0("finalDiffSubT1"))
ggComparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("finalDiffSubT1"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=TRUE)

# for GPS-based taper:
lambdas = Xi %*% (splinePar - splineParGPS)
tvec = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, fault=csz, normalize=TRUE)
plotFault(csz, tvec)

# first plot marginal fit
ggComparePredsToSubs(params, G=G, tvec=tvec, plotNameRoot="Marginal ", nsim = 10000, 
                     subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("finalDiffGPS"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=FALSE)

# now generate predictions for 1700 event
normalPreds = predsGivenSubsidence(params, fault=csz, subDat=T1Dat, niter=10000, G=GT1, prior=FALSE, tvec=tvec, 
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
ggplotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="T1", logScale=FALSE, 
                fileNameRoot=paste0("finalDiffGPST1"))
ggComparePredsToSubs(params, slipPreds=slipPreds, G=GT1, tvec=tvec, plotNameRoot="T1", 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("finalDiffGPST1"), 
                     fault=csz, normalModel=TRUE, useMVNApprox=FALSE, taperedGPSDat=TRUE, 
                     dStar=dStar, normalizeTaper=TRUE)




