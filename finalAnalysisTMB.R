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

fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doTaperDiffPenalty = TRUE, 
                      G=G, debugPlotting=TRUE, logPenaltyPar=log(1), logDiffPenaltyPar=log(1), sharedSpatialProcess=TRUE)
# load("fullFit.RData")
out = load(paste0("~/git/M9/fit_n5_dS25000_diffTRUE_nGPS5_GamTRUE_nGam7_dStarGPS40000_sdTRUE_", 
                  "nVar5_fixInflFALSE_fixedPenTRUE_logPen0_sharePenFALSE_MeanTRUE_nMu5_nMuGPS7_", 
                  "diffPenTRUE_logDiffPen0_fixDiffTRUE_diff0-4.605_pen0-4.605_diffMuFALSE_constrTRUE_", 
                  "useHyperFALSE.RData"))
fullFit = modelInfo

# make parameter table
tab = cbind(fullFit$report$value, fullFit$report$sd)
filterIn = c(1, 3:nrow(tab))
tab = tab[filterIn,]
colnames(tab) = c("Estimate", "SE")

# all the reported parameters are already exponentiated
# parNames = rownames(tab)
# exponentiateI = grepl("log", parNames)

# colnames(tabComb) = c("muZ", "sigmaZ", "gamma", paste0("beta_t", 1:nKnots), paste0("beta_t'", 1:nKnotsGPS), 
#                       paste0("beta_s", 1:nKnotsVar), "phi", "alpha", "psil", "psih", 
#                       )
# rownames(tabComb) = c("MLEs", "SEs")
# keepI = fullFit$report$value != 0
keepI = (nrow(tab)-3):nrow(tab)
library(xtable)
tab = tab[keepI,]
low = tab[,1] + qnorm(0.025, 0, tab[,2])
hi = tab[,1] + qnorm(0.975, 0, tab[,2])
tab = cbind(tab, low, hi)
colnames(tab) = c("Estimate", "SE", "95 CI", "Units")
xtable(tab, digits=3)

# get the Okada matrix and subsidence data for this specific event of interest
isT1 = events=="T1"
T1Dat = dr1[isT1,]
GT1 = G[isT1, ]

# inflate uncertainties as necessary
lowInflate = exp(fullFit$loglowInflateEst)
highInflate = exp(fullFit$loghighInflateEst)
highQual = as.numeric(T1Dat$quality) == 1
lowQual = as.numeric(T1Dat$quality) != 1
T1Dat$Uncertainty[lowQual] = T1Dat$Uncertainty[lowQual]*lowInflate
T1Dat$Uncertainty[highQual] = T1Dat$Uncertainty[highQual]*highInflate
inflateDr1 = dr1
highQual = as.numeric(inflateDr1$quality) == 1
lowQual = as.numeric(inflateDr1$quality) != 1
inflateDr1$Uncertainty[lowQual] = inflateDr1$Uncertainty[lowQual]*lowInflate
inflateDr1$Uncertainty[highQual] = inflateDr1$Uncertainty[highQual]*highInflate

##### plot T1 predictions
normalPreds = predsGivenSubsidenceTMB(fullFit, G=G, gpsDat = threshSlipDat, niter = 5000)

# areal values of zeta
tvec = fullFit$taperVecY
muAreal = normalPreds$zetaEsts * tvec
sdAreal = normalPreds$zetaSD * tvec
medAreal = normalPreds$zetaMed * tvec
l95Areal = normalPreds$zeta025 * tvec
u95Areal = normalPreds$zeta975 * tvec

# get simulations
tab <- normalPreds$predResults
zetaSims = tab$zeta
slipSims = sweep(zetaSims, 1, tvec, "*")
slipPreds = list(meanSlip=muAreal, slipSims=slipSims, Sigma=diag(normalPreds$zetaSD * normalPreds$zetaSD), 
                 muc=normalPreds$zetaEsts)
# list(meanSlip=meanSlip, slipSims=slipSims, Sigma=Sigma, Sigmac=Sigma, muc=muZetaCSZ, 
     # SigmacGPS = SigmaD, mucGPS=muZetaGPS)
# plot results:
ggplotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="T1", plotNameRoot="1700 event ", logScale=FALSE, 
                fileNameRoot=paste0("finalT1TMB"))
ggComparePredsToSubsTMB(fullFit, slipPreds=slipPreds, G=GT1, plotNameRoot="1700 event ", gpsDat=threshSlipDat, 
                     subDat=T1Dat, logScale=FALSE, fileNameRoot=paste0("finalT1TMB"), 
                     fault=csz, nsim = 50000, magLimits = c(9.0, 9.2), binwidth=.0025)

##### plot marginal predictions
normalPreds = predsTMB(fullFit, gpsDat = threshSlipDat, nsim = 50000)

# areal values of zeta
tvec = fullFit$taperVecY
muAreal = normalPreds$meanSlip
sdAreal = sqrt(diag(normalPreds$Sigma)) * tvec
medAreal = normalPreds$meanSlip
l95Areal = normalPreds$muc * tvec + qnorm(0.5, mean=0, sd=sdAreal)
u95Areal = normalPreds$muc * tvec + qnorm(0.95, mean=0, sd=sdAreal)

# get simulations
slipSims = normalPreds$slipSims
slipPreds = normalPreds

# list(meanSlip=meanSlip, slipSims=slipSims, Sigma=Sigma, Sigmac=Sigma, muc=muZetaCSZ, 
# SigmacGPS = SigmaD, mucGPS=muZetaGPS)
# plot results:
ggplotFixedSlip(muAreal, NULL, l95Areal, u95Areal, sdAreal, event="All", plotNameRoot="Marginal ", logScale=FALSE, 
                fileNameRoot=paste0("finalMarginalTMB"))

ggComparePredsToSubsTMB(fullFit, G=G, plotNameRoot="Marginal ", gpsDat=threshSlipDat, 
                        subDat=inflateDr1, logScale=FALSE, fileNameRoot=paste0("finalMarginalTMB"), 
                        fault=csz, nsim = 50000, magLimits = c(8.4, 9.6))

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
# give summaries for slip and simulations for for slip and subs.  All for combined and sub models

# first genereate slip and subsidence simulations:
slipPreds = preds(fullFit, nsim=10000, fault=csz, tvec=fitComb$tvec, 
                      posNormalModel=FALSE, normalModel=TRUE, phiZeta=fitComb$MLEs[length(fitComb$MLEs)-1], 
                      anisotropic=TRUE, taperedGPSDat=TRUE)
subPredsComb = predsToSubsidence(fitComb$MLEs, slipPredsComb, G=G, useMVNApprox=FALSE, subDat=inflateDr1, 
                                 posNormalModel=FALSE, normalModel=TRUE, tvec=fitComb$tvec, dStar=dStar)




