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


