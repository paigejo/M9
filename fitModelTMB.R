library(TMB)
library(alabama)
# source("setup.R")

# fits the combined model with normalized taper using TMB jointly with the variance inflation 
# parameters as well as the gamma spline parameters
fitModelTMB = function(initParams=NULL, gpsDat=slipDatCSZ, gpsDepthThreshold=21000, 
                       G=NULL, subDat=dr1, fault=csz, nKnots=5, 
                       dStar=25000, maxit=500, latRange=c(40, 50), 
                       diffGPSTaper=TRUE, nKnotsGPS=nKnots, reltol=1e-8, 
                       nKnotsVar=5, includeGammaSpline=TRUE, nKnotsGamma=7, 
                       dStarGPS=40000, seed=123, debug=FALSE, doVarSpline=TRUE, 
                       fixInflation=FALSE, maxCount=1, fixedPenalty=TRUE, 
                       logPenaltyPar=log(10), sharedPenalty=FALSE, doMeanSpline=TRUE, 
                       nKnotsMean=5, nKnotsMeanGPS=7, doTaperDiffPenalty=FALSE, logDiffPenaltyPar=log(1), 
                       fixedDiffPenalty=FALSE, debugPlotting=TRUE, saveResults=TRUE, 
                       diffLow=log(1), diffHigh=log(100), penLow=log(1), penHigh=log(100), 
                       diffMean=FALSE, constrainMean=TRUE, muBarrier=1e-04, useAlabama=TRUE, 
                       useHyperpriors=FALSE, sharedSpatialProcess=FALSE, jointShared=FALSE) {
  
  # get all input values
  allInputs = as.list(environment())
  
  # threshold the gps data and multiplied standard deviations by three as in the literature
  threshSlipDat = gpsDat[gpsDat$Depth<gpsDepthThreshold,]
  threshSlipDat$slipErr = threshSlipDat$slipErr*3
  
  # get input parameters and separate them
  if(is.null(initParams))
    initParams = getInitialParameters(nKnots, nKnotsVar, nKnotsGamma, logScale=TRUE, diffGPSTaper=diffGPSTaper)
  
  out = getInputPar(initParams, fault, threshSlipDat, nKnots, diffGPSTaper=diffGPSTaper, nKnotsGPS, taperedGPSDat=TRUE, 
                    anisotropic=TRUE, normalModel=TRUE, nKnotsVar, doVarSpline=TRUE, 
                    includeGammaSpline=TRUE, nKnotsGamma=nKnotsGamma, includeInflation=TRUE)
  logmu = out$muZeta
  betaMean = rep(0, nKnotsMean - 1)
  logMeanGPS = 0
  betaMeanGPS = rep(0, nKnotsMeanGPS - 1)
  betaTaper = out$taperPar
  betaTaperIntercept = betaTaper[1]
  betaTaper = betaTaper[-1]
  betaTaperGPS = out$taperParGPS
  if(is.null(betaTaperGPS))
    betaTaperGPS = rep(0, nKnotsGPS)
  # taperParGPS = out$taperParGPS
  logphi = out$phiZeta
  # lambda0 = out$lambda0
  logalpha = out$alpha
  betasd = out$varPar
  betasdIntercept = betasd[1]
  betasd = betasd[-1]
  betaGamma = out$gammaPar
  betaGammaIntercept = betaGamma[1]
  betaGamma = betaGamma[-1]
  loglowInflate = out$lowInflation
  loghighInflate = out$highInflation
  parscale = out$parscale
  parNames = out$parNames
  parNamesTMB = out$parNamesTMB
  logitOmega = logit(0.25) # initially, assume most of the variation is taken up by iid spatial processes rather than shared
  
  # determine which of the input parameters will be optimized
  includeI = c(TRUE, TRUE, rep(doVarSpline, nKnotsVar - 1), rep(TRUE, nKnots), rep(TRUE, nKnotsGPS * diffGPSTaper), TRUE, 
               rep(includeGammaSpline, nKnotsGamma - 1), rep(!fixInflation, 2), TRUE, TRUE)
  parscale = parscale[includeI]
  parNames = parNames[includeI]
  parNamesTMBShort = parNamesTMB[includeI]
  
  # PARAMETER(logmu);
  # PARAMETER_VECTOR(betaTaperIntercept);
  # PARAMETER_VECTOR(betaTaper);
  # PARAMETER_VECTOR(betasdIntercept);
  # PARAMETER_VECTOR(betasd);
  # PARAMETER_VECTOR(betaGammaIntercept);
  # PARAMETER_VECTOR(betaGamma);
  # PARAMETER(logphi);
  # PARAMETER(logalpha);
  # PARAMETER(loglowInflate);
  # PARAMETER(loghighInflate);
  
  ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(fault)[,3]
  
  ### Rather than training the fault, we redefine an axis to be the strike access in Euclidean space
  ### using a Lambert projection and PCA
  out = straightenFaultLambert()
  faultGeomStraight = out$fault
  scale = out$scale
  parameters = out$projPar
  transformation = out$transformation
  
  ##### compute distance matrices for straightened fault and for straightened gps data
  cszStraight = divideFault2(faultGeomStraight)
  centers = getFaultCenters(csz)[,1:2]
  newCenters = transformation(centers)
  cszStraight$centerX = newCenters[,2]
  cszStraight$centerY = newCenters[,1]
  straightenedGpsCoords = transformation(cbind(threshSlipDat$lon, threshSlipDat$lat))
  
  # calculate along strike and along dip squared distances in kilometers
  strikeCoordsCSZ = cbind(0, cszStraight$centerY)
  dipCoordsCSZ = cbind(cszStraight$centerX, 0)
  squareStrikeDistCsz = rdist(strikeCoordsCSZ)^2
  squareDipDistCsz = rdist(dipCoordsCSZ)^2
  
  # do the same for the gps data
  strikeCoords = cbind(0, straightenedGpsCoords[,1])
  dipCoords = cbind(straightenedGpsCoords[,2], 0)
  squareStrikeDistGps = rdist(strikeCoords)^2
  squareDipDistGps = rdist(dipCoords)^2
  
  # calculate the cross distance matrix (GPS, CSZ)
  squareStrikeDistCross = rdist(strikeCoords, strikeCoordsCSZ)^2
  squareDipDistCross = rdist(dipCoords, dipCoordsCSZ)^2
  # squareStrikeDistJoint = rbind(cbind(squareStrikeDistGps, squareStrikeDistCross), 
  #                               cbind(t(squareStrikeDistCross), squareStrikeDistCsz))
  # squareDipDistJoint = rbind(cbind(squareDipDistGps, squareDipDistCross), 
  #                            cbind(t(squareDipDistCross), squareDipDistCsz))
  
  # compute all required inputs for the TMB cpp code
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)[,-1]
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)[,-1]
  lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
  meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)[,-1]
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)[,-1]
  meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)[,-1]
  DSStrikeCSZ = squareStrikeDistCsz
  DSDipCSZ = squareDipDistCsz
  DSStrikeGPS = squareStrikeDistGps
  DSDipGPS = squareDipDistGps
  zeroMask = eventsEqMask(subDat)
  lowI = as.numeric(as.numeric(subDat$quality) != 1)
  y = -subDat$subsidence
  ysd = subDat$Uncertainty
  x = threshSlipDat$slip
  xsd = threshSlipDat$slipErr
  xDepths = threshSlipDat$Depth
  faultDepths = cszDepths
  dStar = dStar
  dStarGPS = dStarGPS
  sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)[,-1]
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)[,-1]
  gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)[,-1]
  latSeq = seq(latRange[1], latRange[2], l=500)
  deltaPenalty = latSeq[2] - latSeq[1]
  meanBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)[,-1]
  meanBasisGPSPenalty = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)[,-1]
  sdBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnotsVar, lats=latSeq, latRange=latRange)[,-1]
  taperBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)[,-1]
  taperBasisGPSPenalty = getSplineBasis(fault=fault, nKnots=nKnotsGPS, lats=latSeq, latRange=latRange)
  gammaBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)[,-1]
  
  # should choose penalty parameter (kappa) so that 
  # kappa * deltaPenalty * sum((diff(latSeq / 10) / deltaPenalty)^2) = kappa * 0.1
  # is 90% chance between .1, and 10. implying that the penalty factors should be between 1 and 100
  penaltyMean = mean(c(penLow, penHigh))
  penaltySD = (penHigh-penaltyMean)/qnorm(.95)
  
  # should choose penalty parameter (kappa) so that 
  # kappa * integrate(function(x) {1}, 40, 50) = kappa * 10
  # is 90% chance between .1, and 1. implying that the penalty factors should be between 1/100 and 1/10
  diffPenaltyMean = mean(c(diffLow, diffHigh))
  diffPenaltySD = (diffHigh-diffPenaltyMean)/qnorm(.95)
  
  # compile the function and its gradient
  if(!debug)
    compile("fitModelTMB.cpp")
  else
    compile("fitModelTMB.cpp","-O0 -g")
  dyn.load(dynlib("fitModelTMB"))
  set.seed(seed)
  data <- list(y=y, x=x, lambdaBasisY=lambdaBasisY, lambdaBasisX=lambdaBasisX, lambdaBasisXGPS=lambdaBasisXGPS, DSStrikeCSZ=DSStrikeCSZ, 
               DSDipCSZ=DSDipCSZ, DSStrikeGPS=DSStrikeGPS, DSDipGPS=DSDipGPS, zeroMask=zeroMask, lowI=lowI, 
               ysd=ysd, xsd=xsd, xDepths=xDepths, faultDepths=faultDepths, dStar=dStar, sdBasisX=sdBasisX, 
               sdBasisY=sdBasisY, gammaBasis=gammaBasis, G=G, dStarGPS=dStarGPS, deltaPenalty=deltaPenalty, 
               sdBasisPenalty=sdBasisPenalty, taperBasisPenalty=taperBasisPenalty, taperBasisGPSPenalty=taperBasisGPSPenalty, 
               gammaBasisPenalty=gammaBasisPenalty, penaltyMean=penaltyMean, penaltySD=penaltySD, 
               sharedPenalty=as.numeric(sharedPenalty), meanBasisX=meanBasisX, meanBasisXGPS=meanBasisXGPS, meanBasisY=meanBasisY, 
               meanBasisPenalty=meanBasisPenalty, meanBasisGPSPenalty=meanBasisGPSPenalty, doTaperDiffPenalty=as.numeric(doTaperDiffPenalty), 
               diffPenaltyMean=diffPenaltyMean, diffPenaltySD=diffPenaltySD, diffMean=as.numeric(diffMean), 
               useHyperpriors=as.numeric(useHyperpriors), sharedSpatialProcess=as.numeric(sharedSpatialProcess), 
               jointShared=as.numeric(jointShared), DSStrikeCross=squareStrikeDistCross, 
               DSDipCross=squareDipDistCross)
  parameters = list(logmu=logmu, betasdIntercept=betasdIntercept, betasd=betasd, betaMean=betaMean, logMeanGPS=logMeanGPS, betaMeanGPS=betaMeanGPS, 
                    betaTaperIntercept=betaTaperIntercept, betaTaper=betaTaper, betaTaperGPS=betaTaperGPS, 
                    betaGammaIntercept=betaGammaIntercept, betaGamma=betaGamma, 
                    loglowInflate=loglowInflate, loghighInflate=loghighInflate, 
                    logphi=logphi, logalpha=logalpha, betasdPenaltyLogLambda=0, 
                    betaTaperPenaltyLogLambda=0, betaTaperGPSPenaltyLogLambda=0, 
                    betaGammaPenaltyLogLambda=0, betaMeanPenaltyLogLambda=0, 
                    taperDiffPenaltyLogLambda=0, betaMeanGPSPenaltyLogLambda=0, 
                    meanDiffPenaltyLogLambda=0, logitOmega=logitOmega)
  if(fixedPenalty) {
    parameters$betasdPenaltyLogLambda = logPenaltyPar
    parameters$betaTaperPenaltyLogLambda = logPenaltyPar
    parameters$betaTaperGPSPenaltyLogLambda = logPenaltyPar
    parameters$betaGammaPenaltyLogLambda = logPenaltyPar
    parameters$betaMeanPenaltyLogLambda = logPenaltyPar
    betaMeanGPSPenaltyLogLambda = logPenaltyPar
  }
  if(fixedDiffPenalty) {
    taperDiffPenaltyLogLambda = logDiffPenaltyPar
    meanDiffPenaltyLogLambda = logDiffPenaltyPar
  }
  
  # pick which parameters will be fixed
  includeTaperDiffPar =  !fixedDiffPenalty && doTaperDiffPenalty
  map = list()
  if(!includeGammaSpline)
    map = c(map, list(betaGamma=rep(factor(NA), length(betaGamma))))
  if(!doVarSpline)
    map = c(map, list(betasd=rep(factor(NA), length(betasd))))
  if(!doMeanSpline)
    map = c(map, list(betaMean=rep(factor(NA), length(betaMean))))
  if(!diffMean)
    map = c(map, list(logMeanGPS=factor(NA), betaMeanGPS=rep(factor(NA), length(betaMeanGPS))))
  if(fixInflation)
    map = c(map, list(loghighInflate=factor(NA), loglowInflate=factor(NA)))
  if(!diffGPSTaper)
    map = c(map, list(betaTaperGPS=rep(factor(NA), length(betaTaperGPS))))
  if(fixedPenalty) {
    map = c(map, list(betasdPenaltyLogLambda=factor(NA), 
                      betaTaperPenaltyLogLambda=factor(NA), 
                      betaTaperGPSPenaltyLogLambda=factor(NA),
                      betaGammaPenaltyLogLambda=factor(NA), 
                      betaMeanPenaltyLogLambda=factor(NA), 
                      betaMeanGPSPenaltyLogLambda=factor(NA)))
  }
  else if(sharedPenalty) {
    map = c(map, list(betaTaperPenaltyLogLambda=factor(NA), 
                      betaTaperGPSPenaltyLogLambda=factor(NA),
                      betaGammaPenaltyLogLambda=factor(NA), 
                      betaMeanPenaltyLogLambda=factor(NA), 
                      betaMeanGPSPenaltyLogLambda=factor(NA)))
  }
  if(!includeTaperDiffPar) {
    map = c(map, list(taperDiffPenaltyLogLambda=factor(NA)))
    map = c(map, list(meanDiffPenaltyLogLambda=factor(NA)))
  }
  if(!sharedSpatialProcess) {
    map = c(map, list(logitOmega=factor(NA)))
  }
  
  obj <- MakeADFun(data, parameters, DLL="fitModelTMB", map=map)
  obj$hessian <- TRUE
  
  # set constraints on the parameters
  # !lambdaInRange || any(sigmaZeta <= 0) || (muZeta <= 0) || (phiZeta <= 0) || (alpha <= 0.05) || (alpha >= 20))
  # feasible region given by: ui %*% theta - ci >= 0
  # for now, only restrict alpha, which is the last parameter, for numerical stability
  # logalpha - log(0.05) >= 0
  # -logalpha + log(20) >= 0 ==> logalpha <= log(20)
  nPenaltyPar = ifelse(sharedPenalty, 1, 5)
  ui = rbind(c(rep(0, length(obj$par)-1 - (2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar))), 
               1, rep(0, 2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar))), 
             c(rep(0, length(obj$par)-1 - (2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar))), 
               -1, rep(0, 2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar))))
  ci = c(log(0.05), -log(20))
  
  if(constrainMean) {
    # feasible region given by: ui %*% theta - ci >= 0
    # we need meanBasisPenalty %*% betaMean - sdBasisPenalty %*% betas >= 0
    parNames = names(obj$par)
    logMuI = 1
    meanParI = which(grepl("betaMean", parNames))
    logMeanParGPSI = which(grepl("logMeanGPS", parNames))
    meanParGPSI = which(grepl("betaMeanGPS", parNames))
    sdParI = which(grepl("betasd", parNames))
    sdParInterceptI = which(grepl("betasdIntercept", parNames))
    sdParI = sdParI[!(sdParI == sdParInterceptI)]
    meanParI = meanParI[!(meanParI %in% meanParGPSI)]
    
    if(!diffMean) {
      newui = matrix(0, nrow=nrow(meanBasisPenalty), ncol=ncol(ui))
      newui[,logMuI] = 1
      newui[,meanParI] = rbind(meanBasisPenalty)
      newui[,sdParInterceptI] = -1
      newui[,sdParI] = rbind(-sdBasisPenalty)
      ui = rbind(ui, newui)
    } else {
      newui = matrix(0, nrow=2 * nrow(meanBasisPenalty), ncol=ncol(ui))
      newui[,logMuI] = 1
      newui[,logMeanParGPSI] = c(rep(0, nrow(meanBasisPenalty)), rep(1, nrow(meanBasisPenalty)))
      newui[,meanParI] = rbind(meanBasisPenalty, meanBasisPenalty)
      newui[(nrow(meanBasisPenalty) + 1):nrow(newui),meanParGPSI] = meanBasisGPSPenalty
      newui[,sdParInterceptI] = -1
      newui[,sdParI] = rbind(-sdBasisPenalty, -sdBasisPenalty)
      ui = rbind(ui, newui)
    }
    ci = c(ci, rep(0, nrow(newui)))
  }
  obj$ui = ui
  obj$ci = ci
  
  # reparameterize inputs from optim to constrOptim
  obj$f = obj$fn
  obj$grad = obj$gr
  obj$theta = obj$par
  obj$mu=muBarrier
  
  # set optimization control parameters (adding in penalty parameters)
  controls = list(parscale=c(parscale[1:2], rep(1, (nKnotsMean - 1) * doMeanSpline + diffMean * nKnotsMeanGPS), parscale[3:length(parscale)], 
                             rep(1, nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar) + sharedSpatialProcess)))
  obj$control = controls
  finalFitControls = controls
  finalFitControls$parscale = finalFitControls$parscale / 25
  
  keepCalling = TRUE
  lastPar = obj$theta
  count = 0
  allPar = matrix(NA, nrow=maxCount, ncol=length(obj$par))
  objectiveValues = rep(NA, maxCount)
  while(keepCalling) {
    count = count + 1
    
    print("Beginning optimization at: ")
    print(obj$theta)
    
    if(useAlabama) {
      obj$hin = function(x) {c(ui %*% x-ci)}
      obj$hin.jac = function(x) {ui}
      opt = constrOptim.nl(obj$par, obj$fn, obj$gr, obj$hin, obj$hin.jac, control.optim=obj$control)
    } else {
      # opt <- do.call("optim", obj)
      opt <- do.call("constrOptim", obj)
      # opt <- do.call("myConstrOptim", obj)
    }
    
    
    grad = obj$gr()
    if(abs(max(grad, na.rm=TRUE)) < 0.01 || count >= maxCount)
      keepCalling = FALSE
    else {
      obj$par = opt$par
      obj$theta = jitter(opt$par, amount = 0.1 * 0.5^count)
      obj$control = finalFitControls
      lastPar = opt$par
    }
    allPar[count,] = opt$par
    objectiveValues[count] = opt$value
    
    # reorder output parameters
    minPar = opt$par
    optParNames = names(opt$par)
    logmuEst = minPar[which(optParNames == "logmu")]
    betaMeanEst = c(logmuEst, minPar[which(optParNames == "betaMean")])
    logMeanGPSEst = minPar[which(optParNames == "logMeanGPS")]
    betaMeanGPSEst = c(logMeanGPSEst, minPar[which(optParNames == "betaMeanGPS")])
    if(!diffMean)
      logMeanGPSEst = 0
    betasdInterceptEst = minPar[which(optParNames == "betasdIntercept")]
    betasdEst = c(betasdInterceptEst, minPar[which(optParNames == "betasd")])
    betaTaperInterceptEst = minPar[which(optParNames == "betaTaperIntercept")]
    betaTaperEst = c(betaTaperInterceptEst, minPar[which(optParNames == "betaTaper")])
    betaTaperGPSEst = c(minPar[which(optParNames == "betaTaperGPS")])
    betaGammaInterceptEst = minPar[which(optParNames == "betaGammaIntercept")]
    betaGammaEst = c(betaGammaInterceptEst, minPar[which(optParNames == "betaGamma")])
    logphiEst = minPar[which(optParNames == "logphi")]
    logalphaEst = minPar[which(optParNames == "logalpha")]
    loglowInflateEst = minPar[which(optParNames == "loglowInflate")]
    loghighInflateEst = minPar[which(optParNames == "loghighInflate")]
    if(!fixedPenalty) {
      betasdPenaltyLogLambdaEst = minPar[which(optParNames == "betasdPenaltyLogLambda")]
      if(!sharedPenalty) {
        betaTaperPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperPenaltyLogLambda")]
        betaTaperGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperGPSPenaltyLogLambda")]
        betaGammaPenaltyLogLambdaEst = minPar[which(optParNames == "betaGammaPenaltyLogLambda")]
        betaMeanPenaltyLogLambdaEst = minPar[which(optParNames == "betaMeanPenaltyLogLambdaEst")]
      }
      else {
        betaTaperPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaTaperGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaGammaPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaMeanPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      }
    } else {
      betasdPenaltyLogLambdaEst = logPenaltyPar
      betaTaperPenaltyLogLambdaEst = logPenaltyPar
      betaTaperGPSPenaltyLogLambdaEst = logPenaltyPar
      betaGammaPenaltyLogLambdaEst = logPenaltyPar
      betaMeanPenaltyLogLambdaEst = logPenaltyPar
    }
    if(includeTaperDiffPar) {
      taperDiffPenaltyLogLambdaEst = minPar[which(optParNames == "taperDiffPenaltyLogLambdaEst")]
      meanDiffPenaltyLogLambda = minPar[which(optParNames == "meanDiffPenaltyLogLambda")]
    }
    else if(fixedDiffPenalty) {
      taperDiffPenaltyLogLambdaEst = logDiffPenaltyPar
      meanDiffPenaltyLogLambda = logDiffPenaltyPar
    }
    else {
      taperDiffPenaltyLogLambdaEst = NULL
      meanDiffPenaltyLogLambda = NULL
    }
    if(sharedSpatialProcess) {
      logitOmega = minPar[which(optParNames == "logitOmega")]
      omega = expit(logitOmega)
    } else {
      omega = 0
      logitOmega = -Inf
    }
    finalPar = c(logmuEst=logmuEst, betaMeanEst=betaMeanEst, logMeanGPSEst=logMeanGPSEst, betaMeanGPSEst=betaMeanGPSEst, 
                 betasdEst=betasdEst, betaTaperEst=betaTaperEst, betaTaperGPSEst=betaTaperGPSEst, 
                 betaGammaEst=betaGammaEst, logphiEst=logphiEst, logalphaEst=logalphaEst, 
                 loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
                 betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
                 betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
                 betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
                 betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
                 betaMeanPenaltyLogLambdaEst=betaMeanPenaltyLogLambdaEst, 
                 betaMeanGPSPenaltyLogLambda=betaMeanGPSPenaltyLogLambda, 
                 taperDiffPenaltyLogLambdaEst=taperDiffPenaltyLogLambdaEst, 
                 meanDiffPenaltyLogLambda=meanDiffPenaltyLogLambda, omega=omega, 
                 logitOmega=logitOmega)
    
    lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
    lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
    lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
    sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)
    sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
    meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
    meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
    meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)
    gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
    if(!includeGammaSpline)
      gammaBasis = matrix(1, nrow=nrow(threshSlipDat))
    if(!doMeanSpline) {
      meanBasisX = matrix(1, nrow=nrow(threshSlipDat))
      meanBasisXGPS = matrix(1, nrow=nrow(threshSlipDat))
      meanBasisY = matrix(1, nrow=nrow(threshSlipDat))
    }
    
    if(!diffGPSTaper)
      taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst), dStar = dStarGPS)
    else
      taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst + lambdaBasisXGPS %*% betaTaperGPSEst), dStar = dStarGPS)
    taperVecY = taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar)
    sdVecX = exp(sdBasisX %*% betasdEst)
    sdVecY = exp(sdBasisY %*% betasdEst)
    if(diffMean)
      meanVecX = exp(meanBasisX %*% betaMeanEst + diffMean * (meanBasisXGPS %*% betaMeanGPSEst))
    else
      meanVecX = exp(meanBasisX %*% betaMeanEst)
    meanVecY = exp(meanBasisY %*% betaMeanEst)
    gammaVec = exp(gammaBasis %*% betaGammaEst)
    
    if(count == 1) {
      if(debugPlotting)
        browser()
    }
    
    # preview results with plot
    # plot lambda and taper
    latSeq = seq(latRange[1], latRange[2], l=500)
    splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
    splineMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsGPS, lats=latSeq, latRange=latRange)
    lambdaSeq = exp(splineMat %*% betaTaperEst)
    if(diffGPSTaper)
      lambdaSeqGPS = exp(splineMat %*% betaTaperEst + splineMatGPS %*% betaTaperGPSEst)
    else
      lambdaSeqGPS = lambdaSeq
    par(mfrow=c(1,2))
    plot(latSeq, lambdaSeq, type="l", ylim=range(c(lambdaSeq, 0, lambdaSeqGPS)))
    lines(latSeq, lambdaSeqGPS, col="blue")
    tmp = getTaperSpline(betaTaperEst, nKnots=nKnots, dStar=dStar, latRange=latRange, fault=fault, normalize=TRUE, logScale = TRUE)
    plotFault(fault, tmp, varRange=c(0, 1))
    
    # plot taper (on original scale and GPS scale) and standard deviation
    taperSeq = taper(10000, lambdaSeq, dStar = dStar)
    taperSeqGPS = taper(10000, lambdaSeqGPS, dStar = dStarGPS)
    par(mfrow=c(1,2))
    plot(latSeq, taperSeq, type="l", ylim=range(c(taperSeq, 0, taperSeqGPS)))
    lines(latSeq, taperSeqGPS, col="blue")
    legend("bottomright", c("Fault taper", "Locking rate taper"), col=c("black", "blue"), lty=1)
    sdMat = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, lats=latSeq, latRange=latRange)
    sdSeq = exp(sdMat %*% betasdEst)
    plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq)))
    sdSharedSeq = sdSeq * sqrt(omega)
    lines(latSeq, sdSharedSeq, lty=2)
    
    # plot gamma and mean vectors
    latSeq = seq(latRange[1], latRange[2], l=500)
    gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
    if(!includeGammaSpline)
      gammaMat = matrix(1, nrow=length(latSeq))
    gammaSeq = exp(gammaMat %*% betaGammaEst)
    meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
    meanMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)
    if(!doMeanSpline)
      meanMat = matrix(1, nrow=length(latSeq))
    meanSeq = exp(meanMat %*% betaMeanEst)
    if(!doMeanSpline)
      meanMat = matrix(1, nrow=length(latSeq))
    if(diffMean)
      meanSeqGPS = exp(meanMat %*% betaMeanEst + diffMean * (meanMatGPS %*% betaMeanGPSEst))
    else
      meanSeqGPS = exp(meanMat %*% betaMeanEst)
    plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
    plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0, meanSeqGPS)))
    lines(latSeq, meanSeqGPS, col="blue")
    
    par(mfrow=c(1,1))
    meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
    meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
    if(diffMean)
      meanVecX = exp(meanBasisX %*% betaMeanEst + diffMean * (meanBasisXGPS %*% betaMeanGPSEst))
    else
      meanVecX = exp(meanBasisX %*% betaMeanEst)
    sortI = sort(threshSlipDat$lat, index.return=TRUE)$ix
    plot(threshSlipDat$lat, 1000 * meanVecX * taperVecX / threshSlipDat$slip, pch=19, cex=.1,
         ylab="Time (Years)", xlab="Latitude", main="Estimated time between earthquakes", 
         ylim = range(c(1000 * meanVecX * taperVecX / threshSlipDat$slip, 1000 / gammaVec[sortI])))
    lines(threshSlipDat$lat[sortI], 1000 / gammaVec[sortI], col="blue")
    
    plot(threshSlipDat$lat, threshSlipDat$slip / taperVecX, pch=19, cex=.1, 
         ylim=range(c(threshSlipDat$slip / taperVecX, meanSeqGPS * gammaSeq)), 
         main="Estimated untapered locking rate")
    lines(latSeq, meanSeqGPS * gammaSeq, col="blue")
    
    if(count == 1) {
      if(debugPlotting)
        browser()
    }
  }
  minI = which(objectiveValues == min(objectiveValues, na.rm = TRUE))[1]
  minPar = allPar[minI,]
  
  # opt$hessian ## <-- FD hessian from optim
  grad = obj$gr(minPar)
  hess = obj$he(minPar)    ## <-- Analytical hessian
  report = sdreport(obj)
  
  print("optimum parameters for each iteration:")
  print(t(allPar[,1:count]))
  
  print("optimum parameters for best iteration:")
  print(minPar)
  
  print("optimum gradient for best iteration")
  print(grad)
  
  print("sdreport")
  print(report)
  
  # reorder output parameters
  optParNames = names(opt$par)
  logmuEst = minPar[which(optParNames == "logmu")]
  betaMeanEst = c(logmuEst, minPar[which(optParNames == "betaMean")])
  logMeanGPSEst = minPar[which(optParNames == "logMeanGPS")]
  betaMeanGPSEst = c(logMeanGPSEst, minPar[which(optParNames == "betaMeanGPS")])
  if(!diffMean)
    logMeanGPSEst = 0
  betasdInterceptEst = minPar[which(optParNames == "betasdIntercept")]
  betasdEst = c(betasdInterceptEst, minPar[which(optParNames == "betasd")])
  betaTaperInterceptEst = minPar[which(optParNames == "betaTaperIntercept")]
  betaTaperEst = c(betaTaperInterceptEst, minPar[which(optParNames == "betaTaper")])
  betaTaperGPSEst = c(minPar[which(optParNames == "betaTaperGPS")])
  betaGammaInterceptEst = minPar[which(optParNames == "betaGammaIntercept")]
  betaGammaEst = c(betaGammaInterceptEst, minPar[which(optParNames == "betaGamma")])
  logphiEst = minPar[which(optParNames == "logphi")]
  logalphaEst = minPar[which(optParNames == "logalpha")]
  loglowInflateEst = minPar[which(optParNames == "loglowInflate")]
  loghighInflateEst = minPar[which(optParNames == "loghighInflate")]
  if(!fixedPenalty) {
    betasdPenaltyLogLambdaEst = minPar[which(optParNames == "betasdPenaltyLogLambda")]
    if(!sharedPenalty) {
      betaTaperPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperPenaltyLogLambda")]
      betaTaperGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperGPSPenaltyLogLambda")]
      betaGammaPenaltyLogLambdaEst = minPar[which(optParNames == "betaGammaPenaltyLogLambda")]
      betaMeanPenaltyLogLambdaEst = minPar[which(optParNames == "betaMeanPenaltyLogLambdaEst")]
    }
    else {
      betaTaperPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaTaperGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaGammaPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaMeanPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
    }
  } else {
    betasdPenaltyLogLambdaEst = logPenaltyPar
    betaTaperPenaltyLogLambdaEst = logPenaltyPar
    betaTaperGPSPenaltyLogLambdaEst = logPenaltyPar
    betaGammaPenaltyLogLambdaEst = logPenaltyPar
    betaMeanPenaltyLogLambdaEst = logPenaltyPar
  }
  if(includeTaperDiffPar) {
    taperDiffPenaltyLogLambdaEst = minPar[which(optParNames == "taperDiffPenaltyLogLambdaEst")]
    meanDiffPenaltyLogLambda = minPar[which(optParNames == "meanDiffPenaltyLogLambda")]
  }
  else if(fixedDiffPenalty) {
    taperDiffPenaltyLogLambdaEst = logDiffPenaltyPar
    meanDiffPenaltyLogLambda = logDiffPenaltyPar
  }
  else {
    taperDiffPenaltyLogLambdaEst = NULL
    meanDiffPenaltyLogLambda = NULL
  }
  if(sharedSpatialProcess) {
    logitOmega = minPar[which(optParNames == "logitOmega")]
    omega = expit(logitOmega)
  } else {
    omega = 0
    logitOmega = -Inf
  }
  finalPar = c(logmuEst=logmuEst, betaMeanEst=betaMeanEst, logMeanGPSEst=logMeanGPSEst, betaMeanGPSEst=betaMeanGPSEst, 
               betasdEst=betasdEst, betaTaperEst=betaTaperEst, betaTaperGPSEst=betaTaperGPSEst, 
               betaGammaEst=betaGammaEst, logphiEst=logphiEst, logalphaEst=logalphaEst, 
               loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
               betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
               betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
               betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
               betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
               betaMeanPenaltyLogLambdaEst=betaMeanPenaltyLogLambdaEst, 
               betaMeanGPSPenaltyLogLambda=betaMeanGPSPenaltyLogLambda, 
               taperDiffPenaltyLogLambdaEst=taperDiffPenaltyLogLambdaEst, 
               meanDiffPenaltyLogLambda=meanDiffPenaltyLogLambda, omega=omega, 
               logitOmega=logitOmega)
  
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
  lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
  sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
  meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
  meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)
  gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
  if(!includeGammaSpline)
    gammaBasis = matrix(1, nrow=nrow(threshSlipDat))
  if(!doMeanSpline) {
    meanBasisX = matrix(1, nrow=nrow(threshSlipDat))
    meanBasisXGPS = matrix(1, nrow=nrow(threshSlipDat))
    meanBasisY = matrix(1, nrow=nrow(threshSlipDat))
  }
  
  if(!diffGPSTaper)
    taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst), dStar = dStarGPS)
  else
    taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst + lambdaBasisXGPS %*% betaTaperGPSEst), dStar = dStarGPS)
  taperVecY = taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar)
  sdVecX = exp(sdBasisX %*% betasdEst)
  sdVecY = exp(sdBasisY %*% betasdEst)
  if(diffMean)
    meanVecX = exp(meanBasisX %*% betaMeanEst + meanBasisXGPS %*% betaMeanGPSEst)
  else
    meanVecX = exp(meanBasisX %*% betaMeanEst)
  meanVecY = exp(meanBasisY %*% betaMeanEst)
  gammaVec = exp(gammaBasis %*% betaGammaEst)
  
  if(count == 1) {
    if(debugPlotting)
      browser()
  }
  
  # preview results with plot
  # plot lambda and taper
  latSeq = seq(latRange[1], latRange[2], l=500)
  splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
  splineMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsGPS, lats=latSeq, latRange=latRange)
  lambdaSeq = exp(splineMat %*% betaTaperEst)
  if(diffGPSTaper)
    lambdaSeqGPS = exp(splineMat %*% betaTaperEst + splineMatGPS %*% betaTaperGPSEst)
  else
    lambdaSeqGPS = lambdaSeq
  par(mfrow=c(1,2))
  plot(latSeq, lambdaSeq, type="l", ylim=range(c(lambdaSeq, 0, lambdaSeqGPS)))
  lines(latSeq, lambdaSeqGPS, col="blue")
  tmp = getTaperSpline(betaTaperEst, nKnots=nKnots, dStar=dStar, latRange=latRange, fault=fault, normalize=TRUE, logScale = TRUE)
  plotFault(fault, tmp, varRange=c(0, 1))
  
  # plot taper (on original scale and GPS scale) and standard deviation
  taperSeq = taper(10000, lambdaSeq, dStar = dStar)
  taperSeqGPS = taper(10000, lambdaSeqGPS, dStar = dStarGPS)
  par(mfrow=c(1,2))
  plot(latSeq, taperSeq, type="l", ylim=range(c(taperSeq, 0, taperSeqGPS)))
  lines(latSeq, taperSeqGPS, col="blue")
  legend("bottomright", c("Fault taper", "Locking rate taper"), col=c("black", "blue"), lty=1)
  sdMat = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, lats=latSeq, latRange=latRange)
  sdSeq = exp(sdMat %*% betasdEst)
  plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq)))
  sdSharedSeq = sdSeq * sqrt(omega)
  lines(latSeq, sdSharedSeq, lty=2)
  
  # plot gamma and mean vectors
  latSeq = seq(latRange[1], latRange[2], l=500)
  gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
  if(!includeGammaSpline)
    gammaMat = matrix(1, nrow=length(latSeq))
  gammaSeq = exp(gammaMat %*% betaGammaEst)
  meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
  meanMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)
  if(!doMeanSpline)
    meanMat = matrix(1, nrow=length(latSeq))
  meanSeq = exp(meanMat %*% betaMeanEst)
  if(diffMean)
    meanSeqGPS = exp(meanMat %*% betaMeanEst + diffMean * (meanMatGPS %*% betaMeanGPSEst))
  else
    meanSeqGPS = exp(meanMat %*% betaMeanEst)
  plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
  plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0, meanSeqGPS)))
  lines(latSeq, meanSeqGPS, col="blue")
  
  par(mfrow=c(1,1))
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
  meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
  if(diffMean)
    meanVecX = exp(meanBasisX %*% betaMeanEst + diffMean * (meanBasisXGPS %*% betaMeanGPSEst))
  else
    meanVecX = exp(meanBasisX %*% betaMeanEst)
  sortI = sort(threshSlipDat$lat, index.return=TRUE)$ix
  plot(threshSlipDat$lat, 1000 * meanVecX * taperVecX / threshSlipDat$slip, pch=19, cex=.1,
       ylab="Time (Years)", xlab="Latitude", main="Estimated time between earthquakes", 
       ylim = range(c(1000 * meanVecX * taperVecX / threshSlipDat$slip, 1000 / gammaVec[sortI])))
  lines(threshSlipDat$lat[sortI], 1000 / gammaVec[sortI], col="blue")
  
  plot(threshSlipDat$lat, threshSlipDat$slip / taperVecX, pch=19, cex=.1, 
       ylim=range(c(threshSlipDat$slip / taperVecX, meanSeqGPS * gammaSeq)), 
       main="Estimated untapered locking rate")
  lines(latSeq, meanSeqGPS * gammaSeq, col="blue")
  
  # Return results
  # return(list(Ests=Ests, muZetaEst=muZetaEst, sigmaZetaEst=sigmaZetaEst, lambdaEst=NA, 
  #             gammaEst=gammaEst, logLikEst=logLikEst, splineParEst=splinePar, phiEst=phiEst, 
  #             alphaEst=alphaEst, hess=hess, optimTable=optimTable, 
  #             tvec=tvec, tvecGPS=tvecGPS, optPar=opt$par, optGrad=optGrad, 
  #             strikeDistGps = strikeDistGps, dipDistGps = dipDistGps, 
  #             strikeDistCsz = strikeDistCsz, dipDistCsz = dipDistCsz))
  modelInfo = list(obj=obj, data=data, map=map, opt=opt, grad=grad, hess=hess, report=report, logmuEst=logmuEst, betaTaperEst=betaTaperEst, 
                   betaTaperGPSEst=betaTaperGPSEst, logphiEst=logphiEst, logalphaEst=logalphaEst, betasdEst=betasdEst, betasdEst=betasdEst, betaGammaEst=betaGammaEst, 
       betaMeanEst=betaMeanEst, logMeanGPSEst=logMeanGPSEst, betaMeanGPSEst=betaMeanGPSEst, loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
       strikeDistGps = sqrt(squareStrikeDistGps), dipDistGps = sqrt(squareDipDistGps),
       strikeDistCsz = sqrt(squareStrikeDistCsz), dipDistCsz = sqrt(squareDipDistCsz), 
       finalPar=finalPar, betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
       betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
       betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
       betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
       taperDiffPenaltyLogLambdaEst=taperDiffPenaltyLogLambdaEst, 
       taperVecX=taperVecX, taperVecY=taperVecY, sdVecX=sdVecX, sdVecY=sdVecY, gammaVec=gammaVec, 
       gpsDepthThreshold=gpsDepthThreshold, maxCount=maxCount, diffLow=diffLow, diffHigh=diffHigh, 
       penHigh=penHigh, penLow=penLow, transformation=transformation, diffMean=diffMean, 
       allInputs=allInputs)
  
  if(saveResults) {
    fileName = paste0("fit_n", nKnots, "_dS", dStar, "_diff", diffGPSTaper, 
                      "_nGPS", nKnotsGPS, "_Gam", includeGammaSpline, "_nGam", nKnotsGamma, 
                      "_dStarGPS", dStarGPS, "_sd", doVarSpline, "_nVar", nKnotsVar, 
                      "_fixPen", fixedPenalty, "_logPen", 
                      round(logPenaltyPar, 3), "_sharePen", sharedPenalty, "_Mean", doMeanSpline, 
                      "_nMu", nKnotsMean, "_nMuGPS", nKnotsMeanGPS, "_diffPen", doTaperDiffPenalty, 
                      "_logDiffPen", round(logDiffPenaltyPar, 3), "_fixDiff", fixedDiffPenalty, 
                      "_diff", round(diffLow, 3), "-", round(diffHigh, 3), 
                      "_pen", round(round(diffLow, 3), 3), "-", round(penHigh, 3), 
                      "_diffMu", diffMean, "_shared", sharedSpatialProcess, "_joint", jointShared, 
                      "_hyper", useHyperpriors, ".RData")
    save(modelInfo, file=fileName)
  }
  modelInfo
}

getInitialParameters = function(nKnots=5, nKnotsVar=5, nKnotsGamma=7, logScale=FALSE, diffGPSTaper=FALSE, 
                                nKnotsGPS=nKnots) {
  # in order: 
  # mu, five SD parameters, five taper parameters, seven gamma parameters, 
  # low inflation, high inflation, spatial range, alpha/anisotropy parameter
  if(!diffGPSTaper)
    out = c(20, 10, rep(0, nKnotsVar - 1), 1, rep(0, nKnots - 1), 1, rep(0, nKnotsGamma - 1), 1.75, 1.25, 175, 1)
  else
    out = c(20, 10, rep(0, nKnotsVar - 1), 1, rep(0, nKnotsGPS + nKnots - 1), 1, rep(0, nKnotsGamma - 1), 1.75, 1.25, 175, 1)
    
  if(logScale) {
    # mean and standard deviation intercept
    out[1:2] = log(out[1:2])
    
    # intercept for taper
    out[(2 + nKnotsVar)] = log(out[(2 + nKnotsVar)])
    
    # intercept for gamma
    out[(2 + nKnotsVar + nKnots + nKnotsGPS * diffGPSTaper)] = log(out[(2 + nKnotsVar + nKnots + nKnotsGPS * diffGPSTaper)])
    
    # inflation, scale, and anisotropy
    out[(length(out) - 3):length(out)] = log(out[(length(out) - 3):length(out)])
  }
  
  out
}

plotModelInfo = function(modelInfo, latRange=c(40, 50), fault=csz, gpsDat=slipDatCSZ) {
  # modelInfo = list(obj=obj, data=data, map=map, opt=opt, grad=grad, hess=hess, report=report, logmuEst=logmuEst, betaTaperEst=betaTaperEst, 
  #                  logphiEst=logphiEst, logalphaEst=logalphaEst, betasdEst=betasdEst, betasdEst=betasdEst, betaGammaEst=betaGammaEst, 
  #                  loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
  #                  strikeDistGps = sqrt(squareStrikeDistGps), dipDistGps = sqrt(squareDipDistGps),
  #                  strikeDistCsz = sqrt(squareStrikeDistCsz), dipDistCsz = sqrt(squareDipDistCsz), 
  #                  finalPar=finalPar, betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
  #                  betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
  #                  betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
  #                  betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
  #                  taperDiffPenaltyLogLambdaEst=taperDiffPenaltyLogLambdaEst, 
  #                  taperVecX=taperVecX, taperVecY=taperVecY, sdVecX=sdVecX, sdVecY=sdVecY, gammaVec=gammaVec, 
  #                  gpsDepthThreshold=gpsDepthThreshold, maxCount=maxCount)
  parNames = names(modelInfo$opt$par)
  logmuEst = modelInfo$logmuEst
  # betaMeanEst = modelInfo$betaMeanEst
  # betaMeanEst = c(logmuEst, modelInfo$opt$par[which(parNames == "betaMean")])
  betaMeanGPSEst = c(logmuEst, modelInfo$opt$par[which(parNames == "betaMeanGPS")])
  diffMean = length(modelInfo$opt$par[which(parNames == "betaMeanGPS")]) == 0
  diffMean = diffMean && (length(modelInfo$opt$par[which(parNames == "logMeanGPS")]) == 0)
  betaMeanEst = modelInfo$betaMeanEst
  # betaMeanGPSEst = modelInfo$betaMeanGPSEst
  # logMeanGPSEst = modelInfo$opt$par[which(parNames == "logMeanGPS")]
  if(!diffMean)
    logMeanGPSEst = 0
  betasdInterceptEst = modelInfo$betasdInterceptEst
  betasdEst = modelInfo$betasdEst
  betaTaperInterceptEst = modelInfo$betaTaperInterceptEst
  betaTaperEst = modelInfo$betaTaperEst
  betaTaperGPSEst = modelInfo$opt$par[which(parNames == "betaTaperGPS")]
  betaGammaInterceptEst = modelInfo$betaGammaInterceptEst
  betaGammaEst = modelInfo$betaGammaEst
  logphiEst = modelInfo$logphiEst
  logalphaEst = modelInfo$logalphaEst
  loglowInflateEst = modelInfo$loglowInflateEst
  loghighInflateEst = modelInfo$loghighInflateEst
  betasdPenaltyLogLambdaEst = modelInfo$betasdPenaltyLogLambdaEst
  betaTaperPenaltyLogLambdaEst = modelInfo$betaTaperPenaltyLogLambdaEst
  betaTaperGPSPenaltyLogLambdaEst = modelInfo$betaTaperGPSPenaltyLogLambdaEst
  betaGammaPenaltyLogLambdaEst = modelInfo$betaGammaPenaltyLogLambdaEst
  betaMeanPenaltyLogLambdaEst = modelInfo$betaMeanPenaltyLogLambdaEst
  taperDiffPenaltyLogLambdaEst = modelInfo$taperDiffPenaltyLogLambdaEst
  
  gpsDepthThreshold = modelInfo$gpsDepthThreshold
  dStar = modelInfo$data$dStar
  dStarGPS = modelInfo$data$dStarGPS
  
  nKnots = length(betaTaperEst)
  nKnotsGPS = length(betaTaperGPSEst)
  nKnotsGamma = length(betaGammaEst)
  nKnotsVar = length(betasdEst)
  nKnotsMean = length(betaMeanEst)
  nKnotsMeanGPS = length(betaMeanGPSEst)
  
  diffGPSTaper = length(betaTaperGPSEst) != 0
  includeGammaSpline = length(betaGammaEst) != 1
  doMeanSpline = length(betaMeanEst) != 1
  
  # threshold the gps data and multiplied standard deviations by three as in the literature
  threshSlipDat = gpsDat[gpsDat$Depth<gpsDepthThreshold,]
  threshSlipDat$slipErr = threshSlipDat$slipErr*3
  
  taperVecX = modelInfo$taperVecX
  gammaVec = modelInfo$gammaVec
  
  # preview results with plot
  # plot lambda and taper
  latSeq = seq(latRange[1], latRange[2], l=500)
  splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
  splineMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsGPS, lats=latSeq, latRange=latRange)
  lambdaSeq = exp(splineMat %*% betaTaperEst)
  if(diffGPSTaper)
    lambdaSeqGPS = exp(splineMat %*% betaTaperEst + splineMatGPS %*% betaTaperGPSEst)
  else
    lambdaSeqGPS = lambdaSeq
  par(mfrow=c(1,2))
  plot(latSeq, lambdaSeq, type="l", ylim=range(c(lambdaSeq, 0, lambdaSeqGPS)))
  lines(latSeq, lambdaSeqGPS, col="blue")
  tmp = getTaperSpline(betaTaperEst, nKnots=nKnots, dStar=dStar, latRange=latRange, fault=fault, normalize=TRUE, logScale = TRUE)
  plotFault(fault, tmp, varRange=c(0, 1))
  
  # plot taper (on original scale and GPS scale) and standard deviation
  taperSeq = taper(10000, lambdaSeq, dStar = dStar)
  taperSeqGPS = taper(10000, lambdaSeqGPS, dStar = dStarGPS)
  par(mfrow=c(1,2))
  plot(latSeq, taperSeq, type="l", ylim=range(c(taperSeq, 0, taperSeqGPS)))
  lines(latSeq, taperSeqGPS, col="blue")
  legend("bottomright", c("Fault taper", "Locking rate taper"), col=c("black", "blue"), lty=1)
  sdMat = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, lats=latSeq, latRange=latRange)
  sdSeq = exp(sdMat %*% betasdEst)
  plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq)))
  
  # plot gamma and mean vectors
  latSeq = seq(latRange[1], latRange[2], l=500)
  gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
  if(!includeGammaSpline)
    gammaMat = matrix(1, nrow=length(latSeq))
  gammaSeq = exp(gammaMat %*% betaGammaEst)
  meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
  meanMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)
  if(!doMeanSpline)
    meanMat = matrix(1, nrow=length(latSeq))
  meanSeq = exp(meanMat %*% betaMeanEst)
  meanSeqGPS = exp(meanMat %*% betaMeanEst + diffMean * (meanMatGPS %*% betaMeanGPSEst))
  plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
  plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0, meanSeqGPS)))
  lines(latSeq, meanSeqGPS, col="blue")
  
  par(mfrow=c(1,1))
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
  meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
  meanVecX = exp(meanBasisX %*% betaMeanEst + diffMean * (meanBasisXGPS %*% betaMeanGPSEst))
  sortI = sort(threshSlipDat$lat, index.return=TRUE)$ix
  plot(threshSlipDat$lat, 1000 * meanVecX * taperVecX / threshSlipDat$slip, pch=19, cex=.1,
       ylab="Time (Years)", xlab="Latitude", main="Estimated time between earthquakes", 
       ylim = range(c(1000 * meanVecX * taperVecX / threshSlipDat$slip, 1000 / gammaVec[sortI])))
  lines(threshSlipDat$lat[sortI], 1000 / gammaVec[sortI], col="blue")
  
  plot(threshSlipDat$lat, threshSlipDat$slip / taperVecX, pch=19, cex=.1, 
       ylim=range(c(threshSlipDat$slip / taperVecX, meanSeqGPS * gammaSeq)), 
       main="Estimated untapered locking rate")
  lines(latSeq, meanSeqGPS * gammaSeq, col="blue")
}

# fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doTaperDiffPenalty = TRUE, 
#                       G=G, debugPlotting=TRUE, logPenaltyPar=log(1), logDiffPenaltyPar=log(1), 
#                       sharedSpatialProcess=TRUE, jointShared = TRUE, debug=TRUE)


logit <- function(x) {
  log(x/(1-x))
}

expit <- function(x) {
  res = exp(x)/(1+exp(x))
  res[x > 100] = 1
  res[x < -100] = 0
  res
}




