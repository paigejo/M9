library(TMB)
library(alabama)
# source("setup.R")

# fits the combined model with normalized taper using TMB jointly with the variance inflation 
# parameters as well as the gamma spline parameters
fitModelTMB = function(modelInfo=NULL, paramInit=NULL, gpsDat=slipDatCSZ, gpsDepthThreshold=21000, 
                       G=NULL, subDat=dr1, fault=csz, nKnots=5, 
                       dStar=25000, maxit=500, latRange=c(40, 50), 
                       diffGPSTaper=TRUE, nKnotsGPS=nKnots, reltol=1e-8, 
                       nKnotsVar=5, includeGammaSpline=TRUE, nKnotsGamma=7, 
                       dStarGPS=40000, seed=123, debug=FALSE, doVarSpline=TRUE, 
                       fixInflation=FALSE, maxCount=1, fixedPenalty=TRUE, 
                       logPenaltyPar=log(10), sharedPenalty=FALSE, doMeanSpline=TRUE, 
                       nKnotsMean=ifelse(doMeanSpline, 5, 1), nKnotsMeanGPS=ifelse(doMeanSpline, 7, 1), 
                       doSmoothnessPenalty=FALSE, doDiffPenalty=FALSE, logDiffPenaltyPar=log(1), 
                       fixedDiffPenalty=FALSE, debugPlotting=TRUE, saveResults=TRUE, 
                       diffLow=log(1), diffHigh=log(100), penLow=log(1), penHigh=log(100), 
                       diffMean=FALSE, constrainMean=TRUE, muBarrier=1e-04, useAlabama=TRUE, 
                       useHyperpriors=FALSE, sharedSpatialProcess=FALSE, jointShared=FALSE, 
                       inflateVarLocking=FALSE, diffVar=FALSE, nKnotsVarGPS=1, estimateGpsShared=FALSE, 
                       conditionalGamma=FALSE, calcHess=FALSE, recompile=FALSE, varSmoothnessPenalty=FALSE, 
                       varLogLambda=NULL, reparameterizeVar=FALSE) {
  
  if(is.null(modelInfo)) {
    # get all input values
    allInputs = as.list(environment())
    
    # threshold the gps data and multiplied standard deviations by three as in the literature
    threshSlipDat = gpsDat[gpsDat$Depth<gpsDepthThreshold,]
    if(!inflateVarLocking)
      threshSlipDat$slipErr = threshSlipDat$slipErr*3
    
    # get input parameters and separate them
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
    if(!reparameterizeVar)
      betasdInterceptGPS = 0
    else
      betasdInterceptGPS = betasdIntercept
    betasdGPS = rep(0, nKnotsVarGPS - 1)
    betaGamma = out$gammaPar
    betaGammaIntercept = betaGamma[1]
    betaGamma = betaGamma[-1]
    loglowInflate = out$lowInflation
    loghighInflate = out$highInflation
    parscale = out$parscale
    parNames = out$parNames
    parNamesTMB = out$parNamesTMB
    logitOmega = logit(0.25) # initially, assume most of the variation is taken up by iid spatial processes rather than shared
    logitOmega2 = logitOmega
    logLockInflate = log(3^2)
    
    # determine which of the input parameters will be optimized
    includeI = c(TRUE, TRUE, rep(doVarSpline, nKnotsVar - 1), rep(TRUE, nKnots), rep(TRUE, nKnotsGPS * diffGPSTaper), !conditionalGamma, 
                 rep(as.numeric(includeGammaSpline)*as.numeric(!conditionalGamma), nKnotsGamma - 1), rep(!fixInflation, 2), TRUE, TRUE)
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
    sdBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVarGPS, latRange=latRange)[,-1]
    sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)[,-1]
    if(!conditionalGamma) {
      gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)[,-1]
    }
    else {
      gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
    }
    latSeq = seq(latRange[1], latRange[2], l=500)
    deltaPenalty = latSeq[2] - latSeq[1]
    meanBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)[,-1]
    meanBasisGPSPenalty = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)[,-1]
    sdBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnotsVar, lats=latSeq, latRange=latRange)[,-1]
    sdBasisGPSPenalty = getSplineBasis(fault=fault, nKnots=nKnotsVarGPS, lats=latSeq, latRange=latRange)[,-1]
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
    if(!conditionalGamma) {
      if(!debug)
        compile("fitModelTMB.cpp")
      else
        compile("fitModelTMB.cpp","-O0 -g")
      dyn.load(dynlib("fitModelTMB"))
    }
    else {
      if(!debug)
        compile("fitModelTMBconditionalGamma.cpp")
      else
        compile("fitModelTMBconditionalGamma.cpp","-O0 -g")
      dyn.load(dynlib("fitModelTMBconditionalGamma"))
    }
    set.seed(seed)
    data <- list(y=y, x=x, lambdaBasisY=lambdaBasisY, lambdaBasisX=lambdaBasisX, lambdaBasisXGPS=lambdaBasisXGPS, DSStrikeCSZ=DSStrikeCSZ, 
                 DSDipCSZ=DSDipCSZ, DSStrikeGPS=DSStrikeGPS, DSDipGPS=DSDipGPS, zeroMask=zeroMask, lowI=lowI, 
                 ysd=ysd, xsd=xsd, xDepths=xDepths, faultDepths=faultDepths, dStar=dStar, sdBasisX=sdBasisX, 
                 sdBasisY=sdBasisY, gammaBasis=gammaBasis, G=G, dStarGPS=dStarGPS, deltaPenalty=deltaPenalty, 
                 sdBasisPenalty=sdBasisPenalty, sdBasisGPSPenalty=sdBasisGPSPenalty, taperBasisPenalty=taperBasisPenalty, taperBasisGPSPenalty=taperBasisGPSPenalty, 
                 gammaBasisPenalty=gammaBasisPenalty, penaltyMean=penaltyMean, penaltySD=penaltySD, 
                 sharedPenalty=as.numeric(sharedPenalty), meanBasisX=meanBasisX, meanBasisXGPS=meanBasisXGPS, meanBasisY=meanBasisY, 
                 meanBasisPenalty=meanBasisPenalty, meanBasisGPSPenalty=meanBasisGPSPenalty, doDiffPenalty=as.numeric(doDiffPenalty), 
                 diffPenaltyMean=diffPenaltyMean, diffPenaltySD=diffPenaltySD, diffMean=as.numeric(diffMean), 
                 useHyperpriors=as.numeric(useHyperpriors), sharedSpatialProcess=as.numeric(sharedSpatialProcess), 
                 jointShared=as.numeric(jointShared), DSStrikeCross=squareStrikeDistCross, 
                 DSDipCross=squareDipDistCross, diffGPSTaper=as.numeric(diffGPSTaper), sdBasisXGPS=sdBasisXGPS, 
                 diffVar=as.numeric(diffVar), doSmoothnessPenalty=as.numeric(doSmoothnessPenalty), 
                 estimateGpsShared=as.numeric(estimateGpsShared), varSmoothnessPenalty=as.numeric(varSmoothnessPenalty), 
                 reparameterizeVar=as.numeric(reparameterizeVar))
    parameters = list(logmu=logmu, betasdIntercept=betasdIntercept, betasd=betasd, betasdInterceptGPS=betasdInterceptGPS, 
                      betasdGPS=betasdGPS, betaMean=betaMean, logMeanGPS=logMeanGPS, betaMeanGPS=betaMeanGPS, 
                      betaTaperIntercept=betaTaperIntercept, betaTaper=betaTaper, betaTaperGPS=betaTaperGPS, 
                      betaGammaIntercept=betaGammaIntercept, betaGamma=betaGamma, 
                      loglowInflate=loglowInflate, loghighInflate=loghighInflate, 
                      logphi=logphi, logalpha=logalpha, 
                      betasdPenaltyLogLambda=0, betasdGPSPenaltyLogLambda=0, 
                      betaTaperPenaltyLogLambda=0, betaTaperGPSPenaltyLogLambda=0, 
                      betaGammaPenaltyLogLambda=0, betaMeanPenaltyLogLambda=0, 
                      taperDiffPenaltyLogLambda=0, betaMeanGPSPenaltyLogLambda=0, 
                      meanDiffPenaltyLogLambda=0, sdDiffPenaltyLogLambda=0, 
                      logitOmega=logitOmega, logitOmega2=logitOmega2, logLockInflate=logLockInflate)
    if(fixedPenalty) {
      parameters$betasdPenaltyLogLambda = logPenaltyPar
      parameters$betasdGPSPenaltyLogLambda = logPenaltyPar
      parameters$betaTaperPenaltyLogLambda = logPenaltyPar
      parameters$betaTaperGPSPenaltyLogLambda = logPenaltyPar
      parameters$betaGammaPenaltyLogLambda = logPenaltyPar
      parameters$betaMeanPenaltyLogLambda = logPenaltyPar
      parameters$betaMeanGPSPenaltyLogLambda = logPenaltyPar
    }
    if(!is.null(varLogLambda)) {
      parameters$betasdPenaltyLogLambda = varLogLambda
      parameters$betasdGPSPenaltyLogLambda = varLogLambda
    }
    if(fixedDiffPenalty) {
      taperDiffPenaltyLogLambda = logDiffPenaltyPar
      meanDiffPenaltyLogLambda = logDiffPenaltyPar
      sdDiffPenaltyLogLambda = logDiffPenaltyPar
    }
    if(conditionalGamma) {
      parameters$betaGammaIntercept = NULL
      parameters$betaGamma = NULL
      parameters$betaGammaPenaltyLogLambda = NULL
    }
    
    # pick which parameters will be fixed
    includeTaperDiffPar =  !fixedDiffPenalty && doDiffPenalty
    map = list()
    if(!includeGammaSpline) {
      if(!conditionalGamma)
        map = c(map, list(betaGamma=rep(factor(NA), length(betaGamma))))
    }
    if(!doVarSpline)
      map = c(map, list(betasd=rep(factor(NA), length(betasd))))
    if(!doMeanSpline)
      map = c(map, list(betaMean=rep(factor(NA), length(betaMean))))
    if(!diffMean)
      map = c(map, list(logMeanGPS=factor(NA), betaMeanGPS=rep(factor(NA), length(betaMeanGPS))))
    if(!diffVar)
      map = c(map, list(betasdInterceptGPS=factor(NA), betasdGPS=rep(factor(NA), length(betasdGPS))))
    if(fixInflation)
      map = c(map, list(loghighInflate=factor(NA), loglowInflate=factor(NA)))
    if(!diffGPSTaper)
      map = c(map, list(betaTaperGPS=rep(factor(NA), length(betaTaperGPS))))
    if(fixedPenalty || !doSmoothnessPenalty) {
      map = c(map, list(betasdPenaltyLogLambda=factor(NA),
                        betasdGPSPenaltyLogLambda=factor(NA), 
                        betaTaperPenaltyLogLambda=factor(NA), 
                        betaTaperGPSPenaltyLogLambda=factor(NA),
                        betaMeanPenaltyLogLambda=factor(NA), 
                        betaMeanGPSPenaltyLogLambda=factor(NA)))
      if(!conditionalGamma)
        map = c(map, list(betaGammaPenaltyLogLambda=factor(NA)))
    }
    else if(sharedPenalty) {
      map = c(map, list(betasdGPSPenaltyLogLambda=factor(NA), 
                        betaTaperPenaltyLogLambda=factor(NA), 
                        betaTaperGPSPenaltyLogLambda=factor(NA),
                        betaMeanPenaltyLogLambda=factor(NA), 
                        betaMeanGPSPenaltyLogLambda=factor(NA)))
      if(!conditionalGamma)
        map = c(map, list(betaGammaPenaltyLogLambda=factor(NA)))
    }
    if(!includeTaperDiffPar) {
      map = c(map, list(taperDiffPenaltyLogLambda=factor(NA)))
      map = c(map, list(meanDiffPenaltyLogLambda=factor(NA)))
      map = c(map, list(sdDiffPenaltyLogLambda=factor(NA)))
    }
    if(!sharedSpatialProcess) {
      map = c(map, list(logitOmega=factor(NA), logitOmega2=factor(NA)))
    }
    if(!jointShared || !estimateGpsShared) {
      map = c(map, list(logitOmega2=factor(NA)))
    }
    if(!inflateVarLocking) {
      map = c(map, list(logLockInflate=factor(NA)))
    }
    
    if(!conditionalGamma)
      obj <- MakeADFun(data, parameters, DLL="fitModelTMB", map=map)
    else
      obj <- MakeADFun(data, parameters, DLL="fitModelTMBconditionalGamma", map=map)
    obj$hessian <- FALSE
    
    # reparameterize inputs from optim to constrOptim
    obj$f = obj$fn
    obj$grad = obj$gr
    obj$theta = obj$par
    obj$mu=muBarrier
    
    # set optimization control parameters (adding in penalty parameters)
    nPenaltyPar = ifelse(sharedPenalty, 1, 5 - as.numeric(conditionalGamma))
    controls = list(parscale=c(parscale[1:2], rep(1, (nKnotsMean - 1) * as.numeric(doMeanSpline) + as.numeric(diffMean) * nKnotsMeanGPS + as.numeric(diffVar) * nKnotsVarGPS), parscale[3:length(parscale)], 
                               rep(1, nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar) + as.numeric(sharedSpatialProcess) + 
                                     as.numeric(sharedSpatialProcess && jointShared && estimateGpsShared) + as.numeric(inflateVarLocking))))
    obj$control = controls
  }
  else {
    set.seed(seed)
    
    data = modelInfo$data
    parameters = modelInfo$parameters
    map = modelInfo$map
    ui = modelInfo$obj$ui
    ci = modelInfo$obj$ci
    controls = modelInfo$obj$control
    
    # compile the function and its gradient
    if(recompile) {
      if(!conditionalGamma) {
        if(!debug)
          compile("fitModelTMB.cpp")
        else
          compile("fitModelTMB.cpp","-O0 -g")
        dyn.load(dynlib("fitModelTMB"))
      }
      else {
        if(!debug)
          compile("fitModelTMBconditionalGamma.cpp")
        else
          compile("fitModelTMBconditionalGamma.cpp","-O0 -g")
        dyn.load(dynlib("fitModelTMBconditionalGamma"))
      }
      
      if(fixedPenalty) {
        parameters$betasdPenaltyLogLambda = logPenaltyPar
        parameters$betasdGPSPenaltyLogLambda = logPenaltyPar
        parameters$betaTaperPenaltyLogLambda = logPenaltyPar
        parameters$betaTaperGPSPenaltyLogLambda = logPenaltyPar
        parameters$betaGammaPenaltyLogLambda = logPenaltyPar
        parameters$betaMeanPenaltyLogLambda = logPenaltyPar
        parameters$betaMeanGPSPenaltyLogLambda = logPenaltyPar
      }
      if(fixedDiffPenalty) {
        parameters$taperDiffPenaltyLogLambda = logDiffPenaltyPar
        parameters$meanDiffPenaltyLogLambda = logDiffPenaltyPar
        parameters$sdDiffPenaltyLogLambda = logDiffPenaltyPar
      }
      if(!is.null(varLogLambda)) {
        parameters$betasdPenaltyLogLambda = varLogLambda
        parameters$betasdGPSPenaltyLogLambda = varLogLambda
      }
      
      if(!conditionalGamma)
        obj <- MakeADFun(data, parameters, DLL="fitModelTMB", map=map)
      else
        obj <- MakeADFun(data, parameters, DLL="fitModelTMBconditionalGamma", map=map)
      
      obj$hessian <- FALSE
      
      # reparameterize inputs from optim to constrOptim
      obj$f = obj$fn
      obj$grad = obj$gr
      obj$theta = obj$par
      obj$mu=muBarrier
    }
    else {
      obj  = modelInfo$obj
    }
    
    doDiffPenalty = data$doDiffPenalty
    doSmoothnessPenalty = data$doSmoothnessPenalty
    includeTaperDiffPar =  !fixedDiffPenalty && doDiffPenalty
    
    xDepths = data$xDepths
    faultDepths = data$faultDepths
    lambdaBasisX = data$lambdaBasisX
    lambdaBasisY = data$lambdaBasisY
    sdBasisX = data$sdBasisX
    sdBasisXGPS = data$sdBasisXGPS
    sdBasisY = data$sdBasisY
    dStar = data$dStar
    dStarGPS = data$dStarGPS
    gammaBasis = data$gammaBasis
    deltaPenalty = data$deltaPenalty
    sdBasisPenalty = data$sdBasisPenalty
    sdBasisGPSPenalty = data$sdBasisGPSPenalty
    taperBasisPenalty = data$taperBasisPenalty
    taperBasisGPSPenalty = data$taperBasisGPSPenalty
    gammaBasisPenalty = data$gammaBasisPenalty
    meanBasisX = data$meanBasisX
    meanBasisY = data$meanBasisY
    meanBasisXGPS = data$meanBasisXGPS
    meanBasisPenalty = data$meanBasisPenalty
    meanBasisGPSPenalty = data$meanBasisGPSPenalty
    diffGPSTaper = data$diffGPSTaper
    sharedPenalty = data$sharedPenalty
    fixInflation = modelInfo$allInputs$fixInflation
    inflateVarLocking = modelInfo$allInputs$inflateVarLocking
    
    squareStrikeDistGps = data$DSStrikeGPS^2
    squareDipDistGps = data$DSDipGPS^2
    squareStrikeDistCsz = data$DSStrikeCSZ^2
    squareDipDistCsz = data$DSDipCSZ^2
    squareStrikeDistCross = data$strikeDistCross^2
    squareDipDistCross = data$dipDistCross^2
    transformation = modelInfo$transformation
    
    diffMean = modelInfo$diffMean
    diffVar = modelInfo$diffVar
    diffHigh = modelInfo$diffHigh
    penHigh = modelInfo$penHigh
    gpsDepthThreshold = modelInfo$gpsDepthThreshold
    allInputs = modelInfo$allInputs
  }
  
  # set constraints on the parameters
  # !lambdaInRange || any(sigmaZeta <= 0) || (muZeta <= 0) || (phiZeta <= 0) || (alpha <= 0.05) || (alpha >= 20))
  # feasible region given by: ui %*% theta - ci >= 0
  # for now, only restrict alpha, which is the last parameter, for numerical stability
  # logalpha - log(0.05) >= 0
  # -logalpha + log(20) >= 0 ==> logalpha <= log(20)
  nPenaltyPar = ifelse(sharedPenalty, 1, 5)
  ui = rbind(c(rep(0, length(obj$par)-1 - (2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar) + as.numeric(inflateVarLocking))), 
               1, rep(0, 2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar) + as.numeric(inflateVarLocking))), 
             c(rep(0, length(obj$par)-1 - (2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar) + as.numeric(inflateVarLocking))), 
               -1, rep(0, 2 * as.numeric(!fixInflation) + nPenaltyPar * as.numeric(!fixedPenalty) + 2 * as.numeric(includeTaperDiffPar) + as.numeric(inflateVarLocking))))
  ci = c(log(0.05), -log(20))
  
  if(constrainMean) {
    # feasible region given by: ui %*% theta - ci >= 0
    # we need meanBasisPenalty %*% betaMean - sdBasisPenalty %*% betas >= 0 
    # (ensure that mean is larger than standard deviation at all latitudes)
    parNames = names(obj$par)
    
    # get indices of the parameters corresponding to the mean
    logMuI = 1
    meanParI = which(grepl("betaMean", parNames))
    logMeanParGPSI = which(grepl("logMeanGPS", parNames))
    meanParGPSI = which(grepl("betaMeanGPS", parNames))
    meanParI = meanParI[!(meanParI %in% meanParGPSI)]
    
    # get indices of the parameters corresponding to the sd
    sdParI = which(grepl("betasd", parNames))
    sdParInterceptI = which(grepl("betasdIntercept", parNames))
    sdParI = sdParI[!(sdParI == sdParInterceptI)]
    sdParInterceptGPSI = which(grepl("betasdInterceptGPS", parNames))
    sdParInterceptI = sdParInterceptI[!(sdParInterceptI == sdParInterceptGPSI)]
    sdParGPSI = which(grepl("betasdGPS", parNames))
    sdParI = sdParI[!(sdParI %in% sdParGPSI)]
    
    # add in the linear constraints
    if(!(diffMean || diffVar)) {
      newui = matrix(0, nrow=nrow(meanBasisPenalty), ncol=ncol(ui))
      newui[,logMuI] = 1
      newui[,meanParI] = rbind(meanBasisPenalty)
      newui[,sdParInterceptI] = -1
      newui[,sdParI] = rbind(-sdBasisPenalty)
      
      ui = rbind(ui, newui)
    } else {
      newui = matrix(0, nrow=2 * nrow(meanBasisPenalty), ncol=ncol(ui))
      
      # add in constraint for the fault (and part of the constraint for the locking rate data)
      newui[, logMuI] = 1
      newui[, meanParI] = rbind(meanBasisPenalty, meanBasisPenalty)
      
      newui[, sdParInterceptI] = -1
      newui[, sdParI] = rbind(-sdBasisPenalty, -sdBasisPenalty)
      
      # add in constraint for the locking rate data
      if(diffMean) {
        newui[(nrow(meanBasisPenalty) + 1):nrow(newui), logMeanParGPSI] = rep(1, nrow(meanBasisGPSPenalty))
        newui[(nrow(meanBasisPenalty) + 1):nrow(newui), meanParGPSI] = meanBasisGPSPenalty
      }
      
      if(diffVar) {
        newui[(nrow(meanBasisPenalty) + 1):nrow(newui), sdParInterceptGPSI] = -1
        newui[(nrow(meanBasisPenalty) + 1):nrow(newui), sdParGPSI] = -sdBasisGPSPenalty
      }
      
      ui = rbind(ui, newui)
    }
    ci = c(ci, rep(0, nrow(newui)))
  }
  
  obj$ui = ui
  obj$ci = ci
  
  # initialize optimization at the given parameters if necessary
  if(!is.null(paramInit)) {
    obj$par = paramInit
    obj$theta = paramInit
  }
  
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
      obj$par = opt$par* jitter(rep(1, length(opt$par)), amount = 0.05 * 0.5^count)
      obj$theta = obj$par
      obj$control = finalFitControls
      lastPar = opt$par
      obj$control$parscale = finalFitControls$parscale*2^(-count)
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
    betasdInterceptGPSEst = minPar[which(optParNames == "betasdInterceptGPS")]
    betasdGPSEst = c(betasdInterceptGPSEst, minPar[which(optParNames == "betasdGPS")])
    betaTaperInterceptEst = minPar[which(optParNames == "betaTaperIntercept")]
    betaTaperEst = c(betaTaperInterceptEst, minPar[which(optParNames == "betaTaper")])
    betaTaperGPSEst = c(minPar[which(optParNames == "betaTaperGPS")])
    betaGammaInterceptEst = minPar[which(optParNames == "betaGammaIntercept")]
    betaGammaEst = c(betaGammaInterceptEst, minPar[which(optParNames == "betaGamma")])
    if(conditionalGamma) {
      ## use log-transformed OLS estimate
      # calculate the mean
      meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=ifelse(doMeanSpline, nKnotsMean, 1), latRange=latRange)
      meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
      
      if(diffMean)
        meanVecX = exp(meanBasisX %*% betaMeanEst + as.numeric(diffMean) * (meanBasisXGPS %*% betaMeanGPSEst))
      else
        meanVecX = exp(meanBasisX %*% betaMeanEst)
      
      # calculate the taper
      lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
      lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
      
      if(!diffGPSTaper)
        taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst), dStar = dStarGPS)
      else
        taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst + lambdaBasisXGPS %*% betaTaperGPSEst), dStar = dStarGPS)
      
      # now perform the regression
      logX = log(data$x)
      newVarx = (data$xsd / data$x)^2
      offset = log(meanVecX) + log(taperVecX)
      gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
      if(!includeGammaSpline)
        gammaBasis = matrix(1, nrow=nrow(threshSlipDat))
      mod = lm(logX~gammaBasis - 1, offset=offset, weights=1 / newVarx)
      betaGammaEst = coef(mod)
    }
    logphiEst = minPar[which(optParNames == "logphi")]
    logalphaEst = minPar[which(optParNames == "logalpha")]
    loglowInflateEst = minPar[which(optParNames == "loglowInflate")]
    loghighInflateEst = minPar[which(optParNames == "loghighInflate")]
    if(!fixedPenalty && doSmoothnessPenalty) {
      betasdPenaltyLogLambdaEst = minPar[which(optParNames == "betasdPenaltyLogLambda")]
      if(!sharedPenalty) {
        betaTaperPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperPenaltyLogLambda")]
        betaTaperGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperGPSPenaltyLogLambda")]
        betaGammaPenaltyLogLambdaEst = minPar[which(optParNames == "betaGammaPenaltyLogLambda")]
        betaMeanPenaltyLogLambdaEst = minPar[which(optParNames == "betaMeanPenaltyLogLambdaEst")]
        betasdGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betasdGPSPenaltyLogLambda")]
        betaMeanGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betaMeanGPSPenaltyLogLambda")]
      }
      else {
        betaTaperPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaTaperGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaGammaPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaMeanPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaMeanGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betasdGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      }
    } else {
      betasdPenaltyLogLambdaEst = logPenaltyPar
      betaTaperPenaltyLogLambdaEst = logPenaltyPar
      betaTaperGPSPenaltyLogLambdaEst = logPenaltyPar
      betaGammaPenaltyLogLambdaEst = logPenaltyPar
      betaMeanPenaltyLogLambdaEst = logPenaltyPar
      betaMeanGPSPenaltyLogLambdaEst = logPenaltyPar
      betasdGPSPenaltyLogLambdaEst = logPenaltyPar
    }
    if(!is.null(varLogLambda)) {
      betasdPenaltyLogLambdaEst = varLogLambda
      betasdGPSPenaltyLogLambdaEst = varLogLambda
    }
    if(includeTaperDiffPar) {
      taperDiffPenaltyLogLambdaEst = minPar[which(optParNames == "taperDiffPenaltyLogLambdaEst")]
      meanDiffPenaltyLogLambdaEst = minPar[which(optParNames == "meanDiffPenaltyLogLambda")]
      sdDiffPenaltyLogLambdaEst = minPar[which(optParNames == "sdDiffPenaltyLogLambda")]
    }
    else if(fixedDiffPenalty) {
      taperDiffPenaltyLogLambdaEst = logDiffPenaltyPar
      meanDiffPenaltyLogLambdaEst = logDiffPenaltyPar
      sdDiffPenaltyLogLambdaEst = logDiffPenaltyPar
    }
    else {
      taperDiffPenaltyLogLambdaEst = NULL
      meanDiffPenaltyLogLambdaEst = NULL
      sdDiffPenaltyLogLambdaEst = NULL
    }
    if(sharedSpatialProcess) {
      logitOmegaEst = minPar[which(optParNames == "logitOmega")]
      omegaEst = expit(logitOmegaEst)
      
      if(estimateGpsShared) {
        logitOmega2Est = minPar[which(optParNames == "logitOmega2")]
        omega2Est = expit(logitOmega2Est)
      }
      else {
        logitOmega2Est = logitOmegaEst
        omega2Est = omegaEst
      }
    } else {
      omegaEst = 0
      logitOmegaEst = -Inf
      omega2Est = 0
      logitOmega2Est = -Inf
    }
    if(inflateVarLocking) {
      logLockInflateEst = minPar[which(optParNames == "logLockInflate")]
      lockInflateEst = exp(logLockInflateEst)
    } else {
      logLockInflateEst = log(3^2)
      lockInflateEst = exp(logLockInflateEst)
    }
    finalPar = c(logmuEst=logmuEst, betaMeanEst=betaMeanEst, logMeanGPSEst=logMeanGPSEst, betaMeanGPSEst=betaMeanGPSEst, 
                 betasdEst=betasdEst, betasdGPSEst = betasdGPSEst, betaTaperEst=betaTaperEst, betaTaperGPSEst=betaTaperGPSEst, 
                 betaGammaEst=betaGammaEst, logphiEst=logphiEst, logalphaEst=logalphaEst, 
                 loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
                 betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
                 betasdGPSPenaltyLogLambdaEst=betasdGPSPenaltyLogLambdaEst, 
                 betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
                 betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
                 betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
                 betaMeanPenaltyLogLambdaEst=betaMeanPenaltyLogLambdaEst, 
                 betaMeanGPSPenaltyLogLambdaEst=betaMeanGPSPenaltyLogLambdaEst, 
                 taperDiffPenaltyLogLambdaEst=taperDiffPenaltyLogLambdaEst, 
                 meanDiffPenaltyLogLambdaEst=meanDiffPenaltyLogLambdaEst, 
                 sdDiffPenaltyLogLambdaEst = sdDiffPenaltyLogLambdaEst, 
                 omegaEst=omegaEst, logitOmegaEst=logitOmegaEst, omega2Est=omega2Est, logitOmega2Est=logitOmega2Est, 
                 logLockInflateEst=logLockInflateEst, 
                 lockInflateEst=lockInflateEst)
    
    lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
    lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
    lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
    sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)
    sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
    sdBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVarGPS, latRange=latRange)
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
    if(doVarSpline) {
      if(diffVar) {
        if(!reparameterizeVar)
          sdVecX = exp(sdBasisX %*% betasdEst + sdBasisXGPS %*% betasdGPSEst)
        else
          sdVecX = exp(sdBasisXGPS %*% betasdGPSEst)
      }
      else
        sdVecX = exp(sdBasisX %*% betasdEst)
      sdVecY = exp(sdBasisY %*% betasdEst)
    }
    else {
      if(diffVar) {
        if(!reparameterizeVar)
          sdVecX = rep(exp(betasdEst + betasdGPSEst), length(xDepths))
        else
          sdVecX = rep(exp(betasdGPSEst), length(xDepths))
      }
      else
        sdVecX = rep(exp(betasdEst ), length(xDepths))
      sdVecY = rep(exp(betasdEst), length(faultDepths))
    }
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
    sdMatGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVarGPS, lats=latSeq, latRange=latRange)
    if(doVarSpline) {
      sdSeq = exp(sdMat %*% betasdEst)
      if(diffVar) {
        if(!reparameterizeVar)
          sdSeqGPS = exp(sdMat %*% betasdEst + sdMatGPS %*% betasdGPSEst)
        else
          sdSeqGPS = exp(sdMatGPS %*% betasdGPSEst)
      }
      else
        sdSeqGPS = exp(sdMat %*% betasdEst)
    }
    else {
      sdSeq = rep(exp(betasdEst), length(latSeq))
      if(diffVar) {
        if(!reparameterizeVar)
          sdSeqGPS = rep(exp(betasdEst + betasdGPSEst), length(latSeq))
        else
          sdSeqGPS = rep(exp(betasdGPSEst), length(latSeq))
      }
      else
        sdSeqGPS = rep(exp(betasdEst), length(latSeq))
    }
    plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq, sdSeqGPS)))
    sdSharedSeq = sdSeq * sqrt(omegaEst)
    lines(latSeq, sdSharedSeq, lty=2)
    lines(latSeq, sdSeqGPS, col="blue")
    lines(latSeq, sdSeqGPS * sqrt(omegaEst), col="blue", lty=2)
    
    # plot gamma and mean vectors
    latSeq = seq(latRange[1], latRange[2], l=500)
    gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
    if(!includeGammaSpline)
      gammaMat = matrix(1, nrow=length(latSeq))
    gammaSeq = exp(gammaMat %*% betaGammaEst)
    meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
    meanMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)
    if(!doMeanSpline) {
      meanMat = matrix(1, nrow=length(latSeq))
      meanMatGPS = matrix(1, nrow=length(latSeq))
    }
    meanSeq = exp(meanMat %*% betaMeanEst)
    if(diffMean)
      meanSeqGPS = exp(meanMat %*% betaMeanEst + as.numeric(diffMean) * (meanMatGPS %*% betaMeanGPSEst))
    else
      meanSeqGPS = exp(meanMat %*% betaMeanEst)
    plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
    plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0, meanSeqGPS)))
    lines(latSeq, meanSeqGPS, col="blue")
    
    par(mfrow=c(1,1))
    meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
    meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
    if(!doMeanSpline) {
      meanBasisX = matrix(1, nrow=length(threshSlipDat$lat))
      meanBasisXGPS = matrix(1, nrow=length(threshSlipDat$lat))
    }
    if(diffMean)
      meanVecX = exp(meanBasisX %*% betaMeanEst + as.numeric(diffMean) * (meanBasisXGPS %*% betaMeanGPSEst))
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
    
    report = obj$report()
    
    par(mfrow=c(1,2))
    plot(report$xResiduals, threshSlipDat$lat)
    plot(report$yResiduals, subDat$Lat)
    
    if(count == 1) {
      if(debugPlotting)
        browser()
    }
  }
  minI = which(objectiveValues == min(objectiveValues, na.rm = TRUE))[1]
  minPar = allPar[minI,]
  
  # opt$hessian ## <-- FD hessian from optim
  grad = obj$gr(minPar)
  if(calcHess) {
    hess = obj$he(minPar)    ## <-- Analytical hessian
    sdReport = sdreport(obj)
  }
  else {
    # calculate this stuff later
    hess = NULL
    sdReport = NULL
  }
  report = obj$report()
  
  print("optimum parameters for each iteration:")
  print(t(allPar[,1:count]))
  
  print("optimum parameters for best iteration:")
  print(minPar)
  
  print("optimum gradient for best iteration")
  print(grad)
  
  print("sdreport")
  print(report)
  
  optParNames = names(opt$par)
  logmuEst = minPar[which(optParNames == "logmu")]
  betaMeanEst = c(logmuEst, minPar[which(optParNames == "betaMean")])
  logMeanGPSEst = minPar[which(optParNames == "logMeanGPS")]
  betaMeanGPSEst = c(logMeanGPSEst, minPar[which(optParNames == "betaMeanGPS")])
  if(!diffMean)
    logMeanGPSEst = 0
  betasdInterceptEst = minPar[which(optParNames == "betasdIntercept")]
  betasdEst = c(betasdInterceptEst, minPar[which(optParNames == "betasd")])
  betasdInterceptGPSEst = minPar[which(optParNames == "betasdInterceptGPS")]
  betasdGPSEst = c(betasdInterceptGPSEst, minPar[which(optParNames == "betasdGPS")])
  betaTaperInterceptEst = minPar[which(optParNames == "betaTaperIntercept")]
  betaTaperEst = c(betaTaperInterceptEst, minPar[which(optParNames == "betaTaper")])
  betaTaperGPSEst = c(minPar[which(optParNames == "betaTaperGPS")])
  betaGammaInterceptEst = minPar[which(optParNames == "betaGammaIntercept")]
  betaGammaEst = c(betaGammaInterceptEst, minPar[which(optParNames == "betaGamma")])
  if(conditionalGamma) {
    ## use log-transformed OLS estimate
    # calculate the mean
    meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
    meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
    
    if(diffMean)
      meanVecX = exp(meanBasisX %*% betaMeanEst + as.numeric(diffMean) * (meanBasisXGPS %*% betaMeanGPSEst))
    else
      meanVecX = exp(meanBasisX %*% betaMeanEst)
    
    # calculate the taper
    lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
    lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
    
    if(!diffGPSTaper)
      taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst), dStar = dStarGPS)
    else
      taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst + lambdaBasisXGPS %*% betaTaperGPSEst), dStar = dStarGPS)
    
    # now perform the regression
    logX = log(data$x)
    newVarx = (data$xsd / data$x)^2
    offset = log(meanVecX) + log(taperVecX)
    gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
    if(!includeGammaSpline)
      gammaBasis = matrix(1, nrow=nrow(threshSlipDat))
    mod = lm(logX~gammaBasis - 1, offset=offset, weights=1 / newVarx)
    betaGammaEst = coef(mod)
  }
  logphiEst = minPar[which(optParNames == "logphi")]
  logalphaEst = minPar[which(optParNames == "logalpha")]
  loglowInflateEst = minPar[which(optParNames == "loglowInflate")]
  loghighInflateEst = minPar[which(optParNames == "loghighInflate")]
  if(!fixedPenalty && doSmoothnessPenalty) {
    betasdPenaltyLogLambdaEst = minPar[which(optParNames == "betasdPenaltyLogLambda")]
    if(!sharedPenalty) {
      betaTaperPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperPenaltyLogLambda")]
      betaTaperGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betaTaperGPSPenaltyLogLambda")]
      betaGammaPenaltyLogLambdaEst = minPar[which(optParNames == "betaGammaPenaltyLogLambda")]
      betaMeanPenaltyLogLambdaEst = minPar[which(optParNames == "betaMeanPenaltyLogLambdaEst")]
      betasdGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betasdGPSPenaltyLogLambda")]
      betaMeanGPSPenaltyLogLambdaEst = minPar[which(optParNames == "betaMeanGPSPenaltyLogLambda")]
    }
    else {
      betaTaperPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaTaperGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaGammaPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaMeanPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaMeanGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betasdGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
    }
  } else {
    betasdPenaltyLogLambdaEst = logPenaltyPar
    betaTaperPenaltyLogLambdaEst = logPenaltyPar
    betaTaperGPSPenaltyLogLambdaEst = logPenaltyPar
    betaGammaPenaltyLogLambdaEst = logPenaltyPar
    betaMeanPenaltyLogLambdaEst = logPenaltyPar
    betaMeanGPSPenaltyLogLambdaEst = logPenaltyPar
    betasdGPSPenaltyLogLambdaEst = logPenaltyPar
  }
  if(!is.null(varLogLambda)) {
    betasdPenaltyLogLambdaEst = varLogLambda
    betasdGPSPenaltyLogLambdaEst = varLogLambda
  }
  if(includeTaperDiffPar) {
    taperDiffPenaltyLogLambdaEst = minPar[which(optParNames == "taperDiffPenaltyLogLambdaEst")]
    meanDiffPenaltyLogLambdaEst = minPar[which(optParNames == "meanDiffPenaltyLogLambda")]
    sdDiffPenaltyLogLambdaEst = minPar[which(optParNames == "sdDiffPenaltyLogLambda")]
  }
  else if(fixedDiffPenalty) {
    taperDiffPenaltyLogLambdaEst = logDiffPenaltyPar
    meanDiffPenaltyLogLambdaEst = logDiffPenaltyPar
    sdDiffPenaltyLogLambdaEst = logDiffPenaltyPar
  }
  else {
    taperDiffPenaltyLogLambdaEst = NULL
    meanDiffPenaltyLogLambdaEst = NULL
    sdDiffPenaltyLogLambdaEst = NULL
  }
  if(sharedSpatialProcess) {
    logitOmegaEst = minPar[which(optParNames == "logitOmega")]
    omegaEst = expit(logitOmegaEst)
    
    if(estimateGpsShared) {
      logitOmega2Est = minPar[which(optParNames == "logitOmega2")]
      omega2Est = expit(logitOmega2Est)
    }
    else {
      logitOmega2Est = logitOmegaEst
      omega2Est = omegaEst
    }
  } else {
    omegaEst = 0
    logitOmegaEst = -Inf
    omega2Est = 0
    logitOmega2Est = -Inf
  }
  if(inflateVarLocking) {
    logLockInflateEst = minPar[which(optParNames == "logLockInflate")]
    lockInflateEst = exp(logLockInflateEst)
  } else {
    logLockInflateEst = log(3^2)
    lockInflateEst = exp(logLockInflateEst)
  }
  finalPar = c(logmuEst=logmuEst, betaMeanEst=betaMeanEst, logMeanGPSEst=logMeanGPSEst, betaMeanGPSEst=betaMeanGPSEst, 
               betasdEst=betasdEst, betasdGPSEst = betasdGPSEst, betaTaperEst=betaTaperEst, betaTaperGPSEst=betaTaperGPSEst, 
               betaGammaEst=betaGammaEst, logphiEst=logphiEst, logalphaEst=logalphaEst, 
               loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
               betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
               betasdGPSPenaltyLogLambdaEst=betasdGPSPenaltyLogLambdaEst, 
               betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
               betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
               betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
               betaMeanPenaltyLogLambdaEst=betaMeanPenaltyLogLambdaEst, 
               betaMeanGPSPenaltyLogLambdaEst=betaMeanGPSPenaltyLogLambdaEst, 
               taperDiffPenaltyLogLambdaEst=taperDiffPenaltyLogLambdaEst, 
               meanDiffPenaltyLogLambdaEst=meanDiffPenaltyLogLambdaEst, 
               sdDiffPenaltyLogLambdaEst = sdDiffPenaltyLogLambdaEst, 
               omegaEst=omegaEst, logitOmegaEst=logitOmegaEst, omega2Est=omega2Est, logitOmega2Est=logitOmega2Est, 
               logLockInflateEst=logLockInflateEst, 
               lockInflateEst=lockInflateEst)
  
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
  lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
  sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
  sdBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVarGPS, latRange=latRange)
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
  if(doVarSpline) {
    if(diffVar) {
      if(!reparameterizeVar)
        sdVecX = exp(sdBasisX %*% betasdEst + sdBasisXGPS %*% betasdGPSEst)
      else
        sdVecX = exp(sdBasisXGPS %*% betasdGPSEst)
    }
    else
      sdVecX = exp(sdBasisX %*% betasdEst)
    sdVecY = exp(sdBasisY %*% betasdEst)
  }
  else {
    if(diffVar) {
      if(!reparameterizeVar)
        sdVecX = rep(exp(betasdEst + betasdGPSEst), length(xDepths))
      else
        sdVecX = rep(exp(betasdGPSEst), length(xDepths))
    }
    else
      sdVecX = rep(exp(betasdEst ), length(xDepths))
    sdVecY = rep(exp(betasdEst), length(faultDepths))
  }
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
  sdMatGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVarGPS, lats=latSeq, latRange=latRange)
  if(doVarSpline) {
    sdSeq = exp(sdMat %*% betasdEst)
    if(diffVar) {
      if(!reparameterizeVar)
        sdSeqGPS = exp(sdMat %*% betasdEst + sdMatGPS %*% betasdGPSEst)
      else
        sdSeqGPS = exp(sdMatGPS %*% betasdGPSEst)
    }
    else
      sdSeqGPS = exp(sdMat %*% betasdEst)
  }
  else {
    sdSeq = rep(exp(betasdEst), length(latSeq))
    if(diffVar) {
      if(!reparameterizeVar)
        sdSeqGPS = rep(exp(betasdEst + betasdGPSEst), length(latSeq))
      else
        sdSeqGPS = rep(exp(betasdGPSEst), length(latSeq))
    }
    else
      sdSeqGPS = rep(exp(betasdEst), length(latSeq))
  }
  plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq, sdSeqGPS)))
  sdSharedSeq = sdSeq * sqrt(omegaEst)
  lines(latSeq, sdSharedSeq, lty=2)
  lines(latSeq, sdSeqGPS, col="blue")
  lines(latSeq, sdSeqGPS * sqrt(omegaEst), col="blue", lty=2)
  
  # plot gamma and mean vectors
  latSeq = seq(latRange[1], latRange[2], l=500)
  gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
  if(!includeGammaSpline)
    gammaMat = matrix(1, nrow=length(latSeq))
  gammaSeq = exp(gammaMat %*% betaGammaEst)
  meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
  meanMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)
  if(!doMeanSpline) {
    meanMat = matrix(1, nrow=length(latSeq))
    meanMatGPS = matrix(1, nrow=length(latSeq))
  }
  meanSeq = exp(meanMat %*% betaMeanEst)
  if(diffMean)
    meanSeqGPS = exp(meanMat %*% betaMeanEst + as.numeric(diffMean) * (meanMatGPS %*% betaMeanGPSEst))
  else
    meanSeqGPS = exp(meanMat %*% betaMeanEst)
  plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
  plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0, meanSeqGPS)))
  lines(latSeq, meanSeqGPS, col="blue")
  
  par(mfrow=c(1,1))
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
  meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
  if(!doMeanSpline) {
    meanBasisX = matrix(1, nrow=length(threshSlipDat$lat))
    meanBasisXGPS = matrix(1, nrow=length(threshSlipDat$lat))
  }
  if(diffMean)
    meanVecX = exp(meanBasisX %*% betaMeanEst + as.numeric(diffMean) * (meanBasisXGPS %*% betaMeanGPSEst))
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
  
  par(mfrow=c(1,2))
  plot(report$xResiduals, threshSlipDat$lat)
  plot(report$yResiduals, subDat$Lat)
  
  # Return results
  # return(list(Ests=Ests, muZetaEst=muZetaEst, sigmaZetaEst=sigmaZetaEst, lambdaEst=NA, 
  #             gammaEst=gammaEst, logLikEst=logLikEst, splineParEst=splinePar, phiEst=phiEst, 
  #             alphaEst=alphaEst, hess=hess, optimTable=optimTable, 
  #             tvec=tvec, tvecGPS=tvecGPS, optPar=opt$par, optGrad=optGrad, 
  #             strikeDistGps = strikeDistGps, dipDistGps = dipDistGps, 
  #             strikeDistCsz = strikeDistCsz, dipDistCsz = dipDistCsz))
  modelInfo = list(obj=obj, data=data, map=map, parameters=parameters, opt=opt, grad=grad, hess=hess, report=report, sdReport=sdReport, minPar=minPar, logmuEst=logmuEst, betaTaperEst=betaTaperEst, 
                   betaTaperGPSEst=betaTaperGPSEst, logphiEst=logphiEst, logalphaEst=logalphaEst, betasdEst=betasdEst, betasdEst=betasdEst, betasdGPSEst=betasdGPSEst, betaGammaEst=betaGammaEst, 
                   betaMeanEst=betaMeanEst, logMeanGPSEst=logMeanGPSEst, betaMeanGPSEst=betaMeanGPSEst, loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
                   strikeDistGps = sqrt(squareStrikeDistGps), dipDistGps = sqrt(squareDipDistGps),
                   strikeDistCsz = sqrt(squareStrikeDistCsz), dipDistCsz = sqrt(squareDipDistCsz), 
                   strikeDistCross = sqrt(squareStrikeDistCross), dipDistCross = sqrt(squareDipDistCross), 
                   finalPar=finalPar, betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
                   betasdGPSPenaltyLogLambdaEst=betasdGPSPenaltyLogLambdaEst, 
                   betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
                   betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
                   betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
                   taperDiffPenaltyLogLambdaEst=taperDiffPenaltyLogLambdaEst, 
                   omegaEst=omegaEst, omega2Est=omega2Est, diffGPSTaper=diffGPSTaper, 
                   taperVecX=taperVecX, taperVecY=taperVecY, sdVecX=sdVecX, sdVecY=sdVecY, gammaVec=gammaVec, 
                   gpsDepthThreshold=gpsDepthThreshold, maxCount=maxCount, diffLow=diffLow, diffHigh=diffHigh, 
                   penHigh=penHigh, penLow=penLow, transformation=transformation, diffMean=diffMean, 
                   diffVar=diffVar, allInputs=allInputs, allPar=allPar, 
                   doSmoothnessPenalty=doSmoothnessPenalty, doDiffPenalty=doDiffPenalty, 
                   estimateGpsShared=estimateGpsShared, varSmoothnessPenalty=varSmoothnessPenalty, 
                   varLogLambda=varLogLambda, reparameterizeVar=reparameterizeVar)
  
  if(saveResults) {
    shareText = ""
    if(sharedSpatialProcess && jointShared && estimateGpsShared) {
      shareText = "_estGpsShared"
    }
    else if(!estimateGpsShared) {
      shareText = "_jointShare"
    }
    else if(sharedSpatialProcess) {
      shareText = "_share"
    }
    fileName = paste0("fit", "_dS", dStar, "_dSGPS", dStarGPS, 
                      "_n", nKnots, 
                      "_nGPS", ifelse(diffGPSTaper, nKnotsGPS, 0), 
                      "_nGam", ifelse(includeGammaSpline, nKnotsGamma, 1), 
                      "_nVar", ifelse(doVarSpline, nKnotsVar, 1), 
                      "_nVarGPS", ifelse(diffVar, nKnotsVarGPS, 0), 
                      "_nMu", ifelse(doMeanSpline, nKnotsMean, 1), 
                      "_nMuGPS", ifelse(diffMean, nKnotsMeanGPS, 0), 
                      "_smoothP", doSmoothnessPenalty, 
                      "vSmoothP", varSmoothnessPenalty, 
                      "_fixPen", fixedPenalty, "_logPen", round(logPenaltyPar, 3), 
                      "_shareP", sharedPenalty, 
                      "_diffP", doDiffPenalty, 
                      "_logDiffP", round(logDiffPenaltyPar, 3), "_fixDiff", fixedDiffPenalty, 
                      "_varP", is.null(varLogLambda), 
                      "_repar", reparameterizeVar, 
                      "_diff", round(diffLow, 3), "-", round(diffHigh, 3), 
                      "_pen", round(round(diffLow, 3), 3), "-", round(penHigh, 3), 
                      "_hyp", useHyperpriors, 
                      shareText, 
                      "_cnstr", constrainMean, ".RData")
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

plotModelInfo = function(modelInfo, latRange=c(40, 50), fault=csz, gpsDat=slipDatCSZ, subDat=dr1) {
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
  diffMean = modelInfo$diffMean
  diffMean = modelInfo$diffMean
  diffVar = modelInfo$diffVar
  diffGPSTaper = modelInfo$data$diffGPSTaper
  betaMeanEst = modelInfo$betaMeanEst
  betaMeanGPSEst = modelInfo$betaMeanGPSEst
  # betaMeanGPSEst = modelInfo$betaMeanGPSEst
  # logMeanGPSEst = modelInfo$opt$par[which(parNames == "logMeanGPS")]
  if(!diffMean)
    logMeanGPSEst = 0
  betasdInterceptEst = modelInfo$betasdInterceptEst
  betasdEst = modelInfo$betasdEst
  betasdGPSEst = modelInfo$betasdGPSEst
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
  sharedSpatialProcess = modelInfo$sharedSpatialProcess
  jointShared = modelInfo$jointShared
  estimateGpsShared = modelInfo$estimateGpsShared
  omegaEst = modelInfo$omegaEst
  omega2Est = modelInfo$omega2Est
  
  report = modelInfo$report
  
  gpsDepthThreshold = modelInfo$gpsDepthThreshold
  dStar = modelInfo$data$dStar
  dStarGPS = modelInfo$data$dStarGPS
  
  nKnots = length(betaTaperEst)
  nKnotsGPS = length(betaTaperGPSEst)
  nKnotsGamma = length(betaGammaEst)
  nKnotsVar = length(betasdEst)
  nKnotsVarGPS = length(betasdGPSEst)
  nKnotsMean = length(betaMeanEst)
  nKnotsMeanGPS = length(betaMeanGPSEst)
  
  diffGPSTaper = length(betaTaperGPSEst) != 0
  includeGammaSpline = length(betaGammaEst) != 1
  doMeanSpline = length(betaMeanEst) != 1
  doVarSpline = length(modelInfo$betasdEst) != 1
  
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
  sdMatGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVarGPS, lats=latSeq, latRange=latRange)
  if(doVarSpline) {
    sdSeq = exp(sdMat %*% betasdEst)
    if(diffVar)
      sdSeqGPS = exp(sdMat %*% betasdEst + sdMatGPS %*% betasdGPSEst)
    else
      sdSeqGPS = exp(sdMat %*% betasdEst)
  }
  else {
    sdSeq = rep(exp(betasdEst), length(latSeq))
    if(diffVar)
      sdSeqGPS = rep(exp(betasdEst + betasdGPSEst), length(latSeq))
    else
      sdSeqGPS = rep(exp(betasdEst), length(latSeq))
  }
  plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq, sdSeqGPS)))
  sdSharedSeq = sdSeq * sqrt(omegaEst)
  lines(latSeq, sdSharedSeq, lty=2)
  lines(latSeq, sdSeqGPS, col="blue")
  lines(latSeq, sdSeqGPS * sqrt(omegaEst), col="blue", lty=2)
  
  # plot gamma and mean vectors
  latSeq = seq(latRange[1], latRange[2], l=500)
  gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
  if(!includeGammaSpline)
    gammaMat = matrix(1, nrow=length(latSeq))
  gammaSeq = exp(gammaMat %*% betaGammaEst)
  meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
  meanMatGPS = getSplineBasis(fault=fault, nKnots=nKnotsMeanGPS, lats=latSeq, latRange=latRange)
  if(!doMeanSpline) {
    meanMat = matrix(1, nrow=length(latSeq))
    meanMatGPS = matrix(1, nrow=length(latSeq))
  }
  meanSeq = exp(meanMat %*% betaMeanEst)
  if(diffMean)
    meanSeqGPS = exp(meanMat %*% betaMeanEst + as.numeric(diffMean) * (meanMatGPS %*% betaMeanGPSEst))
  else
    meanSeqGPS = exp(meanMat %*% betaMeanEst)
  plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
  plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0, meanSeqGPS)))
  lines(latSeq, meanSeqGPS, col="blue")
  
  par(mfrow=c(1,1))
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
  meanBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMeanGPS, latRange=latRange)
  if(!doMeanSpline) {
    meanBasisX = matrix(1, nrow=length(threshSlipDat$lat))
    meanBasisXGPS = matrix(1, nrow=length(threshSlipDat$lat))
  }
  if(diffMean)
    meanVecX = exp(meanBasisX %*% betaMeanEst + as.numeric(diffMean) * (meanBasisXGPS %*% betaMeanGPSEst))
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
  
  par(mfrow=c(1,2))
  plot(report$xResiduals, threshSlipDat$lat)
  plot(report$yResiduals, subDat$Lat)
  
}

# fullFit = fitModelTMB(fixedPenalty = TRUE, fixedDiffPenalty = TRUE, doDiffPenalty = TRUE, 
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

# return a updated list of parameters based on the optimization parameters
updateParameters = function(optPar, map, parameters) {
  optParNames = unique(names(optPar))
  newParameters = parameters
  for(i in 1:length(optParNames)) {
    thisName = optParNames[i]
    theseParameters = optPar[names(optPar) == thisName]
    parameters[[thisName]] = theseParameters
  }
  parameters
}

updateHessian = function(modelInfo, numerical=FALSE, ...) {
  if(!numerical)
    modelInfo$hess = modelInfo$obj$he()
  else
    modelInfo$hess = hessian(modelInfo$obj$fn, modelInfo$minPar, ...)
  modelInfo
}

updateSDReport = function(modelInfo, recompile=FALSE) {
  if(recompile) {
    if(!conditionalGamma) {
      if(!debug)
        compile("fitModelTMB.cpp")
      else
        compile("fitModelTMB.cpp","-O0 -g")
      dyn.load(dynlib("fitModelTMB"))
    }
    else {
      if(!debug)
        compile("fitModelTMBconditionalGamma.cpp")
      else
        compile("fitModelTMBconditionalGamma.cpp","-O0 -g")
      dyn.load(dynlib("fitModelTMBconditionalGamma"))
    }
  }
    
  # modelInfo$obj$fn(modelInfo$minPar)
  modelInfo$sdReport = sdreport(modelInfo$obj, modelInfo$minPar)
  rownames(modelInfo$sdReport) = names(modelInfo$opt$par)
  modelInfo
}
