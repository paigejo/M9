library(TMB)
library(matrixcalc)
library(MASS)
# source("setup.R")

# fits the combined model with normalized taper using TMB jointly with the variance inflation 
# parameters as well as the gamma spline parameters
fitModelTMB = function(initParams=NULL, gpsDat=slipDatCSZ, gpsDepthThreshold=21000, 
                       G=NULL, subDat=dr1, fault=csz, nKnots=5, 
                       dStar=25000, maxit=500, latRange=c(40, 50), 
                       diffGPSTaper=TRUE, nKnotsGPS=nKnots, reltol=1e-8, 
                       nKnotsVar=5, includeGammaSpline=FALSE, nKnotsGamma=7, 
                       dStarGPS=40000, seed=123, debug=FALSE, doVarSpline=TRUE, 
                       fixInflation=FALSE, maxCount=10, fixedPenalty=TRUE, 
                       logPenaltyPar=log(10), sharedPenalty=FALSE, doMeanSpline=FALSE, 
                       nKnotsMean=5) {
  
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
  betaTaper = out$taperPar
  betaTaperIntercept = betaTaper[1]
  betaTaper = betaTaper[-1]
  betaTaperGPS = out$taperParGPS
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
  
  # determine which of the input parameters will be optimized
  includeI = c(TRUE, TRUE, rep(doVarSpline, nKnotsVar - 1), rep(TRUE, nKnots), rep(TRUE, nKnotsGPS * diffGPSTaper), TRUE, rep(includeGammaSpline, nKnotsGamma - 1), 
               rep(!fixInflation, 2), TRUE, TRUE)
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
  cszStraight$centerX = newCenters[,1]
  cszStraight$centerY = newCenters[,2]
  straightenedGpsCoords = transformation(cbind(threshSlipDat$lon, threshSlipDat$lat))
  
  # calculate along strike and along dip squared distances in kilometers
  strikeCoords = cbind(0, cszStraight$centerY)
  dipCoords = cbind(cszStraight$centerX, 0)
  squareStrikeDistCsz = rdist(strikeCoords)^2
  squareDipDistCsz = rdist(dipCoords)^2
  
  # do the same for the gps data
  strikeCoords = cbind(0, straightenedGpsCoords[,2])
  dipCoords = cbind(straightenedGpsCoords[,1], 0)
  squareStrikeDistGps = rdist(strikeCoords)^2
  squareDipDistGps = rdist(dipCoords)^2
  
  # compute all required inputs for the TMB cpp code
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)[,-1]
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)[,-1]
  lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
  meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)[,-1]
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)[,-1]
  DSStrikeCSZ = squareStrikeDistCsz
  DSDipCSZ = squareDipDistCsz
  DSStrikeGPS = squareStrikeDistGps
  DSDipGPS = squareDipDistGps
  zeroMask = eventsEqMask(subDat)
  lowI = as.numeric(as.numeric(subDat$quality) != 1)
  y = -subDat$subsidence
  ysd = subDat$Uncertainty
  x = threshSlipDat$slip
  xsd = sqrt(log(.5*(sqrt(4*threshSlipDat$slipErr^2/threshSlipDat$slip^2 + 1) + 1)))
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
  sdBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnotsVar, lats=latSeq, latRange=latRange)[,-1]
  taperBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)[,-1]
  taperBasisGPSPenalty = getSplineBasis(fault=fault, nKnots=nKnotsGPS, lats=latSeq, latRange=latRange)
  gammaBasisPenalty = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)[,-1]
  
  # should choose penalty parameter (lambda) so that 
  # lambda * deltaPenalty * sum((diff(latSeq / 10) / deltaPenalty)^2) = lambda * 0.1
  # is 90% chance between 5, and 50. implying that the penalty factors should be between 50 and 500
  penaltyMean = mean(c(log(50), log(500)))
  penaltySD = (log(500)-penaltyMean)/qnorm(.95)
  
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
               sharedPenalty=as.numeric(sharedPenalty), meanBasisX=meanBasisX, meanBasisY=meanBasisY, 
               meanBasisPenalty=meanBasisPenalty)
  parameters = list(logmu=logmu, betasdIntercept=betasdIntercept, betasd=betasd, betaMean=betaMean, 
                    betaTaperIntercept=betaTaperIntercept, betaTaper=betaTaper, betaTaperGPS=betaTaperGPS, 
                    betaGammaIntercept=betaGammaIntercept, betaGamma=betaGamma, 
                    loglowInflate=loglowInflate, loghighInflate=loghighInflate, 
                    logphi=logphi, logalpha=logalpha, betasdPenaltyLogLambda=0, 
                    betaTaperPenaltyLogLambda=0, betaTaperGPSPenaltyLogLambda=0, 
                    betaGammaPenaltyLogLambda=0, betaMeanPenaltyLogLambda=0)
  if(fixedPenalty) {
    parameters$betasdPenaltyLogLambda = logPenaltyPar
    parameters$betaTaperPenaltyLogLambda = logPenaltyPar
    parameters$betaTaperGPSPenaltyLogLambda = logPenaltyPar
    parameters$betaGammaPenaltyLogLambda = logPenaltyPar
    betaMeanPenaltyLogLambda = logPenaltyPar
  }
  
  # pick which parameters will be fixed
  map = list()
  if(!includeGammaSpline)
    map = c(map, list(betaGamma=rep(factor(NA), length(betaGamma))))
  if(!doVarSpline)
    map = c(map, list(betasd=rep(factor(NA), length(betasd))))
  if(!doMeanSpline)
    map = c(map, list(betaMean=rep(factor(NA), length(betaMean))))
  if(fixInflation)
    map = c(map, list(loghighInflate=factor(NA), loglowInflate=factor(NA)))
  if(!diffGPSTaper)
    map = c(map, list(betaTaperGPS=rep(factor(NA), length(betaTaperGPS))))
  if(fixedPenalty) {
    map = c(map, list(betasdPenaltyLogLambda=factor(NA), 
                      betaTaperPenaltyLogLambda=factor(NA), 
                      betaTaperGPSPenaltyLogLambda=factor(NA),
                      betaGammaPenaltyLogLambda=factor(NA), 
                      betaMeanPenaltyLogLambda=factor(NA)))
  }
  else if(sharedPenalty) {
    map = c(map, list(betaTaperPenaltyLogLambda=factor(NA), 
                      betaTaperGPSPenaltyLogLambda=factor(NA),
                      betaGammaPenaltyLogLambda=factor(NA), 
                      betaMeanPenaltyLogLambda=factor(NA)))
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
  ui = rbind(c(rep(0, length(obj$par)-1 - 2 - nPenaltyPar * !fixedPenalty * (!fixInflation)), 1, rep(0, 2 * (!fixInflation) + nPenaltyPar * !fixedPenalty)), 
             c(rep(0, length(obj$par)-1 - 2 - nPenaltyPar * !fixedPenalty * (!fixInflation)), -1, rep(0, 2 * (!fixInflation) + nPenaltyPar * !fixedPenalty)))
  ci = c(log(0.05), -log(20))
  obj$ui = ui
  obj$ci = ci
  
  # reparameterize inputs from optim to constrOptim
  obj$f = obj$fn
  obj$grad = obj$gr
  obj$theta = obj$par
  
  # set optimization control parameters (adding in penalty parameters)
  controls = list(parscale=c(parscale, rep(1, nPenaltyPar * !fixedPenalty + nKnotsMean * doMeanSpline)))
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
    
    # opt <- do.call("optim", obj)
    opt <- do.call("constrOptim", obj)
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
    betasdInterceptEst = minPar[which(optParNames == "betasdIntercept")]
    betasdEst = c(betasdInterceptEst, minPar[which(optParNames == "betasd")])
    betaTaperInterceptEst = minPar[which(optParNames == "betaTaperIntercept")]
    betaTaperEst = c(betaTaperInterceptEst, minPar[which(optParNames == "betaTaper")])
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
      }
      else {
        betaTaperPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaTaperGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
        betaGammaPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      }
    } else {
      betasdPenaltyLogLambdaEst = logPenaltyPar
      betaTaperPenaltyLogLambdaEst = logPenaltyPar
      betaTaperGPSPenaltyLogLambdaEst = logPenaltyPar
      betaGammaPenaltyLogLambdaEst = logPenaltyPar
    }
    
    finalPar = c(logmuEst=logmuEst, betaMeanEst=betaMeanEst, betasdEst=betasdEst, betaTaperEst=betaTaperEst, 
                 betaGammaEst=betaGammaEst, logphiEst=logphiEst, logalphaEst=logalphaEst, 
                 loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
                 betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
                 betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
                 betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
                 betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
                 betaMeanPenaltyLogLambdaEst=betaMeanPenaltyLogLambdaEst)
    
    lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
    lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
    sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)
    sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
    meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
    meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)
    gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
    
    taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst), dStar = dStarGPS)
    taperVecY = taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar)
    sdVecX = exp(sdBasisX %*% betasdEst)
    sdVecY = exp(sdBasisY %*% betasdEst)
    meanVecX = exp(meanBasisX %*% betaMeanEst)
    meanVecY = exp(meanaBasisY %*% betaMeanEst)
    gammaVec = exp(gammaBasis %*% betaGammaEst)
    
    # preview results with plot
    # plot lambda and taper
    latSeq = seq(latRange[1], latRange[2], l=500)
    splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
    lambdaSeq = exp(splineMat %*% betaTaperEst)
    par(mfrow=c(1,2))
    plot(latSeq, lambdaSeq, type="l", ylim=range(c(lambdaSeq, 0)))
    tmp = getTaperSpline(betaTaperEst, nKnots=nKnots, dStar=dStar, latRange=latRange, fault=fault, normalize=TRUE, logScale = TRUE)
    plotFault(fault, tmp, varRange=c(0, 1))
    
    # plot taper (on original scale and GPS scale) and standard deviation
    latSeq = seq(latRange[1], latRange[2], l=500)
    splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
    lambdaSeq = exp(splineMat %*% betaTaperEst)
    taperSeq = taper(10000, lambdaSeq, dStar = dStar)
    taperSeqGPS = taper(10000, lambdaSeq, dStar = dStarGPS)
    par(mfrow=c(1,2))
    plot(latSeq, taperSeq, type="l", ylim=range(c(taperSeq, 0, taperSeqGPS)))
    lines(latSeq, taperSeqGPS, col="blue")
    legend("bottomright", c("Fault taper", "Locking rate taper"), col=c("blue", "black"), lty=1)
    sdMat = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, lats=latSeq, latRange=latRange)
    sdSeq = exp(sdMat %*% betasdEst)
    plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq)))
    
    # plot gamma and mean vectors
    latSeq = seq(latRange[1], latRange[2], l=500)
    gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
    gammaSeq = exp(gammaMat %*% betaGammaEst)
    meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
    meanSeq = exp(meanMat %*% betaMeanEst)
    plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
    plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0)))
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
  betasdInterceptEst = minPar[which(optParNames == "betasdIntercept")]
  betasdEst = c(betasdInterceptEst, minPar[which(optParNames == "betasd")])
  betaTaperInterceptEst = minPar[which(optParNames == "betaTaperIntercept")]
  betaTaperEst = c(betaTaperInterceptEst, minPar[which(optParNames == "betaTaper")])
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
    }
    else {
      betaTaperPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaTaperGPSPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
      betaGammaPenaltyLogLambdaEst = betasdPenaltyLogLambdaEst
    }
  } else {
    betasdPenaltyLogLambdaEst = logPenaltyPar
    betaTaperPenaltyLogLambdaEst = logPenaltyPar
    betaTaperGPSPenaltyLogLambdaEst = logPenaltyPar
    betaGammaPenaltyLogLambdaEst = logPenaltyPar
  }
  
  finalPar = c(logmuEst=logmuEst, betaMeanEst=betaMeanEst, betasdEst=betasdEst, betaTaperEst=betaTaperEst, 
               betaGammaEst=betaGammaEst, logphiEst=logphiEst, logalphaEst=logalphaEst, 
               loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
               betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
               betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
               betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
               betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
               betaMeanPenaltyLogLambdaEst=betaMeanPenaltyLogLambdaEst)
  
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)
  sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
  meanBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsMean, latRange=latRange)
  meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)
  gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
  
  taperVecX = taper(xDepths, exp(lambdaBasisX %*% betaTaperEst), dStar = dStarGPS)
  taperVecY = taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar)
  sdVecX = exp(sdBasisX %*% betasdEst)
  sdVecY = exp(sdBasisY %*% betasdEst)
  meanVecX = exp(meanBasisX %*% betaMeanEst)
  meanVecY = exp(meanaBasisY %*% betaMeanEst)
  gammaVec = exp(gammaBasis %*% betaGammaEst)
  
  # preview results with plot
  # plot lambda and taper
  latSeq = seq(latRange[1], latRange[2], l=500)
  splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
  lambdaSeq = exp(splineMat %*% betaTaperEst)
  par(mfrow=c(1,2))
  plot(latSeq, lambdaSeq, type="l", ylim=range(c(lambdaSeq, 0)))
  tmp = getTaperSpline(betaTaperEst, nKnots=nKnots, dStar=dStar, latRange=latRange, fault=fault, normalize=TRUE, logScale = TRUE)
  plotFault(fault, tmp, varRange=c(0, 1))
  
  # plot taper (on original scale and GPS scale) and standard deviation
  latSeq = seq(latRange[1], latRange[2], l=500)
  splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
  lambdaSeq = exp(splineMat %*% betaTaperEst)
  taperSeq = taper(10000, lambdaSeq, dStar = dStar)
  taperSeqGPS = taper(10000, lambdaSeq, dStar = dStarGPS)
  par(mfrow=c(1,2))
  plot(latSeq, taperSeq, type="l", ylim=range(c(taperSeq, 0, taperSeqGPS)))
  lines(latSeq, taperSeqGPS, col="blue")
  legend("bottomright", c("Fault taper", "Locking rate taper"), col=c("blue", "black"), lty=1)
  sdMat = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, lats=latSeq, latRange=latRange)
  sdSeq = exp(sdMat %*% betasdEst)
  plot(latSeq, sdSeq, type="l", ylim=range(c(0, sdSeq)))
  
  # plot gamma and mean vectors
  latSeq = seq(latRange[1], latRange[2], l=500)
  gammaMat = getSplineBasis(fault=fault, nKnots=nKnotsGamma, lats=latSeq, latRange=latRange)
  gammaSeq = exp(gammaMat %*% betaGammaEst)
  meanMat = getSplineBasis(fault=fault, nKnots=nKnotsMean, lats=latSeq, latRange=latRange)
  meanSeq = exp(meanMat %*% betaMeanEst)
  plot(latSeq, gammaSeq, type="l", ylim=range(c(gammaSeq, 0)))
  plot(latSeq, meanSeq, type="l", ylim=range(c(meanSeq, 0)))
  
  # Return results
  # return(list(Ests=Ests, muZetaEst=muZetaEst, sigmaZetaEst=sigmaZetaEst, lambdaEst=NA, 
  #             gammaEst=gammaEst, logLikEst=logLikEst, splineParEst=splinePar, phiEst=phiEst, 
  #             alphaEst=alphaEst, hess=hess, optimTable=optimTable, 
  #             tvec=tvec, tvecGPS=tvecGPS, optPar=opt$par, optGrad=optGrad, 
  #             strikeDistGps = strikeDistGps, dipDistGps = dipDistGps, 
  #             strikeDistCsz = strikeDistCsz, dipDistCsz = dipDistCsz))
  list(obj=obj, data=data, map=map, opt=opt, grad=grad, hess=hess, report=report, logmuEst=logmuEst, betaTaperEst=betaTaperEst, 
       logphiEst=logphiEst, logalphaEst=logalphaEst, betasdEst=betasdEst, betasdEst=betasdEst, betaGammaEst=betaGammaEst, 
       loglowInflateEst=loglowInflateEst, loghighInflateEst=loghighInflateEst, 
       strikeDistGps = sqrt(squareStrikeDistGps), dipDistGps = sqrt(squareDipDistGps),
       strikeDistCsz = sqrt(squareStrikeDistCsz), dipDistCsz = sqrt(squareDipDistCsz), 
       finalPar=finalPar, betasdPenaltyLogLambdaEst=betasdPenaltyLogLambdaEst, 
       betaTaperPenaltyLogLambdaEst=betaTaperPenaltyLogLambdaEst, 
       betaTaperGPSPenaltyLogLambdaEst=betaTaperGPSPenaltyLogLambdaEst, 
       betaGammaPenaltyLogLambdaEst=betaGammaPenaltyLogLambdaEst, 
       taperVecX=taperVecX, taperVecY=taperVecY, sdVecX=sdVecX, sdVecY=sdVecY, gammaVec=gammaVec)
}

getInitialParameters = function(nKnots=5, nKnotsVar=5, nKnotsGamma=7, logScale=FALSE, diffGPSTaper=FALSE, 
                                nKnotsGPS=nKnots) {
  # in order: 
  # mu, five SD parameters, five taper parameters, seven gamma parameters, 
  # low inflation, high inflation, spatial range, alpha/anisotropy parameter
  if(!diffGPSTaper)
    out = c(20, 15, rep(0, nKnotsVar - 1), 1, rep(0, nKnots - 1), 1, rep(0, nKnotsGamma - 1), 1.75, 1.25, 175, 1)
  else
    out = c(20, 15, rep(0, nKnotsVar - 1), 1, rep(0, nKnotsGPS + nKnots - 1), 1, rep(0, nKnotsGamma - 1), 1.75, 1.25, 175, 1)
    
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

test = fitModelTMB(debug=TRUE, doMeanSpline=TRUE)







