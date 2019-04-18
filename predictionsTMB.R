# this script contains functions for generating predictions using TMB parameterizations

# function for computing predictive distribution given subsidence data.  Note that 
# the subsidence data should only be the data from one earthquake.  The returned 
# values relating to beta correspond to log zeta.  For instance, betaEsts is the 
# vector of estimates of log zeta for the given earthquake.  If prior=TRUE, then 
# Goldfinger 2012 data is used to create a prior for earthquake magnitude.
# NOTE: this function should only be called on data for a single fixed earthquake. 
# modelInfo: an object output by fitModelTMB
# subDat: uncertainties are automatically scaled by the parameters estimated in the modelInfo object
# pts: a set of longitude and latitude coordinates for which we want the predictions. 
#      these points will be transformed using the projection used for the analysis
# event: the earthquake for which we want predictions
predsGivenSubsidenceTMB = function(modelInfo, fault=csz, subDat=dr1, niter=500, gpsDat, 
                                   G=NULL, posNormalModel=FALSE, 
                                   pts=cbind(gpsDat$lon, gpsDat$lat), fastPNSim=TRUE, 
                                   event="T1", xDepths=gpsDat$Depth) {
  
  finalPar = modelInfo$finalPar
  
  # set other relevant parameters
  phiZeta = exp(modelInfo$logphiEst)
  alpha = exp(modelInfo$logalphaEst)
  nuZeta = 3/2
  dStarGPS = modelInfo$data$dStarGPS
  dStar = modelInfo$data$dStar
  
  highInflate = modelInfo$loghighInflateEst
  lowInflate = modelInfo$loglowInflateEst
  
  optParNames = names(modelInfo$opt$par)
  minPar = modelInfo$opt$par
  betaTaperEst = modelInfo$betaTaperEst
  betaTaperGPSEst = minPar[which(optParNames == "betaTaperGPS")]
  betasdEst = modelInfo$betasdEst
  betaMeanEst = modelInfo$betaMeanEst
  betaGammaEst = modelInfo$betaGammaEst
  betaGammaGPSEst = modelInfo$betaGammaGPSEst
  
  nKnots = length(modelInfo$betaTaperEst)
  nKnotsGPS = ncol(modelInfo$data$lambdaBasisXGPS)
  nKnotsVar = length(modelInfo$betasdEst)
  nKnotsGamma = length(modelInfo$betaGammaEst)
  # muZeta = exp(modelInfo$logmu)
  nKnotsMean = length(modelInfo$betaMeanEst)
  nKnotsMeanGPS = length(modelInfo$betaMeanEstGPS)
  
  diffGPSTaper = length(modelInfo$betaTaperGPSEst) != 0
  diffMean = length(betaGammaGPSEst) != 0
  doMeanSpline = modelInfo$allInputs$doMeanSpline
  includeGammaSpline = modelInfo$allInputs$includeGammaSpline
  
  faultDepths = getFaultCenters(fault)[,3]
  xDepths = gpsDat$Depth
  
  # generate spline basis matrix
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnots, latRange=latRange)
  lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsGPS, latRange=latRange)
  if(doMeanSpline) {
    meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)
    meanBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsMean, latRange=latRange)
    if(diffMean)
      meanBasisXGPS = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsMeanGPS, latRange=latRange)
  }
  else {
    meanBasisY = matrix(1, nrow=length(faultDepths), ncol=1)
    meanBasisX = matrix(1, nrow=length(xDepths), ncol=1)
    if(diffMean)
      meanBasisXGPS = matrix(1, nrow=length(xDepths), ncol=1)
  }
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
  sdBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsVar, latRange=latRange)
  if(includeGammaSpline) {
    gammaBasis = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsGamma, latRange=latRange)
  }
  else {
    gammaBasis = matrix(1, nrow=length(xDepths), ncol=1)
  }
  
  # evaluate splines on the fault and for the points of interest
  if(!diffGPSTaper)
    taperVecX = c(taper(xDepths, exp(lambdaBasisX %*% betaTaperEst), dStar = dStarGPS))
  else
    taperVecX = c(taper(xDepths, exp(lambdaBasisX %*% betaTaperEst + lambdaBasisXGPS %*% betaTaperGPSEst), dStar = dStarGPS))
  taperVecY = c(taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar))
  sdVecX = exp(sdBasisX %*% betasdEst)
  sdVecY = exp(sdBasisY %*% betasdEst)
  if(!diffMean)
    meanVecX = exp(meanBasisX %*% betaMeanEst)
  else
    meanVecX = exp(meanBasisX %*% betaMeanEst + meanBasisXGPS %*% betaMeanGPSEst)
  meanVecY = exp(meanBasisY %*% betaMeanEst)
  gammaVec = exp(gammaBasis %*% betaGammaEst)
  
  ## first transform data so it is uncorrelated standard normal
  n = nrow(fault)
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # get the Okada matrix and subsidence data for this specific event of interest
  isT1 = events==event
  T1Dat = subDat[isT1,]
  GT1 = G[isT1, ]
  
  # inflate uncertainties as necessary
  highQual = as.numeric(T1Dat$quality) == 1
  lowQual = as.numeric(T1Dat$quality) != 1
  T1Dat$Uncertainty[lowQual] = T1Dat$Uncertainty[lowQual]*lowInflate
  T1Dat$Uncertainty[highQual] = T1Dat$Uncertainty[highQual]*highInflate
  
  ### Rather than training the fault, we redefine an axis to be the strike access in Euclidean space
  ### using a Lambert projection and PCA
  out = straightenFaultLambert()
  faultGeomStraight = out$fault
  scale = out$scale
  parameters = out$projPar
  transformation = out$transformation
  straightPoints = transformation(pts)
  
  cszStraight = divideFault2(faultGeomStraight)
  centers = getFaultCenters(csz)[,1:2]
  newCenters = transformation(centers)
  cszStraight$centerX = newCenters[,2]
  cszStraight$centerY = newCenters[,1]
  
  # calculate along strike and along dip squared distances in kilometers
  strikeCoordsCSZ = cbind(0, cszStraight$centerX)
  dipCoordsCSZ = cbind(cszStraight$centerY, 0)
  squareStrikeDistCsz = rdist(strikeCoordsCSZ)^2
  squareDipDistCsz = rdist(dipCoordsCSZ)^2
  
  # do the same for the gps data
  strikeCoordsGps = cbind(0, straightPoints[,1])
  dipCoordsGps = cbind(straightPoints[,2], 0)
  squareStrikeDistGps = rdist(strikeCoordsGps)^2
  squareDipDistGps = rdist(dipCoordsGps)^2
  
  # compute the same for the cross distances between the point areal data
  squareDistMatCrossStrike = rdist(strikeCoordsGps, strikeCoordsCSZ)^2
  squareDistMatCrossDip = rdist(dipCoordsGps, dipCoordsCSZ)^2
  
  # compute point, fault, and cross distance matrices
  distMatGPS = sqrt(alpha^2 * squareStrikeDistGps + alpha^(-2) * squareDipDistGps)
  distMatCSZ = sqrt(alpha^2 * squareStrikeDistCsz + alpha^(-2) * squareDipDistCsz)
  distMatCross = sqrt(alpha^2 * squareDistMatCrossStrike + alpha^(-2) * squareDistMatCrossDip)
  
  # now compute the covariances
  # NOTE: the exact points passed to stationary.cov don't matter since the distance matrices are passed
  xs = cbind(fault$longitude, fault$latitude)
  SigmaFault = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                             smoothness=nuZeta, distMat = distMatCSZ)
  SigmaFault = sweep(sweep(SigmaFault, 2, sdVecY, "*"), 1, sdVecY, "*")
  
  # just for point wise predictions, shouldn't be used for locking rate predictions.
  SigmaPointToFault = stationary.cov(x1=pts, x2=strikeCoordsCSZ, Covariance="Matern", theta=phiZeta, 
                                    smoothness=nuZeta, distMat = distMatCross)
  SigmaPointToFault = sweep(sweep(SigmaPointToFault, 2, sdVecY, "*"), 1, sdVecX, "*")
  SigmaPoint = stationary.cov(pts, Covariance="Matern", theta=phiZeta, 
                              smoothness=nuZeta, distMat = distMatGPS)
  SigmaPoint = sweep(sweep(SigmaPoint, 2, sdVecX, "*"), 1, sdVecX, "*")
  
  # combine covariance matrices:
  SigmaZeta = rbind(cbind(SigmaPoint, SigmaPointToFault), 
                    cbind(t(SigmaPointToFault), SigmaFault))
  
  # set up required data variables in the R environment for Stan
  muVec = c(meanVecX, meanVecY)
  X = GT1 %*% diag(taperVecY)
  y = as.vector(-T1Dat$subsidence)
  sigmaY = T1Dat$Uncertainty
  n = length(y)
  pFault = nrow(fault)
  pPoint = nrow(pts)
  
  # need:
  # beta (logZeta areal), zeta (areal), logZetaPoint, seismicMoment, Mw
  # under the normal model, predictions becomes a conditional normal problem
  
  # get the distribution of Y and its covariance to zeta process TODO: include omega parameters here and in the conditional normal
  mvnApprox = estSubsidenceMeanCov(meanVecY, lambda, SigmaFault, GT1, subDat=T1Dat, 
                                   tvec=taperVecY, normalModel=TRUE)
  SigmaY = mvnApprox$Sigma
  diag(SigmaY) = diag(SigmaY) + sigmaY^2
  muY = X %*% muVec[-(1:nrow(pts))]
  SigmaZetaToY = rbind(SigmaPointToFault %*% t(X), 
                       SigmaFault %*% t(X))
  
  # compute the conditional distribution of S * zeta given Y
  # (Xd, muP, muD, SigmaP, SigmaD, SigmaPtoD)
  out = conditionalNormal(y, muVec, muY, SigmaZeta, SigmaY, SigmaZetaToY)
  muc = out$muc
  Sigmac = out$Sigmac
  
  # generate simulations of zeta from predictive distribution
  Lc = t(chol(Sigmac))
  
  if(posNormalModel) {
    notAllPos=TRUE
    nNewSims = niter*2
    zetaSims = matrix(-1, nrow=nrow(Lc), ncol=niter*2) # multiply by two for consistency with Stan MCMC results
    while(notAllPos) {
      # generate simulations until all slips are positive, if necessary
      # can check probability of generating all positive simulation with this code:
      # library(mvtnorm)
      # pmvnorm(upper=rep(0, nrow(csz)), mean=muc[-(1:nrow(gpsDat))], sigma=Sigmac[-(1:nrow(gpsDat)),-(1:nrow(gpsDat))])
      negCol = function(simCol) {
        any(simCol < 0)
      }
      negCols = apply(zetaSims, 2, negCol)
      if(nNewSims != sum(negCols)) {
        nNewSims = sum(negCols)
        print(paste0("number of simulations remaining: ", nNewSims))
      }
      
      if(! fastPNSim) {
        # simulate only for remaining columns
        zSims = matrix(rnorm(nNewSims*nrow(Lc)), nrow=nrow(Lc), ncol=nNewSims)
        thisZetaSims = sweep(Lc %*% zSims, 1, muc, "+") # add muZeta to each zero mean simulation
        zetaSims[,negCols] = thisZetaSims
      }
      else {
        # simulate a bunch and take any sims that are positive
        zSims = matrix(rnorm(niter*2*nrow(Lc)), nrow=nrow(Lc), ncol=niter*2)
        thisZetaSims = sweep(Lc %*% zSims, 1, muc, "+") # add muZeta to each zero mean simulation
        thisPosCols = which(!apply(thisZetaSims, 2, negCol))
        
        if(length(thisPosCols) > nNewSims) {
          zetaSims[,negCols] = thisZetaSims[,thisPosCols[1:nNewSims]]
        }
        else if(length(thisPosCols) > 0) {
          negColsI = which(negCols)
          zetaSims[,negColsI[1:length(thisPosCols)]] = thisZetaSims[,thisPosCols]
        }
        
      }
      
      notAllPos =  any(zetaSims < 0) && posNormalModel
    }
    logZetaSims = log(zetaSims)
  }
  else {
    zSims = matrix(rnorm(2*niter*nrow(Lc)), nrow=nrow(Lc), ncol=2*niter)
    zetaSims = sweep(Lc %*% zSims, 1, muc, "+") # add muZeta to each zero mean simulation
    logZetaSims = matrix(NA, ncol=2, nrow=nrow(zetaSims))
  }
  
  # seperate out point and block averaged zeta values
  zetaPoint = zetaSims[1:nrow(pts),]
  zetaAreal = zetaSims[-(1:nrow(pts)),]
  logZetaPoint = logZetaSims[1:nrow(pts),]
  logZetaAreal = logZetaSims[-(1:nrow(pts)),]
  
  # from zetaSims, generate seismic moments and magnitudes
  slipSims = sweep(zetaAreal, 1, taperVecY, "*")
  mags = apply(slipSims, 2, getMomentFromSlip, fault=fault)
  moments = 10^(mags*1.5 + 9.05)
  
  # compile results into predResults object with:
  # beta (logZeta areal), zeta (areal), logZetaPoint, seismicMoment, Mw (include zetaPoint as well since logZetaPoint is NA)
  predResults = list(beta=logZetaAreal, zeta=zetaAreal, logZetaPoint=logZetaPoint, zetaPoint=zetaPoint, 
                     seismicMoment=moments, Mw=mags)
  
  #   # show results
  #   print(predResults, digits = 1)
  
  # extract the samples of the relevant parameters
  tab = predResults
  betaTab = tab$beta
  zetaTab = tab$zeta
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = exp(logZetaPointTab)
  zetaPointTab = tab$zetaPoint
  M0Tab = tab$seismicMoment
  MwTab = tab$Mw
  
  # compute means, standard deviations, and middle 95% confidence intervals for relevant values
  betaEsts = NA
  betaSD = NA
  betaMed = NA
  beta025 = NA
  beta975 = NA
  zetaEsts = muc[-(1:nrow(pts))]
  zetaSD = sqrt(diag(Sigmac)[-(1:nrow(pts))])
  zetaMed = zetaEsts
  zeta025 = qnorm(.025, zetaEsts, zetaSD)
  zeta975 = qnorm(.975, zetaEsts, zetaSD)
  logZetaPointEsts = NA
  logZetaPointSD = NA
  logZetaPointMed = NA
  logZetaPoint025 = NA
  logZetaPoint975 = NA
  zetaPointEsts = muc[1:nrow(pts)]
  zetaPointSD = sqrt(diag(Sigmac)[1:nrow(pts)])
  zetaPointMed = zetaPointEsts
  zetaPoint025 = qnorm(.025, zetaPointEsts, zetaPointSD)
  zetaPoint975 = qnorm(.975, zetaPointEsts, zetaPointSD)
  M0Est = mean(M0Tab)
  M0SD = sd(M0Tab)
  M0Med = median(M0Tab)
  M0975 = quantile(probs=.975, M0Tab)
  M0025 = quantile(probs=.025, M0Tab)
  # must remove earthquakes with negative seismic moments, if that ever happens, since must take log
  posMags = M0Est > 0
  MwEst = mean(MwTab[posMags])
  MwSD = sd(MwTab[posMags])
  MwMed = median(MwTab[posMags])
  Mw975 = quantile(probs=.975, MwTab[posMags])
  Mw025 = quantile(probs=.025, MwTab[posMags])
  
  # estimate covariance matrix of params using emprical estimate from MCMC samples if necessary
  return(list(predResults=predResults, resultTab=tab, 
                betaEsts=betaEsts, betaSD=betaSD, betaMed=betaMed, beta025=beta025, beta975=beta975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zetaMed=zetaMed, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPointMed=logZetaPointMed, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPointMed=zetaPointMed, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaPointVarMat=NULL, zetaPointVarMat=NULL, 
                M0Est=M0Est, M0SD=M0SD, M0Med=M0Med, M0975=M0975, M0025=M0025, MwEst=MwEst, MwSD=MwSD, Mw975=Mw975, Mw025=Mw025))
}

# generate predictions given only the parameter MLEs (no GPS or subsidence data)
predsTMB = function(modelInfo, nsim=100, fault=csz, gpsDat, 
                    pts=cbind(gpsDat$lon, gpsDat$lat), 
                    posNormalModel=FALSE, fastPNSim=TRUE) {
  finalPar = modelInfo$finalPar
  
  # set other relevant parameters
  phiZeta = exp(modelInfo$logphiEst)
  alpha = exp(modelInfo$logalphaEst)
  nuZeta = 3/2
  dStarGPS = modelInfo$data$dStarGPS
  dStar = modelInfo$data$dStar
  
  highInflate = modelInfo$loghighInflateEst
  lowInflate = modelInfo$loglowInflateEst
  
  optParNames = names(modelInfo$opt$par)
  minPar = modelInfo$opt$par
  betaTaperEst = modelInfo$betaTaperEst
  # betaTaperGPSEst = minPar[which(optParNames == "betaTaperGPS")]
  betaTaperGPSEst = modelInfo$betaTaperGPSEst
  betasdEst = modelInfo$betasdEst
  betaMeanEst = modelInfo$betaMeanEst
  betaMeanGPSEst = modelInfo$betaMeanGPSEst
  betaGammaEst = modelInfo$betaGammaEst
  
  nKnots = length(modelInfo$betaTaperEst)
  nKnotsGPS = ncol(modelInfo$data$lambdaBasisXGPS)
  nKnotsVar = length(modelInfo$betasdEst)
  nKnotsGamma = length(modelInfo$betaGammaEst)
  # muZeta = exp(modelInfo$logmu)
  nKnotsMean = length(modelInfo$betaMeanEst)
  nKnotsMeanGPS = length(modelInfo$betaMeanGPSEst)
  
  diffGPSTaper = length(modelInfo$betaTaperGPSEst) != 0
  diffMean = length(modelInfo$betaMeanGPSEst) != 0
  
  # generate spline basis matrix
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnots, latRange=latRange)
  lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsGPS, latRange=latRange)
  meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)
  meanBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsMean, latRange=latRange)
  meanBasisXGPS = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsMeanGPS, latRange=latRange)
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
  sdBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsVar, latRange=latRange)
  gammaBasis = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsGamma, latRange=latRange)
  
  faultDepths = getFaultCenters(fault)[,3]
  xDepths = gpsDat$Depth
  
  # evaluate splines on the fault and for the points of interest
  taperVecY = c(taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar))
  sdVecX = exp(sdBasisX %*% betasdEst)
  sdVecY = exp(sdBasisY %*% betasdEst)
  if(!diffMean)
    meanVecX = exp(meanBasisX %*% betaMeanEst)
  else
    meanVecX = exp(meanBasisX %*% betaMeanEst + meanBasisXGPS %*% betaMeanGPSEst)
  meanVecY = exp(meanBasisY %*% betaMeanEst)
  gammaVec = exp(gammaBasis %*% betaGammaEst)
  
  # get CSZ prediction coordinates
  xp = cbind(fault$longitude, fault$latitude)
  
  ### Rather than training the fault, we redefine an axis to be the strike access in Euclidean space
  ### using a Lambert projection and PCA
  out = straightenFaultLambert()
  faultGeomStraight = out$fault
  scale = out$scale
  parameters = out$projPar
  transformation = out$transformation
  straightPoints = transformation(pts)
  
  cszStraight = divideFault2(faultGeomStraight)
  centers = getFaultCenters(csz)[,1:2]
  newCenters = transformation(centers)
  cszStraight$centerX = newCenters[,2]
  cszStraight$centerY = newCenters[,1]
  
  # calculate along strike and along dip squared distances in kilometers
  strikeCoordsCSZ = cbind(0, cszStraight$centerY)
  dipCoordsCSZ = cbind(cszStraight$centerX, 0)
  squareStrikeDistCsz = rdist(strikeCoordsCSZ)^2
  squareDipDistCsz = rdist(dipCoordsCSZ)^2
  
  # do the same for the gps data
  strikeCoordsGps = cbind(0, straightPoints[,1])
  dipCoordsGps = cbind(straightPoints[,2], 0)
  squareStrikeDistGps = rdist(strikeCoordsGps)^2
  squareDipDistGps = rdist(dipCoordsGps)^2
  
  # compute point, fault, and cross distance matrices
  distMatGPS = sqrt(alpha^2 * squareStrikeDistGps + alpha^(-2) * squareDipDistGps)
  distMatCSZ = sqrt(alpha^2 * squareStrikeDistCsz + alpha^(-2) * squareDipDistCsz)
  
  # now compute the covariances
  # NOTE: the exact points passed to stationary.cov don't matter since the distance matrices are passed
  SigmaD = stationary.cov(pts, Covariance="Matern", theta=phiZeta, 
                              smoothness=nuZeta, distMat = distMatGPS)
  SigmaD = sweep(sweep(SigmaD, 2, sdVecX, "*"), 1, sdVecX, "*")
  
  xs = cbind(fault$longitude, fault$latitude)
  Sigma = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                              smoothness=nuZeta, distMat = distMatCSZ)
  Sigma = sweep(sweep(Sigma, 2, sdVecY, "*"), 1, sdVecY, "*")
  SigmaL = t(chol(Sigma))
  
  # # generate predictive simulations
  # notAllPos=TRUE
  # zetaSims = matrix(-1, nrow=nrow(xp), ncol=nsim)
  # nNewSims = nsim
  # while(notAllPos) {
  #   # generate simulations until all slips are positive, if necessary
  #   negCol = function(simCol) {
  #     any(simCol < 0)
  #   }
  #   negCols = apply(zetaSims, 2, negCol)
  #   if(nNewSims != sum(negCols)) {
  #     nNewSims = sum(negCols)
  #     print(paste0("number of simulations remaining: ", nNewSims))
  #   }
  #   
  #   zSims = matrix(rnorm(nNewSims*nrow(xp)), nrow=nrow(xp), ncol=nNewSims)
  #   logZetaSims = sweep(SigmaL %*% zSims, 1, muZetaCSZ, "+") # add muZeta to each zero mean simulation
  #   if(!normalModel)
  #     zetaSims[,negCols] = exp(logZetaSims)
  #   else
  #     zetaSims[,negCols] = logZetaSims
  #   
  #   notAllPos =  any(zetaSims < 0) && posNormalModel
  # }
  # slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  muZetaCSZ = meanVecY
  muZetaGPS = meanVecX
  if(posNormalModel) {
    notAllPos=TRUE
    nNewSims = nsim
    zetaSims = matrix(-1, nrow=nrow(SigmaL), ncol=nsim) # multiply by two for consistency with Stan MCMC results
    while(notAllPos) {
      # generate simulations until all slips are positive, if necessary
      # can check probability of generating all positive simulation with this code:
      # library(mvtnorm)
      # pmvnorm(upper=rep(0, nrow(csz)), mean=muc[-(1:nrow(gpsDat))], sigma=Sigmac[-(1:nrow(gpsDat)),-(1:nrow(gpsDat))])
      negCol = function(simCol) {
        any(simCol < 0)
      }
      negCols = apply(zetaSims, 2, negCol)
      if(nNewSims != sum(negCols)) {
        nNewSims = sum(negCols)
        print(paste0("number of simulations remaining: ", nNewSims))
      }
      
      if(! fastPNSim) {
        # simulate only for remaining columns
        zSims = matrix(rnorm(nNewSims*nrow(SigmaL)), nrow=nrow(SigmaL), ncol=nNewSims)
        thisZetaSims = sweep(SigmaL %*% zSims, 1, muZetaCSZ, "+") # add muZeta to each zero mean simulation
        zetaSims[,negCols] = thisZetaSims
      }
      else {
        # simulate a bunch and take any sims that are positive
        zSims = matrix(rnorm(nsim*nrow(SigmaL)), nrow=nrow(SigmaL), ncol=nsim)
        thisZetaSims = sweep(SigmaL %*% zSims, 1, muZetaCSZ, "+") # add muZeta to each zero mean simulation
        thisPosCols = which(!apply(thisZetaSims, 2, negCol))
        
        if(length(thisPosCols) > nNewSims) {
          zetaSims[,negCols] = thisZetaSims[,thisPosCols[1:nNewSims]]
        }
        else if(length(thisPosCols) > 0) {
          negColsI = which(negCols)
          zetaSims[,negColsI[1:length(thisPosCols)]] = thisZetaSims[,thisPosCols]
        }
      }
      
      notAllPos =  any(zetaSims < 0) && posNormalModel
    }
    logZetaSims = log(zetaSims)
  }
  else {
    zSims = matrix(rnorm(nsim*nrow(SigmaL)), nrow=nrow(SigmaL), ncol=nsim)
    zetaSims = sweep(SigmaL %*% zSims, 1, muZetaCSZ, "+") # add muZeta to each zero mean simulation
    logZetaSims = matrix(NA, ncol=2, nrow=nrow(zetaSims))
  }
  slipSims = sweep(zetaSims, 1, taperVecY, FUN="*")
  
  # get mean slip prediction field
  if(!posNormalModel)
    meanSlip = muZetaCSZ * taperVecY
  else {
    meanSlip = apply(slipSims, 1, mean)
    if(nsim < 1000)
      warning("mean slip estimates may be poor with positive normal mode for <1000 simulations")
  }
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, Sigma=Sigma, Sigmac=Sigma, muc=muZetaCSZ, 
              SigmacGPS = SigmaD, mucGPS=muZetaGPS))
}

# Compute subsidence from the prediction simulations using the Okada model (NOTE: 
# returned ``subsidence'' is really uplift here).
# Preds is a list with elements named meanSlip (vector) and slipSims (matrix)
predsToSubsidenceTMB = function(modelInfo, preds, fault=csz, G=NULL, 
                                subDat=dr1, posNormalModel=FALSE) {
  finalPar = modelInfo$finalPar
  
  # set other relevant parameters
  phiZeta = exp(modelInfo$logphiEst)
  alpha = exp(modelInfo$logalphaEst)
  nuZeta = 3/2
  dStarGPS = modelInfo$data$dStarGPS
  dStar = modelInfo$data$dStar
  
  highInflate = modelInfo$loghighInflateEst
  lowInflate = modelInfo$loglowInflateEst
  
  optParNames = names(modelInfo$opt$par)
  minPar = modelInfo$opt$par
  betaTaperEst = modelInfo$betaTaperEst
  # betaTaperGPSEst = minPar[which(optParNames == "betaTaperGPS")]
  betaTaperGPSEst = modelInfo$betaTaperGPSEst
  betasdEst = modelInfo$betasdEst
  betaMeanEst = modelInfo$betaMeanEst
  betaMeanGPSEst = modelInfo$betaMeanGPSEst
  betaGammaEst = modelInfo$betaGammaEst
  
  nKnots = length(modelInfo$betaTaperEst)
  nKnotsGPS = ncol(modelInfo$data$lambdaBasisXGPS)
  nKnotsVar = length(modelInfo$betasdEst)
  nKnotsGamma = length(modelInfo$betaGammaEst)
  # muZeta = exp(modelInfo$logmu)
  nKnotsMean = length(modelInfo$betaMeanEst)
  nKnotsMeanGPS = length(modelInfo$betaMeanGPSEst)
  
  diffGPSTaper = length(modelInfo$betaTaperGPSEst) != 0
  diffMean = length(modelInfo$betaMeanGPSEst) != 0
  
  # generate spline basis matrix
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  
  faultDepths = getFaultCenters(fault)[,3]
  
  # evaluate splines on the fault and for the points of interest
  taperVecY = c(taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar))
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  tvec = taperVecY
  
  # transform slips into subsidences
  meanSlip = preds$meanSlip
  slipSims = preds$slipSims
  meanSub = G %*% cbind(meanSlip)
  subSims = G %*% slipSims
  
  # approximate upper and lower 95% quantiles (with either MVN approximate or simulations)
  # NOTE: use preds$Sigma not preds$Sigmac since Sigmac gives the covariance in mean estimate
  sigmaEps = subDat$Uncertainty
  if(!posNormalModel) {
    subMVN = estSubsidenceMeanCov(preds$muc, lambda, preds$Sigma, G, fault=fault, subDat=subDat, 
                                  normalModel=TRUE, tvec=tvec)
    subMu = subMVN$mu
    subSigma = subMVN$Sigma
    l95 = qnorm(.025, mean=subMu, sd=sqrt(diag(subSigma)))
    u95 = qnorm(.975, mean=subMu, sd=sqrt(diag(subSigma)))
    sigmaNoise = sqrt(diag(subSigma) + sigmaEps^2)
    l95Noise = qnorm(.025, mean=subMu, sd=sigmaNoise)
    u95Noise = qnorm(.975, mean=subMu, sd=sigmaNoise)
    subSimsNoise=NULL
  }
  else if(useMVNApprox) {
    subMVN = estSubsidenceMeanCov(preds$muc, lambda, preds$Sigma, G, subDat=subDat, fault=fault, 
                                  tvec=tvec)
    subMu = subMVN$mu
    subSigma = subMVN$Sigma
    l95 = qnorm(.025, mean=subMu, sd=sqrt(diag(subSigma)))
    u95 = qnorm(.975, mean=subMu, sd=sqrt(diag(subSigma)))
    sigmaNoise = sqrt(diag(subSigma) + sigmaEps^2)
    l95Noise = qnorm(.025, mean=subMu, sd=sigmaNoise)
    u95Noise = qnorm(.975, mean=subMu, sd=sigmaNoise)
    subSimsNoise=NULL
  }
  else {
    l95 = apply(subSims, 1, quantile, probs=.025)
    u95 = apply(subSims, 1, quantile, probs=.975)
    noiseSims = matrix(rnorm(length(subSims), 0, sigmaEps), nrow=nrow(subSims))
    subSimsNoise = subSims + noiseSims
    l95Noise = apply(subSimsNoise, 1, quantile, probs=.025)
    u95Noise = apply(subSimsNoise, 1, quantile, probs=.975)
  }
  # simulate middle 95% interval in observations
  
  return(list(meanSub = meanSub, subSims = subSims, l95=l95, u95=u95, l95Noise=l95Noise, 
              u95Noise=u95Noise, noiseSims=subSimsNoise))
}