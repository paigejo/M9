# for fitting the model assuming the locking rate data product is tapered

fitModel2 = function(initParams=NULL, useMVNApprox=FALSE, gpsDat=slipDatCSZ, 
                     G=NULL, fauxG=NULL, subDat=dr1, fault=csz, nKnots=5, normalizeTaper=TRUE, 
                     dStar=21000, useGrad=FALSE, maxit=500, latRange=c(40, 50), 
                     normalModel=FALSE, doHess=TRUE, corGPS=FALSE, finalFit=FALSE, 
                     diffGPSTaper=FALSE, nKnotsGPS=nKnots, reltol=1e-8, anisotropic=FALSE) {
  
  ##### Set up initial parameter guesses for the MLE optimization if necessary
  # variance of a lognormal (x) is (exp(sigma^2) - 1) * exp(2mu + sigma^2)
  # then sigma^2 = log(var(e^x)) - 2mu - log(e^(sigma^2) - 1)
  # so sigma^2 \approx log(var(e^x)) - 2mu
  # sigmaZetaInit = sqrt(log(var(dr1$subsidence)) - 2*muZetaInit)
  if(is.null(initParams)) {
    lambdaInit=2
    if(!normalModel) {
      muZetaInit = log(20)
      phiInit = 250
      sigmaZetaInit = 1
    }
    else {
      muZetaInit = 20
      phiInit = 175
      sigmaZetaInit = sd(gpsDat$slip)*(muZetaInit/mean(gpsDat$slip))
    }
      
    # muXiInit = log(2.25) # based on plots from exploratory analysis
    # initSplineParams = getInitialSplineEsts(muZetaInit, sigmaZetaInit, lambdaInit, nKnots=nKnots, 
    # normalizeTaper=normalizeTaper, dStar=dStar)
    initSplineParams = getSplineEstsMomentMatch(muZetaInit, sigmaZetaInit, lambdaInit, corMatCSZ, nKnotsMean=nKnots, 
                                                nKnotsVar=nKnots, normalizeTaper=normalizeTaper, dStar=dStar, 
                                                normalModel=normalModel, useGrad=useGrad, fault=fault, 
                                                latRange=latRange, G=G, subDat=subDat)$betaHat
    initParams = c(muZetaInit, sigmaZetaInit, initSplineParams)
    
    # initialize with isotropic model even if we want to fit an anisotropic model eventually
    if(anisotropic) {
      initParams = c(muZetaInit, sigmaZetaInit, initSplineParams, 1)
    }
    
    initParams
  }
  
  # set spatial correlation smoothness
  nuZeta = 3/2
  
  ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
  # NOTE: only the upper triangle is computed for faster computation.  This is 
  #       sufficient for R's Cholesky decomposition
  # since we only want the correlation matrix, set sigmaZeta to 1
  # coords = cbind(csz$longitude, csz$latitude)
  # corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
  #                            onlyUpper=TRUE, Distance="rdist.earth",
  #                            Dist.args=list(miles=FALSE), smoothness=nuZeta)
  # don't compute the areally averaged covariance matrix: takes too long
  # tmpParams = rep(1, 5)
  # corMatCSZ = arealZetaCov(tmpParams, csz, nDown1=9, nStrike1=12)
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  if(is.null(fauxG))
    fauxG = getFauxG()
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(fault)[,3]
  
  ##### compute distance matrices for fault and for gps data
  coordsGPS = cbind(gpsDat$lon, gpsDat$lat)
  coordsCSZ = cbind(fault$longitude, fault$latitude)
  if(!anisotropic) {
    distMatGPS = rdist.earth(coordsGPS, miles=FALSE)
    distMatCSZ = rdist.earth(coordsCSZ, miles=FALSE)
    squareStrikeDistCsz = NULL
    squareDipDistCsz = NULL
    squareStrikeDistGps = NULL
    squareDipDistGps = NULL
  } else {
    distMatGPS = NULL
    distMatCSZ = NULL
    
    # the below code has bin commented out because the resulting covariances are not positive semi-definite
    # ## under anisotropic model, we compute the along strike and along dip distances separately
    # # fault geometry distances
    # out = compute_subfault_distances(csz)
    # squareStrikeDistCsz = out$Dstrike^2
    # squareDipDistCsz = out$Ddip^2
    # 
    # # gps data distances
    # out = compute_dip_strike_distance_gps(slipDatCSZ, csz)
    # squareStrikeDistGps = out$Dstrike^2
    # squareDipDistGps = out$Ddip^2
    
    ## under anisotropic model, we compute the along strike and along dip distances separately
    # straighten fault and gps geometries
    faultGeomStraight = straightenFault3()
    cszStraight = divideFault2(faultGeomStraight)
    cszStraight$centerX = 0.5 * (cszStraight$topLeftX + cszStraight$topRightX)
    cszStraight$centerY = 0.5 * (cszStraight$topLeftY + cszStraight$bottomRightY)
    straightenedGpsCoords = calcStraightGpsCoords(faultGeomStraight, faultGeom, gpsDat)
    
    # convert from new, straightened coordinate system to kilometers
    topI = 13
    bottomI = 1
    originalDist = rdist.earth(faultGeom[c(topI, bottomI), 2:3], miles=FALSE)[1, 2]
    newDist = rdist(cbind(c(faultGeomStraight$topMiddleX[topI], faultGeomStraight$topMiddleX[bottomI]), 
                          c(faultGeomStraight$topMiddleY[topI], faultGeomStraight$topMiddleY[bottomI])))[1, 2]
    kmPerUnit = originalDist / newDist
    
    # calculate along strike and along dip squared distances in kilometers
    strikeCoords = cbind(0, cszStraight$centerY)
    dipCoords = cbind(cszStraight$centerX, 0)
    squareStrikeDistCsz = (rdist(strikeCoords) * kmPerUnit)^2
    squareDipDistCsz = (rdist(dipCoords) * kmPerUnit)^2
    
    # do the same for the gps data
    strikeCoords = cbind(0, straightenedGpsCoords[,2])
    dipCoords = cbind(straightenedGpsCoords[,1], 0)
    squareStrikeDistGps = (rdist(strikeCoords) * kmPerUnit)^2
    squareDipDistGps = (rdist(dipCoords) * kmPerUnit)^2
  }
  
  ##### Do optimization
  optimTable <<- NULL # table to save likelihood optimization steps
  if(is.null(reltol)) {
    if(normalizeTaper)
      reltol=1e-13
    else
      reltol=1e-8
  }
  if(normalizeTaper)
    controls = list(fnscale = -1, reltol=reltol, parscale=c(2, 2, rep(1, length(initParams)-3-anisotropic), 25), maxit=maxit)
  else
    controls = list(fnscale = -1, reltol=reltol, parscale=c(2, 2, rep(1/21000, length(initParams)-3-anisotropic), 25), maxit=maxit)
  if(anisotropic)
    controls$parscale = c(controls$parscale, 2)
  
  if(finalFit) {
    controls$parscale = controls$parscale/25
  }
  
  cszDepthsIn = cszDepths
  corMatGPSIn = NULL
  corMatCSZIn = NULL
  corMatCSZLIn = NULL
  gpsDatIn = gpsDat
  nsimIn = 500 # this doens't matter since it's not used
  dStarIn = dStar
  normalizeTaperIn = normalizeTaper
  useMVNApproxIn = TRUE
  useASLApproxIn = FALSE
  nKnotsIn = nKnots
  GIn = G
  subDatIn = subDat
  faultIn = fault
  useSubPriorIn = FALSE
  useSlipPriorIn = FALSE
  fauxGIn = fauxG
  constrLambdaIn = FALSE
  latRangeIn = latRange
  normalModelIn = normalModel
  taperedGPSDatIn = TRUE
  distMatGPSIn = distMatGPS
  distMatCSZIn = distMatCSZ
  corGPSIn = corGPS
  diffGPSTaperIn = diffGPSTaper
  nKnotsGPSIn = nKnotsGPS
  anisotropicIn = anisotropic
  squareStrikeDistCszIn = squareStrikeDistCsz
  squareDipDistCszIn = squareDipDistCsz
  squareStrikeDistGpsIn = squareStrikeDistGps
  squareDipDistGpsIn = squareDipDistGps
  if(!useGrad)
    opt = optim(initParams, fixedDataLogLik, control=controls, hessian=FALSE, 
                cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                gpsDat=gpsDatIn, nsim=nsimIn, nKnots=nKnotsIn, 
                useMVNApprox=useMVNApproxIn, G=GIn, subDat=subDatIn, fault=faultIn, 
                normalizeTaper=normalizeTaperIn, dStar=dStarIn, useASLApprox=useASLApproxIn, 
                useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
                taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
                corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
                anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
                squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
                squareDipDistGps=squareDipDistGpsIn)
  else {
    lats=seq(latRange[1], latRange[2], l=100)
    Xi = getSplineBasis(latRange=latRange, nKnots=nKnots, lats=lats)
    if(diffGPSTaper) {
      gpsTaperZeroMat = matrix(0, ncol=nKnotsGPS, nrow=100)
      uiIn = cbind(0,0,-Xi, gpsTaperZeroMat, 0)
    }
    else
      uiIn = cbind(0,0,-Xi, 0)
    if(anisotropic)
      uiIn = cbind(uiIn, 0)
    opt = constrOptim(initParams, fixedDataLogLik, fixedDataLogLikGrad, control=controls, hessian=FALSE, 
                cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                gpsDat=gpsDatIn, nKnots=nKnotsIn, 
                useMVNApprox=TRUE, G=GIn, subDat=subDatIn, 
                fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
                useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
                taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
                corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
                anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
                squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
                squareDipDistGps=squareDipDistGpsIn, 
                method="BFGS", ui=uiIn, ci=rep(-10, 100))
    # use BFGS! (constraint makes lambdas not be able to be bigger than 10.  Somehow that works...)
  }
  
  # test = fullDataLogLik(initParams, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, 
  # phiZeta=phiZeta, slipDatCSZ=slipDatCSZ)
  
  # get results of optimization
  muZetaMLE = opt$par[1]
  sigmaZetaMLE = opt$par[2]
  muZetaGPS = muZetaMLE
  muZetaCSZ = muZetaMLE
  if(length(opt$par) > 3 || (nKnots==1)) {
    splinePar = opt$par[(length(opt$par)-(nKnots-1)-1):(length(opt$par) - 1)]
  }
  if(!anisotropic) {
    phiMLE = opt$par[length(opt$par)]
    alphaMLE = NA
  }
  else {
    phiMLE = opt$par[length(opt$par) - 1]
    alphaMLE = opt$par[length(opt$par)]
  }
  logLikMLE = opt$value
  if(doHess) {
    hess.args = list(eps=1e-7, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)
    hess = hessian(fixedDataLogLik, opt$par, method.args=hess.args,
                   cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                   gpsDat=gpsDatIn, nKnots=nKnotsIn, 
                   useMVNApprox=TRUE, G=GIn, subDat=subDatIn, 
                   fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
                   useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                   constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
                   taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
                   corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
                   anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
                   squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
                   squareDipDistGps=squareDipDistGpsIn)
    # hess = opt$hessian
  }
  else {
    hess=NULL
  }
  
  # calculate gradient at optimum to make sure we're at a local maximum
  if(useGrad)
    optGrad = fixedDataLogLikGrad(opt$par, cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                                  gpsDat=gpsDatIn, nKnots=nKnotsIn, 
                                  useMVNApprox=TRUE, G=GIn, subDat=subDatIn, 
                                  fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
                                  useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                                  constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
                                  taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
                                  corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
                                  anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
                                  squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
                                  squareDipDistGps=squareDipDistGpsIn)
  else
    optGrad = jacobian(fixedDataLogLik, opt$par, cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                       gpsDat=gpsDatIn, nKnots=nKnotsIn, 
                       useMVNApprox=TRUE, G=GIn, subDat=subDatIn, 
                       fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
                       useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                       constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
                       taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
                       corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
                       anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
                       squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
                       squareDipDistGps=squareDipDistGpsIn)
  
  # compute spline basis matrix for subsidence data
  Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  XiGPS1 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  
  # now get lambdas
  taperPar = opt$par[3:(length(opt$par)-1)]
  taperParGPS = taperPar[-(1:nKnots)]
  taperPar = taperPar[1:nKnots]
  lambda = Xi %*% taperPar
  lambdaGPS = XiGPS1 %*% taperPar
  if(diffGPSTaper) {
    XiGPS2 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
    lambdaGPS = lambdaGPS - XiGPS2 %*% taperParGPS
  }
  
  # compute taper vector
  tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  tvecGPS = taper(gpsDat$Depth, lambda=lambdaGPS, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  x = log(gpsDat$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  gammaEst = exp(sum((x- log(muZetaMLE) - log(tvecGPS))*ci))
  
  # params is in order: lambda, muZeta, sigmaZeta, lambda0, gamma, splineParams, phi
  MLEs = c(NA, opt$par[1:2], 0.25, gammaEst, opt$par[-c(1:2)])
  
  # compute final anisotropy distances if necessary
  if(anisotropic) {
    strikeDistGps = (1 / alphaMLE) * sqrt(squareStrikeDistGps)
    dipDistGps = alphaMLE * sqrt(squareDipDistGps)
    strikeDistCsz = (1 / alphaMLE) * sqrt(squareStrikeDistCsz)
    dipDistCsz = alphaMLE * sqrt(squareDipDistCsz)
  } else {
    strikeDistGps = NULL
    dipDistGps = NULL
    strikeDistCsz = NULL
    dipDistCsz = NULL
  }
  
  # Return results
  return(list(MLEs=MLEs, muZetaMLE=muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambdaMLE=NA, 
              gammaEst=gammaEst, logLikMLE=logLikMLE, splineParMLE=splinePar, phiMLE=phiMLE, 
              alphaMLE=alphaMLE, hess=hess, optimTable=optimTable, 
              tvec=tvec, tvecGPS=tvecGPS, optPar=opt$par, optGrad=optGrad, 
              strikeDistGps = strikeDistGps, dipDistGps = dipDistGps, 
              strikeDistCsz = strikeDistCsz, dipDistCsz = dipDistCsz))
}


# compute derivative of Matern correlation matrix with respect to phi (assuming 3/2 smoothness)
corrGrad = function(distMat, phiZeta, corMat, 
                    anisotropic=FALSE, alpha=NULL, squareDistMatStrike=NULL, squareDistMatDip=NULL) {
  # distMat * (-1/phiZeta^2 + 1/(-phiZeta^2+(-phiZeta)*distMat)) * corMat
  if(!anisotropic)
    distMat^2/phiZeta^3 * exp(-distMat/phiZeta)
  else {
    ## an anisotropic case, compute derivative with respect to phi and alpha
    # derivative with respect to phi
    dphi = distMat^2/phiZeta^3 * exp(-distMat/phiZeta)
    
    # derivative with respect to alpha
    numerator = alpha * squareDistMatStrike - alpha^(-3) * squareDistMatDip
    denominator = distMat
    # dalpha = (1 - 1 / phiZeta) * corMat * (numerator / denominator)
    dalpha = (1 / (phiZeta + distMat) - 1 / phiZeta) * corMat * (numerator / denominator)
    diag(dalpha) = 0
    
    list(dphi=dphi, dalpha=dalpha)
  }
}

# compute the derivative of the expectation of muX
meanXGrad = function(muZeta, nPar, tvecGPS, tGrad, gpsDat=slipDatCSZ, normalModel=FALSE, 
                     taperedGPSDat=FALSE, anisotropic=FALSE) {
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  x = log(gpsDat$slip)
  
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  if(normalModel) {
    x = gpsDat$slip
    sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
    ci = 1/sigmaXi^2
    ci = ci/sum(ci)
    muXi = sum((log(x)- log(muZeta) - log(tvecGPS))*ci)
    gammaEst = exp(muXi)
  }
  
  # compute gradient of mu_X:
  ExpXGrad = array(0, dim=c(length(x), nPar))
  if(taperedGPSDat) {
    if(!normalModel) {
      ExpXGrad[,1] = 1
      if(taperedGPSDat)
        ExpXGrad[,3:(nPar-1 - anisotropic)] = 1/tGrad
    }
    else {
      # ExpXGrad[,3:(nPar-1)] = gammaEst * muZeta * tGrad + outer(c(tvecGPS * muZeta), c(gammaEst * matrix(ci*tvecGPS, nrow=1) %*% tGrad))
      ExpXGrad[,3:(nPar-1 - anisotropic)] = muZeta * gammaEst * (tvecGPS %*% matrix(-ci/tvecGPS, nrow=1) + diag(nrow=length(tvecGPS))) %*% tGrad
    }
  }
  else {
    # do nothing, since gradient is 0
  }
  
  ExpXGrad
}

# compute the derivative of the covariance matrix of X
covXGrad = function(muZeta, sigmaZeta, tvecGPS, tGrad, corMatGPS, nPar, distMatGPS=NULL, 
                    phiZeta=NULL, gpsDat=slipDatCSZ, normalModel=FALSE, taperedGPSDat=FALSE, 
                    corGPS=FALSE, anisotropic=FALSE, alpha=1, squareDistMatStrike=NULL, 
                    squareDistMatDip=NULL) {
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  x = log(gpsDat$slip)
  
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  if(normalModel) {
    x = gpsDat$slip
    sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
    ci = 1/sigmaXi^2
    ci = ci/sum(ci)
    muXi = sum((log(x)- log(muZeta) - log(tvecGPS))*ci)
    gammaEst = exp(muXi)
    sigmaXi = gpsDat$slipErr
  }
  
  # compute gradient of Sigma_X:
  SigmaGrad = array(0, dim = c(length(x), length(x), nPar))
  if(!normalModel) {
    SigmaGrad[,,2] = 2*sigmaZeta*corMatGPS
    if(taperedGPSDat) {
      SigmaGrad[,,nPar] = sigmaZeta^2 * corrGrad(distMatGPS, phiZeta, corMatGPS, 
                                                 anisotropic, alpha, squareDistMatStrike, squareDistMatDip)
    }
  }
  else {
    dGammaEstDMu = -gammaEst/muZeta
    if(!taperedGPSDat) {
      SigmaGrad[,,1] = 2*gammaEst*sigmaZeta^2*corMatGPS * dGammaEstDMu
      SigmaGrad[,,2] = 2*gammaEst^2*sigmaZeta*corMatGPS
    }
    else {
      SigmaGrad[,,1] = 2*gammaEst* sweep(sweep(sigmaZeta^2*corMatGPS, 1, tvecGPS, "*"), 2, tvecGPS, "*") * dGammaEstDMu
      SigmaGrad[,,2] = 2*gammaEst^2*sigmaZeta*sweep(sweep(corMatGPS, 1, tvecGPS, "*"), 2, tvecGPS, "*")
      
      # use the taper gradient to compute the rest of the covariance gradient
      nKnotsTot = ncol(tGrad)
      for(i in 1:nKnotsTot) {
        pt1 = c(-sigmaZeta^2*2*gammaEst^2*(matrix(ci/tvecGPS, nrow=1)%*%matrix(tGrad[,nKnotsTot-i+1], ncol=1))) * sweep(sweep(corMatGPS, 1, tvecGPS, "*"), 2, tvecGPS, "*")
        # tmp = sweep(sweep(sigmaZeta^2 * corMatGPS, 1, tGrad[,nKnotsTot-i+1], "*"), 2, tvecGPS, "*")
        # pt2 = gammaEst^2 * (tmp + t(tmp))
        tmp = tGrad[,nKnotsTot-i+1] %*% matrix(tvecGPS, nrow=1)
        pt2 = gammaEst^2 * corMatGPS * sigmaZeta^2 * (tmp + t(tmp))
        SigmaGrad[,,nPar-i - anisotropic] = pt1 + pt2
      }
      
      if(!anisotropic) {
        CGrad = corrGrad(distMatGPS, phiZeta, corMatGPS, 
                         anisotropic, alpha, squareDistMatStrike, squareDistMatDip)
        SigmaGrad[,,nPar] = sigmaZeta^2 * gammaEst^2 * sweep(sweep(CGrad, 1, tvecGPS, "*"), 2, tvecGPS, "*")
        
        if(corGPS) {
          extraBit = sweep(sweep(CGrad, 1, sigmaXi, "*"), 2, sigmaXi, "*")
          SigmaGrad[,,nPar] = SigmaGrad[,,nPar] + extraBit
        }
      }
      else {
        out = corrGrad(distMatGPS, phiZeta, corMatGPS, 
                       anisotropic, alpha, squareDistMatStrike, squareDistMatDip)
        CGrad = out$dphi
        alphaGrad = out$dalpha
        
        SigmaGrad[,,nPar - 1] = sigmaZeta^2 * gammaEst^2 * sweep(sweep(CGrad, 1, tvecGPS, "*"), 2, tvecGPS, "*")
        SigmaGrad[,,nPar] = sigmaZeta^2 * gammaEst^2 * sweep(sweep(alphaGrad, 1, tvecGPS, "*"), 2, tvecGPS, "*")
        if(corGPS) {
          # if we assume the same correlation structure in the measurement error of xi, we must account for this
          extraBitPhi = sweep(sweep(CGrad, 1, sigmaXi, "*"), 2, sigmaXi, "*")
          extraBitAlpha = sweep(sweep(alphaGrad, 1, sigmaXi, "*"), 2, sigmaXi, "*")
          SigmaGrad[,,nPar - 1] = SigmaGrad[,,nPar - 1] + extraBitPhi
          SigmaGrad[,,nPar] = SigmaGrad[,,nPar] + extraBitAlpha
        }
      }
    }
  }
  
  SigmaGrad
}

# compute gamma estimate from model fit assuming a normal model
computeGammaEst = function(modelFit, gpsDat, tvecGPS=modelFit$tvecGPS) {
  muZeta = modelFit$muZetaMLE
  
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  x = gpsDat$slip
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  muXi = sum((log(x)- log(muZeta) - log(tvecGPS))*ci)
  gammaEst = exp(muXi)
  
  return(gammaEst)
}