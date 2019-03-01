# for fitting the model assuming the locking rate data product is tapered

fitModel2 = function(initParams=NULL, useMVNApprox=FALSE, gpsDat=slipDatCSZ, 
                     G=NULL, fauxG=NULL, subDat=dr1, fault=csz, nKnots=5, normalizeTaper=TRUE, 
                     dStar=21000, useGrad=FALSE, maxit=500, latRange=c(40, 50), 
                     normalModel=FALSE, doHess=TRUE, corGPS=FALSE, finalFit=FALSE, 
                     diffGPSTaper=FALSE, nKnotsGPS=nKnots, reltol=1e-8, anisotropic=FALSE, 
                     useGradForHess=(finalFit && doHess && useGrad), verbose=TRUE, 
                     nKnotsGamma=7, doGammaSpline=FALSE, nKnotsVar=5, doVarSpline=FALSE, 
                     dStarGPS=dStar) {
  
  ##### Set up initial parameter guesses for the MLE optimization if necessary (OUTDATED INITIALIZATION)
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
  
  # get input parameters
  out = getInputPar(initParams, fault, gpsDat, nKnots, diffGPSTaper, nKnotsGPS, taperedGPSDat, anisotropic, normalModel, nKnotsVar, doVarSpline)
  muZeta = out$muZeta
  muVecCSZ = out$muVecCSZ
  muVecGPS = out$muVecGPS
  sigmaZeta = out$sigmaZeta
  taperPar = out$taperPar
  taperParGPS = out$taperParGPS
  phiZeta = out$phiZeta
  nuZeta = out$nuZeta
  lambda0 = out$lambda0
  alpha = out$alpha
  varPar = out$varPar
  parscale = out$parscale
  parNames = out$parNames
  
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
    
    ### the below code is been commented out because the resulting transformation 
    ### modified the original distances too much
    # ## under anisotropic model, we compute the along strike and along dip distances separately
    # # straighten fault and gps geometries
    # faultGeomStraight = straightenFault3()
    # cszStraight = divideFault2(faultGeomStraight)
    # cszStraight$centerX = 0.5 * (cszStraight$topLeftX + cszStraight$topRightX)
    # cszStraight$centerY = 0.5 * (cszStraight$topLeftY + cszStraight$bottomRightY)
    # straightenedGpsCoords = calcStraightGpsCoords(faultGeomStraight, faultGeom, gpsDat)
    # 
    # # convert from new, straightened coordinate system to kilometers
    # topI = 13
    # bottomI = 1
    # originalDist = rdist.earth(faultGeom[c(topI, bottomI), 2:3], miles=FALSE)[1, 2]
    # newDist = rdist(cbind(c(faultGeomStraight$topMiddleX[topI], faultGeomStraight$topMiddleX[bottomI]),
    #                       c(faultGeomStraight$topMiddleY[topI], faultGeomStraight$topMiddleY[bottomI])))[1, 2]
    # kmPerUnit = originalDist / newDist
    # 
    # # calculate along strike and along dip squared distances in kilometers
    # strikeCoords = cbind(0, cszStraight$centerY * kmPerUnit)
    # dipCoords = cbind(cszStraight$centerX * kmPerUnit, 0)
    # squareStrikeDistCsz = rdist(strikeCoords)^2
    # squareDipDistCsz = rdist(dipCoords)^2
    # 
    # # do the same for the gps data
    # strikeCoords = cbind(0, straightenedGpsCoords[,2] * kmPerUnit)
    # dipCoords = cbind(straightenedGpsCoords[,1] * kmPerUnit, 0)
    # squareStrikeDistGps = rdist(strikeCoords)^2
    # squareDipDistGps = rdist(dipCoords)^2
    
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
    straightenedGpsCoords = transformation(cbind(gpsDat$lon, gpsDat$lat))

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
  }
  
  ##### Do optimization
  optimTable <<- NULL # table to save likelihood optimization steps
  if(is.null(reltol)) {
    if(normalizeTaper)
      reltol=1e-13
    else
      reltol=1e-8
  }
  # if(normalizeTaper)
  #   controls = list(fnscale = -1, reltol=reltol, parscale=c(2, 2, rep(1, length(initParams)-3-anisotropic), 25), maxit=maxit)
  # else
  #   controls = list(fnscale = -1, reltol=reltol, parscale=c(2, 2, rep(1/21000, length(initParams)-3-anisotropic), 25), maxit=maxit)
  # if(anisotropic)
  #   controls$parscale = c(controls$parscale, 2)
  # 
  # if(finalFit) {
  #   controls$parscale = controls$parscale/25
  # }
  controls = list(fnscale = -1, reltol=reltol, parscale=parscale, maxit=maxit)
  
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
  verboseIn = verbose
  doGammaSplineIn = doGammaSpline
  nKnotsGammaIn = nKnotsGamma
  dStarGPSIn=dStarGPS
  nKnotsVarIn=nKnotsVar
  doVarSplineIn=doVarSpline
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
                squareDipDistGps=squareDipDistGpsIn, verbose=verboseIn, doGammaSpline=doGammaSplineIn, 
                nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn, doVarSpline=doVarSplineIn, 
                nKnotsVar=nKnotsVarIn)
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
                squareDipDistGps=squareDipDistGpsIn, verbose=verboseIn, doGammaSpline=doGammaSplineIn, 
                nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn, doVarSpline=doVarSplineIn, 
                nKnotsVar=nKnotsVarIn, 
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
    print("Calculating Hessian...")
    # if(useGradForHess) {
    #   hess.args = list(eps=1e-7, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)
    #   hess = jacobian(fixedDataLogLikGrad, opt$par, method.args=hess.args,
    #                  cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
    #                  gpsDat=gpsDatIn, nKnots=nKnotsIn, 
    #                  useMVNApprox=TRUE, G=GIn, subDat=subDatIn, 
    #                  fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
    #                  useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
    #                  constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
    #                  taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
    #                  corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
    #                  anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
    #                  squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
    #                  squareDipDistGps=squareDipDistGpsIn, verbose=verboseIn)
    # } else {
    #   hess.args = list(eps=1e-7, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)
    #   hess = hessian(fixedDataLogLik, opt$par, method.args=hess.args,
    #                  cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
    #                  gpsDat=gpsDatIn, nKnots=nKnotsIn, 
    #                  useMVNApprox=TRUE, G=GIn, subDat=subDatIn, 
    #                  fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
    #                  useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
    #                  constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
    #                  taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
    #                  corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
    #                  anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
    #                  squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
    #                  squareDipDistGps=squareDipDistGpsIn, verbose=verboseIn)
    #   # hess = opt$hessian
    # }
      hess.args = list(eps=1e-4, d=0.001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE)
      hess = jacobian(fixedDataLogLikGrad, opt$par, method.args=hess.args,
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
                      squareDipDistGps=squareDipDistGpsIn, doGammaSpline=doGammaSplineIn, 
                      nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn, doVarSpline=doVarSplineIn, 
                      nKnotsVar=nKnotsVarIn, verbose=FALSE)
      # numeric hessian is numerically unstable
      # hessNumeric = hessian(fixedDataLogLik, opt$par, method.args=hess.args,
      #                cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
      #                gpsDat=gpsDatIn, nKnots=nKnotsIn, 
      #                useMVNApprox=TRUE, G=GIn, subDat=subDatIn, 
      #                fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
      #                useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
      #                constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
      #                taperedGPSDat=taperedGPSDatIn, distMatGPS=distMatGPSIn, distMatCSZ=distMatCSZIn, 
      #                corGPS=corGPSIn, diffGPSTaper=diffGPSTaperIn, nKnotsGPS=nKnotsGPSIn, 
      #                anisotropic=anisotropicIn, squareStrikeDistCsz=squareStrikeDistCszIn, 
      #                squareDipDistCsz=squareDipDistCszIn, squareStrikeDistGps=squareStrikeDistGpsIn, 
      #                squareDipDistGps=squareDipDistGpsIn, verbose=FALSE)
      # hess = opt$hessian
  }
  else {
    hess=NULL
  }
  
  # calculate gradient at optimum to make sure we're at a local maximum
  print("Calculating Jacobian")
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
                                  squareDipDistGps=squareDipDistGpsIn, doGammaSpline=doGammaSplineIn, 
                                  nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn, doVarSpline=doVarSplineIn, 
                                  nKnotsVar=nKnotsVarIn)
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
                       squareDipDistGps=squareDipDistGpsIn, doGammaSpline=doGammaSplineIn, 
                       nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn, doVarSpline=doVarSplineIn, 
                       nKnotsVar=nKnotsVarIn)
  
  # compute spline basis matrix for subsidence data
  Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  XiGPS1 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  
  # now get lambdas
  taperPar = opt$par[3:(length(opt$par)-1 - anisotropic)]
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
  tvecGPS = taper(gpsDat$Depth, lambda=lambdaGPS, alpha=2, normalize=normalizeTaper, dStar=dStarGPS)
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  x = log(gpsDat$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  
  if(doGammaSpline) {
    ys = x- log(muZetaMLE) - log(tvecGPS)
    X = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
    logModel = lm(ys~X-1, weights=ci)
    gammaEstGps = exp(X %*% logModel$coefficients)
    gammaCoefs = logModel$coefficients
    gammaBasis = X
    gammaEst = list(gammaEstGps=gammaEstGps, gammaCoefs=gammaCoefs, gammaBasis=gammaBasis)
    MLEs = c(NA, opt$par[1:2], 0.25, NA, opt$par[-c(1:2)])
    
    # plot gamma estimate
    sortI=sort(gpsDat$lat, index.return=TRUE)$ix
    plot(exp(ys), gpsDat$lat, pch=".", xlab=TeX("$\\gamma$"), ylab="Latitude", main=TeX("Estimate of $\\gamma$"))
    lines(gammaEstGps[sortI], gpsDat$lat[sortI], col="blue")
    
    quilt.plot(gpsDat$lon, gpsDat$lat, exp(ys), xlab="Longitude", ylab="Latitude", main=TeX("Estimate of $\\gamma$"))
    map("world", "Canada", add=TRUE)
    world(add=TRUE)
    
    # plot of years between earthquakes
    plot(1000/exp(ys), gpsDat$lat, pch=".", xlab="Time (Years)", ylab="Latitude", main="Estimate of Time Between Earthquakes")
    lines(1000/gammaEstGps[sortI], gpsDat$lat[sortI], col="blue")
    
    ticks = c(seq(0, 1000, by=100), seq(1500, max(pretty(1000/exp(ys), n=20)), by=500))
    quilt.plot(gpsDat$lon, gpsDat$lat, log10(1000/exp(ys)), xlab="Longitude", ylab="Latitude", main="Estimate of Years Between Earthquakes", 
               add.legend=FALSE)
    image.plot(add=TRUE, zlim=range(log10(1000/exp(ys))), nlevel=64, legend.only=TRUE, horizontal=FALSE, axis.args=list( at=log10(ticks), labels=ticks))
    map("world", "Canada", add=TRUE)
    world(add=TRUE)
  }
  else {
    gammaEst = exp(sum((x- log(muZetaMLE) - log(tvecGPS))*ci))
    # params is in order: lambda, muZeta, sigmaZeta, lambda0, gamma, splineParams, phi
    MLEs = c(NA, opt$par[1:2], 0.25, gammaEst, opt$par[-c(1:2)])
  }
  
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
                     taperedGPSDat=FALSE, anisotropic=FALSE, doGammaSpline=FALSE, nKnotsGamma=7) {
  
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
    
    if(doGammaSpline) {
      # first compute gamma estimates
      ys = log(x)- log(muZeta) - log(tvecGPS)
      X = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
      logModel = lm(ys~X-1, weights=ci)
      muXi = X %*% logModel$coefficients
      gammaEst = exp(muXi)
      
      # now compute the gradients of gamma estimates
      W = diag(ci)
      H = solve(t(X) %*% W %*% X, t(X) %*% W)
      tempMat = sweep(X %*% H, 1, -gammaEst, "*")
      dGammaMu = rowSums(tempMat) * (1 / muZeta)
      dGammaBeta = sweep(tempMat, 2, 1 / tvecGPS, "*") %*% tGrad
    } else {
      muXi = sum((log(x)- log(muZeta) - log(tvecGPS))*ci)
      gammaEst = exp(muXi)
    }
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
      if(doGammaSpline) {
        ExpXGrad[,3:(nPar-1 - anisotropic)] = sweep(tGrad, 1, gammaEst * muZeta, "*") + sweep(dGammaBeta, 1, tvecGPS * muZeta, "*")
      } else {
        # ExpXGrad[,3:(nPar-1)] = gammaEst * muZeta * tGrad + outer(c(tvecGPS * muZeta), c(gammaEst * matrix(ci*tvecGPS, nrow=1) %*% tGrad))
        ExpXGrad[,3:(nPar-1 - anisotropic)] = muZeta * gammaEst * (tvecGPS %*% matrix(-ci/tvecGPS, nrow=1) + diag(nrow=length(tvecGPS))) %*% tGrad
      }
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
                    squareDistMatDip=NULL, doGammaSpline=FALSE, nKnotsGamma=7) {
  
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
    if(doGammaSpline) {
      ys = log(x) - log(muZeta) - log(tvecGPS)
      X = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
      logModel = lm(ys~X-1, weights=ci)
      muXi = X %*% logModel$coefficients
      gammaEst = exp(muXi)
      
      # now compute the gradients of gamma estimates
      W = diag(ci)
      H = solve(t(X) %*% W %*% X, t(X) %*% W)
      tempMat = sweep(X %*% H, 1, -gammaEst, "*")
      dGammaMu = rowSums(tempMat) * (1 / muZeta)
      dGammaBeta = sweep(tempMat, 2, 1 / tvecGPS, "*") %*% tGrad
    } else {
      muXi = sum((log(x)- log(muZeta) - log(tvecGPS))*ci)
      gammaEst = exp(muXi)
    }
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
    if(!doGammaSpline)
      dGammaEstDMu = -gammaEst/muZeta
    
    if(!taperedGPSDat) {
      if(doGammaSpline)
        stop("Fitting gamma spline while assuming the gps data is not tapered is not yet supported")
      
      SigmaGrad[,,1] = 2*gammaEst*sigmaZeta^2*corMatGPS * dGammaEstDMu
      SigmaGrad[,,2] = 2*gammaEst^2*sigmaZeta*corMatGPS
    }
    else {
      if(!doGammaSpline) {
        SigmaGrad[,,1] = 2*gammaEst* sweep(sweep(sigmaZeta^2*corMatGPS, 1, tvecGPS, "*"), 2, tvecGPS, "*") * dGammaEstDMu
        SigmaGrad[,,2] = 2*gammaEst^2*sigmaZeta*sweep(sweep(corMatGPS, 1, tvecGPS, "*"), 2, tvecGPS, "*")
      }
      else {
        SigmaGrad[,,1] = sweep(sweep(sigmaZeta^2*corMatGPS, 1, tvecGPS * dGammaMu, "*"), 2, tvecGPS * gammaEst, "*")
        SigmaGrad[,,1] = SigmaGrad[,,1] + t(SigmaGrad[,,1])
        SigmaGrad[,,2] = 2*sigmaZeta*sweep(sweep(corMatGPS, 1, tvecGPS * gammaEst, "*"), 2, tvecGPS * gammaEst, "*")
      }
      
      # use the taper gradient to compute the rest of the covariance gradient
      nKnotsTot = ncol(tGrad)
      for(i in 1:nKnotsTot) {
        if(!doGammaSpline) {
          pt1 = c(-sigmaZeta^2*2*gammaEst^2*(matrix(ci/tvecGPS, nrow=1)%*%matrix(tGrad[,nKnotsTot-i+1], ncol=1))) * sweep(sweep(corMatGPS, 1, tvecGPS, "*"), 2, tvecGPS, "*")
          # tmp = sweep(sweep(sigmaZeta^2 * corMatGPS, 1, tGrad[,nKnotsTot-i+1], "*"), 2, tvecGPS, "*")
          # pt2 = gammaEst^2 * (tmp + t(tmp))
          tmp = tGrad[,nKnotsTot-i+1] %*% matrix(tvecGPS, nrow=1)
          pt2 = gammaEst^2 * corMatGPS * sigmaZeta^2 * (tmp + t(tmp))
          SigmaGrad[,,nPar-i - anisotropic] = pt1 + pt2
        }
        else {
          gradVec = dGammaBeta[,nKnotsTot-i+1] * tvecGPS + gammaEst * tGrad[,nKnotsTot-i+1]
          tmp = sweep(sweep(corMatGPS, 1, gradVec, "*"), 2, tvecGPS * gammaEst, "*")
          SigmaGrad[,,nPar-i - anisotropic] = sigmaZeta^2 * (tmp + t(tmp))
        }
      }
      
      if(!anisotropic) {
        CGrad = corrGrad(distMatGPS, phiZeta, corMatGPS, 
                         anisotropic, alpha, squareDistMatStrike, squareDistMatDip)
        
        if(!doGammaSpline)
          SigmaGrad[,,nPar] = sigmaZeta^2 * gammaEst^2 * sweep(sweep(CGrad, 1, tvecGPS, "*"), 2, tvecGPS, "*")
        else
          SigmaGrad[,,nPar] = sigmaZeta^2 * sweep(sweep(CGrad, 1, tvecGPS * gammaEst, "*"), 2, tvecGPS * gammaEst, "*")
        
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
        
        if(!doGammaSpline) {
          SigmaGrad[,,nPar - 1] = sigmaZeta^2 * gammaEst^2 * sweep(sweep(CGrad, 1, tvecGPS, "*"), 2, tvecGPS, "*")
          SigmaGrad[,,nPar] = sigmaZeta^2 * gammaEst^2 * sweep(sweep(alphaGrad, 1, tvecGPS, "*"), 2, tvecGPS, "*")
        }
        else {
          SigmaGrad[,,nPar - 1] = sigmaZeta^2 * sweep(sweep(CGrad, 1, tvecGPS * gammaEst, "*"), 2, tvecGPS * gammaEst, "*")
          SigmaGrad[,,nPar] = sigmaZeta^2 * sweep(sweep(alphaGrad, 1, tvecGPS * gammaEst, "*"), 2, tvecGPS * gammaEst, "*")
        }
        
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

# parametricBootstrap = function(optParams, useMVNApprox=FALSE, gpsDat=slipDatCSZ, 
#                                G=NULL, fauxG=NULL, subDat=dr1, fault=csz, nKnots=5, normalizeTaper=TRUE, 
#                                dStar=21000, useGrad=FALSE, maxit=500, latRange=c(40, 50), 
#                                normalModel=FALSE, doHess=TRUE, corGPS=FALSE, finalFit=FALSE, 
#                                diffGPSTaper=FALSE, nKnotsGPS=nKnots, reltol=1e-8, anisotropic=FALSE, 
#                                useGradForHess=(finalFit && doHess && useGrad), verbose=TRUE) {
#   
#   ##### Set up initial parameter guesses for the MLE optimization if necessary
#   # variance of a lognormal (x) is (exp(sigma^2) - 1) * exp(2mu + sigma^2)
#   # then sigma^2 = log(var(e^x)) - 2mu - log(e^(sigma^2) - 1)
#   # so sigma^2 \approx log(var(e^x)) - 2mu
#   # sigmaZetaInit = sqrt(log(var(dr1$subsidence)) - 2*muZetaInit)
#   if(is.null(optParams)) {
#     lambdaInit=2
#     if(!normalModel) {
#       muZetaInit = log(20)
#       phiInit = 250
#       sigmaZetaInit = 1
#     }
#     else {
#       muZetaInit = 20
#       phiInit = 175
#       sigmaZetaInit = sd(gpsDat$slip)*(muZetaInit/mean(gpsDat$slip))
#     }
#     
#     # muXiInit = log(2.25) # based on plots from exploratory analysis
#     # initSplineParams = getInitialSplineEsts(muZetaInit, sigmaZetaInit, lambdaInit, nKnots=nKnots, 
#     # normalizeTaper=normalizeTaper, dStar=dStar)
#     initSplineParams = getSplineEstsMomentMatch(muZetaInit, sigmaZetaInit, lambdaInit, corMatCSZ, nKnotsMean=nKnots, 
#                                                 nKnotsVar=nKnots, normalizeTaper=normalizeTaper, dStar=dStar, 
#                                                 normalModel=normalModel, useGrad=useGrad, fault=fault, 
#                                                 latRange=latRange, G=G, subDat=subDat)$betaHat
#     optParams = c(muZetaInit, sigmaZetaInit, initSplineParams)
#     
#     # initialize with isotropic model even if we want to fit an anisotropic model eventually
#     if(anisotropic) {
#       optParams = c(muZetaInit, sigmaZetaInit, initSplineParams, 1)
#     }
#     
#     optParams
#   }
#   
#   # set spatial correlation smoothness
#   nuZeta = 3/2
#   
#   ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
#   # NOTE: only the upper triangle is computed for faster computation.  This is 
#   #       sufficient for R's Cholesky decomposition
#   # since we only want the correlation matrix, set sigmaZeta to 1
#   # coords = cbind(csz$longitude, csz$latitude)
#   # corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
#   #                            onlyUpper=TRUE, Distance="rdist.earth",
#   #                            Dist.args=list(miles=FALSE), smoothness=nuZeta)
#   # don't compute the areally averaged covariance matrix: takes too long
#   # tmpParams = rep(1, 5)
#   # corMatCSZ = arealZetaCov(tmpParams, csz, nDown1=9, nStrike1=12)
#   
#   # get Okada linear transformation matrix
#   if(is.null(G)) {
#     nx = 300
#     ny=  900
#     lonGrid = seq(lonRange[1], lonRange[2], l=nx)
#     latGrid = seq(latRange[1], latRange[2], l=ny)
#     G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
#   }
#   if(is.null(fauxG))
#     fauxG = getFauxG()
#   
#   ##### calculate depths of the centers of the CSZ subfaults
#   cszDepths = getFaultCenters(fault)[,3]
#   
#   ##### compute distance matrices for fault and for gps data
#   coordsGPS = cbind(gpsDat$lon, gpsDat$lat)
#   coordsCSZ = cbind(fault$longitude, fault$latitude)
#   if(!anisotropic) {
#     distMatGPS = rdist.earth(coordsGPS, miles=FALSE)
#     distMatCSZ = rdist.earth(coordsCSZ, miles=FALSE)
#     squareStrikeDistCsz = NULL
#     squareDipDistCsz = NULL
#     squareStrikeDistGps = NULL
#     squareDipDistGps = NULL
#   } else {
#     distMatGPS = NULL
#     distMatCSZ = NULL
#     
#     # the below code has bin commented out because the resulting covariances are not positive semi-definite
#     # ## under anisotropic model, we compute the along strike and along dip distances separately
#     # # fault geometry distances
#     # out = compute_subfault_distances(csz)
#     # squareStrikeDistCsz = out$Dstrike^2
#     # squareDipDistCsz = out$Ddip^2
#     # 
#     # # gps data distances
#     # out = compute_dip_strike_distance_gps(slipDatCSZ, csz)
#     # squareStrikeDistGps = out$Dstrike^2
#     # squareDipDistGps = out$Ddip^2
#     
#     ### the below code is been commented out because the resulting transformation 
#     ### modified the original distances too much
#     # ## under anisotropic model, we compute the along strike and along dip distances separately
#     # # straighten fault and gps geometries
#     # faultGeomStraight = straightenFault3()
#     # cszStraight = divideFault2(faultGeomStraight)
#     # cszStraight$centerX = 0.5 * (cszStraight$topLeftX + cszStraight$topRightX)
#     # cszStraight$centerY = 0.5 * (cszStraight$topLeftY + cszStraight$bottomRightY)
#     # straightenedGpsCoords = calcStraightGpsCoords(faultGeomStraight, faultGeom, gpsDat)
#     # 
#     # # convert from new, straightened coordinate system to kilometers
#     # topI = 13
#     # bottomI = 1
#     # originalDist = rdist.earth(faultGeom[c(topI, bottomI), 2:3], miles=FALSE)[1, 2]
#     # newDist = rdist(cbind(c(faultGeomStraight$topMiddleX[topI], faultGeomStraight$topMiddleX[bottomI]),
#     #                       c(faultGeomStraight$topMiddleY[topI], faultGeomStraight$topMiddleY[bottomI])))[1, 2]
#     # kmPerUnit = originalDist / newDist
#     # 
#     # # calculate along strike and along dip squared distances in kilometers
#     # strikeCoords = cbind(0, cszStraight$centerY * kmPerUnit)
#     # dipCoords = cbind(cszStraight$centerX * kmPerUnit, 0)
#     # squareStrikeDistCsz = rdist(strikeCoords)^2
#     # squareDipDistCsz = rdist(dipCoords)^2
#     # 
#     # # do the same for the gps data
#     # strikeCoords = cbind(0, straightenedGpsCoords[,2] * kmPerUnit)
#     # dipCoords = cbind(straightenedGpsCoords[,1] * kmPerUnit, 0)
#     # squareStrikeDistGps = rdist(strikeCoords)^2
#     # squareDipDistGps = rdist(dipCoords)^2
#     
#     ### Rather than training the fault, we redefine an axis to be the strike access in Euclidean space
#     ### using a Lambert projection and PCA
#     out = straightenFaultLambert()
#     faultGeomStraight = out$fault
#     scale = out$scale
#     parameters = out$projPar
#     transformation = out$transformation
#     
#     cszStraight = divideFault2(faultGeomStraight)
#     centers = getFaultCenters(csz)[,1:2]
#     newCenters = transformation(centers)
#     cszStraight$centerX = newCenters[,1]
#     cszStraight$centerY = newCenters[,2]
#     straightenedGpsCoords = transformation(cbind(gpsDat$lon, gpsDat$lat))
#     
#     # calculate along strike and along dip squared distances in kilometers
#     strikeCoords = cbind(0, cszStraight$centerY)
#     dipCoords = cbind(cszStraight$centerX, 0)
#     squareStrikeDistCsz = rdist(strikeCoords)^2
#     squareDipDistCsz = rdist(dipCoords)^2
#     
#     # do the same for the gps data
#     strikeCoords = cbind(0, straightenedGpsCoords[,2] * alpha)
#     dipCoords = cbind(straightenedGpsCoords[,1] / alpha, 0)
#     squareStrikeDistGps = rdist(strikeCoords)^2
#     squareDipDistGps = rdist(dipCoords)^2
#   }
#   
#   # set inputs to the data simulation function
#   cszDepthsIn = cszDepths
#   corMatGPSIn = NULL
#   corMatCSZIn = NULL
#   corMatCSZLIn = NULL
#   gpsDatIn = gpsDat
#   nsimIn = 500 # this doens't matter since it's not used
#   dStarIn = dStar
#   normalizeTaperIn = normalizeTaper
#   useMVNApproxIn = TRUE
#   useASLApproxIn = FALSE
#   nKnotsIn = nKnots
#   GIn = G
#   subDatIn = subDat
#   faultIn = fault
#   useSubPriorIn = FALSE
#   useSlipPriorIn = FALSE
#   fauxGIn = fauxG
#   constrLambdaIn = FALSE
#   latRangeIn = latRange
#   normalModelIn = normalModel
#   taperedGPSDatIn = TRUE
#   distMatGPSIn = distMatGPS
#   distMatCSZIn = distMatCSZ
#   corGPSIn = corGPS
#   diffGPSTaperIn = diffGPSTaper
#   nKnotsGPSIn = nKnotsGPS
#   anisotropicIn = anisotropic
#   squareStrikeDistCszIn = squareStrikeDistCsz
#   squareDipDistCszIn = squareDipDistCsz
#   squareStrikeDistGpsIn = squareStrikeDistGps
#   squareDipDistGpsIn = squareDipDistGps
#   verboseIn = verbose
# }
# 
# generateSimulatedDataSets = function(nsim=100, params, cszDepths, corMatGPS=NULL, corMatCSZL=NULL, gpsDat=slipDatCSZ, 
#                                      subDat=dr1, fault=csz, useMVNApprox=FALSE, verbose=TRUE, muZeta=NULL, 
#                                      G=NULL, nKnots=5, normalizeTaper=TRUE, dStar=21000, useASLApprox=FALSE, 
#                                      useSubPrior=FALSE, useSlipPrior=FALSE, fauxG=NULL, constrLambda=TRUE, latRange=c(40,50), 
#                                      normalModel=FALSE, taperedGPSDat=FALSE, distMatGPS=NULL, distMatCSZ=NULL, 
#                                      corGPS=FALSE, diffGPSTaper=FALSE, nKnotsGPS=nKnots, anisotropic=FALSE, 
#                                      squareStrikeDistCsz=NULL, squareDipDistCsz=NULL, squareStrikeDistGps=NULL, 
#                                      squareDipDistGps=NULL) {
#   ##### get parameters
#   alpha = 1
#   if(!taperedGPSDat) {
#     tvec=NULL
#     if(length(params) <= 3 && (nKnots != 1)) {
#       vecLambda=FALSE
#       if(is.null(muZeta)) {
#         lambda = params[1]
#         muZeta = params[2]
#         sigmaZeta = params[3]
#         muZetaGPS = muZeta
#         muZetaCSZ = muZeta
#         vecMu=FALSE
#       }
#       else {
#         lambda = params[1]
#         sigmaZeta = params[2]
#         muZetaGPS = muZeta[1:nrow(gpsDat)]
#         muZetaCSZ = muZeta[(nrow(gpsDat)+1):length(muZeta)]
#         vecMu=TRUE
#       }
#     }
#     if(length(params) > 3 || (nKnots == 1)) {
#       vecLambda=TRUE
#       if(is.null(muZeta)) {
#         muZeta = params[1]
#         sigmaZeta = params[2]
#         splinePar = params[-c(1:2)]
#         tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar, latRange=latRange, fault=fault)
#         Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
#         lambda = Xi %*% splinePar # lambda is only used for checking if lambda < 0
#         muZetaGPS = muZeta
#         muZetaCSZ = muZeta
#         vecMu=FALSE
#       }
#       else {
#         stop("Code not yet written for spline taper fit along with vector muZeta")
#       }
#     }
#     else
#       tvec = taper(cszDepths, lambda=lambda, dStar=dStar, normalize=normalizeTaper)
#   }
#   else {
#     ## if the GPS data is considered tapered
#     # in this case the parameter vector is the same as otherwise, except it has an additional 
#     # scale parameter, phi, at the end.
#     
#     tvec=NULL
#     if(length(params) <= 3 && (nKnots != 1)) {
#       vecLambda=FALSE
#       if(is.null(muZeta)) {
#         lambda = params[1]
#         muZeta = params[2]
#         sigmaZeta = params[3]
#         muZetaGPS = muZeta
#         muZetaCSZ = muZeta
#         vecMu=FALSE
#       }
#       else {
#         lambda = params[1]
#         sigmaZeta = params[2]
#         muZetaGPS = muZeta[1:nrow(gpsDat)]
#         muZetaCSZ = muZeta[(nrow(gpsDat)+1):length(muZeta)]
#         vecMu=TRUE
#       }
#     }
#     if(length(params) > 3 || (nKnots == 1)) {
#       vecLambda=TRUE
#       if(is.null(muZeta)) {
#         muZeta = params[1]
#         sigmaZeta = params[2]
#         if(!anisotropic)
#           splinePar = params[-c(1:2, length(params))]
#         else
#           splinePar = params[-c(1:2, length(params) - 1, length(params))]
#         if(diffGPSTaper) {
#           splineParGPS = splinePar[-(1:nKnots)]
#           splinePar = splinePar[1:nKnots]
#         }
#         tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar, latRange=latRange, fault=fault)
#         Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
#         lambda = Xi %*% splinePar # lambda is only used for checking if lambda < 0
#         muZetaGPS = muZeta
#         muZetaCSZ = muZeta
#         vecMu=FALSE
#       }
#       else {
#         stop("Code not yet written for spline taper fit along with vector muZeta")
#       }
#     }
#     else
#       tvec = taper(cszDepths, lambda=lambda, dStar=dStar, normalize=normalizeTaper)
#     
#     # now get the scale and smoothness parameters
#     if(taperedGPSDat) {
#       if(!anisotropic)
#         phiZeta = params[length(params)]
#       else {
#         phiZeta = params[length(params) - 1]
#         alpha = params[length(params)]
#       }
#       nuZeta = 3/2
#     }
#     else {
#       out = getCorPar(normalModel)
#       phiZeta = out$phiZeta
#       nuZeta = out$nuZeta
#     }
#     
#     ## compute CSZ and GPS correlation matrices (use onlyUpper option to save time since will just take 
#     ## Cholesky decomposition of upper triangle anyway in both cases)
#     if(alpha >= 0.05 &&  phiZeta >= 0 && alpha <= 20) {
#       if(!anisotropic) {
#         coords = cbind(fault$longitude, fault$latitude)
#         if(phiZeta > 0) {
#           corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
#                                      onlyUpper=FALSE, distMat=distMatCSZ, smoothness=nuZeta)
#           corMatCSZL = t(chol(corMatCSZ))
#           
#           # GPS covariance matrix
#           coords = cbind(gpsDat$lon, gpsDat$lat)
#           corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
#                                      onlyUpper=FALSE, distMat=distMatGPS, smoothness=nuZeta)
#         }
#       }
#       else {
#         # compute distance matrices accounting for anisotropy
#         distMatCsz = sqrt(alpha^2 * squareStrikeDistCsz + (1 / alpha^2) * squareDipDistCsz)
#         distMatGps = sqrt(alpha^2 * squareStrikeDistGps + (1 / alpha^2) * squareDipDistGps)
#         
#         # now generate correlation matrices and the Cholesky decomposition if necessary
#         corMatCSZ = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
#                                    onlyUpper=FALSE, distMat=distMatCsz, smoothness=nuZeta)
#         corMatCSZL = t(chol(corMatCSZ))
#         
#         corMatGPS = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
#                                    onlyUpper=FALSE, distMat=distMatGps, smoothness=nuZeta)
#         corMatGPSL = t(corMatGPS)
#       }
#     }
#   }
#   lambda0 = 0.25
#   
#   # compute unit slip Okada seadef if necessary
#   if(is.null(G)) {
#     nx = 300
#     ny=  900
#     lonGrid = seq(lonRange[1], lonRange[2], l=nx)
#     latGrid = seq(latRange[1], latRange[2], l=ny)
#     G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
#   }
#   if(is.null(fauxG))
#     fauxG = getFauxG()
#   
#   # compute covariance of LOG of Zeta and its Cholesky decomp
#   SigmaZetaGPSL = sigmaZeta * corMatGPSL
#   SigmaZetaCSZL = sigmaZeta * corMatCSZL
#   
#   # compute gps taper
#   if(taperedGPSDat) {
#     XiGPS1 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
#     lambda = XiGPS1 %*% splinePar
#     if(diffGPSTaper) {
#       XiGPS2 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
#       lambda = lambda - XiGPS2 %*% splineParGPS
#     }
#     tvecGPS = taper(gpsDat$Depth, lambda, alpha=2, dStar=dStar, normalize=normalizeTaper)
#   }
#   else {
#     tvecGPS = rep(1, length(tvec))
#   }
#   
#   # Given distribution of zeta, generate subsidence data
#   generateSubsidenceData = function(nsim=1, muZeta, lambda, sigmaZeta, SigmaZeta, G, subDat=dr1, 
#                                     tvec=taper(csz$depth, lambda=lambda, normalize=normalizeTaper), 
#                                     dStar=21000, normalizeTaper=TRUE, normalModel=FALSE) {
#     
#     # This is the key step: approximate G %*% T %*% Zeta with a MVN
#     mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G, subDat=subDat, 
#                                      tvec=tvec, normalModel=normalModel)
#     mu = mvnApprox$mu
#     Sigma = mvnApprox$Sigma
#     
#     # add data noise to covariance matrix
#     diag(Sigma) = diag(Sigma) + subDat$Uncertainty^2
#     
#     # # get log likelihood
#     # y = -subDat$subsidence# MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
#     # lnLik = logLikGP(y - mu, chol(Sigma))
#     
#     # now generate the data
#     SigmaL = t(chol(Sigma))
#     eps = matrix(rnorm(nrow(Sigma) * nsim), nrow=nrow(Sigma))
#     uplift = sweep(matrix(SigmaL %*% eps, ncol=nsim), 1, mu, "+")
#     subsidence = -uplift
#     
#     subsidence
#   }
#   
#   # Do the same for the distribution of the gps data
#   generateGpsData = function(nsim=1, muZeta, SigmaZeta, gpsDat=slipDatCSZ, tvec = rep(1, length(sigmaXi)), 
#                       normalModel=FALSE, corGPS=FALSE) {
#     
#     ##### get conditional MLE of muXi
#     # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
#     # Transformation from additive error to  multiplicative lognormal model with asympototic 
#     # median and variance matching.
#     sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
#     
#     # estimate muXi MLE with inverse variance weighting
#     x = log(gpsDat$slip)
#     
#     ci = 1/sigmaXi^2
#     ci = ci/sum(ci)
#     if(normalModel) {
#       x = gpsDat$slip
#       sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
#       ci = 1/sigmaXi^2
#       ci = ci/sum(ci)
#       muXi = sum((log(x)- log(muZeta) - log(tvec))*ci)
#       gammaEst = exp(muXi)
#       
#       xCntr = x - gammaEst * tvec * muZeta
#       if(any(tvec != 1))
#         SigmaZeta = sweep(sweep(SigmaZeta, 1, tvec, "*"), 2, tvec, "*")
#       SigmaZeta = gammaEst^2 * SigmaZeta
#       sigmaXi = gpsDat$slipErr
#     }
#     else {
#       muXi = sum((x- muZeta - log(tvec))*ci)
#       
#       xCntr = x - muZeta - muXi - log(tvec)
#     }
#     
#     # if GPS measurement error is correlated, make it correlated in the same way as the latent process
#     if(corGPS) {
#       sigmaTZetas = sqrt(diag(SigmaZeta))
#       CXi = sweep(sweep(SigmaZeta, 1, 1/sigmaTZetas, "*"), 2, 1/sigmaTZetas, "*")
#       SigmaXi = sweep(sweep(CXi, 1, sigmaXi, "*"), 2, sigmaXi, "*")
#       Sigma = SigmaZeta + SigmaXi
#     }
#     else
#       Sigma = SigmaZeta + diag(sigmaXi^2)
#     
#     # # simulate data
#     # return(logLikGP(xCntr, chol(Sigma)))
#     # now generate the data
#     SigmaL = t(chol(Sigma))
#     eps = matrix(rnorm(nrow(Sigma) * nsim), nrow=nrow(Sigma))
#     lockrate = sweep(matrix(SigmaL %*% eps, ncol=nsim), 1, xCntr, "+")
#     lockrate
#   }
#   
#   ##### simulate the data sets
#   SigmaZetaCSZ = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
#   SigmaZetaGPS = SigmaZetaGPSL %*% t(SigmaZetaGPSL)
#   simsSubsidence = generateSubsidenceData(nsim, muZetaCSZ, lambda, sigmaZeta, SigmaZetaCSZ, G, subDat=subDat, 
#                             tvec=tvec, dStar=dStar, normalizeTaper=normalizeTaper, normalModel=normalModel)
#   simsGps = generateGpsData(nsim, muZetaGPS, SigmaZetaGPS, gpsDat, normalModel=normalModel, tvec=tvecGPS, corGPS=corGPS)
#   
#   ##### fit the model to the simulated data sets
#   fitParameters = matrix(nrow=length(params), ncol=nsim)
#   for(i in 1:nsim) {
#     # get the simulated data sets
#     thisSubDat = subDat
#     thisSubDat$subsidence = simsSubsidence[,i]
#     thisGpsDat = gpsDat
#     thisGpsDat = simsGps[,i]
#     
#     # fit the model (start at the known MLEs to save time)
#     fit = fitModel2(params, cszDepths, corMatGPS, corMatCSZL, gpsDat=thisGpsDat, 
#                     subDat=thisSubDat, fault, 500, useMVNApprox, verbose=FALSE, muZeta, 
#                     G, nKnots, normalizeTaper, dStar, useASLApprox, 
#                     useSubPrior, useSlipPrior, fauxG, constrLambda, latRange, 
#                     normalModel, taperedGPSDat, distMatGPS, distMatCSZ, 
#                     corGPS, diffGPSTaper, nKnotsGPS, anisotropic, 
#                     squareStrikeDistCsz, squareDipDistCsz, squareStrikeDistGps, 
#                     squareDipDistGps)
#     
#     # store the results
#     fitParameters[,i] = fit$optPar
#   }
#   
#   list(fitParameters=fitParameters, mean=apply(fitParameters, 1, mean), 
#        se=apply(fitParameters, 1, sd))
# }