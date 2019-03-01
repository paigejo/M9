# Archived functions that are no longer used, but held for reference.


##### fit GPS covariance range parameter (phiZeta) in km units
# zeta ~ exp(GP(muZeta, sigmaZeta * rhoZeta(.)))
# inputs:
# sigmaInit: initial guess for sigmaZeta
# phiInit: initial guess for correlation range parameters phi
# NOTE: assumes nuZeta, Matern smoothness parameter, is 3/2.  If phiInit
# is left as NULL, initial guess is set to max distance between locations 
# divided by 3
fitGPSCovariance = function(phiInit=NULL, doLog=TRUE) {
  ##### subset GPS data so it's only over the fault geometry
  slipDatCSZ = getFaultGPSDat()
  
  # guess sigma
  sigmaInit=sd(slipDatCSZ$slip)
  
  # set up and transform data
  coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  if(doLog)
    slip = log(slipDatCSZ$slip)
  else
    slip = slipDatCSZ$slip
  #   coords = cbind(slipDat$lon, slipDat$lat)
  #   logSlip = log(slipDat$slip)
  
  # NOTE: may want to further transform logSlip by subtracting
  #       log taper
  
  ##### guess lambda (noise/signal)
  # use delta method:
  # sqrt(n)(hat(theta) - theta) -> N(0, sigma^2)
  # f(x) = log(x), f'(x) = 1/x
  # sqrt(n)(log(hat(theta)) - log(theta)) -> N(0, sigma^2/theta^2)
  # take minimum because noise appears small, values are smooth
  if(doLog)
    noiseEst = min(slipDatCSZ$slipErr^2/slipDatCSZ$slip^2)
  else
    noiseEst = min(slipDatCSZ$slip^2)
  signalEst=  var(slip)
  lambdaGuess = noiseEst/(signalEst - noiseEst)
  
  #   optim.args = list(method = "BFGS", control = list(fnscale = -1, 
  #                                                     ndeps = rep(log(1.1), length(cov.params.guess) + 
  #                                                                   1), reltol = 1e-04, maxit = 10))
  
  # fit MLE range parameter with chordal distance
  distFun <<- function(x) { rdist.earth(x, miles=FALSE) }
  MLE = myMLESpatialProcess(coords, slip, lambda.start = lambdaGuess, 
                            cov.args=list(Covariance="Matern", Distance="rdist.earth", smoothness=1.5, 
                                          Dist.args=list(miles=FALSE)), 
                            Distance="distFun")
  phiMLE = MLE$summary[7]
  hess = MLE$hess
  lambda = MLE$lambda.fixed
  # mu = MLE$mu
  colnames(hess) = c("llambda", "lphi")
  
  # compute SE phi from hessian of log phi
  ## get variance this way rather than inverting full hessian since hess not full rank
  if(abs(hess[1,1]) < 1e-4) {
    phiVar = -1/hess[2,2]
  }
  else {
    parVar = solve(-hess)
    phiVar = parVar[2,2]
  }
  phiVar = phiVar*phiMLE^2
  phiSE = sqrt(phiVar)
  
  # return MLE
  list(phiMLE=phiMLE, lPhiHess=hess, phiSE=phiSE, lambda=lambda, mKrigObj=MLE)
}
# Results:
## fitGPSCovariance()
# $phiMLE
# [1] 232.5806
# 
# $lPhiHess
# llambda      lphi
# 0  0.000000
# theta       0 -9.248716
# 
# $phiSE
# [1] 76.47732
# 
# $lambda
# 1.22189e-28 

## for normal model:
# fitGPSCovariance(doLog=FALSE)
# $phiMLE
# [1] 174.3635
# 
# $lPhiHess
# llambda      lphi
# -3.507766  -9.19334
# theta -9.193340 -36.82721
# 
# $phiSE
# [1] 48.86449
# 
# $lambda
# 7.074392e-06 

## using all slipDat we get:
# $phiMLE
# [1] 28.19848
# 
# $lPhiHess
# llambda       lphi
# -24.48357  -11.78139
# theta -11.78139 -316.05643
# 
# $phiSE
# [1] 1.586147

getCorPar = function(normalModel=FALSE) {
  if(!normalModel) {
    phiZeta = 232.5722 # MLE based on fitGPSCovariance result
    nuZeta = 3/2 # we assume this is the Matern smoothness parameter
    lambda = 1.22189e-28 # ratio of nugget variance to sill
  }
  else if(normalModel) {
    phiZeta = 174.3635 # MLE based on fitGPSCovariance result
    nuZeta = 3/2 # we assume this is the Matern smoothness parameter
    lambda = 7.074392e-06 # ratio of nugget variance to sill
  }
  
  list(phiZeta=phiZeta, nuZeta=nuZeta, lambda=lambda)
}

##### function for performing full MLE fit of model
# params[1]: lambda
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: lambda0
# NOTE: params[5]: muXi is not required since it's estimated conditionally
# NOTE: currently the Okada model is the most time-consuming step in the optimization.
#       It may help to do block optimization, fixing lambda0 for a while before changing it.
#       Alternatively, we could assume lambda0=0.25 as is often done.
doFullFit = function(initParams=NULL, nsim=500, useMVNApprox=FALSE) {
  ##### Set up initial parameter guesses for the MLE optimization if necessary
  if(is.null(initParams)) {
    lambdaInit=2
    muZetaInit = log(20)
    # variance of a lognormal (x) is (exp(sigma^2) - 1) * exp(2mu + sigma^2)
    # then sigma^2 = log(var(e^x)) - 2mu - log(e^(sigma^2) - 1)
    # so sigma^2 \approx log(var(e^x)) - 2mu
    # sigmaZetaInit = sqrt(log(var(dr1$subsidence)) - 2*muZetaInit)
    sigmaZetaInit = 1
    lambda0Init = 0.25
    # muXiInit = log(2.25) # based on plots from exploratory analysis
    initParams = c(lambdaInit, muZetaInit, sigmaZetaInit, lambda0Init)
  }
  
  # get spatial correlation parameters (computed based on fitGPSCovariance)
  corPar = getCorPar()
  phiZeta = corPar$phiZeta
  nuZeta = corPar$nuZeta
  
  ##### subset GPS data so it's only over the fault geometry
  slipDatCSZ = getFaultGPSDat()
  
  ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
  # NOTE: only the upper triangle is computed for faster computation.  This is 
  #       sufficient for R's Cholesky decomposition
  # since we only want the correlation matrix, set sigmaZeta to 1
  #   coords = cbind(csz$longitude, csz$latitude)
  #   corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
  #                              onlyUpper=TRUE, Distance="rdist.earth", 
  #                              Dist.args=list(miles=FALSE), smoothness=nuZeta)
  tmpParams = rep(1, 5)
  # corMatCSZ = arealZetaCov(tmpParams, csz, nDown1=9, nStrike1=12)
  # load the precomputed correlation matrix
  arealCSZCor = getArealCorMat(fault)
  corMatCSZ = arealCSZCor
  corMatCSZL = t(chol(corMatCSZ))
  coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                             onlyUpper=TRUE, smoothness=nuZeta, 
                             Distance="rdist.earth", Dist.args=list(miles=FALSE))
  
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(csz)[,3]
  
  ##### Do optimization
  lastParams <<- NULL # used for updating Okada matrix only when necessary
  optimTable <<- NULL # table to save likelihood optimization steps
  controls = list(fnscale = -1, reltol=10^-5)
  opt = optim(initParams, fullDataLogLik, control=controls, hessian=TRUE, 
              cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, 
              phiZeta=phiZeta, slipDatCSZ=slipDatCSZ, nsim=nsim, 
              useMVNApprox=useMVNApprox)
  
  # test = fullDataLogLik(initParams, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, 
  # phiZeta=phiZeta, slipDatCSZ=slipDatCSZ)
  
  # get results of optimization
  lambdaMLE = opt$par[1]
  muZetaMLE = opt$par[2]
  sigmaZetaMLE = opt$par[3]
  lambda0MLE = opt$par[4]
  logLikMLE = opt$value
  hess = opt$hessian
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  logX = log(slipDatCSZ$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  muXiMLE = sum(logX*ci) - muZetaMLE
  
  # Return results
  list(MLEs=c(opt$par, muXiMLE), lambdaMLE=lambdaMLE, muZetaMLE = muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambda0MLE=lambda0MLE, 
       muXiMLE=muXiMLE, logLikMLE=logLikMLE, hess=hess, optimTable=optimTable)
}

##### Computes the full data log likelihood given the parameters.
## inputs:
# if muZeta is constant and passed as parameter:
# params[1]: lambda
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: lambda0
# if muZeta is vector with muZeta = c(muZetaPoint, muZetaAreal):
# params[1]: lambda
# params[2]: sigmaZeta
# params[3]: lambda0
fullDataLogLik = function(params, cszDepths, corMatGPS, corMatCSZL, phiZeta, gpsDat=slipDatCSZ, 
                          nsim=500, useMVNApprox=FALSE, muZeta=NULL) {
  ##### get parameters
  if(is.null(muZeta)) {
    lambda = params[1]
    muZeta = params[2]
    sigmaZeta = params[3]
    lambda0 = params[4]
    muZetaGPS = muZeta
    muZetaCSZ = muZeta
    vecMu=FALSE
  }
  else {
    lambda = params[1]
    sigmaZeta = params[2]
    lambda0 = params[3]
    muZetaGPS = muZeta[1:nrow(gpsDat)]
    muZetaCSZ = muZeta[(nrow(gpsDat)+1):length(muZeta)]
    vecMu=TRUE
  }
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(lastParams) || lastParams[4] != lambda0) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G <<- okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slip=1, poisson=lambda0)
  }
  
  ##### update lastParams
  lastParams <<- params
  
  ##### make sure params are in correct range and update optimTable
  if(lambda < 0 || sigmaZeta < 0 || lambda0 < 0) {
    newRow = rbind(c(params, -10000, -5000, -5000, NA))
    colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "lambda0", "lnLik", 
                         "subLnLik", "GPSLnLik", "lnLikSE")
    print(newRow)
    optimTable <<- rbind(optimTable, newRow)
    return(-10000) # return a very small likelihood
  }
  
  ##### compute covariance of EXPONENT of Zeta and its Cholesky decomp
  SigmaZetaGPS = sigmaZeta^2 * corMatGPS
  SigmaZetaCSZL = sigmaZeta * corMatCSZL
  
  ##### Compute Likelihood of subsidence and GPS data
  # sub = subsidenceLnLik(G, cszDepths, SigmaZetaCSZL, lambda, muZeta, nsim=nsim)
  if(!useMVNApprox)
    sub = subsidenceLnLikMod(G, cszDepths, SigmaZetaCSZL, lambda, muZetaCSZ, nsim=nsim)
  else {
    SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
    sub = subsidenceLnLikMod2(muZetaCSZ, lambda, sigmaZeta, SigmaZeta, G)
  }
  GPS = GPSLnLik(muZetaGPS, SigmaZetaGPS, gpsDat)
  lnLik = sub[1] + GPS
  lnLikSE = sub[2]
  
  ##### update parameter optimization table
  newRow = rbind(c(params, lnLik, sub[1], GPS, lnLikSE))
  if(!vecMu) {
    colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "lambda0", "lnLik", 
                         "subLnLik", "GPSLnLik", "lnLikSE")
  }
  else {
    colnames(newRow) = c("lambda", "sigmaZeta", "lambda0", "lnLik", 
                         "subLnLik", "GPSLnLik", "lnLikSE")
  }
  print(newRow)
  optimTable <<- rbind(optimTable, newRow)
  
  ##### return full data likelihood
  lnLik
}

# the newer version of this function only includes inputs relevant for the final analysis
# for fitting the model assuming the locking rate data product is tapered
fitModel2 = function(initParams=NULL, useMVNApprox=FALSE, gpsDat=slipDatCSZ, 
                     G=NULL, fauxG=NULL, subDat=dr1, fault=csz, nKnots=5, normalizeTaper=TRUE, 
                     dStar=21000, useGrad=FALSE, maxit=500, latRange=c(40, 50), 
                     normalModel=FALSE, doHess=TRUE, corGPS=FALSE, finalFit=FALSE, 
                     diffGPSTaper=FALSE, nKnotsGPS=nKnots, reltol=1e-8, anisotropic=FALSE, 
                     useGradForHess=(finalFit && doHess && useGrad), verbose=TRUE, 
                     nKnotsGamma=7, doGammaSpline=FALSE, nKnotsSD=5, doSDSpline=FALSE, 
                     dStarGPS=dStar) {
  
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
  verboseIn = verbose
  doGammaSplineIn = doGammaSpline
  nKnotsGammaIn = nKnotsGamma
  dStarGPSIn=dStarGPS
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
                nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn)
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
                      nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn, 
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
                    nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn, verbose=FALSE)
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
                                  nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn)
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
                       nKnotsGamma=nKnotsGammaIn, dStarGPS=dStarGPSIn)
  
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

##### Compute subsidence data likelihood
# Get likelihood of Y in:
# S[i] = t(d; lambda) * zeta[i],         zeta ~ exp(GP(muZeta, sigmaZeta * rhoZeta(.)))
# Y[i] = g(S, lambda0)[i] + eps[i],      eps  ~ MVN(0, diag(sigma))
#      = (G %*% T %*% Zeta)[i] + eps[i]
## inputs:
# G: okada matrix from okadaAll() function
# cszDepths: vector of depths of the csz subfault centers (the CSZ fault geometry)
# SigmaZetaL: lower triangular part of Cholesky decomp of covariance matrix GP in exponent of Zeta
# lambda: taper decay
# muZeta: mean of GP in exponent of Zeta
# nsim: number of Monte Carlo simulations of G %*% T %*% Zeta for likelihood estimation
# assumed values:
# dr1: subsidence data
# csz: CSZ fault geometry
# NOTE: lnLik SE is roughly 2.5 at ~500 simulations and it takes ~1sec/500sim
# NOTE: this is O(nsim) due to matrix multiply as currently implemented.  If 
#       coordinates are rounded to grid, could use circulant embedding to get 
#       faster?
subsidenceLnLik = function(G, cszDepths, SigmaZetaL, lambda, muZeta, nsim=3000) {
  # m is number of csz subfaults
  m = nrow(SigmaZetaL)
  
  # get log taper values
  logt = log(taper(cszDepths, lambda=lambda))
  
  # get updated mean vector of GP in exponent of Zeta
  meanZeta = logt + muZeta
  
  # simulate G %*% T %*% Zeta (each col of sims is a simulation.  nsim simulations for each event)
  Zs = matrix(rnorm(m*nsim*length(uniqueEvents)), nrow=m, ncol=nsim*length(uniqueEvents))
  GPs0 = SigmaZetaL %*% Zs # mean zero
  GPs = sweep(GPs0, 1, STATS = meanZeta, FUN="+")
  Ss = exp(GPs)
  sims = G %*% Ss
  
  # y is the data, uncert is its standard deviation
  y = -dr1$subsidence # MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  uncert = dr1$Uncertainty
  
  # Y[i] - (G %*% T %*% Zeta)[i] = eps[i]
  simEps = sweep(-sims, 1, y, "+")
  
  # calculate the likelihood for the given event
  eventLnLik = function(i) {
    eventName = uniqueEvents[i]
    inds = as.character(dr1$event) == eventName
    thisUncert = uncert[inds]
    
    # calculate simulation likelihoods
    # NOTE: due to lack of numerical instability must estimate each individual observation's
    #       likelihood, then log them, then sum them to get total log lik estimate.
    simLiks = function(simVals) {
      dnorm(simVals, sd=thisUncert)
    }
    colStart = (i-1)*nsim + 1
    colEnd = i*nsim
    
    # get the individual likelihoods for each simulated observation
    liks = apply(matrix(simEps[inds,colStart:colEnd], ncol=colEnd-colStart+1), 2, simLiks)
    
    # take mean and get SE for estimate of each observation's lik
    if(!is.null(dim(liks))) {
      obsLiks = apply(liks, 1, mean)
      obsLikSEs = apply(liks, 1, sd)/sqrt(nsim)
    }
    else{
      obsLiks = mean(liks)
      obsLikSEs = sd(liks)/sqrt(nsim)
    }
    
    # convert to log like estimates
    obsLnLiks = log(obsLiks)
    obsLnLikSEs = obsLikSEs/obsLiks
    
    # get full event log-lik estimate
    lnLik = sum(obsLnLiks)
    lnLikSE = sqrt(sum(obsLnLikSEs^2))
    return(c(lnLik=lnLik, lnLikSE=lnLikSE))
  }
  
  # get all event log likelihoods and standard errors
  lnLiks = sapply(1:length(uniqueEvents), eventLnLik)
  lnLikSEs = lnLiks[2,]
  lnLiks = lnLiks[1,]
  
  # now estimate total log likelihood.  Get SE as well
  lnLik = sum(lnLiks)
  lnLikSE = sqrt(sum(lnLikSEs^2))
  return(c(lnLik=lnLik, lnLikSE=lnLikSE))
  
  #   Don't do this anymore because it's not numerically stable.  Instead 
  #   use above log likelihood estimation
  #   # recursive function for computing expectation and variance of product of independent likelihood estimates
  #   prodLik = function(liks, likSEs) {
  #     # base case 1
  #     if(length(liks) == 1) {
  #       return(c(lik=liks, likSE = likSEs))
  #     }
  #     # calculate product update
  #     # E(XY) = EX EY
  #     # Var(XY) = [EX]^2 VarY + [EY]^2 VarX + Var(X) Var(Y)
  #     newLik = liks[1]*liks[2]
  #     likVars = likSEs^2
  #     newLikSE = sqrt(liks[1]^2 * likVars[2] + liks[2]^2 * likVars[1] + likVars[1]*likVars[2])
  #     
  #     # base case 2
  #     if(length(liks) == 2)
  #       return(c(lik=newLik, likSE=newLikSE))
  #     
  #     # general case
  #     prodLik(c(newLik, liks[3:length(liks)]), c(newLikSE, likSEs[3:length(likSEs)]))
  #   }
  #   
  #   return(prodLik(liks, likSEs))
}


# functions for fitting sdInflation MLE:

##### Computes the full data log likelihood given the parameters.
## inputs:
# params[1]: lambda
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: lambda0
# params[5]: sdInflation
# NOTE: params[5] was muXi but there is an easy conditional estimate of this
fullDataLogLik = function(params, cszDepths, corMatGPS, corMatCSZL, phiZeta, slipDatCSZ, 
                          nsim=500, useMVNApprox=FALSE) {
  ##### get parameters
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  sdInflation = params[5]
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(lastParams) || lastParams[4] != lambda0) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G <<- okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slip=1, poisson=lambda0)
  }
  
  ##### update lastParams
  lastParams <<- params
  
  ##### make sure params are in correct range and update optimTable
  if(lambda < 0 || sigmaZeta < 0 || lambda0 < 0 || sdInflation < 0) {
    newRow = rbind(c(params, -10000, -5000, -5000, NA))
    colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "lambda0", "lnLik", 
                         "subLnLik", "GPSLnLik", "lnLikSE")
    print(newRow)
    optimTable <<- rbind(optimTable, newRow)
    return(-10000) # return a very small likelihood
  }
  
  ##### compute covariance of exponent of Zeta and its Cholesky decomp
  SigmaZetaGPS = sigmaZeta^2 * corMatGPS
  SigmaZetaCSZL = sigmaZeta * corMatCSZL
  
  ##### Compute Likelihood of subsidence and GPS data
  # sub = subsidenceLnLik(G, cszDepths, SigmaZetaCSZL, lambda, muZeta, nsim=nsim)
  if(!useMVNApprox)
    sub = subsidenceLnLikMod(G, cszDepths, SigmaZetaCSZL, lambda, muZeta, nsim=nsim)
  else {
    SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
    sub = subsidenceLnLikMod2(muZeta, lambda, sigmaZeta, SigmaZeta, G)
  }
  GPS = GPSLnLik(muZeta, SigmaZetaGPS, slipDatCSZ, sdInflation=sdInflation)
  lnLik = sub[1] + GPS
  lnLikSE = sub[2]
  
  ##### update parameter optimization table
  newRow = rbind(c(params, lnLik, sub[1], GPS, lnLikSE))
  colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "lambda0", "sdInflation", "lnLik", 
                       "subLnLik", "GPSLnLik", "lnLikSE")
  print(newRow)
  optimTable <<- rbind(optimTable, newRow)
  
  ##### return full data likelihood
  lnLik
}

##### compute likelihood of GPS data
# Get likelihood of X in:
# X[i] = zeta[i] * xi[i], 
# ( or X[i]/t(d, lambda) = ... )
# where
# zeta ~ exp(GP(muZeta, sigmaZeta * rhoZeta(.)))
# xi   ~ exp(MVN(muXi, diag(sigmaXi)))
# and hence,
# logX ~ MVN(muXi + muZeta + log(t(d; lambda)), SigmaZeta + diag(sigmaXi))
# ( or logX ~ MVN(muXi + muZeta, SigmaZeta + diag(sigmaXi)) )
# The later seems preferable, since non-physical aspects of the GPS data 
# imply it has not yet been tapered, but this might warrant looking into.
## inputs:
# muXi: mean of MVN in exponent of xi (vector of GPS slip values).  NOTE:
# this input is no longer required.  It is estimated w/in the fxn
# muZeta: mean of GP in exponent of Zeta (constant)
# SigmaZeta: Covariance of exponent of zeta (of dimension |sigmaXi|)
# slipDatCSZ: GPS data restricted to CSZ fault geometry
# logt: log taper values evaluated at the GPS data
# sdInflation: the amount to add to the standard deviation of the GPS data
GPSLnLik = function(muZeta, SigmaZeta, slipDatCSZ, logt = rep(0, length(sigmaXi)), sdInflation=0) {
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*(slipDatCSZ$slipErr+sdInflation)^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  logX = log(slipDatCSZ$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  muXi = sum(logX*ci) - muZeta
  
  # center log GPS data so its mean is 0
  logXCntr = logX - muXi - muZeta - logt
  
  # get likelihood
  Sigma = SigmaZeta + diag(sigmaXi^2)
  logLikGP(logXCntr, chol(Sigma))
}


############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# fitting functions assuming lambda0 (okada model elasticity) is fixed at 0.25

##### function for performing full MLE fit of model
# params[1]: lambda
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: lambda0
# NOTE: params[5]: muXi is not required since it's estimated conditionally
# NOTE: currently the Okada model is the most time-consuming step in the optimization.
#       It may help to do block optimization, fixing lambda0 for a while before changing it.
#       Alternatively, we could assume lambda0=0.25 as is often done.
doFixedFit = function(initParams=NULL, nsim=500, useMVNApprox=FALSE, fitSDInflation=FALSE) {
  ##### Set up initial parameter guesses for the MLE optimization if necessary
  if(is.null(initParams)) {
    lambdaInit=2
    muZetaInit = log(20)
    # variance of a lognormal (x) is (exp(sigma^2) - 1) * exp(2mu + sigma^2)
    # then sigma^2 = log(var(e^x)) - 2mu - log(e^(sigma^2) - 1)
    # so sigma^2 \approx log(var(e^x)) - 2mu
    # sigmaZetaInit = sqrt(log(var(dr1$subsidence)) - 2*muZetaInit)
    sigmaZetaInit = 1
    # muXiInit = log(2.25) # based on plots from exploratory analysis
    sdInflationInit=0
    if(fitSDInflation)
      initParams = c(lambdaInit, muZetaInit, sigmaZetaInit, sdInflationInit)
    else
      initParams = c(lambdaInit, muZetaInit, sigmaZetaInit)
  }
  
  phiZeta = 232.5722 # MLE based on fitGPSCovariance result
  nuZeta = 3/2 # we assume this is the Matern smoothness parameter
  
  ##### subset GPS data so it's only over the fault geometry
  slipDatCSZ = getFaultGPSDat()
  
  ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
  # NOTE: only the upper triangle is computed for faster computation.  This is 
  #       sufficient for R's Cholesky decomposition
  # since we only want the correlation matrix, set sigmaZeta to 1
  #   coords = cbind(csz$longitude, csz$latitude)
  #   corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
  #                              onlyUpper=TRUE, Distance="rdist.earth", 
  #                              Dist.args=list(miles=FALSE), smoothness=nuZeta)
  tmpParams = rep(1, 5)
  # corMatCSZ = arealZetaCov(tmpParams, csz, nDown1=9, nStrike1=12)
  # load the precomputed correlation matrix
  load("arealCSZCor.RData")
  corMatCSZ = arealCSZCor
  corMatCSZL = t(chol(corMatCSZ))
  coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                             onlyUpper=TRUE, smoothness=nuZeta, 
                             Distance="rdist.earth", Dist.args=list(miles=FALSE))
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(csz)[,3]
  
  ##### Do optimization
  lastParams <<- NULL # used for updating Okada matrix only when necessary
  optimTable <<- NULL # table to save likelihood optimization steps
  controls = list(fnscale = -1, reltol=10^-5)
  opt = optim(initParams, fixedDataLogLik, control=controls, hessian=TRUE, 
              cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, 
              phiZeta=phiZeta, slipDatCSZ=slipDatCSZ, nsim=nsim, 
              useMVNApprox=useMVNApprox)
  
  # test = fullDataLogLik(initParams, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, 
  # phiZeta=phiZeta, slipDatCSZ=slipDatCSZ)
  
  # get results of optimization
  lambdaMLE = opt$par[1]
  muZetaMLE = opt$par[2]
  sigmaZetaMLE = opt$par[3]
  if(fitSDInflation) {
    SDInflationMLE = opt$par[4]
    MLEs = opt$par
  }
  else {
    SDInflationMLE = 0
    MLEs = c(opt$par, 0)
  }
  logLikMLE = opt$value
  hess = opt$hessian
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*(slipDatCSZ$slipErr+SDInflationMLE)^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  logX = log(slipDatCSZ$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  muXiMLE = sum(logX*ci) - muZetaMLE
  MLEs = c(MLEs, 0.25, muXiMLE)
  names(MLEs) = c("lambda", "muZeta", "sigmaZeta", "SDInflation", "lambda0", "muXi")
  
  # Return results
  list(MLEs=MLEs, lambdaMLE=lambdaMLE, muZetaMLE = muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambda0MLE=0.25, 
       muXiMLE=muXiMLE, logLikMLE=logLikMLE, hess=hess, optimTable=optimTable)
}

##### Computes the full data log likelihood given the parameters.
## inputs:
# params[1]: lambda
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: sdInflation
# NOTE: params[5] was muXi but there is an easy conditional estimate of this
fixedDataLogLik = function(params, cszDepths, corMatGPS, corMatCSZL, phiZeta, slipDatCSZ, 
                           nsim=500, useMVNApprox=FALSE) {
  ##### get parameters
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  sdInflation = params[4]
  if(is.na(sdInflation)) {
    sdInflation = 0
    params = c(params, sdInflation)
  }
  lambda0 = 0.25
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(lastParams)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G <<- okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slip=1, poisson=lambda0)
  }
  
  ##### update lastParams
  lastParams <<- params
  
  ##### make sure params are in correct range and update optimTable
  if(lambda < 0 || sigmaZeta < 0 || sdInflation < 0) {
    newRow = rbind(c(params, -10000, -5000, -5000, NA))
    colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "sdInflation", "lnLik", 
                         "subLnLik", "GPSLnLik", "lnLikSE")
    print(newRow)
    optimTable <<- rbind(optimTable, newRow)
    return(-10000) # return a very small likelihood
  }
  
  ##### compute covariance of exponent of Zeta and its Cholesky decomp
  SigmaZetaGPS = sigmaZeta^2 * corMatGPS
  SigmaZetaCSZL = sigmaZeta * corMatCSZL
  
  ##### Compute Likelihood of subsidence and GPS data
  # sub = subsidenceLnLik(G, cszDepths, SigmaZetaCSZL, lambda, muZeta, nsim=nsim)
  if(!useMVNApprox)
    sub = subsidenceLnLikMod(G, cszDepths, SigmaZetaCSZL, lambda=0.25, muZeta, nsim=nsim)
  else {
    SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
    sub = subsidenceLnLikMod2(muZeta, lambda, sigmaZeta, SigmaZeta, G)
  }
  
  GPS = GPSLnLik(muZeta, SigmaZetaGPS, slipDatCSZ, sdInflation=sdInflation)
  lnLik = sub[1] + GPS
  lnLikSE = sub[2]
  
  ##### update parameter optimization table
  newRow = rbind(c(params, lnLik, sub[1], GPS, lnLikSE))
  colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "sdInflation", "lnLik", 
                       "subLnLik", "GPSLnLik", "lnLikSE")
  print(newRow)
  optimTable <<- rbind(optimTable, newRow)
  
  ##### return full data likelihood
  lnLik
}

# this version of the function was cleaned up to only allow for the input 
# arguments that ended up being included in the final analysis
##### Computes the full data log likelihood given the parameters.
## inputs:
# if muZeta not provided directly, then:
# params[1]: lambda
# params[2]: muZeta
# params[3]: sigmaZeta
# else:
# params[1]: lambda
# params[2]: sigmaZeta
# NOTE: muXi was provided but there is an easy conditional estimate of this
fixedDataLogLik = function(params, cszDepths, corMatGPS=NULL, corMatCSZL=NULL, gpsDat=slipDatCSZ, 
                           subDat=dr1, fault=csz, nsim=500, useMVNApprox=FALSE, verbose=TRUE, muZeta=NULL, 
                           G=NULL, nKnots=5, normalizeTaper=TRUE, dStar=21000, useASLApprox=FALSE, 
                           useSubPrior=FALSE, useSlipPrior=FALSE, fauxG=NULL, constrLambda=TRUE, latRange=c(40,50), 
                           normalModel=FALSE, taperedGPSDat=FALSE, distMatGPS=NULL, distMatCSZ=NULL, 
                           corGPS=FALSE, diffGPSTaper=FALSE, nKnotsGPS=nKnots, anisotropic=FALSE, 
                           squareStrikeDistCsz=NULL, squareDipDistCsz=NULL, squareStrikeDistGps=NULL, 
                           squareDipDistGps=NULL, nKnotsGamma=7, doGammaSpline=FALSE, dStarGPS=dStar) {
  ##### get parameters
  alpha = 1
  if(!taperedGPSDat) {
    tvec=NULL
    if(length(params) <= 3 && (nKnots != 1)) {
      vecLambda=FALSE
      if(is.null(muZeta)) {
        lambda = params[1]
        muZeta = params[2]
        sigmaZeta = params[3]
        muZetaGPS = muZeta
        muZetaCSZ = muZeta
        vecMu=FALSE
      }
      else {
        lambda = params[1]
        sigmaZeta = params[2]
        muZetaGPS = muZeta[1:nrow(gpsDat)]
        muZetaCSZ = muZeta[(nrow(gpsDat)+1):length(muZeta)]
        vecMu=TRUE
      }
    }
    if(length(params) > 3 || (nKnots == 1)) {
      vecLambda=TRUE
      if(is.null(muZeta)) {
        muZeta = params[1]
        sigmaZeta = params[2]
        splinePar = params[-c(1:2)]
        tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar, latRange=latRange, fault=fault)
        Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
        lambda = Xi %*% splinePar # lambda is only used for checking if lambda < 0
        muZetaGPS = muZeta
        muZetaCSZ = muZeta
        vecMu=FALSE
      }
      else {
        stop("Code not yet written for spline taper fit along with vector muZeta")
      }
    }
    else
      tvec = taper(cszDepths, lambda=lambda, dStar=dStar, normalize=normalizeTaper)
  }
  else {
    ## if the GPS data is considered tapered
    # in this case the parameter vector is the same as otherwise, except it has an additional 
    # scale parameter, phi, at the end.
    
    tvec=NULL
    if(length(params) <= 3 && (nKnots != 1)) {
      vecLambda=FALSE
      if(is.null(muZeta)) {
        lambda = params[1]
        muZeta = params[2]
        sigmaZeta = params[3]
        muZetaGPS = muZeta
        muZetaCSZ = muZeta
        vecMu=FALSE
      }
      else {
        lambda = params[1]
        sigmaZeta = params[2]
        muZetaGPS = muZeta[1:nrow(gpsDat)]
        muZetaCSZ = muZeta[(nrow(gpsDat)+1):length(muZeta)]
        vecMu=TRUE
      }
    }
    if(length(params) > 3 || (nKnots == 1)) {
      vecLambda=TRUE
      if(is.null(muZeta)) {
        muZeta = params[1]
        sigmaZeta = params[2]
        if(!anisotropic)
          splinePar = params[-c(1:2, length(params))]
        else
          splinePar = params[-c(1:2, length(params) - 1, length(params))]
        if(diffGPSTaper) {
          splineParGPS = splinePar[-(1:nKnots)]
          splinePar = splinePar[1:nKnots]
        }
        tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar, latRange=latRange, fault=fault)
        Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
        lambda = Xi %*% splinePar # lambda is only used for checking if lambda < 0
        muZetaGPS = muZeta
        muZetaCSZ = muZeta
        vecMu=FALSE
      }
      else {
        stop("Code not yet written for spline taper fit along with vector muZeta")
      }
    }
    else
      tvec = taper(cszDepths, lambda=lambda, dStar=dStar, normalize=normalizeTaper)
    
    # now get the scale and smoothness parameters
    if(taperedGPSDat) {
      if(!anisotropic)
        phiZeta = params[length(params)]
      else {
        phiZeta = params[length(params) - 1]
        alpha = params[length(params)]
      }
      nuZeta = 3/2
    }
    else {
      out = getCorPar(normalModel)
      phiZeta = out$phiZeta
      nuZeta = out$nuZeta
    }
    
    ## compute CSZ and GPS correlation matrices (use onlyUpper option to save time since will just take 
    ## Cholesky decomposition of upper triangle anyway in both cases)
    if(alpha >= 0.05 &&  phiZeta >= 0 && alpha <= 20) {
      if(!anisotropic) {
        coords = cbind(fault$longitude, fault$latitude)
        if(phiZeta > 0) {
          corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                                     onlyUpper=FALSE, distMat=distMatCSZ, smoothness=nuZeta)
          corMatCSZL = t(chol(corMatCSZ))
          
          # GPS covariance matrix
          coords = cbind(gpsDat$lon, gpsDat$lat)
          corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                                     onlyUpper=FALSE, distMat=distMatGPS, smoothness=nuZeta)
        }
      }
      else {
        # compute distance matrices accounting for anisotropy
        distMatCsz = sqrt(alpha^2 * squareStrikeDistCsz + (1 / alpha^2) * squareDipDistCsz)
        distMatGps = sqrt(alpha^2 * squareStrikeDistGps + (1 / alpha^2) * squareDipDistGps)
        
        # now generate correlation matrices and the Cholesky decomposition if necessary
        corMatCSZ = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
                                   onlyUpper=FALSE, distMat=distMatCsz, smoothness=nuZeta)
        corMatCSZL = t(chol(corMatCSZ))
        
        corMatGPS = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
                                   onlyUpper=FALSE, distMat=distMatGps, smoothness=nuZeta)
      }
    }
  }
  lambda0 = 0.25
  
  ##### set up column names for parameter optimization table
  if(!diffGPSTaper)
    knotNames = paste0("beta", 1:(nKnots))
  else
    knotNames = c(paste0("beta", 1:(nKnots)), paste0("beta'", 1:(nKnots)))
  if(!vecMu) {
    if(length(params) > 3 || (nKnots == 1)) {
      if(!taperedGPSDat) {
        cNames = c("muZeta", "sigmaZeta", knotNames, 
                   "LL", "subLL", "GPSLL", "priorLL", "LLSE")
      }
      else {
        if(!anisotropic)
          cNames = c("muZeta", "sigmaZeta", knotNames, "phi", 
                     "LL", "subLL", "GPSLL", "priorLL", "LLSE")
        else
          cNames = c("muZeta", "sigmaZeta", knotNames, "phi", "alpha", 
                     "LL", "subLL", "GPSLL", "priorLL", "LLSE")
      }
    }
    else {
      if(!taperedGPSDat) {
        cNames = c("lambda", "muZeta", "sigmaZeta", "LL", 
                   "subLL", "GPSLL", "priorLL", "LLSE")
      }
      else {
        cNames = c("lambda", "muZeta", "sigmaZeta", "phiZeta", "LL", 
                   "subLL", "GPSLL", "priorLL", "LLSE")
      }
    }
  }
  else {
    if(length(params) > 2 || (nKnots == 1)) {
      if(!taperedGPSDat) {
        cNames = c("sigmaZeta", knotNames, 
                   "LL", "subLL", "GPSLL", "priorLL", "LLSE")
      }
      else {
        cNames = c("sigmaZeta", knotNames, "phiZeta", 
                   "LL", "subLL", "GPSLL", "priorLL", "LLSE")
      }
    }
    else {
      if(!taperedGPSDat) {
        cNames = c("lambda", "sigmaZeta", "LL", 
                   "subLL", "GPSLL", "priorLL", "LLSE")
      }
      else {
        cNames = c("lambda", "sigmaZeta", "phiZeta", "LL", 
                   "subLL", "GPSLL", "priorLL", "LLSE")
      }
    }
  }
  if(useSlipPrior)
    cNames = c(cNames, "slipPriorLL")
  if(useSubPrior)
    cNames = c(cNames, "subPriorLL")
  if(useSlipPrior || useSlipPrior)
    cNames = c(cNames, "priorLJac")
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  if(is.null(fauxG))
    fauxG = getFauxG()
  
  # plot taper (on original scale and GPS scale)
  latSeq = seq(latRange[1], latRange[2], l=500)
  splineMat = getSplineBasis(fault=fault, nKnots=nKnots, lats=latSeq, latRange=latRange)
  lambdaSeq = splineMat %*% splinePar
  lambdaSeqGPS = 0
  if(diffGPSTaper) {
    # add red line for GPS taper
    XiGPS2 = getSplineBasis(fault=fault, nKnots=nKnotsGPS, lats=latSeq, latRange=latRange)
    lambdaSeqGPS = lambdaSeq - XiGPS2 %*% splineParGPS
  }
  par(mfrow=c(1,2))
  plot(latSeq, lambdaSeq, type="l", ylim=range(c(lambdaSeq, lambdaSeqGPS, 0)))
  if(diffGPSTaper)
    lines(latSeq, lambdaSeqGPS, col="red")
  tmp = getTaperSpline(splinePar, nKnots=nKnots, dStar=dStar, latRange=latRange, fault=fault, normalize=normalizeTaper)
  plotFault(fault, tmp, varRange=c(0, 1))
  
  ##### make sure params are in correct range and update optimTable
  lambdaInRange=TRUE
  if(!vecLambda)
    lambdaInRange = (lambda>0) && all(lambda < 15)
  else if(constrLambda && any(lambdaSeq < 0) && any(lambdaSeqGPS < 0))
    lambdaInRange = FALSE
  else if(any(abs(lambdaSeq) > 15) && any(abs(lambdaSeqGPS) > 15))
    lambdaInRange=FALSE
  
  if(!lambdaInRange || (sigmaZeta <= 0) || (muZeta <= 0) || (phiZeta <= 0) || (alpha <= 0.05) || (alpha >= 20)) {
    if(useSlipPrior && useSubPrior) {
      if(!vecMu)
        newRow = rbind(c(params, -10000, -5000, -5000, 0, NA, NA, NA, NA))
      else
        newRow = rbind(c(params[1:2], splinePar, -10000, -5000, -5000, 0, NA, NA, NA, NA))
    }
    else {
      if(!vecMu)
        newRow = rbind(c(params, -10000, -5000, -5000, 0, NA))
      else
        newRow = rbind(c(params[1:2], splinePar, -10000, -5000, -5000, 0, NA))
    }
    colnames(newRow) = cNames
    if(verbose)
      print(newRow)
    optimTable <<- rbind(optimTable, newRow)
    return(-10000) # return a very small likelihood
  }
  
  ##### compute covariance of LOG of Zeta and its Cholesky decomp
  SigmaZetaGPS = sigmaZeta^2 * corMatGPS
  SigmaZetaCSZL = sigmaZeta * corMatCSZL
  
  ##### Compute Likelihood of subsidence and GPS data
  # sub = subsidenceLnLik(G, cszDepths, SigmaZetaCSZL, lambda, muZeta, nsim=nsim)
  SigmaZeta=NULL
  if(!useMVNApprox & !useASLApprox)
    sub = subsidenceLnLikMod(G, cszDepths, SigmaZetaCSZL, lambda=0.25, muZetaCSZ, nsim=nsim, 
                             tvec=tvec)
  else if(useASLApprox) {
    SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
    sub = subsidenceLnLikASL(muZetaCSZ, lambda, sigmaZeta, SigmaZeta, G, subDat=subDat, 
                             tvec=tvec, nsim=nsim)
  }
  else {
    SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
    sub = subsidenceLnLikMod2(muZetaCSZ, lambda, sigmaZeta, SigmaZeta, G, subDat=subDat, 
                              tvec=tvec, dStar=dStar, normalizeTaper=normalizeTaper, normalModel=normalModel)
  }
  
  if(taperedGPSDat) {
    XiGPS1 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
    lambda = XiGPS1 %*% splinePar
    if(diffGPSTaper) {
      XiGPS2 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
      lambda = lambda - XiGPS2 %*% splineParGPS
    }
    tvecGPS = taper(gpsDat$Depth, lambda, alpha=2, dStar=dStarGPS, normalize=normalizeTaper)
  }
  else {
    tvecGPS = rep(1, length(tvec))
  }
  
  GPS = GPSLnLik(muZetaGPS, SigmaZetaGPS, gpsDat, normalModel=normalModel, tvec=tvecGPS, corGPS=corGPS, 
                 doGammaSpline=doGammaSpline, nKnotsGamma=nKnotsGamma)
  lnLik = sub[1] + GPS
  lnLikSE = sub[2]
  priorLnLik = 0
  if(useSubPrior && useSlipPrior) {
    # if necessary, compute and add prior log likelihood as well
    if(is.null(SigmaZeta))
      SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
    if(is.null(tvec))
      tvec= taper(cszDepths, lambda, normalize=normalizeTaper, dStar=dStar)
    subPriorLik = priorSubLikSFrechet(muZeta, sigmaZeta, tvec, G, fauxG, subDat)
    slipPriorLik = priorSlipLikSFrechet(muZeta, sigmaZeta, tvec)
    priorLnLik = ifelse(useSubPrior, subPriorLik, 0) + ifelse(useSlipPrior, slipPriorLik, 0)
    lnLik = lnLik + priorLnLik
    
    # add the log jacobian factor
    Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
    priorJac = priorJacobian(muZeta, sigmaZeta, lambda, Xi, tvec, G, fauxG, subDat, 
                             fauxObs=getFauxObs(), fault, normalizeTaper, dStar)
    lnLik = lnLik + priorJac
  }
  else if(useSubPrior || useSlipPrior) {
    stop("Jacobian not derived in this case")
  }
  
  ##### update parameter optimization table
  if(!vecMu)
    newRow = c(params)
  else
    newRow = c(params[1:2], splinePar)
  newRow = c(newRow, lnLik, sub[1], GPS, priorLnLik, lnLikSE)
  if(useSlipPrior && useSubPrior)
    newRow = c(newRow, slipPriorLik, subPriorLik, priorJac)
  newRow = rbind(newRow)
  colnames(newRow) = cNames
  if(verbose)
    print(newRow, digits=5)
  optimTable <<- rbind(optimTable, newRow)
  
  ##### return full data likelihood
  lnLik
}

# an older version of this function. The newer version only includes inputs 
# that are actually incorporated in the final analysis
# compute the full gradient of the likelihood
# For now, I assume constant mean, so params is of the form:
# params[1]: muZeta
# params[2]: sigmaZeta
# params[3:(3+nKnots-1)]: taperPar
# NOTE: the last three inputs (phiZeta, useMVNApprox, and muZeta) are ignored, but 
#       are necessary since they are also inputs to fixedDataLogLikGrad.
fixedDataLogLikGrad = function(params, cszDepths, corMatGPS=NULL, corMatCSZL=NULL, gpsDat=slipDatCSZ, 
                               subDat=dr1, fault=csz, nsim=500, useMVNApprox=FALSE, verbose=TRUE, muZeta=NULL, 
                               G=NULL, nKnots=5, normalizeTaper=TRUE, dStar=21000, useASLApprox=FALSE, 
                               useSubPrior=FALSE, useSlipPrior=FALSE, fauxG=NULL, constrLambda=TRUE, latRange=c(40,50), 
                               normalModel=FALSE, taperedGPSDat=FALSE, distMatGPS=NULL, distMatCSZ=NULL, 
                               corGPS=FALSE, diffGPSTaper=FALSE, nKnotsGPS=nKnots, anisotropic=FALSE, 
                               squareDipDistGps=NULL, squareStrikeDistGps=NULL, 
                               squareDipDistCsz=NULL, squareStrikeDistCsz=NULL, returnRow=TRUE, 
                               nKnotsGamma=7, doGammaSpline=FALSE, dStarGPS=dStar, 
                               southernFun=nKnotsGPS == 1) {
  
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  muVecGPS = rep(muZeta, nrow(gpsDat))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots + diffGPSTaper*nKnotsGPS)]
  if(diffGPSTaper) {
    taperParGPS = taperPar[-(1:nKnots)]
    taperPar = taperPar[1:nKnots]
  }
  if(taperedGPSDat) {
    if(!anisotropic)
      phiZeta = params[length(params)]
    else
      phiZeta = params[length(params) - 1]
  }
  else{
    out = getCorPar(normalModel)
    phiZeta = out$phiZeta
  }
  nuZeta=3/2
  lambda0 = 0.25
  if(anisotropic)
    alpha = params[length(params)]
  else
    alpha = 1
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # compute spline basis matrix for subsidence data
  Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  XiGPS = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  lambdaGPS = XiGPS %*% taperPar
  
  # if GPS data has adjusted taper, update lambdas here:
  if(diffGPSTaper) {
    XiGPS2 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
    lambdaGPS = lambdaGPS - XiGPS2 %*% taperParGPS
  }
  
  # compute taper vector
  tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  tvecGPS = NULL
  if(taperedGPSDat)
    tvecGPS = taper(gpsDat$Depth, lambda=lambdaGPS, alpha=2, normalize=normalizeTaper, dStar=dStarGPS)
  
  # compute covariance matrices for latent process on fault (CSZ) and locking data (GPS)
  if(taperedGPSDat) {
    if(!anisotropic) {
      coords = cbind(fault$longitude, fault$latitude)
      corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                                 onlyUpper=FALSE, distMat=distMatCSZ, smoothness=nuZeta)
      corMatCSZL = t(chol(corMatCSZ))
      covMatCSZ = sigmaZeta^2 * corMatCSZ
      
      # GPS covariance matrix
      coords = cbind(gpsDat$lon, gpsDat$lat)
      corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                                 onlyUpper=FALSE, distMat=distMatGPS, smoothness=nuZeta)
    } else {
      # compute distance matrices accounting for anisotropy
      distMatCSZ = sqrt(alpha^2 * squareStrikeDistCsz + (1 / alpha^2) * squareDipDistCsz)
      distMatGPS = sqrt(alpha^2 * squareStrikeDistGps + (1 / alpha^2) * squareDipDistGps)
      
      # now generate correlation matrices and the Cholesky decomposition if necessary
      corMatCSZ = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
                                 onlyUpper=FALSE, distMat=distMatCSZ, smoothness=nuZeta)
      covMatCSZ = sigmaZeta^2 * corMatCSZ
      corMatCSZL = t(chol(corMatCSZ))
      
      corMatGPS = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
                                 onlyUpper=FALSE, distMat=distMatGPS, smoothness=nuZeta)
    }
  }
  else {
    arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
    corMatCSZ = arealCSZCor
    covMatCSZ = sigmaZeta^2 * corMatCSZ
  }
  
  # get zeta mean vector
  if(!normalModel)
    expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
  else
    expectZeta = muVecCSZ
  
  # get zeta covariance matrix (covariance of zeta, NOT log(zeta))
  if(!normalModel) {
    covZeta = exp(covMatCSZ) - 1
    covZeta = sweep(covZeta, 1, expectZeta, "*")
    covZeta = sweep(covZeta, 2, expectZeta, "*")
  }
  else
    covZeta = covMatCSZ
  
  # make sure tvec is in reasonable range
  # if(max(tvec < .075))
  #   return(rep(0, length(params)))
  
  # compute expected uplift
  expectY = G %*% (tvec * expectZeta)
  
  # get covariance gradient
  SigmaYGrad = covGrad(covZeta, lambda, tvec, Xi, G, fault, normalizeTaper, dStar, 
                       muZeta, sigmaZeta, subDat=subDat, normalModel=normalModel, 
                       distMat=distMatCSZ, phiZeta=phiZeta, corMatCSZ=corMatCSZ, 
                       anisotropic=anisotropic, alpha=alpha, 
                       squareStrikeDist=squareStrikeDistCsz, 
                       squareDipDist=squareDipDistCsz)
  
  subGrad = subLnLikGrad(expectZeta, expectY, lambda, Xi, G, covZeta, tvec=tvec, 
                         fault, subDat, normalizeTaper, dStar, pow=2, 
                         muZeta, sigmaZeta, SigmaYGrad=SigmaYGrad, normalModel=normalModel, 
                         phiZeta=phiZeta, distMat=distMatCSZ, corMatCSZ=corMatCSZ, 
                         anisotropic=anisotropic, alpha=alpha, 
                         squareStrikeDist=squareStrikeDistCsz, 
                         squareDipDist=squareDipDistCsz)
  GPSGrad = GPSLnLikGrad(muZeta, sigmaZeta, gpsDat, nKnots=nKnots, normalModel=normalModel, 
                         corMatGPS=corMatGPS, nPar=length(params), phiZeta=phiZeta,
                         distMatGPS=distMatGPS, taperedGPSDat=taperedGPSDat, tvecGPS=tvecGPS, 
                         dStar=dStarGPS, latRange=latRange, fault=fault, normalizeTaper=normalizeTaper, 
                         lambda=lambdaGPS, corGPS=corGPS, diffGPSTaper=diffGPSTaper, 
                         anisotropic=anisotropic, alpha=alpha, 
                         squareStrikeDist=squareStrikeDistGps, 
                         squareDipDist=squareDipDistGps, doGammaSpline=doGammaSpline, 
                         nKnotsGamma=nKnotsGamma)
  
  # adjust subsidence gradient to include zeros at the indices of the GPS taper parameters if necessary
  if(diffGPSTaper) {
    subGrad = matrix(c(subGrad[1:(2+nKnots)], rep(0, nKnotsGPS), subGrad[-(1:(2+nKnots))]), nrow=1)
  }
  
  # calculate prior gradient if necessary (DEPRECATED)
  priorGrad = 0
  if(useSlipPrior || useSubPrior) {
    subPriorGrad = 0
    slipPriorGrad = 0
    
    if(normalModel)
      stop("priors not yet derived for normal model")
    
    # get faux data varY gradient
    SigmaYGradFaux = covGrad(covZeta, lambda, tvec, Xi, fauxG, fault, normalizeTaper, dStar, 
                             muZeta, sigmaZeta, subDat=getFauxDat(subDat))
    
    if(useSubPrior) {
      if(is.null(fauxG))
        fauxG = getFauxG()
      subPriorGrad = priorSubLikGradSFrechet(muZeta, sigmaZeta, lambda, Xi, tvec, SigmaYGrad, G, fauxG, subDat, 
                                             fault=fault, normalizeTaper=normalizeTaper, dStar=dStar)
    }
    if(useSlipPrior)
      slipPriorGrad = priorSlipLikGradSFrechet(muZeta, sigmaZeta, tvec, lambda, Xi, fault, normalizeTaper, 
                                               dStar)
    if(useSubPrior && useSlipPrior) {
      fullG = rbind(G, fauxG)
      EYGrad = -meanGrad(expectZeta, lambda, Xi, fullG, fault, normalizeTaper, dStar, muZeta, sigmaZeta) # TAKE NEGATIVE FOR SUBSIDENCE
      priorJacGrad = priorJacobianGrad(muZeta, sigmaZeta, lambda, Xi, tvec, EYGrad, SigmaYGrad, SigmaYGradFaux, G, 
                                       fauxG, subDat, fauxObs=getFauxObs(), fault, normalizeTaper, dStar)
    }
    else if(useSubPrior || useSlipPrior) {
      stop("Jacobian not derived in this case")
    }
    
    priorGrad = priorJacGrad + useSubPrior*subPriorGrad + useSlipPrior*slipPriorGrad
  }
  
  if(returnRow)
    subGrad + GPSGrad + priorGrad
  else
    c(subGrad + GPSGrad + priorGrad)
}
# corGradNumeric(params, distMatCSZ, 1, anisotropic, squareStrikeDistCsz, squareDipDistCsz)
# corrGrad(distMatCSZ, phiZeta, corMatCSZ, anisotropic, alpha, squareStrikeDistCsz, squareDipDistCsz)
# jacobian(fixedDataLogLik, params, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, fault=fault, verbose=FALSE,
#          useMVNApprox=TRUE, G=G, nKnots=nKnots, dStar=dStar, useSubPrior=useSubPrior, useSlipPrior=useSlipPrior,
#          fauxG=fauxG, constrLambda=constrLambda, latRange=latRange, normalModel=normalModel, gpsDat=gpsDat, subDat=subDat,
#          anisotropic=anisotropic, squareDipDistGps=squareDipDistGps, squareStrikeDistGps=squareStrikeDistGps,
#          squareDipDistCsz=squareDipDistCsz, squareStrikeDistCsz=squareStrikeDistCsz, taperedGPSDat=taperedGPSDat,
#          normalizeTaper=normalizeTaper, corGPS=corGPS, diffGPSTaper=diffGPSTaper, nKnotsGPS=nKnotsGPS, nKnotsGamma=nKnotsGamma,
#          doGammaSpline=doGammaSpline)
# [,1]      [,2]      [,3]    [,4]      [,5]      [,6]      [,7]       [,8]      [,9]
# [1,] -0.2456049 0.3572869 -26.66681 3.43431 -14.42633 -19.16979 -9.013843 0.09438024 -22.38458
# jacobian(fixedDataLogLik, params, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, verbose=FALSE, method.args=list(d=.0000001), 
# phiZeta=phiZeta, useMVNApprox=TRUE, G=G, nKnots=nKnots, dStar=dStar, useSubPrior=useSubPrior, useSlipPrior=useSlipPrior, fauxG=fauxG)
# test = covYGradNumeric(params, G, nKnots, normalizeTaper, subDat, fault, colNum=1, dStar, latRange, normalModel,
#                        taperedGPSDat, distMatCSZ, cszDepths, anisotropic, squareDipDistCsz, squareStrikeDistCsz)
# test = covXGradNumeric(params, nKnots=nKnotsGPS, normalizeTaper=normalizeTaper, gpsDat=gpsDat,
#                        colNum=1, dStar=dStar, latRange=latRange, normalModel=normalModel,
#                        taperedGPSDat=taperedGPSDat, distMatGPS=distMatGPS, corGPS=corGPS, anisotropic=anisotropic,
#                        squareDipDist=squareDipDistGps, squareStrikeDist=squareStrikeDistGps)
# for new model:
# jacobian(fixedDataLogLik, params, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, fault=fault, verbose=FALSE,
#          useMVNApprox=TRUE, G=G, nKnots=nKnots, dStar=dStar, useSubPrior=useSubPrior, useSlipPrior=useSlipPrior,
#          fauxG=fauxG, constrLambda=constrLambda, latRange=latRange, normalModel=normalModel, gpsDat=gpsDat, subDat=subDat,
#          taperedGPSDat=taperedGPSDat, distMatGPS=distMatGPS, distMatCSZ=distMatCSZ, normalizeTaper=normalizeTaper,
#          corGPS=corGPS, diffGPSTaper=diffGPSTaper, nKnotsGPS=nKnotsGPS)

# for non-normalized taper:
# jacobian(fixedDataLogLik, params, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, fault=fault, verbose=FALSE,
#          useMVNApprox=TRUE, G=G, nKnots=nKnots, dStar=dStar, useSubPrior=useSubPrior, useSlipPrior=useSlipPrior,
#          fauxG=fauxG, constrLambda=constrLambda, latRange=latRange, normalModel=normalModel, gpsDat=gpsDat, subDat=subDat,
#          taperedGPSDat=taperedGPSDat, distMatGPS=distMatGPS, distMatCSZ=distMatCSZ, normalizeTaper=normalizeTaper, 
#          method.args = list(eps=1/2000000))

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# calculate subsidence log likelihood assuming multivariate normal approximation 
# to the lognormal distribution

source("approxDistn.R")

# Approximate the distribution of G %*% T %*% Zeta with a MVN to get likelihood
subsidenceLnLikMod2 = function(muZeta, lambda, sigmaZeta, SigmaZeta, G) {
  # get vector of taper values
  tvec = taper(csz$depth, lambda = lambda)
  
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G)
  mu = mvnApprox$mu
  Sigma = mvnApprox$Sigma
  
  # add data noise to covariance matrix
  diag(Sigma) = diag(Sigma) + dr1$Uncertainty^2
  
  # get log likelihood
  y = -dr1$subsidence# MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  lnLik = logLikGP(y - mu, chol(Sigma))
  
  return(c(lnLik=lnLik, lnLikSE=0))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions for computing best GPS data sd inflation

# This function fits the optimal additive standard deviation inflation for 
# the GPS data conditional on the rest of the parameters.
getBestSDInflation = function(params, initGuess=30) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # compute Okada matrix G
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda0)
  
  # get coordinates for GPS data
  xs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  
  # compute relevant covariances (don't add inflation factor to SigmaD yet)
  load("arealCSZCor.RData")
  SigmaZ = arealCSZCor * sigmaZeta^2
  SigmaZtoD = t(pointArealZetaCov(params, xs, csz, nDown=9, nStrike=12))
  SigmaD = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                          smoothness=nuZeta, Distance="rdist.earth", 
                          Dist.args=list(miles=FALSE)) * sigmaZeta^2
  
  ##### Do optimization
  # use Brent method since it's a 1-D optimization problem
  lastParams <<- NULL # used for updating Okada matrix only when necessary
  optimTable <<- NULL # table to save likelihood optimization steps
  controls = list(fnscale = -1, reltol=10^-5)
  opt = optim(initGuess, subLikGivenGPS, control=controls, hessian=TRUE, 
              params=params, SigmaZtoD=SigmaZtoD, SigmaD=SigmaD, G=G, 
              method="Brent", lower=0, upper=1e6)
  
  # get results of optimization
  bestInflation = opt$par
  hess = opt$hessian
  
  return(list(sdInflation=bestInflation, hessian=hess, optimTable=optimTable))
}

# function for computing likelihood of subsidence data given GPS data.  The 
# areal covariance matrix for zeta is unnecessary to input because the 
# correlation matrix is loaded from memory.
# NOTE: This function assumes subsidence data is multivariate normal!
# SigmaD is covariance matrix for the GPS data (without any desired added variance)
# SigmaZtoD is cross covariance matrix between zeta and GPS data
subLikGivenGPS = function(sdInflation, params, SigmaZtoD, SigmaD, G=NULL) {
  ##### make sure params are in correct range and update optimTable
  if(sdInflation < 0) {
    newRow = rbind(c(sdInflation, -10000))
    colnames(newRow) = c("sdInflation", "lnLik")
    print(newRow)
    optimTable <<- rbind(optimTable, newRow)
    return(-10000) # return a very small likelihood
  }
  
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # get data
  logX = log(slipDatCSZ$slip)
  Y = -dr1$subsidence
  
  # compute inflated variance sigmaXi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*(slipDatCSZ$slipErr+sdInflation)^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # update SigmaD with added variance
  diag(SigmaD) = diag(SigmaD) + sigmaXi^2
  
  # load covariance matrix of sigma over CSZ fault cells
  load("arealCSZCor.RData")
  SigmaZeta = arealCSZCor * sigmaZeta^2
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda0)
  }
  
  # get taper vector
  tvec = taper(csz$depth, lambda=lambda)
  
  # compute distribution of zeta conditional on GPS data
  predDistn = conditionalNormal(logX, muZeta, muZeta + muXi, SigmaZeta, SigmaD, SigmaZtoD)
  muc = predDistn$muc
  Sigmac = predDistn$Sigmac
  
  ## compute f(Y|X)
  # for now, we'll assume normality of Y even though that's not a great 
  # assumption.
  moments = estSubsidenceMeanCov(muc, lambda, sigmaZeta, Sigmac, G, tvec, TRUE, csz)
  muEst = moments$mu
  SigmaEstU = chol(moments$Sigma + diag(dr1$Uncertainty^2))
  logFYGivenX = logLikGP(Y - muEst, SigmaEstU)
  
  ##### update parameter optimization table
  newRow = rbind(c(sdInflation, logFYGivenX))
  colnames(newRow) = c("sdInflation", "logFYGivenX")
  print(newRow)
  optimTable <<- rbind(optimTable, newRow)
  
  logFYGivenX
}



################# exponential priors
##### exponential priors and gradients
# get prior on subsidence 95th quantile (exponential with 95% < 2m)
priorSubLikExp = function(muZeta, sigmaZeta, tvec, G=NULL, fauxG=NULL, subDat=dr1, fauxObs=getFauxObs(), rate=qexp(.95)/3) {
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  if(is.null(fauxG))
    fauxG = getFauxG(fauxObs)
  fullG = rbind(G, fauxG)
  
  # make our faux observations look like actual dr1 data (the event doesn't really matter in our case 
  # since we just want marginal variance and expectation)
  fauxDat = data.frame(list(Lon=fauxObs[,1], Lat=fauxObs[,2], Uncertainty=mean(subDat$Uncertainty), 
                            event="T1"))
  
  # combine relevant part of real and faux observations
  tmp = data.frame(list(Lon=subDat$Lon, Lat=subDat$Lat, Uncertainty=subDat$Uncertainty, event=subDat$event))
  fullDat = rbind(tmp, fauxDat)
  
  load("arealCSZCor.RData")
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, covMatCSZ, fullG, subDat=fullDat, 
                                   tvec=tvec)
  mu = -mvnApprox$mu # MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  Sigma = mvnApprox$Sigma
  
  ##### get 95th quantiles of subsidence
  sigmas = sqrt(diag(Sigma) + fullDat$Uncertainty^2)
  quants = qnorm(.95, mu, sigmas)
  
  # get maximum expected subsidence
  maxSub = max(quants)
  
  # return log prior likelihood
  dexp(maxSub, rate, log=TRUE)
}

# get prior subsidence likelihood gradient
priorSubLikGradExp = function(muZeta, sigmaZeta, lambda, Xi, tvec, SigmaYGrad, G=NULL, fauxG=NULL, subDat=dr1, 
                              fauxObs=getFauxObs(), rate=qexp(.95)/3, fault=csz, normalizeTaper=TRUE, 
                              dStar=21000) {
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  if(is.null(fauxG))
    fauxG = getFauxG(fauxObs)
  fullG = rbind(G, fauxG)
  
  ##### generate faux data
  # make our faux observations look like actual dr1 data (the event doesn't really matter in our case 
  # since we just want marginal variance and expectation)
  fauxDat = data.frame(list(Lon=fauxObs[,1], Lat=fauxObs[,2], Uncertainty=mean(subDat$Uncertainty), 
                            event="T1"))
  
  # combine relevant part of real and faux observations
  tmp = data.frame(list(Lon=subDat$Lon, Lat=subDat$Lat, Uncertainty=subDat$Uncertainty, event=subDat$event))
  fullDat = rbind(tmp, fauxDat)
  
  # get zeta mean vector and log zeta covariance
  load("arealCSZCor.RData")
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  expectZeta = exp(muZeta + diag(covMatCSZ)/2)
  
  # get zeta covariance matrix (covariance of zeta, NOT log(zeta))
  covZeta = exp(covMatCSZ) - 1
  covZeta = sweep(covZeta, 1, expectZeta, "*")
  covZeta = sweep(covZeta, 2, expectZeta, "*")
  
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, covMatCSZ, fullG, subDat=fullDat, 
                                   tvec=tvec)
  mu = -mvnApprox$mu # MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  Sigma = mvnApprox$Sigma
  
  ##### get 95th quantiles of subsidence
  sigmas = sqrt(diag(Sigma) + fullDat$Uncertainty^2)
  quants = qnorm(.95, mu, sigmas)
  
  # get maximum expected subsidence
  maxSubI = which.max(quants)
  
  ##### compute gradient
  # gradient of EY. MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  expectGrad = -meanGrad(expectZeta, lambda, Xi, fullG, fault, normalizeTaper, dStar, muZeta, sigmaZeta, index=maxSubI)
  
  # get gradient of SigmaY (for faux observations, we already have for true observations)
  SigmaYGradFaux = covGrad(covZeta, lambda, tvec, Xi, fauxG, fault, normalizeTaper, dStar, 
                           muZeta, sigmaZeta, subDat=fauxDat)
  
  # get only parts of SigmaY gradient we need (gradient of cov(Y_i, Y_i), no cross-covariances needed)
  getVarGrad = function(SigmaGrad) {
    n = dim(SigmaGrad)[1]
    p = dim(SigmaGrad)[3]
    inds1 = rep(1:n, p)
    inds2 = rep(1:p, each=n)
    matrix(SigmaGrad[cbind(inds1, inds1, inds2)], ncol=p)
  }
  varGradObs = getVarGrad(SigmaYGrad)
  varGradFaux = getVarGrad(SigmaYGradFaux)
  varGrad = rbind(varGradObs, varGradFaux)
  sdGrad = sweep(varGrad, 1, 1/2 * 1/sigmas, "*")
  
  # compute the gradient we're actually interested in: gradient of log prior likelihood
  -rate * (expectGrad + sdGrad[maxSubI,] * qnorm(.95))
}

# prior on expected earthquake slip.  95% chance that max slip < 40m
priorSlipLikExp = function(muZeta, sigmaZeta, tvec, rate=qexp(.95)/40) {
  
  # get log zeta variances
  load("arealCSZCor.RData")
  corMatCSZ = arealCSZCor
  varCSZ = sigmaZeta^2 * diag(corMatCSZ)
  sigmaCSZ = sqrt(varCSZ)
  
  # compute slip log scale means and standard deviations
  muSlip = muZeta + log(tvec)
  sigmaSlip = sigmaCSZ
  
  # compute 95th percentiles of distributions
  quants = qlnorm(.95, muSlip, sigmaSlip)
  
  # get log likelihood
  maxSlip = max(quants)
  dexp(maxSlip, rate, log=TRUE)
}

# prior on expected earthquake slip.  95% chance that max slip < 40m
priorSlipLikGradExp = function(muZeta, sigmaZeta, tvec, lambda, Xi, fault=csz, normalizeTaper=TRUE, 
                               dStar=21000, rate=qexp(.95)/40) {
  
  # get log zeta variances
  load("arealCSZCor.RData")
  corMatCSZ = arealCSZCor
  varCSZ = sigmaZeta^2 * diag(corMatCSZ)
  sigmaCSZ = sqrt(varCSZ)
  d = sqrt(diag(corMatCSZ))
  
  # compute slip log scale means and standard deviations
  muSlip = muZeta + log(tvec)
  sigmaSlip = sigmaCSZ
  
  # compute 95th percentiles of slip distributions
  quants = qlnorm(.95, muSlip, sigmaSlip)
  
  # get maximum slip
  maxSlipI = which.max(quants)
  q95 = quants[maxSlipI]
  
  ##### compute gradients
  # get taper gradient
  tGrad = taperGrad(fault$depth, Xi, lambda, dStar, normalizeTaper)
  
  # get 95th percentile gradient
  di = d[maxSlipI]
  ti = tvec[maxSlipI]
  tGradI = tGrad[maxSlipI,]
  q95Grad = matrix(c(q95, di * qnorm(.95) * q95, q95/ti * tGradI), nrow=1)
  
  # now get log slip prior gradient
  -rate * q95Grad
}
