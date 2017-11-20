


# Function for compute conditional mean and variance for normal distribution given data.
# Xp | Xd has (conditional) mean:
# muc = muP + SigmaPtoD %*% SigmaD^(-1) %*% (Xd - muD)
# and (conditional) variance:
# Sigmac = SigmaP - SigmaPtoD %*% SigmaD^(-1) %*% SigmaDtoP
conditionalNormal = function(Xd, muP, muD, SigmaP, SigmaD, SigmaPtoD) {
  
  # SigmaDInv = solve(SigmaD) # NOTE: luckily we only need to do this once.
  # muc = muP + SigmaPtoD %*% SigmaDTildeInv %*% (Xd - muD)
  # Sigmac = SigmaP - SigmaPtoD %*% SigmaDInv %*% SigmaDtoP
  
  # compute conditional mean and variance of zeta
  muc = muP + SigmaPtoD %*% solve(SigmaD, Xd - muD)
  Sigmac = SigmaP - SigmaPtoD %*% solve(SigmaD, t(SigmaPtoD))
  
  return(list(muc=muc, Sigmac=Sigmac))
}


# Function for generating simulations from the model predictive distribution given 
# the GPS data.  Same is predsGivenGPS but also generates predictions and 
# simulations at GPS coordinates
predsGivenGPSFull = function(params, nsim=100, muVec=NULL, gpsDat=slipDatCSZ, fault=csz, normalizeTaper=FALSE) {
  # get fit MLEs
  if(is.null(muVec)) {
    lambda = params[1]
    muZeta = params[2]
    sigmaZeta = params[3]
    lambda0 = params[4]
    muXi = params[5]
    muZetaGPS = muZeta
  }
  else {
    if(length(muVec) == 1)
      muVec = rep(muVec, nrow(gpsDat) + nrow(fault))
    
    lambda = params[1]
    sigmaZeta = params[3]
    lambda0 = params[4]
    muXi = params[5]
    muZeta=muVec
    muZetaGPS = muVec[1:nrow(gpsDat)]
  }
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # get log GPS data
  logX = log(gpsDat$slip)
  
  # get GPS data
  xs = cbind(gpsDat$lon, gpsDat$lat)
  
  # compute relevant covariances
  arealCSZCor = getArealCorMat(fault)
  SigmaB = arealCSZCor * sigmaZeta^2
  SigmaSB = pointArealZetaCov(params, xs, fault, nDown=9, nStrike=12)
  SigmaS = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                          smoothness=nuZeta, Distance="rdist.earth", 
                          Dist.args=list(miles=FALSE)) * sigmaZeta^2
  
  # now block them into predictions and data covariance matrices
  SigmaP = cbind(rbind(SigmaS, t(SigmaSB)), rbind(SigmaSB, SigmaB))
  SigmaD = SigmaS + diag(sigmaXi^2)
  SigmaPtoD = rbind(SigmaS, t(SigmaSB))
  
  # get block means
  muP = muZeta
  if(is.null(muVec))
    muD = rep(muZetaGPS + muXi, length(logX))
  else
    muD = muZetaGPS + muXi
  
  # compute conditional normal mean and standard error in mean estimate
  Xd = logX
  condDistn = conditionalNormal(Xd, muP, muD, SigmaP, SigmaD, SigmaPtoD)
  muc = condDistn$muc
  Sigmac = condDistn$Sigmac
  
  ##### now generate the conditional simulations.  Note that we still use the
  ##### marginal covariance structure, but we've conditionally updated the mean.
  # get prediction locations
  # get CSZ prediction coordinates
  xd = cbind(gpsDat$lon, gpsDat$lat)  # d is GPS locations
  xp = cbind(fault$longitude, fault$latitude) # p is fault areal locations
  nd = nrow(xd)
  np = nrow(xp)
  point = 1:nd
  areal = (nd+1):(np+nd)
  
  # Cholesky decomp used for simulations (don't add Sigmac because that 
  # can lead to non-physical results.  Here we just treat the mean as 
  # constant and use the marginal variability)
  # NOTE: deflate variance to account for increased variation due to
  #       conditional mean (using MSE = Var + bias^2 formula)
  # SigmaPred = SigmaP + Sigmac
  # varDeflation = 1 - mean((muc[point] - muZeta)^2)/sigmaZeta^2
  # SigmaPred = SigmaP * varDeflation
  # SigmaPred = SigmaP
  SigmaPred = Sigmac
  SigmaL = t(chol(SigmaPred))
  
  # generate predictive simulations
  zSims = matrix(rnorm(nsim*(nrow(xp)+nrow(xd))), nrow=nrow(xp)+nrow(xd), ncol=nsim)
  logZetaSims0 = SigmaL %*% zSims # each column is a zero mean simulation
  logZetaSims = sweep(logZetaSims0, 1, muc, "+") # add conditional mean to simulations
  zetaSims = exp(logZetaSims)
  tvec = taper(c(getFaultCenters(fault)[,3], gpsDat$Depth), lambda = lambda, normalize=normalizeTaper)
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # seperate areal average sims from GPS point location sims
  slipSimsGPS = slipSims[point,]
  slipSims = slipSims[areal,]
  
  # get mean slip prediction field
  meanSlip = exp(muc[areal] + diag(SigmaPred[areal,areal])/2) * tvec[areal]
  meanSlipGPS = exp(muc[point] + diag(SigmaPred[point,point])/2) * tvec[point]
  
  return(list(meanSlip=meanSlip, meanSlipGPS=meanSlipGPS, slipSims=slipSims, slipSimsGPS=slipSimsGPS, 
              Sigmac=Sigmac[areal, areal], muc=muc[areal], SigmacGPS = Sigmac[point, point], 
              mucGPS=muc[point], Sigma=SigmaPred[areal,areal], SigmaGPS=SigmaPred[point,point]))
}

# Function for generating simulations from the model predictive distribution given 
# the GPS data.  The predictive model is given in the presentation for the week of
# 08/30/17:
# log(zeta) | X has (conditional) mean:
# muc = muZeta + SigmaPtoD %*% (SigmaD + diag(sigmaXi^2))^(-1) %*% (log(X) - muZeta - muXi)
# and (conditional) variance"
# Sigmac = SigmaP - SigmaPtoD %*% (SigmaD + diag(sigmaXi^2))^(-1) %*% SigmaDtoP
predsGivenGPS = function(params, nsim=100, muVec=NULL, tvec=NULL, fault=csz, posNormalModel=FALSE, 
                         normalModel=posNormalModel, normalizeTaper=FALSE, dStar=28000) {
  # get fit MLEs
  if(is.null(muVec)) {
    lambda = params[1]
    muVecGPS = params[2]
    muVecCSZ = params[2]
    sigmaZeta = params[3]
    lambda0 = params[4]
    muXi = params[5]
  }
  else {
    lambda = params[1]
    sigmaZeta = params[3]
    lambda0 = params[4]
    muXi = params[5]
    muVecGPS = muVec[1:nrow(slipDatCSZ)]
    muVecCSZ = muVec[(nrow(slipDatCSZ)+1):length(muVec)]
  }
  
  # get the taper
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  
  # set other relevant parameters
  corPar = getCorPar(normalModel=normalModel)
  phiZeta = corPar$phiZeta
  nuZeta = corPar$nuZeta
  
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  if(!normalModel)
    sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  else
    sigmaXi = slipDatCSZ$slipErr
  
  # get log GPS data
  if(!normalModel)
    x = log(slipDatCSZ$slip)
  else
    x = slipDatCSZ$slip
  
  # get GPS data and CSZ prediction coordinates
  xd = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  xp = cbind(fault$longitude, fault$latitude)
  
  # compute relevant covariance matrices
  SigmaD = stationary.cov(xd, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                          theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  if(normalModel)
    SigmaD = muXi^2 * SigmaD
  SigmaD = SigmaD + diag(sigmaXi^2)
  SigmaDInv = solve(SigmaD) # NOTE: luckily we only need to do this once.
  
  # These lines have been replaced with the appropriate code for computing 
  # covariance between areal averages of zeta and points or other averages:
  #   SigmaPtoD = stationary.cov(xp, xd, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
  #                              theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  #   SigmaDtoP = t(SigmaPtoD)
  #   SigmaP = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
  #                           theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  arealCSZCor = getArealCorMat(fault, normalModel = normalModel)
  SigmaP = arealCSZCor * sigmaZeta^2
  SigmaDtoP = pointArealZetaCov(params, xd, fault, nDown=9, nStrike=12, normalModel=normalModel)
  if(normalModel) {
    SigmaP = muXi^2 * SigmaP
    SigmaDtoP = muXi^2 * SigmaDtoP
  }
  SigmaPtoD = t(SigmaDtoP)
  
  # compute conditional mean and variance of zeta and take Cholesky decomp
  # NOTE: use total variance, conditional covariance gives the SEs for the 
  #       conditional mean estimate.
  if(!normalModel)
    muX = muVecGPS + muXi
  else
    muX = muVecGPS * muXi
  muc = muVecCSZ + SigmaPtoD %*% SigmaDInv %*% (x - muX)
  Sigmac = SigmaP - SigmaPtoD %*% SigmaDInv %*% SigmaDtoP
  #   varDeflation = 1 - mean((muc - muZeta)^2)/sigmaZeta^2
  #   SigmaPred = SigmaP * varDeflation # total variance is conditional variance
  SigmaPred = Sigmac
  SigmaPredL = t(chol(SigmaPred))
  
  # # generate predictive simulations
  # zSims = matrix(rnorm(nsim*nrow(xp)), nrow=nrow(xp), ncol=nsim)
  # logZetaSims = sweep(SigmaPredL %*% zSims, 1, muc, FUN="+") # add muc to each zero mean simulation
  # if(!normalModel)
  #   zetaSims = exp(logZetaSims)
  # else
  #   zetaSims = logZetaSims
  # slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # generate predictive simulations
  notAllPos=TRUE
  zetaSims = matrix(-1, nrow=nrow(xp), ncol=nsim)
  while(notAllPos) {
    # generate simulations until all slips are positive, if necessary
    negCol = function(simCol) {
      any(simCol < 0)
    }
    negCols = apply(zetaSims, 2, negCol)
    nNewSims = sum(negCols)
    
    zSims = matrix(rnorm(nNewSims*nrow(xp)), nrow=nrow(xp), ncol=nNewSims)
    logZetaSims = sweep(SigmaPredL %*% zSims, 1, muc, "+") # add muZeta to each zero mean simulation
    if(!normalModel)
      zetaSims[,negCols] = exp(logZetaSims)
    else
      zetaSims[,negCols] = logZetaSims
    
    notAllPos =  any(zetaSims < 0) && posNormalModel
  }
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # get mean slip prediction field
  if(!normalModel)
    meanSlip = exp(muc + diag(SigmaPred)/2) * tvec
  else if(!posNormalModel)
    meanSlip = muc * tvec
  else {
    meanSlip = apply(slipSims, 1, mean)
    if(nsim < 1000)
      warning("mean slip estimates may be poor with positive normal mode for <1000 simulations")
  }
  
  ##### generate predictions at GPS locations
  mucGPS = muVecGPS + SigmaD %*% SigmaDInv %*% (x - muX)
  SigmacGPS = SigmaD - (SigmaD - 2*diag(sigmaXi^2) + diag(sigmaXi^4) %*% SigmaDInv)
  SigmaPredGPS = SigmacGPS + SigmaD
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, tvec=tvec, muc=muc, Sigmac=Sigmac, 
              mucGPS=mucGPS, SigmacDiagGPS = diag(SigmacGPS), Sigma=SigmaPred, SigmaGPS=SigmaPredGPS))
}

# generate predictions given only the parameter MLEs (no GPS or subsidence data)
preds = function(params, nsim=100, fault=csz, muVec=NULL, tvec=rep(params[1], nrow(fault)), 
                 posNormalModel=FALSE, normalModel=posNormalModel, phiZeta=NULL, 
                 taperedGPSDat=FALSE) {
  # get parameters
  if(is.null(muVec)) {
    lambda = params[1]
    muZeta = params[2]
    sigmaZeta = params[3]
    muXi = params[5]
    muZetaGPS = rep(muZeta, nrow(slipDatCSZ))
    muZetaCSZ = rep(muZeta, nrow(fault))
  }
  else {
    lambda = params[1]
    sigmaZeta = params[3]
    muXi = params[5]
    muZeta = muVec
    muZetaGPS = muVec[1:nrow(slipDatCSZ)]
    muZetaCSZ = muVec[(nrow(slipDatCSZ)+1):length(muVec)]
  }
  
  # set other relevant parameters
  if(is.null(phiZeta)) {
    corPar = getCorPar(normalModel=normalModel)
    phiZeta = corPar$phiZeta
    nuZeta = corPar$nuZeta
  }
  else {
    nuZeta = 3/2
  }
  
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  if(!normalModel)
    sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  else
    sigmaXi = slipDatCSZ$slipErr
  
  # get CSZ prediction coordinates
  xp = cbind(fault$longitude, fault$latitude)
  
  # compute relevant covariance matrices
  # NOTE: previous code replaced with code calculating covariance between 
  #       areal averages of zeta
#   Sigma = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
#                          theta=phiZeta, smoothness=nuZeta, onlyUpper=TRUE) * sigmaZeta^2
  # Load the precomputed correlation matrix.
  if(!taperedGPSDat) {
    arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  }
  else {
    phiZ = params[length(params)]
    coordsZ = cbind(fault$longitude, fault$latitude)
    distMatZ = rdist.earth(coordsZ, miles=FALSE)
    arealCSZCor = stationary.cov(coordsZ, Covariance="Matern", theta=phiZ,
                             onlyUpper=FALSE, distMat=distMatZ, smoothness=3/2)
  }
  Sigma = arealCSZCor * sigmaZeta^2
  SigmaL = t(chol(Sigma))
  
  # generate predictive simulations
  notAllPos=TRUE
  zetaSims = matrix(-1, nrow=nrow(xp), ncol=nsim)
  nNewSims = nsim
  while(notAllPos) {
    # generate simulations until all slips are positive, if necessary
    negCol = function(simCol) {
      any(simCol < 0)
    }
    negCols = apply(zetaSims, 2, negCol)
    if(nNewSims != sum(negCols)) {
      nNewSims = sum(negCols)
      print(paste0("number of simulations remaining: ", nNewSims))
    }
    
    zSims = matrix(rnorm(nNewSims*nrow(xp)), nrow=nrow(xp), ncol=nNewSims)
    logZetaSims = sweep(SigmaL %*% zSims, 1, muZetaCSZ, "+") # add muZeta to each zero mean simulation
    if(!normalModel)
      zetaSims[,negCols] = exp(logZetaSims)
    else
      zetaSims[,negCols] = logZetaSims
    
    notAllPos =  any(zetaSims < 0) && posNormalModel
  }
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # get mean slip prediction field
  if(!normalModel)
    meanSlip = exp(muZetaCSZ + diag(Sigma)/2) * tvec
  else if(!posNormalModel)
    meanSlip = muZetaCSZ * tvec
  else {
    meanSlip = apply(slipSims, 1, mean)
    if(nsim < 1000)
      warning("mean slip estimates may be poor with positive normal mode for <1000 simulations")
  }
  
  # compute covariance at GPS locations
  xd = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  SigmaD = stationary.cov(xd, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                          theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, Sigma=Sigma, Sigmac=Sigma, muc=muZetaCSZ, 
              SigmacGPS = SigmaD, mucGPS=muZetaGPS))
}

# generate pointwise predictions over a grid given only the parameter MLEs (no GPS or subsidence data)
predsPoint = function(params, nsim=100, fault=csz, muVec=NULL, lonLatGrid=NULL, dStar=26000, nKnots=5, 
                      normalizeTaper=FALSE) {
  # get parameters
  if(is.null(muVec)) {
    lambda = params[1]
    muZeta = params[2]
    sigmaZeta = params[3]
    muXi = params[5]
    muZetaGPS = rep(muZeta, nrow(slipDatCSZ))
    muZetaCSZ = rep(muZeta, nrow(fault))
  }
  else {
    lambda = params[1]
    sigmaZeta = params[3]
    muXi = params[5]
    muZeta = muVec
    muZetaGPS = muVec[1:nrow(slipDatCSZ)]
    muZetaCSZ = muVec[(nrow(slipDatCSZ)+1):length(muVec)]
  }
  # get tvec
  if(is.na(lambda))
    taperPar = params[6:(5+nKnots)]
  
  # make the grid over which to make pointwise estimates (~7040 points)
  if(is.null(lonLatGrid)) {
    latRange=range(slipDatCSZ$lat)
    lonRange=range(slipDatCSZ$lon)
    nx = 80
    ny = 240
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    lonLatGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
    lonLatGrid = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2]))
    
    # make sure we only generate predictions over our fault geometry (which has some gaps in it)
    lonLatGrid = as.matrix(getFaultGPSDat(lonLatGrid))
  }
  
  muZetaCSZ = rep(muZeta, nrow(lonLatGrid))
  
  # get predicted depth at prediction locations
  phiZeta = 232.5722
  out = fastTps(cbind(slipDat$lon, slipDat$lat), slipDat$Depth, m=3, theta=phiZeta, lon.lat=TRUE, 
                Dist.args=list(miles=FALSE, method="greatcircle"))
  depths = predict(out, lonLatGrid)
  negDepth = depths < 0
  depths[negDepth] = 0
  
  # get tvec
  predData = data.frame(list(latitude=lonLatGrid[,2], depth=depths))
  tvec = getTaperSpline(taperPar, fault=predData, dStar=dStar, normalize=normalizeTaper)
  
  # compute covariance at prediction locations
  nuZeta=3/2
  Sigma = stationary.cov(lonLatGrid, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                          theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  
  # generate predictive simulations
  SigmaL = t(chol(Sigma))
  zSims = matrix(rnorm(nsim*nrow(lonLatGrid)), nrow=nrow(lonLatGrid), ncol=nsim)
  logZetaSims = sweep(SigmaL %*% zSims, 1, muZetaCSZ, "+") # add muZeta to each zero mean simulation
  zetaSims = exp(logZetaSims)
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # get mean slip prediction field
  meanSlip = exp(muZeta + diag(Sigma)/2) * tvec
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, Sigma=Sigma, lonLatGrid=lonLatGrid))
}

# generate predictions given only the parameter MLEs (no GPS or subsidence data) at 
# GPS locations as well as areal averages over the fault grid cells.
predsArealAndLoc = function(params, nsim=100, fault=csz, gpsDat=slipDatCSZ, muVec=params[2], 
                            normalizeTaper=FALSE, dStar=28000) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # get vector mean
  if(length(muVec) == 1)
    muVec = rep(muVec, nrow(fault)+nrow(gpsDat))
  muVecGPS = muVec[1:nrow(gpsDat)]
  muVecCSZ = muVec[(nrow(gpsDat)+1):length(muVec)]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # get CSZ prediction coordinates
  xd = cbind(gpsDat$lon, gpsDat$lat)
  xp = cbind(fault$longitude, fault$latitude)
  nd = nrow(xd)
  np = nrow(xp)
  
  # compute relevant covariance matrices
  # NOTE: previous code replaced with code calculating covariance between 
  #       areal averages of zeta
  #   Sigma = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
  #                          theta=phiZeta, smoothness=nuZeta, onlyUpper=TRUE) * sigmaZeta^2
  if(! identical(fault, csz))
    SigmaCSZ = arealZetaCov(params, fault, nDown1=9, nStrike1=12)
  else {
    # in this case, it takes too long to compute.  Load the precomputed correlation matrix.
    arealCSZCor = getArealCorMat(fault)
    SigmaCSZ = arealCSZCor * sigmaZeta^2
  }
  SigmaD = stationary.cov(xd, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                          theta=phiZeta, smoothness=nuZeta, onlyUpper=TRUE) * sigmaZeta^2
  SigmaDtoP = pointArealZetaCov(params, xd, fault, nDown=9, nStrike=12)
  SigmaPtoD = t(SigmaDtoP)
  Sigma = cbind(rbind(SigmaD, SigmaPtoD), rbind(SigmaDtoP, SigmaCSZ))
  SigmaL = t(chol(Sigma))
  
  # generate predictive simulations
  zSims = matrix(rnorm(nsim*(nrow(xp)+nrow(xd))), nrow=nrow(xp)+nrow(xd), ncol=nsim)
  logZetaSims = sweep(SigmaL %*% zSims, 1, muVec, "+") # add muZeta (vector) to each zero mean simulation
  zetaSims = exp(logZetaSims)
  tvec = taper(c(gpsDat$Depth, getFaultCenters(fault)[,3]), lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # seperate areal average sims from GPS point location sims
  point = 1:nd
  areal = (nd+1):(nd+np)
  slipSimsGPS = slipSims[point,]
  slipSims = slipSims[areal,]
  
  # get mean slip prediction field
  meanSlip = rep(exp(muZeta) + SigmaCSZ^2/2, nrow(xp)) * tvec[areal]
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, slipSimsGPS=slipSimsGPS, Sigmac=Sigma, 
              muc=rep(muZeta, nrow(fault)), SigmacGPS = SigmaD, mucGPS=rep(muZeta, nrow(xd))))
}

# Compute subsidence from the prediction simulations using the Okada model (NOTE: 
# returned ``subsidence'' is really uplift here).
# Preds is a list with elements named meanSlip (vector) and slipSims (matrix)
predsToSubsidence = function(params, preds, fault=csz, useMVNApprox=TRUE, G=NULL, 
                             subDat=dr1, posNormalModel=FALSE, normalModel=posNormalModel, tvec=NULL, 
                             normalizeTaper=FALSE, dStar=28000) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # get predictions from input list
  meanSlip = preds$meanSlip
  slipSims = preds$slipSims
  
  # get taper if necessary
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # transform slips into subsidences
  meanSub = G %*% cbind(meanSlip)
  subSims = G %*% slipSims
  
  # approximate upper and lower 95% quantiles (with either MVN approximate or simulations)
  # NOTE: use preds$Sigma not preds$Sigmac since Sigmac gives the covariance in mean estimate
  sigmaEps = subDat$Uncertainty
  if(normalModel && !posNormalModel && useMVNApprox) {
    subMVN = estSubsidenceMeanCov(preds$muc, lambda, sigmaZeta, preds$Sigma, G, fault=fault, subDat=subDat, 
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
    subMVN = estSubsidenceMeanCov(preds$muc, lambda, sigmaZeta, preds$Sigma, G, subDat=subDat, fault=fault, 
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

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
##### functions for bootstrapping

# Function for computing standard errors for nonnegative weighted least squares.
# Assumes regression coefficients follow a folded normal distribution, which is 
# why traditional WLS estimates are needed in addition to the NNWLS estimates.
# - beta is the vector of nnls coefficient estimates
# - betaWLS is the vector of WLS coefficient estimates
# - V is the diagonal of the data variance matrix.  Here we assume the 
#   variance matrix is itself diagonal as in weighted least squares.
# - A is the design matrix.  In our case this is G %*% T
getNNLSSE = function(beta, betaWLS, V, A) {
  # make the variance matrix
  Vinv = diag(V^(-1))
  Sigma = pseudoinverse(A %*% Vinv %*% t(A), 10^(-7))
  
  # compute the parameter standard errors under normal assumption
  SEs = sqrt(diag(Sigma))
  
  # Under our conditional normal distribution (positive normal) 
  # the the distribution will be cut off at 0 so we need to modify 
  # the SEs to account for this:
  
  
  return(SEs)
}
# NOTE: might want to also make a function for confidence intervals

# bootstrap residuals to get standard errors for NNLS fit
nnlsBootstrapSE = function(A, y, yFitted, nSamples=10000) {
  resids = y - yFitted
  n = length(y)
  
  # store regression estimates (beta) in a matrix with nSamples rows
  betas = matrix(nrow=nSamples, ncol=ncol(A))
  for(i in 1:nSamples) {
    # print progress
    if(i %% 500 == 0)
      print(paste0("iteration ", i, "/", nSamples))
    
    # resample residuals and get nnls parameter estimates
    epsStar = sample(resids, n, replace = TRUE)
    yStar = yFitted + epsStar
    nnlsMod = nnls(A, yStar)
    betas[i,] = nnlsMod$x
  }
  
  # get boostrap marginal standard errors
  SEs = apply(betas, 2, sd)
  
  return(SEs)
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
##### functions for generating predictions when conditioning on all the data

##### Generate samples from the predictive distribution conditional on all the data 
##### assuming log zeta, log X, and Y is all multivariate normal.
genFullPredsMVN = function(params, nsim=1000, sdAdd=0, normalizeTaper=FALSE, dStar=28000) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*(slipDatCSZ$slipErr+sdAdd)^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # get data
  logX = log(slipDatCSZ$slip)
  Y = -dr1$subsidence
  
  # get Okada linear transformation matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda0)
  
  # get taper vector
  tvec = taper(csz$depth, lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  
  # compute G %*% T
  GT = sweep(G, 2, tvec, "*")
  
  # get coordinates for GPS data
  xs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  
  # compute relevant covariances
  arealCSZCor = getArealCorMat(fault)
  SigmaB = arealCSZCor * sigmaZeta^2
  SigmaSB = pointArealZetaCov(params, xs, csz, nDown=9, nStrike=12)
  SigmaS = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                          smoothness=nuZeta, Distance="rdist.earth", 
                          Dist.args=list(miles=FALSE)) * sigmaZeta^2
  SigmaYMod = diag(exp(muZeta + diag(SigmaB)/2)) %*% t(GT)
  SigmaBY = SigmaB %*% SigmaYMod
  SigmaSY = SigmaSB %*% SigmaYMod
  subDistn = getSubsidenceVarianceMat(params, fault = csz, G = G)
  SigmaY = subDistn$Sigma
  
  # now block them into predictions and data covariance matrices
  SigmaP = cbind(rbind(SigmaB, SigmaSB), rbind(t(SigmaSB), SigmaS))
  SigmaD = cbind(rbind(SigmaS + diag(sigmaXi^2), t(SigmaSY)), rbind(SigmaSY, SigmaY))
  SigmaPtoD = cbind(rbind(t(SigmaSB), SigmaS), rbind(SigmaBY, SigmaSY))
  
  # get block means
  muP = muZeta
  muD = c(rep(muZeta + muXi, length(logX)), GT %*% exp(muZeta + diag(SigmaB)/2))
  
  # compute conditional normal mean and standard error in mean estimate
  Xd = c(logX, Y)
  condDistn = conditionalNormal(Xd, muP, muD, SigmaP, SigmaD, SigmaPtoD)
  muc = condDistn$muc
  Sigmac = condDistn$Sigmac
  
  ##### now generate the conditional simulations.  Note that we still use the
  ##### marginal covariance structure, but we've conditionally updated the mean.
  # get prediction locations
  # get CSZ prediction coordinates
  xd = cbind(slipDatCSZ$lon, slipDatCSZ$lat)  # d is GPS locations
  xp = cbind(csz$longitude, csz$latitude) # p is fault areal locations
  nd = nrow(xd)
  np = nrow(xp)
  areal = 1:np
  point = (np+1):(np+nd)
  
  # Cholesky decomp used for simulations (don't add Sigmac because that 
  # can lead to non-physical results.  Here we just treat the mean as 
  # constant and use the marginal variability)
  # NOTE: deflate variance to account for increased variation due to
  #       conditional mean (using MSE = Var + bias^2 formula)
  # SigmaPred = SigmaP + Sigmac
  # varDeflation = 1 - mean((muc[point] - muZeta)^2)/sigmaZeta^2
  # SigmaPred = SigmaP * varDeflation
  # SigmaPred = SigmaP
  SigmaPred = Sigmac
  SigmaL = t(chol(SigmaPred))
  
  # generate predictive simulations
  zSims = matrix(rnorm(nsim*(nrow(xp)+nrow(xd))), nrow=nrow(xp)+nrow(xd), ncol=nsim)
  logZetaSims0 = SigmaL %*% zSims # each column is a zero mean simulation
  logZetaSims = sweep(logZetaSims0, 1, muc, "+") # add conditional mean to simulations
  zetaSims = exp(logZetaSims)
  tvec = taper(c(csz$depth, slipDatCSZ$Depth), lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # seperate areal average sims from GPS point location sims
  slipSimsGPS = slipSims[point,]
  slipSims = slipSims[areal,]
  
  # get mean slip prediction field
  meanSlip = exp(muc[areal] + diag(SigmaPred[areal,areal])/2) * tvec[areal]
  meanSlipGPS = exp(muc[point] + diag(SigmaPred[point,point])/2) * tvec[point]
  
  return(list(meanSlip=meanSlip, meanSlipGPS=meanSlipGPS, slipSims=slipSims, slipSimsGPS=slipSimsGPS, 
              Sigmac=Sigmac[areal, areal], muc=muc[areal], SigmacGPS = Sigmac[point, point], 
              mucGPS=muc[point], Sigma=SigmaPred[areal,areal], SigmaGPS=SigmaPred[point,point]))
}

##### Try importance sampling for approximate predictive distribution.  Sample from 
##### log zeta(B:) conditional on GPS and subsidence data under a MVN approximation 
##### to the true predictive distribution, and weight by
##### f(X | log zeta(B:)) f(Y | log zeta(B:)) f(zeta) / [f(X) f(Y) f*(zeta)], 
##### where f*(zeta) is the MVN approximation to the predictive distribution.  
##### Once all samples are taken, the probability of accepting the sample is the 
##### weight divided by the max weight of all the samples.  Note that the hope is 
##### that the distribution sampled from has a decent portion of the density 
##### aligning with the high density regions of the weighting function.  Otherwise 
##### the estimate will be very poor.
genFullPreds = function(params, nsim=1000, normalizeTaper=FALSE, dStar=28000) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # get Okada linear transformation matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda0)
  
  # get taper vector
  tvec = taper(csz$depth, lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  
  # generate simulations of zeta marginally
  zeta = predsArealAndLoc(params, nsim)
  zetaGPS = zeta$slipSimsGPS
  zetaAreal = zeta$slipSims
  
  # convert to subsidence simulations (compute G %*% T %*% zeta)
  subs = G %*% sweep(zetaAreal, 1, tvec, "*")
  
  ## compute f(Y | zeta) for each simulation of zeta
  eps = dr1$Uncertainty
  diffs = sweep(subs, 1, -dr1$subsidence, "-")
  logFYGivenZeta = apply(dnorm(diffs, 0, eps, log = TRUE), 2, sum)
  
  ## compute f(X | zeta) for each simulation of zeta
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  logX = log(slipDatCSZ$slip)
  sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  diffs = sweep(zetaGPS, 1, logX, "-")
  logFXGivenZeta = apply(dnorm(diffs, 0, sigmaXi, log = TRUE), 2, sum)
  
  ## compute f(Y)
  # for now, we'll assume normality of Y even though that's not a great 
  # assumption.
  arealCSZCor = getArealCorMat(fault)
  SigmaZeta = arealCSZCor * sigmaZeta^2
  moments = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G, tvec, TRUE, csz)
  muEst = moments$mu
  SigmaEstU = chol(moments$Sigma + diag(dr1$Uncertainty^2))
  logFY = logLikGP(-dr1$subsidence - muEst, SigmaEstU)
  
  ## compute f(X)
  logXCntr = logX - muXi - muZeta
  coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                             onlyUpper=TRUE, smoothness=nuZeta, 
                             Distance="rdist.earth", Dist.args=list(miles=FALSE))
  SigmaX = corMatGPS * sigmaZeta^2 + diag(sigmaXi^2)
  logFX = logLikGP(logXCntr, chol(SigmaX))
  
  # get sample weights
  logWeights = logFXGivenZeta + logFYGivenZeta - logFY - logFX
  weights = exp(logWeights)
  
  # compute weighted samples
  weightedZeta = exp(sweep(log(zeta), 2, logWeights, "+"))
  
  # estimate mean, standard deviation, standard error
  muEst = apply(weightedZeta, 1, mean)
  sdEst = apply(weightedZeta, 1, sd)
  seEst = sdEst/nsim
  
  return(list(muEst=muEst, sdEst=sdEst, seEst=seEst, zetaGivenXSims=zetaGivenX, weightedZetaSims=weightedZetaSims))
}

##### Try importance sampling for full predictive distribution.  Sample from 
##### log zeta(B:) | X and weight by f(Y | log zeta(B:)) / f(Y).  Alternatively, 
##### one could sample from log zeta(B:) marginally and weight by
##### f(X | log zeta(B:)) f(Y | log zeta(B:)) / [f(X) f(Y)].  Once all samples 
##### are taken, the probability of accepting the sample is the weight divided 
##### by the max weight of all the samples.  Note that the hope is that the 
##### distribution sampled from is heavy tailed, according the rejection sampling 
##### algorithm.
genFullPredsGPS = function(params, nsim=25000, normalizeTaper=FALSE, dStar=28000) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # get Okada linear transformation matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda0)
  
  # get taper vector
  tvec = taper(csz$depth, lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  
  # generate simulations of zeta given the GPS data
  zetaGivenX = predsGivenGPS(params, nsim)$slipSims
  
  # convert to subsidence simulations (compute G %*% T %*% zeta | X)
  subsGivenX = G %*% sweep(zetaGivenX, 1, tvec, "*")
  
  # compute f(Y | zeta) for each simulation of zetaGivenX
  eps = dr1$Uncertainty
  diffs = sweep(subsGivenX, 1, -dr1$subsidence, "-")
  logFYGivenZeta = apply(dnorm(diffs, 0, eps, log = TRUE), 2, sum)
  
  # compute f(Y)
  # for now, we'll assume normality of Y even though that's not a great 
  # assumption.
  arealCSZCor = getArealCorMat(fault)
  SigmaZeta = arealCSZCor * sigmaZeta^2
  moments = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G, tvec, TRUE, csz)
  muEst = moments$mu
  SigmaEstU = chol(moments$Sigma + diag(dr1$Uncertainty^2))
  logFY = logLikGP(-dr1$subsidence - muEst, SigmaEstU)
  
  # get sample weights
  logWeights = logFYGivenZeta - logFY
  weights = exp(logWeights)
  
  # compute weighted samples
  weightedZeta = exp(sweep(log(zetaGivenX), 2, logWeights, "+"))
  
  # estimate mean, standard deviation, standard error
  muEst = apply(weightedZeta, 1, mean)
  sdEst = apply(weightedZeta, 1, sd)
  seEst = sdEst/nsim
  
  return(list(muEst=muEst, sdEst=sdEst, seEst=seEst, zetaGivenXSims=zetaGivenX, weightedZetaSims=weightedZetaSims))
}
# max weight so far: -1253.168

##### Try importance sampling for full predictive distribution.  Sample from 
##### log zeta(B:) marginally and weight by
##### f(X | log zeta(B:)) f(Y | log zeta(B:)) / [f(X) f(Y)].  Once all samples 
##### are taken, the probability of accepting the sample is the weight divided 
##### by the max weight of all the samples.  Note that the hope is that the 
##### distribution sampled from has a decent portion of the density aligning with 
##### the high density regions of the weighting function.  Ohterwise, the 
##### estimate will be very poor.
genFullPredsMarginal = function(params, nsim=1000) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # get Okada linear transformation matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda0)
  
  # get taper vector
  tvec = taper(csz$depth, lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  
  # generate simulations of zeta marginally
  zeta = predsArealAndLoc(params, nsim)
  zetaGPS = zeta$slipSimsGPS
  zetaAreal = zeta$slipSims
  
  # convert to subsidence simulations (compute G %*% T %*% zeta)
  subs = G %*% sweep(zetaAreal, 1, tvec, "*")
  
  ## compute f(Y | zeta) for each simulation of zeta
  eps = dr1$Uncertainty
  diffs = sweep(subs, 1, -dr1$subsidence, "-")
  logFYGivenZeta = apply(dnorm(diffs, 0, eps, log = TRUE), 2, sum)
  
  ## compute f(X | zeta) for each simulation of zeta
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  logX = log(slipDatCSZ$slip)
  sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  diffs = sweep(zetaGPS, 1, logX, "-")
  logFXGivenZeta = apply(dnorm(diffs, 0, sigmaXi, log = TRUE), 2, sum)
  
  ## compute f(Y)
  # for now, we'll assume normality of Y even though that's not a great 
  # assumption.
  arealCSZCor = getArealCorMat(fault)
  SigmaZeta = arealCSZCor * sigmaZeta^2
  moments = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G, tvec, TRUE, csz)
  muEst = moments$mu
  SigmaEstU = chol(moments$Sigma + diag(dr1$Uncertainty^2))
  logFY = logLikGP(-dr1$subsidence - muEst, SigmaEstU)
  
  ## compute f(X)
  logXCntr = logX - muXi - muZeta
  coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                             onlyUpper=TRUE, smoothness=nuZeta, 
                             Distance="rdist.earth", Dist.args=list(miles=FALSE))
  SigmaX = corMatGPS * sigmaZeta + diag(sigmaXi)
  logFX = logLikGP(logXCntr, chol(SigmaX))
  
  # get sample weights
  logWeights = logFXGivenZeta + logFYGivenZeta - logFY - logFX
  weights = exp(logWeights)
  
  # compute weighted samples
  weightedZeta = exp(sweep(log(zeta), 2, logWeights, "+"))
  
  # estimate mean, standard deviation, standard error
  muEst = apply(weightedZeta, 1, mean)
  sdEst = apply(weightedZeta, 1, sd)
  seEst = sdEst/nsim
  
  return(list(muEst=muEst, sdEst=sdEst, seEst=seEst, zetaGivenXSims=zetaGivenX, weightedZetaSims=weightedZetaSims))
}

# function for computing predictive distribution given subsidence data.  Note that 
# the subsidence data should only be the data from one earthquake.  The returned 
# values relating to beta correspond to log zeta.  For instance, betaEsts is the 
# vector of estimates of log zeta for the given earthquake.  If prior=TRUE, then 
# Goldfinger 2012 data is used to create a prior for earthquake magnitude.
# NOTE: this function should only be called on data for a single fixed earthquake.  
predsGivenSubsidence = function(params, muVec=params[2], fault=csz, subDat=dr1, niter=500, estVarMat=TRUE, 
                                G=NULL, prior=FALSE, normalizeTaper=FALSE, dStar=28000, 
                                tvec=taper(getFaultCenters(fault)[,3], lambda=params[1], 
                                           normalize=normalizeTaper, dStar=dStar), priorMaxSlip=NA, 
                                normalModel=FALSE, posNormalModel=FALSE, taperedGPSDat=FALSE, 
                                gpsDat=slipDatCSZ, fastPNSim=TRUE) {
  if(posNormalModel && !normalModel)
    stop("normalModel must be set to TRUE if posNormalModel is set to TRUE")
  
  # first rows correspond to GPS data means, last rows correspond to csz 
  # grid cell means
  if(length(muVec) == 1)
    muVec = rep(muVec, nrow(gpsDat)+nrow(fault))
  
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  if(!taperedGPSDat) {
    corPar = getCorPar(normalModel=normalModel)
    phiZeta = corPar$phiZeta
    nuZeta = corPar$nuZeta
  }
  else {
    phiZeta = params[length(params)]
    nuZeta = 3/2
  }
  
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
  
  # compute areal log zeta covariance matrix
  if(!taperedGPSDat) {
    if(identical(fault, csz)) {
      arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
      SigmaArea = arealCSZCor * sigmaZeta^2
    }
    else {
      SigmaArea = arealZetaCov(params, fault, nDown1=9, nStrike1=12, normalModel=normalModel)
    }
    
    # compute point log zeta covariance matrix and cross-covariance
    xs = cbind(gpsDat$lon, gpsDat$lat)
    SigmaPointToArea = pointArealZetaCov(params, xs, fault, nDown=9, nStrike=12, normalModel=normalModel)
    SigmaPoint = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                                smoothness=nuZeta, Distance="rdist.earth", 
                                Dist.args=list(miles=FALSE)) * sigmaZeta^2
  }
  else {
    # under tapered GPS data model, don't compute the areal covariance matrix.  Approximate it 
    # with the pointwise covariance matrix
    
    xs = cbind(fault$longitude, fault$latitude)
    SigmaArea = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                               smoothness=nuZeta, Distance="rdist.earth", 
                               Dist.args=list(miles=FALSE)) * sigmaZeta^2
    
    xsPt = cbind(gpsDat$lon, gpsDat$lat)
    SigmaPointToArea = stationary.cov(x1=xsPt, x2=xs, Covariance="Matern", theta=phiZeta, 
                                      smoothness=nuZeta, Distance="rdist.earth", 
                                      Dist.args=list(miles=FALSE)) * sigmaZeta^2
    SigmaPoint = stationary.cov(xsPt, Covariance="Matern", theta=phiZeta, 
                                smoothness=nuZeta, Distance="rdist.earth", 
                                Dist.args=list(miles=FALSE)) * sigmaZeta^2
  }
  
  # combine covariance matrices:
  SigmaZeta = cbind(rbind(SigmaPoint, t(SigmaPointToArea)), 
                    rbind(SigmaPointToArea, SigmaArea))
  
  # set up required data variables in the R environment for Stan
  priorSigmaL = t(chol(SigmaZeta))
  priorMu = as.vector(muVec)
  X = G %*% diag(tvec)
  y = as.vector(-subDat$subsidence)
  sigmaY = subDat$Uncertainty
  n = length(y)
  pAreal = nrow(fault)
  pPoint = nrow(xs)
  
  if(!normalModel) {
    # load rstan and precompile parallel MC code for faster fitting in parallel.  Also use O3 optimization
    set_cppo(mode = "fast")
    rstan_options(auto_write = TRUE)
    numCores = parallel::detectCores()
    options(mc.cores = numCores)
  }
  
  if(normalModel || posNormalModel) {
    # need:
    # beta (logZeta areal), zeta (areal), logZetaPoint, seismicMoment, Mw
    # under the normal model, predictions becomes a conditional normal problem
    
    # get the distribution of Y and its covariance to zeta process
    mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaArea, G, subDat=subDat, 
                                     tvec=tvec, normalModel=normalModel)
    SigmaY = mvnApprox$Sigma
    diag(SigmaY) = diag(SigmaY) + sigmaY^2
    muY = X %*% muVec[1:nrow(fault)]
    SigmaZetaToY = rbind(SigmaPointToArea %*% t(X), 
                          SigmaArea %*% t(X))
    
    # compute the conditional distribution of zeta given Y
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
          else {
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
    zetaPoint = zetaSims[1:nrow(gpsDat),]
    zetaAreal = zetaSims[-(1:nrow(gpsDat)),]
    logZetaPoint = logZetaSims[1:nrow(gpsDat),]
    logZetaAreal = logZetaSims[-(1:nrow(gpsDat)),]
    
    # from zetaSims, generate seismic moments and magnitudes
    slipSims = sweep(zetaAreal, 1, tvec, "*")
    mags = apply(slipSims, 2, getMomentFromSlip, fault=fault)
    moments = 10^(mags*1.5 + 9.05)
    
    # compile results into predResults object with:
    # beta (logZeta areal), zeta (areal), logZetaPoint, seismicMoment, Mw (include zetaPoint as well since logZetaPoint is NA)
    predResults = list(beta=logZetaAreal, zeta=zetaAreal, logZetaPoint=logZetaPoint, zetaPoint=zetaPoint, 
                       seismicMoment=moments, Mw=mags)
  }
  else if(!prior && is.na(priorMaxSlip)) {
    # fit the regression model with Stan
    # don't use a prior on the EQ magnitude
    areas = fault$length * fault$width
    mu = 4e10 # rigidty (value as assumed by Randy)
    
    predResults <- stan(file = "predsGivenSub.stan", iter=niter, 
                        data=list(n=n, pAreal=pAreal, pPoint=pPoint, X=X, y=y, 
                                  sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL, 
                                  tvec=tvec, areas=areas, rigidity=mu), 
                        pars=c("alpha"), include=FALSE)
  }
  else if(is.na(priorMaxSlip)) {
    ## use prior on EQ magnitude from Goldfinger 2012 data
    # read in Goldfinger data
    goldfinger = read.csv("goldfinger2012Table8.csv", header=TRUE)
    dr1Events = c(1:10, 12:17, 19:20, 22, 24, 28:29)
    dr1Events = 1:31 # only keep events after and including T13
    goldfinger = goldfinger[dr1Events,]
    M0 = goldfinger$seismicMoment
    
    # convert from dyne cm to N m
    M0 = M0 * 10^(-7)
    
    # get exponential rate parameter for fault grid cells
    # NOTE: this is before tapering, since zeta itself is not tapered
    areas = fault$length * fault$width
    mu = 4e10 # rigidty (value as assumed by Randy)
    tau = 1/mean(M0, na.rm=TRUE)
    
    predResults <- stan(file = "predsGivenSubExp.stan", iter=niter, 
                        data=list(n=n, pAreal=pAreal, pPoint=pPoint, X=X, y=y, 
                                  sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL, 
                                  tvec=tvec, areas=areas, rigidity=mu, tau=tau), 
                        pars=c("alpha"), include=FALSE)
  }
  else {
    # priorMaxSlip is the 95th percentile of our exponential distribution (maybe 500m by default)
    matchFun = function(tau) {
      abs(qexp(.95, rate=tau) - priorMaxSlip)
    }
    tau = optimize(matchFun, c(0, 1))$minimum
    
    # get other necessary values
    areas = fault$length * fault$width
    mu = 4e10 # rigidty (value as assumed by Randy)
    
    predResults <- stan(file = "predsGivenSubMaxExp.stan", iter=niter, 
                        data=list(n=n, pAreal=pAreal, pPoint=pPoint, X=X, y=y, 
                                  sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL, 
                                  tvec=tvec, areas=areas, rigidity=mu, tau=tau), 
                        pars=c("alpha"), include=FALSE)
  }
#   # show results
#   print(predResults, digits = 1)
  
  # extract the samples of the relevant parameters
  if(!normalModel)
    tab <- extract(predResults, permuted = TRUE) # return a list of arrays 
  else
    tab = predResults
  betaTab = tab$beta
  zetaTab = tab$zeta
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = exp(logZetaPointTab)
  if(normalModel)
    zetaPointTab = tab$zetaPoint
  M0Tab = tab$seismicMoment
  MwTab = tab$Mw
  
  # compute means, standard deviations, and middle 95% confidence intervals for relevant values
  if(normalModel && !posNormalModel) {
    betaEsts = NA
    betaSD = NA
    betaMed = NA
    beta025 = NA
    beta975 = NA
    zetaEsts = muc[-(1:nrow(gpsDat))]
    zetaSD = sqrt(diag(Sigmac)[-(1:nrow(gpsDat))])
    zetaMed = zetaEsts
    zeta025 = qnorm(.025, zetaEsts, zetaSD)
    zeta975 = qnorm(.975, zetaEsts, zetaSD)
    logZetaPointEsts = NA
    logZetaPointSD = NA
    logZetaPointMed = NA
    logZetaPoint025 = NA
    logZetaPoint975 = NA
    zetaPointEsts = muc[1:nrow(gpsDat)]
    zetaPointSD = sqrt(diag(Sigmac)[1:nrow(gpsDat)])
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
  }
  else{
    if(!normalModel && !posNormalModel){
      betaTab = t(betaTab)
      zetaTab = t(zetaTab)
      logZetaPointTab = t(logZetaPointTab)
      zetaPointTab = t(zetaPointTab)
    }
    betaEsts = rowMeans(betaTab)
    betaSD = apply(betaTab, 1, sd)
    betaMed = apply(betaTab, 1, median)
    beta025 = apply(betaTab, 1, quantile, probs=0.025)
    beta975 = apply(betaTab, 1, quantile, probs=0.975)
    zetaEsts = rowMeans(zetaTab)
    zetaSD = apply(zetaTab, 1, sd)
    zetaMed = apply(zetaTab, 1, median)
    zeta025 = apply(zetaTab, 1, quantile, probs=0.025)
    zeta975 = apply(zetaTab, 1, quantile, probs=0.975)
    logZetaPointEsts = rowMeans(logZetaPointTab)
    logZetaPointSD = apply(logZetaPointTab, 1, sd)
    logZetaPointMed = apply(logZetaPointTab, 1, median)
    logZetaPoint025 = apply(logZetaPointTab, 1, quantile, probs=0.025)
    logZetaPoint975 = apply(logZetaPointTab, 1, quantile, probs=0.975)
    zetaPointEsts = rowMeans(zetaPointTab)
    zetaPointSD = apply(zetaPointTab, 1, sd)
    zetaPointMed = apply(zetaPointTab, 1, median)
    zetaPoint025 = apply(zetaPointTab, 1, quantile, probs=0.025)
    zetaPoint975 = apply(zetaPointTab, 1, quantile, probs=0.975)
    M0Est = mean(M0Tab)
    M0SD = sd(M0Tab)
    M0Med = median(M0Tab)
    M0975 = quantile(probs=.975, M0Tab)
    M0025 = quantile(probs=.025, M0Tab)
    MwEst = mean(MwTab)
    MwSD = sd(MwTab)
    MwMed = median(MwTab)
    Mw975 = quantile(probs=.975, MwTab)
    Mw025 = quantile(probs=.025, MwTab)
  }
  
  # estimate covariance matrix of params using emprical estimate from MCMC samples if necessary
  if(!estVarMat)
    return(list(predResults=predResults, resultTab=tab, 
                betaEsts=betaEsts, betaSD=betaSD, betaMed=betaMed, beta025=beta025, beta975=beta975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zetaMed=zetaMed, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPointMed=logZetaPointMed, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPointMed=zetaPointMed, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaPointVarMat=NULL, zetaPointVarMat=NULL, 
                M0Est=M0Est, M0SD=M0SD, M0Med=M0Med, M0975=M0975, M0025=M0025, MwEst=MwEst, MwSD=MwSD, Mw975=Mw975, Mw025=Mw025))
  else {
    if(normalModel) {
      # for a normal model we have the exact variance matrix of the estimates
      betaVarMat = NA
      zetaVarMat = Sigmac[-(1:nrow(gpsDat)), -(1:nrow(gpsDat))]
      logZetaPointVarMat = NA
      zetaPointVarMat = Sigmac[1:nrow(gpsDat), 1:nrow(gpsDat)]
    }
    else {
      # for a non-normal model we estimate the empirical variance matrix
      betaTabCntr = sweep(betaTab, 2, betaEsts, "-")
      zetaTabCntr = sweep(zetaTab, 2, zetaEsts, "-")
      logZetaPointTabCntr = sweep(logZetaPointTab, 2, logZetaPointEsts, "-")
      zetaPointTabCntr = sweep(zetaPointTab, 2, zetaPointEsts, "-")
      betaVarMat = t(betaTabCntr) %*% betaTabCntr/nrow(betaTabCntr)
      zetaVarMat = t(zetaTabCntr) %*% zetaTabCntr/nrow(zetaTabCntr)
      logZetaPointVarMat = t(logZetaPointTabCntr) %*% logZetaPointTabCntr/nrow(logZetaPointTabCntr)
      zetaPointVarMat = t(zetaPointTabCntr) %*% zetaPointTabCntr/nrow(zetaPointTabCntr)
    }
    return(list(predResults=predResults, resultTab=tab, 
                betaEsts=betaEsts, betaSD=betaSD, betaMed=betaMed, beta025=beta025, beta975=beta975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zetaMed=zetaMed, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPointMed=logZetaPointMed, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPointMed=zetaPointMed, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                betaVarMat=betaVarMat, zetaVarMat=zetaVarMat, 
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat, 
                M0Est=M0Est, M0SD=M0SD, M0Med=M0Med, M0975=M0975, M0025=M0025, MwEst=MwEst, MwSD=MwSD, Mw975=Mw975, Mw025=Mw025))
  }
}

# given the output from stan, computes mean, sd, and middle 95% interval estimators along with
# (if the user requests) variance/covariance matrices for the means
extractStanPredictions = function(stanResults, estVarMat=TRUE) {
  # extract the samples of the relevant parameters
  tab <- extract(stanResults, permuted = TRUE) # return a list of arrays 
  betaTab = tab$beta
  zetaTab = tab$zeta
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = exp(logZetaPointTab)
  M0Tab = tab$seismicMoment
  MwTab = tab$Mw
  
  # compute means, standard deviations, and middle 95% confidence intervals for relevant values
  betaEsts = colMeans(betaTab)
  betaSD = apply(betaTab, 2, sd)
  beta025 = apply(betaTab, 2, quantile, probs=0.025)
  beta975 = apply(betaTab, 2, quantile, probs=0.975)
  zetaEsts = colMeans(zetaTab)
  zetaSD = apply(zetaTab, 2, sd)
  zeta025 = apply(zetaTab, 2, quantile, probs=0.025)
  zeta975 = apply(zetaTab, 2, quantile, probs=0.975)
  logZetaPointEsts = colMeans(logZetaPointTab)
  logZetaPointSD = apply(logZetaPointTab, 2, sd)
  logZetaPoint025 = apply(logZetaPointTab, 2, quantile, probs=0.025)
  logZetaPoint975 = apply(logZetaPointTab, 2, quantile, probs=0.975)
  zetaPointEsts = colMeans(zetaPointTab)
  zetaPointSD = apply(zetaPointTab, 2, sd)
  zetaPoint025 = apply(zetaPointTab, 2, quantile, probs=0.025)
  zetaPoint975 = apply(zetaPointTab, 2, quantile, probs=0.975)
  M0Est = mean(M0Tab)
  M0SD = sd(M0Tab)
  M0975 = quantile(probs=.975, M0Tab)
  M0025 = quantile(probs=.025, M0Tab)
  MwEst = mean(MwTab)
  MwSD = sd(MwTab)
  Mw975 = quantile(probs=.975, MwTab)
  Mw025 = quantile(probs=.025, MwTab)
  
  # estimate covariance matrix of params using emprical estimate from MCMC samples if necessary
  if(!estVarMat)
    return(list(predResults, 
                betaEsts, betaSD, beta025, beta975, 
                zetaEsts, zetaSD, zeta025, zeta975, 
                logZetaPointEsts, logZetaPointSD, logZetaPoint025, logZetaPoint975, 
                zetaPointEsts, zetaPointSD, zetaPoint025, zetaPoint975, 
                M0Est, M0SD, M0975, M0025, MwEst, MwSD, Mw975, Mw025))
  else {
    betaTabCntr = sweep(betaTab, 2, betaEsts, "-")
    zetaTabCntr = sweep(zetaTab, 2, zetaEsts, "-")
    logZetaPointTabCntr = sweep(logZetaPointTab, 2, logZetaPointEsts, "-")
    zetaPointTabCntr = sweep(zetaPointTab, 2, zetaPointEsts, "-")
    betaVarMat = t(betaTabCntr) %*% betaTabCntr/nrow(betaTabCntr)
    zetaVarMat = t(zetaTabCntr) %*% zetaTabCntr/nrow(zetaTabCntr)
    logZetaPointVarMat = t(logZetaPointTabCntr) %*% logZetaPointTabCntr/nrow(logZetaPointTabCntr)
    zetaPointVarMat = t(zetaPointTabCntr) %*% zetaPointTabCntr/nrow(zetaPointTabCntr)
    return(list(betaEsts=betaEsts, betaSD=betaSD, beta025=beta025, beta975=beta975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                betaVarMat=betaVarMat, zetaVarMat=zetaVarMat, 
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat, 
                M0Est, M0SD, M0975, M0025, MwEst, MwSD, Mw975, Mw025))
  }
}

# this function updates the estimate of the mean vector of log zeta by weighting the 
# estimates of each earthquake by their variances.
updateMu = function(params, muVec=params[2], fault=csz, niter=1000, gpsDat=slipDatCSZ, G=NULL, 
                    usePrior=FALSE, subDat=dr1, normalizeTaper=FALSE, dStar=28000, 
                    tvec=taper(getFaultCenters(fault)[,3], lambda=params[1], dStar=dStar, 
                               normalize=normalizeTaper)) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  if(length(muVec) == 1)
    muVec = rep(muVec, nrow(gpsDat) + nrow(fault))
  
  # precompute G if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # Only use events that have at least 5 observations (leaves out T12, T11, and T9a)
  minObs = 5
  numObs = table(subDat$event)
  threshUniqueEvents = uniqueEvents[numObs >= minObs]
  
  sampleSize = length(threshUniqueEvents)+1
  ## for each event, get estimates for mean of log zeta and covariance matrices
  ## (for areal averages and point estimates over GPS data support)
  muMat = matrix(nrow=nrow(fault), ncol=sampleSize)
  sdMat = matrix(nrow=nrow(fault), ncol=sampleSize)
  muMatPoint = matrix(nrow=nrow(gpsDat), ncol=sampleSize)
  sdMatPoint = matrix(nrow=nrow(gpsDat), ncol=sampleSize)
  Sigmas = list()
  SigmasPoint = list()
  W = matrix(nrow=nrow(fault), ncol=sampleSize)
  WPoint = matrix(nrow=nrow(gpsDat), ncol=sampleSize)
  allStanResults = list()
  
  # first get estimates of log zeta from subsidence data attributed to each event
  for(e in 1:length(threshUniqueEvents)) {
    print(paste0("estimating log zeta for event ", threshUniqueEvents[e]))
    
    # all get event MCMC data
    eventInds = events == threshUniqueEvents[e]
    thisSubDat = subDat[eventInds,]
    thisG = G[eventInds,]
    eventPreds = predsGivenSubsidence(params, muVec=muVec, fault=fault, subDat=thisSubDat, niter=niter, 
                                      G=thisG, prior=usePrior, tvec=tvec)
    # areal values
    muMat[,e] = eventPreds$betaEsts
    sdMat[,e] = eventPreds$betaSD
    Sigmas[[e]] = eventPreds$betaVarMat
    W[,e] = (eventPreds$betaSD)^(-1)
    # pointwise values
    muMatPoint[,e] = eventPreds$logZetaPointEsts
    sdMatPoint[,e] = eventPreds$logZetaPointSD
    SigmasPoint[[e]] = eventPreds$logZetaPointVarMat
    WPoint[,e] = (eventPreds$logZetaPointSD)^(-1)
    # full Stan MCMC dataset
    allStanResults[[e]] = eventPreds$predResults
  }
  
  # now get predictions for log zeta given GPS data
  print("estimating log zeta for GPS data")
  eventPreds = predsGivenGPSFull(params, nsim=2, muVec=muVec)
  # areal estimates
  muMat[,e+1] = eventPreds$muc
  sdMat[,e+1] = diag(eventPreds$Sigmac)
  Sigmas[[e+1]] = eventPreds$Sigmac
  W[,e+1] = diag(eventPreds$Sigmac)^(-1)
  # pointwise estimates
  muMatPoint[,e+1] = eventPreds$mucGPS
  sdMatPoint[,e+1] = diag(eventPreds$SigmacGPS)
  SigmasPoint[[e+1]] = eventPreds$SigmacGPS
  WPoint[,e+1] = diag(eventPreds$SigmacGPS)^(-1)
  
  # convert W matrix from related weights to absolute weights (summing to 1)
  W = sweep(W, 1, rowSums(W), "/")
  WPoint = sweep(WPoint, 1, rowSums(WPoint), "/")
  
  ## calculate new mean estimate and its variance matrix
  newMu = rowSums(W*muMat)
  varMat = matrix(0, nrow=nrow(fault), ncol=nrow(fault))
  newMuPoint = rowSums(WPoint*muMatPoint)
  varMatPoint = matrix(0, nrow=nrow(gpsDat), ncol=nrow(gpsDat))
  for(k in 1:sampleSize) {
    # areal values
    thisVar = Sigmas[[k]]
    thisWeights = W[,k]
    thisVar = sweep(sweep(thisVar, 1, thisWeights, "*"), 2, thisWeights, "*")
    varMat = varMat + thisVar
    
    #pointwise values
    thisVarPoint = SigmasPoint[[k]]
    thisWeightsPoint = WPoint[,k]
    thisVarPoint = sweep(sweep(thisVarPoint, 1, thisWeightsPoint, "*"), 2, thisWeightsPoint, "*")
    varMatPoint = varMatPoint + thisVarPoint
  }
  
  ## calculate new mean estimate based only on subsidence and its variance matrix
  WSub = W[,1:(sampleSize-1)]
  WSub = sweep(WSub, 1, rowSums(WSub), "/")
  WPointSub = WPoint[,1:(sampleSize-1)]
  WPointSub = sweep(WPointSub, 1, rowSums(WPointSub), "/")
  newMuSub = rowSums(WSub*muMat[,1:(sampleSize-1)])
  varMatSub = matrix(0, nrow=nrow(fault), ncol=nrow(fault))
  newMuSubPoint = rowSums(WPointSub*muMatPoint[,1:(sampleSize-1)])
  varMatSubPoint = matrix(0, nrow=nrow(gpsDat), ncol=nrow(gpsDat))
  for(k in 1:(sampleSize-1)) {
    #areal values
    thisVar = Sigmas[[k]]
    thisWeights = WSub[,k]
    thisVar = sweep(sweep(thisVar, 1, thisWeights, "*"), 2, thisWeights, "*")
    varMatSub = varMatSub + thisVar
    # pointwise values
    thisVarPoint = SigmasPoint[[k]]
    thisWeightsPoint = WPointSub[,k]
    thisVarPoint = sweep(sweep(thisVarPoint, 1, thisWeightsPoint, "*"), 2, thisWeightsPoint, "*")
    varMatSubPoint = varMatSubPoint + thisVarPoint
  }
  
  return(list(newMu=newMu, varMat=varMat, muMat=muMat, Sigmas=Sigmas, 
              newMuSub=newMuSub, newMuPoint=newMuPoint, 
              varMatPoint=varMatPoint, muMatPoint=muMatPoint, SigmasPoint=SigmasPoint, 
              newMuSubPoint=newMuSubPoint, allStanResults=allStanResults))
}

updateMuGivenStan = function(params, stanResults, muVec=params[2], fault=csz, gpsDat=slipDatCSZ, 
                             G=NULL, subDat=dr1) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  if(length(muVec) == 1)
    muVec = rep(muVec, nrow(gpsDat) + nrow(fault))
  
  # precompute G
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # Only use events that have at least 5 observations (leaves out T12, T11, and T9a)
  minObs = 5
  numObs = table(subDat$event)
  threshUniqueEvents = uniqueEvents[numObs >= minObs]
  
  sampleSize = length(threshUniqueEvents)+1
  ## for each event, get estimates for mean of log zeta and covariance matrices
  ## (for areal averages and point estimates over GPS data support)
  muMat = matrix(nrow=nrow(fault), ncol=sampleSize)
  sdMat = matrix(nrow=nrow(fault), ncol=sampleSize)
  muMatPoint = matrix(nrow=nrow(gpsDat), ncol=sampleSize)
  sdMatPoint = matrix(nrow=nrow(gpsDat), ncol=sampleSize)
  Sigmas = list()
  SigmasPoint = list()
  W = matrix(nrow=nrow(fault), ncol=sampleSize)
  WPoint = matrix(nrow=nrow(gpsDat), ncol=sampleSize)
  allStanResults = list()
  
  # first get estimates of log zeta from subsidence data attributed to each event
  for(e in 1:length(threshUniqueEvents)) {
    print(paste0("estimating log zeta for event ", threshUniqueEvents[e]))
    
    # all get event MCMC data
    eventInds = events == threshUniqueEvents[e]
    thisSubDat = subDat[eventInds,]
    thisG = G[eventInds,]
    eventPreds = extractStanPredictions(stanResults[[e]])
    # areal values
    muMat[,e] = eventPreds$betaEsts
    sdMat[,e] = eventPreds$betaSD
    Sigmas[[e]] = eventPreds$betaVarMat
    W[,e] = (eventPreds$betaSD)^(-1)
    # pointwise values
    muMatPoint[,e] = eventPreds$logZetaPointEsts
    sdMatPoint[,e] = eventPreds$logZetaPointSD
    SigmasPoint[[e]] = eventPreds$logZetaPointVarMat
    WPoint[,e] = (eventPreds$logZetaPointSD)^(-1)
    # full Stan MCMC dataset
    allStanResults[[e]] = eventPreds$predResults
  }
  
  # now get predictions for log zeta given GPS data
  print("estimating log zeta for GPS data")
  eventPreds = predsGivenGPSFull(params, nsim=2, muVec=muVec)
  # areal estimates
  muMat[,e+1] = eventPreds$muc
  sdMat[,e+1] = diag(eventPreds$Sigmac)
  Sigmas[[e+1]] = eventPreds$Sigmac
  W[,e+1] = diag(eventPreds$Sigmac)^(-1)
  # pointwise estimates
  muMatPoint[,e+1] = eventPreds$mucGPS
  sdMatPoint[,e+1] = diag(eventPreds$SigmacGPS)
  SigmasPoint[[e+1]] = eventPreds$SigmacGPS
  WPoint[,e+1] = diag(eventPreds$SigmacGPS)^(-1)
  
  # convert W matrix from related weights to absolute weights (summing to 1)
  W = sweep(W, 1, rowSums(W), "/")
  WPoint = sweep(WPoint, 1, rowSums(WPoint), "/")
  
  ## calculate new mean estimate and its variance matrix
  newMu = rowSums(W*muMat)
  varMat = matrix(0, nrow=nrow(fault), ncol=nrow(fault))
  newMuPoint = rowSums(WPoint*muMatPoint)
  varMatPoint = matrix(0, nrow=nrow(gpsDat), ncol=nrow(gpsDat))
  for(k in 1:sampleSize) {
    # areal values
    thisVar = Sigmas[[k]]
    thisWeights = W[,k]
    thisVar = sweep(sweep(thisVar, 1, thisWeights, "*"), 2, thisWeights, "*")
    varMat = varMat + thisVar
    
    #pointwise values
    thisVarPoint = SigmasPoint[[k]]
    thisWeightsPoint = WPoint[,k]
    thisVarPoint = sweep(sweep(thisVarPoint, 1, thisWeightsPoint, "*"), 2, thisWeightsPoint, "*")
    varMatPoint = varMatPoint + thisVarPoint
  }
  
  ## calculate new mean estimate based only on subsidence and its variance matrix
  WSub = W[,1:(sampleSize-1)]
  WSub = sweep(WSub, 1, rowSums(WSub), "/")
  WPointSub = WPoint[,1:(sampleSize-1)]
  WPointSub = sweep(WPointSub, 1, rowSums(WPointSub), "/")
  newMuSub = rowSums(WSub*muMat[,1:(sampleSize-1)])
  varMatSub = matrix(0, nrow=nrow(fault), ncol=nrow(fault))
  newMuSubPoint = rowSums(WPointSub*muMatPoint[,1:(sampleSize-1)])
  varMatSubPoint = matrix(0, nrow=nrow(gpsDat), ncol=nrow(gpsDat))
  for(k in 1:(sampleSize-1)) {
    #areal values
    thisVar = Sigmas[[k]]
    thisWeights = WSub[,k]
    thisVar = sweep(sweep(thisVar, 1, thisWeights, "*"), 2, thisWeights, "*")
    varMatSub = varMatSub + thisVar
    # pointwise values
    thisVarPoint = SigmasPoint[[k]]
    thisWeightsPoint = WPointSub[,k]
    thisVarPoint = sweep(sweep(thisVarPoint, 1, thisWeightsPoint, "*"), 2, thisWeightsPoint, "*")
    varMatSubPoint = varMatSubPoint + thisVarPoint
  }
  
  return(list(newMu=newMu, varMat=varMat, muMat=muMat, Sigmas=Sigmas, 
              newMuSub=newMuSub, newMuPoint=newMuPoint, 
              varMatPoint=varMatPoint, muMatPoint=muMatPoint, SigmasPoint=SigmasPoint, 
              newMuSubPoint=newMuSubPoint, allStanResults=allStanResults))
}


###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
# functions for computing predictive distribution given subsidence data using maximum slip and 
# subsidence priors

# Note that the subsidence data should only be the data from one earthquake.  The returned 
# values relating to beta correspond to log zeta.  For instance, betaEsts is the 
# vector of estimates of log zeta for the given earthquake.  If ev="all", uses all the subsidence 
# data.  Otherwise uses the data from the specified event.
# NOTE: this function should only be called on data for a single fixed earthquake.  
predsGivenSubsidence2 = function(params, muVec=params[2], fault=csz, subDat=dr1, niter=500, estVarMat=TRUE, 
                                G=NULL, fauxG=NULL, dStar=28000, normalizeTaper=FALSE, 
                                tvec=taper(getFaultCenters(fault)[,3], lambda=params[1], dStar=dStar, 
                                           normalize=normalizeTaper), 
                                fauxDat=getFauxDat(), ev="T1", FrechetFGumbelT=TRUE, slipThresh=Inf) {
  # first rows correspond to GPS data means, last rows correspond to csz 
  # grid cell means
  if(length(muVec) == 1)
    muVec = rep(muVec, nrow(slipDatCSZ)+nrow(fault))
  
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
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
  if(is.null(fauxG)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    fauxG = okadaAll(fault, lonGrid, latGrid, cbind(fauxDat$Lon, fauxDat$Lat), slip=1, poisson=lambda0)
  }
  
  # compute areal log zeta covariance matrix
  if(identical(fault, csz)) {
    arealCSZCor = getArealCorMat(fault)
    SigmaArea = arealCSZCor * sigmaZeta^2
  }
  else {
    SigmaArea = arealZetaCov(params, fault, nDown1=9, nStrike1=12)
  }
  
  # compute point log zeta covariance matrix and cross-covariance
  xs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  SigmaPointToArea = pointArealZetaCov(params, xs, fault, nDown=9, nStrike=12)
  SigmaPoint = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                              smoothness=nuZeta, Distance="rdist.earth", 
                              Dist.args=list(miles=FALSE)) * sigmaZeta^2
  
  # combine covariance matrices:
  SigmaZeta = cbind(rbind(SigmaPoint, t(SigmaPointToArea)), 
                    rbind(SigmaPointToArea, SigmaArea))
  
  # combine subDat and fauxDat
  tmp = data.frame(list(Lon=subDat$Lon, Lat=subDat$Lat, Uncertainty=subDat$Uncertainty, event=subDat$event, 
                        subsidence=subDat$subsidence))
  fullDat = rbind(tmp, fauxDat)
  
  # get relevant variables for the specified event (and make sure we don't base likelihood on faux data)
  fullG = rbind(G, fauxG)
  if(ev == "all")
    goodEvent = c(rep(TRUE, nrow(G)), rep(FALSE, nrow(fauxG)))
  else
    goodEvent = c((events == ev), rep(FALSE, nrow(fauxG)))
  goodG = fullG[goodEvent,]
  GRest = fullG[!goodEvent,]
  
  # set up required data variables in the R environment for Stan
  priorSigmaL = t(chol(SigmaZeta))
  priorMu = as.vector(muVec)
  X = goodG %*% diag(tvec)
  XRest = GRest %*% diag(tvec)
  y = as.vector(-fullDat$subsidence[goodEvent])
  sigmaY = fullDat$Uncertainty[goodEvent]
  n = length(y)
  N = nrow(G)
  nFaux = nrow(fauxDat)
  pAreal = nrow(fault)
  pPoint = nrow(xs)
  
  # get other necessary values
  areas = fault$length * fault$width
  mu = 4e10 # rigidty (value as assumed by Randy)
  
  # load rstan and precompile parallel MC code for faster fitting in parallel.  Also use O3 optimization
  set_cppo(mode = "fast")
  rstan_options(auto_write = FALSE)
  numCores = parallel::detectCores()
  options(mc.cores = numCores)
  
  # do the MCMC with Stan
  if(FrechetFGumbelT) {
    # get Gumbel prior parameters
    slipPriorPar = getGumbelParSlip()
    subPriorPar = getGumbelParSub()
    
    muSlip = slipPriorPar[1]
    scaleSlip = slipPriorPar[2]
    muSub = subPriorPar[1]
    scaleSub = subPriorPar[2]
    
    if(is.infinite(slipThresh))
      predResults <- stan(file = "predsGivenSubMaxAllGumbel.stan", iter=niter, 
                          data=list(n=n, N=N, nFaux=nFaux, pAreal=pAreal, pPoint=pPoint, X=X, XRest=XRest, y=y, 
                                    sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL, tvec=tvec, 
                                    areas=areas, rigidity=mu, muSlip=muSlip, scaleSlip=scaleSlip, 
                                    muSub=muSub, scaleSub=scaleSub), 
                          pars=c("alpha", "subs", "subsRest"), include=FALSE)
    else
      predResults <- stan(file = "predsGivenSubMaxAllGumbelSlipThresh.stan", iter=niter, 
                          data=list(n=n, N=N, nFaux=nFaux, pAreal=pAreal, pPoint=pPoint, X=X, XRest=XRest, y=y, 
                                    sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL, tvec=tvec, 
                                    areas=areas, rigidity=mu, muSlip=muSlip, scaleSlip=scaleSlip, 
                                    muSub=muSub, scaleSub=scaleSub, slipThresh=slipThresh), 
                          pars=c("alpha", "subs", "subsRest"), include=FALSE)
  }
  else {
    # get Frechet prior parameters
    slipPriorPar = getSFrechetParSlip()
    subPriorPar = getSFrechetParSub()
    
    mSlip = slipPriorPar[1]
    scaleSlip = slipPriorPar[2]
    shapeSlip = slipPriorPar[3]
    mSub = subPriorPar[1]
    scaleSub = subPriorPar[2]
    shapeSub = subPriorPar[3]
    
    # do MCMC
    predResults <- stan(file = "predsGivenSubMaxAllFrechet.stan", iter=niter, 
                        data=list(n=n, N=N, nFaux=nFaux, pAreal=pAreal, pPoint=pPoint, X=X, XRest=XRest, y=y, 
                                  sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL, tvec=tvec, 
                                  areas=areas, rigidity=mu, mSlip=mSlip, scaleSlip=scaleSlip, shapeSlip=shapeSlip, 
                                  mSub=mSub, scaleSub=scaleSub, shapeSub=shapeSub), 
                        pars=c("alpha", "subs", "subsRest"), include=FALSE)
  }
  #   # show results
  #   print(predResults, digits = 1)
  
  # extract the samples of the relevant parameters
  tab <- extract(predResults, permuted = TRUE) # return a list of arrays 
  logZetaArealTab = tab$logZetaAreal
  zetaTab = tab$zeta # areal
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = exp(logZetaPointTab)
  M0Tab = tab$seismicMoment
  MwTab = tab$Mw
  
  # compute means, standard deviations, and middle 95% confidence intervals for relevant values
  logZetaArealEsts = colMeans(logZetaArealTab)
  logZetaArealSD = apply(logZetaArealTab, 2, sd)
  logZetaArealMed = apply(logZetaArealTab, 2, median)
  logZetaAreal025 = apply(logZetaArealTab, 2, quantile, probs=0.025)
  logZetaAreal975 = apply(logZetaArealTab, 2, quantile, probs=0.975)
  zetaEsts = colMeans(zetaTab)
  zetaSD = apply(zetaTab, 2, sd)
  zetaMed = apply(zetaTab, 2, median)
  zeta025 = apply(zetaTab, 2, quantile, probs=0.025)
  zeta975 = apply(zetaTab, 2, quantile, probs=0.975)
  logZetaPointEsts = colMeans(logZetaPointTab)
  logZetaPointSD = apply(logZetaPointTab, 2, sd)
  logZetaPointMed = apply(logZetaPointTab, 2, median)
  logZetaPoint025 = apply(logZetaPointTab, 2, quantile, probs=0.025)
  logZetaPoint975 = apply(logZetaPointTab, 2, quantile, probs=0.975)
  zetaPointEsts = colMeans(zetaPointTab)
  zetaPointSD = apply(zetaPointTab, 2, sd)
  zetaPointMed = apply(zetaPointTab, 2, median)
  zetaPoint025 = apply(zetaPointTab, 2, quantile, probs=0.025)
  zetaPoint975 = apply(zetaPointTab, 2, quantile, probs=0.975)
  M0Est = mean(M0Tab)
  M0SD = sd(M0Tab)
  M0Med = median(M0Tab)
  M0975 = quantile(probs=.975, M0Tab)
  M0025 = quantile(probs=.025, M0Tab)
  MwEst = mean(MwTab)
  MwSD = sd(MwTab)
  MwMed = median(MwTab)
  Mw975 = quantile(probs=.975, MwTab)
  Mw025 = quantile(probs=.025, MwTab)
  
  # estimate covariance matrix of params using emprical estimate from MCMC samples if necessary
  if(!estVarMat)
    return(list(predResults=predResults, 
                logZetaArealEsts=logZetaArealEsts, logZetaArealSD=logZetaArealSD, logZetaArealMed=logZetaArealMed, logZetaAreal025=logZetaAreal025, logZetaAreal975=logZetaAreal975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zetaMed=zetaMed, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPointMed=logZetaPointMed, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPointMed=zetaPointMed, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaPointVarMat=NULL, zetaPointVarMat=NULL, 
                M0Est=M0Est, M0SD=M0SD, M0Med=M0Med, M0975=M0975, M0025=M0025, M0s=M0Tab, MwEst=MwEst, MwSD=MwSD, Mw975=Mw975, Mw025=Mw025, Mws=MwTab))
  else {
    logZetaArealTabCntr = sweep(logZetaArealTab, 2, logZetaArealEsts, "-")
    zetaTabCntr = sweep(zetaTab, 2, zetaEsts, "-")
    logZetaPointTabCntr = sweep(logZetaPointTab, 2, logZetaPointEsts, "-")
    zetaPointTabCntr = sweep(zetaPointTab, 2, zetaPointEsts, "-")
    logZetaArealVarMat = t(logZetaArealTabCntr) %*% logZetaArealTabCntr/nrow(logZetaArealTabCntr)
    zetaVarMat = t(zetaTabCntr) %*% zetaTabCntr/nrow(zetaTabCntr)
    logZetaPointVarMat = t(logZetaPointTabCntr) %*% logZetaPointTabCntr/nrow(logZetaPointTabCntr)
    zetaPointVarMat = t(zetaPointTabCntr) %*% zetaPointTabCntr/nrow(zetaPointTabCntr)
    return(list(predResults=predResults, 
                logZetaArealEsts=logZetaArealEsts, logZetaArealSD=logZetaArealSD, logZetaArealMed=logZetaArealMed, logZetaAreal025=logZetaAreal025, logZetaAreal975=logZetaAreal975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zetaMed=zetaMed, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPointMed=logZetaPointMed, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPointMed=zetaPointMed, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaArealVarMat=logZetaArealVarMat, zetaVarMat=zetaVarMat, 
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat, 
                M0Est=M0Est, M0SD=M0SD, M0Med=M0Med, M0975=M0975, M0025=M0025, M0s=M0Tab, MwEst=MwEst, MwSD=MwSD, Mw975=Mw975, Mw025=Mw025, Mws=MwTab))
  }
}

# given the output from stan, computes mean, sd, and middle 95% interval estimators along with
# (if the user requests) variance/covariance matrices for the means
extractStanPredictions2 = function(stanResults, estVarMat=TRUE) {
  # extract the samples of the relevant parameters
  tab <- extract(stanResults, permuted = TRUE) # return a list of arrays 
  logZetaArealTab = tab$logZetaAreal
  zetaTab = tab$zeta
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = exp(logZetaPointTab)
  M0Tab = tab$seismicMoment
  MwTab = tab$Mw
  
  # compute means, standard deviations, and middle 95% confidence intervals for relevant values
  logZetaArealEsts = colMeans(logZetaArealTab)
  logZetaArealSD = apply(logZetaArealTab, 2, sd)
  logZetaAreal025 = apply(logZetaArealTab, 2, quantile, probs=0.025)
  logZetaAreal975 = apply(logZetaArealTab, 2, quantile, probs=0.975)
  zetaEsts = colMeans(zetaTab)
  zetaSD = apply(zetaTab, 2, sd)
  zeta025 = apply(zetaTab, 2, quantile, probs=0.025)
  zeta975 = apply(zetaTab, 2, quantile, probs=0.975)
  logZetaPointEsts = colMeans(logZetaPointTab)
  logZetaPointSD = apply(logZetaPointTab, 2, sd)
  logZetaPoint025 = apply(logZetaPointTab, 2, quantile, probs=0.025)
  logZetaPoint975 = apply(logZetaPointTab, 2, quantile, probs=0.975)
  zetaPointEsts = colMeans(zetaPointTab)
  zetaPointSD = apply(zetaPointTab, 2, sd)
  zetaPoint025 = apply(zetaPointTab, 2, quantile, probs=0.025)
  zetaPoint975 = apply(zetaPointTab, 2, quantile, probs=0.975)
  M0Est = mean(M0Tab)
  M0SD = sd(M0Tab)
  M0975 = quantile(probs=.975, M0Tab)
  M0025 = quantile(probs=.025, M0Tab)
  MwEst = mean(MwTab)
  MwSD = sd(MwTab)
  Mw975 = quantile(probs=.975, MwTab)
  Mw025 = quantile(probs=.025, MwTab)
  
  # estimate covariance matrix of params using emprical estimate from MCMC samples if necessary
  if(!estVarMat)
    return(list(predResults, 
                logZetaArealEsts, logZetaArealSD, logZetaAreal025, logZetaAreal975, 
                zetaEsts, zetaSD, zeta025, zeta975, 
                logZetaPointEsts, logZetaPointSD, logZetaPoint025, logZetaPoint975, 
                zetaPointEsts, zetaPointSD, zetaPoint025, zetaPoint975, 
                M0Est, M0SD, M0975, M0025, MwEst, MwSD, Mw975, Mw025))
  else {
    logZetaArealTabCntr = sweep(logZetaArealTab, 2, logZetaArealEsts, "-")
    zetaTabCntr = sweep(zetaTab, 2, zetaEsts, "-")
    logZetaPointTabCntr = sweep(logZetaPointTab, 2, logZetaPointEsts, "-")
    zetaPointTabCntr = sweep(zetaPointTab, 2, zetaPointEsts, "-")
    logZetaArealVarMat = t(logZetaArealTabCntr) %*% logZetaArealTabCntr/nrow(logZetaArealTabCntr)
    zetaVarMat = t(zetaTabCntr) %*% zetaTabCntr/nrow(zetaTabCntr)
    logZetaPointVarMat = t(logZetaPointTabCntr) %*% logZetaPointTabCntr/nrow(logZetaPointTabCntr)
    zetaPointVarMat = t(zetaPointTabCntr) %*% zetaPointTabCntr/nrow(zetaPointTabCntr)
    return(list(logZetaArealEsts=logZetaArealEsts, logZetaArealSD=logZetaArealSD, logZetaAreal025=logZetaAreal025, logZetaAreal975=logZetaAreal975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaArealVarMat=logZetaArealVarMat, zetaVarMat=zetaVarMat, 
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat, 
                M0Est, M0SD, M0975, M0025, MwEst, MwSD, Mw975, Mw025))
  }
}

###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
# make predictions that threshold slip

# Note that the subsidence data should only be the data from one earthquake.  The returned 
# values relating to beta correspond to log zeta.  For instance, betaEsts is the 
# vector of estimates of log zeta for the given earthquake.  If ev="all", uses all the subsidence 
# data.  Otherwise uses the data from the specified event.
# NOTE: this function should only be called on data for a single fixed earthquake.  
predsGivenSubsidence3 = function(params, muVec=params[2], fault=csz, subDat=dr1, niter=500, estVarMat=TRUE, 
                                 G=NULL, fauxG=NULL, dStar=28000, normalizeTaper=FALSE, 
                                 tvec=taper(getFaultCenters(fault)[,3], lambda=params[1], dStar=dStar, 
                                            normalize=normalizeTaper), 
                                 fauxDat=getFauxDat(), ev="T1", FrechetFGumbelT=TRUE, slipThresh=80, adaptDelta=.99, alg="NUTS") {
  # first rows correspond to GPS data means, last rows correspond to csz 
  # grid cell means
  if(length(muVec) == 1)
    muVec = rep(muVec, nrow(slipDatCSZ)+nrow(fault))
  
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
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
  if(is.null(fauxG)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    fauxG = okadaAll(fault, lonGrid, latGrid, cbind(fauxDat$Lon, fauxDat$Lat), slip=1, poisson=lambda0)
  }
  
  # compute areal log zeta covariance matrix
  if(identical(fault, csz)) {
    arealCSZCor = getArealCorMat(fault)
    SigmaArea = arealCSZCor * sigmaZeta^2
  }
  else {
    SigmaArea = arealZetaCov(params, fault, nDown1=9, nStrike1=12)
  }
  
  # compute point log zeta covariance matrix and cross-covariance
  xs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  SigmaPointToArea = pointArealZetaCov(params, xs, fault, nDown=9, nStrike=12)
  SigmaPoint = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                              smoothness=nuZeta, Distance="rdist.earth", 
                              Dist.args=list(miles=FALSE)) * sigmaZeta^2
  
  # combine covariance matrices:
  SigmaZeta = cbind(rbind(SigmaPoint, t(SigmaPointToArea)), 
                    rbind(SigmaPointToArea, SigmaArea))
  
  # combine subDat and fauxDat
  tmp = data.frame(list(Lon=subDat$Lon, Lat=subDat$Lat, Uncertainty=subDat$Uncertainty, event=subDat$event, 
                        subsidence=subDat$subsidence))
  fullDat = rbind(tmp, fauxDat)
  
  # get relevant variables for the specified event (and make sure we don't base likelihood on faux data)
  fullG = rbind(G, fauxG)
  if(ev == "all")
    goodEvent = c(rep(TRUE, nrow(G)), rep(FALSE, nrow(fauxG)))
  else
    goodEvent = c((events == ev), rep(FALSE, nrow(fauxG)))
  goodG = fullG[goodEvent,]
  GRest = fullG[!goodEvent,]
  
  # set up required data variables in the R environment for Stan
  priorSigmaL = t(chol(SigmaZeta))
  priorMu = as.vector(muVec)
  X = goodG %*% diag(tvec)
  XRest = GRest %*% diag(tvec)
  y = as.vector(-fullDat$subsidence[goodEvent])
  sigmaY = fullDat$Uncertainty[goodEvent]
  n = length(y)
  N = nrow(G)
  nFaux = nrow(fauxDat)
  pAreal = nrow(fault)
  pPoint = nrow(xs)
  
  # get other necessary values
  areas = fault$length * fault$width
  mu = 4e10 # rigidty (value as assumed by Randy)
  
  # load rstan and precompile parallel MC code for faster fitting in parallel.  Also use O3 optimization
  set_cppo(mode = "fast")
  rstan_options(auto_write = FALSE)
  numCores = parallel::detectCores()
  options(mc.cores = numCores)
  
  # do the MCMC with Stan
  predResults <- stan(file = "predsGivenSubThreshSlip.stan", iter=niter, algorithm=alg, 
                      data=list(n=n, N=N, nFaux=nFaux, pAreal=pAreal, pPoint=pPoint, G=goodG, y=y, 
                                sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL, tvec=tvec, 
                                areas=areas, rigidity=mu, slipThresh=slipThresh), 
                      pars=c("alpha", "slips", "standardSubs"), include=FALSE, control=list(adapt_delta=adaptDelta))
  #   # show results
  #   print(predResults, digits = 1)
  
  # extract the samples of the relevant parameters
  tab <- extract(predResults, permuted = TRUE) # return a list of arrays 
  logZetaArealTab = tab$logZetaAreal
  zetaTab = tab$zeta # areal
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = exp(logZetaPointTab)
  M0Tab = tab$seismicMoment
  MwTab = tab$Mw
  
  # compute means, standard deviations, and middle 95% confidence intervals for relevant values
  logZetaArealEsts = colMeans(logZetaArealTab)
  logZetaArealSD = apply(logZetaArealTab, 2, sd)
  logZetaArealMed = apply(logZetaArealTab, 2, median)
  logZetaAreal025 = apply(logZetaArealTab, 2, quantile, probs=0.025)
  logZetaAreal975 = apply(logZetaArealTab, 2, quantile, probs=0.975)
  zetaEsts = colMeans(zetaTab)
  zetaSD = apply(zetaTab, 2, sd)
  zetaMed = apply(zetaTab, 2, median)
  zeta025 = apply(zetaTab, 2, quantile, probs=0.025)
  zeta975 = apply(zetaTab, 2, quantile, probs=0.975)
  logZetaPointEsts = colMeans(logZetaPointTab)
  logZetaPointSD = apply(logZetaPointTab, 2, sd)
  logZetaPointMed = apply(logZetaPointTab, 2, median)
  logZetaPoint025 = apply(logZetaPointTab, 2, quantile, probs=0.025)
  logZetaPoint975 = apply(logZetaPointTab, 2, quantile, probs=0.975)
  zetaPointEsts = colMeans(zetaPointTab)
  zetaPointSD = apply(zetaPointTab, 2, sd)
  zetaPointMed = apply(zetaPointTab, 2, median)
  zetaPoint025 = apply(zetaPointTab, 2, quantile, probs=0.025)
  zetaPoint975 = apply(zetaPointTab, 2, quantile, probs=0.975)
  M0Est = mean(M0Tab)
  M0SD = sd(M0Tab)
  M0Med = median(M0Tab)
  M0975 = quantile(probs=.975, M0Tab)
  M0025 = quantile(probs=.025, M0Tab)
  MwEst = mean(MwTab)
  MwSD = sd(MwTab)
  MwMed = median(MwTab)
  Mw975 = quantile(probs=.975, MwTab)
  Mw025 = quantile(probs=.025, MwTab)
  
  # estimate covariance matrix of params using emprical estimate from MCMC samples if necessary
  if(!estVarMat)
    return(list(predResults=predResults, 
                logZetaArealEsts=logZetaArealEsts, logZetaArealSD=logZetaArealSD, logZetaArealMed=logZetaArealMed, logZetaAreal025=logZetaAreal025, logZetaAreal975=logZetaAreal975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zetaMed=zetaMed, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPointMed=logZetaPointMed, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPointMed=zetaPointMed, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat, 
                M0Est=M0Est, M0SD=M0SD, M0Med=M0Med, M0975=M0975, M0025=M0025, M0s=M0Tab, MwEst=MwEst, MwSD=MwSD, Mw975=Mw975, Mw025=Mw025, Mws=MwTab))
  else {
    logZetaArealTabCntr = sweep(logZetaArealTab, 2, logZetaArealEsts, "-")
    zetaTabCntr = sweep(zetaTab, 2, zetaEsts, "-")
    logZetaPointTabCntr = sweep(logZetaPointTab, 2, logZetaPointEsts, "-")
    zetaPointTabCntr = sweep(zetaPointTab, 2, zetaPointEsts, "-")
    logZetaArealVarMat = t(logZetaArealTabCntr) %*% logZetaArealTabCntr/nrow(logZetaArealTabCntr)
    zetaVarMat = t(zetaTabCntr) %*% zetaTabCntr/nrow(zetaTabCntr)
    logZetaPointVarMat = t(logZetaPointTabCntr) %*% logZetaPointTabCntr/nrow(logZetaPointTabCntr)
    zetaPointVarMat = t(zetaPointTabCntr) %*% zetaPointTabCntr/nrow(zetaPointTabCntr)
    return(list(predResults=predResults, 
                logZetaArealEsts=logZetaArealEsts, logZetaArealSD=logZetaArealSD, logZetaArealMed=logZetaArealMed, logZetaAreal025=logZetaAreal025, logZetaAreal975=logZetaAreal975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zetaMed=zetaMed, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPointMed=logZetaPointMed, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPointMed=zetaPointMed, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaArealVarMat=logZetaArealVarMat, zetaVarMat=zetaVarMat, 
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat, 
                M0Est=M0Est, M0SD=M0SD, M0Med=M0Med, M0975=M0975, M0025=M0025, M0s=M0Tab, MwEst=MwEst, MwSD=MwSD, Mw975=Mw975, Mw025=Mw025, Mws=MwTab))
  }
}

# given the output from stan, computes mean, sd, and middle 95% interval estimators along with
# (if the user requests) variance/covariance matrices for the means
extractStanPredictions3 = function(stanResults, estVarMat=TRUE) {
  # extract the samples of the relevant parameters
  tab <- extract(stanResults, permuted = TRUE) # return a list of arrays 
  logZetaArealTab = tab$logZetaAreal
  zetaTab = tab$zeta
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = exp(logZetaPointTab)
  M0Tab = tab$seismicMoment
  MwTab = tab$Mw
  
  # compute means, standard deviations, and middle 95% confidence intervals for relevant values
  logZetaArealEsts = colMeans(logZetaArealTab)
  logZetaArealSD = apply(logZetaArealTab, 2, sd)
  logZetaAreal025 = apply(logZetaArealTab, 2, quantile, probs=0.025)
  logZetaAreal975 = apply(logZetaArealTab, 2, quantile, probs=0.975)
  zetaEsts = colMeans(zetaTab)
  zetaSD = apply(zetaTab, 2, sd)
  zeta025 = apply(zetaTab, 2, quantile, probs=0.025)
  zeta975 = apply(zetaTab, 2, quantile, probs=0.975)
  logZetaPointEsts = colMeans(logZetaPointTab)
  logZetaPointSD = apply(logZetaPointTab, 2, sd)
  logZetaPoint025 = apply(logZetaPointTab, 2, quantile, probs=0.025)
  logZetaPoint975 = apply(logZetaPointTab, 2, quantile, probs=0.975)
  zetaPointEsts = colMeans(zetaPointTab)
  zetaPointSD = apply(zetaPointTab, 2, sd)
  zetaPoint025 = apply(zetaPointTab, 2, quantile, probs=0.025)
  zetaPoint975 = apply(zetaPointTab, 2, quantile, probs=0.975)
  M0Est = mean(M0Tab)
  M0SD = sd(M0Tab)
  M0975 = quantile(probs=.975, M0Tab)
  M0025 = quantile(probs=.025, M0Tab)
  MwEst = mean(MwTab)
  MwSD = sd(MwTab)
  Mw975 = quantile(probs=.975, MwTab)
  Mw025 = quantile(probs=.025, MwTab)
  
  # estimate covariance matrix of params using emprical estimate from MCMC samples if necessary
  if(!estVarMat)
    return(list(predResults, 
                logZetaArealEsts, logZetaArealSD, logZetaAreal025, logZetaAreal975, 
                zetaEsts, zetaSD, zeta025, zeta975, 
                logZetaPointEsts, logZetaPointSD, logZetaPoint025, logZetaPoint975, 
                zetaPointEsts, zetaPointSD, zetaPoint025, zetaPoint975, 
                M0Est, M0SD, M0975, M0025, MwEst, MwSD, Mw975, Mw025))
  else {
    logZetaArealTabCntr = sweep(logZetaArealTab, 2, logZetaArealEsts, "-")
    zetaTabCntr = sweep(zetaTab, 2, zetaEsts, "-")
    logZetaPointTabCntr = sweep(logZetaPointTab, 2, logZetaPointEsts, "-")
    zetaPointTabCntr = sweep(zetaPointTab, 2, zetaPointEsts, "-")
    logZetaArealVarMat = t(logZetaArealTabCntr) %*% logZetaArealTabCntr/nrow(logZetaArealTabCntr)
    zetaVarMat = t(zetaTabCntr) %*% zetaTabCntr/nrow(zetaTabCntr)
    logZetaPointVarMat = t(logZetaPointTabCntr) %*% logZetaPointTabCntr/nrow(logZetaPointTabCntr)
    zetaPointVarMat = t(zetaPointTabCntr) %*% zetaPointTabCntr/nrow(zetaPointTabCntr)
    return(list(logZetaArealEsts=logZetaArealEsts, logZetaArealSD=logZetaArealSD, logZetaAreal025=logZetaAreal025, logZetaAreal975=logZetaAreal975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                logZetaArealVarMat=logZetaArealVarMat, zetaVarMat=zetaVarMat, 
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat, 
                M0Est, M0SD, M0975, M0025, MwEst, MwSD, Mw975, Mw025))
  }
}


#####
# fit evaluation via cross-validation

# perform cross-validation for the given model and event.  Only try to predict subsidence data for the given event using predictive distribution
doCVSub = function(params, muVec=params[2], fault=csz, subDat=dr1, testEvent="T1", niter=500, 
                  G=NULL, prior=FALSE, normalizeTaper=FALSE, dStar=28000, 
                  tvec=taper(getFaultCenters(fault)[,3], lambda=params[1], normalize=normalizeTaper), 
                  priorMaxSlip=NA, normalModel=FALSE, posNormalModel=FALSE, nFold=20, seed=123, 
                  bySite=FALSE, taperedGPSDat=FALSE, gpsDat=slipDatCSZ, inPar=TRUE, fastPNSim=FALSE) {
  
  # set the seed so that data is scrambled the same way every time, no matter the model
  set.seed(seed)
  
  # first get the event dataset and scramble it
  eventDat = subDat[as.character(subDat$event) == testEvent,]
  N = nrow(eventDat)
  scrambleI = sample(1:N, N, replace=FALSE)
  eventDat = eventDat[scrambleI,]
  
  # if CV is performed by site, nFold argument is ignored
  sites = unique(eventDat$Site)
  siteLats = aggregate(eventDat$Lat, list(eventDat$Site), mean)[,2]
  sortI = sort(siteLats, index.return=TRUE)$ix
  sites = sites[sortI]
  if(bySite)
    nFold = length(sites)
  
  # separate G into components for event of interest and for everything else
  GEvent = G[as.character(subDat$event) == testEvent,]
  
  # now do CV
  indices = seq(1, N, l=nFold+1)
  startIs = indices[1:nFold]
  endIs = indices[2:(nFold+1)]
  MSEs = c()
  biases = c()
  variances = c()
  cols = rainbow(nFold)
  plot(dr1$Lon, dr1$Lat, main="CV split", xlab="Lon", ylab="Lat", type="n")
  if(!inPar) {
    for(k in 1:nFold) {
      # print progress
      print(paste0("iteration ", k, "/", nFold))
      
      # set train and test datasets
      startI = startIs[k]
      endI = endIs[k]
      if(!bySite) {
        testDat = eventDat[startI:endI,]
        trainDat = eventDat[-(startI:endI),]
      }
      else {
        site = sites[k]
        testDat = eventDat[eventDat$Site == site,]
        trainDat = eventDat[eventDat$Site != site,]
      }
      # add test dat to plot
      text(testDat$Lon, testDat$Lat, col=cols[k], labels=as.character(k))
      
      # set G for predicting at training data locations
      if(!bySite)
        thisG = GEvent[-(startI:endI),]
      else
        thisG = GEvent[eventDat$Site != site,]
      
      # generate simulations from prediction areal zeta field fastPNSim: 951.112, else: 1207.864
      system.time(preds <- predsGivenSubsidence(params, muVec, fault, trainDat, niter, FALSE, thisG, prior, 
                                                normalizeTaper=normalizeTaper, dStar=dStar, tvec, priorMaxSlip, 
                                                normalModel, posNormalModel, taperedGPSDat=taperedGPSDat, 
                                                gpsDat=gpsDat, fastPNSim=FALSE))
      zetaSims = preds$resultTab$zeta
      
      # make sure all prediction simulations have the same format
      if(!normalModel)
        zetaSims = t(zetaSims)
      
      # get resulting subsidence predictions
      if(!bySite)
        GPred = GEvent[startI:endI,]
      else
        GPred = matrix(GEvent[eventDat$Site == site,], ncol=ncol(G))
      GTPred = sweep(GPred, 2, tvec, "*")
      subPreds = GTPred %*% zetaSims
      
      # compute MSE
      yTest = -testDat$subsidence
      errs = sweep(subPreds, 1, yTest, "-")
      mse = mean(errs^2)
      
      MSEs[k] = mse
      biases[k] = mean(errs)
      variances[k] = mean(apply(errs, 1, var))
      
      print(paste0("iteration MSE: ", mse))
    }
  }
  else {
    # in this case we do things in parallel
    
    # first plot the CV partition
    for(i in 1:nFold) {
      # set train and test datasets
      startI = startIs[i]
      endI = endIs[i]
      if(!bySite) {
        testDat = eventDat[startI:endI,]
        trainDat = eventDat[-(startI:endI),]
      }
      else {
        site = sites[i]
        testDat = eventDat[eventDat$Site == site,]
        trainDat = eventDat[eventDat$Site != site,]
      }
      
      # add test dat to plot
      text(testDat$Lon, testDat$Lat, col=cols[i], labels=as.character(i))
    }
    
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1]-1) #not to overload your computer
    registerDoParallel(cl)
    
    writeLines(c(""), "log.txt")
    
    # now in parallel compute the CV MSE, bias, and variance
    results = foreach(i= 1:nFold, .combine=cbind, .verbose=TRUE) %dopar% {
      setwd("~/git/M9/")
      
      # write output to the log.txt file
      sink("log.txt", append=TRUE)
      
      # first source everything
      library(fields)
      source("taper.R")
      source('predictions.R')
      source('fitModel.R')
      source('plotSubfault.R')
      source("loadFloodDat.R")
      source("exploratoryAnalysisFuns.R") # -418.9, 319
      source("splines.R")
      
      
      # print progress
      print(paste0("iteration ", i, "/", nFold))
      
      # set train and test datasets
      startI = startIs[i]
      endI = endIs[i]
      if(!bySite) {
        testDat = eventDat[startI:endI,]
        trainDat = eventDat[-(startI:endI),]
      }
      else {
        site = sites[i]
        testDat = eventDat[eventDat$Site == site,]
        trainDat = eventDat[eventDat$Site != site,]
      }
      
      # set G for predicting at training data locations
      if(!bySite)
        thisG = GEvent[-(startI:endI),]
      else
        thisG = GEvent[eventDat$Site != site,]
      
      # generate simulations from prediction areal zeta field
      preds = predsGivenSubsidence(params, muVec, fault, trainDat, niter, FALSE, thisG, prior, 
                                   normalizeTaper=normalizeTaper, dStar=dStar, tvec, priorMaxSlip, 
                                   normalModel, posNormalModel, taperedGPSDat=taperedGPSDat, 
                                   gpsDat=gpsDat)
      zetaSims = preds$resultTab$zeta
      
      # make sure all prediction simulations have the same format
      if(!normalModel)
        zetaSims = t(zetaSims)
      
      # get resulting subsidence predictions
      if(!bySite)
        GPred = GEvent[startI:endI,]
      else
        GPred = matrix(GEvent[eventDat$Site == site,], ncol=ncol(G))
      GTPred = sweep(GPred, 2, tvec, "*")
      subPreds = GTPred %*% zetaSims
      
      # compute MSE
      yTest = -testDat$subsidence
      errs = sweep(subPreds, 1, yTest, "-")
      mse = mean(errs^2)
      
      MSEs[i] = mse
      biases[i] = mean(errs)
      variances[i] = mean(apply(errs, 1, var))
      
      print(paste0("iteration MSE: ", mse))
      
      return(c(mse, mean(errs), mean(apply(errs, 1, var))))
    }
    MSEs = results[1,]
    biases = results[2,]
    variances = results[3,]
    
    #stop cluster
    stopCluster(cl)
  }
  
  return(list(MSE=mean(MSEs), bias=mean(biases), variance=mean(variances), MSEs=MSEs, 
              biases=biases, variances=variances, nFold=nFold, nObs=nrow(eventDat), 
              weight=sum(1/eventDat$Uncertainty^2)))
}


# perform cross-validation for the given model and event.  Only try to predict subsidence data for the given 
# event using marginal distribution.  niter only matters for posnormal model, where mean and var too difficult 
# to compute analytically
getEventMSE = function(params, fault=csz, subDat=dr1, testEvent="T1", niter=500, G=NULL, fauxG=NULL, 
                       normalModel=FALSE, posNormalModel=FALSE, seed=123, normalizeTaper=FALSE, 
                       tvec=taper(getFaultCenters(fault)[,3], lambda=params[1], normalize=normalizeTaper), 
                       dStar=28000, taperedGPSDat=FALSE, gpsDat=slipDatCSZ) {
  
  # set the seed so that data is scrambled the same way every time, no matter the model
  set.seed(seed)
  
  # first get the event dataset and scramble it
  eventDat = subDat[as.character(subDat$event) == testEvent,]
  otherDat = subDat[as.character(subDat$event) != testEvent,]
  N = nrow(eventDat)
  
  # separate G into components for event of interest and for everything else
  GEvent = G[as.character(subDat$event) == testEvent,]
  GOther = G[as.character(subDat$event) != testEvent,]
  
  # plot data
  plot(subDat$Lon, subDat$Lat, main="Event data", xlab="Lon", ylab="Lat", type="n")
  points(eventDat$Lon, eventDat$Lat, pch="+", col="red")
  map("world", "Canada", add=TRUE)
  US(add=TRUE)
  
  muVec = NULL
  
  ## now get MSE, bias, variance
    
  # unless we use a truncated normal model, there is a analytical form for bias, variance.  
  # NOTE: need to check analytical solution first
  # if(!posNormalModel)
  #   niter = 2
  
  # get marginal subsidence distribution (`subsidence` predictions are really uplift)
  if(taperedGPSDat)
    phiZeta = params[length(params)]
  else 
    phiZeta=NULL
  preds = preds(params, niter, fault, muVec, tvec, posNormalModel, normalModel, phiZeta, taperedGPSDat)
  subPreds = predsToSubsidence(params, preds, fault, FALSE, GEvent, eventDat, 
                               posNormalModel, normalModel, tvec, normalizeTaper, dStar)
  
  # under positive normal model, must estimate bias and variance from simulations.  Otherwise 
  # MSE is determined analytically.  NOTE: need to verify analytical solution
  subSims = subPreds$subSims
  err = sweep(subSims, 1, -eventDat$subsidence)
  bias = mean(err)
  variances = apply(err, 1, var)
  variance = mean(variances)
  mse = mean(err^2)
  
  return(list(MSE=mse, bias=bias, variance=variance, weight=sum(1/eventDat$Uncertainty^2), nObs = nrow(eventDat)))
}