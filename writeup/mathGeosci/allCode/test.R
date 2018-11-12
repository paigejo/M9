# function for testing normality assumption of G T zeta
testNormality = function(params, fault=faultGeom, nsim=10000) {
  # draw from marginal distribution of slip
  out = preds(params, fault=fault, nsim=nsim)
  
  # compute subsidences
  out = predsToSubsidence(params, out, fault=fault)
  subSims = out$subSims
  
  # make histogram and add normal approximation
  hist(subSims, freq=F, breaks=550, xlim=c(-2,1), main="Simulated Marginal Subsidences")
  subMu = mean(subSims)
  subSigma = sd(subSims)
  xs = seq(-2, 1, l=200)
  lines(xs, dnorm(xs, subMu, subSigma), col="red")
  # hist(subSims, freq=F, breaks=150, xlim=c(-3, 2))
}
# NOTE: For the 20 or 240 size faults, definitely not normal.  
#       True distribution is much more peaked in the center and heavy tailed

# functions for testing the gradient calculations
sqErrM = function(params, args, normalizeTaper=TRUE, makeLambdaPos=FALSE) {
  XiLat = args[[1]]
  Xi = args[[2]]
  dStar = args[[3]]
  Xmat = args[[4]]
  y = args[[5]]
  varVec = args[[6]]
  taperPar = params
  mu = NULL
  sigma = NULL
  if(length(params) > 5) {
    # compute gradient wrt mu and sigma as well
    mu = params[1]
    sigma = params[2]
    taperPar = params[3:length(params)]
    G = args[[7]]
    arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
    corMatCSZ = arealCSZCor
    covMatCSZ = sigma^2 * corMatCSZ
    expectZeta = exp(mu + diag(covMatCSZ)/2)
    Xmat = G %*% diag(expectZeta)
  }
  
  # compute taper parameters over CSZ grid
  lambdaLats = XiLat %*% taperPar
  lambdaCSZ = Xi %*% taperPar
  
  if(makeLambdaPos && any(c(lambdaLats, lambdaCSZ) < 0))
    return(1e15)
  tvec = taper(csz$depth, lambda=lambdaCSZ, normalize=normalizeTaper, dStar=dStar)
  
  # check how well we match the spline mean estimate
  mean((y - Xmat %*% tvec)^2/varVec)
}
sqErrV = function(params, args, normalizeTaper=TRUE, makeLambdaPos=FALSE) {
  XiLat = args[[1]]
  Xi = args[[2]]
  dStar = args[[3]]
  covZeta = args[[4]]
  G = args[[5]]
  varVec = args[[6]]
  if(length(params) > 5) {
    # compute gradient wrt mu and sigma as well
    mu = params[1]
    sigma = params[2]
    taperPar = params[3:length(params)]
    corZeta = args[[7]]
    covZeta = corZeta * (exp(sigma^2) - 1) * exp(2*mu + sigma^2/2)
  }
  
  # compute taper parameters over CSZ grid
  lambdaLats = XiLat %*% taperPar
  lambdaCSZ = Xi %*% taperPar
  
  if(makeLambdaPos && any(c(lambdaLats, lambdaCSZ) < 0))
    return(1e15)
  tvec = taper(csz$depth, lambda=lambdaCSZ, normalize=normalizeTaper, dStar=dStar)
  
  # check how well we match the spline variance estimate
  GT = sweep(G, 2, tvec, "*")
  VSubsidence = GT %*% covZeta %*% t(GT)
  mean((varVec - diag(VSubsidence))^2)/(2*mean(varVec^2))
}
varY = function(params, args, G, normalizeTaper=TRUE, ind=494) {
  XiLat = args[[1]]
  Xi = args[[2]]
  dStar = args[[3]]
  covZeta = args[[4]]
  G = args[[5]]
  taperPar = params
  if(length(args) > 5) {
    # if pickMuSigma
    mu = params[1]
    sigma = params[2]
    taperPar = params[3:length(params)]
    corZeta = args[[6]]
    covZeta = corZeta * (exp(sigma^2) - 1) * exp(2*mu + sigma^2/2)
    arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
    corMatCSZ = arealCSZCor
    covMatCSZ = sigma^2 * corMatCSZ
    expectZeta = exp(mu + diag(covMatCSZ)/2)
  }
  
  # compute taper parameters and values over CSZ grid
  lambdaLats = XiLat %*% taperPar
  lambdaCSZ = Xi %*% taperPar
  tvec = taper(csz$depth, lambda=lambdaCSZ, normalize=normalizeTaper, dStar=dStar)
  
  # get varY (excluding var epsilon, which is constant w.r.t. taper params)
  GT = sweep(G, 2, tvec, "*")
  GTi = GT[ind,]
  varY = diag(GTi %*% covZeta %*% GTi)
}
taperVals = function(params, depths, Xi, dStar, normalizeTaper=TRUE) {
  # meanGrad = function(expectZeta, lambda, tvec, Xi, G, fault=csz, 
  #                     normalizeTaper=TRUE, dStar=21000) {
  
  #tGrad = taperGrad(csz$depth, Xi, lambda, dStar=dStar, normalize=normalizeTaper)
  
  # compute taper parameters over CSZ grid
  lambdaCSZ = Xi %*% params[3:7]
  taper(depths, lambda=lambdaCSZ, normalize=normalizeTaper, dStar=dStar)
}
# test = jacobian(taperVals, params, depths=cszDepths, Xi=Xi, dStar=dStar, normalizeTaper=normalizeTaper)

# jacobian(subsidenceLnLikMod2Test, params, distMatCSZ=distMatCSZ, gpsDat=gpsDat, subDat=subDat,
#          G=G, nKnots=nKnots, normalizeTaper=normalizeTaper, dStar=dStar, fault=fault, latRange=latRange,
#          normalModel=normalModel, taperedGPSDat=taperedGPSDat)
subsidenceLnLikMod2Test = function(params, distMatCSZ, gpsDat=slipDatCSZ, subDat=dr1, G=NULL, nKnots=5, normalizeTaper=TRUE, 
                                   dStar=28000, fault=csz, latRange = c(40, 50), normalModel=FALSE, taperedGPSDat=FALSE) {
  cszDepths = getFaultCenters(fault)[,3]
  
  ##### get parameters
  muZeta = params[1]
  sigmaZeta = params[2]
  if(!taperedGPSDat)
    splinePar = params[-c(1:2)]
  else
    splinePar = params[-c(1:2, length(params))]
  
  # compute correlation matrix for GPS data
  if(!taperedGPSDat) {
    corPar = getCorPar(normalModel=normalModel)
    phiZeta = corPar$phiZeta
    nuZeta = corPar$nuZeta
  }
  else {
    phiZeta = params[length(params)]
    nuZeta=3/2
  }
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # get taper
  Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambda = Xi %*% splinePar
  tvec = taper(cszDepths, lambda, alpha=2, dStar=dStar, normalize=normalizeTaper)
  
  ##### compute covariance of LOG of Zeta and its Cholesky decomp
  coords = cbind(fault$longitude, fault$latitude)
  corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                             onlyUpper=FALSE, distMat=distMatCSZ, smoothness=nuZeta)
  SigmaZetaCSZ = sigmaZeta^2 * corMatCSZ
  
  ##### Compute Likelihood of subsidence and GPS data
  subsidenceLnLikMod2(rep(muZeta, nrow(fault)), lambda, sigmaZeta, SigmaZetaCSZ, G, subDat=subDat, tvec=tvec, 
                      dStar=dStar, normalModel=normalModel)[1]
}
# jacobian(subsidenceLnLikMod2Test, params, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, gpsDat=gpdDat,
#          subDat=subDat, muZeta=muZeta, G=G, nKnots=nKnots, normalizeTaper=normalizeTaper, dStar=dStar,
#          fault=fault, latRange=latRange, normalModel=normalModel)

# jacobian(GPSLnLikTest, params, gpsDat=gpsDat, fault=fault, muZeta=muZeta, normalizeTaper=normalizeTaper, G=G,
#          dStar=dStar, latRange=latRange, normalModel=normalModel, taperedGPSDat=taperedGPSDat)
# GPSLnLikTest(params, gpsDat=gpsDat, fault=fault, muZeta=muZeta, normalizeTaper=normalizeTaper, G=G,
#              dStar=dStar, latRange=latRange, normalModel=normalModel, taperedGPSDat=taperedGPSDat)
GPSLnLikTest = function(params, gpsDat=slipDatCSZ, fault=csz, muZeta=NULL, normalizeTaper=TRUE, G=NULL, 
                        dStar=28000, latRange=c(40,50), normalModel=FALSE, taperedGPSDat=FALSE, corGPS=FALSE) {
  ##### get parameters
  muZeta = params[1]
  sigmaZeta = params[2]
  if(!taperedGPSDat)
    splinePar = params[-c(1:2)]
  else
    splinePar = params[-c(1:2, length(params))]
  
  # compute correlation matrix for GPS data
  if(!taperedGPSDat) {
    corPar = getCorPar(normalModel=normalModel)
    phiZeta = corPar$phiZeta
    nuZeta = corPar$nuZeta
  }
  else {
    phiZeta = params[length(params)]
    nuZeta=3/2
  }
  coords = cbind(gpsDat$lon, gpsDat$lat)
  corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                             onlyUpper=FALSE, smoothness=nuZeta, 
                             Distance="rdist.earth", Dist.args=list(miles=FALSE))
  
  ##### compute covariance of EXPONENT of Zeta and its Cholesky decomp
  SigmaZetaGPS = sigmaZeta^2 * corMatGPS
  
  if(taperedGPSDat) {
    XiGPS = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
    lambda = XiGPS %*% splinePar
    tvecGPS= taper(gpsDat$Depth, lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  }
  else {
    tvecGPS = rep(1, nrow(gpsDat))
  }
  
  GPSLnLik(muZeta, SigmaZetaGPS, gpsDat, normalModel=normalModel, tvec=tvecGPS, corGPS=corGPS)
}

subLnLikGradTest = function(params, corMatCSZL, G=NULL, corZeta=NULL, normalizeTaper=TRUE, fault=csz, subDat=dr1, 
                            nKnots=5, dStar=26000) {
  
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # compute covariance of LOG of Zeta
  SigmaZetaCSZ = sigmaZeta^2 * corMatCSZL %*% t(corMatCSZL)
  
  # get zeta mean vector
  expectZeta = exp(muVecCSZ + diag(SigmaZetaCSZ)/2)
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  covZeta = exp(sigmaZeta^2 * corMatCSZ) - 1
  covZeta = sweep(covZeta, 1, expectZeta, "*")
  covZeta = sweep(covZeta, 2, expectZeta, "*")
  
  # compute spline basis matrix for subsidence data
  latRange = c(40, 50)
  Xi = bs(fault$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # compute taper vector
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # compute expected uplift
  expectY = G %*% (tvec * expectZeta)
  
  subLnLikGrad(expectZeta, expectY, lambda, Xi, G, covZeta, tvec=tvec, 
               fault, subDat, normalizeTaper, dStar, pow=2, 
               muZeta, sigmaZeta)
}

# test = jacobian(meanYTest, params, G=G, fault=fault, subDat=subDat, nKnots=nKnots,
#                 normalizeTaper=normalizeTaper, dStar=dStar, latRange=latRange, normalModel=normalModel)
meanYTest = function(params, G=NULL, fault=csz, subDat=dr1, nKnots=5, 
                     normalizeTaper=TRUE, dStar=26000, latRange=c(40,50), normalModel=FALSE) {
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  SigmaZetaCSZ = sigmaZeta^2 * corMatCSZ
  
  # get zeta mean vector
  if(!normalModel)
    expectZeta = exp(muVecCSZ + diag(SigmaZetaCSZ)/2)
  else
    expectZeta = muVecCSZ
  
  # compute spline basis matrix for subsidence data
  Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # compute taper vector
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # get G %*% T and expectation of Y
  GT = sweep(G, 2, tvec, "*")
  GT %*% expectZeta
}

getMeanYi = function(params, G, muZeta=NULL, ind=494, nKnots=5, normalizeTaper=TRUE, dStar=21000, fault=csz) {
  
  Xi = getSplineBasis(fault, nKnots=nKnots)
  
  ##### get parameters
  if(length(params) <= 3) {
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
  if(length(params) > 3) {
    if(is.null(muZeta)) {
      muZeta = params[1]
      sigmaZeta = params[2]
      splinePar = params[-c(1:2)]
      tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar)
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
    tvec = taper(cszDepths, lambda=lambda, dStar=dStar, alpha=2)
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  SigmaZetaCSZ = sigmaZeta^2 * corMatCSZ
  
  # get zeta mean vector
  expectZeta = exp(muZetaCSZ + diag(SigmaZetaCSZ)/2)
  
  # get mean Y
  GT = sweep(G, 2, tvec, "*")
  expectYs = GT %*% expectZeta
  
  expectYs[ind]
}

meanGradTest = function(params, G=NULL, fault=csz, subDat=dr1, dStar=26000, latRange=c(40,50), normalModel=FALSE) {
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  
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
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  SigmaZetaCSZ = sigmaZeta^2 * corMatCSZ
  
  # get zeta mean vector
  if(!normalModel)
    expectZeta = exp(muVecCSZ + diag(SigmaZetaCSZ)/2)
  else
    expectZeta = muVecCSZ
  
  meanGrad(expectZeta, lambda, Xi, G, muZeta=muZeta, sigmaZeta=sigmaZeta, 
           dStar=dStar, fault=fault, normalModel=normalModel)
}

covGradTest = function(params, corMatCSZL, G=NULL, normalizeTaper=TRUE, fault=csz, subDat=dr1, 
                       nKnots=5, dStar=26000) {
  
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # compute covariance of LOG of Zeta and its Cholesky decomp
  SigmaZetaCSZ = sigmaZeta^2 * corMatCSZL %*% t(corMatCSZL)
  
  # get zeta mean vector
  expectZeta = exp(muVecCSZ + diag(SigmaZetaCSZ)/2)
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  covZeta = exp(SigmaZetaCSZ) - 1
  covZeta = sweep(covZeta, 1, expectZeta, "*")
  covZeta = sweep(covZeta, 2, expectZeta, "*")
  
  # compute spline basis matrix for subsidence data
  latRange = c(40, 50)
  Xi = bs(fault$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # compute taper vector
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  covGrad(covZeta, lambda, tvec, Xi, G, fault, normalizeTaper, dStar, 
          muZeta, sigmaZeta)
}

# compute SigmaY gradient numerically. Should match with covGrad function.
# If colNum is NA, then return entire thing, otherwise return a specific column 
# of the gradient tensor (c, :, :)
# covYGradNumeric(params, G, nKnots, normalizeTaper, subDat, fault, colNum=1, dStar,
                # latRange, normalModel, taperedGPSDat, distMatCSZ, cszDepths)
covYGradNumeric = function(params, G=NULL, nKnots=5, normalizeTaper=TRUE, subDat=dr1, 
                           fault=csz, colNum=NA, dStar=28000, latRange=c(40,50), normalModel=FALSE, 
                           taperedGPSDat=FALSE, distMatCSZ, cszDepths) {
  
  # return the given column of SigmaY for the input parameters
  getSigmaYCol = function(params, colI=1) {
    SigmaY = getSigmaY(params, G, nKnots, normalizeTaper, subDat, fault, dStar, latRange, normalModel, 
                       taperedGPSDat, distMatCSZ, cszDepths)
    
    return(SigmaY[,colI])
  }
  
  if(is.na(colNum)) {
    # now compute the jacobian for each row of SigmaY and combine into the full gradient
    fullGrad = array(dim=c(nrow(G), length(params), nrow(G)))
    for(i in 1:nrow(G)) {
      # print progress
      if((i %% 25) == 0)
        print(paste0("current col is ", i, "/", nrow(G)))
      
      # compute column gradient
      colGrad = jacobian(getSigmaYCol, params, colI=i)
      fullGrad[,,i] = colGrad
    }
    
    # permute dimensions of array so the third dimension corrsponds to the parameter
    fullGrad = aperm(fullGrad, c(1, 3, 2))
  }
  else {
    # only take gradient of the given column of SigmaY
    fullGrad = jacobian(getSigmaYCol, params, colI=colNum)
  }
  
  return(fullGrad)
}
covYiHessNumeric = function(params, corMatCSZL, G=NULL, nKnots=5, normalizeTaper=TRUE, subDat=dr1, 
                            fault=csz, corZeta=NULL, ind=494, dStar=26000) {
  
  # return the given element of SigmaY for the input parameters
  getSigmaYi = function(params, i=1) {
    SigmaY = getSigmaY(params, G, nKnots, normalizeTaper, subDat, fault, dStar)
    
    return(SigmaY[i,i])
  }
  
  # only take gradient of the given column of SigmaY
  fullHess = hessian(getSigmaYi, params, i=ind)
  
  return(fullHess)
}

# return SigmaY for the given input parameters
getSigmaY = function(params, G=NULL, nKnots=5, normalizeTaper=TRUE, subDat=dr1, 
                     fault=csz, dStar=28000, latRange=c(40, 50), normalModel=FALSE, 
                     taperedGPSDat=FALSE, distMatCSZ, cszDepths) {
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  if(taperedGPSDat)
    phiZeta=params[length(params)]
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  if(!taperedGPSDat) {
    arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
    corMatCSZ = arealCSZCor
    covMatCSZ = sigmaZeta^2 * corMatCSZ
  }
  else {
    coords = cbind(fault$longitude, fault$latitude)
    nuZeta = 3/2
    corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                               onlyUpper=FALSE, distMat=distMatCSZ, smoothness=nuZeta)
    covMatCSZ = sigmaZeta^2 * corMatCSZ
  }
  # expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
  # covZeta = exp(covMatCSZ) - 1
  # covZeta = sweep(covZeta, 1, expectZeta, "*")
  # covZeta = sweep(covZeta, 2, expectZeta, "*")
  
  # compute spline basis matrix for subsidence data
  Xi = getSplineBasis(nKnots=nKnots, fault=fault, latRange=latRange)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # compute taper vector
  
  tvec = taper(cszDepths, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # now ready to compute SigmaY
  out = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, covMatCSZ, G, tvec, 
                             fault=fault, subDat=subDat, normalModel=normalModel)
  SigmaY = out$Sigma
  diag(SigmaY) = diag(SigmaY) + subDat$Uncertainty^2
  
  return(SigmaY)
}
# compute correlation matrix gradient numerically. Should match with corGrad function.
# If colNum is NA, then return entire thing, otherwise return a specific column 
# of the gradient tensor (c, :, :)
# corGradNumeric(params, distMatCSZ, 1)
corGradNumeric = function(params, distMat, colNum=NA) {
  
  # return the given column of SigmaY for the input parameters
  getCorCol = function(params, colI=1) {
    corMat = getCorZeta(params, distMat)
    
    return(corMat[,colI])
  }
  
  if(is.na(colNum)) {
    # now compute the jacobian for each row of SigmaY and combine into the full gradient
    fullGrad = matrix(nrow=nrow(distMat), ncol=nrow(distMat))
    for(i in 1:nrow(distMat)) {
      # print progress
      if((i %% 25) == 0)
        print(paste0("current col is ", i, "/", nrow(distMat)))
      
      # compute column gradient
      colGrad = jacobian(getCorCol, params, colI=i)
      fullGrad[,i] = colGrad
    }
  }
  else {
    # only take gradient of the given column of SigmaY
    fullGrad = jacobian(getCorCol, params, colI=colNum)
  }
  
  return(fullGrad)
}
# get arbitrary correlation matrix for a given distMat (for testing taperedGPS data model)
getCorZeta = function(params, distMat) {
  phiZeta = params[length(params)]
  nuZeta = 3/2
  corMat = stationary.cov(cbind(1,1), Covariance="Matern", theta=phiZeta,
                             onlyUpper=FALSE, distMat=distMat, smoothness=nuZeta)
  
  corMat
}

###### test subsidence log likelihood (each of its parts)
logLikGPYQuad = function(params, G=NULL, nKnots=5, normalizeTaper=TRUE, subDat=dr1, 
                         fault=csz, dStar=26000) {
  
  # first get Variance matrix of Y and invert
  SigmaY = getSigmaY(params, G, nKnots, normalizeTaper, subDat, fault, dStar)
  SigmaYInv = qr.solve(SigmaY)
  
  # get Y (uplift, so negative subsidence)
  Y = -subDat$subsidence
  
  #compute relevant part of log likelihood:
  t(Y) %*% SigmaYInv %*% Y
}

logLikGPYMuQuad = function(params, G=NULL, nKnots=5, normalizeTaper=TRUE, subDat=dr1, 
                           fault=csz, dStar=26000) {
  
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  
  # compute spline basis matrix for subsidence data
  latRange = c(40, 50)
  Xi = bs(fault$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # now get lambdas and taper vector
  lambda = Xi %*% taperPar
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # get log zeta covariance matrix
  
  
  # get muY
  GT = sweep(G, 2, tvec, "*")
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
  muY = GT %*% expectZeta
  
  # first get Variance matrix of Y and invert
  SigmaY = getSigmaY(params, G, nKnots, normalizeTaper, subDat, fault, dStar)
  SigmaYInv = qr.solve(SigmaY)
  
  # get Y (uplift, so negative subsidence)
  Y = -subDat$subsidence
  
  #compute relevant part of log likelihood:
  -2*t(Y) %*% SigmaYInv %*% muY
}

logLikGPMuQuad = function(params, G=NULL, nKnots=5, normalizeTaper=TRUE, subDat=dr1, 
                          fault=csz, dStar=26000) {
  
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  
  # compute spline basis matrix for subsidence data
  latRange = c(40, 50)
  Xi = bs(fault$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # now get lambdas and taper vector
  lambda = Xi %*% taperPar
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # get muY
  GT = sweep(G, 2, tvec, "*")
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
  muY = GT %*% expectZeta
  
  # first get Variance matrix of Y and invert
  SigmaY = getSigmaY(params, G, nKnots, normalizeTaper, subDat, fault, dStar)
  SigmaYInv = qr.solve(SigmaY)
  
  #compute relevant part of log likelihood:
  t(muY) %*% SigmaYInv %*% muY
}

logLikGPLogDetQuad = function(params, G=NULL, nKnots=5, normalizeTaper=TRUE, subDat=dr1, 
                              fault=csz, dStar=26000) {
  
  # first get Variance matrix of Y and invert
  SigmaY = getSigmaY(params, G, nKnots, normalizeTaper, subDat, fault, dStar)
  # method 1: Cholesky decomp
  SigmaYU = chol(SigmaY)
  2*sum(log(diag(SigmaYU)))
}


subLnLikGradTest = function(params, G=NULL, tvec=NULL, 
                            fault=csz, subDat=dr1, normalizeTaper=TRUE, dStar=26000, pow=2, 
                            corZeta=NULL, normalModel=FALSE) {
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # get zeta covariance matrix (covariance of zeta, NOT log(zeta))
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  if(!normalModel) {
    expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
    covZeta = exp(sigmaZeta^2 * corMatCSZ) - 1
    covZeta = sweep(covZeta, 1, expectZeta, "*")
    covZeta = sweep(covZeta, 2, expectZeta, "*")
  }
  else {
    expectZeta = muVecCSZ
    covZeta = sigmaZeta^2 * corMatCSZ
  }
  
  # compute correlation matrix (of zeta, NOT log(zeta))
  if(is.null(corZeta)) {
    sigmas = sqrt(diag(covZeta))
    corZeta = sweep(sweep(covZeta, 1, sigmas, "/"), 2, sigmas, "/")
  }
  
  # compute spline basis matrix for subsidence data
  latRange = c(40, 50)
  Xi = bs(fault$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # compute taper vector
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # get Y MVN approximation
  out = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, covMatCSZ, G, tvec, normalModel=normalModel)
  SigmaY = out$Sigma
  diag(SigmaY) = diag(SigmaY) + subDat$Uncertainty^2
  expectY = out$mu
  
  ## do some precomputations:
  
  # get covariance and mean gradients (covariance and mean of Y)
  SigmaYGrad = covGrad(covZeta, lambda, tvec, Xi, G, fault, normalizeTaper, dStar, 
                       muZeta, sigmaZeta, normalModel=normalModel)
  ExpYGrad = meanGrad(expectZeta, lambda, Xi, G, fault, normalizeTaper, dStar, muZeta, sigmaZeta, 
                      normalModel=normalModel)
  
  nPar = dim(SigmaYGrad)[3]
  
  # compute Sigma_Y^-1
  SigmaYInv = solve(SigmaY)
  
  # compute gradient of Sigma_Y^-1, which is Sigma_Y^-1 dSigma_Y^-1/dtheta_i Sigma_Y^-1
  dSigmaYInv = array(dim = c(dim(SigmaYInv), nPar))
  for(i in 1:nPar) {
    dSigmaYInv[,,i] = SigmaYInv %*% SigmaYGrad[,,i] %*% SigmaYInv
  }
  
  # get Y (which is uplift, so negative of subsidence)
  Y = -subDat$subsidence
  
  ## Now compute likelihood gradient:
  # four parts: Y^T Sigma^-1 Y, -2 Y^T Sigma^-1 mu_Y, mu_y^T Sigma^-1 mu_y, and log|Sigma|
  
  yQuad = 1:nPar
  yMuQuad = 1:nPar
  muQuad = 1:nPar
  logDet = 1:nPar
  for(i in 1:nPar) {
    # gradient of Y^T Sigma^-1 Y
    yQuad[i] = -t(Y) %*% dSigmaYInv[,,i] %*% Y
    
    # gradient of -2 Y^T Sigma^-1 mu_Y
    yMuQuad[i] = -2*(t(Y) %*% SigmaYInv %*% ExpYGrad[,i] - t(Y) %*% dSigmaYInv[,,i] %*% expectY)
    
    # gradient of mu_Y^T Sigma_Y^-1 mu_Y
    muQuad[i] = t(expectY) %*% SigmaYInv %*% ExpYGrad[,i] - t(expectY) %*% dSigmaYInv[,,i] %*% expectY +
      t(ExpYGrad[,i]) %*% SigmaYInv %*% expectY
    
    # gradient of log|Sigma_Y|
    logDet[i] = sum(diag(SigmaYInv %*% SigmaYGrad[,,i]))
  }
  
  # sum to get the full gradient (and multiply by -0.5)
  fullGrad = -0.5 * (yQuad + yMuQuad + muQuad + logDet)
  
  return(fullGrad)
}

GPSLnLikGradTest = function(params) {
  muZeta = params[1]
  sigmaZeta = params[2]
  GPSLnLikGrad(muZeta, sigmaZeta)
}

priorSubLikTest = function(params, G=NULL, fauxObs=getFauxObs(), fauxG=getFauxG(), 
                           subDat=dr1, fault=csz, nKnots=15, normalizeTaper=TRUE, dStar=21000) {
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  # compute spline basis matrix for subsidence data
  latRange = c(40, 50)
  Xi = bs(fault$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # compute taper vector
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # priorSubLik(muZeta, sigmaZeta, tvec, G, fauxG, subDat, fauxObs)
  priorSubLikSFrechet(muZeta, sigmaZeta, tvec, G, fauxG, subDat, fauxObs)
}

priorSlipLikTest = function(params, nKnots=15, dStar=21000, fault=csz, normalizeTaper=TRUE) {
  ##### get parameters
  muZeta = params[1]
  muVecCSZ = rep(muZeta, nrow(fault))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  
  # compute spline basis matrix for subsidence data
  latRange = c(40, 50)
  Xi = bs(fault$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # now get lambdas
  lambda = Xi %*% taperPar
  
  # compute taper vector
  tvec = taper(fault$depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # priorSlipLik(muZeta, sigmaZeta, tvec)
  priorSlipLikSFrechet(muZeta, sigmaZeta, tvec)
}



##### prior gradient testing functions
getQY = function(params, Xi, G, fauxG, subDat=dr1, fauxObs=getFauxObs(), fault=csz, 
                 normalizeTaper=TRUE, dStar=21000, muZeta=NULL) {
  ##### get parameters
  if(length(params) <= 3) {
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
  if(length(params) > 3) {
    if(is.null(muZeta)) {
      muZeta = params[1]
      sigmaZeta = params[2]
      splinePar = params[-c(1:2)]
      tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar)
      Xi = getSplineBasis(fault, nKnots=nKnots)
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
    tvec = taper(cszDepths, lambda=lambda, dStar=dStar, alpha=2)
  
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
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
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
  qY = max(quants)
  qY
}

getQS = function(params, nKnots=5, muZeta=NULL, normalizeTaper=TRUE, dStar=21000, fault=csz, cszDepths=fault$depth) {
  ##### get parameters
  if(length(params) <= 3) {
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
  if(length(params) > 3) {
    if(is.null(muZeta)) {
      muZeta = params[1]
      sigmaZeta = params[2]
      splinePar = params[-c(1:2)]
      tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar)
      Xi = getSplineBasis(fault, nKnots=nKnots)
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
    tvec = taper(cszDepths, lambda=lambda, dStar=dStar, alpha=2)
  
  # get log zeta variances
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
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
  max(quants)
}

getPriorJacobian = function(params, nKnots=5, G, fauxG, subDat=dr1, fauxObs=getFauxObs(), 
                            fault=csz, normalizeTaper=TRUE, dStar=21000, muZeta=NULL) {
  
  Xi = getSplineBasis(fault, nKnots=nKnots)
  
  ##### get parameters
  if(length(params) <= 3) {
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
  if(length(params) > 3) {
    if(is.null(muZeta)) {
      muZeta = params[1]
      sigmaZeta = params[2]
      splinePar = params[-c(1:2)]
      tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar)
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
    tvec = taper(cszDepths, lambda=lambda, dStar=dStar, alpha=2)
  
  priorJacobian(muZeta, sigmaZeta, lambda, Xi, tvec, G, fauxG, subDat, 
                fauxObs=fauxObs, fault, normalizeTaper, dStar)
}

# test if log absolute determinant of Jacobian is correct numerically
getDetJacobian = function(params, G, fauxG) {
  
  Xi = getSplineBasis()
  g = function(pars) {
    qS = getQS(pars)
    qY = getQY(pars, Xi, G, fauxG)
    return(c(pars[-c(1,2)], qS, qY))
  }
  
  log(abs(det(jacobian(g, params))))
}


#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
#########################################################################################################
# functions for testing whether the prior Jacobian is correct.  Do MCMC with constant 1 likelihood and 
# try to recreate the desired distribution on our quantity of interest

plotTestPriorJacobian = function(resultsMCMC) {
  ##### first print acceptance probability
  print(paste0("acceptance probability is: ", mean(resultsMCMC[,ncol(resultsMCMC)])))
  
  ##### now plot traceplots of all parameters
  for(i in 1:(ncol(resultsMCMC-1))) {
    vals = resultsMCMC[,i]
    plot(vals, main=paste0("par ", i, " traceplot"), type="l")
  }
  
  ##### now plot histograms of qS and qY versus prior densities
  qSVals = resultsMCMC[,8]
  qYVals = resultsMCMC[,9]
  
  # qY histogram
  hist(qSVals, main="qS Histogram Vs. Prior Density", freq=FALSE, xlim=c(0, 100), breaks=50)
  xs = seq(0, 100, length=200)
  sfrechetPar = getSFrechetParSlip()
  ys = dfrechet(xs, location=sfrechetPar[1], scale=sfrechetPar[2], shape=sfrechetPar[3])
  lines(xs, ys, col="blue")
  
  # qY histogram
  hist(qYVals, main="qY Histogram Vs. Prior Density", freq=FALSE, xlim=c(0, 7), breaks=50)
  xs = seq(0, 7, length=200)
  sfrechetPar = getSFrechetPar()
  ys = dfrechet(xs, location=sfrechetPar[1], scale=sfrechetPar[2], shape=sfrechetPar[3])
  lines(xs, ys, col="blue")
}

testPriorJacobian = function(initParams=NULL, niter=500, G=NULL, fauxG=NULL, muSigma=.1, sigmaSigma=.02, betaSigma=muSigma) {
  # set fixed constants for test
  normalizeTaper=TRUE
  nKnots=5
  dStar=21000
  muZeta=NULL
  subDat=dr1
  fault=csz
  
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
    # initSplineParams = getInitialSplineEsts(muZetaInit, sigmaZetaInit, lambdaInit, nKnots=nKnots, 
    # normalizeTaper=normalizeTaper, dStar=dStar)
    initSplineParams = getSplineEstsMomentMatch(muZetaInit, sigmaZetaInit, lambdaInit, corMatCSZ, nKnotsMean=nKnots, 
                                                nKnotsVar=nKnots, normalizeTaper=normalizeTaper, dStar=dStar)$betaHat
    
    initParams = c(muZetaInit, sigmaZetaInit, initSplineParams)
  }
  
  # get fixed constants
  phiZeta = 232.5722 # MLE based on fitGPSCovariance result
  nuZeta = 3/2 # we assume this is the Matern smoothness parameter
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  if(is.null(fauxG))
    fauxG = getFauxG()
  
  ##### now that we're done with the precomputations, we're ready to begin MCMC
  tabMCMC = testPriorJacobianMCMC(initParams, niter, G, fauxG, muSigma, sigmaSigma, betaSigma)
  return(tabMCMC)
}

# do the MCMC.  Return a table of the results
testPriorJacobianMCMC = function(initParams, niter, G, fauxG, muSigma=.1, sigmaSigma=.01, betaSigma=muSigma) {
  
  ##### do some precomputations
  Xi = getSplineBasis()
  
  ##### do MCMC
  # table will have model parameters as well as qS, qY, and whether the most recent proposal was accepted (1) 
  # or not (0)
  parTab = matrix(nrow=niter, ncol=length(initParams)+3)
  colnames(parTab) = c("muZeta", "sigmaZeta", paste0("beta", 1:5), "qS", "qY", "accept")
  currPar = initParams
  currLogLik = testPriorLogLik(currPar, G, fauxG)
  accept=1
  for(it in 1:niter) {
    # print progress
    if(it %% 25 == 1)
      print(paste0("Iteration ", it, "/", niter))
    
    # generate current parameter sample
    if(it != 1) {
      # generate proposed parameters
      propPar = rProp(currPar, muSigma, sigmaSigma, betaSigma)
      
      # decide whether or not to accept them
      propLogLik = testPriorLogLik(propPar, G, fauxG)
      numer = propLogLik + dProp(propPar, currPar, sigmaSigma)
      denom = currLogLik + dProp(currPar, propPar, sigmaSigma)
      acceptProb = exp(numer - denom)
      accept = runif(1) < acceptProb
      
      if(accept) {
        currPar = propPar
        currLogLik = propLogLik
      }
    }
    
    # get qS and qY if necessary
    if(accept) {
      qS = getQS(currPar)
      qY = getQY(currPar, Xi, G, fauxG)
    }
    
    # store sample
    parTab[it,] = c(currPar, qS, qY, accept)
  }
  print("MCMC complete.")
  
  return(parTab)
}

# generate random MCMC proposal parameters.  muZeta and betas have normal proposal distribution with 
# mean currPar and standard deviations muSigma and betaSigma, while sigmaZeta has lognormal proposal 
# distribution with median sigmaZeta and standard deviation (on log scale) sigmaSigma.
rProp = function(currPar, muSigma=1, sigmaSigma=1, betaSigma=muSigma) {
  # get parts of currPar
  muZeta = currPar[1]
  sigmaZeta = currPar[2]
  splinePar = currPar[-c(1,2)]
  
  # generate proposals
  rMu = rnorm(1, muZeta, sigmaZeta)
  rSigma = rlnorm(1, log(sigmaZeta), sigmaSigma)
  rBetas = rnorm(length(splinePar), splinePar, betaSigma)
  
  # return proposals
  return(c(rMu, rSigma, rBetas))
}

# only accounts for sigmaZeta probability since that is the only non-symmetric proposal distribution
dProp = function(fromPar, toPar, sigmaSigma=1, logProb=TRUE) {
  # get relevant parts of parameter vectors
  sigmaFrom = fromPar[2]
  sigmaTo = toPar[2]
  
  # return proposal probability
  dlnorm(sigmaTo, log(sigmaFrom), sigmaSigma, log=logProb)
}

testPriorLogLik = function(params, G, fauxG) {
  # set fixed constants for test
  normalizeTaper=TRUE
  nKnots=5
  dStar=21000
  muZeta=NULL
  subDat=dr1
  fault=csz
  
  ##### get parameters
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
      tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar)
      Xi = getSplineBasis(fault, nKnots=nKnots)
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
  
  lnLik = 0
  subPriorLik = priorSubLikSFrechet(muZeta, sigmaZeta, tvec, G, fauxG, subDat)
  slipPriorLik = priorSlipLikSFrechet(muZeta, sigmaZeta, tvec)
  priorLnLik = subPriorLik + slipPriorLik
  lnLik = lnLik + priorLnLik
  
  # add the log jacobian factor
  Xi = getSplineBasis(fault, nKnots=nKnots)
  priorJac = priorJacobian(muZeta, sigmaZeta, lambda, Xi, tvec, G, fauxG, subDat, 
                           fauxObs=getFauxObs(), fault, normalizeTaper, dStar)
  lnLik = lnLik + priorJac
  
  return(lnLik)
}

## test gradient in mean and covariance of X

# covariance of X
# test = covXGradNumeric(params, nKnots, normalizeTaper, gpsDat, 1, dStar, latRange, normalModel,
#                        taperedGPSDat, distMatGPS, corGPS=corGPS)
covXGradNumeric = function(params, nKnots=5, normalizeTaper=TRUE, gpsDat=slipDatCSZ, 
                           colNum=NA, dStar=28000, latRange=c(40,50), normalModel=FALSE, 
                           taperedGPSDat=FALSE, distMatGPS, corGPS=FALSE) {
  
  # return the given column of SigmaY for the input parameters
  getSigmaXCol = function(params, colI=1) {
    SigmaX = getSigmaX(params, nKnots, normalizeTaper, gpsDat, dStar, latRange, normalModel, 
                       taperedGPSDat, distMatGPS, corGPS)
    
    return(SigmaX[,colI])
  }
  
  if(is.na(colNum)) {
    # now compute the jacobian for each row of SigmaY and combine into the full gradient
    fullGrad = array(dim=c(nrow(gpsDat), length(params), nrow(gpsDat)))
    for(i in 1:nrow(gpsDat)) {
      # print progress
      if((i %% 25) == 0)
        print(paste0("current col is ", i, "/", nrow(gpsDat)))
      
      # compute column gradient
      colGrad = jacobian(getSigmaXCol, params, colI=i)
      fullGrad[,,i] = colGrad
    }
    
    # permute dimensions of array so the third dimension corrsponds to the parameter
    fullGrad = aperm(fullGrad, c(1, 3, 2))
  }
  else {
    # only take gradient of the given column of SigmaY
    fullGrad = jacobian(getSigmaXCol, params, colI=colNum)
  }
  
  return(fullGrad)
}

# return SigmaY for the given input parameters
getSigmaX = function(params, nKnots=5, normalizeTaper=TRUE, gpsDat=slipDatCSZ, 
                     dStar=28000, latRange=c(40, 50), normalModel=FALSE, 
                     taperedGPSDat=FALSE, distMatGPS, corGPS=FALSE) {
  ##### get parameters
  muZeta = params[1]
  muVecGPS = rep(muZeta, nrow(gpsDat))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  if(taperedGPSDat)
    phiZeta=params[length(params)]
  
  # compute spline basis matrix for subsidence data
  XiGPS = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  
  # now get lambdas
  lambda = XiGPS %*% taperPar
  
  # compute taper vector
  tvec = taper(gpsDat$Depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # estimate muXi MLE with inverse variance weighting
  x = log(gpsDat$slip)
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  if(normalModel) {
    x = gpsDat$slip
    
    muXi = sum((log(x)- log(muZeta) - log(tvec))*ci)
    gammaEst = exp(muXi)
    sigmaXi = gpsDat$slipErr
  }
  else {
    muXi = sum((x- muZeta - log(tvec))*ci)
  }
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  coords = cbind(gpsDat$lon, gpsDat$lat)
  if(!taperedGPSDat) {
    corPar = getCorPar(normalModel)
    phiZeta=corPar$phiZeta
  }
  nuZeta = 3/2
  corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                             onlyUpper=FALSE, distMat=distMatGPS, smoothness=nuZeta)
  covMatGPS = sigmaZeta^2 * corMatGPS
  
  # now ready to compute SigmaX
  if(taperedGPSDat)
    SigmaX = gammaEst^2 * sweep(sweep(covMatGPS, 1, tvec, "*"), 2, tvec, "*")
  else
    SigmaX = covMatGP
  
  if(corGPS) {
    sigmaX = sqrt(diag(SigmaX))
    Cx = sweep(sweep(SigmaX, 1, 1/sigmaX, "*"), 2, 1/sigmaX, "*")
    SigmaXi = sweep(sweep(Cx, 1, sigmaXi, "*"), 2, sigmaXi, "*")
    SigmaX = SigmaX + SigmaXi
  }
  else
    diag(SigmaX) = diag(SigmaX) + sigmaXi^2
  
  return(SigmaX)
}

# covariance of X inverted
# test = covXInvGradNumeric(params, nKnots, normalizeTaper, gpsDat, 1, dStar, latRange, normalModel, taperedGPSDat, distMatGPS)
covXInvGradNumeric = function(params, nKnots=5, normalizeTaper=TRUE, gpsDat=slipDatCSZ, 
                              colNum=NA, dStar=28000, latRange=c(40,50), normalModel=FALSE, 
                              taperedGPSDat=FALSE, distMatGPS) {
  
  # return the given column of SigmaY for the input parameters
  getSigmaXInvCol = function(params, colI=1) {
    SigmaX = getSigmaX(params, nKnots, normalizeTaper, gpsDat, dStar, latRange, normalModel, 
                       taperedGPSDat, distMatGPS)
    
    return(solve(SigmaX)[,colI])
  }
  
  if(is.na(colNum)) {
    # now compute the jacobian for each row of SigmaY and combine into the full gradient
    fullGrad = array(dim=c(nrow(gpsDat), length(params), nrow(gpsDat)))
    for(i in 1:nrow(gpsDat)) {
      # print progress
      if((i %% 25) == 0)
        print(paste0("current col is ", i, "/", nrow(gpsDat)))
      
      # compute column gradient
      colGrad = jacobian(getSigmaXInvCol, params, colI=i)
      fullGrad[,,i] = colGrad
    }
    
    # permute dimensions of array so the third dimension corrsponds to the parameter
    fullGrad = aperm(fullGrad, c(1, 3, 2))
  }
  else {
    # only take gradient of the given column of SigmaY
    fullGrad = jacobian(getSigmaXInvCol, params, colI=colNum)
  }
  
  return(fullGrad)
}

# mean of X
# jacobian(expXNumeric, params, gpsDat=gpsDat, nKnots=nKnots, latRange=latRange,
#          dStar=dStar, normalizeTaper=normalizeTaper, taperedGPSDat=taperedGPSDat,
#          normalModel=normalModel)

expXNumeric = function(params, gpsDat=slipDatCSZ, nKnots=5, latRange=c(40,50), dStar=28000, 
                       normalizeTaper=TRUE, taperedGPSDat=FALSE, normalModel=FALSE) {
  ##### get parameters
  muZeta = params[1]
  muVecGPS = rep(muZeta, nrow(gpsDat))
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  lambda0 = 0.25
  if(taperedGPSDat)
    phiZeta=params[length(params)]
  
  # compute spline basis matrix for subsidence data
  XiGPS = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  
  # now get lambdas
  lambda = XiGPS %*% taperPar
  
  # compute taper vector
  tvec = taper(gpsDat$Depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # estimate muXi MLE with inverse variance weighting
  x = log(gpsDat$slip)
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  if(normalModel) {
    x = gpsDat$slip
    
    muXi = sum((log(x)- log(muZeta) - log(tvec))*ci)
    gammaEst = exp(muXi)
    sigmaXi = gpsDat$slipErr
  }
  else {
    muXi = sum((x- muZeta - log(tvec))*ci)
  }
  
  # compute MVN expectation
  if(normalModel)
    gammaEst * tvec * muZeta
  else
    muXi + muZeta + log(tvec)
}

###### test GPS log likelihood (each of its parts)
# jacobian(logLikGPXQuad, params, distMatGPS=distMatGPS, nKnots=nKnots, normalizeTaper=normalizeTaper,
#          gpsDat=gpsDat, dStar=dStar, normalModel=normalModel, taperedGPSDat=taperedGPSDat, latRange=latRange,
#          method.args=list(eps=1e-5, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
logLikGPXQuad = function(params, distMatGPS, nKnots=5, normalizeTaper=TRUE, gpsDat=slipDatCSZ, 
                         dStar=28000, normalModel=FALSE, taperedGPSDat=FALSE, latRange=c(40,50)) {
  
  # first get Variance matrix of X and invert
  SigmaX = getSigmaX(params, nKnots, normalizeTaper, gpsDat, dStar, latRange, normalModel, 
                     taperedGPSDat, distMatGPS)
  SigmaXInv = qr.solve(SigmaX)
  
  # get Y (uplift, so negative subsidence)
  x = gpsDat$slip
  if(!normalModel)
    x = log(x)
  
  #compute relevant part of log likelihood:
  t(x) %*% SigmaXInv %*% x
}

# jacobian(logLikGPXMuQuad, params, distMatGPS=distMatGPS, nKnots=nKnots, normalizeTaper=normalizeTaper,
#          gpsDat=gpsDat, dStar=dStar, latRange=latRange, normalModel=normalModel, taperedGPSDat=taperedGPSDat)
logLikGPXMuQuad = function(params, distMatGPS, nKnots=5, normalizeTaper=TRUE, gpsDat=slipDatCSZ, dStar=28000, 
                           latRange=c(40,50), normalModel=FALSE, taperedGPSDat=FALSE) {
  
  ##### get parameters
  muZeta = params[1]
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  
  # get locking data
  x = gpsDat$slip
  if(!normalModel)
    x = log(x)
  
  # first get Variance matrix of X and invert
  SigmaX = getSigmaX(params, nKnots, normalizeTaper, gpsDat, dStar, latRange, normalModel, 
                     taperedGPSDat, distMatGPS)
  SigmaXInv = qr.solve(SigmaX)
  
  # compute spline basis matrix for subsidence data
  XiGPS = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  
  # now get lambdas and taper vector
  lambda = XiGPS %*% taperPar
  tvec = taper(gpsDat$Depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # get scaling parameter conditional estimate
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  if(normalModel) {
    x = gpsDat$slip
    
    muXi = sum((log(x)- log(muZeta) - log(tvec))*ci)
    gammaEst = exp(muXi)
  }
  else {
    muXi = sum((x- muZeta - log(tvec))*ci)
  }
  
  # get muX
  if(normalModel) {
    muX = gammaEst * tvec * muZeta
  }
  else {
    muX = muXi + log(tvec) + muZeta
  }
  
  #compute relevant part of log likelihood:
  -2*t(x) %*% SigmaXInv %*% muX
}

# jacobian(logLikGPMuXQuad, params, distMatGPS=distMatGPS, nKnots=nKnots, normalizeTaper=normalizeTaper,
#          gpsDat=gpsDat, dStar=dStar, latRange=latRange, normalModel=normalModel, taperedGPSDat=taperedGPSDat)
logLikGPMuXQuad = function(params, distMatGPS, nKnots=5, normalizeTaper=TRUE, gpsDat=slipDatCSZ, dStar=28000, 
                          latRange=c(40,50), normalModel=FALSE, taperedGPSDat=FALSE) {
  
  ##### get parameters
  muZeta = params[1]
  sigmaZeta = params[2]
  taperPar = params[3:(2+nKnots)]
  
  # get locking data
  x = gpsDat$slip
  if(!normalModel)
    x = log(x)
  
  # first get Variance matrix of X and invert
  SigmaX = getSigmaX(params, nKnots, normalizeTaper, gpsDat, dStar, latRange, normalModel, 
                     taperedGPSDat, distMatGPS)
  SigmaXInv = qr.solve(SigmaX)
  
  # compute spline basis matrix for subsidence data
  XiGPS = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  
  # now get lambdas and taper vector
  lambda = XiGPS %*% taperPar
  tvec = taper(gpsDat$Depth, lambda=lambda, alpha=2, normalize=normalizeTaper, dStar=dStar)
  
  # get scaling parameter conditional estimate
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  if(normalModel) {
    x = gpsDat$slip
    
    muXi = sum((log(x)- log(muZeta) - log(tvec))*ci)
    gammaEst = exp(muXi)
  }
  else {
    muXi = sum((x- muZeta - log(tvec))*ci)
  }
  
  # get muX
  if(normalModel) {
    muX = gammaEst * tvec * muZeta
  }
  else {
    muX = muXi + log(tvec) + muZeta
  }
  
  #compute relevant part of log likelihood:
  t(muX) %*% SigmaXInv %*% muX
}

# jacobian(logLikGPXLogDetQuad, params, distMatGPS=distMatGPS, nKnots=nKnots, normalizeTaper=normalizeTaper, gpsDat=gpsDat,
#          dStar=dStar, latRange=latRange, normalModel=normalModel, taperedGPSDat=taperedGPSDat)
logLikGPXLogDetQuad = function(params, distMatGPS, nKnots=5, normalizeTaper=TRUE, gpsDat=slipDatCSZ, dStar=28000, 
                               latRange=c(40,50), normalModel=FALSE, taperedGPSDat=FALSE) {
  
  # first get Variance matrix of X and invert
  SigmaX = getSigmaX(params, nKnots, normalizeTaper, gpsDat, dStar, latRange, normalModel, 
                     taperedGPSDat, distMatGPS)
  SigmaXInv = qr.solve(SigmaX)
  
  # method 1: Cholesky decomp
  SigmaXU = chol(SigmaX)
  2*sum(log(diag(SigmaXU)))
}