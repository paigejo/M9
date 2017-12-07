# test approximation of normal by lognormal

genLogMuSigma = function(mu, Sigma) {
  # estimate mean using median estimate
  logMu = log(mu)
  
  # estimate individual variances using plug in logMu estimate
  logSigmas = sqrt(log(0.5*(sqrt(4*diag(Sigma) * exp(-2*logMu) + 1) + 1)))
  
  # estimate the rest of the variance matrix
  logSigma = sweep(Sigma, 1, exp(-logMu - logSigmas^2/2), "*")
  logSigma = sweep(logSigma, 2, exp(-logMu - logSigmas^2/2), "*")
  logSigma = log(logSigma + 1)
  
  # NOTE: this line theoretically shouldn't change anything, right?
  diag(logSigma) = logSigmas^2
  
  return(list(mu=logMu, Sigma=logSigma))
}

# get the multivariate normal approximating the linear transformation of the lognormal 
# process from the Okada model and tapering matrices.  i.e. get a MVN approximation to 
# the distribution of G %*% T %*% Zeta.  The mean and covariance returned are the true 
# mean and covariance of G %*% T %*% Zeta.
# Here SigmaZeta is for the csz fault geometry GP areal covariance matrix (for LOG 
# of zeta). setEventsIndep sets covariance elements from different events to be 
# independent.  subDat is the subsidence dataset used to generate G
estSubsidenceMeanCov = function(muZeta, lambda, sigmaZeta=NULL, SigmaZeta, G, tvec=NULL, 
                                setEventsIndep=TRUE, fault=csz, subDat=dr1, normalModel=FALSE, 
                                dStar=28000, normalizeTaper=TRUE) {
  
  # get vector of taper values
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], lambda = lambda, dStar=dStar, normalize=normalize, alpha=2)
  
  # get G T
  GT = sweep(G, 2, tvec, "*")
  
  # compute MVN mean and variance
  if(length(muZeta) == 1)
    muZeta = rep(muZeta, nrow(SigmaZeta))
  sigmaZeta = sqrt(diag(SigmaZeta))
  if(!normalModel)
    meanVec = exp(muZeta + diag(SigmaZeta)/2)
  else
    meanVec = muZeta
  
  ##### get the full mean vector and covariance matrix
  # mean
  mu = GT %*% cbind(meanVec)
  
  # covariance matrix
  if(!normalModel) {
    covZeta = exp(SigmaZeta) - 1
    covZeta = sweep(covZeta, 1, meanVec, "*")
    covZeta = sweep(covZeta, 2, meanVec, "*")
  }
  else
    covZeta = SigmaZeta
  Sigma = GT %*% covZeta %*% t(GT)
  
  ##### set covariance of observations from different events to 0
  if(setEventsIndep) {
    eventEq = eventsEqMask(subDat)
    Sigma = Sigma * eventEq # set inter-event covariance to 0 here
  }
  
  return(list(mu=mu, Sigma=Sigma))
}

# this function returns a n x n matrix (for n subDat observations) with 0 or 1 
# elements.  The ijth element is 1 if the event of the ith observation is the 
# same is the event of the jth observation, and 0 otherwise.
eventsEqMask = function(subDat=dr1) {
  subDatEvents = as.character(subDat$event)
  
  eventEq = matrix(0, nrow=nrow(subDat), ncol=nrow(subDat))
  thisUniqueEvents = as.character(sort(factor(unique(subDatEvents), levels=uniqueEvents)))
  for(i in 1:length(thisUniqueEvents)) {
    # get event data
    e = thisUniqueEvents[i]
    eventInds = (subDatEvents == e)
    eventEq = eventEq + (eventInds %o% eventInds)
  }
  return(eventEq)
}

# Function for computing Variance matrix of subsidence data (sorted by event
# so that it is block diagonal) along with its eigendecomposition
# - params: e.g. fixedFitMVN$MLEs.  The Parameters of the model
# - nDown,nStrike: parameters for approximation of areal average covariance 
#                  computations.
getSubsidenceVarianceMat = function(params, fault = faultGeom, G = NULL, nDown=9, 
                                    nStrike=12, subDat=dr1) {
  subDatEvents = as.character(subDat$event)
  thisUniqueEvents = as.character(sort(factor(unique(subDatEvents), levels=uniqueEvents)))
  
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  n = nrow(fault)
  
  # get Okada linear transformation matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  if(is.null(G))
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  
  # get taper values
  tvec = taper(getFaultCenters(fault)[,3], lambda=lambda)
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # get CSZ EQ slip covariance matrix (accounting for areal averages 
  # rather than assuming point observations)
  xp = cbind(fault$longitude, fault$latitude)
  #   SigmaZeta = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
  #                              theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  # load the precomputed areal correlation matrix
  if(identical(fault,csz) && nDown==9 && nStrike==12) {
    load("arealCSZCor.RData")
    SigmaZeta = arealCSZCor * sigmaZeta^2
  }
  else
    SigmaZeta = arealZetaCov(params, fault, nDown1=nDown, nStrike1=nStrike)
  slipParams = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G, tvec, subDat=subDat)
  Sigma = slipParams$Sigma
  
  # add measurement variance
  diag(Sigma) = diag(Sigma) + subDat$Uncertainty^2
  
  # get eigendecompositions of each submatrix
  decomp = list(values=c(), vectors=matrix(0, nrow=nrow(Sigma), ncol=ncol(Sigma)))
  for(i in 1:length(thisUniqueEvents)) {
    # get event data
    e = thisUniqueEvents[i]
    eventInds = (subDatEvents == e)
    
    # get eigendecompositions of submatrix
    subMat = Sigma[eventInds, eventInds]
    subDecomp = eigen(subMat)
    decomp$values = c(decomp$values, subDecomp$values)
    decomp$vectors[eventInds,eventInds] = subDecomp$vectors
  }
  
  return(list(Sigma=Sigma, decomp=decomp))
}

# 2nd function for approximating the areal covariance matrix of Zeta.  
# For areal covariances, just take the average covariance between 
# all sets of points in each region of the input fault.  nDown and 
# nStrike control over which points from each region to average the 
# covariance over.  Uses the centers of the subdivisions of the 
# fault.
# NOTE: The difference between this function and the last one is that 
#       the faults are subdivided all at once to avoid recomputing the 
#       subdivisions (since that takes a ton of time).  
arealZetaCov = function(params=rep(1, 3), fault1=csz, fault2=NULL, nDown1=3, nStrike1=4, 
                        nDown2=nDown1, nStrike2=nStrike1, corPar=NULL, normalModel=FALSE) {
  if(is.null(fault2))
    fault2 = fault1
  
  # set fixed covariance parameters if necessary
  sigmaZetaMLE = params[3]
  if(is.null(corPar)) {
    corPar = getCorPar(normalModel=normalModel)
    nuZeta = corPar$nuZeta
    phiZeta = corPar$phiZeta
  }
  
  # divide the faults
  subfaults1 = divideFault(fault1, nDown=nDown1, nStrike=nStrike1)
  centers1 = matrix(getFaultCenters(subfaults1)[,1:2], ncol=2)
  n1 = nDown1*nStrike1
  subfaults2 = divideFault(fault2, nDown=nDown2, nStrike=nStrike2)
  centers2 = matrix(getFaultCenters(subfaults2)[,1:2], ncol=2)
  n2 = nDown2*nStrike2
  
  # function for computing one row of covariance matrix
  covMatRow = function(subfaultI, fault) {
    # get which rows in coords2 are in the relevant region in fault2
    startI = (subfaultI - 1)*n2 + 1
    endI = startI + n2 - 1
    thisCenters2 = centers2[startI:endI,]
    
    # compute cross covariance matrix
    covMat = stationary.cov(thisCenters2, centers1, "Matern", "rdist.earth", 
                            list(miles=FALSE), theta=phiZeta, smoothness=nuZeta)
    
    # average over all n1*n2 covariances in blocks
    covMat = matrix(covMat, nrow=n1*n2)
    covRow = colMeans(covMat)
    
    return(covRow)
  }
  
  # get full covariance matrix
  covMat = sapply(1:nrow(fault2), covMatRow, fault=fault1) * sigmaZetaMLE^2
  
  return(covMat)
}

# function for computing cross covariance between point value of zeta and 
# areal average of zeta.  Resulting cov mat has dimension (npts x nRegions)
pointArealZetaCov = function(params, coords, fault, nDown=3, nStrike=4, normalModel=FALSE) {
  # we treat the coords like regions with zero length and width
  fakeFault = cbind(rep(1, nrow(coords)), coords, matrix(0, nrow=nrow(coords), ncol=5))
  fakeFault = data.frame(fakeFault)
  names(fakeFault) = names(faultGeom)
  
  # now we just call the usual areal zeta covariance matrix function
  return(arealZetaCov(params, fakeFault, fault, nDown1=1, nStrike1=1, 
                      nDown2=nDown, nStrike2=nStrike, normalModel=normalModel))
}
# test = pointArealZetaCov(MLEs, coords, csz, nDown=9, nStrike=12)

#####
# functions for double lognormal approximation for distribution of subsidences
# take in matrix of subsidence observations at each location.  Columns represent a set of 
# observations in a simulation
# return pPos, pNeg, muPos, sdPos, muNeg, sdNeg
approxDoubleLN = function(simMat) {
  # first get probability subsidence will be positive or negative at each location
  simMatPos = simMat > 0
  pPos = apply(simMatPos, 1, mean)
  pNeg = 1-pPos
  
  # compute parameters for positive lognormal
  logMat = log(simMat)
  logMat[!simMatPos] = NA
  logMat[!is.finite(logMat)] = NA
  muPos = apply(logMat, 1, mean, na.rm=TRUE)
  sdPos = apply(logMat, 1, sd, na.rm=TRUE)
  
  # compute parameters for negative lognormal
  logMat = log(-simMat)
  logMat[simMatPos] = NA
  logMat[!is.finite(logMat)] = NA
  muNeg = apply(logMat, 1, mean, na.rm=TRUE)
  sdNeg = apply(logMat, 1, sd, na.rm=TRUE)
  
  # make sure to put NAs when subsidence will only be positive or negative at a location
  muPos[!is.finite(muPos)] = NA
  sdPos[!is.finite(sdPos)] = NA
  muNeg[!is.finite(muNeg)] = NA
  sdNeg[!is.finite(sdNeg)] = NA
  
  return(list(pPos=pPos, pNeg=pNeg, muPos=muPos, sdPos=sdPos, muNeg=muNeg, sdNeg=sdNeg))
}

# uses Asymmetric shifted Laplace distribution approximation 
# for distribution of subsidences
# take in matrix of subsidence observations at each location.  Columns represent a set of 
# observations in a simulation.  Breaks is the number of histogram breaks 
# used to estimate the mode
# return mHats, lambdaHats, and kappaHats
approxASL = function(simMat, breaks=500) {
  getMode = function(vals) {
    # only consider simulation values in middle 95% of data to ensure 
    # breaks are high enough resolution
    inRange = function(x, low, hi) {
      tmp = x[x > low]
      tmp[tmp < hi]
    }
    vals = inRange(vals, quantile(vals, probs=0.025), quantile(vals, probs=0.975))
    
    # compute histogram of data to estimate mode
    empCDF = hist(vals, plot=F, breaks=breaks)
    modeI = which.max(empCDF$density)
    mode = (empCDF$breaks[modeI] + empCDF$breaks[modeI+1])*0.5
  }
  
  # use the modes of the distributions to estimate shift parameter m
  mHats = apply(simMat, 1, getMode)
  cntrSims = sweep(simMat, 1, mHats, "-")
  simMatPos = cntrSims > 0
  pPos = apply(simMatPos, 1, mean)
  pNeg = 1-pPos
  
  # get asymmetry parameter, kappa
  kappaHats = sqrt(1/pPos - 1)
  
  # get rate parameter lambda
  # lambdaHats = (1 - kappaHats^2)/(kappaHats*(rowMeans(cntrSims)))
  lambdaHats = sqrt((1+kappaHats^4)/(apply(simMat, 1, var)*kappaHats^2))
  
  # alternative estimators of lambda and kappa:
#   meanDiff = mean(subSims[i,])-mi
#   V = var(subSims[i,])
#   lamiTest = sqrt(2/(V - meanDiff^2))
#   y = lamiTest*meanDiff
#   kappaiTest = (-y+sqrt(y^2 + 4))/2
  
  # now reparameterize to agree with ald package
  # sigmaHats = 1/lambdaHats
  # tauHats = 0.5 + sign(pPos - pNeg) * 0.5 * sqrt(1 - 4*kappaHats/(kappaHats^2+1))
  # tauHats = 0.5 + sign(pPos - pNeg) * 0.5 * sqrt(1 - 4/kappaHats)
  
  return(list(mHats=mHats, kappaHats=kappaHats, lambdaHats=lambdaHats))
}

# ASL distribution functions.  ASL = Asymmetric Laplace Distribution
dASL = function(x, m=0, lambda=1, kappa=1, log=FALSE) {
  fracFac = lambda/(kappa + 1/kappa)
  logExpFac = -(x-m)*lambda*sign(x-m)*kappa^(sign(x-m))
  
  if(log)
    return(log(fracFac) + logExpFac)
  else
    return(fracFac * exp(logExpFac))
}

pASL = function(x, m=0, lambda=1, kappa=1) {
  ifelse(x <= m, 
    kappa^2/(1+kappa^2) * exp((x-m)*lambda/kappa), 
    1 - 1/(1+kappa^2) * exp(-lambda * kappa * (x-m)))
}

qASL = function(x, m=0, lambda=1, kappa=1) {
  ifelse(x <= kappa^2/(1+kappa^2), 
         m + kappa/lambda * log(x*(1+kappa^2)/kappa^2), 
         m - log((1+kappa^2)*(1-x))/(lambda*kappa))
}



# DEPRECATED function for approximating the areal covariances between the 
# average of zeta over a subregion and the average over another 
# subregion.  Helper function for compArealSigmaZetaCov (and 
# possibly other functions to be made soon)
arealZetaCovSubfault = function(params, subfault1, subfault2, 
                                      nDown1=3, nStrike1=4, 
                                      nDown2=nDown1, nStrike2=nStrike2) {
  # set fixed covariance parameters
  sigmaZetaMLE = params[3]
  nuZeta = 3/2 # assumed Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # divide fault for approximation
  subfaults1 = divideSubfault(subfault1, nDown1, nStrike1)
  subfaults2 = divideSubfault(subfault2, nDown2, nStrike2)
  coords1 = matrix(getFaultCenters(subfaults1)[,1:2], ncol=2)
  coords2 = matrix(getFaultCenters(subfaults2)[,1:2], ncol=2)
  
  # calculate average covariance
  covMat = stationary.cov(coords1, coords2, Covariance="Matern", Distance="rdist.earth", 
                          theta=phiZeta, smoothness=nuZeta, Dist.args = list(miles=FALSE))
  avgCov = mean(covMat) * sigmaZetaMLE^2
  
  return(avgCov)
}

# DEPRECATED function for approximating the areal covariance matrix of Zeta.  
# For areal covariances, just take the average covariance between 
# all sets of points in each region of the input fault.  nDown and 
# nStrike control over which points from each region to average the 
# covariance over.  Uses the centers of the subdivisions of the 
# fault.
arealZetaCov2 = function(params, fault1, fault2=NULL, nDown1=3, nStrike1=4, 
                        nDown2=nDown1, nStrike2=nStrike1) {
  if(is.null(fault2))
    fault2 = fault1
  
  # set fixed covariance parameters
  sigmaZetaMLE = params[3]
  nuZeta = 3/2 # assumed Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # function for computing one row of covariance matrix
  covMatRow = function(subfault, fault) {
    apply(fault, 1, arealZetaCovSubfault, params=params, subfault1=subfault, 
          nDown1=nDown1, nStrike1=nStrike1, nDown2=nDown2, nStrike2=nStrike2)
  }
  
  # get full covariance matrix
  covMat = apply(fault2, 1, covMatRow, fault=fault1)
  
  return(covMat)
}

getArealCorMat = function(fault=csz, params=NULL, normalModel=FALSE) {
  if(is.null(params)) {
    corPar = getCorPar(normalModel=normalModel)
    nuZeta = corPar$nuZeta
    phiZeta = corPar$phiZeta
    params = c(nuZeta=nuZeta, phiZeta=phiZeta)
  }
  else {
    nuZeta = params[1]
    phiZeta = params[2]
  }
  
  if(identical(fault, csz)) {
    if(max(abs(params - c(nuZeta=3/2, phiZeta=232.5722))) < 1e-3)
      load("arealCSZCor.RData")
    else if(max(abs(params - c(nuZeta=3/2, phiZeta=174.3635))) < 1e-3)
      load("arealCSZCorNormal.RData")
    else
      stop("areal correlation matrix for zeta not computed for the given correlation parameters")
  }
  else {
    # in this case we assume the given fault is a subset of csz with the given latitude range
    hiLatFault = 47.5
    lowLatFault = 42.5
    inRangeFault = (csz$latitude < hiLatFault) & (csz$latitude > lowLatFault)
    
    # NOTE: this is temporary for testing what happens when fitting only part of the data:
    if(max(abs(params - c(nuZeta=3/2, phiZeta=232.5722))) < 1e-3)
      load("arealCSZCor.RData")
    else if(max(abs(params - c(nuZeta=3/2, phiZeta=174.3635))) < 1e-3)
      load("arealCSZCorNormal.RData")
    else
      stop("areal correlation matrix for zeta not computed for the given correlation parameters")
    
    arealCSZCor = arealCSZCor[inRangeFault, inRangeFault]
    # stop("not implemented yet")
  }
  
  return(arealCSZCor)
}

# for naively adjusting mu for positive normal.  Returns what mu should be to get an actual mean of 
# muZeta.  Requires the (constant) muZeta and the covariance matrix covMatCSZ of the GP.
# NOTE: this fit might take a while

getPosNormMu = function(muZeta, covMatCSZ, startN=20, initNewMu=mean(muZeta)) {
  require(tmvtnorm)
  n = startN
  maxN = nrow(covMatCSZ)
  newMu = initNewMu
  
  while(n < maxN) {
    print(paste0("N: ", n))
    newMu = getPosNormMuN(muZeta, covMatCSZ, n, newMuInit=newMu)
    n = n*2
    if(n > maxN)
      n = maxN
  }
  
  # one last fit at maxN:
  print(paste0("N: ", maxN))
  newMu = getPosNormMuN(muZeta, covMatCSZ, maxN, newMuInit=newMu)
  
  # compute the probability all slips are positive
  print(paste0("Probability all positive: ", pmvnorm(upper=rep(0, maxN), mean=rep(-newMu, maxN), sigma=covMatCSZ)))
  
  newMu
}

# helper function for getPosNormMuN
getPosNormMuN = function(muZeta, covMatCSZ, n=10, newMuInit=mean(muZeta), extraFac=1) {
  sigmaTest = covMatCSZ[1:n,1:n]
  
  newMu = newMuInit
  muDiff = Inf
  lastMuDiff=NULL
  while(abs(muDiff) > .01) {
    meanTest = rep(newMu, n)
    out <- mtmvnorm(mean=meanTest, sigma=sigmaTest, lower=rep(0, n), doComputeVariance = FALSE)
    adjustedMu = mean(out$tmean)
    muDiff = adjustedMu - mean(muZeta)
    newMu = newMu - muDiff*extraFac
    
    # for linear extrapolation, extraFac=1, otherwise, usually must be larger for faster convergence
    if(!is.null(lastMuDiff))
      extraFac = extraFac + muDiff/lastMuDiff
    lastMuDiff=muDiff
    
    print(paste0("muDiff: ", muDiff, "; newMu: ", newMu))
  }
  newMu
}