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
# Here SigmaZeta is for the csz fault geometry GP areal covariance matrix.
# setEventsIndep sets covariance elements from different events to be independent.
# subDat is the subsidence dataset used to generate G
estSubsidenceMeanCov = function(muZeta, lambda, sigmaZeta, SigmaZeta, G, tvec=NULL, 
                                setEventsIndep=TRUE, fault=csz, subDat=dr1) {
  subDatEvents = as.character(subDat$event)
  
  # get vector of taper values
  if(is.null(tvec))
    tvec = taper(fault$depth, lambda = lambda)
  
  # get G T^2 G^T
  GT = sweep(G, 2, tvec, "*")
  
  # compute MVN mean and variance
  if(length(muZeta) == 1)
    muZeta = rep(muZeta, nrow(SigmaZeta))
  if(length(muZeta) == 1)
    sigmaZeta = rep(sigmaZeta, nrow(SigmaZeta))
  else
    sigmaZeta = sqrt(diag(SigmaZeta))
  meanVec = exp(muZeta + sigmaZeta^2/2)
  
  ##### get the full mean vector and covariance matrix
  # mean
  mu = GT %*% cbind(meanVec)
  
  # covariance matrix
  covZeta = exp(SigmaZeta) - 1
  covZeta = sweep(covZeta, 1, meanVec, "*")
  covZeta = sweep(covZeta, 2, meanVec, "*")
  Sigma = GT %*% covZeta %*% t(GT)
  
  ##### set covariance of observations from different events to 0
  if(setEventsIndep) {
    eventEq = matrix(0, nrow=nrow(Sigma), ncol=ncol(Sigma))
    thisUniqueEvents = as.character(sort(factor(unique(subDatEvents), levels=uniqueEvents)))
    for(i in 1:length(thisUniqueEvents)) {
      # get event data
      e = thisUniqueEvents[i]
      eventInds = (subDatEvents == e)
      eventEq = eventEq + (eventInds %o% eventInds)
    }
    Sigma = Sigma * eventEq # set inter-event covariance to 0 here
  }
  
  return(list(mu=mu, Sigma=Sigma))
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
  tvec = taper(fault$depth, lambda=lambda)
  
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
arealZetaCov = function(params, fault1, fault2=NULL, nDown1=3, nStrike1=4, 
                        nDown2=nDown1, nStrike2=nStrike1) {
  if(is.null(fault2))
    fault2 = fault1
  
  # set fixed covariance parameters
  sigmaZetaMLE = params[3]
  nuZeta = 3/2 # assumed Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
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
pointArealZetaCov = function(params, coords, fault, nDown=3, nStrike=4) {
  # we treat the coords like regions with zero length and width
  fakeFault = cbind(rep(1, nrow(coords)), coords, matrix(0, nrow=nrow(coords), ncol=5))
  fakeFault = data.frame(fakeFault)
  names(fakeFault) = names(faultGeom)
  
  # now we just call the usual areal zeta covariance matrix function
  return(arealZetaCov(params, fakeFault, fault, nDown1=1, nStrike1=1, 
                      nDown2=nDown, nStrike2=nStrike))
}
# test = pointArealZetaCov(MLEs, coords, csz, nDown=9, nStrike=12)










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
