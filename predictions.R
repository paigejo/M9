


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
predsGivenGPSFull = function(params, nsim=100, muVec=NULL, gpsDat=slipDatCSZ, fault=csz) {
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
  load("arealCSZCor.RData")
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
  tvec = taper(c(fault$depth, gpsDat$Depth), lambda = lambda)
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
predsGivenGPS = function(params, nsim=100) {
  # get fit MLEs
  if(is.null(muVec)) {
    lambda = params[1]
    muZeta = params[2]
    sigmaZeta = params[3]
    lambda0 = params[4]
    muXi = params[5]
  }
  else {
    lambda = params[1]
    sigmaZeta = params[2]
    lambda0 = params[3]
    muXi = params[4]
  }
  
  # set other relevant parameters
  nuZeta = 3/2 # Matern smoothness
  phiZeta = 232.5722 # fit from fitGPSCovariance()
  
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # get log GPS data
  logX = log(slipDatCSZ$slip)
  
  # get GPS data and CSZ prediction coordinates
  xd = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  xp = cbind(csz$longitude, csz$latitude)
  
  # compute relevant covariance matrices
  SigmaD = stationary.cov(xd, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                          theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  SigmaD = SigmaD + diag(sigmaXi^2)
  SigmaDTildeInv = solve(SigmaDTilde) # NOTE: luckily we only need to do this once.
  
  # These lines have been replaced with the appropriate code for computing 
  # covariance between areal averages of zeta and points or other averages:
  #   SigmaPtoD = stationary.cov(xp, xd, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
  #                              theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  #   SigmaDtoP = t(SigmaPtoD)
  #   SigmaP = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
  #                           theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  SigmaDtoP = pointArealZetaCov(params, xd, csz, nDown=9, nStrike=12)
  SigmaPtoD = t(SigmaDtoP)
  load("arealCSZCor.RData")
  SigmaP = arealCSZCor * sigmaZeta^2
  
  # compute conditional mean and variance of zeta and take Cholesky decomp
  # NOTE: use total variance, conditional covariance gives the SEs for the 
  #       conditional mean estimate.
  muc = muZeta + SigmaPtoD %*% SigmaDTildeInv %*% (logX - muZeta - muXi)
  Sigmac = SigmaP - SigmaPtoD %*% SigmaDTildeInv %*% SigmaDtoP
  #   varDeflation = 1 - mean((muc - muZeta)^2)/sigmaZeta^2
  #   SigmaPred = SigmaP * varDeflation # total variance is conditional variance
  SigmaPred = Sigmac
  SigmaPredL = t(chol(SigmaPred))
  
  # generate predictive simulations
  zSims = matrix(rnorm(nsim*nrow(xp)), nrow=nrow(xp), ncol=nsim)
  logZetaSims = sweep(SigmaPredL %*% zSims, 1, muc, FUN="+") # add muc to each zero mean simulation
  zetaSims = exp(logZetaSims)
  tvec = taper(csz$depth, lambda = lambda)
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # get mean slip prediction field
  meanSlip = exp(muc + diag(SigmaPred)/2) * tvec
  
  ##### generate predictions at GPS locations
  mucGPS = muZeta + SigmaD %*% SigmaDTildeInv %*% (logX - muZeta - muXi)
  SigmacGPS = SigmaD - (SigmaDTilde - 2*diag(sigmaXi^2) + diag(sigmaXi^4) %*% SigmaDTildeInv)
  SigmaPredGPS = SigmacGPS + SigmaD
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, tvec=tvec, muc=muc, Sigmac=Sigmac, 
              mucGPS=mucGPS, SigmacDiagGPS = diag(SigmacGPS), Sigma=SigmaPred, SigmaGPS=SigmaPredGPS))
}

# generate predictions given only the parameter MLEs (no GPS or subsidence data)
preds = function(params, nsim=100, fault=csz) {
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
  sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # get CSZ prediction coordinates
  xp = cbind(fault$longitude, fault$latitude)
  
  # compute relevant covariance matrices
  # NOTE: previous code replaced with code calculating covariance between 
  #       areal averages of zeta
#   Sigma = stationary.cov(xp, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
#                          theta=phiZeta, smoothness=nuZeta, onlyUpper=TRUE) * sigmaZeta^2
  if(! identical(fault, csz))
    Sigma = arealZetaCov(params, fault, nDown1=9, nStrike1=12)
  else {
    # in this case, it takes too long to compute.  Load the precomputed correlation matrix.
    load("arealCSZCor.RData")
    Sigma = arealCSZCor * sigmaZeta^2
  }
  SigmaL = t(chol(Sigma))
  
  # generate predictive simulations
  zSims = matrix(rnorm(nsim*nrow(xp)), nrow=nrow(xp), ncol=nsim)
  logZetaSims = SigmaL %*% zSims + muZeta # add muZeta to each zero mean simulation
  zetaSims = exp(logZetaSims)
  tvec = taper(fault$depth, lambda = lambda)
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # get mean slip prediction field
  meanSlip = rep(exp(muZeta) + sigmaZeta^2/2, nrow(xp)) * tvec
  
  # compute covariance at GPS locations
  xd = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  SigmaD = stationary.cov(xd, Covariance="Matern", Distance="rdist.earth", Dist.args=list(miles=FALSE), 
                          theta=phiZeta, smoothness=nuZeta) * sigmaZeta^2
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, Sigma=Sigma, Sigmac=Sigma, muc=rep(muZeta, nrow(fault)), 
              SigmacGPS = SigmaD, mucGPS=rep(muZeta, nrow(xd))))
}

# generate predictions given only the parameter MLEs (no GPS or subsidence data) at 
# GPS locations as well as areal averages over the fault grid cells.
predsArealAndLoc = function(params, nsim=100, fault=csz, gpsDat=slipDatCSZ, muVec=params[2]) {
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
    load("arealCSZCor.RData")
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
  tvec = taper(c(gpsDat$Depth, fault$depth), lambda = lambda)
  slipSims = sweep(zetaSims, 1, tvec, FUN="*")
  
  # seperate areal average sims from GPS point location sims
  point = 1:nd
  areal = (nd+1):(nd+np)
  slipSimsGPS = slipSims[point,]
  slipSims = slipSims[areal,]
  
  # get mean slip prediction field
  meanSlip = rep(exp(muZeta) + sigmaZeta^2/2, nrow(xp)) * tvec[areal]
  
  return(list(meanSlip=meanSlip, slipSims=slipSims, slipSimsGPS=slipSimsGPS, Sigmac=Sigma, 
              muc=rep(muZeta, nrow(fault)), SigmacGPS = SigmaD, mucGPS=rep(muZeta, nrow(xd))))
}

# Compute subsidence from the prediction simulations using the Okada model.
# Preds is a list with elements named meanSlip (vector) and slipSims (matrix)
predsToSubsidence = function(params, preds, fault=csz, useMVNApprox=TRUE) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  # get predictions from input list
  meanSlip = preds$meanSlip
  slipSims = preds$slipSims
  
  # get Okada linear transformation matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda0)
  
  # transform slips into subsidences
  meanSub = G %*% cbind(meanSlip)
  subSims = G %*% slipSims
  
  # approximate upper and lower 95% quantiles (with either MVN approximate or simulations)
  # NOTE: use preds$Sigma not preds$Sigmac since Sigmac gives the covariance in mean estimate
  if(useMVNApprox) {
    subMVN = estSubsidenceMeanCov(preds$muc, lambda, sigmaZeta, preds$Sigma, G, fault=fault)
    subMu = subMVN$mu
    subSigma = subMVN$Sigma
    l95 = qnorm(.025, mean=subMu, sd=sqrt(diag(subSigma)))
    u95 = qnorm(.975, mean=subMu, sd=sqrt(diag(subSigma)))
  }
  else {
    l95 = apply(subSims, 1, quantile, probs=.025)
    u95 = apply(subSims, 1, quantile, probs=.975)
  }
  return(list(meanSub = meanSub, subSims = subSims, l95=l95, u95=u95))
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
genFullPredsMVN = function(params, nsim=1000, sdAdd=0) {
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
  tvec = taper(csz$depth, lambda=lambda)
  
  # compute G %*% T
  GT = sweep(G, 2, tvec, "*")
  
  # get coordinates for GPS data
  xs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  
  # compute relevant covariances
  load("arealCSZCor.RData")
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
  tvec = taper(c(csz$depth, slipDatCSZ$Depth), lambda = lambda)
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
genFullPreds = function(params, nsim=1000) {
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
  tvec = taper(csz$depth, lambda=lambda)
  
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
  load("arealCSZCor.RData")
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
genFullPredsGPS = function(params, nsim=25000) {
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
  tvec = taper(csz$depth, lambda=lambda)
  
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
  load("arealCSZCor.RData")
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
  tvec = taper(csz$depth, lambda=lambda)
  
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
  load("arealCSZCor.RData")
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
# vector of estimates of log zeta for the given earthquake.
# NOTE: this function should only be called on data for a single fixed earthquake.
predsGivenSubsidence = function(params, muVec=params[2], fault=csz, subDat=dr1, niter=500, estVarMat=TRUE, 
                                G=NULL) {
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
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda)
  }
  
  # get taper values
  tvec = taper(fault$depth, lambda=lambda)
  
  # compute areal log zeta covariance matrix
  if(identical(fault, csz)) {
    load("arealCSZCor.RData")
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
  
  # set up required data variables in the R environment for Stan
  priorSigmaL = t(chol(SigmaZeta))
  priorMu = as.vector(muVec)
  X = G %*% diag(tvec)
  y = as.vector(-subDat$subsidence)
  sigmaY = subDat$Uncertainty
  n = length(y)
  pAreal = nrow(fault)
  pPoint = nrow(xs)
  
  # load rstan and precompile parallel MC code for faster fitting in parallel.  Also use O3 optimization
  set_cppo(mode = "fast")
  rstan_options(auto_write = TRUE)
  numCores = parallel::detectCores()
  options(mc.cores = numCores)
  
  # fit the regression model with Stan
  predResults <- stan(file = "predsGivenSub.stan", iter=niter, 
                      data=list(n=n, pAreal=pAreal, pPoint=pPoint, X=X, y=y, 
                                sigmaY=sigmaY, priorMu=priorMu, priorSigmaL=priorSigmaL))
  
#   # show results
#   print(predResults, digits = 1)
  
  # extract the samples of the relevant parameters
  tab <- extract(predResults, permuted = TRUE) # return a list of arrays 
  betaTab = tab$beta
  zetaTab = tab$zeta
  logZetaPointTab = tab$logZetaPoint
  zetaPointTab = tab$zetaPoint
  
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
  
  if(!estVarMat)
    return(list(predResults, 
                betaEsts, betaSD, beta025, beta975, 
                zetaEsts, zetaSD, zeta025, zeta975, 
                logZetaPointEsts, logZetaPointSD, logZetaPoint025, logZetaPoint975, 
                zetaPointEsts, zetaPointSD, zetaPoint025, zetaPoint975))
  else {
    # compute covariance matrix of estimates using emprical estimate from MCMC samples
    betaTabCntr = sweep(betaTab, 2, betaEsts, "-")
    zetaTabCntr = sweep(zetaTab, 2, zetaEsts, "-")
    logZetaPointTabCntr = sweep(logZetaPointTab, 2, logZetaPointEsts, "-")
    zetaPointTabCntr = sweep(zetaPointTab, 2, zetaPointEsts, "-")
    betaVarMat = t(betaTabCntr) %*% betaTabCntr/nrow(betaTabCntr)
    zetaVarMat = t(zetaTabCntr) %*% zetaTabCntr/nrow(zetaTabCntr)
    logZetaPointVarMat = t(logZetaPointTabCntr) %*% logZetaPointTabCntr/nrow(logZetaPointTabCntr)
    zetaPointVarMat = t(zetaPointTabCntr) %*% zetaPointTabCntr/nrow(zetaPointTabCntr)
    return(list(predResults=predResults, 
                betaEsts=betaEsts, betaSD=betaSD, beta025=beta025, beta975=beta975, 
                zetaEsts=zetaEsts, zetaSD=zetaSD, zeta025=zeta025, zeta975=zeta975,
                logZetaPointEsts=logZetaPointEsts, logZetaPointSD=logZetaPointSD, logZetaPoint025=logZetaPoint025, logZetaPoint975=logZetaPoint975,
                zetaPointEsts=zetaPointEsts, zetaPointSD=zetaPointSD, zetaPoint025=zetaPoint025, zetaPoint975=zetaPoint975,
                betaVarMat=betaVarMat, zetaVarMat=zetaVarMat, 
                logZetaPointVarMat=logZetaPointVarMat, zetaPointVarMat=zetaPointVarMat))
  }
}

# this function updates the estimate of the mean vector of log zeta by weighting the 
# estimates of each earthquake by their variances.
updateMu = function(params, muVec=params[2], fault=csz, niter=1000, gpsDat=slipDatCSZ) {
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
  if(length(muVec) == 1)
    muVec = rep(muVec, nrow(gpsDat) + nrow(fault))
  
  # precompute G
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(fault, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=lambda)
  
  # Only use events that have at least 5 observations (leaves out T12, T11, and T9a)
  minObs = 5
  numObs = table(dr1$event)
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
    subDat = dr1[eventInds,]
    thisG = G[eventInds,]
    eventPreds = predsGivenSubsidence(params, muVec=muVec, fault=fault, subDat=subDat, niter=niter, G=thisG)
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
