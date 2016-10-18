# fit model:
# S[i] = t(d; lambda) * zeta[i],         zeta ~ exp(GP(muZeta, sigmaZeta * rhoZeta(.)))
# Y[i] = g(S, lambda0)[i] + eps[i],      eps  ~ MVN(0, diag(sigma))
# X[i] = zeta[i] * xi[i],                xi   ~ exp(MVN(muXi, diag(sigmaXi)))
# ( or X[i]/t(d, lambda) = ... )
# Here t(d, lambda) is the normalized double exponential taper function 
# from taper.R with parameter lambda. 
# 

##### load in data and the necessary functions
library(fields)
setwd("~/git/M9/")
# source("loadFloodDat.R")
source("taper.R")
source("okada.R")

##### fit GPS covariance range parameter (phiZeta) in km units
# zeta ~ exp(GP(muZeta, sigmaZeta * rhoZeta(.)))
# inputs:
# sigmaInit: initial guess for sigmaZeta
# phiInit: initial guess for correlation range parameters phi
# NOTE: assumes nuZeta, Matern smoothness parameter, is 3/2.  If phiInit
# is left as NULL, initial guess is set to max distance between locations 
# divided by 3
fitGPSCovariance = function(phiInit=NULL) {
  ##### subset GPS data so it's only over the fault geometry
  slipDatCSZ = getFaultGPSDat()
  
  # guess sigma
  sigmaInit=sd(slipDatCSZ$slip)
  
  # set up and transform data
  coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  logSlip = log(slipDatCSZ$slip)
  
  # NOTE: may want to further transform logSlip by subtracting
  #       log taper
  
  ##### guess lambda (noise/signal)
  # use delta method:
  # sqrt(n)(hat(theta) - theta) -> N(0, sigma^2)
  # f(x) = log(x), f'(x) = 1/x
  # sqrt(n)(log(hat(theta)) - log(theta)) -> N(0, sigma^2/theta^2)
  # take minimum because noise appears small, values are smooth
  noiseEst = min(slipDatCSZ$slipErr^2/slipDatCSZ$slip^2)
  signalEst=  var(logSlip)
  lambdaGuess = noiseEst/(signalEst - noiseEst)
  
  # fit MLE range parameter with chordal distance
  distFun <<- function(x) { rdist.earth(x, miles=FALSE) }
  MLE = MLESpatialProcess(coords, logSlip, lambda.start = lambdaGuess, 
                          cov.args=list(Covariance="Matern", Distance="rdist.earth", smoothness=1.5, 
                                        Dist.args=list(miles=FALSE)), 
                          Distance="distFun")
  phiMLE = MLE$summary[7]
  
  # return MLE
  list(phiMLE=phiMLE)
}
# Results:
# fitGPSCovariance()
# phiZeta = 232.5722

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
  
  ##### compute covariance of exponent of Zeta and its Cholesky decomp
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
# muZeta: mean of GP in exponent of Zeta (constant).  If vector then must be 
#         same length as GPS data
# SigmaZeta: Covariance of exponent of zeta (of dimension |sigmaXi|)
# gpsDat: GPS data restricted to CSZ fault geometry
# logt: log taper values evaluated at the GPS data
GPSLnLik = function(muZeta, SigmaZeta, gpsDat=slipDatCSZ, logt = rep(0, length(sigmaXi))) {
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  logX = log(gpsDat$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  muXi = sum((logX-muZeta)*ci)
  
  # center log GPS data so its mean is 0
  logXCntr = logX - muXi - muZeta - logt
  
  # get likelihood
  Sigma = SigmaZeta + diag(sigmaXi^2)
  logLikGP(logXCntr, chol(Sigma))
}

##### MVN log likelihood
## inputs:
# dat: observations.  If multiple realizations, dat is a matrix with each column being a 
#      different set of observations for one realization
# SimgaU: upper triangular Cholesky decomp of covariance matrix
logLikGP = function(dat, SigmaU) {
  if(!is.matrix(dat))
    dat = matrix(dat, nrow=length(dat))
  n = nrow(dat)
  
  #calcualte GP log likelihood (log likelihood of S0j)
  # (-1/2) t(y) %*% Sigma0^-1 y - (1/2) log det Sigma0 - (n/2) log(2pi)
  # define z = U^-1 %*% y, and x = L^-1 %*% U^-1 %*% y
  # then:
  z = backsolve(SigmaU, dat, transpose=TRUE)
  x = backsolve(SigmaU, z) #make sure x is a column
  # get full loglik.  First sum() only used for multiple realizations
  log.likGP = sum(-(1/2) * colSums(dat * x) - sum(log(diag(SigmaU))) - (n/2)*log(2*pi))
  
  return(log.likGP)
}

getGPSSubfaultDepths = function() {
  
  # helper function for determining if GPS data is within a specific subfault geometry
  getSubfaultGPSDat = function(i) {
    row = csz[i,]
    geom = calcGeom(row)
    corners = geom$corners[,1:2]
    in.poly(cbind(slipDat$lon, slipDat$lat), corners)
  }
  
  # construct logical matrix, where each column is the result of getSubfaultGPSDat(j)
  # Hence, row represents data index, column represents subfault index.
  inSubfaults = sapply(1:nrow(csz), getSubfaultGPSDat)
  
  # get the depth of each subfault by taking mean
  subfaultDepth = function(inds) {
    subfaultDat = slipDat[inds,]
    mean(subfaultDat$Depth)
  }
  apply(inSubfaults, 2, subfaultDepth)
}

##### Compute subsidence data likelihood
# same and subsidenceLnLik but cheats a tad to be faster.  Downside is lnLik SE 
# estimate is anti-conservative.  Uses MC Likelihood estimation, but doesn't 
# assume normality of subsidence data.
# muZeta: mean of log zeta field areal average values.  If vector must be 
#         same length as nrow(fault)
subsidenceLnLikMod = function(G, cszDepths, SigmaZetaL, lambda, muZeta, nsim=3000, 
                              subDat=dr1) {
  # m is number of csz subfaults
  m = nrow(SigmaZetaL)
  
  # get log taper values
  logt = log(taper(cszDepths, lambda=lambda))
  
  # get updated mean vector of GP in exponent of Zeta
  meanZeta = logt + muZeta
  
  # simulate G %*% T %*% Zeta (each col of sims is a simulation.  nsim simulations for each event)
  Zs = matrix(rnorm(m*nsim), nrow=m, ncol=nsim)
  GPs0 = SigmaZetaL %*% Zs # mean zero
  GPs = sweep(GPs0, 1, STATS = meanZeta, FUN="+")
  Ss = exp(GPs)
  sims = G %*% Ss
  
  # y is the data, uncert is its standard deviation
  y = -subDat$subsidence# MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  uncert = subDat$Uncertainty
  
  # Y[i] - (G %*% T %*% Zeta)[i] = eps[i]
  simEps = sweep(-sims, 1, y, "+")
  
#   # calculate the likelihood for the given event
#   eventLnLik = function(i) {
#     eventName = uniqueEvents[i]
#     inds = as.character(dr1$event) == eventName
#     thisUncert = uncert[inds]
#     
#     # calculate simulation likelihoods
#     # NOTE: due to lack of numerical instability must estimate each individual observation's
#     #       likelihood, then log them, then sum them to get total log lik estimate.
#     simLiks = function(simVals) {
#       dnorm(simVals, sd=thisUncert)
#     }
#     colStart = (i-1)*nsim + 1
#     colEnd = i*nsim
#     
#     # get the individual likelihoods for each simulated observation
#     liks = apply(matrix(simEps[inds,], ncol=nsim), 2, simLiks)
#     
#     # take mean and get SE for estimate of each observation's lik
#     if(!is.null(dim(liks))) {
#       obsLiks = apply(liks, 1, mean)
#       obsLikSEs = apply(liks, 1, sd)/sqrt(nsim)
#     }
#     else{
#       obsLiks = mean(liks)
#       obsLikSEs = sd(liks)/sqrt(nsim)
#     }
#     
#     # convert to log like estimates
#     obsLnLiks = log(obsLiks)
#     obsLnLikSEs = obsLikSEs/obsLiks
#     
#     # get full event log-lik estimate
#     lnLik = sum(obsLnLiks)
#     lnLikSE = sqrt(sum(obsLnLikSEs^2))
#     return(c(lnLik=lnLik, lnLikSE=lnLikSE))
#   }
  
  # get the individual likelihoods for each observation
  likMat = apply(simEps, 2, dnorm, sd=dr1$Uncertainty)
  obsLiks = apply(likMat, 1, mean)
  obsLikSEs = apply(likMat, 1, sd)/sqrt(nsim)
  
  # conver to log likelihoods
  obsLnLiks = log(obsLiks)
  obsLnLikSEs = obsLikSEs/obsLiks
  
  #get full data log likelihood estimate
  lnLik = sum(obsLnLiks)
  lnLikSE = sqrt(sum(obsLnLikSEs^2))
  
#   # get all event log likelihoods and standard errors
#   lnLiks = sapply(1:length(uniqueEvents), eventLnLik)
#   lnLikSEs = lnLiks[2,]
#   lnLiks = lnLiks[1,]
#   
#   # now estimate total log likelihood.  Get SE as well
#   lnLik = sum(lnLiks)
#   lnLikSE = sqrt(sum(lnLikSEs^2))
  return(c(lnLik=lnLik, lnLikSE=lnLikSE))
}


############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################

# fitting functions assuming lambda0 (okada model elasticity) is fixed at 0.25

##### function for performing full MLE fit of model
# if muVec not included, params are of form:
# params[1]: lambda
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: lambda0
# otherwise, they are of form:
# params[1]: lambda
# params[2]: sigmaZeta
# params[3]: lambda0
# NOTE: muXi is not required since it's estimated conditionally
# NOTE: currently the Okada model is the most time-consuming step in the optimization.
#       It may help to do block optimization, fixing lambda0 for a while before changing it.
#       Alternatively, we could assume lambda0=0.25 as is often done.
doFixedFit = function(initParams=NULL, nsim=500, useMVNApprox=FALSE, gpsDat=slipDatCSZ, 
                      corMatGPS=NULL, muVec=NULL) {
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
    if(is.null(muVec))
      initParams = c(lambdaInit, muZetaInit, sigmaZetaInit)
    else
      initParams = c(lambdaInit, sigmaZetaInit)
  }
  
  # make sure correct parameter vector format is given if muZeta is a vector
  if(!is.null(muVec)) {
    if(length(initParams) > 2)
      initParams = c(initParams[1], initParams[3]) # lambda and sigmaZeta
  }
  
  # set GPS and CSZ mean fields
  
  
  phiZeta = 232.5722 # MLE based on fitGPSCovariance result
  nuZeta = 3/2 # we assume this is the Matern smoothness parameter
  
  ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
  # NOTE: only the upper triangle is computed for faster computation.  This is 
  #       sufficient for R's Cholesky decomposition
  # since we only want the correlation matrix, set sigmaZeta to 1
  #   coords = cbind(csz$longitude, csz$latitude)
  #   corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
  #                              onlyUpper=TRUE, Distance="rdist.earth", 
  #                              Dist.args=list(miles=FALSE), smoothness=nuZeta)
  # tmpParams = rep(1, 5)
  # corMatCSZ = arealZetaCov(tmpParams, csz, nDown1=9, nStrike1=12)
  # load the precomputed correlation matrix
  load("arealCSZCor.RData")
  corMatCSZ = arealCSZCor
  corMatCSZL = t(chol(corMatCSZ))
  coords = cbind(gpsDat$lon, gpsDat$lat)
  if(is.null(corMatGPS)) {
    corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                               onlyUpper=TRUE, smoothness=nuZeta, 
                               Distance="rdist.earth", Dist.args=list(miles=FALSE))
  }
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(csz)[,3]
  
  ##### Do optimization
  lastParams <<- NULL # used for updating Okada matrix only when necessary
  optimTable <<- NULL # table to save likelihood optimization steps
  controls = list(fnscale = -1, reltol=10^-5)
  opt = optim(initParams, fixedDataLogLik, control=controls, hessian=TRUE, 
              cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, 
              phiZeta=phiZeta, gpsDat=gpsDat, nsim=nsim, 
              useMVNApprox=useMVNApprox, muZeta=muVec)
  
  # test = fullDataLogLik(initParams, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, 
  # phiZeta=phiZeta, slipDatCSZ=slipDatCSZ)
  
  # get results of optimization
  if(is.null(muVec)) {
    lambdaMLE = opt$par[1]
    muZetaMLE = opt$par[2]
    sigmaZetaMLE = opt$par[3]
    muZetaGPS = muZetaMLE
    muZetaCSZ = muZetaMLE
  }
  else {
    lambdaMLE = opt$par[1]
    sigmaZetaMLE = opt$par[2]
    muZetaMLE=muVec
    muZetaGPS = muVec[1:nrow(gpsDat)]
    muZetaCSZ = muVec[(nrow(gpsDat)+1):length(muVec)]
  }
  logLikMLE = opt$value
  hess = opt$hessian
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  logX = log(gpsDat$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  muXiMLE = sum((logX-muZetaGPS)*ci)
  
  # params is in order: lambda, muZeta, sigmaZeta, lambda0, muXi
  if(is.null(muVec))
    MLEs = c(opt$par, 0.25, muXiMLE)
  else
    MLEs = c(opt$par[1], NA, opt$par[2], 0.25, muXiMLE)
  
  # Return results
  list(MLEs=MLEs, lambdaMLE=lambdaMLE, muZetaMLE=muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambda0MLE=0.25, 
       muXiMLE=muXiMLE, logLikMLE=logLikMLE, hess=hess, optimTable=optimTable)
}

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
fixedDataLogLik = function(params, cszDepths, corMatGPS, corMatCSZL, phiZeta, gpsDat=slipDatCSZ, 
                           subDat=dr1, nsim=500, useMVNApprox=FALSE, verbose=TRUE, muZeta=NULL) {
  ##### get parameters
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
  lambda0 = 0.25
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(lastParams)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G <<- okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
  ##### update lastParams
  lastParams <<- params
  
  ##### make sure params are in correct range and update optimTable
  if(lambda < 0 || sigmaZeta < 0) {
    if(vecMu)
      newRow = rbind(c(params[1:2], -10000, -5000, -5000, NA))
    else
      newRow = rbind(c(params, -10000, -5000, -5000, NA))
    if(!vecMu) {
      colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "lnLik", 
                           "subLnLik", "GPSLnLik", "lnLikSE")
    }
    else {
      colnames(newRow) = c("lambda", "sigmaZeta", "lnLik", 
                           "subLnLik", "GPSLnLik", "lnLikSE")
    }
    if(verbose)
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
    sub = subsidenceLnLikMod(G, cszDepths, SigmaZetaCSZL, lambda=0.25, muZetaCSZ, nsim=nsim)
  else {
    SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
    sub = subsidenceLnLikMod2(muZetaCSZ, lambda, sigmaZeta, SigmaZeta, G, subDat=subDat)
  }
    
  GPS = GPSLnLik(muZetaGPS, SigmaZetaGPS, gpsDat)
  lnLik = sub[1] + GPS
  lnLikSE = sub[2]
  
  ##### update parameter optimization table
  if(!vecMu) {
    newRow = rbind(c(params, lnLik, sub[1], GPS, lnLikSE))
    colnames(newRow) = c("lambda", "muZeta", "sigmaZeta", "lnLik", 
                         "subLnLik", "GPSLnLik", "lnLikSE")
  }
  else {
    newRow = rbind(c(params[1:2], lnLik, sub[1], GPS, lnLikSE))
    colnames(newRow) = c("lambda", "sigmaZeta", "lnLik", 
                         "subLnLik", "GPSLnLik", "lnLikSE")
  }
  if(verbose)
    print(newRow)
  optimTable <<- rbind(optimTable, newRow)
  
  ##### return full data likelihood
  lnLik
}

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
# here muZeta can be either a constant or a vector of length ncol(G)
subsidenceLnLikMod2 = function(muZeta, lambda, sigmaZeta, SigmaZeta, G, subDat=dr1) {
  # get vector of taper values
  tvec = taper(csz$depth, lambda = lambda)
  
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, SigmaZeta, G, subDat=subDat)
  mu = mvnApprox$mu
  Sigma = mvnApprox$Sigma
  
  # add data noise to covariance matrix
  diag(Sigma) = diag(Sigma) + subDat$Uncertainty^2
  
  # get log likelihood
  y = -subDat$subsidence# MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  lnLik = logLikGP(y - mu, chol(Sigma))
  
  return(c(lnLik=lnLik, lnLikSE=0))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions relating to variance inflation of GPS data:

# Do n-fold cross-validation over the subsidence data to determine optimal 
# sd inflation additive factor.  In other words, use all GPS data and only 
# part of the subsidence data to fit the model, then compute the predictive 
# likelihood of the left out data.
sdInflationCV = function(initParams=NULL, nfold=10, sdInflationGrid=NULL, 
                         nsim=500, useMVNApprox=TRUE) {
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
    initParams = c(lambdaInit, muZetaInit, sigmaZetaInit)
  }
  
  phiZeta = 232.5722 # MLE based on fitGPSCovariance result
  nuZeta = 3/2 # we assume this is the Matern smoothness parameter
  
  ##### set up default grid sequence of values for sdInflationGrid
  if(is.null(sdInflationGrid)) {
    maxVal = 100
    minVal = .05
    sdInflationGrid = c(0, 10^seq(log10(minVal), log10(maxVal), l=19))
  }
  
  ##### subset GPS data so it's only over the fault geometry
  slipDatCSZ = getFaultGPSDat()
  logX = log(slipDatCSZ$slip)
  
  ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
  # NOTE: only the upper triangle is computed for faster computation.  This is 
  #       sufficient for R's Cholesky decomposition
  # since we only want the correlation matrix, set sigmaZeta to 1
  #   coords = cbind(csz$longitude, csz$latitude)
  #   corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
  #                              onlyUpper=TRUE, Distance="rdist.earth", 
  #                              Dist.args=list(miles=FALSE), smoothness=nuZeta)
  tmpParams = rep(1, 5)
  
  # get coordinates for GPS data
  xs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  
  # calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(csz)[,3]
  
  # precompute areal and point to areal correlation matrices
  params = c(NA,NA,1,NA,NA)
  CorSB = pointArealZetaCov(params, xs, csz, nDown=9, nStrike=12)
  CorS = stationary.cov(xs, Covariance="Matern", theta=phiZeta, 
                          smoothness=nuZeta, Distance="rdist.earth", 
                          Dist.args=list(miles=FALSE))
  load("arealCSZCor.RData")
  CorB = arealCSZCor
  CorBL = t(chol(CorB))
  
  # precompute Okada matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  GFull = okadaAll(csz, lonGrid, latGrid, cbind(dr1$Lon, dr1$Lat), slip=1, poisson=0.25)
  
  ##### Do cross-validation
  # scramle the data (and save how to unscramble it)
  scrambleInds = sample(1:nrow(dr1), nrow(dr1))
  dr1Scramble = dr1[scrambleInds,]
  eventsScramble = events[scrambleInds]
  unscrambleInds = sort(scrambleInds, index.return=TRUE)$ix
  
  # preallocate variables to save during the CV and begin CV
  lnLikGrid = matrix(NA, nrow=length(sdInflationGrid), ncol=nfold)
  maxLnLik = -Inf
  for(i in 1:length(sdInflationGrid)) {
    sdInflation = sdInflationGrid[i]
    print(paste0("sdInflation: ", sdInflation))
    
    for(j in 1:nfold) {
      
      # set subsidence data to be used for fitting
      # (make sure data is still sorted by event)
      predInds = ((1:nrow(dr1)) <= nrow(dr1)*j/nfold) & 
        (((1:nrow(dr1)) > nrow(dr1)*(j-1)/nfold))
      fitInds = !predInds
      scrambledFitInds = scrambleInds[fitInds]
      thisEvents = events[scrambledFitInds]
      ord = order(factor(thisEvents, levels=uniqueEvents))
      YFitOrder = scrambledFitInds[ord]
      subDatFit = dr1[YFitOrder,]
      
      # add to GPS data's variance
      gpsDatInflate = slipDatCSZ
      gpsDatInflate$slipErr = gpsDatInflate$slipErr + sdInflation
      
      # Do optimization
      lastParams <<- NULL # used for updating Okada matrix only when necessary
      optimTable <<- NULL # table to save likelihood optimization steps
      controls = list(fnscale = -1, reltol=10^-5)
      opt = optim(initParams, fixedDataLogLik, control=controls, hessian=TRUE, 
                  cszDepths=cszDepths, corMatGPS=CorS, corMatCSZL=CorBL, 
                  phiZeta=phiZeta, gpsDat=gpsDatInflate, subDat=subDatFit, nsim=nsim, 
                  useMVNApprox=useMVNApprox, verbose=FALSE)
      params = opt$par
      
      # get fit MLEs
      lambda = params[1]
      muZeta = params[2]
      sigmaZeta = params[3]
      lambda0 = 0.25
      
      # estimate muXi MLE with inverse variance weighting
      # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
      # Transformation from additive error to  multiplicative lognormal model with asympototic 
      # median and variance matching.
      sigmaXi = sqrt(log(.5*(sqrt(4*gpsDatInflate$slipErr^2/gpsDatInflate$slip^2 + 1) + 1)))
      logX = log(slipDatCSZ$slip)
      ci = 1/sigmaXi^2
      ci = ci/sum(ci)
      muXi = sum(logX*ci) - muZeta
      
      # set all parameters optimized in fit
      params = c(params, lambda0, muXi)
      
      # get taper values
      tvec = taper(cszDepths, lambda=lambda)
      
      # compute G %*% T
      GFit = GFull[YFitOrder,]
      GTFit = sweep(GFit, 2, tvec, "*")
      
      # get data for prediction
      scrambledPredInds = scrambleInds[predInds]
      thisEvents = events[scrambledPredInds]
      ord = order(factor(thisEvents, levels=uniqueEvents))
      YPredOrder = scrambledPredInds[ord]
      subDatPred = dr1[YPredOrder,]
      
      # compute G %*% T
      GPred = GFull[YPredOrder,]
      GTPred = sweep(GPred, 2, tvec, "*")
      
      # get predictive distribution
      SigmaSB = CorSB * sigmaZeta^2
      SigmaS = CorS * sigmaZeta^2
      SigmaB = CorB * sigmaZeta^2
      SigmaYModFit = diag(exp(muZeta + diag(SigmaB)/2)) %*% t(GTFit)
      SigmaYModPred = diag(exp(muZeta + diag(SigmaB)/2)) %*% t(GTPred)
      # SigmaBYFit = SigmaB %*% SigmaYModFit
      SigmaSYFit = SigmaSB %*% SigmaYModFit
      SigmaSYPred = SigmaSB %*% SigmaYModPred
      subDistn = getSubsidenceVarianceMat(params, fault=csz, G=GFull, subDat=dr1)
      SigmaY = subDistn$Sigma
      SigmaYFit = SigmaY[YFitOrder,YFitOrder]
      SigmaYPred = SigmaY[YPredOrder,YPredOrder]
      SigmaYPredToFit = SigmaY[YPredOrder, YFitOrder]
      
      # now block them into predictions and data covariance matrices
      SigmaP = SigmaYPred
      SigmaD = cbind(rbind(SigmaS + diag(sigmaXi^2), t(SigmaSYFit)), rbind(SigmaSYFit, SigmaYFit))
      SigmaPtoD = cbind(t(SigmaSYPred), SigmaYPredToFit)
      
      # now get predictive likelihood of the left out subsidence data
      fYPred = subLikGivenAll(params, SigmaPtoD, SigmaD, GPred=GPred, GFit=GFit, 
                              subDatPred=subDatPred, subDatFit=subDatFit)
      
      # get log likelihood
      lnLik = fYPred
      lnLikGrid[i,j] = lnLik
      
      # check if this is the best run so far.  If so store results
      if(lnLik > maxLnLik) {
        maxLnLik = lnLik
        bestOpt = opt
        bestOptimTable = optimTable
        bestI = i
      }
      
      print(paste0("Rep ", j, "/", nfold, " complete with lnLik ", lnLik))
    }
    print(paste0("Mean lnLik for inflation factor of ", sdInflation, " is: ", mean(lnLikGrid[i,])))
  }
  
  # get results of optimization for the best CV value of sdInflation
  bestSDInflation = sdInflationGrid[bestI]
  lambdaMLE = bestOpt$par[1]
  muZetaMLE = bestOpt$par[2]
  sigmaZetaMLE = bestOpt$par[3]
  muXiMLE = bestOpt$par[5]
  logLikMLE = bestOpt$value
  hess = bestOpt$hessian
  
  ##### get conditional MLE of muXi
  # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to  multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*(slipDatCSZ$slipErr+bestSDInflation)^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # estimate muXi MLE with inverse variance weighting
  logX = log(slipDatCSZ$slip)
  ci = 1/sigmaXi^2
  ci = ci/sum(ci)
  muXiMLE = sum(logX*ci) - muZetaMLE
  
  # Return results
  list(sdInflationGrid=sdInflationGrid, lnLikGrid=lnLikGrid, bestSDInflation=bestSDInflation, 
       MLEs=c(bestOpt$par, 0.25, muXiMLE), lambdaMLE=lambdaMLE, muZetaMLE = muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambda0MLE=0.25, 
       muXiMLE=muXiMLE, logLikMLE=logLikMLE, hess=hess, optimTable=bestOptimTable)
}

# function for computing likelihood of subsidence data given GPS and subsidence data.  The 
# areal covariance matrix for zeta is unnecessary to input because the 
# correlation matrix is loaded from memory.
# NOTE: This function assumes subsidence data is multivariate normal!
# SigmaP is the covariance matrix of the subsidence prediction data
# SigmaD is covariance matrix for the GPS and subsidence data
# SigmaPtoD is cross covariance matrix between Ypred and GPS and Yfit data
subLikGivenAll = function(params, SigmaPtoD, SigmaD, GPred=NULL, GFit=NULL, gpsDat=slipDatCSZ, 
                          subDatPred=dr1, subDatFit=subDatPred) {
  
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
  logX = log(gpsDat$slip)
  YPred = -subDatPred$subsidence
  YFit = -subDatFit$subsidence
  
  # compute inflated variance sigmaXi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # load covariance matrix of sigma over CSZ fault cells
  load("arealCSZCor.RData")
  SigmaZeta = arealCSZCor * sigmaZeta^2
  
  # get Okada linear transformation matrix
  if(is.null(GFit) || is.null(GPred)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    if(is.null(GPred))
      GPred = okadaAll(csz, lonGrid, latGrid, cbind(subDatPred$Lon, subDatPred$Lat), slip=1, poisson=lambda0)
    if(is.null(GFit))
      GFit = okadaAll(csz, lonGrid, latGrid, cbind(subDatFit$Lon, subDatFit$Lat), slip=1, poisson=lambda0)
  }
  
  # compute taper values and G %*% T
  tvec = taper(csz$depth, lambda=lambda)
  GTPred = sweep(GPred, 2, tvec, "*")
  GTFit = sweep(GFit, 2, tvec, "*")
  
  # compute covariance of prediction subsidence data
  subPredDistn = getSubsidenceVarianceMat(params, csz, GPred, subDat=subDatPred)
  SigmaP = subPredDistn$Sigma
  
  # calculation data and prediction means
  muD = c(rep(muZeta+muXi, length(logX)), GTFit%*%rep(exp(muZeta + sigmaZeta^2/2), length(tvec)))
  muP = GTPred%*%rep(exp(muZeta + sigmaZeta^2/2), length(tvec))
  
  # compute predictive distribution of subsidence given fit data assuming 
  # multivariate normality
  predDistn = conditionalNormal(Xd=c(logX, YFit), muD=muD, muP=muP, SigmaP=SigmaP, 
                                SigmaD=SigmaD, SigmaPtoD=SigmaPtoD)
  muc = predDistn$muc
  #Sigmac = predDistn$Sigmac
  
  # NOTE: the predictive covariance should not be the conditional 
  #       covariance, because that represents the covariance in the mean
  #       rather than covariance in each individual observation.  Instead, 
  #       I use the deflated marginal variance, where the variance is 
  #       deflated using the formula (MSE = var + bias^2).  Ideally, 
  #       the model var should be refit after assuming this conditional mean 
  #       possibly iteratively instead of using this approximation.
  varDeflation = 1 - mean((muc - muZeta)^2)/sigmaZeta^2
  
  ## compute f(Y_Pred|X,Y_Fit) (still assuming multivariate normality)
  logLikPred = logLikGP(YPred - muc, chol(SigmaP * varDeflation))
  
  logLikPred
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions for performing iterative mean-parameter estimation

# fit model with nonstationary mean using iterative method and data imputation.  
# params are updated (fitModelParams is called) maxIter times, and the mean is 
# updated (fitModelMean is called) maxIter - 1 times.  loadDat is the name of 
# an RData file containing the progress of this function, and if it it included 
# this function will start again from where it left off when loadDat was saved.  
# saveFile is the name of the RData file to save this function's progress to.  
# other inputs are inputs to doFixedFit() and updateMu() functions (niterMCMC 
# is the niter input to updateMu()).
fitModelIterative = function(initParams=NULL, nsim=500, useMVNApprox=TRUE, gpsDat=slipDatCSZ, 
                             corMatGPS=NULL, muVec=NULL, maxIter=5, fault=csz, niterMCMC=250, 
                             loadDat=NULL, saveFile="iterFitProgress.RData") {
  funIns = list(initParams=initParams, nsim=nsim, useMVNApprox=useMVNApprox, gpsDat=gpsDat, 
                corMatGPS=corMatGPS, muVec=muVec, maxIter=maxIter, fault=fault, niterMCMC=niterMCMC)
  
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
    initParams = c(lambdaInit, muZetaInit, sigmaZetaInit)
  }
  
  # store each updated mean and MLE sets
  env = environment()
  assign("muMatCSZ", matrix(nrow=nrow(fault), ncol=0), envir=env)
  assign("muMatGPS", matrix(nrow=nrow(slipDatCSZ), ncol=0), envir=env)
  assign("parMat", matrix(nrow=5, ncol=0), envir=env)
  assign("logLiks", c(), envir=env)
  
  # store optimization tables and other important things
  assign("optimTables", list(), envir=env)
  assign("lastStanResults", list(), envir=env)
  assign("lastHess", c(), envir=env)
  
  #called maxIter times
  fitModelParams = function(params, muVec=NULL, currIter=1) {
    # save current progress
    muMatCSZ=get("muMatCSZ", envir=env)
    muMatGPS=get("muMatGPS", envir=env)
    parMat=get("parMat", envir=env)
    logLiks=get("logLiks", envir=env)
    optimTables=get("optimTables", envir=env)
    lastHess=get("lastHess", envir=env)
    lastStanResults=get("lastStanResults", envir=env)
    
    state = list(muMatCSZ=muMatCSZ, muMatGPS=muMatGPS, parMat=parMat, logLiks=logLiks, 
                 optimTables=optimTables, lastHess=lastHess, lastStanResults=lastStanResults)
    subFunName = "fitModelParams"
    subFunIns = list(params=params, muVec=muVec, currIter=currIter)
    save(funIns, state, subFunName, subFunIns, file=saveFile)
    
    # base case: if past maxIter, resutrn results
    if(currIter > maxIter)
      return(list(params=params, muVec=muVec))
    
    print(paste0("fitting model parameters, current iteration is: ", currIter))
    
    # fit the parameters
    parFit = doFixedFit(params, nsim, useMVNApprox, gpsDat, corMatGPS, muVec)
#     list(MLEs=c(opt$par, 0.25, muXiMLE), lambdaMLE=lambdaMLE, muZetaMLE=muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambda0MLE=0.25, 
#          muXiMLE=muXiMLE, logLikMLE=logLikMLE, hess=hess, optimTable=optimTable)
    
    # update parMat with new parameters
    newPar = parFit$MLEs
    parMat = get("parMat", envir=env)
    if(length(newPar) < 5)
      assign("parMat", cbind(parMat, c(newPar[1], NA, newPar[2:length(newPar)])), envir=env)
    else
      assign("parMat", cbind(parMat, newPar), envir=env)
    
    # update optimTables, hessian, logLik
    optimTable = get("optimTables", envir=env)
    logLiks = get("logLiks", envir=env)
    assign("optimTables", c(optimTables, list(parFit$optimTable)), envir=env)
    assign("lastHess", parFit$hess, envir=env)
    assign("logLiks", c(logLiks, parFit$logLikMLE), envir=env)
    
    if(is.null(muVec))
      return(fitModelMean(newPar, currIter=currIter+0.5))
    else
      return(fitModelMean(newPar, muVec, currIter=currIter+0.5))
  }
  
  #called maxIter - 1  times
  fitModelMean = function(params, muVec=params[2], currIter=1) {
    # save current progress
    muMatCSZ=get("muMatCSZ", envir=env)
    muMatGPS=get("muMatGPS", envir=env)
    parMat=get("parMat", envir=env)
    logLiks=get("logLiks", envir=env)
    optimTables=get("optimTables", envir=env)
    lastHess=get("lastHess", envir=env)
    lastStanResults=get("lastStanResults", envir=env)
    
    state = list(muMatCSZ=muMatCSZ, muMatGPS=muMatGPS, parMat=parMat, logLiks=logLiks, 
                 optimTables=optimTables, lastHess=lastHess, lastStanResults=lastStanResults)
    subFunName = "fitModelMean"
    subFunIns = list(params=params, muVec=muVec, currIter=currIter)
    save(funIns, state, subFunName, subFunIns, file=saveFile)
    
    # base case: if past maxIter, resutrn results
    if(currIter > maxIter) {
      return(list(params=params, muVec=muVec))
    }
    print(paste0("fitting model mean, current iteration is: ", currIter))
    
    newMu = updateMu(params, muVec=muVec, fault=csz, niter=niterMCMC)
#     list(newMu=newMu, varMat=varMat, muMat=muMat, Sigmas=Sigmas, 
#          newMuSub=newMuSub, newMuPoint=newMuPoint, 
#          varMatPoint=varMatPoint, muMatPoint=muMatPoint, SigmasPoint=SigmasPoint, 
#          newMuSubPoint=newMuSubPoint, allStanResults=allStanResults)
    
    # update muMat
    muMatCSZ = get("muMatCSZ", envir=env)
    muMatGPS = get("muMatGPS", envir=env)
    lastStanResults = get("lastStanResults", envir=env)
    assign("muMatCSZ", cbind(muMatCSZ, newMu$newMu), envir=env)
    assign("muMatGPS", cbind(muMatGPS, newMu$newMuPoint), envir=env)
    assign("lastStanResults", newMu$allStanResults, envir=env)
    
    return(fitModelParams(params, c(newMu$newMuPoint, newMu$newMu), currIter=currIter+0.5))
  }
  
  # run the iterative model (either from where it left off or the beginning)
  if(!is.null(loadDat)) {
    # since user input previous progress, start from there
    
    # load the data 
    # loadDat contains: funIns, state, subFunName, and subFunIns
    # (lastIns should now be overwritten.)
    load(loadDat)
    
    # get state where last left off
    assign("muMatCSZ", state$muMatCSZ, envir=env)
    assign("muMatGPS", state$muMatGPS, envir=env)
    assign("parMat", state$parMat, envir=env)
    assign("logLiks", state$logLiks, envir=env)
    assign("optimTables", state$optimTables, envir=env)
    assign("lastHess", state$lastHess, envir=env)
    assign("lastStanResults", state$lastStanResults, envir=env)
    
    # get function inputs where last left off
    initParams=funIns$initParams
    nsim=funIns$nsim
    useMVNApprox=funIns$useMVNApprox
    gpsDat=funIns$gpsDat
    corMatGPS=funIns$corMatGPS
    muVec=funIns$muVec
#     maxIter=funIns$maxIter
    funIns$maxIter = maxIter # in case we want to add more iterations later.
    fault=funIns$fault
    niterMCMC=funIns$niterMCMC
    
    # run iterative method from where we left off
    out = do.call(subFunName, subFunIns)
  }
  else{
    # run the iterative method from the beginning
    out = fitModelParams(initParams, currIter=1)
  }
  
  # get results
  parMat = get("parMat", envir=env)
  optimTable = get("optimTables", envir=env)
  lastHess = get("lastHess", envir=env)
  logLiks = get("logLiks", envir=env)
  muMatCSZ = get("muMatCSZ", envir=env)
  muMatGPS = get("muMatGPS", envir=env)
  lastStanResults=get("lastStanResults", envir=env)
  MLEs = out$params
  muVec = out$muVec
  muVecGPS = muVec[1:nrow(gpsDat)]
  muVecCSZ = muVec[(nrow(gpsDat)+1):length(muVec)]
  
  # return results
  return(list(MLEs=MLEs, muVecGPS=muVecGPS, muVecCSZ=muVecCSZ, muMatGPS=muMatGPS, muMatCSZ=muMatCSZ, 
              parMat=parMat, logLiks=logLiks, optimTables=optimTables, hess=lastHess, 
              lastStanResults=lastStanResults))
}










