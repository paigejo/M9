# Archived functions that are no longer used, but held for reference.

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
