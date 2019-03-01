##### load in data and the necessary functions
library(fields)
setwd("~/git/M9/")
# source("loadFloodDat.R")
source("taper.R")
source("okada.R")

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
# muXi: mean of MVN in Log of xi (vector of GPS slip values).  NOTE:
# this input is no longer required.  It is estimated w/in the fxn
# muZeta: mean of GP in LOG of Zeta (constant).  If vector then must be 
#         same length as GPS data
# SigmaZeta: Covariance of LOG of zeta (of dimension |sigmaXi|)
# gpsDat: GPS data restricted to CSZ fault geometry
GPSLnLik = function(muZeta, SigmaZeta, gpsDat=slipDatCSZ, tvec = rep(1, length(sigmaXi)), 
                    normalModel=FALSE, corGPS=FALSE, doGammaSpline=FALSE, nKnotsGamma=7) {
  
  if(any(is.na(tvec)) || any(tvec < 1e-4))
    return(-10000)
  
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
      # generate conditional estimate of gamma as a function of latitude using weighted least squares on log scale
      ys = log(x)- log(muZeta) - log(tvec)
      X = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
      logModel = lm(ys~X-1, weights=ci)
      gammaEst = exp(X %*% logModel$coefficients)
      
      xCntr = x - gammaEst * tvec * muZeta
      SigmaZeta = sweep(sweep(SigmaZeta, 1, tvec * gammaEst, "*"), 2, tvec * gammaEst, "*")
    } else {
      muXi = sum((log(x)- log(muZeta) - log(tvec))*ci)
      gammaEst = exp(muXi)
      
      xCntr = x - gammaEst * tvec * muZeta
      if(any(tvec != 1))
        SigmaZeta = sweep(sweep(SigmaZeta, 1, tvec, "*"), 2, tvec, "*")
      SigmaZeta = gammaEst^2 * SigmaZeta
    }
    sigmaXi = gpsDat$slipErr
  }
  else {
    if(doGammaSpline)
      stop("fitting gamma not supported for lognormal model")
    
    muXi = sum((x- muZeta - log(tvec))*ci)
    
    xCntr = x - muZeta - muXi - log(tvec)
  }
  
  # if GPS measurement error is correlated, make it correlated in the same way as the latent process
  if(corGPS) {
    sigmaTZetas = sqrt(diag(SigmaZeta))
    CXi = sweep(sweep(SigmaZeta, 1, 1/sigmaTZetas, "*"), 2, 1/sigmaTZetas, "*")
    SigmaXi = sweep(sweep(CXi, 1, sigmaXi, "*"), 2, sigmaXi, "*")
    Sigma = SigmaZeta + SigmaXi
  }
  else
    Sigma = SigmaZeta + diag(sigmaXi^2)
  
  # get likelihood
  return(logLikGP(xCntr, chol(Sigma)))
}
# hist(solve(t(chol(Sigma))) %*% xCntr, freq=F)
# xs = seq(-4, 4, l=100)
# lines(xs, dnorm(xs),col="blue")

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
                              subDat=dr1, tvec=taper(cszDepths, lambda=lambda)) {
  # m is number of csz subfaults
  m = nrow(SigmaZetaL)
  
  # get log taper values
  logt = log(tvec)
  
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
  # uncert = subDat$Uncertainty
  
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
  likMat = apply(simEps, 2, dnorm, sd=subDat$Uncertainty)
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
                      corMatGPS=NULL, muVec=NULL, G=NULL, subDat=dr1, dStar=21000, normalizeTaper=TRUE) {
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
      initParams = c(lambdaInit, NA, sigmaZetaInit)
  }
  
  # get spatial correlation parameters (computed based on fitGPSCovariance)
  corPar = getCorPar()
  phiZeta = corPar$phiZeta
  nuZeta = corPar$nuZeta
  
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
  arealCSZCor = getArealCorMat(fault)
  corMatCSZ = arealCSZCor
  corMatCSZL = t(chol(corMatCSZ))
  coords = cbind(gpsDat$lon, gpsDat$lat)
  if(is.null(corMatGPS)) {
    corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                               onlyUpper=TRUE, smoothness=nuZeta, 
                               Distance="rdist.earth", Dist.args=list(miles=FALSE))
  }
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(csz)[,3]
  
  ##### Do optimization
  optimTable <<- NULL # table to save likelihood optimization steps
  controls = list(fnscale = -1, reltol=10^-5)
  
  # must rename variables for hessian calculation
  cszDepthsIn = cszDepthsIn
  corMatGPSIn = corMatGPS
  corMatCSZIn = corMatCSZ
  corMatCSZLIn = corMatCSZL
  phiZetaIn = phiZeta
  gpsDatIn = gpsDat
  nsimIn = nsim
  dStarIn = dStar
  normalizeTaperIn = normalizeTaper
  useMVNApproxIn = useMVNApprox
  muVecIn = muVec
  GIn = G
  subDatIn = subDat
  opt = optim(initParams, fixedDataLogLik, control=controls, hessian=FALSE, 
              cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
              phiZeta=phiZetaIn, gpsDat=gpsDatIn, nsim=nsimIn, dStar=dStarIn, normalizeTaper=normalizeTaperIn, 
              useMVNApprox=useMVNApproxIn, muZeta=muVecIn, G=GIn, subDat=subDatIn)
  
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
  hess.args = list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)
  hess = hessian(fixedDataLogLik, opt$par, 
                 cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                 phiZeta=phiZetaIn, gpsDat=gpsDatIn, nsim=nsimIn, dStar=dStarIn, normalizeTaper=normalizeTaperIn, 
                 useMVNApprox=useMVNApproxIn, muZeta=muVecIn, G=GIn, subDat=subDatIn)
  
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
fixedDataLogLik = function(params, cszDepths, gpsDat=slipDatCSZ, 
                           subDat=dr1, fault=csz, verbose=TRUE, 
                           G=NULL, nKnots=5, normalizeTaper=TRUE, dStar=21000, 
                           constrLambda=TRUE, latRange=c(40,50), 
                           distMatGPS=NULL, distMatCSZ=NULL, 
                           corGPS=FALSE, diffGPSTaper=FALSE, nKnotsGPS=nKnots, anisotropic=FALSE, 
                           squareStrikeDistCsz=NULL, squareDipDistCsz=NULL, squareStrikeDistGps=NULL, 
                           squareDipDistGps=NULL, nKnotsGamma=7, doGammaSpline=FALSE, dStarGPS=dStar, 
                           nKnotsVar=5, doVarSpline=FALSE) {
  ##### get parameters
  out = getInputPar(params, fault, gpsDat, nKnots, diffGPSTaper, nKnotsGPS, taperedGPSDat=TRUE, anisotropic, normalModel=TRUE, 
                    nKnotsVar, doVarSpline)
  muZeta = out$muZeta
  muVecCSZ = out$muVecCSZ
  muVecGPS = out$muVecGPS
  sigmaZeta = out$sigmaZeta
  splinePar = out$taperPar
  splineParGPS = out$taperParGPS
  phiZeta = out$phiZeta
  nuZeta = out$nuZeta
  lambda0 = out$lambda0
  alpha = out$alpha
  varPar = out$varPar
  parNames = out$parNames
  
  muZetaGPS = muZeta
  muZetaCSZ = muZeta
  vecMu=FALSE
  nuZeta = 3/2
  
  # generate taper lambda parameter as a function of latitude for the subsidence data
  tvec = getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar, latRange=latRange, fault=fault)
  Xi = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambda = Xi %*% splinePar # lambda is only used for checking if lambda < 0
  
  # do the same for the gps data
  XiGPS1 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
  lambda = XiGPS1 %*% splinePar
  if(diffGPSTaper) {
    XiGPS2 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
    lambda = lambda - XiGPS2 %*% splineParGPS
  }
  
  # calculate latent field standard deviation as a function of latitude for the subsidence data
  if(doVarSpline) {
    XiVar = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
    sigmaZeta = XiVar %*% varPar
    
    XiVarGPS = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsVar, latRange=latRange)
    sigmaZetaGPS = XiVarGPS %*% varPar
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
  lambda0 = 0.25
  
  ##### set up column names for parameter optimization table
  cNames = c(parNames, "LL", "subLL", "GPSLL", "priorLL", "LLSE")
  
  ##### compute unit slip Okada seadef if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=lambda0)
  }
  
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
  if(constrLambda && any(lambdaSeq < 0) && any(lambdaSeqGPS < 0))
    lambdaInRange = FALSE
  else if(any(abs(lambdaSeq) > 15) && any(abs(lambdaSeqGPS) > 15))
    lambdaInRange=FALSE
  
  if(!lambdaInRange || any(sigmaZeta <= 0) || (muZeta <= 0) || (phiZeta <= 0) || (alpha <= 0.05) || (alpha >= 20)) {
    newRow = rbind(c(params, -10000, -5000, -5000, 0, NA))
    colnames(newRow) = cNames
    if(verbose)
      print(newRow)
    optimTable <<- rbind(optimTable, newRow)
    return(-10000) # return a very small likelihood
  }
  
  ##### compute covariance of LOG of Zeta and its Cholesky decomp
  if(doVarSpline) {
    SigmaZetaGPS = sweep(sweep(corMatGPS, 1, sigmaZeta, "*"), 2, sigmaZeta, "*")
    SigmaZetaCSZL = sweep(corMatCSZL, 1, sigmaZetaGPS, "*")
  } else {
    SigmaZetaGPS = sigmaZeta^2 * corMatGPS
    SigmaZetaCSZL = sigmaZeta * corMatCSZL
  }
  
  ##### Compute Likelihood of subsidence and GPS data
  SigmaZeta = SigmaZetaCSZL %*% t(SigmaZetaCSZL)
  sub = subsidenceLnLikMod2(muZetaCSZ, lambda, sigmaZeta, SigmaZeta, G, subDat=subDat, 
                            tvec=tvec, dStar=dStar, normalizeTaper=normalizeTaper, normalModel=TRUE)
  
  XiGPS1 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  lambda = XiGPS1 %*% splinePar
  if(diffGPSTaper) {
    XiGPS2 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
    lambda = lambda - XiGPS2 %*% splineParGPS
  }
  tvecGPS = taper(gpsDat$Depth, lambda, alpha=2, dStar=dStarGPS, normalize=normalizeTaper)
  
  GPS = GPSLnLik(muZetaGPS, SigmaZetaGPS, gpsDat, normalModel=TRUE, tvec=tvecGPS, corGPS=corGPS, 
                 doGammaSpline=doGammaSpline, nKnotsGamma=nKnotsGamma)
  lnLik = sub[1] + GPS
  lnLikSE = sub[2]
  priorLnLik = 0
  
  ##### update parameter optimization table
  newRow = c(params, lnLik, sub[1], GPS, priorLnLik, lnLikSE)
  newRow = rbind(newRow)
  colnames(newRow) = cNames
  if(verbose)
    print(newRow, digits=5)
  optimTable <<- rbind(optimTable, newRow)
  
  ##### return full data likelihood
  lnLik
}
# jacobian(fixedDataLogLik, params, cszDepths=cszDepths, corMatGPS=corMatGPS, corMatCSZL=corMatCSZL, phiZeta=phiZeta, 
# useMVNApprox=TRUE, G=G, nKnots=nKnots, dStar=dStar, useSubPrior=useSubPrior, useSlipPrior=useSlipPrior, fauxG=fauxG)
# starting quantile values:
# Browse[2]> getQS(params)
# [1] 238.1374
# Browse[2]> getQY(params, Xi=getSplineBasis(), G, fauxG)
# [1] 36.08274
# pars = splineFit21k5$MLEs
# pars = c(pars[2], pars[3], pars[6:length(pars)])
# getQS(pars)
# getQY(pars, getSplineBasis(), G, fauxG)

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
# here muZeta and sigmaZeta can be either a constant or a vector of length ncol(G)
subsidenceLnLikMod2 = function(muZeta, lambda, sigmaZeta, SigmaZeta, G, subDat=dr1, 
                               tvec=taper(csz$depth, lambda=lambda, normalize=normalizeTaper), 
                               dStar=21000, normalizeTaper=TRUE, normalModel=FALSE) {
  
  # if sigmaZeta gets too large, variance explodes, getting lnLik is numerically infeasible
  if(!normalModel && (any(sigmaZeta > 4) || any(sigmaZeta < 0)))
    return(c(lnLik=-5000, lnLikSE=NA))
  
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, SigmaZeta, G, subDat=subDat, 
                                   tvec=tvec, normalModel=normalModel)
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
                         nsim=500, useMVNApprox=TRUE, subDat=dr1) {
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
  
  # get spatial correlation parameters (computed based on fitGPSCovariance)
  corPar = getCorPar()
  phiZeta = corPar$phiZeta
  nuZeta = corPar$nuZeta
  
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
  arealCSZCor = getArealCorMat(fault)
  CorB = arealCSZCor
  CorBL = t(chol(CorB))
  
  # precompute Okada matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  GFull = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  
  ##### Do cross-validation
  # scramle the data (and save how to unscramble it)
  scrambleInds = sample(1:nrow(subDat), nrow(subDat))
  subDatScramble = subDat[scrambleInds,]
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
      predInds = ((1:nrow(subDat)) <= nrow(subDat)*j/nfold) & 
        (((1:nrow(subDat)) > nrow(subDat)*(j-1)/nfold))
      fitInds = !predInds
      scrambledFitInds = scrambleInds[fitInds]
      thisEvents = events[scrambledFitInds]
      ord = order(factor(thisEvents, levels=uniqueEvents))
      YFitOrder = scrambledFitInds[ord]
      subDatFit = subDat[YFitOrder,]
      
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
      subDatPred = subDat[YPredOrder,]
      
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
      subDistn = getSubsidenceVarianceMat(params, fault=csz, G=GFull, subDat=subDat)
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
  # get spatial correlation parameters (computed based on fitGPSCovariance)
  corPar = getCorPar()
  phiZeta = corPar$phiZeta
  nuZeta = corPar$nuZeta
  
  # get data
  logX = log(gpsDat$slip)
  YPred = -subDatPred$subsidence
  YFit = -subDatFit$subsidence
  
  # compute inflated variance sigmaXi. Derivation from week 08_30_17.Rmd presentation.  
  # Transformation from additive error to multiplicative lognormal model with asympototic 
  # median and variance matching.
  sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  
  # load covariance matrix of sigma over CSZ fault cells
  arealCSZCor = getArealCorMat(fault)
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
                             loadDat=NULL, saveFile="iterFitProgress.RData", usePrior=FALSE, 
                             subDat=dr1) {
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
  
  # get Okada linear transformation matrix
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  
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
    
    # base case: if past maxIter, return results
    if(currIter > maxIter)
      return(list(params=params, muVec=muVec))
    
    print(paste0("fitting model parameters, current iteration is: ", currIter))
    
    # fit the parameters
    parFit = doFixedFit(params, nsim, useMVNApprox, gpsDat, corMatGPS, muVec, G=G, subDat=subDat)
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
    
    newMu = updateMu(params, muVec=muVec, fault=csz, niter=niterMCMC, G=G, usePrior=usePrior, subDat=subDat)
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


############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions for fitting model assuming spline basis for lambda in latitude

##### function for performing full MLE fit of model
# if muVec not included, params are of form:
# params[1]: lambda (NA, use spline basis instead)
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: lambda0 (0.25)
# params[5]: muXi
# params[6:(6+nKnots-1)]: lambda spline coefficient basis coefficients
# otherwise, they are of form:
# params[1]: lambda
# params[2]: sigmaZeta
# params[3]: lambda0
# initParams is of the form: c(muZeta, sigmaZeta, splineParams)
# NOTE: fixing params[6:(6+nKnots-1)] = lambda for one lambda corresponds to 
#       the previous model, lambda is constant throughout space
# NOTE2: only set useGrad to TRUE when using a MVN approximation (when useMVNApprox=TRUE)
doFitSpline = function(initParams=NULL, nsim=500, useMVNApprox=FALSE, gpsDat=slipDatCSZ, 
                      corMatGPS=NULL, muVec=NULL, G=NULL, fauxG=NULL, subDat=dr1, fault=csz, nKnots=5, normalizeTaper=TRUE, 
                      dStar=21000, useASLApprox=FALSE, useGrad=FALSE, maxit=500, useSubPrior=FALSE, 
                      useSlipPrior=TRUE, constrLambda=TRUE, latRange=c(40, 50), normalModel=FALSE) {
  
  # load the precomputed correlation matrix
  arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  corMatCSZ = arealCSZCor
  
  ##### Set up initial parameter guesses for the MLE optimization if necessary
  if(is.null(initParams)) {
    lambdaInit=2
    if(!normalModel)
      muZetaInit = log(20)
    else
      muZetaInit = 20
    # variance of a lognormal (x) is (exp(sigma^2) - 1) * exp(2mu + sigma^2)
    # then sigma^2 = log(var(e^x)) - 2mu - log(e^(sigma^2) - 1)
    # so sigma^2 \approx log(var(e^x)) - 2mu
    # sigmaZetaInit = sqrt(log(var(dr1$subsidence)) - 2*muZetaInit)
    if(!normalModel)
      sigmaZetaInit = 1
    else
      sigmaZetaInit = sd(gpsDat$slip)*(muZetaInit/mean(gpsDat$slip))
    # muXiInit = log(2.25) # based on plots from exploratory analysis
    # initSplineParams = getInitialSplineEsts(muZetaInit, sigmaZetaInit, lambdaInit, nKnots=nKnots, 
                                            # normalizeTaper=normalizeTaper, dStar=dStar)
    initSplineParams = getSplineEstsMomentMatch(muZetaInit, sigmaZetaInit, lambdaInit, corMatCSZ, nKnotsMean=nKnots, 
                                                nKnotsVar=nKnots, normalizeTaper=normalizeTaper, dStar=dStar, 
                                                normalModel=normalModel, useGrad=useGrad, fault=fault, 
                                                latRange=latRange, G=G, subDat=subDat)$betaHat
    
    initParams = c(muZetaInit, sigmaZetaInit, initSplineParams)
  }
  
  # get spatial correlation parameters (computed based on fitGPSCovariance)
  corPar = getCorPar(normalModel=normalModel)
  phiZeta = corPar$phiZeta
  nuZeta = corPar$nuZeta
  
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
  
  # compute correlation matrix Cholesky decomposition
  corMatCSZL = t(chol(corMatCSZ))
  coords = cbind(gpsDat$lon, gpsDat$lat)
  if(is.null(corMatGPS)) {
    corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                               onlyUpper=TRUE, smoothness=nuZeta, 
                               Distance="rdist.earth", Dist.args=list(miles=FALSE))
  }
  
  # compute correlation matrix for GPS locations if necessary (only necessary for gradient)
  if(useGrad) {
    coords = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
    corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                               onlyUpper=FALSE, smoothness=nuZeta, 
                               Distance="rdist.earth", Dist.args=list(miles=FALSE))
  }
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  if(is.null(fauxG)) {
    if(useSubPrior)
      fauxG = getFauxG()
    else
      fauxG = NULL
  }
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(csz)[,3]
  
  ##### Do optimization
  optimTable <<- NULL # table to save likelihood optimization steps
  controls = list(fnscale = -1, reltol=1e-7, parscale=rep(.1, length(initParams)), maxit=maxit)
  
  cszDepthsIn = cszDepths
  corMatGPSIn = corMatGPS
  corMatCSZIn = corMatCSZ
  corMatCSZLIn = corMatCSZL
  phiZetaIn = phiZeta
  gpsDatIn = gpsDat
  nsimIn = nsim
  dStarIn = dStar
  normalizeTaperIn = normalizeTaper
  useMVNApproxIn = useMVNApprox
  useASLApproxIn = useASLApprox
  muVecIn = muVec
  nKnotsIn = nKnots
  GIn = G
  subDatIn = subDat
  faultIn = fault
  useSubPriorIn = useSubPrior
  useSlipPriorIn = useSlipPrior
  fauxGIn = fauxG
  constrLambdaIn = constrLambda
  latRangeIn = latRange
  normalModelIn = normalModel
  if(!useGrad)
    opt = optim(initParams, fixedDataLogLik, control=controls, hessian=FALSE, 
                cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                phiZeta=phiZetaIn, gpsDat=gpsDatIn, nsim=nsimIn, nKnots=nKnotsIn, 
                useMVNApprox=useMVNApproxIn, muZeta=muVecIn, G=GIn, subDat=subDatIn, fault=faultIn, 
                normalizeTaper=normalizeTaperIn, dStar=dStarIn, useASLApprox=useASLApproxIn, 
                useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn)
  else {
    # opt = optim(initParams, fixedDataLogLik, fixedDataLogLikGrad, control=controls, hessian=FALSE, 
    #             cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
    #             phiZeta=phiZetaIn, gpsDat=gpsDatIn, nKnots=nKnotsIn, 
    #             useMVNApprox=TRUE, muZeta=muVecIn, G=GIn, subDat=subDatIn, 
    #             normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
    #             method="BFGS")
    
    if(constrLambda) {
      latSeq = seq(latRange[1], latRange[2], l=100)
      U = cbind(0, 0, getSplineBasis(nKnots=nKnots, lats=latSeq, latRange=latRange))
      C = rep(0, nrow(U))
      opt = constrOptim(initParams, fixedDataLogLik, fixedDataLogLikGrad, control=controls, hessian=TRUE, 
                        cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                        phiZeta=phiZetaIn, gpsDat=gpsDatIn, nKnots=nKnotsIn, 
                        useMVNApprox=TRUE, muZeta=muVecIn, G=GIn, subDat=subDatIn, 
                        fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
                        useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                        constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
                        method="BFGS", ui=U, ci=C)
    }
    else {
      opt = optim(initParams, fixedDataLogLik, fixedDataLogLikGrad, control=controls, hessian=TRUE, 
                  cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                  phiZeta=phiZetaIn, gpsDat=gpsDatIn, nKnots=nKnotsIn, 
                  useMVNApprox=TRUE, muZeta=muVecIn, G=GIn, subDat=subDatIn, 
                  fault=faultIn, normalizeTaper=normalizeTaperIn, dStar=dStarIn, 
                  useSubPrior=useSubPriorIn, useSlipPrior=useSlipPriorIn, fauxG=fauxGIn, 
                  constrLambda=constrLambdaIn, latRange=latRangeIn, normalModel=normalModelIn, 
                  method="BFGS")
    }
  }
  
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
  if(length(opt$par) > 3 || (nKnots==1)) {
    splinePar = opt$par[(length(opt$par)-(nKnots-1)):length(opt$par)]
  }
  logLikMLE = opt$value
  # hess.args = list(eps=1e-5, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)
  # hess = hessian(fixedDataLogLik, opt$par, method.args=hess.args, 
  #                cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
  #                phiZeta=phiZetaIn, gpsDat=gpsDatIn, nsim=nsimIn, nKnots=nKnotsIn, 
  #                useMVNApprox=useMVNApproxIn, muZeta=muVecIn, G=GIn, subDat=subDatIn, 
  #                normalizeTaper=normalizeTaperIn, dStar=dStarIn, useASLApprox=useASLApproxIn)
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
  if(!normalModel)
    muXiMLE = sum((logX-muZetaGPS)*ci)
  if(normalModel) {
    muXiMLE = sum((logX-log(muZetaGPS))*ci)
    muXiMLE = exp(muXiMLE)
  }
  
  # params is in order: lambda, muZeta, sigmaZeta, lambda0, muXi, splineParams
  if(is.null(muVec))
    MLEs = c(NA, opt$par[1:2], 0.25, muXiMLE, opt$par[-c(1:2)])
  else
    MLEs = c(NA, NA, opt$par[1], 0.25, muXiMLE, opt$par[-1])
  
  # Return results
  list(MLEs=MLEs, lambdaMLE=lambdaMLE, muZetaMLE=muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambda0MLE=NA, 
       muXiMLE=muXiMLE, logLikMLE=logLikMLE, splineParMLE=splinePar, hess=hess, optimTable=optimTable, 
       tvec=getTaperSpline(splinePar, fault=fault, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar, latRange=latRange))
}

# Note that although lambda estimates can be negative, technically since it is 
# squared in the double exponential that is not a serious problem.
getInitialSplineEsts = function(muVecCSZ, sigmaZeta, lambda, G=NULL, nKnots=5, subDat=dr1, 
                                normalizeTaper=TRUE, dStar=21000, initPar = c(lambda, rep(0, nKnots))) {
  # compute spline basis coefficients for subsidence data
  latRange = c(40, 50)
  Xi = bs(subDat$Lat, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  mod = lm(subDat$subsidence ~ Xi - 1)
  preds = Xi %*% coef(mod)
  
  if(length(muVecCSZ)==1)
    muVecCSZ = rep(muVecCSZ, nrow(csz))
  arealCSZCor = getArealCorMat(fault)
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
  
  # compute Okada matrix if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  
  # recompute spline basis over csz fault geometry
  Xi = bs(csz$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # compute initial estimates for spline basis coefficients for CSZ taper
  # Xi betaHat = G diag(E[zeta]) Xi betaPrimeHat
  y = -preds
  Xmat = G %*% diag(expectZeta)
  
  lats = seq(40, 50, l=100)
  XiLat = bs(lats, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  XiLat = cbind(rep(1, nrow(XiLat)), XiLat)
  sqErr = function(params) {
    lambdaLats = XiLat %*% params
    lambdaCSZ = Xi %*% params
    if(any(c(lambdaLats, lambdaCSZ) < 0))
      return(1e15)
    mean(((y - Xmat %*% taper(csz$depth, lambda=Xi%*%params, normalize=normalizeTaper, dStar=dStar))/subDat$Uncertainty)^2)
  }
  opt = optim(initPar, sqErr)
  betaPrimeHat = opt$par
  
  lambdaVecCSZ = Xi %*% betaPrimeHat
  lambdaVecLat = XiLat %*% betaPrimeHat
  
  return(list(betaHat = betaPrimeHat, lambdaVecCSZ=lambdaVecCSZ, lats=lats, lambdaVecLat=lambdaVecLat))
}

# this function is like getInitialSplineEsts except it matches mean as well 
# as variance using spline estimates of fitSplinesToSubDat.  Here we assume 
# that there can be a different number of knots for mean and variance, but 
# the taper function must have the same number of knots as the mean spline 
# fit.  If initPar has length longer than nKnotsMean+1, it is assumed that 
# the user included guesses for muZeta and sigmaZeta as well.
getSplineEstsMomentMatch = function(muVecCSZ, sigmaZeta, lambda, corMatCSZ=NULL, splineParMean=NULL, 
                                    splineParVar=NULL, G=NULL, nKnotsMean=5, nKnotsVar=5, 
                                    subDat=dr1, normalizeTaper=TRUE, dStar=21000, 
                                    initPar = c(lambda, rep(.01, nKnotsMean-1)), useGrad=TRUE, normalModel=FALSE, 
                                    fault=csz, latRange=c(40,50)) {
  
  if(useGrad && normalModel)
    useGrad = FALSE
  
  if(length(initPar) > nKnotsMean)
    pickMuSigma = TRUE
  else
    pickMuSigma = FALSE
  
  # first compute mean and variance spline parameters, if necessary
  if(is.null(splineParMean) || is.null(splineParVar)) {
    splineVals = fitSplinesToSubDat(nKnotsMean=nKnotsMean, nKnotsVar=nKnotsVar, subDat=subDat)
    splineParMean = splineVals$meanSplinePar
    splineParVar = splineVals$varSplinePar
  }
  
  # compute spline mean
  XiMean = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange, lats=subDat$Lat)
  meanVec = XiMean %*% splineParMean # meanVec is average subsidence (not uplift) vs lat
  
  # compute spline variance
  XiVar = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange, lats=subDat$Lat)
  varVec = XiVar %*% splineParVar
  
  # get zeta mean vector
  if(length(muVecCSZ)==1)
    muVecCSZ = rep(muVecCSZ, nrow(fault))
  arealCSZCor = getArealCorMat(fault)
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  if(!normalModel)
    expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
  else
    expectZeta = muVecCSZ
  
  # get zeta covariance matrix (covariance of zeta, not log(zeta))
  if(!normalModel) {
    covZeta = exp(covMatCSZ) - 1
    covZeta = sweep(covZeta, 1, expectZeta, "*")
    covZeta = sweep(covZeta, 2, expectZeta, "*")
  }
  else
    covZeta = covMatCSZ
  
  # compute correlation matrix (of zeta, not log(zeta))
  if(pickMuSigma) {
    sigmas = sqrt(diag(covZeta))
    corZeta = sweep(sweep(covZeta, 1, sigmas, "/"), 2, sigmas, "/")
  }
  
  # compute Okada matrix if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  
  # compute estimates for spline basis coefficients for CSZ taper
  # Xi betaHat = G diag(E[zeta]) taper(Xi betaPrimeHat)
  # Xmat is G diag(E[zeta])
  y = -meanVec # y is average uplift (not subsidence) vs latitude fit to sub dat
  Xmat = G %*% diag(expectZeta) # multiply this by tvec to get expected uplift
  
  # get basis matrix for the taper parameter over CSZ geometry and latitudes
  Xi = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange, lats=fault$latitude)
  lats = seq(latRange[1], latRange[2], l=100)
  XiLat = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange, lats=lats)
  
  # function measuring how well the moments are matched by the 
  # taper function parameters.  Only checks how well it matches at the data locations.
  pow=2 # may change later.  Default is 2
  sqErr = function(params, makeLambdaPos=FALSE) {
    if(pickMuSigma) {
      mu = params[1]
      sigma = params[2]
      taperPar = params[3:length(params)]
    }
    else
      taperPar = params
    
    # compute taper parameters over CSZ grid
    lambdaLats = XiLat %*% taperPar
    lambdaCSZ = Xi %*% taperPar
    
    if(makeLambdaPos && any(c(lambdaLats, lambdaCSZ) < 0))
      return(1e15)
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaCSZ, normalize=normalizeTaper, dStar=dStar)
    
    # if necessary, update covZeta and Xmat (variance and expectation of true field)
    if(pickMuSigma) {
      if(!normalModel)
        covZeta = corZeta * (exp(sigma^2) - 1) * exp(2*mu + sigma^2/2)
      else
        covZeta = corZeta * sigma^2
      covMatCSZ = sigma^2 * corMatCSZ
      if(!normalModel)
        expectZeta = exp(mu + diag(covMatCSZ)/2)
      else
        expectZeta = rep(mu, nrow(covMatCSZ))
      Xmat = G %*% diag(expectZeta)
    }
    
    # check how well we match the spline variance estimate
    GT = sweep(G, 2, tvec, "*")
    VSubsidence = GT %*% covZeta %*% t(GT)
    if(pow == 2)
      varMatch = mean((varVec - diag(VSubsidence))^pow)/(2*mean(varVec^2))
    else
      varMatch = mean((varVec - diag(VSubsidence))^pow)
      
    # check how well we match the spline mean estimate
    meanMatch = mean(((y - Xmat %*% tvec)/sqrt(varVec))^pow)
    
    meanMatch + varMatch
  }
  sqErrGrad = function(params) {
    if(pickMuSigma) {
      mu = params[1]
      sigma = params[2]
      taperPar = params[3:length(params)]
    }
    else
      taperPar = params
    
    # compute taper parameters over CSZ grid
    lambdaLats = XiLat %*% taperPar
    lambdaCSZ = Xi %*% taperPar
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaCSZ, normalize=normalizeTaper, dStar=dStar)
    
    # if necessary, update covZeta and Xmat (variance and expectation of true field)
    # NOTE: corZeta and covZeta or the correlation and covariance of zeta, not log(zeta)
    if(pickMuSigma) {
      covZeta = corZeta * (exp(sigma^2) - 1) * exp(2*mu + sigma^2/2)
      covMatCSZ = sigma^2 * corMatCSZ
      expectZeta = exp(mu + diag(covMatCSZ)/2)
      Xmat = G %*% diag(expectZeta)
    }
    
    # get expected uplift
    expectY = Xmat %*% tvec
    
    # get subsidence variance
    GT = sweep(G, 2, tvec, "*")
    VarVecY = c(diag(GT %*% covZeta %*% t(GT)))
    
    if(!pickMuSigma) {
      MG = MGrad(expectZeta, expectY, lambdaCSZ, Xi, G, y, varVec, fault=csz, 
                 normalizeTaper=normalizeTaper, dStar=dStar, pow=pow)
      VG = VGrad(covZeta, VarVecY, lambdaCSZ, tvec, Xi, G, varVec, fault=csz, 
                 subDat=subDat, normalizeTaper=normalizeTaper, dStar=dStar, pow=pow)
    }
    else {
      MG = MGrad(expectZeta, expectY, lambdaCSZ, Xi, G, y, varVec, fault=csz, 
                 normalizeTaper=normalizeTaper, dStar=dStar, pow=pow, muZeta=mu, 
                 sigmaZeta=sigma)
      VG = VGrad(covZeta, VarVecY, lambdaCSZ, tvec, Xi, G, varVec, fault=csz, 
                 subDat=subDat, normalizeTaper=normalizeTaper, dStar=dStar, pow=pow, 
                 muZeta=mu, sigmaZeta=sigma, corZeta=corZeta)
    }
    
    c(MG + VG)
  }
  # use robust Nelder-Mead then use use BFGS for gradient optimization
  controls = list(reltol=1e-7)
  opt = optim(initPar, sqErr, control=controls, method="Nelder-Mead")
  if(useGrad)
    opt = optim(opt$par, sqErr, sqErrGrad, control=controls, method="BFGS")
  
  if(!pickMuSigma) {
    muZeta = muVecCSZ
    betaPrimeHat = opt$par
  }
  else {
    muZeta = opt$par[1]
    sigmaZeta = opt$par[2]
    betaPrimeHat = opt$par[3:length(opt$par)]
  }
  lambdaVecCSZ = Xi %*% betaPrimeHat
  lambdaVecLat = XiLat %*% betaPrimeHat
  
  return(list(betaHat = betaPrimeHat, lambdaVecCSZ=lambdaVecCSZ, lats=lats, lambdaVecLat=lambdaVecLat, 
              muZeta=muZeta, sigmaZeta=sigmaZeta))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# function for estimating subsidence likelihood using copula and asymmetric 
# shifted laplace distribution
# Approximate the distribution of G %*% T %*% Zeta with a MVN to get likelihood
# here muZeta can be either a constant or a vector of length ncol(G)
subsidenceLnLikASL = function(muZeta, lambda, sigmaZeta, SigmaZeta, G, subDat=dr1, 
                               tvec=taper(csz$depth, lambda=lambda), nsim=10000) {
  # .9578
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, SigmaZeta, G, subDat=subDat, 
                                   tvec=tvec)
  mu = mvnApprox$mu
  Sigma = mvnApprox$Sigma
  
  # add data noise to covariance matrix
  diag(Sigma) = diag(Sigma) + subDat$Uncertainty^2
  
  # generate 10000 subsidence simulations
  params=c(lambda, muZeta, sigmaZeta, 0.25, NA)
  reals = preds(params, nsim=nsim, fault=csz, tvec=tvec)
  subs = predsToSubsidence(params, reals, useMVNApprox=FALSE, G=G, subDat=subDat)
  subSims = subs$noiseSims
  
  # compute the parameters of the ASL
  parASL = approxASL(subSims, breaks=1000)
  mHats = parASL$mHats
  lambdaHats = parASL$lambdaHats
  kappaHats = parASL$kappaHats
  
  # backtransform data to multivariate uniform
  unifSims = pASL(-subDat$subsidence, mHats, lambdaHats, kappaHats)
  
  # transform data to have correct marginal Gaussian distribution
  sds = sqrt(diag(Sigma))
  normSims = qnorm(unifSims, mean=mu, sd=sds)
  
  # get log likelihood
  y = normSims # MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  lnLik = logLikGP(y - mu, chol(Sigma))
  
  if(is.nan(lnLik)) {
    print(NaN)
  }
  
  return(c(lnLik=lnLik, lnLikSE=0))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions for performing likelihood coordinate ascent

#
# Note that although lambda estimates can be negative, technically since it is 
# squared in the double exponential that is not a serious problem.
splineAscent = function(muVecCSZ, sigmaZeta, lambda, G=NULL, nKnots=5, subDat=dr1, 
                        normalizeTaper=TRUE, dStar=21000, initPar = c(lambda, rep(0, nKnots))) {
  # compute spline basis coefficients for subsidence data
  latRange = c(40, 50)
  Xi = bs(subDat$Lat, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  mod = lm(subDat$subsidence ~ Xi - 1)
  preds = Xi %*% coef(mod)
  
  if(length(muVecCSZ)==1)
    muVecCSZ = rep(muVecCSZ, nrow(csz))
  arealCSZCor = getArealCorMat(fault)
  corMatCSZ = arealCSZCor
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  expectZeta = exp(muVecCSZ + diag(covMatCSZ)/2)
  
  # compute Okada matrix if necessary
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  
  # recompute spline basis over csz fault geometry
  Xi = bs(csz$latitude, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  Xi = cbind(rep(1, nrow(Xi)), Xi)
  
  # compute initial estimates for spline basis coefficients for CSZ taper
  # Xi betaHat = G diag(E[zeta]) Xi betaPrimeHat
  y = -preds
  Xmat = G %*% diag(expectZeta)
  
  lats = seq(40, 50, l=100)
  XiLat = bs(lats, df=nKnots-1, intercept=FALSE, Boundary.knots=latRange)
  XiLat = cbind(rep(1, nrow(XiLat)), XiLat)
  sqErr = function(params) {
    lambdaLats = XiLat %*% params
    lambdaCSZ = Xi %*% params
    if(any(c(lambdaLats, lambdaCSZ) < 0))
      return(1e15)
    mean(((y - Xmat %*% taper(csz$depth, lambda=Xi%*%params, normalize=normalizeTaper, dStar=dStar))/subDat$Uncertainty)^2)
  }
  opt = optim(initPar, sqErr)
  betaPrimeHat = opt$par
  
  lambdaVecCSZ = Xi %*% betaPrimeHat
  lambdaVecLat = XiLat %*% betaPrimeHat
  
  return(list(betaHat = betaPrimeHat, lambdaVecCSZ=lambdaVecCSZ, lats=lats, lambdaVecLat=lambdaVecLat))
}

##### function for performing full MLE fit of model using COORDINATE ASCENT
# if muVec not included, params are of form:
# params[1]: lambda (NA, use spline basis instead)
# params[2]: muZeta
# params[3]: sigmaZeta
# params[4]: lambda0 (0.25)
# params[5]: muXi
# params[6:(6+nKnots-1)]: lambda spline coefficient basis coefficients
# otherwise, they are of form:
# params[1]: lambda
# params[2]: sigmaZeta
# params[3]: lambda0
# initParams is of the form: c(muZeta, sigmaZeta, splineParams)
# NOTE: fixing params[6:(6+nKnots-1)] = lambda for one lambda corresponds to 
#       the previous model, lambda is constant throughout space
doFitSplineCoord = function(initParams=NULL, nsim=500, useMVNApprox=FALSE, gpsDat=slipDatCSZ, 
                            corMatGPS=NULL, muVec=NULL, G=NULL, subDat=dr1, nKnots=5, normalizeTaper=TRUE, 
                            dStar=21000, useASLApprox=FALSE, latRange=c(40,50)) {
  # load the precomputed correlation matrix
  arealCSZCor = getArealCorMat(fault)
  corMatCSZ = arealCSZCor
  
  ##### first get spline estimates of mean and variance in subsidence data space
  subSplineFit = fitSplinesToSubDat(nKnotsMean=nKnots, nKnotsVar=nKnots)
  meanSplinePar = subSplineFit$meanSplinePar
  varSplinePar = subSplineFit$varSplinePar
  
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
    initSplineParams = getSplineEstsMomentMatch(muZetaInit, sigmaZetaInit, lambdaInit, corMatCSZ, nKnotsMean=nKnots, 
                                                nKnotsVar=nKnots, normalizeTaper=normalizeTaper, dStar=dStar, 
                                                splineParMean=meanSplinePar, splineParVar=varSplinePar)$betaHat
    initParams = c(muZetaInit, sigmaZetaInit, initSplineParams)
  }
  else{
    initSplineParams = initParams[3:length(initParams)]
  }
  
  # get spatial correlation parameters (computed based on fitGPSCovariance)
  corPar = getCorPar()
  phiZeta = corPar$phiZeta
  nuZeta = corPar$nuZeta
  
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
  
  # compute Cholesky decomposition of corMatCSZ
  corMatCSZL = t(chol(corMatCSZ))
  coords = cbind(gpsDat$lon, gpsDat$lat)
  if(is.null(corMatGPS)) {
    corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta, 
                               onlyUpper=TRUE, smoothness=nuZeta, 
                               Distance="rdist.earth", Dist.args=list(miles=FALSE))
  }
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(csz, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(csz)[,3]
  
  ##### make wrapper around optimization function to set the spline parameters
  # input parameters: 
  # muZeta, sigmaZeta
  #
  splinePar = initSplineParams
  splineLogLik = function(par) {
    # extract relevant parameters
    muZeta = par[1]
    sigmaZeta = par[2]
    
    # do crude estimation of taper parameters (use last estimate as initial guess 
    # for minimizer)
    splinePar = getSplineEstsMomentMatch(muZeta, sigmaZeta, NULL, corMatCSZ, nKnotsMean=nKnots, 
                                         nKnotsVar=nKnots, normalizeTaper=normalizeTaper, 
                                         dStar=dStar, initPar = splinePar, G=G, 
                                         splineParMean=meanSplinePar, splineParVar=varSplinePar)$betaHat
    
    # now get the likelihood given the estimate of the splinePar
    allPar = c(par, splinePar)
    fixedDataLogLik(allPar, cszDepths=cszDepths, corMatGPS=corMatGPS, 
                    corMatCSZL=corMatCSZL, phiZeta=phiZeta, gpsDat=gpsDat, 
                    nsim=nsim, nKnots=nKnots, useMVNApprox=useMVNApprox, 
                    muZeta=muVec, G=G, subDat=subDat, 
                    normalizeTaper=normalizeTaper, dStar=dStar, 
                    useASLApprox=useASLApprox)
  }
  
  ##### Do initial optimization
  print("performing initial optimization with coordinate descent")
  optimTable <<- NULL # table to save likelihood optimization steps
  controls = list(fnscale = -1, reltol=1e-7)
  initOpt = optim(initParams[1:2], splineLogLik, control=controls, hessian=FALSE)
  initMLEs = initOpt$par
  initSplineEsts = getSplineEstsMomentMatch(initMLEs[1], initMLEs[2], NULL, corMatCSZ, nKnotsMean=nKnots, 
                                            nKnotsVar=nKnots, normalizeTaper=normalizeTaper, dStar=dStar, 
                                            initPar = splinePar, G=G, 
                                            splineParMean=meanSplinePar, splineParVar=varSplinePar)$betaHat
  
  ##### now use that output from that optimization as initialization
  initParams = c(initMLEs, initSplineEsts)
  
  cszDepthsIn = cszDepths
  corMatGPSIn = corMatGPS
  corMatCSZIn = corMatCSZ
  corMatCSZLIn = corMatCSZL
  phiZetaIn = phiZeta
  gpsDatIn = gpsDat
  nsimIn = nsim
  dStarIn = dStar
  normalizeTaperIn = normalizeTaper
  useMVNApproxIn = useMVNApprox
  useASLApproxIn = useASLApprox
  muVecIn = muVec
  nKnotsIn = nKnots
  GIn = G
  subDatIn = subDat
  print("beginning final optimization")
  opt = optim(initParams, fixedDataLogLik, control=controls, hessian=FALSE, 
              cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
              phiZeta=phiZetaIn, gpsDat=gpsDatIn, nsim=nsimIn, nKnots=nKnotsIn, 
              useMVNApprox=useMVNApproxIn, muZeta=muVecIn, G=GIn, subDat=subDatIn, 
              normalizeTaper=normalizeTaperIn, dStar=dStarIn, useASLApprox=useASLApproxIn)
  
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
  if(length(opt$par) > 3) {
    splinePar = opt$par[(length(opt$par)-(nKnots-1)):length(opt$par)]
  }
  logLikMLE = opt$value
  print("calculating hessian")
  hess.args = list(eps=1e-5, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)
  hess = hessian(fixedDataLogLik, opt$par, method.args=hess.args, 
                 cszDepths=cszDepthsIn, corMatGPS=corMatGPSIn, corMatCSZL=corMatCSZLIn, 
                 phiZeta=phiZetaIn, gpsDat=gpsDatIn, nsim=nsimIn, nKnots=nKnotsIn, 
                 useMVNApprox=useMVNApproxIn, muZeta=muVecIn, G=GIn, subDat=subDatIn, 
                 normalizeTaper=normalizeTaperIn, dStar=dStarIn, useASLApprox=useASLApproxIn)
  
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
  
  # params is in order: lambda, muZeta, sigmaZeta, lambda0, muXi, splineParams
  if(is.null(muVec))
    MLEs = c(NA, opt$par[1:2], 0.25, muXiMLE, opt$par[-c(1:2)])
  else
    MLEs = c(NA, NA, opt$par[1], 0.25, muXiMLE, opt$par[-1])
  
  # Return results
  list(MLEs=MLEs, lambdaMLE=lambdaMLE, muZetaMLE=muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambda0MLE=NA, 
       muXiMLE=muXiMLE, logLikMLE=logLikMLE, splineParMLE=splinePar, hess=hess, optimTable=optimTable, 
       tvec=getTaperSpline(splinePar, nKnots=nKnots, normalize=normalizeTaper, dStar=dStar, latRange=latRange))
}

### gradient functions for coordinate descent
# gradient of E[Y] with respect to beta vec, the taper parameters.  Returns an n x p matrix 
# of derivatives.  If muZeta and sigmaZeta are given, then derivatives w.r.t. them will be 
# returned along with derivatives w.r.t. taper par in the order: (mu, sigma, taper par)
# NOTE: techically this gives the uplift gradient.  In order to really have the gradient of 
#       E[Y] you must take the negative of this
meanGrad = function(expectZeta, lambda, Xi, G, fault=csz, 
                    normalizeTaper=TRUE, dStar=21000, muZeta=NULL, sigmaZeta=NULL, 
                    index=1:nrow(G), normalModel=FALSE, phiZeta=NULL, anisotropic=FALSE) {
  
  # first get the indices we actually want (defaults to all indices) and make sure we still have a matrix
  G = matrix(G[index,], ncol=ncol(G))
  
  # get taper gradient
  depths = getFaultCenters(fault)[,3]
  tGrad = taperGrad(depths, Xi, lambda, dStar=dStar, normalize=normalizeTaper)
  
  # G %*% diag(E[zeta]) %*% dtvec/dbetavec
  taperPart = sweep(G, 2, expectZeta, "*") %*% tGrad
  
  # if user wants mu and sigma gradients as well:
  muPart = NULL
  sigmaPart = NULL
  phiPart = NULL
  if(!is.null(muZeta)) {
    # compute G %*% T
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, dStar=dStar, normalize=normalizeTaper)
    GT = sweep(G, 2, tvec, "*")
    
    # get d E[Y]/ d mu
    if(!normalModel)
      muPart = GT %*% expectZeta
    else
      muPart = rowSums(GT)
  }
  if(!is.null(sigmaZeta)) {
    # sigmaPart = sigmaZeta * muPart
    
    if(!normalModel) {
      # compute G %*% T if necessary
      if(is.null(GT)) {
        tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, dStar=dStar, normalize=normalizeTaper)
        GT = sweep(G, 2, tvec, "*")
      }
      
      # get d vector
      arealCSZCor = getArealCorMat(fault)
      DCD = arealCSZCor
      dLog = sqrt(diag(DCD))
      sigmas = sqrt(diag(DCD * sigmaZeta^2))
      
      # get d E[Y]/ d sigma
      sigmaPart = GT %*% (dLog * sigmas * expectZeta)
    }
    else
      sigmaPart = 0
  }
  
  if(!is.null(phiZeta)) {
    phiPart = 0
  }
  
  if(!anisotropic)
    cbind(muPart, sigmaPart, taperPart, phiPart)
  else
    cbind(muPart, sigmaPart, taperPart, phiPart, 0) # alpha component of gradient is 0
}

# gradient of Var(Y) with respect to beta vec, the taper parameters.  Returns an n x p matrix 
# of derivatives for p basis elements and n fault cells.  If muZeta and sigmaZeta are given, 
# then derivatives w.r.t. them will be returned along with derivatives w.r.t. taper par in 
# the order: (mu, sigma, taper par)
varGrad = function(varMatCSZ, lambda, tvec, Xi, G, fault=csz, normalizeTaper=TRUE, dStar=21000, 
                   muZeta=NULL, sigmaZeta=NULL, corZeta=NULL, index=1:nrow(G)) {
  # we only get the relevant index of the variance gradient (which is by default everything)
  G = matrix(G[index,], ncol=ncol(G))
  
  # get taper gradient
  depths = getFaultCenters(fault)[,3]
  tGrad = taperGrad(depths, Xi, lambda, dStar=dStar, normalize=normalizeTaper)
  
  ## we want (dt/dbeta_j) t^T + t (dt/dbeta_j)^T for each j
  # first make tensor of dimension (n, n, p) where (., ., j) is 
  # (dt/dbeta_j) t^T + t (dt/dbeta_j)^T
  symmGradTensor = array(dim=c(nrow(tGrad), nrow(tGrad), ncol(tGrad)))
  for(i in 1:dim(symmGradTensor)[3]) {
    outerMat = matrix(tGrad[,i] %o% tvec, nrow=length(tvec), ncol=length(tvec))
    symmGradTensor[,,i] = varMatCSZ * (outerMat + t(outerMat))
  }
  
  taperPart = matrix(nrow=nrow(G), ncol=ncol(tGrad))
  for(j in 1:dim(symmGradTensor)[3]) {
    taperPart[,j] = c(diag(G %*% symmGradTensor[,,j] %*% t(G)))
  }
  
  # now compute muPart and sigmaPart of gradient matrix if necessary
  muPart = NULL
  sigmaPart = NULL
  if(!is.null(muZeta)) {
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, dStar=dStar, normalize=normalizeTaper)
    GT = sweep(G, 2, tvec, "*")
    quadForm = diag(GT %*% corZeta %*% t(GT))
    
    muPart = quadForm * 2 * (exp(sigmaZeta^2) - 1) * exp(2*muZeta + sigmaZeta^2) # divide last sigmaZeta^2 by 2?
  }
  if(!is.null(sigmaZeta))
    sigmaPart = quadForm * 2 * sigmaZeta * exp(2*muZeta + sigmaZeta^2) * (2*exp(sigmaZeta^2) - 1)
  
  cbind(muPart, sigmaPart, taperPart)
}

# compute gradient of mean matching index, M = (1/n) sum ((mu_i - EY_i)/sigma(lat_i))^2
# with respect to beta.  Returns a p dimensional vector for p = #(beta)
MGrad = function(expectZeta, expectY, lambda, Xi, G, splineMeanVec, splineVarVec, fault=csz, 
                 normalizeTaper=TRUE, dStar=21000, pow=2, muZeta=NULL, sigmaZeta=NULL) {
  
  # get dEY/dbeta gradient matrix (n x p)
  EGradMat = meanGrad(expectZeta, lambda, Xi, G, fault, normalizeTaper, dStar, 
                      muZeta=muZeta, sigmaZeta=sigmaZeta)
  
  # compute gradient vector
  n = nrow(EGradMat)
  (-pow/n)*(matrix(abs((splineMeanVec - expectY)/splineVarVec)^(pow-1)*sign(splineMeanVec - expectY), nrow=1) %*% EGradMat)
}

# compute gradient of variance matching index, 
# V = (1/2n\bar{sigma}^4) sum (sigma(lat_i) - Var(Y_i) - Var(eps_i))^2
# with respect to beta.  Returns a p dimensional vector for p = #(beta)
VGrad = function(varMatCSZ, VarVecY, lambda, tvec, Xi, G, splineVarVec, fault=csz, 
                 subDat=dr1, normalizeTaper=TRUE, dStar=21000, pow=2, muZeta=NULL, 
                 sigmaZeta=NULL, corZeta=NULL) {
  
  # get dVarY/dbeta gradient matrix (n x p)
  VarGradMat = varGrad(varMatCSZ, lambda, tvec, Xi, G, fault=csz, normalizeTaper, dStar, 
                       muZeta=muZeta, sigmaZeta=sigmaZeta, corZeta=corZeta)
  
  # compute gradient vector
  n = nrow(VarGradMat)
  sigmaBar4 = mean(splineVarVec^2)
  if(pow == 2)
    return((-pow/(2*n*sigmaBar4)) * (matrix(abs(splineVarVec - VarVecY)^(pow-1)*sign(splineVarVec-VarVecY), nrow=1) %*% VarGradMat))
  else
    return((-pow/n) * (matrix(abs(splineVarVec - VarVecY)^(pow-1)*sign(splineVarVec-VarVecY), nrow=1) %*% VarGradMat))
}

# fit splines for mean and variance to subsidence data (not the negative of the subsidence data).  
# NOTE: this does not depend on taper and Zeta parameters
fitSplinesToSubDat = function(varSplineParInit=NULL, nKnotsMean=5, nKnotsVar=4, subDat=dr1) {
  # compute spline basis coefficients for mean of subsidence data
  latRange = c(40, 50)
  latSeq = seq(latRange[1], latRange[2], l=1000)
  XiMean = getSplineBasis(csz, nKnots=nKnotsVar, latRange=latRange, lats=subDat$Lat)
  XiMeanLat = getSplineBasis(csz, nKnots=nKnotsVar, latRange=latRange, lats=latSeq)
  mod = lm(subDat$subsidence ~ XiMean - 1)
  meanSplinePar = coef(mod)
  subMean = XiMean %*% meanSplinePar
  
  # center data so it's zero mean
  subCntr = subDat$subsidence - subMean
  
  # compute spline basis for variance of subsidence data
  XiVarLat = getSplineBasis(csz, nKnots=nKnotsVar, latRange=latRange, lats=latSeq)
  XiVar = getSplineBasis(csz, nKnots=nKnotsVar, latRange=latRange, lats=subDat$Lat)
  
  # compute log likelihood assuming independence
  datLogLik = function(splinePar) {
    # first make sure there's no negative variances
    varsLat = XiVarLat %*% splinePar
    if(any(varsLat < .05^2))
      return(-10000)
    
    # now compute log likelihood
    sds = sqrt(XiVar %*% splinePar + subDat$Uncertainty^2)
    sum(dnorm(subCntr, sd=sds, log=TRUE))
  }
  
  # set initial splinePar for variance spline if necessary
  if(is.null(varSplineParInit)) {
    # set initial variance to be constant (and ensure it's positive)
    meanVar = max(var(subDat$subsidence) - mean(subDat$Uncertainty^2), .01)
    varSplineParInit = c(meanVar, rep(0, nKnotsVar-1))
  }
  
  # choose optimal splinePar
  controls = list(fnscale = -1, reltol=1e-6)
  opt = optim(varSplineParInit, datLogLik, control=controls)
  varSplinePar = opt$par
  
  return(list(meanSplinePar=meanSplinePar, varSplinePar=varSplinePar))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions for computing the gradient of the likelihood (uses some 
# moment matching gradient functions)

# gradient of Cov(Y_i, Y_j) with respect to beta vec, the taper parameters, mu, and sigma.  
# Returns an n x p matrix of derivatives for p basis elements and n fault cells.  If 
# muZeta and sigmaZeta are given, then derivatives w.r.t. them will be returned along 
# with derivatives w.r.t. taper par in the order: (mu, sigma, taper par).  The 
# returned tensor is 3-dimensional with dimensions: (Y_i, Y_j, theta_k)
# NOTE: varMatCSZ is SigmaZeta (variance of zeta not of log-zeta)
covGrad = function(varMatCSZ, lambda, tvec=NULL, Xi, G, fault=csz, normalizeTaper=TRUE, dStar=21000, 
                   muZeta=NULL, sigmaZeta=NULL, corZeta=NULL, subDat=dr1, normalModel=FALSE, phiZeta=NULL, 
                   corMatCSZ=NULL, distMat=NULL, anisotropic=FALSE, alpha=1, squareStrikeDist=NULL, 
                   squareDipDist=NULL, nKnotsVar=5, doVarSpline=FALSE) {
  # get taper gradient
  depths = getFaultCenters(fault)[,3]
  tGrad = taperGrad(depths, Xi, lambda, dStar=dStar, normalize=normalizeTaper)
  
  # compute taper vector if necessary
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], alpha=2, lambda=lambda, dStar=dStar, normalize=normalizeTaper)
  
  ## we want (dt/dbeta_j) t^T + t (dt/dbeta_j)^T for each j
  # first make tensor of dimension (n, n, p) where (., ., j) is 
  # (dt/dbeta_j) t^T + t (dt/dbeta_j)^T
  symmGradTensor = array(dim=c(nrow(tGrad), nrow(tGrad), ncol(tGrad)))
  for(i in 1:dim(symmGradTensor)[3]) {
    outerMat = matrix(tGrad[,i] %o% tvec, nrow=length(tvec), ncol=length(tvec))
    symmGradTensor[,,i] = varMatCSZ * (outerMat + t(outerMat))
  }
  
  taperPart = array(dim=c(nrow(G), nrow(G), ncol(tGrad)))
  for(j in 1:dim(symmGradTensor)[3]) {
    taperPart[,,j] = G %*% symmGradTensor[,,j] %*% t(G)
  }
  
  # get areal correlation matrix (technically not correlation, rather covariance for sigma^2=1)
  if(is.null(phiZeta)) {
    arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
  }
  else {
    arealCSZCor = corMatCSZ
  }
  
  # compute corZeta if necessary
  if(!is.null(muZeta) || !is.null(sigmaZeta) && !normalModel) {
    warning("lognormal model is deprecated, and code is not guaranteed to work correctly in this case")
    if(doVarSpline)
      stop("variance spline with lognormal model is not supported")
    sigmas = sqrt(diag(varMatCSZ))
    corZeta = sweep(sweep(varMatCSZ, 1, sigmas, "/"), 2, sigmas, "/")
    theorVar = (exp(sigmaZeta^2) - 1) * exp(2*muZeta + sigmaZeta^2)
    d = sigmas/sqrt(theorVar)
    
    # also compute d and dLog vectors for mu and sigma gradients
    DCD = arealCSZCor
    dLog = sqrt(diag(DCD))
  }
  
  # now compute muPart alphaPart and phiPart and sigmaPart of gradient matrix if necessary
  muPart = NULL
  sigmaPart = NULL
  phiPart = NULL
  alphaPart = NULL
  if(!is.null(muZeta)) {
    if(!normalModel) {
      tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, dStar=dStar, normalize=normalizeTaper)
      GTD = sweep(G, 2, tvec*d, "*")
      quadForm = GTD %*% corZeta %*% t(GTD)
      
      muPart = quadForm * 2 * (exp(sigmaZeta^2) - 1) * exp(2*muZeta + sigmaZeta^2) # divide last sigmaZeta^2 by 2?
    }
    else
      muPart = matrix(0, nrow=nrow(G), ncol=nrow(G))
  }
  
  GT = NULL
  if(!is.null(sigmaZeta)) {
    # compute derivative of SigmaZeta w.r.t. sigmaZeta 
    
    if(!normalModel) {
      # (SigmaZeta is variance of zeta, not of log zeta)
      # first get C from DCD
      C = sweep(sweep(DCD, 1, dLog, "/"), 2, dLog, "/")
      crossTerm = 2 * C * outer(dLog, dLog)
      Dp2 = outer(dLog^2, dLog^2, "+")
      fullMat = Dp2 + crossTerm
      dSigmaZeta = fullMat * exp(sigmaZeta^2/2 * fullMat) - Dp2 * exp(sigmaZeta^2/2 * Dp2)
      dSigmaZeta = dSigmaZeta * (exp(2*muZeta) * sigmaZeta)
      
      GT = sweep(G, 2, tvec, "*")
      sigmaPart = GT %*% dSigmaZeta %*% t(GT)
    }
    else {
      GT = sweep(G, 2, tvec, "*")
      sigmaPart = 2*sigmaZeta * GT %*% arealCSZCor %*% t(GT)
    }
  }
  
  # phi and if necessary alpha component of gradient
  if(!is.null(phiZeta)) {
    
    out = corrGrad(distMat, phiZeta, arealCSZCor, anisotropic, alpha, squareStrikeDist, squareDipDist)
    if(!anisotropic)
      Cgrad = out
    else {
      Cgrad = out$dphi
      alphaGrad = out$dalpha
    }
    
    if(!normalModel) {
      phiPart = (varMatCSZ + exp(2*muZeta + sigmaZeta^2)) * sigmaZeta^2 * arealCSZCor * Cgrad
    }
    else {
      phiPart = sigmaZeta^2 * Cgrad
    }
    if(is.null(GT))
      GT = sweep(G, 2, tvec, "*")
    phiPart = GT %*% phiPart %*% t(GT)
    
    if(anisotropic) {
      if(!normalModel)
        stop("non-normality is not supported under the anisotropic assumption")
      
      alphaPart = sigmaZeta^2 * alphaGrad
      alphaPart = GT %*% alphaPart %*% t(GT)
    }
  }
  
  # return a 3-dimensional array, where the derivative w.r.t. each parameter is 
  # a slice with the third dimension fixed.
  if(!anisotropic)
    fullGrad = abind(muPart, sigmaPart, taperPart, phiPart, along=3)
  else
    fullGrad = abind(muPart, sigmaPart, taperPart, phiPart, alphaPart, along=3)
  
  # set different events to be independent (ijth element of gradient should be 0 if 
  # ith event not the same as jth event, since ijth element of covariance is 0)
  eventMask = eventsEqMask(subDat)
  fullGrad = sweep(fullGrad, c(1,2), eventMask, "*")
  
  return(fullGrad)
}

# compute the normal log-likelihood gradient for the subsidence data
# varMatCSZ is covariance of Zeta not of log(Zeta)
subLnLikGrad = function(expectZeta, expectY, lambda, Xi, G, varMatCSZ, tvec=NULL, 
                        fault=csz, subDat=dr1, normalizeTaper=TRUE, dStar=21000, pow=2, 
                        muZeta=NULL, sigmaZeta, SigmaYGrad=NULL, normalModel=FALSE, 
                        phiZeta=NULL, distMat=NULL, corMatCSZ=NULL, 
                        anisotropic=FALSE, alpha=1, squareStrikeDist=NULL, squareDipDist=NULL) {
  ## do some precomputations:
  
  # get taper vector
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambda, alpha=pow, dStar=dStar, normalize=normalizeTaper)
  
  # get covariance and mean gradients (covariance and mean of Y)
  if(is.null(SigmaYGrad))
    SigmaYGrad = covGrad(varMatCSZ, lambda, tvec, Xi, G, fault, normalizeTaper, dStar, 
                         muZeta, sigmaZeta, normalModel=normalModel, subDat=subDat, 
                         phiZeta=phiZeta, corMatCSZ=corMatCSZ, distMat=distMat, 
                         anisotropic=anisotropic, alpha=alpha, 
                         squareStrikeDist=squareStrikeDist, 
                         squareDipDist=squareDipDist)
  ExpYGrad = meanGrad(expectZeta, lambda, Xi, G, fault=fault, normalizeTaper, dStar, muZeta, 
                      sigmaZeta, normalModel=normalModel, phiZeta=phiZeta, anisotropic=anisotropic)
  
  nPar = dim(SigmaYGrad)[3]
  
  # now get Sigma_Y = GT Sigma_zeta T G^T + observation noise
  # NOTE: here Sigma_zeta is not the same as Sigma_{log(zeta)}.  This is the correlation matrix of log(zeta)
  if(is.null(corMatCSZ)) {
    arealCSZCor = getArealCorMat(fault, normalModel=normalModel)
    corMatCSZ = arealCSZCor
  }
  covMatCSZ = sigmaZeta^2 * corMatCSZ # get covariance of log(Zeta)
  out = estSubsidenceMeanCov(muZeta, lambda, covMatCSZ, G, tvec, fault=fault, 
                             subDat=subDat, normalModel=normalModel)
  SigmaY = out$Sigma
  diag(SigmaY) = diag(SigmaY) + subDat$Uncertainty^2
  
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

# # compute the gradient of the GPS data log likelihood for:
# # params[1] = muZeta
# # params[2] = sigmaZeta
# # params[3:7] = taperPar
# GPSLnLikGrad = function(muZeta, sigmaZeta, gpsDat=slipDatCSZ, corMatGPS=NULL, nPar=7, normalModel=FALSE) {
#   
#   # get gps data
#   logX = log(gpsDat$slip)
#   if(!normalModel)
#     x = log(gpsDat$slip)
#   else
#     x = gpsDat$slip
# 
#   if(is.null(corMatGPS)) {
#     # get spatial correlation parameters (computed based on fitGPSCovariance)
#     corPar = getCorPar(normalModel=normalModel)
#     phiZeta = corPar$phiZeta
#     nuZeta = corPar$nuZeta
#     
#     # compute correlation matrix for GPS data
#     coords = cbind(gpsDat$lon, gpsDat$lat)
#     corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
#                                onlyUpper=FALSE, smoothness=nuZeta,
#                                Distance="rdist.earth", Dist.args=list(miles=FALSE))
#   }
# 
#   ## compute expectation of logX
#   # Calculate the standard error vector of xi. Derivation from week 08_30_17.Rmd presentation.
#   # Transformation from additive error to  multiplicative lognormal model with asympototic
#   # median and variance matching.
#   sigmaXi = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
# 
#   # estimate muXi MLE with inverse variance weighting
#   ci = 1/sigmaXi^2
#   ci = ci/sum(ci)
#   
#   if(!normalModel) {
#     muXi = sum(logX*ci) - muZeta
#     muX = muXi + muZeta
#   }
#   else {
#     gammaEst = exp(sum(logX*ci))/muZeta
#     muX = gammaEst * muZeta
#   }
#   ## precompute some things before calculating the lnLik gradient
# 
#   # get Sigma_X gradient
#   if(!normalModel)
#     SigmaGrad = 2*sigmaZeta*corMatGPS
#   else {
#     dGammaEstDMu = -gammaEst/muZeta
#     SigmaGradMu = 2*gammaEst*sigmaZeta^2*corMatGPS * dGammaEstDMu
#     SigmaGradSigma = 2*gammaEst^2*sigmaZeta*corMatGPS
#   }
#   
#   # get Sigma_X (variance of log X)
#   if(!normalModel)
#     SigmaX = sigmaZeta^2 * corMatGPS
#   else {
#     SigmaX = gammaEst^2 * sigmaZeta^2 * corMatGPS
#     sigmaXi = gpsDat$slipErr
#   }
#   diag(SigmaX) = diag(SigmaX) + sigmaXi^2
# 
#   # compute Sigma_X^-1
#   SigmaXInv = solve(SigmaX)
# 
#   # compute gradient of Sigma_X^-1, which is Sigma_X^-1 dSigma_X^-1/dtheta_i Sigma_X^-1
#   if(!normalModel)
#     dSigmaXInv = SigmaXInv %*% SigmaGrad %*% SigmaXInv
#   else {
#     # compute gradient of Sigma_X^-1, which is Sigma_X^-1 dSigma_X^-1/dtheta_i Sigma_X^-1
#     dSigmaXInvMu = SigmaXInv %*% SigmaGradMu %*% SigmaXInv
#     dSigmaXInvSigma = SigmaXInv %*% SigmaGradSigma %*% SigmaXInv
#   }
#   
#   # turn muX into vector form
#   if(length(muX) == 1)
#     muX = rep(muX, nrow(gpsDat))
# 
#   ## Now compute likelihood gradient:
#   # four parts: X^T Sigma^-1 X, -2 X^T Sigma^-1 mu_X, mu_y^T Sigma^-1 mu_y, and log|Sigma|
# 
#   xQuad = rep(0, nPar)
#   xMuQuad = rep(0, nPar)
#   muQuad = rep(0, nPar)
#   logDet = rep(0, nPar)
#   # gradient of X^T Sigma^-1 X
#   if(!normalModel)
#     xQuad[2] = -t(x) %*% dSigmaXInv %*% x
#   else {
#     xQuad[1] = -t(x) %*% dSigmaXInvMu %*% x
#     xQuad[2] = -t(x) %*% dSigmaXInvSigma %*% x
#   }
#   
#   # gradient of -2 X^T Sigma^-1 mu_X
#   if(!normalModel)
#     xMuQuad[2] = 2*(t(x) %*% dSigmaXInv %*% muX)
#   else {
#     xMuQuad[1] = 2*(t(x) %*% dSigmaXInvMu %*% muX)
#     xMuQuad[2] = 2*(t(x) %*% dSigmaXInvSigma %*% muX)
#   }  
#   # gradient of mu_X^T Sigma_X^-1 mu_X
#   if(!normalModel)
#     muQuad[2] = - t(muX) %*% dSigmaXInv %*% muX
#   else {
#     muQuad[1] = - t(muX) %*% dSigmaXInvMu %*% muX
#     muQuad[2] = - t(muX) %*% dSigmaXInvSigma %*% muX
#   }  
#   # gradient of log|Sigma_X|
#   if(!normalModel)
#     logDet[2] = sum(diag(SigmaXInv %*% SigmaGrad))
#   else {
#     logDet[1] = sum(diag(SigmaXInv %*% SigmaGradMu))
#     logDet[2] = sum(diag(SigmaXInv %*% SigmaGradSigma))
#   }
#   # sum to get the full gradient (and multiply by -0.5)
#   fullGrad = -0.5 * (xQuad + xMuQuad + muQuad + logDet)
#   return(fullGrad)
#   # return(0)
# }

# compute the gradient of the GPS data log likelihood for:
# params[1] = muZeta
# params[2] = sigmaZeta
# params[3:7] = taperPar
GPSLnLikGrad = function(muZeta, sigmaZeta, gpsDat=slipDatCSZ, corMatGPS=NULL, nPar=7, normalModel=FALSE, 
                        phiZeta=NULL, distMatGPS=NULL, 
                        taperedGPSDat=FALSE, tvecGPS=NULL, nKnots=5, dStar=21000, latRange=c(40,50), fault=csz, 
                        normalizeTaper=TRUE, lambda=NULL, corGPS=FALSE, diffGPSTaper=FALSE, 
                        anisotropic=FALSE, alpha=1, squareStrikeDist=NULL, squareDipDist=NULL, 
                        doGammaSpline=FALSE, nKnotsGamma=7) {
  
  # get gps data
  logX = log(gpsDat$slip)
  if(!normalModel)
    x = log(gpsDat$slip)
  else
    x = gpsDat$slip
  
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
      ys = log(x)- log(muZeta) - log(tvecGPS)
      X = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
      logModel = lm(ys~X-1, weights=ci)
      muXi = X %*% logModel$coefficients
      gammaEst = exp(muXi)
    } else {
      muXi = sum((log(x)- log(muZeta) - log(tvecGPS))*ci)
      gammaEst = exp(muXi)
    }
  }
  else {
    muXi = sum((x- muZeta - log(tvecGPS))*ci)
  }
  
  if(is.null(corMatGPS)) {
    # get spatial correlation parameters (computed based on fitGPSCovariance)
    corPar = getCorPar(normalModel=normalModel)
    phiZeta = corPar$phiZeta
    nuZeta = corPar$nuZeta
    
    # compute correlation matrix for GPS data
    coords = cbind(gpsDat$lon, gpsDat$lat)
    corMatGPS = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                               onlyUpper=FALSE, smoothness=nuZeta,
                               Distance="rdist.earth", Dist.args=list(miles=FALSE))
  }
  
  ## precompute some things before calculating the lnLik gradient
  
  # get taper gradient if necessary
  # computing taper gradient is necessary under this model
  # compute spline basis matrix for subsidence data, get taper gradient
  if(taperedGPSDat) {
    Xi = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
    Xi2 = Xi
    if(diffGPSTaper)
      Xi2 = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGPS, latRange=latRange)
    depths = gpsDat$Depth
    tGrad = taperGrad(depths, Xi, lambda, dStar=dStar, normalize=normalizeTaper, diffGPSTaper=diffGPSTaper, 
                      nKnots==nKnotsGPS, Xi2)
  }
  
  # compute gradient of mu_X:
  ExpXGrad = meanXGrad(muZeta, nPar, tvecGPS, tGrad, gpsDat=gpsDat, normalModel=TRUE, 
                       taperedGPSDat=taperedGPSDat, anisotropic=anisotropic, 
                       doGammaSpline=doGammaSpline, nKnotsGamma=nKnotsGamma)
  
  # compute gradient of Sigma_X:
  SigmaGrad = covXGrad(muZeta, sigmaZeta, tvecGPS, tGrad, corMatGPS, nPar, distMatGPS, 
                       phiZeta, gpsDat, normalModel, taperedGPSDat, corGPS, 
                       anisotropic, alpha, squareStrikeDist, squareDipDist, 
                       doGammaSpline=doGammaSpline, nKnotsGamma=nKnotsGamma)
  
  # get Sigma_X (variance of log X)
  if(!normalModel) {
    SigmaZeta = sigmaZeta^2 * corMatGPS
    sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  }
  else {
    SigmaZeta = sigmaZeta^2 * corMatGPS
    sigmaXi = gpsDat$slipErr
  }
  if(taperedGPSDat && normalModel) {
    if(!doGammaSpline)
      SigmaX = gammaEst^2 * sweep(sweep(SigmaZeta, 1, tvecGPS, "*"), 2, tvecGPS, "*")
    else
      SigmaX = sweep(sweep(SigmaZeta, 1, tvecGPS * gammaEst, "*"), 2, tvecGPS * gammaEst, "*")
  }
  else
    SigmaX = SigmaZeta
  
  # add in measurement noise to SigmaX
  if(corGPS) {
    sigmaX = sqrt(diag(SigmaX))
    Cx = sweep(sweep(SigmaX, 1, 1/sigmaX, "*"), 2, 1/sigmaX, "*")
    SigmaXi = sweep(sweep(Cx, 1, sigmaXi, "*"), 2, sigmaXi, "*")
  }
  else
    SigmaXi = diag(sigmaXi)
  SigmaX = SigmaX + SigmaXi
  
  # compute Sigma_X^-1
  SigmaXInv = solve(SigmaX)
  
  # compute gradient of Sigma_X^-1, which is Sigma_X^-1 dSigma_X^-1/dtheta_i Sigma_X^-1
  dSigmaXInv = array(0, dim=c(length(x), length(x), nPar))
  for(i in 1:nPar) {
    dSigmaXInv[,,i] = SigmaXInv %*% SigmaGrad[,,i] %*% SigmaXInv
  }
  
  # compute expectation of X
  if(normalModel)
    expectX = gammaEst * muZeta * tvecGPS
  else
    expectX = muXi + muZeta + log(tvecGPS)
  
  ## Now compute likelihood gradient:
  # four parts: X^T Sigma^-1 X, -2 X^T Sigma^-1 mu_X, mu_y^T Sigma^-1 mu_y, and log|Sigma|
  xQuad = 1:nPar
  xMuQuad = 1:nPar
  muQuad = 1:nPar
  logDet = 1:nPar
  for(i in 1:nPar) {
    # gradient of X^T Sigma^-1 X
    xQuad[i] = -t(x) %*% dSigmaXInv[,,i] %*% x
    
    # gradient of -2 X^T Sigma^-1 mu_X
    xMuQuad[i] = -2*(t(x) %*% SigmaXInv %*% ExpXGrad[,i] - t(x) %*% dSigmaXInv[,,i] %*% expectX)
    
    # gradient of mu_X^T Sigma_X^-1 mu_X
    muQuad[i] = t(expectX) %*% SigmaXInv %*% ExpXGrad[,i] - t(expectX) %*% dSigmaXInv[,,i] %*% expectX +
      t(ExpXGrad[,i]) %*% SigmaXInv %*% expectX
    
    # gradient of log|Sigma_X|
    logDet[i] = sum(diag(SigmaXInv %*% SigmaGrad[,,i]))
  }
  
  # sum to get the full gradient (and multiply by -0.5)
  fullGrad = -0.5 * (xQuad + xMuQuad + muQuad + logDet)
  return(fullGrad)
  # return(0)
}

# compute the full gradient of the likelihood
# For now, I assume constant mean, so params is of the form:
# params[1]: muZeta
# params[2]: sigmaZeta
# params[3:(3+nKnots-1)]: taperPar
# NOTE: the last three inputs (phiZeta, useMVNApprox, and muZeta) are ignored, but 
#       are necessary since they are also inputs to fixedDataLogLikGrad.
fixedDataLogLikGrad = function(params, cszDepths, gpsDat=slipDatCSZ, 
                               subDat=dr1, fault=csz, verbose=TRUE, 
                               G=NULL, nKnots=5, normalizeTaper=TRUE, dStar=21000, 
                               constrLambda=TRUE, latRange=c(40,50), 
                               taperedGPSDat=FALSE, distMatGPS=NULL, distMatCSZ=NULL, 
                               corGPS=FALSE, diffGPSTaper=FALSE, nKnotsGPS=nKnots, anisotropic=FALSE, 
                               squareDipDistGps=NULL, squareStrikeDistGps=NULL, 
                               squareDipDistCsz=NULL, squareStrikeDistCsz=NULL, returnRow=TRUE, 
                               nKnotsGamma=7, doGammaSpline=FALSE, dStarGPS=dStar, 
                               southernFun=nKnotsGPS == 1, nKnotsVar=5, doVarSpline=FALSE) {
  
  ##### get parameters
  # muZeta = params[1]
  # muVecCSZ = rep(muZeta, nrow(fault))
  # muVecGPS = rep(muZeta, nrow(gpsDat))
  # sigmaZeta = params[2]
  # taperPar = params[3:(2+nKnots + diffGPSTaper*nKnotsGPS)]
  # if(diffGPSTaper) {
  #   taperParGPS = taperPar[-(1:nKnots)]
  #   taperPar = taperPar[1:nKnots]
  # }
  # if(!anisotropic)
  #   phiZeta = params[length(params)]
  # else
  #   phiZeta = params[length(params) - 1]
  # nuZeta=3/2
  # lambda0 = 0.25
  # if(anisotropic)
  #   alpha = params[length(params)]
  # else
  #   alpha = 1
  
  ##### get parameters
  out = getInputPar(params, fault, gpsDat, nKnots, diffGPSTaper, nKnotsGPS, taperedGPSDat=TRUE, anisotropic, normalModel=TRUE, 
                    nKnotsVar, doVarSpline)
  muZeta = out$muZeta
  muVecCSZ = out$muVecCSZ
  muVecGPS = out$muVecGPS
  sigmaZeta = out$sigmaZeta
  splinePar = out$taperPar
  splineParGPS = out$taperParGPS
  phiZeta = out$phiZeta
  nuZeta = out$nuZeta
  lambda0 = out$lambda0
  alpha = out$alpha
  varPar = out$varPar
  
  muZetaGPS = muZeta
  muZetaCSZ = muZeta
  vecMu=FALSE
  nuZeta = 3/2
  
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
  tvecGPS = taper(gpsDat$Depth, lambda=lambdaGPS, alpha=2, normalize=normalizeTaper, dStar=dStarGPS)
  
  # compute covariance matrices for latent process on fault (CSZ) and locking data (GPS)
  if(!anisotropic) {
    coords = cbind(fault$longitude, fault$latitude)
    corMatCSZ = stationary.cov(coords, Covariance="Matern", theta=phiZeta,
                               onlyUpper=FALSE, distMat=distMatCSZ, smoothness=nuZeta)
    corMatCSZL = t(chol(corMatCSZ))
    
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
    corMatCSZL = t(chol(corMatCSZ))
    
    corMatGPS = stationary.cov(NA, Covariance="Matern", theta=phiZeta,
                               onlyUpper=FALSE, distMat=distMatGPS, smoothness=nuZeta)
  }
  
  if(!doVarSpline) {
    XiVarCSZ = getSplineBasis(data.frame(list(latitude=subDat$lat)), nKnots=nKnotsVar, latRange=latRange)
    sigmaZetaCSZ = XiVarCSZ %*% varPar
    
    covMatCSZ = sweep(sweep(corMatCSZ, 1, sigmaZetaCSZ, "*"), 2, sigmaZetaCSZ, "*")
  } else {
    covMatCSZ = sigmaZeta^2 * corMatCSZ
  }
  
  # get zeta mean vector
  expectZeta = muVecCSZ
  
  # get zeta covariance matrix (covariance of zeta, NOT log(zeta))
  covZeta = covMatCSZ
  
  # make sure tvec is in reasonable range
  # if(max(tvec < .075))
  #   return(rep(0, length(params)))
  
  # compute expected uplift
  expectY = G %*% (tvec * expectZeta)
  
  # get covariance gradient
  SigmaYGrad = covGrad(covZeta, lambda, tvec, Xi, G, fault, normalizeTaper, dStar, 
                       muZeta, sigmaZeta, subDat=subDat, normalModel=TRUE, 
                       distMat=distMatCSZ, phiZeta=phiZeta, corMatCSZ=corMatCSZ, 
                       anisotropic=anisotropic, alpha=alpha, 
                       squareStrikeDist=squareStrikeDistCsz, 
                       squareDipDist=squareDipDistCsz, nKnotsVar=nKnotsVar, doVarSpline=doVarSpline)
  
  subGrad = subLnLikGrad(expectZeta, expectY, lambda, Xi, G, covZeta, tvec=tvec, 
                         fault, subDat, normalizeTaper, dStar, pow=2, 
                         muZeta, sigmaZeta, SigmaYGrad=SigmaYGrad, normalModel=normalModel, 
                         phiZeta=phiZeta, distMat=distMatCSZ, corMatCSZ=corMatCSZ, 
                         anisotropic=anisotropic, alpha=alpha, 
                         squareStrikeDist=squareStrikeDistCsz, 
                         squareDipDist=squareDipDistCsz)
  GPSGrad = GPSLnLikGrad(muZeta, sigmaZeta, gpsDat, nKnots=nKnots, normalModel=TRUE, 
                         corMatGPS=corMatGPS, nPar=length(params), phiZeta=phiZeta,
                         distMatGPS=distMatGPS, taperedGPSDat=TRUE, tvecGPS=tvecGPS, 
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

# separate input parameters into their respective categories, and also generate the respective scales for optimization (parscale for optim)
getInputPar = function(params, fault=csz, gpsDat=slipDatCSZ, nKnots=5, diffGPSTaper=FALSE, nKnotsGPS=5, 
                       taperedGPSDat=TRUE, anisotropic=FALSE, normalModel=TRUE, nKnotsVar=nKnots, doVarSpline=FALSE, 
                       finalFit=FALSE, includeGammaSpline=FALSE, nKnotsGamma=5, includeInflation=FALSE) {
  parscale = rep(0, length(params))
  parNames = as.character(parscale)
  
  ##### get parameters
  muZeta = params[1]
  parscale[1] = 2
  parNames[1] = "muZeta"
  muVecCSZ = rep(muZeta, nrow(fault))
  muVecGPS = rep(muZeta, nrow(gpsDat))
  
  startI = 2
  
  if(!doVarSpline) {
    sigmaZeta = params[2]
    varPar = NULL
    parscale[startI] = 2
    parNames[startI] = "sigmaZeta"
    taperPar = params[(startI + 1):(startI + nKnots + diffGPSTaper*nKnotsGPS)]
    parscale[(startI + 1):(startI+nKnots + diffGPSTaper*nKnotsGPS)] = rep(1, nKnots + diffGPSTaper*nKnotsGPS)
    parNames[(startI + 1):(startI+nKnots + diffGPSTaper*nKnotsGPS)] = paste0("beta^t_", 1:nKnots)
    
    startI = startI+nKnots + diffGPSTaper*nKnotsGPS + 1
  }
  else {
    sigmaZeta = NULL
    varPar = params[startI:(startI - 1 + nKnotsVar)]
    parscale[startI:(startI - 1 + nKnotsVar)] = 1
    parNames[startI:(startI - 1 + nKnotsVar)] = paste0("beta^s_", 1:nKnotsVar)
    startI = nKnotsVar + startI
    endI = startI - 1 + nKnots + diffGPSTaper*nKnotsGPS
    
    taperPar = params[startI:endI]
    parscale[startI:endI] = 1
    parNames[startI:(startI - 1 + nKnots)] = paste0("beta^t_", 1:nKnots)
    if(diffGPSTaper)
      parNames[(startI + nKnots):endI] = paste0("beta^t'_", 1:nKnotsGPS)
    
    startI = endI + 1
  }
  
  if(includeGammaSpline) {
    gammaPar = params[startI:(startI - 1 + nKnotsGamma)]
    parscale[startI:(startI - 1 + nKnotsGamma)] = 1
    parNames[startI:(startI - 1 + nKnotsGamma)] = paste0("beta^g_", 1:nKnotsGamma)
    startI = nKnotsGamma + startI
  } else {
    gammaPar = NULL
  }
  
  if(includeInflation) {
    lowInflation = params[startI]
    highInflation = params[startI + 1]
    parNames[startI] = "phi_l"
    parNames[startI + 1] = "phi_h"
    startI = startI + 2
  } else {
    lowInflation = NULL
    highInflation = NULL
  }
  
  if(diffGPSTaper) {
    taperParGPS = taperPar[-(1:nKnots)]
    taperPar = taperPar[1:nKnots]
  } else {
    taperParGPS = NULL
  }
  if(taperedGPSDat) {
    if(!anisotropic) {
      parscale[length(params)] = 25
      parNames[length(params)] = "phiZeta"
      phiZeta = params[length(params)]
      
      alpha = 1
    }
    else {
      parscale[length(params) - 1] = 25
      parNames[length(params) - 1] = "phiZeta"
      phiZeta = params[length(params) - 1]
      
      parscale[length(params)] = 2
      parNames[length(params)] = "alpha"
      alpha = params[length(params)]
    }
  }
  else{
    out = getCorPar(normalModel)
    phiZeta = out$phiZeta
    
    stop("untapered gps data no longer supported")
  }
  
  nuZeta=3/2
  lambda0 = 0.25
  
  if(finalFit)
    parscale = parscale / 25
  
  list(muZeta=muZeta, muVecCSZ=muVecCSZ, muVecGPS=muVecGPS, sigmaZeta=sigmaZeta, taperPar=taperPar, 
       taperParGPS=taperParGPS, phiZeta=phiZeta, nuZeta=nuZeta, lambda0=lambda0, alpha=alpha, varPar=varPar, 
       gammaPar=gammaPar, parscale=parscale, parNames=parNames, lowInflation=lowInflation, highInflation=highInflation)
}


