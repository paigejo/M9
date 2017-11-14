# this script contains functions related to priors, prior gradients, and faux observations


# return the set of points along the coast we want to ensure probably have <2m subsidence
# (aside from the true observations)
getFauxObs = function() {
  fauxObs = matrix(c(
    -126.998236, 49.818024, 
    -126.658771, 49.580435, 
    -126.395718, 49.411433, 
    -126.144274, 49.273292, 
    -125.541669, 48.920783, 
    -125.116728, 48.731546, 
    -124.698432, 48.590345, 
    -124.726965, 48.386018, 
    -124.717353, 48.143782, 
    -124.637417, 47.909926, 
    -124.444246, 47.754555, 
    -124.362284, 47.564145, 
    -124.413923, 42.510306, 
    -124.363353, 42.185168, 
    -124.211893, 41.873892, 
    -124.080846, 41.547055, 
    -124.098906, 41.242677, 
    -124.363634, 40.275861, 
    -124.102275, 40.093722
  ), byrow = TRUE, ncol=2)
  
  return(fauxObs)
}
getFauxDat = function(subDat=dr1) {
  fauxObs = getFauxObs()
  fauxDat = data.frame(list(Lon=fauxObs[,1], Lat=fauxObs[,2], Uncertainty=mean(subDat$Uncertainty), 
                            event="T1", subsidence=NA))
  return(fauxDat)
}

# get the Okada matrix for the faux observations
getFauxG = function(fauxObs = getFauxObs()) {
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  fauxG = okadaAll(csz, lonGrid, latGrid, fauxObs, slip=1, poisson=0.25)
  
  return(fauxG)
}

##### Shifted Frechet priors and gradients

# get prior on subsidence 95th quantile (exponential with 95% < 2m)
priorSubLikSFrechet = function(muZeta, sigmaZeta, tvec, G=NULL, fauxG=NULL, subDat=dr1, fauxObs=getFauxObs(), frechetPar=getSFrechetParSub()) {
  
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
  
  # compute the log likelihood
  dfrechet(maxSub, location=frechetPar[1], scale=frechetPar[2], shape=frechetPar[3], log=TRUE)
}

# get prior subsidence likelihood gradient
priorSubLikGradSFrechet = function(muZeta, sigmaZeta, lambda, Xi, tvec, SigmaYGrad, G=NULL, fauxG=NULL, subDat=dr1, 
                                   fauxObs=getFauxObs(), frechetPar=getSFrechetParSub(), fault=csz, normalizeTaper=TRUE, 
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
  maxSub = quants[maxSubI]
  
  ##### compute gradient
  # gradient of EY.  MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
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
  
  # compute the quantile gradient wrt thetaVec
  qGrad = matrix((expectGrad + sdGrad[maxSubI,] * qnorm(.95)), nrow=1)
  
  # compute the gradient we actually care about: the gradient of the log prior density
  mu = frechetPar[1]
  sigma = frechetPar[2]
  alpha = frechetPar[3]
  -(1+alpha)/(maxSub-mu) * qGrad + alpha*sigma^alpha * (maxSub-mu)^(-alpha - 1) * qGrad
}

# prior on expected earthquake slip.  95% chance that max slip < 40m
priorSlipLikSFrechet = function(muZeta, sigmaZeta, tvec, frechetPar=getSFrechetParSlip()) {
  
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
  
  dfrechet(maxSlip, location=frechetPar[1], scale=frechetPar[2], shape=frechetPar[3], log=TRUE)
}

# prior on expected earthquake slip.  95% chance that max slip < 40m
priorSlipLikGradSFrechet = function(muZeta, sigmaZeta, tvec, lambda, Xi, fault=csz, normalizeTaper=TRUE, 
                                    dStar=21000, frechetPar=getSFrechetParSlip()) {
  
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
  
  # compute the gradient we actually care about: the gradient of the log prior density
  mu = frechetPar[1]
  sigma = frechetPar[2]
  alpha = frechetPar[3]
  -(1+alpha)/(q95-mu) * q95Grad + alpha*sigma^alpha * (q95-mu)^(-alpha - 1) * q95Grad
}

##### prior Jacobian functions

# get prior log jacobian 
priorJacobian = function(muZeta, sigmaZeta, lambda, Xi, tvec, G, fauxG, subDat=dr1, 
                         fauxObs=getFauxObs(), fault=csz, normalizeTaper=TRUE, dStar=21000) {
  fullG = rbind(G, fauxG)
  
  ##### generate faux data
  # make our faux observations look like actual dr1 data (the event doesn't really matter in our case 
  # since we just want marginal variance and expectation)
  fauxDat = data.frame(list(Lon=fauxObs[,1], Lat=fauxObs[,2], Uncertainty=mean(subDat$Uncertainty), 
                            event="T1"))
  
  # combine relevant part of real and faux observations
  tmp = data.frame(list(Lon=subDat$Lon, Lat=subDat$Lat, Uncertainty=subDat$Uncertainty, event=subDat$event))
  fullDat = rbind(tmp, fauxDat)
  
  # get log zeta variances
  load("arealCSZCor.RData")
  corMatCSZ = arealCSZCor
  ds = sqrt(diag(corMatCSZ))
  varCSZ = sigmaZeta^2 * diag(corMatCSZ)
  sigmaCSZ = sqrt(varCSZ)
  covMatCSZ = sigmaZeta^2 * corMatCSZ
  expectZeta = exp(muZeta + diag(covMatCSZ)/2)
  
  # get zeta covariance matrix (covariance of zeta, NOT log(zeta))
  covZeta = exp(covMatCSZ) - 1
  covZeta = sweep(covZeta, 1, expectZeta, "*")
  covZeta = sweep(covZeta, 2, expectZeta, "*")
  
  # compute SigmaYGrad and SigmaYGradFaux
  SigmaYGrad = covGrad(covZeta, lambda, tvec, Xi, G, fault, normalizeTaper, dStar, muZeta, 
                       sigmaZeta, subDat=subDat)
  SigmaYGradFaux = covGrad(covZeta, lambda, tvec, Xi, fauxG, fault, normalizeTaper, dStar, muZeta, 
                           sigmaZeta, subDat=fauxDat)
  
  # compute slip log scale means and standard deviations
  muSlip = muZeta + log(tvec)
  sigmaSlip = sigmaCSZ
  
  # compute 95th percentiles of distributions
  quants = qlnorm(.95, muSlip, sigmaSlip)
  
  # get qS (maximum slip) and di
  maxSI = which.max(quants)
  qS = quants[maxSI]
  di = ds[maxSI]
  
  # get gradient of qY (max subsidence)
  qYG = qYGrad(muZeta, sigmaZeta, lambda, Xi, tvec, SigmaYGrad, SigmaYGradFaux, G, fauxG, subDat, 
               fauxObs, fault, normalizeTaper, dStar)
  qSG = qSGrad(qS, maxSI, tvec, Xi, lambda, dStar, normalizeTaper, fault)
  
  # return prior log absolute value jacobian
  jac = abs(qYG[,1] * qSG[2] - qYG[,2] * qSG[1])
  return(log(jac))
}

# inputs done
qYGrad = function(muZeta, sigmaZeta, lambda, Xi, tvec, SigmaYGrad, SigmaYGradFaux, G, fauxG, 
                  subDat=dr1, fauxObs=getFauxObs(), fault=csz, normalizeTaper=TRUE, dStar=21000) {
  
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
  # gradient of EYi. MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  expectGrad = -meanGrad(expectZeta, lambda, Xi, fullG, fault, normalizeTaper, dStar, 
                         muZeta, sigmaZeta, index=maxSubI)
  
  # get only parts of SigmaY gradient we need (gradient of cov(Y_i, Y_i) for each i, no cross-covariances needed)
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
  
  # now compute qY gradient
  qYGrad = expectGrad + qnorm(.95) * sdGrad[maxSubI,]
  return(qYGrad)
}

# inputs done
qSGrad = function(qS, qSI, tvec, Xi, lambda, dStar=21000, normalizeTaper=TRUE, fault=csz, computeBetaGrad=FALSE) {
  # get di
  load("arealCSZCor.RData")
  DCD = arealCSZCor
  ds = sqrt(diag(DCD))
  di = ds[qSI]
  
  # compute the gradient of qS w.r.t. muZeta and sigmaZeta
  qGrad = c(qS, di*qnorm(.95)*qS)
  
  # add beta gradient if necessary
  if(computeBetaGrad) {
    depths = fault$depth
    tGrad = taperGrad(depths, Xi, lambda, dStar, normalizeTaper)
    ti = tvec[qSI]
    tGradI = tGrad[qSI,]
    
    qGrad = c(qGrad, (tGradI * qS)/ti)
  }
  
  return(qGrad)
}

# inputs done
# dEYi/d(thetaVec, (mu, sigma))
# note that this truly returns the hessian in E[Y].  Hence, this gives the hessian in subsidence, not uplift.
EYHess = function(muZeta, sigmaZeta, tvec, fullG, EYGrad, tGrad, fullDat) {
  
  ##### Get EY
  # also compute d and dLog vectors for mu and sigma gradients
  load("arealCSZCor.RData")
  DCD = arealCSZCor
  dLog = sqrt(diag(DCD))
  covMatCSZ = sigmaZeta^2 * arealCSZCor
  expectZeta = exp(muZeta + diag(covMatCSZ)/2)
  
  ##### get 95th quantiles of subsidence and largest of them
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, covMatCSZ, fullG, subDat=fullDat, 
                                   tvec=tvec)
  mu = -mvnApprox$mu # MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  Sigma = mvnApprox$Sigma
  
  # get 95th quantiles of subsidence
  sigmas = sqrt(diag(Sigma) + fullDat$Uncertainty^2)
  quants = qnorm(.95, mu, sigmas)
  
  # get maximum expected subsidence
  maxSubI = which.max(quants)
  
  ##### compute hessian of EYi wrt mu and sigma
  # hessian wrt mu
  muEYHess = EYGrad[maxSubI,]
  
  # hessian wrt sigma
  Gi = fullG[maxSubI,]
  GTi = Gi * tvec
  EYHessBeta = ((dLog^2 * sigmaZeta) * expectZeta) %*% sweep(tGrad, 1, Gi, "*")
  EYHessMu = sum(GTi * (dLog^2 * sigmaZeta) * expectZeta)
  EYHessSigma = sum(GTi * (dLog^2 * expectZeta + dLog^4*sigmaZeta^2 * expectZeta))
  sigmaEYHess = c(EYHessMu, EYHessSigma, EYHessBeta)
  
  ##### return results
  EYHess = rbind(muEYHess, -sigmaEYHess)
  return(EYHess)
}

# inputs done
# dVarYi/d(thetaVec, (mu, sigma))
varYHess = function(muZeta, sigmaZeta, SigmaYGrad, SigmaYGradFaux, fullG, fullDat, tGrad, tvec) {
  
  ##### compute covZeta gradient
  # also compute d and dLog vectors for mu and sigma gradients
  load("arealCSZCor.RData")
  DCD = arealCSZCor
  dLog = sqrt(diag(DCD))
  covMatCSZ = sigmaZeta^2 * arealCSZCor
  expectZeta = exp(muZeta + diag(covMatCSZ)/2)
  
  # compute derivative of SigmaZeta w.r.t. sigmaZeta 
  # (SigmaZeta is variance of zeta, not of log zeta)
  # first get C from DCD
  C = sweep(sweep(DCD, 1, dLog, "/"), 2, dLog, "/")
  crossTerm = 2 * C * outer(dLog, dLog)
  Dp2 = outer(dLog^2, dLog^2, "+")
  fullMat = Dp2 + crossTerm
  dSigmaZeta = fullMat * exp(sigmaZeta^2/2 * fullMat) - Dp2 * exp(sigmaZeta^2/2 * Dp2)
  dSigmaZeta = dSigmaZeta * (exp(2*muZeta) * sigmaZeta)
  
  ##### get maximum subsidence index
  # get zeta covariance matrix (covariance of zeta, NOT log(zeta))
  covZeta = exp(covMatCSZ) - 1
  covZeta = sweep(covZeta, 1, expectZeta, "*")
  covZeta = sweep(covZeta, 2, expectZeta, "*")
  
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, covMatCSZ, fullG, subDat=fullDat, 
                                   tvec=tvec)
  mu = -mvnApprox$mu # MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  Sigma = mvnApprox$Sigma
  
  # get 95th quantiles of subsidence
  sigmas = sqrt(diag(Sigma) + fullDat$Uncertainty^2)
  quants = qnorm(.95, mu, sigmas)
  
  # get maximum expected subsidence
  maxSubI = which.max(quants)
  
  ########## compute dVarYi/d(thetaVec, sigma)
  ##### compute hess of VarYi wrt betaj *****************
  Gi = c(fullG[maxSubI,])
  GTi = c(Gi * tvec)
  hessVarYiBeta = 2 * GTi %*% dSigmaZeta %*% sweep(tGrad, 1, Gi, "*")
  
  ##### compute hess of VarYi wrt muZeta *****************
  hessVarYiMu = 2 * GTi %*% dSigmaZeta %*% GTi
  
  ##### compute second derivative of covZeta wrt sigmaZeta
  dSigmaZetaMod = fullMat^2 * exp(sigmaZeta^2/2 * fullMat) - Dp2^2 * exp(sigmaZeta^2/2 * Dp2)
  dSigmaZetaMod = dSigmaZetaMod * (exp(2*muZeta) * sigmaZeta^2)
  d2SigmaZeta = dSigmaZeta * (1/sigmaZeta) + dSigmaZetaMod
  
  ##### compute hess of VarYi wrt sigmaZeta *****************
  hessVarYiSigma = GTi %*% d2SigmaZeta %*% GTi
  
  ##### concatenate results
  sigmaHessVarYi = matrix(c(hessVarYiMu, hessVarYiSigma, hessVarYiBeta), nrow=1)
  
  ########## compute dVarYi/d(thetaVec, mu)
  ##### get only parts of SigmaY gradient we need (gradient of cov(Y_i, Y_i), no cross-covariances needed)
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
  
  ##### now get the hessian
  muHessVarYi = 2 * varGrad[maxSubI,]
  
  ########## concatenate everything and return results
  hessVarYi = rbind(muHessVarYi, sigmaHessVarYi)
  return(hessVarYi)
}

# inputs done
# compute hessian of qY: d^2 qY/d(thetaVec, mu/sigma)
qYHess = function(muZeta, sigmaZeta, lambda, Xi, tvec, EYGrad, SigmaYGrad, SigmaYGradFaux, 
                  G=NULL, fauxG=NULL, subDat=dr1, fauxObs=getFauxObs(), fault=csz, 
                  normalizeTaper=TRUE, dStar=21000) {
  
  fullG = rbind(G, fauxG)
  
  ##### generate faux data
  # make our faux observations look like actual dr1 data (the event doesn't really matter in our case 
  # since we just want marginal variance and expectation)
  fauxDat = data.frame(list(Lon=fauxObs[,1], Lat=fauxObs[,2], Uncertainty=mean(subDat$Uncertainty), 
                            event="T1"))
  
  # combine relevant part of real and faux observations
  tmp = data.frame(list(Lon=subDat$Lon, Lat=subDat$Lat, Uncertainty=subDat$Uncertainty, event=subDat$event))
  fullDat = rbind(tmp, fauxDat)
  
  ##### compute taper gradient
  depths = fault$depth
  tGrad = taperGrad(depths, Xi, lambda, dStar, normalizeTaper)
  
  # get hessians of EY and VarY
  EHess = EYHess(muZeta, sigmaZeta, tvec, fullG, EYGrad, tGrad, fullDat)
  VarHess = varYHess(muZeta, sigmaZeta, SigmaYGrad, SigmaYGradFaux, fullG, fullDat, tGrad, tvec)
  
  # also compute d and dLog vectors for mu and sigma gradients
  load("arealCSZCor.RData")
  covMatCSZ = sigmaZeta^2 * arealCSZCor
  
  ##### get maximum subsidence index
  
  # This is the key step: approximate G %*% T %*% Zeta with a MVN
  mvnApprox = estSubsidenceMeanCov(muZeta, lambda, sigmaZeta, covMatCSZ, fullG, subDat=fullDat, 
                                   tvec=tvec)
  mu = -mvnApprox$mu # MUST TAKE NEGATIVE HERE, SINCE SUBSIDENCE NOT UPLIFT!!!!
  Sigma = mvnApprox$Sigma
  
  # get 95th quantiles of subsidence
  sigmas = sqrt(diag(Sigma) + fullDat$Uncertainty^2)
  quants = qnorm(.95, mu, sigmas)
  
  # get maximum expected subsidence index and variance
  maxSubI = which.max(quants)
  varYi = Sigma[maxSubI, maxSubI]
  if( maxSubI > nrow(G)) {
    maxSubI = maxSubI - nrow(G) + 1
    SigmaYiGrad = SigmaYGradFaux[maxSubI, maxSubI,]
  }
  else
    SigmaYiGrad = SigmaYGrad[maxSubI, maxSubI,]
  
  ##### compute qYHess (1st derivative is row, 2nd is col)
  term1 = -0.5 * (varYi)^(-3/2) * outer(SigmaYiGrad[1:2], SigmaYiGrad)
  term2 = varYi^(-1/2) * VarHess
  hessQY = EHess + (qnorm(.95)/2) * (term1 + term2)
  return(hessQY)
}

# inputs done
# compute the gradient of the log jacobian factor for the prior
priorJacobianGrad = function(muZeta, sigmaZeta, lambda, Xi, tvec, EYGrad, SigmaYGrad, SigmaYGradFaux, 
                             G=NULL, fauxG=NULL, subDat=dr1, fauxObs=getFauxObs(), fault=csz, 
                             normalizeTaper=TRUE, dStar=21000) {
  
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
  ds = sqrt(diag(corMatCSZ))
  
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
  
  # get log zeta variances
  varCSZ = sigmaZeta^2 * diag(corMatCSZ)
  sigmaCSZ = sqrt(varCSZ)
  
  # compute slip log scale means and standard deviations
  muSlip = muZeta + log(tvec)
  sigmaSlip = sigmaCSZ
  
  # compute 95th percentiles of distributions
  quants = qlnorm(.95, muSlip, sigmaSlip)
  
  # get max slip
  maxSlipI = which.max(quants)
  qS = quants[maxSlipI]
  di = ds[maxSlipI]
  
  # get maximum expected observed subsidence
  maxSubI = which.max(quants)
  qY = quants[maxSubI]
  
  ##### get relevant gradients/hessians of qS and qY
  hessQY = qYHess(muZeta, sigmaZeta, lambda, Xi, tvec, EYGrad, SigmaYGrad, SigmaYGradFaux, 
                  G, fauxG, subDat, fauxObs, fault, normalizeTaper, dStar)
  
  gradQY = qYGrad(muZeta, sigmaZeta, lambda, Xi, tvec, SigmaYGrad, SigmaYGradFaux, G, fauxG, 
                  subDat, fauxObs, fault, normalizeTaper, dStar)
  
  gradQS = qSGrad(qS, maxSlipI, tvec, Xi, lambda, dStar, normalizeTaper, fault, computeBetaGrad=TRUE)
  
  ##### compute jacobian and its gradient
  jac = gradQY[1] * gradQS[2] - gradQY[2] * gradQS[1]
  gradJac = hessQY[1,] * gradQS[2] + gradQY[1] * di * qnorm(.95) * gradQS - hessQY[2,] * qS - gradQY[2] * gradQS
  gradLogJac = sign(jac)*gradJac/abs(jac)
  
  return(gradLogJac)
}


##### functions for getting and plotting prior parameters
getFrechetPar = function(ps=c(.05, .95), qs=c(.75, 5)) {
  p1 = ps[1]
  p2 = ps[2]
  q1 = qs[1]
  q2 = qs[2]
  
  # get alpha (shape)
  numer = log(log(p1)/log(p2))
  denom = log(q2) - log(q1)
  alpha = numer/denom
  
  # get sigma (scale)
  sigma=q1 * (-log(p1))^(1/alpha)
  
  return(c(sigma=sigma, alpha=alpha))
}
plotFrechetPrior = function(ps=c(.05, .95), qs=c(.5, 3), ylim=NULL) {
  xs = seq(0, 5, l=200)
  frechetPar = getFrechetPar(ps, qs)
  ys = dfrechet(xs, scale=frechetPar[1], shape=frechetPar[2])
  if(is.null(ylim))
    ylim = c(0, max(ys))
  plot(xs, ys, type="l", main="Frechet prior", ylim=ylim, ylab="Density", xlab="x")
}
getSFrechetPar = function(ps=c(.05, .95), qs=c(.5, 5), m=-.5) {
  p1 = ps[1]
  p2 = ps[2]
  q1 = qs[1]
  q2 = qs[2]
  
  # get alpha (shape)
  numer = log(log(p1)/log(p2))
  denom = log(q2-m) - log(q1-m)
  alpha = numer/denom
  
  # get sigma (scale)
  sigma=(q1-m) * (-log(p1))^(1/alpha)
  
  return(c(m=m, sigma=sigma, alpha=alpha))
}
getSFrechetParSub = function() {
  getSFrechetPar(ps=c(.05, .95), qs=c(.5, 3), m=-1)
}
getSFrechetParSlip = function(ps=c(.05, .95), qs=c(5, 40)) {
  getSFrechetPar(ps=ps, qs=qs, m=0)
}
plotSFrechetPrior = function(ps=c(.05, .95), qs=c(.5, 3), m=-1, ylim=NULL) {
  xs = seq(0, 5, l=200)
  sfrechetPar = getSFrechetPar(ps, qs, m)
  ys = dfrechet(xs, location=m, scale=sfrechetPar[2], shape=sfrechetPar[3])
  if(is.null(ylim))
    ylim = c(0, max(ys))
  plot(xs, ys, type="l", main="Shifted Frechet prior", ylim=ylim, ylab="Density", xlab="x")
}
plotSFrechetPriorSub = function(ylim=NULL) {
  xs = seq(0, 5, l=200)
  sfrechetPar = getSFrechetParSub()
  ys = dfrechet(xs, location=m, scale=sfrechetPar[2], shape=sfrechetPar[3])
  if(is.null(ylim))
    ylim = c(0, max(ys))
  plot(xs, ys, type="l", main="Shifted Frechet prior", ylim=ylim, ylab="Density", xlab="x")
}
plotSFrechetPriorSlip = function(ps=c(.05, .95), qs=c(5, 40), ylim=NULL) {
  xs = seq(0, 60, l=200)
  sfrechetPar = getSFrechetParSlip(ps=ps, qs=qs)
  ys = dfrechet(xs, scale=sfrechetPar[2], shape=sfrechetPar[3])
  if(is.null(ylim))
    ylim = c(0, max(ys))
  plot(xs, ys, type="l", main="Shifted Frechet slip prior", ylim=ylim, ylab="Density", xlab="x")
}

getGumbelPar = function(ps=c(.05, .95), qs=c(.5, 3)) {
  p1 = ps[1]
  p2 = ps[2]
  q1 = qs[1]
  q2 = qs[2]
  
  # get kappa (not a parameter but used in my derivations)
  kappa = log(-log(p1)) / log(-log(p2))
  
  # get mu (location)
  mu = (q1 - kappa * q2)/(1 - kappa)
  
  # get beta (scale)
  beta = -(q1 - mu)/(log(-log(p1)))
  
  return(c(mu=mu, beta=beta))
}
getGumbelParSub = function(ps=c(.05, .95), qs=c(.5, 3)) {
  getGumbelPar(ps=ps, qs=qs)
}
getGumbelParSlip = function(ps=c(.05, .95), qs=c(5, 20)) {
  getGumbelPar(ps=ps, qs=qs)
}
plotGumbelPrior = function(ps=c(.05, .95), qs=c(.5, 3), ylim=NULL) {
  xs = seq(0, 5, l=200)
  gumbelPar = getGumbelPar(ps=ps, qs=qs)
  ys = dgumbel(xs, location=gumbelPar[1], scale=gumbelPar[2])
  if(is.null(ylim))
    ylim = c(0, max(ys))
  plot(xs, ys, type="l", main="Gumbel prior", ylim=ylim, ylab="Density", xlab="x")
}
plotGumbelPriorSlip = function(ps=c(.05, .95), qs=c(5, 20), ylim=NULL) {
  xs = seq(0, 60, l=200)
  gumbelPar = getGumbelParSlip(ps=ps, qs=qs)
  ys = dgumbel(xs, location=gumbelPar[1], scale=gumbelPar[2])
  if(is.null(ylim))
    ylim = c(0, max(ys))
  plot(xs, ys, type="l", main="Gumbel prior slip", ylim=ylim, ylab="Density", xlab="x")
}
