library(TMB)
# source("setup.R")

# fits the combined model with normalized taper using TMB jointly with the variance inflation 
# parameters as well as the gamma spline parameters
fitModelTMB = function(initParams=NULL, gpsDat=slipDatCSZ, gpsDepthThreshold=21000, 
                       G=NULL, subDat=dr1, fault=csz, nKnots=5, 
                       dStar=25000, maxit=500, latRange=c(40, 50), 
                       finalFit=FALSE, 
                       diffGPSTaper=FALSE, nKnotsGPS=nKnots, reltol=1e-8, 
                       nKnotsGamma=7, nKnotsVar=5, includeGammaSpline=FALSE, 
                       dStarGPS=40000, seed=123, debug=FALSE, doVarSpline=TRUE, 
                       fixInflation=TRUE) {
  
  # threshold the gps data and multiplied standard deviations by three as in the literature
  threshSlipDat = gpsDat[gpsDat$Depth<gpsDepthThreshold,]
  threshSlipDat$slipErr = threshSlipDat$slipErr*3
  
  # get input parameters and separate them
  if(is.null(initParams))
    initParams = getInitialParameters(nKnots, nKnotsVar, nKnotsGamma, logScale=TRUE)
  
  out = getInputPar(initParams, fault, threshSlipDat, nKnots, diffGPSTaper=FALSE, nKnotsGPS, taperedGPSDat=TRUE, 
                    anisotropic=TRUE, normalModel=TRUE, nKnotsVar, doVarSpline=doVarSpline, 
                    includeGammaSpline=includeGammaSpline, nKnotsGamma=nKnotsGamma, includeInflation=TRUE)
  logmu = out$muZeta
  betaTaper = out$taperPar
  betaTaperIntercept = betaTaper[1]
  betaTaper = betaTaper[-1]
  # taperParGPS = out$taperParGPS
  logphi = out$phiZeta
  # lambda0 = out$lambda0
  logalpha = out$alpha
  betasd = out$varPar
  betasdIntercept = betasd[1]
  betasd = betasd[-1]
  betaGamma = out$gammaPar
  betaGammaIntercept = betaGamma[1]
  betaGamma = betaGamma[-1]
  loglowInflate = out$lowInflation
  loghighInflate = out$highInflation
  parscale = out$parscale
  parNames = out$parNames
  
  # PARAMETER(logmu);
  # PARAMETER_VECTOR(betaTaperIntercept);
  # PARAMETER_VECTOR(betaTaper);
  # PARAMETER_VECTOR(betasdIntercept);
  # PARAMETER_VECTOR(betasd);
  # PARAMETER_VECTOR(betaGammaIntercept);
  # PARAMETER_VECTOR(betaGamma);
  # PARAMETER(logphi);
  # PARAMETER(logalpha);
  # PARAMETER(loglowInflate);
  # PARAMETER(loghighInflate);
  
  ##### calculate correlation matrices for Zeta in km (for CSZ grid and the GPS data)
  
  # get Okada linear transformation matrix
  if(is.null(G)) {
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    G = okadaAll(fault, lonGrid, latGrid, cbind(subDat$Lon, subDat$Lat), slip=1, poisson=0.25)
  }
  
  ##### calculate depths of the centers of the CSZ subfaults
  cszDepths = getFaultCenters(fault)[,3]
  
  ### Rather than training the fault, we redefine an axis to be the strike access in Euclidean space
  ### using a Lambert projection and PCA
  out = straightenFaultLambert()
  faultGeomStraight = out$fault
  scale = out$scale
  parameters = out$projPar
  transformation = out$transformation
  
  ##### compute distance matrices for straightened fault and for straightened gps data
  cszStraight = divideFault2(faultGeomStraight)
  centers = getFaultCenters(csz)[,1:2]
  newCenters = transformation(centers)
  cszStraight$centerX = newCenters[,1]
  cszStraight$centerY = newCenters[,2]
  straightenedGpsCoords = transformation(cbind(threshSlipDat$lon, threshSlipDat$lat))
  
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
  
  # compute all required inputs for the TMB cpp code
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)[,-1]
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnots, latRange=latRange)[,-1]
  DSStrikeCSZ = squareStrikeDistCsz
  DSDipCSZ = squareDipDistCsz
  DSStrikeGPS = squareStrikeDistGps
  DSDipGPS = squareDipDistGps
  zeroMask = eventsEqMask(subDat)
  lowI = as.numeric(as.numeric(subDat$quality) != 1)
  y = -subDat$subsidence
  ysd = subDat$Uncertainty
  x = threshSlipDat$slip
  xsd = sqrt(log(.5*(sqrt(4*threshSlipDat$slipErr^2/threshSlipDat$slip^2 + 1) + 1)))
  xDepths = threshSlipDat$Depth
  faultDepths = cszDepths
  dStar = dStar
  dStarGPS = dStarGPS
  sdBasisX = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsVar, latRange=latRange)[,-1]
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)[,-1]
  gammaBasis = getSplineBasis(data.frame(list(latitude=threshSlipDat$lat)), nKnots=nKnotsGamma, latRange=latRange)[,-1]
  
  # combine the initial parameter guesses
  
  # compile the function and its gradient
  if(!debug)
    compile("fitModelTMB.cpp")
  else
    compile("fitModelTMB.cpp","-O0 -g")
  dyn.load(dynlib("fitModelTMB"))
  set.seed(seed)
  data <- list(y=y, x=x, lambdaBasisY=lambdaBasisY, lambdaBasisX=lambdaBasisX, DSStrikeCSZ=DSStrikeCSZ, 
               DSDipCSZ=DSDipCSZ, DSStrikeGPS=DSStrikeGPS, DSDipGPS=DSDipGPS, zeroMask=zeroMask, lowI=lowI, 
               ysd=ysd, xsd=xsd, xDepths=xDepths, faultDepths=faultDepths, dStar=dStar, sdBasisX=sdBasisX, 
               sdBasisY=sdBasisY, gammaBasis=gammaBasis, G=G, dStarGPS=dStarGPS)
  parameters = list(logmu=logmu, betaTaper=betaTaper, betasd=betasd, betaGamma=betaGamma, logphi=logphi, 
                    logalpha=logalpha, loglowInflate=loglowInflate, loghighInflate=loghighInflate, 
                    betaTaperIntercept = betaTaperIntercept, betasdIntercept=betasdIntercept, 
                    betaGammaIntercept=betaGammaIntercept)
  
  # before we got going, pick which parameters will be fixed
  map = list()
  if(!includeGammaSpline)
    map = c(map, list(betaGamma=NA))
  if(!doVarSpline)
    map = c(map, list(betasd=NA))
  if(!doVarSpline)
    map = c(map, list(betasd=NA))
  if(fixInflation)
    map = c(map, list(loghighInflate=NA, loglowInflate=NA))
  
  obj <- MakeADFun(data, parameters, DLL="fitModelTMB", map=map)
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  obj$he()    ## <-- Analytical hessian
  sdreport(obj)
  
  # Return results
  # return(list(MLEs=MLEs, muZetaMLE=muZetaMLE, sigmaZetaMLE=sigmaZetaMLE, lambdaMLE=NA, 
  #             gammaEst=gammaEst, logLikMLE=logLikMLE, splineParMLE=splinePar, phiMLE=phiMLE, 
  #             alphaMLE=alphaMLE, hess=hess, optimTable=optimTable, 
  #             tvec=tvec, tvecGPS=tvecGPS, optPar=opt$par, optGrad=optGrad, 
  #             strikeDistGps = strikeDistGps, dipDistGps = dipDistGps, 
  #             strikeDistCsz = strikeDistCsz, dipDistCsz = dipDistCsz))
  list(obj=obj, opt=opt)
}

getInitialParameters = function(nKnots=5, nKnotsVar=5, nKnotsGamma=7, logScale=FALSE) {
  # in order: 
  # mu, five SD parameters, five taper parameters, seven gamma parameters, 
  # low inflation, high inflation, spatial range, alpha/anisotropy parameter
  out = c(20, 15, rep(0, nKnotsVar - 1), 1, rep(0, nKnots - 1), 1, rep(0, nKnotsGamma - 1), 1.75, 1.25, 175, 1)
  
  if(logScale) {
    # mean and standard deviation intercept
    out[1:2] = log(out[1:2])
    
    # intercept for gamma
    out[(2 + nKnotsVar + nKnots)] = log(out[(2 + nKnotsVar + nKnots)])
    
    # inflation, scale, and anisotropy
    out[(length(out) - 3):length(out)] = log(out[(length(out) - 3):length(out)])
  }
  
  out
}

# test = fitModelTMB(getInitialParameters(), debug=TRUE)







