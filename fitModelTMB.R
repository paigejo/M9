library(TMB)
source("setup.R")

# fits the combined model with normalized taper using TMB jointly with the variance inflation 
# parameters as well as the gamma spline parameters
fitModelTMB = function(initParams, gpsDat=slipDatCSZ, 
                       G=NULL, subDat=dr1, fault=csz, nKnots=5, 
                       dStar=21000, maxit=500, latRange=c(40, 50), 
                       normalModel=TRUE, doHess=TRUE, corGPS=FALSE, finalFit=FALSE, 
                       diffGPSTaper=FALSE, nKnotsGPS=nKnots, reltol=1e-8, 
                       nKnotsGamma=7, nKnotsVar=5, 
                       dStarGPS=dStar, seed=123, debug=FALSE) {
  
  # get input parameters
  out = getInputPar(initParams, fault, gpsDat, nKnots, diffGPSTaper=FALSE, nKnotsGPS, taperedGPSDat=TRUE, 
                    anisotropic=TRUE, normalModel=TRUE, nKnotsVar, doVarSpline=TRUE, 
                    includeGammaSpline=TRUE, nKnotsGamma=nKnotsGamma, includeInflation=TRUE)
  mu = out$muZeta
  betaTaper = out$taperPar
  # taperParGPS = out$taperParGPS
  phi = out$phiZeta
  # lambda0 = out$lambda0
  alpha = out$alpha
  betasd = out$varPar
  betaGamma = out$gammaPar
  lowInflate = out$lowInflation
  highInflate = out$highInflation
  parscale = out$parscale
  parNames = out$parNames
  
  # PARAMETER(mu);
  # PARAMETER_VECTOR(betaTaper);
  # PARAMETER_VECTOR(betasd);
  # PARAMETER_VECTOR(betaGamma);
  # PARAMETER(phi);
  # PARAMETER(alpha);
  # PARAMETER(lowInflate);
  # PARAMETER(highInflate);
  
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
  straightenedGpsCoords = transformation(cbind(gpsDat$lon, gpsDat$lat))
  
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
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnots, latRange=latRange)
  DSStrikeCSZ = squareStrikeDistCsz
  DSDipCSZ = squareDipDistCsz
  DSStrikeGPS = squareStrikeDistGps
  DSDipGPS = squareDipDistGps
  zeroMask = eventsEqMask(subDat)
  lowI = as.numeric(as.numeric(subDat$quality) != 1)
  y = -subDat$subsidence
  ysd = subDat$Uncertainty
  x = gpsDat$slip
  xsd = sqrt(log(.5*(sqrt(4*gpsDat$slipErr^2/gpsDat$slip^2 + 1) + 1)))
  xDepths = gpsDat$Depth
  faultDepths = cszDepths
  dStar = dStar
  sdBasisX = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
  sdBasisY = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsVar, latRange=latRange)
  gammaBasis = getSplineBasis(data.frame(list(latitude=gpsDat$lat)), nKnots=nKnotsGamma, latRange=latRange)
  
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
               sdBasisY=sdBasisY, gammaBasis=gammaBasis, G=G)
  parameters = list(mu=mu, betaTaper=betaTaper, betasd=betasd, betaGamma=betaGamma, phi=phi, 
                    alpha=alpha, lowInflate=lowInflate, highInflate=highInflate)
  obj <- MakeADFun(data, parameters, DLL="fitModelTMB")
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

getInitialParameters = function() {
  # in order: 
  # mu, five SD parameters, five taper parameters, seven gamma parameters, 
  # low inflation, high inflation, spatial range, alpha/anisotropy parameter
  c(20, 15, rep(0, 4), 1, rep(0, 4), 1, rep(0, 6), 1, 1, 175, 1)
}

test = fitModelTMB(getInitialParameters(), debug=TRUE)







