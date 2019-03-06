##### functions for spline bases

# return the spline basis matrix
getSplineBasis = function(fault=csz, latRange=c(40,50), nKnots=5, 
                          lats=fault$latitude, southernFun=FALSE, ...) {
  if(!southernFun)
    splineMat = do.call("getBSplineBasis", c(list(fault=fault, latRange=latRange, nKnots=nKnots, lats=lats), list(...)))
  else
    splineMat = do.call("getBSplineBasis", c(list(fault=fault, latRange=latRange, nKnots=nKnots-1, lats=lats), list(...)))
  # do.call("getNormalBasis", c(list(fault=fault, latRange=latRange, nKnots=nKnots, lats=lats, intercept=FALSE), list(...)))
  
  if(southernFun) {
    # make an extra basis function just for the south:
    splineMat2 = do.call("getBSplineBasis", c(list(fault=fault, latRange=latRange, nKnots=nKnots+4, lats=lats), list(...)))
    splineMat = cbind(splineMat, splineMat2[,ncol(splineMat2)])
  }
  
  splineMat
}

# gets taper vector for the given fault and spline basis parameters
getTaperSpline = function(splinePar, fault=csz, latRange=c(40, 50), nKnots=5, 
                          normalize=TRUE, dStar=21000, lats=fault$latitude, 
                          logScale=FALSE) {
  splineMat = getSplineBasis(fault, latRange, nKnots, lats)
  
  lambdas = splineMat %*% splinePar
  if(logScale)
    lambdas = exp(lambdas)
  c(taper(getFaultCenters(fault)[,3], lambda=lambdas, dStar=dStar, normalize=normalize))
}

# get b-spline basis matrix
getBSplineBasis = function(fault=csz, latRange=range(fault$latitude), nKnots=5, 
                           lats=fault$latitude, intercept=FALSE) {
  # if only one knot, return the constant 1 basis
  if(nKnots == 1) {
    return(matrix(rep(1, length(lats)), ncol=1))
  }
  
  # splineMat = bs(lats, df=nKnots-1, intercept=intercept, Boundary.knots=latRange)
  splineMat = bs(latRange[1] + latRange[2] - lats, df=nKnots-1, intercept=intercept, Boundary.knots=latRange)
  if(!intercept)
    return(cbind(rep(1, nrow(splineMat)), splineMat))
  else
    return(splineMat)
}

# make gaussian spline basis matrix
getNormalBasis = function(fault=csz, latRange=range(fault$latitude), nKnots=5, 
                          lats=fault$latitude, sds=NULL, intercept=FALSE) {
  
  if(!intercept) {
    nKnots = nKnots+1
  }
  
  if(is.null(sds))
    sds = (latRange[2] - latRange[1])/(2*nKnots-2) # so each knot is 1 sd away
  nSDAway = 2/3
  sds = sds/nSDAway
  means = seq(latRange[1], latRange[2], l=nKnots-1)
  
  mydnorm = function(m) {
    dnorm(lats, mean=m, sd=sds)
  }
  Xi = sapply(means, mydnorm)
  
  if(intercept)
    return(cbind(1, Xi))
  else
    return(Xi)
}




