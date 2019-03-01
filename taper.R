
# ideas for tapers
# h = function(x) {ifelse(x <= 0, 0, exp(-1/x))}
# g = function(x, a=1) { a*h(x)/(a*h(x) + h(1-x))}
# 
# gInv = function(x, a=1) {
#   ans = x
#   ans[x == 0] = 0
#   ans[x == 1] = 1
#   inds = (x != 0) & (x != 1)
#   
#   term1 = (log(-(a*x - a)/x) + 2)/(2*log(-(a*x - a)/x))
#   term2 = sqrt(log(-(a*x - a)/x)^2 + 4)/(2*log(-(a*x - a)/x))
#   ans[inds] = term1[inds] - term2[inds]
#   
#   ans
# }

# Randy's exponential taper, renormalized
# taper = function(d, lam=1, dStar=1) {
#   ans = (1 - exp(-lam*(dStar-d)/dStar))/(1 - exp(-lam))
#   ans[x > dStar] = 0
#   return(ans)
# }

# power exponential taper, renormalized
taper = function(d, lambda=1, alpha=2, dStar=21000, normalize=TRUE, approxNear0=TRUE) {
  scaledD = abs(d/dStar)^alpha * lambda^alpha
  if(normalize) {
    ans = 1 - (1 - exp(-scaledD))/(1 - exp(-lambda^alpha))
    ans[d > dStar] = 0
  }
  else
    ans = exp(-lambda^2*d^2)
  # ans[abs(scaledD) < 1e-17] = 0 # account for numerical instability near 0 (usually only occurs if something's wrong)
  ans[d <= 0] = 1
  if(approxNear0 && normalize) {
    smallLam = abs(lambda) < .0000005
    ans[smallLam] = 1 - (d[smallLam]/dStar)^2 # numerical instability so take the limit as lambda to 0
  }
  return(ans)
}

# only for alpha=2.  Returns an n x p matrix of derivatives for n depths and p basis coefficients.
# Xi is the n x p taper basis matrix
taperGrad = function(d, Xi, lambda=1, dStar=21000, normalize=TRUE, diffGPSTaper=FALSE, diffXi2=FALSE, Xi2=Xi) {
  dRatio = d/dStar
  scaledD = dRatio * lambda
  expD = exp(-scaledD^2) # expD is the unnormalized taper function
  
  if(normalize) {
    expLam = exp(-lambda^2)
    pVec = ((1 - expD)*expLam*2*lambda)/(1-expLam)^2    -    expD*2*scaledD*dRatio/(1 - expLam)
    pVec[d > dStar] = 0
  }
  else {
    pVec = -2*d^2*lambda*exp(-d^2*lambda^2)
  }
  
  # compute diag(pVec) * Xi
  derivMat = sweep(Xi, 1, pVec, "*")
  
  # if GPS data has an adjusted taper with extra parameters, modify the gradient here
  if(diffGPSTaper && !diffXi2)
    return(cbind(derivMat, -derivMat))
  else if(diffGPSTaper)
    return(cbind(derivMat, -sweep(Xi2, 1, pVec, "*")))
  else
    return(derivMat)
}

# solves the taper equation for d, the depth, given the taper value
getFracTaperDepth = function(lambdas, depthFrac, alpha=2, dStar=21000, normalize=TRUE, approxNear0=TRUE) {
  if(!normalize || alpha != 2) {
    stop("non-normalized taper and non-standard power not yet implemented")
  }
  else {
    ans = sqrt(-dStar^2 * log((1 - exp(-lambdas^2)) * depthFrac + exp(-lambdas^2)) / lambdas^2)
    
    if(approxNear0) {
      smallLam = abs(lambdas) < .0000005
      ans[smallLam] = sqrt(-dStar^2 * (depthFrac - 1))
    }
  }
  
  ans
}


