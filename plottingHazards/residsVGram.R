getVGram = function(seaDef, pctRes, N=1500) {
  resids = getPctResDat(seaDef)
  
  #get coordinates of pctRes in meters
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  X = coordDat$X
  Y = coordDat$Y
  gridXY = make.surface.grid(list(x=X, y=Y))
  
  #must be able to remove points outside of region of interest with mask
  mask = resids$mask
  
  #sample data associated with locations, residuals, and realizations
  gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
  maskGridXYI = gridXYI[mask,]
  maskSampleI = sample(1:nrow(maskGridXYI), N)
  sampleXYI = maskGridXYI[maskSampleI,]
  sampleRealI = sample(1:dim(pctRes)[3], N, replace=TRUE)
  resSelect = cbind(sampleXYI, sampleRealI)
  
  sampleX = X[sampleXYI[,1]]
  sampleY = Y[sampleXYI[,2]]
  sampleRes = pctRes[resSelect]
  sampleReal = sampleRealI
  sampleGrid = matrix(c(sampleX, sampleY), ncol=2)
  
  #concatenate variogram results, only comparing samples from same realization 
  #to each other
  for(r in 1:dim(pctRes)[3]) {
    thisReal = sampleReal == r
    thisGridXY = sampleGrid[thisReal,]
    thisRes = sampleRes[thisReal]
    
    #only compare consecutive points so our point cloud isn't too dense
    n = sum(thisReal)
    id = matrix(c(1:n, 2:n, 1), ncol=2)
    
    if(r == 1) {
      VG = vgram(thisGridXY, thisRes, id=id)
    }
    else {
      tmp = vgram(thisGridXY, thisRes, id=id)
      VG$d = c(VG$d, tmp$d)
      VG$vgram = c(VG$vgram, tmp$vgram)
    }
  }
  
  return(VG)
}

#nugget is tausq
#sill is sigmasq
#phi input is equivalent to range in Matern?

#fields:
# Exponential:
#   
#   phi* exp( -d/range)
# 
# Matern:
#   nu = smoothness, con is normalizing constant, phi = marginal variance
#   phi*con*(d\^nu) * besselK(d , nu )

#geoR:
#Matern correlation:
#ρ(h) = (1/(2^(κ-1) * Γ(κ))) * ((h/φ)^κ) * K_{κ}(h/φ)ρ(h) = (1/(2^(κ-1) * Γ(κ))) * ((h/φ)^κ) * K_{κ}(h/φ)
#note: this is the same is fields with nu=κ, d=h/φ?
plotVGram = function(VG, kappa=3, phi, sill, nugget=0, ylim=NULL) {
  #f is our variogram fit to plot
  f = function(d) {
    d = d/phi
    nu = kappa
    sill - (Matern(d, nu=nu, phi=(sill - nugget)) + nugget)
  }
  
  xs = seq(from=min(VG$d), to=max(VG$d), length=100)
  if(is.null(ylim))
    plot(VG, pch=19, cex=.5)
  else
    plotVGMean(VG, pch=19, cex=.5, ylim=ylim)
  lines(xs, f(xs), col="blue")
}

# sillInit: sd(apply(pctRes, 3, median, na.rm=TRUE))
# rangeInit = max(VG$d)/2
# nuggetInit: 0
# 
fitVGram = function(VG, sillInit, rangeInit = max(VG$d), nuggetInit=0, kappa=3) {
  #get data to fit and initial parameter guesses
  ys = VG$vgram
  ds = VG$d
  
  #make Matern variance function
  f = function(d, phi=rangeInit, kap=kappa, s=sillInit, n=nuggetInit) {
    #if parameters don't make sense, return bad number
    if(phi <= 0 || s < n || n < 0) {
      return(rep(-10, length(d)))
    }
    d = d/phi
    s - (Matern(d, nu=kap, phi=(s - n)) + n)
  }
  
#   fit = nls(ys ~ f(ds, phi=range, s=sill, n=nugget), 
#       start=list(sill=sillInit, nugget=nuggetInit, range=rangeInit), 
#       lower=list(sill=0.000001, nugget=0, range=0.000001), algorithm="port",
#       upper=list(sill=5, nugget=1, range=1500000), 
#       control=nls.control(maxiter=100, tol=10, printEval=TRUE, 
#                           warnOnly=TRUE), 
#       trace=TRUE)
  fit = nls(ys ~ f(ds, phi=range, s=sill, n=nugget), 
            start=list(sill=sillInit, nugget=nuggetInit, range=rangeInit), 
            control=nls.control(maxiter=100, tol=1, printEval=TRUE, 
                                warnOnly=TRUE), trace=TRUE)
  
  summary(fit)
  coefs = coef(fit)
  sill = coefs[1]
  nugget = coefs[2]
  range = coefs[3]
  
  #plot variogram
  xs = seq(0, max(VG$d), length=500)
  plotVGMean(VG, ylim=c(0, .1), 
             main="Empirical Matern Variogram and Fit")
  lines(xs, sqrt(f(xs, phi=range, s=sill, n=nugget)), col="green")
  
  return(fit)
}

#compute least squares using grid search
gridLS = function(xs, ys, f, paramGridList) {
  
  MSE = c()
  for(r in 1:nrow(paramGridList)) {
    params = as.list(paramGridList[r,])
    preds = do.call(f, c(list(xs), params))
    MSE[r] = mean((preds - ys)^2)
  }
  minI = which.min(MSE)
  summTable = cbind(MSE, paramGridList)
  bestParams = paramGridList[minI,]
  bestMSE = MSE[minI]
  
  return(list(summTable = summTable, OLSE=bestParams, OLSE.MSE=bestMSE))
}

#compute maximum likelihood using grid search
#NOTE: paramGridList should have the same number of rows as the length of input 
#sill.grid and nugget.grid
gridML = function(xs, ys, realizations, sill.grid, nugget.grid, paramGridList, ...) {
  if(!is.matrix(paramGridList)) {
    name = names(paramGridList)
    paramGridList = matrix(paramGridList, nrow=length(paramGridList))
    colnames(paramGridList) = name
  }
  
  updateIters = round(quantile(1:nrow(paramGridList), probs=seq(.01, .99, by=.01)))
  count = 0
  logLik = c()
  for(r in 1:nrow(paramGridList)) {
    sill = sill.grid[r]
    nugget = nugget.grid[r]
    params = as.list(paramGridList[r,])
    logLik[r] = do.call("likeEval.multi", c(list(x1=xs, y=ys, realizations=realizations, 
                                                 nugget=nugget, sill=sill), params, list(...)))
    
    #give progress updates
    if(any(updateIters == r)) {
      count = count + 1
      print(paste0(count, " percent complete, ", r, "/", nrow(paramGridList), " iterations complete."))
      print(paste0("Current MLE:"))
      maxI = which.max(logLik[1:r])
      print(paste0("logLik: ", logLik[maxI], " sill: ", sill.grid[maxI], " nugget: ", 
                   nugget.grid[maxI], " range: ", paramGridList[maxI,]))
    }
  }
  maxI = which.max(logLik)
  summTable = cbind(logLik, sill.grid, nugget.grid, paramGridList)
  bestParams = c(sill=sill.grid[maxI], nugget=nugget.grid[maxI], paramGridList[maxI,])
  bestLogLik = logLik[maxI]
  
  return(list(summTable = summTable, MLE=bestParams, MLE.logLik=bestLogLik))
}
