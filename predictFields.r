#N: number of pctRes observations to sample
#M: number of predictions
predictFields = function(seaDef=NULL, pctRes=NULL, model=NULL, krigeMLE=NULL, N=1000, 
                         M=3, predHull=NULL, real=4, mu=0, expand=1+1e-07) {
  resids = getPctResDat(seaDef)
  if(is.null(pctRes)) {
    pctRes = resids$pctRes
    pctResIn = FALSE
  }
  else {
    pctResIn = TRUE
  }
  
  #get parameters from model or krigeMLE:
  if(!is.null(model)) {
    aniso.pars=model$aniso.pars
    phi=model$phi
    lambda = 29606.39 #based on Krig function with some poor guesses at other parameters
    smoothness=model$kappa
  }
  else if(!is.null(krigeMLE)) {
    aniso.pars = c(pi/2, krigeMLE$cov.args.MLE$theta)
    phi = krigeMLE$cov.args.MLE$phi
    lambda = krigeMLE$lambda.MLE
    smoothness=3
  }
  else{
    stop("must either stupply model or krigeMLE")
  }
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  
  pctRes = resids$pctRes
  mask = resids$mask
  Mw = resids$Mw
  type = resids$type
  probs = resids$probs
  outThresh = quantile(abs(pctRes), .99, na.rm=TRUE)
  
  X = coordDat$X
  Y = coordDat$Y
  dx = coordDat$distPerCellX
  dy = coordDat$distPerCellY
  
  ##### make prediction points (approximately equidistant)
  #make prediction grid
  gridXY = make.surface.grid(list(x=X, y=Y))
  maskXY = gridXY[mask,]
  #   gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  #   gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewXY = make.surface.grid(list(x=gridNewX, y=gridNewY))
  
  #make concave hull around prediction mask to get final prediction points
  if(is.null(predHull))
    ashape = ahull(jitter(maskXY[,1], factor=.025), jitter(maskXY[,2], factor=.025), alpha=3000)
  else
    ashape = predHull
  indx=ashape$arcs[,"end1"]  
  hullPts <- maskXY[indx,]                  # extract the boundary points from maskXY
  maskNewXY = in.poly(gridNewXY, hullPts, convex.hull=FALSE)
  predXY = gridNewXY[maskNewXY,]
  predXY[,1] = predXY[,1]/aniso.pars[2] #scale for anisotropy
  
  #plot prediction points
  plotSubPoly(rbind(hullPts, hullPts[1,]), maskXY)
  
  ##### subsample data for training for prediction
  gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
  maskGridXYI = gridXYI[mask,]
  maskSampleI = sample(1:nrow(maskGridXYI), N)
  sampleXYI = maskGridXYI[maskSampleI,]
  #sampleRealI = sample(1:dim(pctRes)[3], N, replace=TRUE)
  sampleRealI = rep(real, N) # take all data from same realization
  resSelect = cbind(sampleXYI, sampleRealI)
  
  sampleX = X[sampleXYI[,1]]
  sampleY = Y[sampleXYI[,2]]
  sampleRes = pctRes[resSelect]
  
  #remove outliers
  if(!pctResIn) {
    outThresh = quantile(abs(pctRes), .99, na.rm=TRUE)
    outI = abs(sampleRes) > outThresh
    sampleX = sampleX[!outI]
    sampleY = sampleY[!outI]
    sampleRes = sampleRes[!outI]
  }
  else {
    outI = is.finite(sampleRes)
    sampleX = sampleX[outI]
    sampleY = sampleY[outI]
    sampleRes = sampleRes[outI]
  }
  
  sampleGrid = list(x=sampleX, y=sampleY)
  quilt.plot(sampleX, sampleY, sampleRes)
  
  #make geodata object with subsampled data
  reals = sampleRealI
  res = sampleRes
  Ys = sampleY
  Xs = sampleX/aniso.pars[2] #scale coordinates for anisotropy
  gridXY = cbind(Xs, Ys)
  #geoDat = as.geodata(cbind(Xs, Ys, res), realizations=reals)
  
  ##### make predictions
  cov.args=list(Covariance="Matern", smoothness=smoothness, 
                theta=phi, Distance="rdist")
  out = mKrig(gridXY, res, cov.args=cov.args, m=1, lambda=lambda)
  
  #make plots
#   surface(out)
#   
#   pred1 = predict(out, predXY)
#   quilt.plot(predXY, pred1)
#   
#   pred2 = predictSurface(out, predXY)
#   predSE = predictSurfaceSE(out, grid.list=predXY, chull.mask=rbind(hullPts, hullPts[1,]))
#   surface(pred2)
  
  test = sim.mKrig.approx(out, predictionPoints=predXY, M=M, 
                          gridExpansion=expand)
  test$Ensemble = test$Ensemble + mu
  
  plot.new()
  par(mfrow=c(2,2))
  minRes = min(min(sampleRes), min(test$Ensemble))
  maxRes = max(max(sampleRes), max(test$Ensemble))
  quilt.plot(sampleX, sampleY, sampleRes, main="Observations", 
             zlim=c(minRes, maxRes))
  for(i in 1:min(c(3, M))) {
    quilt.plot(test$predictionPoints, test$Ensemble[,i], 
               main=paste0("Prediction ", i), zlim=c(minRes, maxRes))
  }
  par(mfrow=c(1,1))
  
  # return(list(krigMod = out, pred=pred, predSE=predSE))
  
  return(list(krigMod = out, pred=test))
}

fitFields = function(seaDef=NULL, model, pctRes=NULL, N=500, predHull=NULL, optim.args=NULL) {
  resids = getPctResDat(seaDef)
  if(is.null(pctRes)) {
    pctRes = resids$pctRes
    pctResIn = FALSE
  }
  else {
    pctResIn = TRUE
  }
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  files = seaDef$files
  defDat = seaDef$dat
  
  SS = resids$SS
  mask = resids$mask
  Mw = resids$Mw
  type = resids$type
  probs = resids$probs
  
  files = files[-SS]
  defDat = defDat[,,-SS]
  
  X = coordDat$X
  Y = coordDat$Y
  dx = coordDat$distPerCellX
  dy = coordDat$distPerCellY
  
  ##### make prediction points (approximately equidistant)
  #make prediction grid
  gridXY = make.surface.grid(list(x=X, y=Y))
  maskXY = gridXY[mask,]
  #   gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  #   gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewXY = make.surface.grid(list(x=gridNewX, y=gridNewY))
  
  #make concave hull around prediction mask to get final prediction points
  if(is.null(predHull))
    ashape = ahull(jitter(maskXY[,1], factor=.025), jitter(maskXY[,2], factor=.025), alpha=3000)
  else
    ashape = predHull
  indx=ashape$arcs[,"end1"]  
  hullPts <- maskXY[indx,]                  # extract the boundary points from maskXY
  maskNewXY = in.poly(gridNewXY, hullPts, convex.hull=FALSE)
  predXY = gridNewXY[maskNewXY,]
  predXY[,1] = predXY[,1]
  
  #plot prediction points
  plotSubPoly(rbind(hullPts, hullPts[1,]), maskXY)
  
  ##### subsample data for training for prediction
  gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
  maskGridXYI = gridXYI[mask,]
  maskSampleI = sample(1:nrow(maskGridXYI), N)
  sampleXYI = maskGridXYI[maskSampleI,]
  #sampleRealI = sample(1:dim(pctRes)[3], N, replace=TRUE)
  sampleRealI = rep(4, N) # take all data from 
  resSelect = cbind(sampleXYI, sampleRealI)
  
  sampleX = X[sampleXYI[,1]]
  sampleY = Y[sampleXYI[,2]]
  sampleRes = pctRes[resSelect]
  
  #remove outliers
  if(!pctResIn) {
    outThresh = quantile(abs(pctRes), .99, na.rm=TRUE)
    outI = abs(sampleRes) > outThresh
    sampleX = sampleX[!outI]
    sampleY = sampleY[!outI]
    sampleRes = sampleRes[!outI]
  }
  else {
    outI = is.finite(sampleRes)
    sampleX = sampleX[outI]
    sampleY = sampleY[outI]
    sampleRes = sampleRes[outI]
  }
  
  sampleGrid = list(x=sampleX, y=sampleY)
  quilt.plot(sampleX, sampleY, sampleRes)
  
  #make geodata object with subsampled data
  reals = sampleRealI
  res = sampleRes
  Ys = sampleY
  Xs = sampleX
  gridXY = cbind(Xs, Ys)
  #geoDat = as.geodata(cbind(Xs, Ys, res), realizations=reals)
  
  #same as stationary.cov but theta controls only scaling for x points
  myStatCov <<- function(x1, x2=NULL, theta = 1, ...) {
    x1 = x1*theta
    return(do.call("stationary.cov", c(list(x1=x1, x2=x2), list(...))))
  }
  
  ##### make predictions
  aniso.pars = model$aniso.pars
  lambda = model$nugget/model$sigmasq # lambda is nugget/sill ?
  phi=model$phi
  smoothness = model$kappa
  cov.args=list(Covariance="Matern", smoothness=smoothness, 
                Distance="rdist")
  cov.args.guess=list(phi=phi, theta=model$aniso.pars[2])
  
  if(is.null(optim.args)) {
    optim.args = list(method = "BFGS", 
                      control=list(fnscale = -1, 
                                   ndeps = rep(log(1.1), length(cov.args.guess)+1), 
                                   reltol=1e-08, maxit=50))
    obj1 = mKrig.MLE.joint(gridXY, res, cov.fun="myStatCov", 
                          cov.args=cov.args, optim.args=optim.args, 
                          cov.params.guess=cov.args.guess, 
                          lambda.guess=lambda, verbose=TRUE)
  }
  else {
    obj1 = mKrig.MLE.joint(gridXY, res, cov.fun="myStatCov", 
                          cov.args=cov.args, 
                          cov.params.guess=cov.args.guess, 
                          lambda.guess=lambda, verbose=TRUE)
  }
  
  lambda.grid=10^seq(-4,4, length=9)
  theta.grid=seq(1, 1.5, length=4)
  theta.grid=rep(theta.grid, length(lambda.grid))
  phi.grid = seq(phi/2, phi*2, length=4)
  par.grid=make.surface.grid(list(theta=theta.grid, phi=phi.grid, lambda=lambda.grid))
  lambda.grid=par.grid[,3]
  par.grid=par.grid[,1:2]
  obj2 = mKrig.MLE(gridXY, res, cov.fun="myStatCov", cov.args=cov.args, 
                   par.grid=par.grid, lambda=lambda.grid, 
                   lambda.profile=FALSE, verbose=TRUE)
  
  return(list(obj1=obj1, obj2=obj2))
}

fitFields.multi = function(seaDef, model, pctRes=NULL, N=500, predHull=NULL, optim.args=NULL) {
  #NOTE: assumes data has constant mean field, no lat/lon regression terms
  
  resids = getPctResDat(seaDef)
  if(is.null(pctRes)) {
    pctRes = resids$pctRes
    pctResIn = FALSE
  }
  else {
    pctResIn = TRUE
  }
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  files = seaDef$files
  defDat = seaDef$dat
  
  SS = resids$SS
  mask = resids$mask
  Mw = resids$Mw
  type = resids$type
  probs = resids$probs
  
  files = files[-SS]
  defDat = defDat[,,-SS]
  
  X = coordDat$X
  Y = coordDat$Y
  dx = coordDat$distPerCellX
  dy = coordDat$distPerCellY
  
  ##### make prediction points (approximately equidistant)
  #make prediction grid
  gridXY = make.surface.grid(list(x=X, y=Y))
  maskXY = gridXY[mask,]
  #   gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  #   gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewXY = make.surface.grid(list(x=gridNewX, y=gridNewY))
  
  #make concave hull around prediction mask to get final prediction points
  if(is.null(predHull))
    ashape = ahull(jitter(maskXY[,1], factor=.025), jitter(maskXY[,2], factor=.025), alpha=3000)
  else
    ashape = predHull
  indx=ashape$arcs[,"end1"]  
  hullPts <- maskXY[indx,]                  # extract the boundary points from maskXY
  maskNewXY = in.poly(gridNewXY, hullPts, convex.hull=FALSE)
  predXY = gridNewXY[maskNewXY,]
  predXY[,1] = predXY[,1]
  
  #plot prediction points
  plotSubPoly(rbind(hullPts, hullPts[1,]), maskXY)
  
  ##### subsample data for training for prediction
  gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
  maskGridXYI = gridXYI[mask,]
  maskSampleI = sample(1:nrow(maskGridXYI), N)
  sampleXYI = maskGridXYI[maskSampleI,]
  sampleRealI = sample(1:dim(pctRes)[3], N, replace=TRUE)
  #sampleRealI = rep(4, N) # take all data from 
  resSelect = cbind(sampleXYI, sampleRealI)
  
  sampleX = X[sampleXYI[,1]]
  sampleY = Y[sampleXYI[,2]]
  sampleRes = pctRes[resSelect]
  
  #remove outliers
  if(!pctResIn) {
    outThresh = quantile(abs(pctRes), .99, na.rm=TRUE)
    outI = abs(sampleRes) > outThresh
    sampleX = sampleX[!outI]
    sampleY = sampleY[!outI]
    sampleRes = sampleRes[!outI]
    sampleReals = sampleRealI[!outI]
  }
  else {
    outI = is.finite(sampleRes)
    sampleX = sampleX[outI]
    sampleY = sampleY[outI]
    sampleRes = sampleRes[outI]
    sampleReals = sampleRealI[outI]
  }
  
  sampleGrid = list(x=sampleX, y=sampleY)
  quilt.plot(sampleX, sampleY, sampleRes)
  
  #make geodata object with subsampled data
  reals = sampleReals
  res = sampleRes
  Ys = sampleY
  Xs = sampleX
  gridXY = cbind(Xs, Ys)
  #geoDat = as.geodata(cbind(Xs, Ys, res), realizations=reals)
  
  #same as stationary.cov but theta controls only scaling for x points
  myStatCov <<- function(x1, x2=NULL, theta = 1, ...) {
    x1 = x1*theta
    return(do.call("stationary.cov", c(list(x1=x1, x2=x2), list(...))))
  }
  
  ##### make predictions (commented out portions would be used in joint
  ##### optimization)
  aniso.pars = model$aniso.pars
  #lambda = model$nugget/model$sigmasq # lambda is nugget/sill ?
  phi=model$phi
  smoothness = model$kappa
  cov.args=list(Covariance="Matern", smoothness=smoothness, 
                Distance="rdist")
  #cov.args.guess=list(phi=phi, theta=model$aniso.pars[2])
  
  lambda.grid=10^seq(-4,4, length=9)
  theta.grid=seq(1, 1.5, length=4)
  theta.grid=rep(theta.grid, length(lambda.grid))
  phi.grid = seq(phi/2, phi*4, length=8)
  par.grid=make.surface.grid(list(theta=theta.grid, phi=phi.grid, lambda=lambda.grid))
  lambda.grid=par.grid[,3]
  par.grid=par.grid[,1:2]
  obj = mKrig.MLE.multi(gridXY, res, reals, cov.fun="myStatCov", cov.args=cov.args, 
                   par.grid=par.grid, lambda=lambda.grid, m=1)
  
  return(obj)
}





