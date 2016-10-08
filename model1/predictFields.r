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

#N: number of pctRes observations to sample
#M: number of predictions
myPredictFields = function(seaDef=NULL, pctRes=NULL, model=NULL, krigeMLE=NULL, N=1000, 
                         M=3, predHull=NULL, real=4, mu=0) {
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
    # aniso.pars = c(pi/2, krigeMLE$cov.args.MLE$theta)
    aniso.pars = c(pi/2, 1)
    phi = krigeMLE$cov.args.MLE[[1]]
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
  
  ##### make predictions using Cholesky decomposition
  #compute covariance matrix
  cov.args=list(Covariance="Matern", smoothness=smoothness, 
                theta=phi, Distance="rdist")
  out = mKrig(gridXY, res, cov.args=cov.args, m=1, lambda=lambda)
  marginalVar = out$rho.MLE
  covMat = stationary.cov(gridXY, Covariance = "Matern", 
                          smoothness=smoothness, theta=phi, 
                          Distance="rdist")
  covMat = covMat*marginalVar
  
  #do cholesky decomposition
  L = t(chol(covMat))
  
  #make predictions
  preds = matrix(NA, nrow=nrow(L), ncol=M)
  for(i in 1:M) {
    z = rnorm(ncol(L))
    preds[,i] = L %*% z
  }
  preds = preds + mu
  
  #plot predictions
  plot.new()
  par(mfrow=c(2,2))
  minRes = min(min(sampleRes), min(preds))
  maxRes = max(max(sampleRes), max(preds))
  quilt.plot(sampleX, sampleY, sampleRes, main="Observations", 
             zlim=c(minRes, maxRes))
  for(i in 1:min(c(3, M))) {
    quilt.plot(gridXY, preds[,i], 
               main=paste0("Prediction ", i), zlim=c(minRes, maxRes))
  }
  par(mfrow=c(1,1))
  
  #return results
  # return(list(krigMod = out, pred=pred, predSE=predSE))
  
  return(pred=list(gridXY=gridXY, preds))
}

#N: number of pctRes observations to sample
#M: number of predictions
#model: a list with phi (range), kappa (smoothness), tausq (nugget), and sigmasq (sill)
myPredictFields2 = function(seaDef=NULL, pctRes=NULL, model=NULL, N=1000, 
                           M=3, predHull=NULL, real=4, mu=0) {
  resids = getPctResDat(seaDef)
  if(is.null(pctRes)) {
    pctRes = resids$pctRes
    pctResIn = FALSE
  }
  else {
    pctResIn = TRUE
  }
  
  #get parameters from model
  phi=model$phi
  smoothness=model$kappa
  marginalVar = model$sigmasq
  nugget = model$tausq
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  
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
  predXY[,1] = predXY[,1]
  
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
  
  #make geodata object with subsampled data
  reals = sampleRealI
  res = sampleRes
  Ys = sampleY
  Xs = sampleX
  gridXY = cbind(Xs, Ys)
  #geoDat = as.geodata(cbind(Xs, Ys, res), realizations=reals)
  
  ##### make predictions using Cholesky decomposition
  #compute covariance matrix
  covMat = stationary.cov(gridXY, Covariance = "Matern", 
                          smoothness=smoothness, theta=phi, 
                          Distance="rdist")
  covMat = covMat*marginalVar + diag(nugget, nrow=nrow(covMat))
  
  #do cholesky decomposition
  L = t(chol(covMat))
  
  #make predictions
  preds = matrix(NA, nrow=nrow(L), ncol=M)
  for(i in 1:M) {
    z = rnorm(ncol(L))
    preds[,i] = L %*% z
  }
  preds = preds + mu
  
  #plot predictions
  par(mfrow=c(2,2))
  minRes = min(min(sampleRes), min(preds))
  maxRes = max(max(sampleRes), max(preds))
  quilt.plot(sampleX, sampleY, sampleRes, main="Observations", 
             zlim=c(minRes, maxRes))
  for(i in 1:min(c(3, M))) {
    quilt.plot(gridXY, preds[,i], 
               main=paste0("Prediction ", i), zlim=c(minRes, maxRes))
  }
  par(mfrow=c(1,1))
  
  #return results
  # return(list(krigMod = out, pred=pred, predSE=predSE))
  
  return(pred=list(gridXY=gridXY, preds))
}

######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
######################################################################
#fitting functions


fitFields = function(seaDef=NULL, model, pctRes=NULL, N=500, predHull=NULL, 
                     optim.args=NULL, real=4) {
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
  sampleRealI = rep(real, N) # take all data from realization "real"
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
  
#   if(is.null(optim.args)) {
#     optim.args = list(method = "BFGS", 
#                       control=list(fnscale = -1, 
#                                    ndeps = rep(log(1.1), length(cov.args.guess)+1), 
#                                    reltol=1e-08, maxit=50))
#     obj1 = mKrig.MLE.joint(gridXY, res, cov.fun="myStatCov", 
#                           cov.args=cov.args, optim.args=optim.args, 
#                           cov.params.guess=cov.args.guess, 
#                           lambda.guess=lambda, verbose=TRUE)
#   }
#   else {
#     obj1 = mKrig.MLE.joint(gridXY, res, cov.fun="myStatCov", 
#                           cov.args=cov.args, 
#                           cov.params.guess=cov.args.guess, 
#                           lambda.guess=lambda, verbose=TRUE)
#   }
  
  lambda.grid=10^seq(-4,4, length=31)
  theta.grid = seq(phi/2, phi*2, length=20)
  par.grid=make.surface.grid(list(theta=theta.grid, lambda=lambda.grid))
  lambda.grid=par.grid[,2]
  par.grid=make.surface.grid(list(theta=par.grid[,1]))
  obj2 = mKrig.MLE(gridXY, res, cov.args=cov.args, 
                   par.grid=par.grid, lambda=lambda.grid, 
                   lambda.profile=FALSE)
  
  lambda.grid=10^seq(-4,4, length=31)
  theta.grid = seq(phi/2, phi*2, length=20)
  par.grid=make.surface.grid(list(theta=theta.grid, lambda=log10(lambda.grid)))
  quilt.plot(par.grid, obj2$summary[,2], main="Negative Log Likelihood", 
             xlab="Theta", ylab="log10(Lambda)")
  
  return(obj2)
}

fitFields.multi = function(seaDef, model, pctRes=NULL, N=500, 
                           predHull=NULL, optim.args=NULL) {
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
                Distance="rdist", realizations=reals)
  #cov.args.guess=list(phi=phi, theta=model$aniso.pars[2])
  
#   lambda.grid=10^seq(-4,4, length=9)
#   theta.grid=seq(1, 1.5, length=4)
#   # theta.grid=rep(theta.grid, length(lambda.grid))
#   phi.grid = seq(phi/2, phi*4, length=8)
#   par.grid=make.surface.grid(list(theta=theta.grid, range=phi.grid, lambda=lambda.grid))
#   lambda.grid=par.grid[,3]
#   par.grid=par.grid[,1:2]
#   obj = mKrig.MLE.multi(gridXY, res, reals, cov.fun="myStatCov", cov.args=cov.args, 
#                    par.grid=par.grid, lambda=lambda.grid, m=1)
  
  
  
  lambda.grid=10^seq(-4,4, length=31)
  theta.grid = seq(phi/2, phi*4, length=20)
  par.grid=make.surface.grid(list(theta=theta.grid, lambda=lambda.grid))
  lambda.grid=par.grid[,2]
  par.grid=make.surface.grid(list(range=par.grid[,1]))
#   obj = mKrig.MLE.multi(gridXY, res, reals, cov.fun="stationary.cov", cov.args=cov.args, 
#                         par.grid=par.grid, lambda=lambda.grid, m=1)
  obj = mKrig.MLE(gridXY, res, cov.fun="stationary.cov.multi", cov.args=cov.args, 
                        par.grid=par.grid, lambda=lambda.grid, m=1, find.trA=FALSE)
  
  lambda.grid=10^seq(-4,4, length=31)
  theta.grid = seq(phi/2, phi*10, length=20)
  par.grid=make.surface.grid(list(theta=theta.grid, lambda=log10(lambda.grid)))
  quilt.plot(par.grid, obj$summary[,2], main="Log Likelihood", 
             xlab="Theta", ylab="log10(Lambda)")
  
  return(obj)
}

#model: a list with phi (range), kappa (smoothness), tausq (nugget), and sigmasq (sill)
plotSample = function(seaDef=NULL, pctRes=NULL, model=NULL, N=1000, 
                      predHull=NULL, mu=0) {
  regressOut = getPctResDat(seaDef)
  resids = getResDat(seaDef)
  if(is.null(pctRes)) {
    pctRes = regressOut$pctRes
    pctResIn = FALSE
  }
  else {
    pctResIn = TRUE
  }
  
  #get parameters from model
  phi=model$phi
  smoothness=model$kappa
  sill = model$sigmasq
  nugget = model$tausq
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  
  mask = regressOut$mask
  outThresh = quantile(abs(pctRes), .99, na.rm=TRUE)
  
  X = coordDat$X
  Y = coordDat$Y
  
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
  
  #####get N samples from each realization and all 15*5 simulations
  allGridXY = array(dim=c(N, 2, 90))
  allResSelect = array(dim=c(N, 3, 90))
  allPctRes = matrix(nrow=N, ncol=90)
  for(i in 1:90) {
    #print progress
    if(i == 1)
      print("Sampling realization data")
    if(i == 16) 
      print("Generating simulations")
    
    #sample data from realizations
    par(mfrow=c(1,1))
    out = samplePctRes(seaDef, pctRes, N=N, real=((i-1) %% 15)+1, thresh=Inf)
    
    #get XY data, resI (which will allow us to get regression predictions
    #at sample locations)
    allGridXY[,,i] = out$gridXY
    allResSelect[,,i] = out$resSelect
    
    #res data will depend on whether this is a simulation or a realization
    #(15 of each)
    if(i <= 15) {
      #if realization, just use data from samplePctRes
      allPctRes[,i] = out$res
    }
    else {
      #if simulation, generate data
      covMat = stationary.cov(out$gridXY, Covariance = "Matern", 
                              smoothness=smoothness, theta=phi, 
                              Distance="rdist")
      covMat = covMat*sill + diag(nugget, nrow=nrow(covMat))
      
      #do cholesky decomposition
      L = t(chol(covMat))
      
      #make predictions
      z = rnorm(ncol(L))
      preds = L %*% z
      allPctRes[,i] = preds
    }
  }
  
  #plot predictions on percent difference from regression scale
  #(6x5 grid, top 15 are observations, bottom 15 are simulations)
  par(mfrow=c(6, 5))
  minRes = max(c(min(allPctRes), -1))
  maxRes = min(c(max(allPctRes), 1))
  #plot observations
  for(i in 1:15) {
    real = i
    
    quilt.plot(allGridXY[,,i], allPctRes[,i], main=paste0("Observation ", i), 
               zlim=c(minRes, maxRes))
  }
  #plot simulations
  for(i in 16:30) {
    quilt.plot(allGridXY[,,i], allPctRes[,i], 
               main=paste0("Prediction ", i-15), zlim=c(minRes, maxRes))
  }
  
  #####now plot actual data versus predictions
  
  #restructure pctRes so we can mask certain points and make predictions
  pctRes = array(pctRes, dim=c(dim(pctRes)[1]*dim(pctRes)[2], dim(pctRes)[3]))
  
  #get regression model and variables
  mod = regressOut$mod
  Mw = regressOut$Mw
  type = regressOut$type
  probs = regressOut$probs
  Mw102 = 10^Mw
  
  #get regression predictions
  testIn = data.frame(Mw102=Mw102, type2=type)
  testOut = predict(mod, testIn)
  preds = array(0, dim=dim(pctRes))
  preds[mask, ] = t(testOut)
  
  #restructure predictions so we can get predictions for each simulation and
  #realization
  preds = array(preds, dim=dim(regressOut$pctRes))
  
  #get all regression predictions for each sample location in 
  #simulations and realizations
  allRegress = matrix(nrow=N, ncol=90)
  for(i in 1:90) {
    allRegress[,i] = preds[allResSelect[,,i]]
  }
  
  #get all final predictions and observations for each sample point
  allSeaDefs = allRegress + allRegress*(allPctRes + mu)
  
  #calculate data ranges (make sure paired observations and simulations 
  #have same color range)
  ranges = matrix(ncol=2, nrow=15)
  allRanges = t(apply(allSeaDefs, 2, range))
  for(i in 1:15) {
    ranges[i, 1] = min(c(allRanges[i, 1], allRanges[i+15, 1]))
    ranges[i, 2] = max(c(allRanges[i, 2], allRanges[i+15, 2]))
  }
  
  #plot observations
  for(i in 1:15) {
    real = i
    
    quilt.plot(allGridXY[,,i], allSeaDefs[,i], main=paste0("Observed SeaDef ", i), zlim=ranges[i,])
  }
  #plot simulations
  for(i in 16:30) {
    quilt.plot(allGridXY[,,i], allSeaDefs[,i], zlim=ranges[i-15,], 
               main=paste0("Predicted SeaDef ", i-15))
  }
  
  ##### now plot same as above but 5 simulations per observation
  par(mfrow=c(2, 3))
  
  #calculate data ranges (make sure paired observations and simulations 
  #have same color range)
  ranges = matrix(ncol=2, nrow=15)
  for(i in 1:15) {
    group = seq(i, i+15*5, by=15)
    ranges[i, 1] = min(allRanges[group, 1])
    ranges[i, 2] = max(allRanges[group, 2])
  }
  
  #plot observations versus corresponding simulations
  for(i in 1:15) {
    group = seq(i, i+15*5, by=15)
    
    quilt.plot(allGridXY[,,i], allSeaDefs[,i], main=paste0("Observed SeaDef ", i), zlim=ranges[i,])
    
    for(j in 2:length(group)) {
      quilt.plot(allGridXY[,,group[j]], allSeaDefs[,group[j]], zlim=ranges[i,], 
                 main=paste0("Predicted SeaDef"))
    }
  }
  
  par(mfrow=c(1,1))
  
  invisible(NULL)
}

#sample data from realization
samplePctRes = function(seaDef, pctRes, N=1000, real=1, thresh=Inf) {
  #oversample since we end up getting rid of some points
  actualN = N
  N = N*1.5
  
  #get X and Y coordinates
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  X = coordDat$X
  Y = coordDat$Y
  
  # get domain mask
  resids = getPctResDat(seaDef)
  mask = resids$mask
  
  ##### subsample data for this realization
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
  
  #remove outliers, non-finite data
  outI = !(is.na(sampleRes) | (abs(sampleRes) > thresh))
  sampleX = sampleX[outI]
  sampleY = sampleY[outI]
  sampleRes = sampleRes[outI]
  resSelect=resSelect[outI,]
  
  #final subsampled data
  res = sampleRes[1:actualN]
  Ys = sampleY[1:actualN]
  Xs = sampleX[1:actualN]
  gridXY = cbind(Xs, Ys)
  resSelect=resSelect[1:actualN,]
  
  return(list(gridXY = gridXY, res=res, resSelect=resSelect))
}

