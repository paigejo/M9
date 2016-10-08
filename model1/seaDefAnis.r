library(geoR)
library(alphahull)

##### function goal: plot anisotropy in variograms
plotAnisVarios = function(seaDef, N=1000) {
  # load seafloor deformations, residuals from regression, coord distance information
  resids = getResDat(seaDef)
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  files = seaDef$files
  defDat = seaDef$dat
  
  SS = resids$SS
  resDat = resids$resids
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
  
  pctRes = resDat/defDat
  pctRes[defDat == 0] = 0 #both resDat and defDat are 0, so set error to
  
  ##### fit anisotropic model
  
  #subset data (how to deal with multiple fields? subset random data from random fields?)
  #For now, I will estimate a single residual field, no matter the realization
  
  # subsample residual data from each realization
  
  for(r in 1:dim(resDat)[3]) {
#     gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
#     maskGridXYI = gridXYI[mask,]
#     maskSampleI = sample(1:nrow(maskGridXYI), N)
#     sampleXYI = maskGridXYI[maskSampleI,]
#     sampleRealI = sample(1:dim(resDat)[3], N, replace=TRUE)
#     resSelect = cbind(sampleXYI, sampleRealI)
#     
#     sampleX = X[sampleXYI[,1]]
#     sampleY = Y[sampleXYI[,2]]
#     sampleRes = resDat[resSelect]
#     
#     sampleGrid = list(x=sampleX, y=sampleY)
#     quilt.plot(sampleX, sampleY, sampleRes)
    
    #now try same thing for but pct error
    #(this looks much better.  Surprisingly, correlation length is much longer north to south)
    
    gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
    maskGridXYI = gridXYI[mask,]
    maskSampleI = sample(1:nrow(maskGridXYI), N)
    sampleXYI = maskGridXYI[maskSampleI,]
    sampleRealI = rep(r, N)
    resSelect = cbind(sampleXYI, sampleRealI)
    
    sampleX = X[sampleXYI[,1]]
    sampleY = Y[sampleXYI[,2]]
    sampleRes = pctRes[resSelect]
    
    sampleGrid = list(x=sampleX, y=sampleY)
    quilt.plot(sampleX, sampleY, sampleRes)
    
    # fit the model
    reals = sampleRealI
    Xs = sampleX
    Ys = sampleY
    mdist = max(Ys) - min(Ys)
    res = sampleRes
    geoDat = as.geodata(cbind(Xs, Ys, res), realizations=reals)
    
#     vario = variog(geoDat)
#     v.env <- variog.mc.env(geoDat, obj.variog=vario)
#     vario4 = variog4(geoDat, direction=c(0, pi/2))
#     plot(vario4, env=v.env, main=paste0("Omni-Directional Variogram For Realization ", r))
    
    par(mfrow=c(2,2))
    varioOmni = variog(geoDat, max.dist=mdist, estimator.type="classical")
    varioNS = variog(geoDat, max.dist=mdist, direction=0)
    varioEW = variog(geoDat, max.dist=mdist, direction=pi/2)
    plot(varioOmni, main=paste0("Omni-Directional Variogram For Realization ", r))
    par(mfrow=c(1,1))
  }
}

##### function goal: fit anisotropic model to seadef data
fitAnisModel = function(seaDef, pctRes = NULL, N=100) {
  # load seafloor deformations, residuals from regression, coord distance information
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
  
  gridXY = make.surface.grid(list(x=X, y=Y))
  maskXY = gridXY[mask,]
  
  ##### fit anisotropic model
  
  #subset data (how to deal with multiple fields? subset random data from random fields?)
  #For now, I will estimate a single residual field, no matter the realization
  
  # subsample residual data from all realizations
#   gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
#   maskGridXYI = gridXYI[mask,]
#   maskSampleI = sample(1:nrow(maskGridXYI), N)
#   sampleXYI = maskGridXYI[maskSampleI,]
#   sampleRealI = sample(1:dim(resDat)[3], N, replace=TRUE)
#   resSelect = cbind(sampleXYI, sampleRealI)
#   
#   sampleX = X[sampleXYI[,1]]
#   sampleY = Y[sampleXYI[,2]]
#   sampleRes = resDat[resSelect]
#   
#   sampleGrid = list(x=sampleX, y=sampleY)
#   quilt.plot(sampleX, sampleY, sampleRes, zlim=quantile(sampleRes, c(.05, .95)), nx=200, ny=200)
#   quilt.plot(sampleX, sampleY, sampleRes, nx=200, ny=200)
  
  #now try same thing for but pct error
  #(this looks much better.  Surprisingly, correlation length is much longer north to south)
  
  gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
  maskGridXYI = gridXYI[mask,]
  maskSampleI = sample(1:nrow(maskGridXYI), N)
  sampleXYI = maskGridXYI[maskSampleI,]
  sampleRealI = sample(1:dim(pctRes)[3], N, replace=TRUE)
  #sampleRealI = rep(5, N)
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
    sampleRealI = sampleRealI[!outI]
  }
  else {
    outI = is.finite(sampleRes)
    sampleX = sampleX[outI]
    sampleY = sampleY[outI]
    sampleRes = sampleRes[outI]
    sampleRealI = sampleRealI[outI]
  }
  
  sampleGrid = list(x=sampleX, y=sampleY)
  #quilt.plot(sampleX, sampleY, sampleRes, zlim=quantile(sampleRes, c(.05, .95)))
  quilt.plot(sampleX, sampleY, sampleRes)
  
  # fit the model
  reals = sampleRealI
  Xs = sampleX
  Ys = sampleY
  res = sampleRes
  geoDat = as.geodata(cbind(Xs, Ys, res), realizations=reals)
  
  kappa=3 #Matern smoothness parameter
  #estimate process variance, ignoring outliers
  #sig2 = var(pctRes[abs(pctRes) < quantile(abs(pctRes), .99)])
  sig2 = var(pctRes, na.rm=TRUE)
  #phi=10000
  #phis = c((max(Y) - min(Y))/20, (max(Y) - min(Y))/10, (max(Y) - min(Y))/7)
  phis = c((max(Y) - min(Y))/3)
  ini.cov.pars = cbind(rep(sig2, length(phis)), phis)
  limits = pars.limits(phi=c(lower=0, upper=(max(Y) - min(Y))))
  out = likfit(geoDat, kappa=kappa, ini.cov.pars=ini.cov.pars, psiA=pi/2, fix.psiR=TRUE, 
               realisations=reals, limits=limits, lik.method="REML")
  
  print(paste0("Final phi is: ", out$phi))
  print(paste0("Final psiR is: ", out$aniso.pars[2]))
  print(paste0("Final phi*psiR is: ", out$phi * out$aniso.pars[2]))
  print(paste0("Final sigma2 is: ", out$sigmasq))
  print(paste0("Final sigma is: ", sqrt(out$sigmasq)))
  
  return(out)
}

##### goal: make predictions
# N: number of data subsamples for predictor training
# numPredictions
krigeModel = function(seaDef, model, N = 100, numPredictions = 2, predHull=NULL) {
  ##### load data:
  # load seafloor deformations, residuals from regression, coord distance information
  resids = getPctResDat(seaDef)
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  files = seaDef$files
  defDat = seaDef$dat
  
  SS = resids$SS
  pctRes = resids$pctRes
  mask = resids$mask
  Mw = resids$Mw
  type = resids$type
  probs = resids$probs
  outThresh = quantile(abs(pctRes), .99, na.rm=TRUE)
  
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
  
  #plot prediction points
  plotSubPoly(rbind(hullPts, hullPts[1,]), maskXY)
  
  ##### subsample data for training for prediction
  gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
  maskGridXYI = gridXYI[mask,]
  maskSampleI = sample(1:nrow(maskGridXYI), N)
  sampleXYI = maskGridXYI[maskSampleI,]
  sampleRealI = sample(1:dim(pctRes)[3], N, replace=TRUE)
  #sampleRealI = rep(15, N)
  resSelect = cbind(sampleXYI, sampleRealI)
  
  sampleX = X[sampleXYI[,1]]
  sampleY = Y[sampleXYI[,2]]
  sampleRes = pctRes[resSelect]
  
  #remove outliers
  outI = abs(sampleRes) > outThresh
  sampleX = sampleX[!outI]
  sampleY = sampleY[!outI]
  sampleRes = sampleRes[!outI]
  
  sampleGrid = list(x=sampleX, y=sampleY)
  quilt.plot(sampleX, sampleY, sampleRes)
  
  #make geodata object with subsampled data
  reals = sampleRealI
  Xs = sampleX
  Ys = sampleY
  res = sampleRes
  geoDat = as.geodata(cbind(Xs, Ys, res), realizations=reals)
  
  ##### make predictions
  kappa=model$kappa #Matern smoothness parameter
  aniso.pars = model$aniso.pars
  lambda = model$lambda
  mod.cont = model.control(kappa=kappa, aniso.pars=aniso.pars, lambda=lambda)
  nu = model$cov.pars[2]
  phi.discrete = seq(min(dx, dy), model$phi*2, l=51)
  probs = dexp(phi.discrete, rate=1/nu)
  probs=probs/sum(probs)
  # prior.cont = prior.control(phi.prior="exponential", phi=model$cov.pars[2])
  # prior.cont = prior.control(phi.prior="reciprocal")
  # prior.cont = prior.control(phi.prior=probs, phi.discrete=phi.discrete)
  prior.cont = prior.control(phi.prior="exponential", phi=nu, phi.discrete=phi.discrete)
  # prior.cont = prior.control(phi.prior="rec", phi.discrete=phi.discrete)
  out.cont  = output.control(simulations=TRUE, n.post=11)
  krigeMod1 = krige.bayes(geoDat, locations=gridNewXY, borders=rbind(hullPts,hullPts[1,]), 
                         model=mod.cont, prior=prior.cont, output=out.cont)
  
  kc = krige.control(obj.m=model)
  oc = output.control(simulations=TRUE, n.pred=numPredictions)
  krigeMod2 = krige.conv(geoDat, locations=gridNewXY, borders=rbind(hullPts,hullPts[1,]), 
                         output = oc, krige=kc)
  
  par(mfrow=c(2,1), mar=c(3,3,0.5, 0.5))
  i=4
  image(krigeMod2, val=krigeMod2$simulation[,i], zlim=c(-.3, .3))
  legend.krige(x.leg=c(22400000, 22600000), y.leg=c(4600000, 5200000), 
               values=krigeMod2$simulation[,i], vertical=TRUE, offset.leg=-2)
}

makePredHull = function(seaDef, alpha=4000) {
  ##### load data:
  # load seafloor deformations, residuals from regression, coord distance information
  resids = getResDat(seaDef)
  mask = resids$mask
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  X = coordDat$X
  Y = coordDat$Y
  
  ##### make prediction points (approximately equidistant)
  #make prediction grid
  gridXY = make.surface.grid(list(x=X, y=Y))
  maskXY = gridXY[mask,]
  gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewXY = make.surface.grid(list(x=gridNewX, y=gridNewY))
  
  #make concave hull around prediction mask to get final prediction points
  ashape = ahull(jitter(maskXY[,1], factor=.025), jitter(maskXY[,2], factor=.025), alpha=alpha)
  
  return(ashape)
}

plotPredHull = function(seaDef, ashape, toMeters=TRUE) {
  resids = getPctResDat(seaDef)
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  
  mask = resids$mask
  
  if(toMeters) {
  X = coordDat$X
  Y = coordDat$Y
  dx = coordDat$distPerCellX
  dy = coordDat$distPerCellY
  }
  else {
    X = lon
    Y = lat
    dx = lon[2] - lon[1]
    dy = lat[2] - lat[1]
  }
  
  ##### make prediction points (approximately equidistant)
  #make prediction grid
  gridXY = make.surface.grid(list(x=X, y=Y))
  maskXY = gridXY[mask,]
  #   gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  #   gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewX = seq(min(maskXY[,1]), max(maskXY[,1]), l=30)
  gridNewY = seq(min(maskXY[,2]), max(maskXY[,2]), l=100)
  gridNewXY = make.surface.grid(list(x=gridNewX, y=gridNewY))
  
  indx=ashape$arcs[,"end1"]  
  hullPts <- maskXY[indx,]                  # extract the boundary points from maskXY
  
  #plot prediction points
  plotSubPoly(rbind(hullPts, hullPts[1,]), maskXY)
  
  invisible(NULL)
}

plotSubPoly = function(bds, maskXY, range=NULL) {
  if(! is.null(range))
    bds = rbind(bds[range,], bds[range[1],])
  quilt.plot(maskXY, rep(1, nrow(maskXY)))
  polygon(x=bds[,1], y=bds[,2], col=rgb(.5, .1, .1, .5))
}

#N is num samples
plotSurfXY = function(seaDef, N = 3000) {
  ##### load and prepare data:
  # load seafloor deformations, residuals from regression, coord distance information
  resids = getResDat(seaDef)
  
  R = 3959 # in miles
  R = 6.371*10^6 #in meters
  coordDat = calcLonLatDist(seaDef$lon, seaDef$lat, R = R)
  
  lon = seaDef$lon
  lat = seaDef$lat
  files = seaDef$files
  defDat = seaDef$dat
  
  SS = resids$SS
  resDat = resids$resids
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
  gridXY = make.surface.grid(list(x=X, y=Y))
  
  #calculate percent error data
  pctRes = resDat/defDat
  pctRes[defDat == 0] = 0 #both resDat and defDat are 0, so set error to
  
  #mask data and grid
  maskXY = gridXY[mask,]
  
  ##### subsample data for training for prediction
  gridXYI = make.surface.grid(list(xi=1:length(X), yi=1:length(Y)))
  maskGridXYI = gridXYI[mask,]
  maskSampleI = sample(1:nrow(maskGridXYI), N)
  sampleXYI = maskGridXYI[maskSampleI,]
  sampleRealI = sample(1:dim(resDat)[3], N, replace=TRUE)
  #sampleRealI = rep(15, N)
  resSelect = cbind(sampleXYI, sampleRealI)
  
  sampleX = X[sampleXYI[,1]]
  sampleY = Y[sampleXYI[,2]]
  sampleRes = pctRes[resSelect]
  
  for(r in min(sampleRealI):max(sampleRealI)) {
    # get data from realization r
    sampleRI = sampleRealI == r
    sampleXR = sampleX[sampleRI]
    sampleYR = sampleY[sampleRI]
    sampleResR = sampleRes[sampleRI]
    #plot full surface
    sampleGrid = list(x=sampleXR, y=sampleYR)
    quilt.plot(sampleXR, sampleYR, sampleResR, main=paste0("2D Surface Realization ", r))
    
    # plot x and y components
    plot(sampleXR, sampleResR, main=paste0("X Surface Realization ", r), pch=19, cex=.5)
    plot(sampleYR, sampleResR, main=paste0("Y Surface Realization ", r), pch=19, cex=.5)
  }
  
}
