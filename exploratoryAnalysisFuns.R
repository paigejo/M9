# miscellaneous functions used for exploratory analysis

# makes plots comparing the predicted subsidence levels to the subsidence data.
# if prediction inputs are left NULL, they are computed using params
comparePredsToSubs = function(params, slipPreds=NULL, slipPredsGPS=NULL, subPreds=NULL, 
                              subPredsGPS=NULL, nsim=100, plotNameRoot="full", 
                              savePlots=TRUE, G=NULL, fileNameRoot=plotNameRoot, 
                              muVec=NULL, useGPS=FALSE, tvec=NULL) {
  # get parameters
  if(is.null(muVec)) {
    lambdaMLE = params[1]
    muZetaMLE = params[2]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
    muZetaCSZ = rep(muZetaMLE, nrow(csz))
  }
  else {
    lambdaMLE = params[1]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaMLE = muVec
    muZetaGPS = muVec[1:nrow(slipDatCSZ)]
    muZetaCSZ = muVec[(nrow(slipDatCSZ)+1):length(muVec)]
  }
  
  # get taper if necessary
  if(is.null(tvec))
    tvec = taper(csz$depth, lambda=lambdaMLE)
  
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = preds(params, nsim=nsim, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec)
  if(is.null(slipPredsGPS))
    slipPredsGPS = predsGivenGPS(params, nsim=nsim, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec)
  if(is.null(subPreds)) {
    if(is.null(G))
      subPreds = predsToSubsidence(params, slipPreds, useMVNApprox = FALSE)
    else
      subPreds = predsToSubsidence(params, slipPreds, G=G, useMVNApprox = FALSE)
  }
  if(is.null(subPredsGPS)) {
    if(is.null(G))
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, useMVNApprox = FALSE)
    else
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, G=G, useMVNApprox = FALSE)
  }
  meanSlip = slipPreds$meanSlip
  meanSub = subPreds$meanSub
  u95 = subPreds$u95
  l95 = subPreds$l95
  u95Noise = subPreds$u95Noise
  l95Noise = subPreds$l95Noise
  meanSlipGPS = slipPredsGPS$meanSlip
  meanSubGPS = subPredsGPS$meanSub
  u95GPS = subPredsGPS$u95
  l95GPS = subPredsGPS$l95
  u95NoiseGPS = subPredsGPS$u95Noise
  l95NoiseGPS = subPredsGPS$l95Noise
  
  ##### generate mean seaDef field from Okada model
  # set Okada subsidence grid
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  coordGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
  meanSeaDef = okada(csz, lonGrid, latGrid, slips=meanSlip)
  
  ##### plot subsidence predictions
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubsidencePredictions.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  plotFault(csz, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
            xlim=lonRange)
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(meanSeaDef)), nx=120, ny=120, 
           main=paste0(plotNameRoot, "Mean Uplift (m)"), xlab="Longitude")
  map("world", "Canada", add=TRUE, lwd=1)
  US(add=TRUE, lwd=1)
  points(dr1$Lon, dr1$Lat, cex=1, pch=3, col="red")
  plotFault(csz, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, dr1$subsidence[inds]))
  # simulated subsidence data from Okada model using marginal distribution
  plot(dr1$subsidence, dr1$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, "95% Confidence Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, dr1$Lat, pch=19, cex=.3, col="blue")
  ord = order(dr1$Lat)
  lines(-l95[ord], dr1$Lat[ord], col="blue")
  lines(-u95[ord], dr1$Lat[ord], col="blue")
  # simulated subsidence data from Okada model using marginal distribution with observation noise
  plot(dr1$subsidence, dr1$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, "95% Prediction Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, dr1$Lat, pch=19, cex=.3, col="blue")
  ord = order(dr1$Lat)
  lines(-l95Noise[ord], dr1$Lat[ord])
  lines(-u95Noise[ord], dr1$Lat[ord])
  # simulated subsidence data from Okada model using GPS locking rates
#   plot(dr1$subsidence, dr1$Lat, pch="+", col="red", ylim=latRange, 
#        xlim=subRange, main=paste0(fileNameRoot, " Simulated Subsidence Given GPS Data"), ylab="", 
#        xlab="Subsidence (m)")
#   points(-meanSubGPS, dr1$Lat, pch=19, cex=.3, col="blue")
#   ord = order(dr1$Lat)
#   lines(-l95GPS[ord], dr1$Lat[ord], col="blue")
#   lines(-u95GPS[ord], dr1$Lat[ord], col="blue")
#   lines(-l95NoiseGPS[ord], dr1$Lat[ord])
#   lines(-u95NoiseGPS[ord], dr1$Lat[ord])
  if(savePlots)
    dev.off()
  
  ##### plot GPS predictions
  # first generate the GPS predictions
#   mucGPS = slipPreds$mucGPS
#   if(!is.null(slipPredsGPS$SigmacDiagGPS))
#     SigmaDiagGPS = slipPredsGPS$SigmacDiagGPS
#   else
#     SigmaDiagGPS = diag(slipPredsGPS$SigmacGPS)
#   sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
#   meanPredsGPS = exp(mucGPS + muXi + SigmaDiagGPS/2 + sigmaXi^2/2)
#   
#   # Now make the plot
#   if(savePlots)
#     pdf(file=paste0(fileNameRoot, "GPSPredictions.pdf"), width=7, height=5)
#   par(mfrow=c(1,2))
#   # plot mean slip field times Xi (mean locking rate)
#   datRange = range(c(meanPredsGPS, slipDatCSZ$slip))
#   quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, meanPredsGPS, 
#              main=paste0(plotNameRoot, " Predicted Mean Locking Rate (mm/yr)"), 
#              xlim=lonRange, ylim=latRange, zlim=datRange)
#   plotFault(csz, new=FALSE, plotData = FALSE)
#   map("world", "Canada", lwd=1, add=TRUE)
#   US(add=TRUE, lwd=1)
#   # plot GPS data
#   quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="Observed Locking Rate (mm/yr)", 
#              xlim=lonRange, ylim=latRange, zlim=datRange)
#   plotFault(csz, new=FALSE, plotData = FALSE)
#   map("world", "Canada", lwd=1, add=TRUE)
#   US(add=TRUE, lwd=1)
#   if(savePlots)
#     dev.off()
  
  invisible(NULL)
}

comparePredsToGPS = function(params, muVec=NULL) {
  # get data
  logX = log(slipDatCSZ$slip)
  sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  if(is.null(muVec))
    muVec = rep(muZeta, length(logX))
  else
    muVec = muVec[1:nrow(slipDatCSZ)]
  
  par(mfrow=c(1,2))
  hist(logX, freq=FALSE, main="log GPS data with Gaussian fit", xlab="Log Locking Rate")
  xs = seq(0, 5, l=200)
  lines(xs, dnorm(xs, mean(logX, sd(logX))), col="blue")
  hist(rnorm(1000*length(sigmaXi), muXi+muVec, sqrt(sigmaXi^2+sigmaZeta^2)), xlim=range(logX), 
       breaks=500, freq=F, ylim=c(0, 1.2), main="Fit log GPS data distribution", 
       xlab="Predicted Log Locking Rate")
  lines(xs, dnorm(xs, mean(logX, sd(logX))), col="blue")
}

getAvgSlipFromMoment = function(mag=9.0, rigidity=4*10^10, lambda=NULL) {
  # compute taper values
  if(is.null(lambda))
    lambda = fixedFitMVN$lambdaMLE
  tvec = taper(csz$depth, lambda=lambda)
  
  # get fault total area
  areas = csz$length*csz$width
  totArea = sum(areas)
  totWeightedArea = sum(areas/tvec)
  
  Mo_desired = 10^(1.5*mag + 9.05)
  
  meanSlipPreTaper = Mo_desired / (rigidity * totArea)
  meanSlipPostTaper = Mo_desired / (rigidity * totWeightedArea)
  
  return(c(pre=meanSlipPreTaper, post=meanSlipPostTaper))
}

getMomentFromAvgSlip = function(slip=3, rigidity=4*10^10, lambda=NULL) {
  # compute taper values
  if(is.null(lambda))
    lambda = fixedFitMVN$lambdaMLE
  tvec = taper(csz$depth, lambda=lambda)
  
  # get fault total area
  areas = csz$length*csz$width
  totArea = sum(areas)
  totWeightedArea = sum(areas/tvec)
  
  MoNoTaper = slip*rigidity*totArea
  MoTaper = slip*rigidity*totWeightedArea
  
  magNoTaper = (log10(MoNoTaper) - 9.05)/1.5
  magTaper = (log10(MoTaper) - 9.05)/1.5
  
  return(c(taper=magTaper, noTaper=magNoTaper))
}

getMomentFromSlip = function(slips, rigidity=4*10^10, doTaper=FALSE, lambda=1, depths=csz$depth) {
  # get fault total area
  areas = csz$length*csz$width
  totArea = sum(areas)
  
  if(doTaper)
    slips = slips*taper(depths, lambda)
  
  Mo = sum(slips*rigidity*areas)
  
  (log10(Mo) - 9.05)/1.5
}

#based on http://onlinelibrary.wiley.com/doi/10.1002/9780470824566.app1/pdf and 
#http://blogs.sas.com/content/iml/2010/12/10/converting-between-correlation-and-covariance-matrices.html
hessianToCorrMat = function(hessMat) {
  covMat = solve(-hessMat)
  D = diag(1/sqrt(diag(covMat)))
  corrMat = D %*% covMat %*% D
  return(corrMat)
}
hessianToCovMat = function(hessMat) {
  solve(-hessMat)
}