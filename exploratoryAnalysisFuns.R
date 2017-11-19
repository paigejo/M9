# miscellaneous functions used for exploratory analysis

# makes plots comparing the predicted subsidence levels to the subsidence data.
# if prediction inputs are left NULL, they are computed using params
comparePredsToSubs = function(params, slipPreds=NULL, slipPredsGPS=NULL, subPreds=NULL, 
                              subPredsGPS=NULL, nsim=100, plotNameRoot="full", 
                              savePlots=TRUE, G=NULL, fileNameRoot=plotNameRoot, 
                              muVec=NULL, useGPS=FALSE, tvec=NULL, subDat=dr1, 
                              logScale=FALSE, fault=csz, latRange=c(40, 50), 
                              posNormalModel=FALSE, normalModel=posNormalModel, doGPSPred=FALSE, 
                              useMVNApprox=FALSE, taperedGPSDat=FALSE, dStar=28000, normalizeTaper=FALSE) {
  # get parameters
  if(is.null(muVec)) {
    lambdaMLE = params[1]
    muZetaMLE = params[2]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
    muZetaCSZ = rep(muZetaMLE, nrow(fault))
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
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaMLE, dStar=dStar, normalize=normalizeTaper)
  
  #
  if(taperedGPSDat)
    phiZeta = params[length(params)]
  else
    phiZeta = NULL
  
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = preds(params, nsim=nsim, fault=fault, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec, 
                      posNormalModel=posNormalModel, normalModel=normalModel, phiZeta=phiZeta)
  if(is.null(slipPredsGPS) && doGPSPred)
    slipPredsGPS = predsGivenGPS(params, nsim=nsim, muVec=c(muZetaGPS, muZetaCSZ), fault=fault, tvec=tvec, 
                                 posNormalModel=posNormalModel, normalModel=normalModel)
  if(is.null(subPreds)) {
    if(is.null(G))
      subPreds = predsToSubsidence(params, slipPreds, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPreds = predsToSubsidence(params, slipPreds, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  if(is.null(subPredsGPS) && doGPSPred) {
    if(is.null(G))
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  meanSlip = slipPreds$meanSlip
  meanSub = subPreds$meanSub
  u95 = subPreds$u95
  l95 = subPreds$l95
  u95Noise = subPreds$u95Noise
  l95Noise = subPreds$l95Noise
  slipSD = apply(slipPreds$slipSims, 1, sd)
  if(doGPSPred) {
    meanSlipGPS = slipPredsGPS$meanSlip
    meanSubGPS = subPredsGPS$meanSub
    u95GPS = subPredsGPS$u95
    l95GPS = subPredsGPS$l95
    u95NoiseGPS = subPredsGPS$u95Noise
    l95NoiseGPS = subPredsGPS$l95Noise
  }
  
  ##### generate mean seaDef field from Okada model
  # set Okada subsidence grid
  lonRange=c(-128, -122.5)
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  coordGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
  meanSeaDef = okada(fault, lonGrid, latGrid, slips=meanSlip)
  
  ##### plot subsidence predictions
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubsidencePredictions.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  if(!logScale) {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, " Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  }
  else {
    plotFault(fault, plotVar=meanSlip, legend.mar=6, main=paste0(plotNameRoot, " Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, varRange=c(.1, max(meanSlip)), logScale=TRUE)
  }
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(meanSeaDef)), nx=120, ny=120, 
           main=paste0(plotNameRoot, " Mean Uplift (m)"), xlab="Longitude")
  map("world", "Canada", add=TRUE, lwd=1)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  plotFault(fault, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, subDat$subsidence))
  # simulated subsidence data from Okada model using marginal distribution
  plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, " 95% Confidence Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, subDat$Lat, pch=19, cex=.3, col="blue")
  ord = order(subDat$Lat)
  lines(-l95[ord], subDat$Lat[ord], col="blue")
  lines(-u95[ord], subDat$Lat[ord], col="blue")
  # simulated subsidence data from Okada model using marginal distribution with observation noise
  plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, " 95% Prediction Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, subDat$Lat, pch=19, cex=.3, col="blue")
  ord = order(subDat$Lat)
  lines(-l95Noise[ord], subDat$Lat[ord])
  lines(-u95Noise[ord], subDat$Lat[ord])
  if(savePlots)
    dev.off()
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubsidencePredictions3.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  if(!logScale) {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  }
  else {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, varRange=c(.1, max(meanSlip)), logScale=TRUE)
  }
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(meanSeaDef)), nx=120, ny=120, 
           main=paste0(plotNameRoot, " Mean Uplift (m)"), xlab="Longitude")
  map("world", "Canada", add=TRUE, lwd=1)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  plotFault(csz, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, subDat$subsidence))
  # simulated subsidence data from Okada model using marginal distribution
  plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, "95% Confidence Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, subDat$Lat, pch=19, cex=.3, col="blue")
  ord = order(subDat$Lat)
  lines(-l95[ord], subDat$Lat[ord], col="blue")
  lines(-u95[ord], subDat$Lat[ord], col="blue")
  # plot normalized residuals distribution
  resids = subDat$subsidence/(-meanSub) - 1
  hist(resids, breaks=50, main="Histogram of normalized residuals", xlab="Normalized Residuals", 
       freq=F)
  xs = seq(0, max(resids), l=100)-1
  cleanResids = resids[is.finite(resids)]
  lines(xs, dlnorm(xs+1, 0, var(log(cleanResids+1), na.rm = TRUE)), col="blue")
  
  if(savePlots)
    dev.off()
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubsidencePredictions2.pdf"), width=10, height=5)
  par(mfrow=c(1,3))
  if(!logScale) {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  }
  else {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, varRange=c(.1, max(meanSlip)), logScale=TRUE)
  }
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(meanSeaDef)), nx=120, ny=120, 
           main=paste0(plotNameRoot, " Mean Uplift (m)"), xlab="Longitude")
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  plotFault(fault, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, subDat$subsidence))
  # simulated subsidence data from Okada model using marginal distribution
  plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, "95% Confidence Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, subDat$Lat, pch=19, cex=.5, col="blue")
  ord = order(subDat$Lat)
  lines(-l95[ord], subDat$Lat[ord], col="blue")
  lines(-u95[ord], subDat$Lat[ord], col="blue")
  if(savePlots)
    dev.off()
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubsidencePredictions4.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  if(!logScale) {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  }
  else {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, varRange=c(.1, max(meanSlip)), logScale=TRUE)
  }
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(meanSeaDef)), nx=120, ny=120, 
           main=paste0(plotNameRoot, " Mean Uplift (m)"), xlab="Longitude")
  map("world", "Canada", add=TRUE, lwd=1)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  plotFault(csz, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, subDat$subsidence))
  # simulated subsidence data from Okada model using marginal distribution
  plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, "95% Confidence Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, subDat$Lat, pch=19, cex=.3, col="blue")
  ord = order(subDat$Lat)
  lines(-l95[ord], subDat$Lat[ord], col="blue")
  lines(-u95[ord], subDat$Lat[ord], col="blue")
  # plot magnitude distribution
  mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault, dStar=dStar, normalizeTaper=normalizeTaper)
  cleanMags = mags[is.finite(mags)]
  hist(cleanMags, breaks=50, main="Histogram of earthquake magnitudes", xlab="Magnitude", 
       freq=F)
  xs = seq(min(cleanMags), max(cleanMags), l=100)
  lines(xs, dnorm(xs, mean(cleanMags), sd(cleanMags)), col="blue")
  abline(v=mean(cleanMags), col="black")
  abline(v=quantile(cleanMags, probs=.975), col="purple", lty=2)
  abline(v=quantile(cleanMags, probs=.025), col="purple", lty=2)
  
  if(savePlots)
    dev.off()
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubsidencePredictionsFinal.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  if(!logScale) {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  }
  else {
    plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, varRange=c(.1, max(meanSlip)), logScale=TRUE)
  }
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotFault(fault, plotVar = slipSD, legend.mar=6, main=paste0(plotNameRoot, "Slip SD (m)"), 
            xlim=lonRange, ylim=latRange)
  map("world", "Canada", add=TRUE, lwd=1)
  US(add=TRUE, lwd=1)
  points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
  plotFault(csz, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, subDat$subsidence))
  # simulated subsidence data from Okada model using marginal distribution
  plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
       xlim=subRange, main=paste0(plotNameRoot, "95% Prediction Band"), ylab="", 
       xlab="Subsidence (m)")
  points(-meanSub, subDat$Lat, pch=19, cex=.3, col="blue")
  ord = order(subDat$Lat)
  lines(-l95Noise[ord], subDat$Lat[ord], col="blue")
  lines(-u95Noise[ord], subDat$Lat[ord], col="blue")
  # plot magnitude distribution
  mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault, dStar=dStar, normalizeTaper=normalizeTaper)
  cleanMags = mags[is.finite(mags)]
  hist(cleanMags, breaks=50, main="Histogram of earthquake magnitudes", xlab="Magnitude", 
       freq=F)
  xs = seq(min(cleanMags), max(cleanMags), l=100)
  lines(xs, dnorm(xs, mean(cleanMags), sd(cleanMags)), col="blue")
  abline(v=mean(cleanMags), col="black")
  abline(v=quantile(cleanMags, probs=.975), col="purple", lty=2)
  abline(v=quantile(cleanMags, probs=.025), col="purple", lty=2)
  
  if(savePlots)
    dev.off()
  
  ##### plot 6 slip simulations
  if(ncol(slipPreds$slipSims) >= 6) {
    if(savePlots)
      pdf(file=paste0(fileNameRoot, "SlipSimulations.pdf"), width=10, height=10)
    par(mfrow=c(2,3))
    for(i in 1:6) {
      if(!logScale) {
        plotFault(fault, plotVar = slipPreds$slipSims[,i], legend.mar=6, main=paste0(plotNameRoot, "Coseismic Slips ", i, " (m)"), 
                  xlim=lonRange, ylim=latRange)
      }
      else {
        plotFault(fault, plotVar = slipPreds$slipSims[,i], legend.mar=6, main=paste0(plotNameRoot, "Coseismic Slips ", i, " (m)"), 
                  xlim=lonRange, ylim=latRange, varRange=c(.1, max(slipPreds$slipSims[,i])), logScale=TRUE)
      }
      map("world", "Canada", lwd=1, add=TRUE)
      US(add=TRUE, lwd=1)
      points(subDat$Lon, subDat$Lat, cex=1, pch=3, col="red")
    }
    
    if(savePlots)
      dev.off()
  }
  
  ##### plot 6 subsidence simulations
  if(ncol(subPreds$subSims) >= 6) {
    if(savePlots)
      pdf(file=paste0(fileNameRoot, "SubsidenceSimulations.pdf"), width=8, height=10)
    par(mfrow=c(2,3))
    for(i in 1:6) {
      plot(subDat$subsidence, subDat$Lat, pch="+", col="red", ylim=latRange, 
           xlim=subRange, main=paste0(plotNameRoot, "Coseismic Subsidences ", i), ylab="", 
           xlab="Subsidence (m)")
      points(-subPreds$subSims[,i], subDat$Lat, pch=19, cex=.5, col="blue")
    }
    
    if(savePlots)
      dev.off()
  }
  
  ##### plot slip quantiles
  slipSims = slipPreds$slipSims
  u95 = apply(slipSims, 1, quantile, probs=.975)
  l95 = apply(slipSims, 1, quantile, probs=.025)
  medSlip = apply(slipSims, 1, median)
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "slipDistn.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  # mean
  plotFault(fault, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
            xlim=lonRange, ylim=latRange, logScale = logScale, varRange=c(.1, max(meanSlip)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  # median
  plotFault(fault, plotVar = medSlip, legend.mar=6, main=paste0(plotNameRoot, "Median Slip (m)"), 
            xlim=lonRange, ylim=latRange, ylab="", logScale = logScale, varRange=c(.1, max(medSlip)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  # 2.5th percentile
  plotFault(fault, plotVar = l95, legend.mar=6, main=paste0(plotNameRoot, "2.5th Percentile Slip (m)"), 
            xlim=lonRange, ylim=latRange, logScale = logScale, varRange=c(min(.01, c(l95)), max(l95)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  ## 97.5th percentile
  plotFault(fault, plotVar = u95, legend.mar=6, main=paste0(plotNameRoot, "97.5th Percentile Slip (m)"), 
            xlim=lonRange, ylim=latRange, ylab="", logScale = logScale, varRange=c(.1, max(u95)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  
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

# plot 2x2 grid of plots with mean, 50th percentile, 2.5th percentile, and 97.5th percentile slips
plotSlipDistribution = function(params, slipPreds=NULL, slipPredsGPS=NULL, subPreds=NULL, 
                                subPredsGPS=NULL, nsim=100, plotNameRoot="full", 
                                savePlots=TRUE, fileNameRoot=plotNameRoot, 
                                muVec=NULL, useGPS=FALSE, tvec=NULL, logScale=FALSE, taperedSlipDat=FALSE, 
                                fault=csz, normalizeTaper=FALSE, dStar=28000) {
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
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaMLE, normalize=normalizeTaper, dStar=dStar)
  
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = preds(params, nsim=nsim, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec, taperedSlipDat=taperedSlipDat)
  meanSlip = slipPreds$meanSlip
  
  slipSims = slipPreds$slipSims
  u95 = apply(slipSims, 1, quantile, probs=.975)
  l95 = apply(slipSims, 1, quantile, probs=.025)
  medSlip = apply(slipSims, 1, median)
  
  ##### First plot areal slips
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "slipDistn.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  # mean
  plotFault(csz, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, "Mean Slip (m)"), 
            xlim=lonRange, ylim=latRange, logScale = logScale, varRange=c(.1, max(meanSlip)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  # median
  plotFault(csz, plotVar = medSlip, legend.mar=6, main=paste0(plotNameRoot, "Median Slip (m)"), 
  xlim=lonRange, ylim=latRange, ylab="", logScale = logScale, varRange=c(.1, max(medSlip)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  # 2.5th percentile
  plotFault(csz, plotVar = l95, legend.mar=6, main=paste0(plotNameRoot, "2.5th Percentile Slip (m)"), 
            xlim=lonRange, ylim=latRange, logScale = logScale, varRange=c(.01, max(l95)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  ## 97.5th percentile
  plotFault(csz, plotVar = u95, legend.mar=6, main=paste0(plotNameRoot, "97.5th Percentile Slip (m)"), 
            xlim=lonRange, ylim=latRange, ylab="", logScale = logScale, varRange=c(.1, max(u95)))
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  
  if(savePlots)
    dev.off()
}

plotFixedSlip = function(meanSlip, medSlip=NULL, l95, u95, slipSD=NULL, plotNameRoot="full", 
                         savePlots=TRUE, fileNameRoot=plotNameRoot, logScale=FALSE, 
                         event="All", maxVal=Inf) {
  if(event == "All")
    inds = 1:nrow(dr1)
  else
    inds = events == event
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "slipDistn.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  # mean
  if(logScale)
    plotFault(csz, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, " Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, logScale = logScale, varRange=c(.1, min(max(meanSlip), maxVal)))
  else
    plotFault(csz, plotVar = meanSlip, legend.mar=6, main=paste0(plotNameRoot, " Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  # median/sd
  if(is.null(slipSD)) {
    if(logScale)
      plotFault(csz, plotVar = medSlip, legend.mar=6, main=paste0(plotNameRoot, " Median Slip (m)"), 
                xlim=lonRange, ylim=latRange, ylab="", logScale = logScale, varRange=c(.1, min(max(medSlip), maxVal)))
    else
      plotFault(csz, plotVar = medSlip, legend.mar=6, main=paste0(plotNameRoot, " Median Slip (m)"), 
                xlim=lonRange, ylim=latRange, ylab="")
    map("world", "Canada", lwd=1, add=TRUE)
    US(add=TRUE, lwd=1)
    points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  }
  else if(is.null(medSlip)) {
    if(logScale)
      plotFault(csz, plotVar = slipSD, legend.mar=6, main=paste0(plotNameRoot, " Slip SD (m)"), 
                xlim=lonRange, ylim=latRange, ylab="", logScale = logScale, varRange=c(.1, min(max(slipSD), maxVal)))
    else
      plotFault(csz, plotVar = slipSD, legend.mar=6, main=paste0(plotNameRoot, " Slip SD (m)"), 
                xlim=lonRange, ylim=latRange, ylab="")
    map("world", "Canada", lwd=1, add=TRUE)
    US(add=TRUE, lwd=1)
    points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  }
  # 2.5th percentile
  if(logScale)
    plotFault(csz, plotVar = l95, legend.mar=6, main=paste0(plotNameRoot, " 2.5th Percentile Slip (m)"), 
              xlim=lonRange, ylim=latRange, logScale = logScale, varRange=c(.01, min(max(l95), maxVal)))
  else
    plotFault(csz, plotVar = l95, legend.mar=6, main=paste0(plotNameRoot, " 2.5th Percentile Slip (m)"), 
              xlim=lonRange, ylim=latRange)
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  ## 97.5th percentile
  if(logScale)
    plotFault(csz, plotVar = u95, legend.mar=6, main=paste0(plotNameRoot, " 97.5th Percentile Slip (m)"), 
              xlim=lonRange, ylim=latRange, ylab="", logScale = logScale, varRange=c(.1, min(max(u95), maxVal)))
  else
    plotFault(csz, plotVar = u95, legend.mar=6, main=paste0(plotNameRoot, " 97.5th Percentile Slip (m)"), 
              xlim=lonRange, ylim=latRange, ylab="")
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  
  if(savePlots)
    dev.off()
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

plotSplineUncertainty = function(splinePar, covMat, nKnots=5, niter=1000, latRange=c(40,50), 
                                 diffGPSTaper=FALSE, nKnotsGPS=5, latsOnX=TRUE, 
                                 main=TeX("$\\lambda$ 95% Confidence Band")) {
  L = t(chol(covMat))
  lats = seq(latRange[1], latRange[2], l=100)
  Xi = getSplineBasis(csz, latRange, nKnots, lats)
  if(diffGPSTaper) {
    XiGPS = getSplineBasis(csz, latRange, nKnotsGPS, lats)
    Xi = cbind(Xi, -XiGPS)
  }
  
  Zs = matrix(rnorm(ncol(L)*niter), nrow=ncol(L), ncol=niter)
  sims = Xi %*% L %*% Zs
  lows = apply(sims, 1, quantile, probs=0.025)
  his = apply(sims, 1, quantile, probs=0.975)
  cntr = Xi %*% splinePar
  lows = cntr + lows
  his = cntr + his
  
  if(latsOnX) {
    yRange=c(min(lows), max(his))
    plot(lats, cntr, main=main, ylim=yRange, type="l", 
         col="blue", ylab=TeX("$\\lambda$"), xlab="Latitude")
    lines(lats, lows, col="black")
    lines(lats, his, col="black")
  }
  else {
    xRange=c(min(lows), max(his))
    plot(cntr, lats, main=main, xlim=xRange, type="l", 
         col="blue", xlab=TeX("$\\lambda$"), ylab="Latitude")
    lines(lows, lats, col="black")
    lines(his, lats, col="black")
  }
  
  invisible(list(lows=lows, his=his, cntr=cntr))
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions for computing earhtquake magnitude

getAvgSlipFromMoment = function(mag=9.0, rigidity=4*10^10, lambda=NULL, dStar=28000, normalizeTaper=FALSE) {
  # compute taper values
  if(is.null(lambda))
    lambda = fixedFitMVN$lambdaMLE
  tvec = taper(getFaultCenters(csz)[,3], lambda=lambda, normalize=normalizeTaper, dStar=dStar)
  
  # get fault total area
  areas = csz$length*csz$width
  totArea = sum(areas)
  totWeightedArea = sum(areas/tvec)
  
  Mo_desired = 10^(1.5*mag + 9.05)
  
  meanSlipPreTaper = Mo_desired / (rigidity * totArea)
  meanSlipPostTaper = Mo_desired / (rigidity * totWeightedArea)
  
  return(c(pre=meanSlipPreTaper, post=meanSlipPostTaper))
}

getMomentFromAvgSlip = function(slip=3, rigidity=4*10^10, lambda=NULL, normalizeTaper=FALSE, dStar=28000) {
  # compute taper values
  if(is.null(lambda))
    lambda = fixedFitMVN$lambdaMLE
  tvec = taper(getFaultCenters(csz)[,3], lambda=lambda, dStar=dStar, normalize=normalizeTaper)
  
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

getMomentFromSlip = function(slips, rigidity=4*10^10, doTaper=FALSE, lambda=1, fault=csz, depths=getFaultCenters(fault)[,3], 
                             normalizeTaper=FALSE, dStar=28000) {
  # get fault total area
  areas = fault$length*fault$width
  totArea = sum(areas)
  
  if(doTaper)
    slips = slips*taper(depths, lambda, dStar=dStar, normalize=normalizeTaper)
  
  Mo = sum(slips*rigidity*areas)
  
  if(Mo <= 0) {
    warning("negative seismic moment observed")
    return(0)
  }
  
  (log10(Mo) - 9.05)/1.5
}

############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
############################################################################
# functions for computing asymptotic correlation and covariance matrices
# using hessian of log-likleihood

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

##### for variogram fitting
# make 95% envelope using this correlation model.  Takes a couple minutes, but it's 
# recommended to do ~500 simulations
plotCorrelationEnvelope = function(nsim=101, vgtype="correlogram") {
  # first fit the mKrig object
  print("Fitting covariance...")
  out = fitGPSCovariance(doLog=FALSE)
  MLE = out$mKrigObj
  betas = MLE$d
  Xmat = MLE$Tmatrix
  slipDatCntr = slipDatCSZ$slip - Xmat %*% betas
  
  # compute distances
  print("computing distances, correlations, and simulations...")
  locs = cbind(slipDatCSZ$lon, slipDatCSZ$lat)
  distMat = rdist.earth(locs, miles=FALSE)
  n <- nrow(locs)
  is = rep(1:n, n)
  js = rep(1:n, rep(n, n))
  ind <- is > js
  id <- cbind(is, js)[ind, ]
  distVec <- rdist.earth.vec(locs[id[, 1], ], locs[id[, 2],], miles=FALSE)
  
  # get correlation parameters
  corPar = getCorPar(TRUE)
  nu = corPar$nuZeta
  phi = corPar$phiZeta
  lambda = corPar$lambda # nugget to sill ratio
  rho = out$mKrigObj$rho.MLE
  
  # compute original correlogram
  ogvg = myvgram(locs, slipDatCntr, lon.lat=TRUE, d=distVec, type=vgtype, breaks=seq(0, max(distVec), l=60), 
                 colMeans=0, sigma=sqrt(rho))
  # ogvg = vgram(locs, slipDatCntr, lon.lat=TRUE, d=distVec, type=vgtype, breaks=seq(0, max(distVec), l=60))
  cntrs = ogvg$centers
  
  # get correlation matrix (for sigma=1)
  corMat = stationary.cov(locs, Covariance = "Matern", distMat=distMat, onlyUpper = TRUE, theta=phi, 
                          smoothness=nu) * (1-lambda)
  diag(corMat) = 1
  
  # convert to covariance matrix
  covMat = corMat * rho
  
  # simulate spatial field
  L = t(chol(covMat))
  fieldSims = L %*% matrix(rnorm(nrow(locs)*nsim), ncol=nsim)
  
  # compute correlograms of the simulated fields
  print("generating variograms")
  getCors = function(sim) {
    # vg = vgram(locs, sim, lon.lat=TRUE, d=distVec, type=vgtype, breaks=seq(0, max(distVec), l=60))
    vg = myvgram(locs, sim, lon.lat=TRUE, d=distVec, type=vgtype, breaks=seq(0, max(distVec), l=60), 
                 colMeans = 0, sigma=sqrt(rho))
    vg$stats[2,]
  }
  allCors = apply(fieldSims, 2, getCors)
  
  # get simulation envelope
  print("finishing up...")
  lohi = apply(allCors, 1, quantile, probs=c(.025, .975))
  
  # compute theoretical correlogram/variogram
  xs = seq(0, max(distVec), l=100)
  cors = Matern(xs, range=phi, smoothness = 3/2) * (1 - lambda)
  covs = cors * rho
  vars = rho - covs
  if(vgtype == "correlogram")
    vgs = cors
  else
    vgs = vars
  
  ## now plot the results
  # plot original data correlogram
  yrange = range(c(lohi, ogvg$stats[2,]))
  plot(cntrs, ogvg$stats[2,], ylim=yrange, type="o", xlab="Distance (km)", 
       ylab="Variogram Estimate", main=paste0("95% Monte Carlo ", vgtype, " envelope"))
  
  # add the simulation envelope
  lines(cntrs, lohi[1,], col="blue", lty=2)
  lines(cntrs, lohi[2,], col="blue", lty=2)
  
  # add theoretical correlogram/variogram
  lines(xs, vgs, col="green")
  
  # plot original data correlogram with mean of the simulations
  yrange = range(c(lohi, ogvg$stats[2,]))
  plot(cntrs, ogvg$stats[2,], ylim=yrange, type="o", xlab="Distance (km)", 
       ylab="Variogram Estimate", main=paste0("95% Monte Carlo ", vgtype, " envelope"))
  
  # add the simulation envelope
  lines(cntrs, lohi[1,], col="blue", lty=2)
  lines(cntrs, lohi[2,], col="blue", lty=2)
  
  # add theoretical correlogram/variogram
  lines(xs, vgs, col="green")
  
  # add the mean of the simulations
  means = apply(allCors, 1, mean)
  lines(cntrs, means, col="purple")
}


myvgram = function (loc, y, id = NULL, d = NULL, lon.lat = FALSE, dmax = NULL, 
                    N = NULL, breaks = NULL, type = c("variogram", "covariogram", 
                                                      "correlogram"), 
                    colMeans=NULL, sigma=NULL) 
{
  type = match.arg(type)
  y <- cbind(y)
  if (is.null(id)) {
    n <- nrow(loc)
    is = rep(1:n, n)
    js = rep(1:n, rep(n, n))
    ind <- is > js
    id <- cbind(is, js)[ind, ]
  }
  if (is.null(d)) {
    loc <- as.matrix(loc)
    if (lon.lat) {
      d <- rdist.earth.vec(loc[id[, 1], ], loc[id[, 2], 
                                               ])
    }
    else {
      d <- rdist.vec(loc[id[, 1], ], loc[id[, 2], ])
    }
  }
  if(is.null(colMeans))
    colMeans <- apply(y, 2, mean, na.rm = TRUE)
  yCntr = sweep(y, 2, colMeans)
  if (type == "correlogram") {
    if(is.null(sigma))
      sigma = apply(yCntr, 2, sd, na.rm = TRUE)
    yCntr = sweep(yCntr, 2, (1/sigma), FUN = "*")
  }
  y1Cntr = yCntr[id[, 1], ]
  y2Cntr = yCntr[id[, 2], ]
  
  if (type == "variogram") {
    vg <- 0.5 * rowMeans(cbind((y1Cntr - y2Cntr)^2), na.rm = TRUE)
  }
  else {
    vg <- rowMeans(cbind(y1Cntr * y2Cntr), na.rm = TRUE)
  }
  call <- match.call()
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  ind <- d <= dmax & !is.na(vg)
  out <- list(d = d[ind], vgram = vg[ind], call = call, type = type)
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("vgram", class(out))
  out
}

myvgram2 = function (loc, y, id = NULL, d = NULL, lon.lat = FALSE, dmax = NULL, 
                    N = NULL, breaks = NULL, type = c("variogram", "covariogram", 
                                                      "correlogram"), 
                    colMeans=NULL, sigmas=NULL) 
{
  type = match.arg(type)
  y <- cbind(y)
  if (is.null(id)) {
    n <- nrow(loc)
    is = rep(1:n, n)
    js = rep(1:n, rep(n, n))
    ind <- is > js
    id <- cbind(is, js)[ind, ]
  }
  if (is.null(d)) {
    loc <- as.matrix(loc)
    if (lon.lat) {
      d <- rdist.earth.vec(loc[id[, 1], ], loc[id[, 2], 
                                               ])
    }
    else {
      d <- rdist.vec(loc[id[, 1], ], loc[id[, 2], ])
    }
  }
  if(is.null(colMeans))
    colMeans <- apply(y, 2, mean, na.rm = TRUE)
  yCntr = sweep(y, 2, colMeans)
  if (type == "correlogram") {
    if(is.null(sigma))
      sigmas = rep(apply(yCntr, 2, sd, na.rm = TRUE), ncol(yCntr))
    yCntr = sweep(yCntr, 2, (1/sigmas), FUN = "*")
  }
  y1Cntr = yCntr[id[, 1], ]
  y2Cntr = yCntr[id[, 2], ]
  
  if (type == "variogram") {
    vg <- 0.5 * rowMeans(cbind((y1Cntr - y2Cntr)^2), na.rm = TRUE)
  }
  else {
    vg <- rowMeans(cbind(y1Cntr * y2Cntr), na.rm = TRUE)
  }
  call <- match.call()
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  ind <- d <= dmax & !is.na(vg)
  out <- list(d = d[ind], vgram = vg[ind], call = call, type = type)
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("vgram", class(out))
  out
}

##### for plotting cdf of minimum of MVN slips
plotMinMVN = function(mu, Sigma, minSeq=seq(min(mu)-2*max(sqrt(diag(Sigma))), mean(mu), l=50), fault=csz, add0=TRUE) {
  require(mvtnorm)
  
  if(add0)
    minSeq = sort(c(0, minSeq))
  
  # make mu into vector if necessary
  if(length(mu) == 1)
    mu = rep(mu, nrow(Sigma))
  
  # compute probabilities
  getPGreaterThanT = function(t) {
    pmvnorm(upper=rep(-t, nrow(fault)), mean=-mu, sigma=Sigma)
  }
  probs = sapply(minSeq, getPGreaterThanT)
  
  # plot results
  plot(minSeq, probs, type="l", main="Probability all slips > t", xlab="t", ylab="Probability")
  
  if(add0)
    abline(v=0, col="red", lty=2)
  
  # print probability at 0
  if(add0)
    print(paste0("prob > 0: ", probs[minSeq == 0]))
}

####################################################################################################
####################################################################################################
####################################################################################################
##### ggplot plotting functions ####################################################################
####################################################################################################
####################################################################################################
####################################################################################################

# plot a bunch of slip simulations on the fault geometry over a map
plotFaultSims = function(sims, fault=csz, varRange=NULL, 
                         logScale=FALSE, xlim=c(-128, -122), ylim=c(39.5, 50.5), 
                         xlab="Longitude", ylab="Latitude", main="Slip Simulations", 
                         clab="Slip (m)", nrow=2, ncol=8) {
  
  fullDF = do.call("rbind", replicate(ncol(sims), fault, simplify = FALSE))
  fullDF$sim = factor(rep(1:ncol(sims), each=nrow(sims)))
  fullDF$plotVar = c(sims)
  plotVar = "plotVar"
  
  # get relevant map data
  states <- map_data("state")
  west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
  canada = map_data("world", "Canada")
  
  if(!is.data.frame(fullDF))
    fullDF = data.frame(fullDF)
  
  # rename row$Fault so that each fault gets unique name
  fullDF$Fault=1:nrow(fullDF)
  
  # make fault polygons
  faultDat = getFaultPolygons(fullDF)
  faultDat = merge(faultDat, fullDF[,c("Fault", plotVar)], by=c("Fault"))
  faultDat$plotVar=faultDat[,plotVar]
  faultDat$sim = rep(fullDF$sim, each=4)
  
  # grey maps plot:
  bg = ggplot(faultDat, aes(x=longitude, y=latitude)) + 
    geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
    coord_fixed(xlim = xlim,  ylim = ylim, ratio = 1.3, expand=FALSE) + 
    labs(x=xlab, y=ylab) + scale_x_continuous("Longitude", c(-127, -125, -123), labels=c("-127", "", "-123"), limits=c(-360, 360)) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='lightblue1'), strip.background = element_blank(),
          strip.text.x = element_blank())
  
  # generate fault polygons portion of plot
  faultPoly = geom_polygon(aes(fill=plotVar, group=factor(Fault)), color="black")
  
  if(logScale)
    pl = bg + faultPoly + facet_wrap(~sim, nrow=nrow, ncol=ncol) + scale_fill_distiller(clab, palette = "Spectral", direction=-1, trans="log") + 
    ggtitle(main)
  else
    pl = bg + faultPoly + facet_wrap(~sim, nrow=nrow, ncol=ncol) + scale_fill_distiller(clab, palette = "Spectral", direction=-1) +
    ggtitle(main)
  # ii <- cut(values, breaks = seq(min(values), max(values), len = 100), 
  #           include.lowest = TRUE)
  # ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
  # colors <- colorRampPalette(c("lightblue", "blue"))(99)[ii]
  # if(logScale)
  #   faultPoly = faultPoly + scale_fill_manual("", palette = "Spectral", direction=-1, trans="log")
  # else
  #   faultPoly = faultPoly + scale_fill_distiller("", palette = "Spectral", direction=-1)
  
  pl + guides(color=FALSE) 
}

# a simplified ggplot version of comparePredsToSubs
ggComparePredsToSubs = function(params, slipPreds=NULL, slipPredsGPS=NULL, subPreds=NULL, 
                                subPredsGPS=NULL, nsim=100, plotNameRoot="full", 
                                savePlots=TRUE, G=NULL, fileNameRoot=plotNameRoot, 
                                muVec=NULL, useGPS=FALSE, tvec=NULL, subDat=dr1, 
                                logScale=FALSE, fault=csz, latRange=c(40, 50), 
                                posNormalModel=FALSE, normalModel=posNormalModel, doGPSPred=FALSE, 
                                useMVNApprox=FALSE, taperedGPSDat=FALSE, dStar=28000, normalizeTaper=FALSE, 
                                noTitle=FALSE, lwd=.5) {
  # get parameters
  if(is.null(muVec)) {
    lambdaMLE = params[1]
    muZetaMLE = params[2]
    sigmaZetaMLE = params[3]
    muXi = params[5]
    muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
    muZetaCSZ = rep(muZetaMLE, nrow(fault))
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
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaMLE, dStar=dStar, normalize=normalizeTaper)
  
  #
  if(taperedGPSDat)
    phiZeta = params[length(params)]
  else
    phiZeta = NULL
  
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = preds(params, nsim=nsim, fault=fault, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec, 
                      posNormalModel=posNormalModel, normalModel=normalModel, phiZeta=phiZeta)
  if(is.null(slipPredsGPS) && doGPSPred)
    slipPredsGPS = predsGivenGPS(params, nsim=nsim, muVec=c(muZetaGPS, muZetaCSZ), fault=fault, tvec=tvec, 
                                 posNormalModel=posNormalModel, normalModel=normalModel)
  if(is.null(subPreds)) {
    if(is.null(G))
      subPreds = predsToSubsidence(params, slipPreds, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPreds = predsToSubsidence(params, slipPreds, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                   posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  if(is.null(subPredsGPS) && doGPSPred) {
    if(is.null(G))
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    else
      subPredsGPS = predsToSubsidence(params, slipPredsGPS, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                      posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
  }
  meanSlip = slipPreds$meanSlip
  meanSub = subPreds$meanSub
  u95 = subPreds$u95
  l95 = subPreds$l95
  u95Noise = subPreds$u95Noise
  l95Noise = subPreds$l95Noise
  slipSD = apply(slipPreds$slipSims, 1, sd)
  if(doGPSPred) {
    meanSlipGPS = slipPredsGPS$meanSlip
    meanSubGPS = subPredsGPS$meanSub
    u95GPS = subPredsGPS$u95
    l95GPS = subPredsGPS$l95
    u95NoiseGPS = subPredsGPS$u95Noise
    l95NoiseGPS = subPredsGPS$l95Noise
  }
  
  ##### generate mean seaDef field from Okada model
  # set Okada subsidence grid
  lonRange=c(-128, -122.5)
  nx = 300
  ny=  900
  lonGrid = seq(lonRange[1], lonRange[2], l=nx)
  latGrid = seq(latRange[1], latRange[2], l=ny)
  coordGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
  
  ##### plot subsidence predictions
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "ggSubsidencePredictions.pdf"), width=8, height=10)
  
  # slip mean
  if(!logScale) {
    pl1 = ggPlotFaultDat(fault, plotVar=meanSlip, main=paste0(plotNameRoot, "Mean Slip (m)"), 
                   xlim=lonRange, ylim=latRange, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
  }
  else {
    pl1 = ggPlotFaultDat(fault, plotVar=meanSlip, main=paste0(plotNameRoot, "Mean Slip (m)"), 
              xlim=lonRange, ylim=latRange, logScale=TRUE, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
  }
  
  # slip SE
  pl2 = ggPlotFaultDat(fault, plotVar = slipSD, main=paste0(plotNameRoot, "Slip SD (m)"), 
            xlim=lonRange, ylim=latRange, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
  
  # simulated subsidence data from Okada model using marginal distribution
  subRange = range(c(-meanSub, -l95Noise, -u95Noise, subDat$subsidence))
  ord = order(subDat$Lat)
  ordDat = subDat[ord,]
  ordL95 = l95Noise[ord]
  ordU95 = u95Noise[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl3 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  # plot magnitude distribution
  mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault, dStar=dStar, normalizeTaper=normalizeTaper)
  cleanMags = mags[is.finite(mags)]
  pl4 = qplot(cleanMags) + labs(x="Magnitudes", y="Frequency") + 
    geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.975), linetype=2) + 
    geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.025), linetype=2) + 
    geom_vline(col="purple", xintercept=mean(cleanMags)) + 
    ggtitle(paste0(plotNameRoot, "Histogram of earthquake magnitudes")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  if(noTitle) {
    pl1 = pl1 + ggtitle(NULL)
    pl2 = pl2 + ggtitle(NULL)
    pl3 = pl3 + ggtitle(NULL)
    pl4 = pl4 + ggtitle(NULL)
  }
  
  # combine plots into one
  multiplot(pl1, pl3, pl2, pl4, layout=matrix(1:4, ncol=2))
  
  dev.off()
  # tmp = arrangeGrob(pl1, pl3, pl2, pl4, layout_matrix=matrix(1:4, ncol=2), heights=rep(41,4), widths=rep(4,4))
  # 
  # grid.arrange(pl1, pl3, pl2, pl4, layout_matrix=matrix(1:4, ncol=2), heights=1:4, widths=1:4)
  
  invisible(NULL)
}

# plots 2x2 grid of fault plots
ggplotFixedSlip = function(meanSlip, medSlip=NULL, l95, u95, slipSD=NULL, plotNameRoot="full", 
                         savePlots=TRUE, fileNameRoot=plotNameRoot, logScale=FALSE, 
                         event="All", subDat=dr1) {
  if(event == "All")
    inds = 1:nrow(subDat)
  else
    inds = as.character(subDat$event) == event
  sortDat = subDat[inds,]
  
  # mean
  pl1 = ggPlotFaultDat(csz, meanSlip, logScale=logScale, xlim=lonRange, ylim=latRange, 
                   main=paste0(plotNameRoot, " Mean Slip (m)"), clab="") + 
    geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  
  # median/sd
  if(is.null(slipSD)) {
    pl2 = ggPlotFaultDat(csz, medSlip, logScale=logScale, xlim=lonRange, ylim=latRange, 
                         main=paste0(plotNameRoot, " Median Slip (m)"), clab="") + 
      geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  }
  else if(is.null(medSlip)) {
    pl2 = ggPlotFaultDat(csz, slipSD, logScale=logScale, xlim=lonRange, ylim=latRange, 
                         main=paste0(plotNameRoot, " Slip SD (m)"), clab="") + 
      geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  }
  # 2.5th percentile
  pl3 = ggPlotFaultDat(csz, l95, logScale=logScale, xlim=lonRange, ylim=latRange, 
                       main=paste0(plotNameRoot, " 2.5th Percentile Slip (m)"), clab="") + 
    geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  
  ## 97.5th percentile
  pl4 = ggPlotFaultDat(csz, u95, logScale=logScale, xlim=lonRange, ylim=latRange, 
                       main=paste0(plotNameRoot, " 97.5th Percentile Slip (m)"), clab="") + 
    geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "ggSlipDistn.pdf"), width=8, height=10)
  
  ## now put all plots together onto grid
  multiplot(pl1, pl3, pl2, pl4, cols=2)
  
  if(savePlots)
    dev.off()
}

# plot subsidence predictions against each other
ggCompareSubs = function(params, 
                         subPreds1, subPreds2, subPreds3, subPreds4, 
                         subDat1=dr1, subDat2=subDat1, subDat3=subDat1, subDat4=subDat3, 
                         tvec=NULL, 
                         plotNameRoot1="full", plotNameRoot2="full", plotNameRoot3="full", plotNameRoot4="full", 
                         savePlots=TRUE, fileNameRoot="", 
                         logScale=FALSE, fault=csz, latRange=c(40, 50), 
                         posNormalModel=FALSE, normalModel=posNormalModel, 
                         useMVNApprox=FALSE, taperedGPSDat=FALSE, dStar=25000, 
                         normalizeTaper=FALSE, noTitle=FALSE) {
  
  # get parameters
  lambdaMLE = params[1]
  muZetaMLE = params[2]
  sigmaZetaMLE = params[3]
  muXi = params[5]
  muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
  muZetaCSZ = rep(muZetaMLE, nrow(fault))
  
  # get taper if necessary
  if(is.null(tvec))
    tvec = taper(getFaultCenters(fault)[,3], lambda=lambdaMLE, dStar=dStar, normalize=normalizeTaper)
  
  #
  if(taperedGPSDat)
    phiZeta = params[length(params)]
  else
    phiZeta = NULL
  
  meanSub1 = subPreds1$meanSub
  u95Noise1 = subPreds1$u95Noise
  l95Noise1 = subPreds1$l95Noise
  meanSub2 = subPreds2$meanSub
  u95Noise2 = subPreds2$u95Noise
  l95Noise2 = subPreds2$l95Noise
  meanSub3 = subPreds3$meanSub
  u95Noise3 = subPreds3$u95Noise
  l95Noise3 = subPreds3$l95Noise
  meanSub4 = subPreds4$meanSub
  u95Noise4 = subPreds4$u95Noise
  l95Noise4 = subPreds4$l95Noise
  
  ## Make plots
  
  #plot 1
  subRange = range(c(-meanSub1, -l95Noise1, -u95Noise1, subDat1$subsidence, 
                     -meanSub2, -l95Noise2, -u95Noise2, subDat2$subsidence, 
                     -meanSub3, -l95Noise3, -u95Noise3, subDat3$subsidence, 
                     -meanSub4, -l95Noise4, -u95Noise4, subDat4$subsidence))
  ord = order(subDat1$Lat)
  ordDat = subDat1[ord,]
  ordL95 = l95Noise1[ord]
  ordU95 = u95Noise1[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl1 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot1, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  #plot 2
  ord = order(subDat2$Lat)
  ordDat = subDat2[ord,]
  ordL95 = l95Noise2[ord]
  ordU95 = u95Noise2[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl2 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot2, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  #plot 3
  ord = order(subDat3$Lat)
  ordDat = subDat3[ord,]
  ordL95 = l95Noise3[ord]
  ordU95 = u95Noise3[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl3 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot3, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  #plot 4
  ord = order(subDat4$Lat)
  ordDat = subDat4[ord,]
  ordL95 = l95Noise4[ord]
  ordU95 = u95Noise4[ord]
  tmp = cbind(ordDat, ordL95, ordU95)
  pl4 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
    labs(x="Subsidence (m)", y="") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    ggtitle(paste0(plotNameRoot4, "95% Prediction Band")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
  
  # remove titles if necessary
  if(noTitle) {
    pl1 = pl1 + ggtitle(NULL)
    pl2 = pl2 + ggtitle(NULL)
    pl3 = pl3 + ggtitle(NULL)
    pl4 = pl4 + ggtitle(NULL)
  }
  
  ## save plots
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "ggCompareSubs.pdf"), width=8, height=10)
  
  # put all plots together onto grid
  multiplot(pl1, pl3, pl2, pl4, cols=2)
  
  if(savePlots)
    dev.off()
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, byrow=FALSE) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols), 
                     byrow=byrow)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

##### functions for plotting multiple fields on a grid

ggplotSubsidenceGrid = function(allSubs, allPlotNames=NULL, savePlots=TRUE, fileNameRoot="", 
                                event="All", subDat=dr1, nr=NULL, nc=2, byrow=TRUE, 
                                latRange=c(40,50)) {
  if(is.null(allPlotNames)) {
    for(i in 1:ncol(allSubs))
      allPlotNames[i] = list(NULL)
  }
  
  if(event == "All")
    inds = 1:nrow(subDat)
  else
    inds = as.character(subDat$event) == event
  sortDat = subDat[inds,]
  
  # mean
  plots = list()
  subRange = range(c(allSubs, sortDat$subsidence))
  for(i in 1:ncol(allSubs)) {
    sortDat$theseSubs = allSubs[,i]
    pl = ggplot() + 
      geom_point(aes(x=subsidence, y=Lat), col="red", shape=3, data=sortDat) +
      geom_point(aes(x=theseSubs, y=Lat), col="blue", size=.3, data=sortDat) + 
      coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='white')) + 
      guides(color=FALSE) + 
      labs(x="Subsidence (m)", y="Latitude")
    plots = c(plots, list(pl))
  }
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SubGrid.pdf"), width=8, height=10)
  
  ## now put all plots together onto grid
  multiplot(plotlist=plots, cols=nc, byrow=byrow)
  
  if(savePlots)
    dev.off()
  
}

# plots nr by nc grid of fault plots
# allSlips is a list of slips to put on csz fault
# nr and nc is number of rows and columns of grid
# byrow is whether to put plots in row-major or column-major order
ggplotSlipGrid = function(allSlips, allPlotNames=NULL, savePlots=TRUE, 
                          fileNameRoot="", logScale=FALSE, 
                          event="All", subDat=dr1, nr=NULL, nc=2, byrow=TRUE, 
                          lwd=.5) {
  if(is.null(allPlotNames)) {
    for(i in 1:ncol(allSlips))
      allPlotNames[i] = list(NULL)
  }
  
  if(event == "All")
    inds = 1:nrow(subDat)
  else
    inds = as.character(subDat$event) == event
  sortDat = subDat[inds,]
  
  # mean
  plots = list()
  slipRange = range(allSlips)
  for(i in 1:ncol(allSlips)) {
    pl = ggPlotFaultDat(csz, allSlips[,i], varRange=slipRange, logScale=logScale, 
                        xlim=lonRange, ylim=latRange, 
                        main=allPlotNames[[i]], clab="", lwd=lwd) + 
      geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=sortDat)
    plots = c(plots, list(pl))
  }
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "SlipGrid.pdf"), width=8, height=10)
  
  ## now put all plots together onto grid
  multiplot(plotlist=plots, cols=nc, byrow=byrow)
  
  if(savePlots)
    dev.off()
  else
    return(plots)
}


ggplotSplineUncertainty = function(splinePar, covMat, nKnots=5, niter=1000, latRange=c(40,50), 
                                 diffGPSTaper=FALSE, nKnotsGPS=5, latsOnX=TRUE, 
                                 main=TeX("$\\lambda$ 95% Confidence Band"), uncertaintyBands=TRUE) {
  
  lats = seq(latRange[1], latRange[2], l=100)
  Xi = getSplineBasis(csz, latRange, nKnots, lats)
  if(diffGPSTaper) {
    XiGPS = getSplineBasis(csz, latRange, nKnotsGPS, lats)
    Xi = cbind(Xi, -XiGPS)
  }
  
  # draw simulations of the taper function
  cntr = Xi %*% splinePar
  if(uncertaintyBands) {
    L = t(chol(covMat))
    Zs = matrix(rnorm(ncol(L)*niter), nrow=ncol(L), ncol=niter)
    sims = Xi %*% L %*% Zs
    lows = apply(sims, 1, quantile, probs=0.025)
    his = apply(sims, 1, quantile, probs=0.975)
    lows = cntr + lows
    his = cntr + his
  }
  
  if(uncertaintyBands)
    lambdaRange =c(min(lows), max(his))
  else
    lambdaRange = range(cntr)
  
  # plot results
  if(latsOnX) {
    pl= ggplot() + 
      geom_path(aes(y=cntr, x=lats), col="blue") +
      geom_hline(col="black", yintercept = 0) +
      guides(col=FALSE) + ggtitle(main) + labs(x=TeX("$\\lambda$"), y="Latitude") + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='white'))
    
    if(uncertaintyBands)
      pl = pl + geom_path(aes(y=his, x=lats, linetype=2), col="blue") + 
      geom_path(aes(y=lows, x=lats, linetype=2), col="blue")
  }
  else {
    pl = ggplot() + 
      geom_path(aes(x=cntr, y=lats), col="blue") +
      coord_cartesian(xlim=lambdaRange, ylim=latRange) + 
      geom_vline(col="black", xintercept = 0) +
      guides(col=FALSE) + ggtitle(main) + labs(x=TeX("$\\lambda$"), y="Latitude") + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='white'))
    
    if(uncertaintyBands)
      pl = pl + geom_path(aes(x=lows, y=lats), col="blue", linetype=2L) +
      geom_path(aes(x=his, y=lats), col="blue", linetype=2L)
  }
  
  pl
}


# compare mean slip, sub preds, mags of different models
ggCompareModels = function(modelFitList, 
                           nsim=100, plotNameRoot="", savePlots=TRUE, 
                           G=NULL, fileNameRoot=plotNameRoot, muVec=NULL, subDat=dr1, 
                           logScale=FALSE, fault=csz, latRange=c(40, 50), 
                           posNormalModelVec=rep(FALSE, length(modelFitList)), 
                           normalModelVec=rep(TRUE, length(modelFitList)), 
                           useMVNApprox=FALSE, noTitle=TRUE, taperedGPSDat=TRUE, 
                           magRange=NULL, lwd=.5) {
  
  plots = list()
  for(i in 1:length(modelFitList)) {
    fit = modelFitList[[i]]
    params = fit$MLEs
    normalModel = normalModelVec[i]
    posNormalModel = posNormalModelVec[i]
    tvec = fit$tvec
    slipPreds = NULL
    subPreds = NULL
    
    # get parameters
    if(is.null(muVec)) {
      lambdaMLE = params[1]
      muZetaMLE = params[2]
      sigmaZetaMLE = params[3]
      muXi = params[5]
      muZetaGPS = rep(muZetaMLE, nrow(slipDatCSZ))
      muZetaCSZ = rep(muZetaMLE, nrow(fault))
    }
    else {
      lambdaMLE = params[1]
      sigmaZetaMLE = params[3]
      muXi = params[5]
      muZetaMLE = muVec
      muZetaGPS = muVec[1:nrow(slipDatCSZ)]
      muZetaCSZ = muVec[(nrow(slipDatCSZ)+1):length(muVec)]
    }
    
    #
    if(taperedGPSDat)
      phiZeta = params[length(params)]
    else
      phiZeta = NULL
    
    # generate predictions if they are left NULL by the user
    if(is.null(slipPreds))
      slipPreds = preds(params, nsim=nsim, fault=fault, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec, 
                        posNormalModel=posNormalModel, normalModel=normalModel, phiZeta=phiZeta)
    if(is.null(subPreds)) {
      if(is.null(G))
        subPreds = predsToSubsidence(params, slipPreds, useMVNApprox = useMVNApprox, subDat=subDat, 
                                     posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
      else
        subPreds = predsToSubsidence(params, slipPreds, G=G, useMVNApprox = useMVNApprox, subDat=subDat, 
                                     posNormalModel=posNormalModel, normalModel=normalModel, tvec=tvec)
    }
    meanSlip = slipPreds$meanSlip
    meanSub = subPreds$meanSub
    u95 = subPreds$u95
    l95 = subPreds$l95
    u95Noise = subPreds$u95Noise
    l95Noise = subPreds$l95Noise
    slipSD = apply(slipPreds$slipSims, 1, sd)
    myQuant = function(xs) {
      obs = xs[1]
      xs = xs[-1]
      mean(xs <= obs)
    }
    subQuant = apply(cbind(subDat$subsidence, -subPreds$noiseSims), 1, myQuant)
    outOfBounds = (subQuant < .025) | (subQuant > .975)
    normResids = qnorm(p=subQuant)
    normResids[!is.finite(normResids)] = NA
    
    ##### generate mean seaDef field from Okada model
    # set Okada subsidence grid
    lonRange=c(-128, -122.5)
    nx = 300
    ny=  900
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    coordGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
    
    ##### plot subsidence predictions
    
    # slip mean
    if(!logScale) {
      pl1 = ggPlotFaultDat(fault, plotVar=meanSlip, main=paste0(plotNameRoot, "Mean Slip (m)"), 
                           xlim=lonRange, ylim=latRange, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
    }
    else {
      pl1 = ggPlotFaultDat(fault, plotVar=meanSlip, main=paste0(plotNameRoot, "Mean Slip (m)"), 
                           xlim=lonRange, ylim=latRange, logScale=TRUE, clab="", lwd=lwd) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
    }
    plots = c(plots, list(pl1))
    
    # simulated subsidence data from Okada model using marginal distribution
    subRange = range(c(-meanSub, -l95Noise, -u95Noise, subDat$subsidence))
    ord = order(subDat$Lat)
    ordDat = subDat[ord,]
    ordL95 = l95Noise[ord]
    ordU95 = u95Noise[ord]
    # tmp = cbind(ordDat, ordL95, ordU95)
    # pl2 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    #   scale_y_continuous("Latitude", limits=latRange) + 
    #   labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    #   geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    #   ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
    #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(), 
    #         panel.background = element_rect(fill='white'))
    allShapes = rep(19, length(normResids))
    allShapes[outOfBounds] = 18
    allShapes = allShapes[ord]
    normResids = normResids[ord]
    aboveBounds = subQuant[ord] > .975
    belowBounds = subQuant[ord] < .025
    outOfBounds = outOfBounds[ord]
    tmp = cbind(ordDat, ordL95, ordU95, normResids, outOfBounds, aboveBounds, belowBounds)
    pl2 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col=normResids, fill=normResids), shape=19, size=.3, data=tmp[!outOfBounds,]) + 
      geom_point(aes(x=subsidence, y=Lat), shape=17, col="purple", size=.5, data=tmp[aboveBounds,], inherit.aes=FALSE) +
      geom_point(aes(x=subsidence, y=Lat), shape=17, col="green", size=.5, data=tmp[belowBounds,], inherit.aes=FALSE) +
      scale_y_continuous("Latitude", limits=latRange) + 
      labs(x="Subsidence (m)", y="Latitude") + guides(shape=FALSE, fill=FALSE) + 
      scale_color_distiller("", palette = "RdBu", direction=-1, limits=c(-qnorm(.975), qnorm(.975))) + 
      ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='white'))
    # geom_point(aes(x=subsidence, y=Lat, fill=normResids), shape=17, col="black", size=.3, data=tmp[outOfBounds,], inherit.aes=FALSE) +
    # scale_shape_manual(values=c(19, 17)) + 
    plots = c(plots, list(pl2))
    
    
    # plot magnitude distribution
    mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault, dStar=dStar, normalizeTaper=normalizeTaper)
    cleanMags = mags[is.finite(mags)]
    if(is.null(magRange))
      magRange=range(cleanMags)
    pl3 = qplot(cleanMags, xlim=magRange) + labs(x="Magnitudes", y="Frequency") + 
      geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.975), linetype=2) + 
      geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.025), linetype=2) + 
      geom_vline(col="purple", xintercept=mean(cleanMags)) + 
      ggtitle(paste0(plotNameRoot, "Histogram of earthquake magnitudes")) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='white'))
    plots = c(plots, list(pl3))
  }
  
  if(noTitle) {
    for(i in 1:length(plots)) {
      plots[[i]] = plots[[i]] + ggtitle(NULL)
    }
  }
  
  if(savePlots)
    pdf(file=paste0(fileNameRoot, "ModelComparison.pdf"), width=8, height=10)
  
  # combine plots into one
  multiplot(plotlist = plots, cols=3, byrow=TRUE)
  
  if(savePlots)
    dev.off()
  # tmp = arrangeGrob(pl1, pl3, pl2, pl4, layout_matrix=matrix(1:4, ncol=2), heights=rep(41,4), widths=rep(4,4))
  # 
  # grid.arrange(pl1, pl3, pl2, pl4, layout_matrix=matrix(1:4, ncol=2), heights=1:4, widths=1:4)
  
  invisible(NULL)
}

# function for plotting model locking normalized residuals versus latitude
ggplotLockingResiduals = function(modelFit, tvecGPS, gpsDat, latRange=c(40,50), 
                                  main=NULL) {
  
  params = modelFit$optPar
  muZ = params[1]
  sigmaZ = params[2]
  gamma = modelFit$gammaEst
  preds = (muZ * gamma) * tvecGPS
  resids = gpsDat$slip - preds
  normalizedResids = resids/sqrt(gamma^2*sigmaZ^2*tvecGPS^2 + gpsDat$slipErr^2)
  lats = gpsDat$lat
  ggplot() + geom_point(aes(x=normalizedResids, y=lats), col="blue", size=.5) + 
    labs(x="Normalized Residuals", y="Latitude") + 
    geom_vline(col="black", xintercept=0) + 
    ggtitle(main) + coord_fixed(xlim=range(normalizedResids), ylim=latRange, expand=FALSE) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
}

# function for plotting model locking normalized residuals versus latitude
ggplotSubsidenceResiduals = function(modelFit, tvec, subDat, G, latRange=c(40,50), 
                                  main=NULL, fault=csz) {
  # get model parameters
  params = modelFit$optPar
  muZ = params[1]
  sigmaZ = params[2]
  phiZ = params[length(params)]
  
  # compute covariance matrix of T %*% Z
  coordsZ = cbind(fault$longitude, fault$latitude)
  distMatZ = rdist.earth(coordsZ, miles=FALSE)
  corMatZ = stationary.cov(coordsZ, Covariance="Matern", theta=phiZ,
                             onlyUpper=FALSE, distMat=distMatZ, smoothness=3/2)
  covMatZ = sigmaZ^2 * corMatZ
  covMatTZ = sweep(sweep(covMatZ, 1, tvec, "*"), 2, tvec, "*")
  
  # compute marginal variances of subsidences
  covMatSubs = G %*% covMatTZ %*% t(G) + diag(subDat$Uncertainty^2)
  sigmaSubs = sqrt(diag(covMatSubs))
  
  # compute predicted subsidences
  slipPreds = muZ * tvec
  subPreds = -(G %*% slipPreds)
  resids = subDat$subsidence - subPreds
  normalizedResids = resids/sigmaSubs
  lats = subDat$Lat
  ggplot() + geom_point(aes(x=normalizedResids, y=lats), col="blue", size=.5) + 
    labs(x="Normalized Residuals", y="Latitude") + 
    geom_vline(col="black", xintercept=0) + 
    ggtitle(main) + coord_fixed(xlim=range(normalizedResids), ylim=latRange, expand=FALSE) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white'))
}


