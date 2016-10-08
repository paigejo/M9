# miscellaneous functions used for exploratory analysis

# makes plots comparing the predicted subsidence levels to the subsidence data.
# if prediction inputs are left NULL, they are computed using params
comparePredsToSubs = function(params, slipPreds=NULL, subPreds=NULL, nsim=100, 
                              plotNameRoot="full", savePlots=TRUE) {
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = predsGivenGPS(params, nsim=nsim)
  if(is.null(subPreds))
    subPreds = predsToSubsidence(params, slipPreds)
  meanSlip = slipPreds$meanSlip
  meanSub = subPreds$meanSub
  u95 = subPreds$u95
  l95 = subPreds$l95
  
  # get fit MLEs
  lambda = params[1]
  muZeta = params[2]
  sigmaZeta = params[3]
  lambda0 = params[4]
  muXi = params[5]
  
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
    pdf(file=paste0(plotNameRoot, "SubsidencePredictions.pdf"), width=8, height=10)
  par(mfrow=c(2,2))
  plotFault(csz, plotVar = meanSlip, legend.mar=6, main="Mean Slip (m)", 
            xlim=lonRange)
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  inds = event == "T1"
  points(dr1$Lon[inds], dr1$Lat[inds], cex=1, pch=3, col="red")
  # seadef from Okada model generated from the slips
  plotTopo(list(lon=lonGrid, lat=latGrid, dat=t(meanSeaDef)), nx=120, ny=120, main="Seafloor Deformation (m)")
  map("world", "Canada", add=TRUE, lwd=1)
  US(add=TRUE, lwd=1)
  points(dr1$Lon, dr1$Lat, cex=1, pch=3, col="red")
  plotFault(csz, new=FALSE, lwd=1, plotData=FALSE)
  # T1 subsidence data
  subRange = range(c(-meanSub, -l95, -u95, dr1$subsidence[inds]))
  plot(dr1$subsidence, dr1$Lat, pch=19, cex=.3, ylim=latRange, 
       xlim=subRange, main="Observed Subsidence", ylab="", xlab="Subsidence (m)")
  # simulated subsidence data from Okada model using GPS locking rates
  plot(-meanSub, dr1$Lat, pch=19, cex=.3, ylim=latRange, col="blue", 
       xlim=subRange, main="Simulated Subsidence", ylab="", xlab="Subsidence (m)")
  ord = order(dr1$Lat)
  lines(-l95[ord], dr1$Lat[ord])
  lines(-u95[ord], dr1$Lat[ord])
  if(savePlots)
    dev.off()
  
  ##### plot GPS predictions
  # first generate the GPS predictions
  mucGPS= slipPreds$mucGPS
  SigmaDiagGPS = diag(slipPreds$SigmacGPS)
  sigmaXi = sqrt(log(.5*(sqrt(4*slipDatCSZ$slipErr^2/slipDatCSZ$slip^2 + 1) + 1)))
  meanPredsGPS = exp(mucGPS + muXi + SigmaDiagGPS/2 + sigmaXi^2/2)
  
  # Now make the plot
  if(savePlots)
    pdf(file=paste0(plotNameRoot, "GPSPredictions.pdf"), width=7, height=5)
  par(mfrow=c(1,2))
  # plot mean slip field times Xi (mean locking rate)
  datRange = range(c(meanPredsGPS, slipDatCSZ$slip))
  quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, meanPredsGPS, main="Predicted Mean Locking Rate (mm/yr)", 
             xlim=lonRange, ylim=latRange, zlim=datRange)
  plotFault(csz, new=FALSE, plotData = FALSE)
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  # plot GPS data
  quilt.plot(slipDatCSZ$lon, slipDatCSZ$lat, slipDatCSZ$slip, main="Observed Locking Rate (mm/yr)", 
             xlim=lonRange, ylim=latRange, zlim=datRange)
  plotFault(csz, new=FALSE, plotData = FALSE)
  map("world", "Canada", lwd=1, add=TRUE)
  US(add=TRUE, lwd=1)
  if(savePlots)
    dev.off()
  
  invisible(NULL)
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

getMomentFromSlip = function(slips, rigidity=4*10^10) {
  # get fault total area
  areas = csz$length*csz$width
  totArea = sum(areas)
  
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