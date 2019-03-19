# plotters for TMB fits
# a simplified ggplot version of comparePredsToSubs
ggComparePredsToSubsTMB = function(modelInfo, slipPreds=NULL, subPreds=NULL, 
                                nsim=100, plotNameRoot="fullTMB", 
                                savePlots=TRUE, G=NULL, fileNameRoot=plotNameRoot, 
                                muVec=NULL, subDat=dr1, 
                                gpsDat=slipDatCSZ, logScale=FALSE, fault=csz, latRange=c(40, 50), 
                                posNormalModel=FALSE, pts=cbind(gpsDat$lon, gpsDat$lat), 
                                noTitle=FALSE, lwd=.5, magLimits=c(8.5, 9.5), binwidth=NULL) {
  
  finalPar = modelInfo$finalPar
  
  # set other relevant parameters
  phiZeta = exp(modelInfo$logphiEst)
  alpha = exp(modelInfo$logalphaEst)
  nuZeta = 3/2
  dStarGPS = modelInfo$data$dStarGPS
  dStar = modelInfo$data$dStar
  
  highInflate = modelInfo$loghighInflateEst
  lowInflate = modelInfo$loglowInflateEst
  
  optParNames = names(modelInfo$opt$par)
  minPar = modelInfo$opt$par
  betaTaperEst = modelInfo$betaTaperEst
  betaTaperGPSEst = minPar[which(optParNames == "betaTaperGPS")]
  betasdEst = modelInfo$betasdEst
  betaMeanEst = modelInfo$betaMeanEst
  betaGammaEst = modelInfo$betaGammaEst
  betaGammaGPSEst = modelInfo$betaGammaGPSEst
  
  nKnots = length(modelInfo$betaTaperEst)
  nKnotsGPS = ncol(modelInfo$data$lambdaBasisXGPS)
  nKnotsVar = length(modelInfo$betasdEst)
  nKnotsGamma = length(modelInfo$betaGammaEst)
  # muZeta = exp(modelInfo$logmu)
  nKnotsMean = length(modelInfo$betaMeanEst)
  
  diffGPSTaper = length(modelInfo$betaTaperGPSEst)
  diffMean = length(betaGammaGPSEst) != 0
  
  # generate spline basis matrix
  lambdaBasisY = getSplineBasis(fault, nKnots=nKnots, latRange=latRange)
  lambdaBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnots, latRange=latRange)
  lambdaBasisXGPS = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsGPS, latRange=latRange)
  meanBasisY = getSplineBasis(fault, nKnots=nKnotsMean, latRange=latRange)
  meanBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsMean, latRange=latRange)
  if(diffMean)
    meanBasisXGPS = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsMeanGPS, latRange=latRange)
  sdBasisY = getSplineBasis(fault, nKnots=nKnotsVar, latRange=latRange)
  sdBasisX = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsVar, latRange=latRange)
  gammaBasis = getSplineBasis(data.frame(list(latitude=pts[,2])), nKnots=nKnotsGamma, latRange=latRange)
  
  faultDepths = getFaultCenters(fault)[,3]
  xDepths = gpsDat$Depth
  
  # evaluate splines on the fault and for the points of interest
  taperVecY = c(taper(faultDepths, exp(lambdaBasisY %*% betaTaperEst), dStar = dStar))
  sdVecX = exp(sdBasisX %*% betasdEst)
  sdVecY = exp(sdBasisY %*% betasdEst)
  if(!diffMean)
    meanVecX = exp(meanBasisX %*% betaMeanEst)
  else
    meanVecX = exp(meanBasisX %*% betaMeanEst + meanBasisXGPS %*% betaMeanGPSEst)
  meanVecY = exp(meanBasisY %*% betaMeanEst)
  gammaVec = exp(gammaBasis %*% betaGammaEst)
  
  # generate predictions if they are left NULL by the user
  if(is.null(slipPreds))
    slipPreds = predsTMB(modelInfo, nsim=nsim, fault=fault, gpsDat=gpsDat, posNormalModel=posNormalModel)
  if(is.null(subPreds)) {
    subPreds = predsToSubsidenceTMB(modelInfo, slipPreds, G=G, subDat=subDat, 
                                   posNormalModel=posNormalModel)
  }
  meanSlip = slipPreds$meanSlip
  meanSub = subPreds$meanSub
  u95 = subPreds$u95
  l95 = subPreds$l95
  u95Noise = subPreds$u95Noise
  l95Noise = subPreds$l95Noise
  slipSD = sqrt(diag(slipPreds$Sigma)) * taperVecY
  
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
  # subRange = range(c(-meanSub, -l95Noise, -u95Noise, subDat$subsidence))
  ord = order(subDat$Lat)
  ordDat = subDat[ord,]
  ordL95 = l95Noise[ord]
  ordU95 = u95Noise[ord]
  ordLat = subDat$Lat[ord]
  ordMeanSub = subPreds$meanSub[ord]
  # tmp = cbind(ordDat, ordL95, ordU95)
  # pl3 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
  #   coord_fixed(xlim=subRange, ylim=latRange, expand=FALSE) + 
  #   labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
  #   geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
  #   ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
  #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(), 
  #         panel.background = element_rect(fill='white'))
  tmp = data.frame(meanSub=ordMeanSub, Lat=ordLat)
  pl3 = ggplot() + geom_point(aes(x=subsidence, y=Lat), col="red", shape=3, size=.3, data=ordDat) + 
    geom_point(aes(x=-meanSub, y=Lat), col="blue", shape=19, size=.3, data=tmp) + 
    ggtitle(paste0(plotNameRoot, "Subsidence Predictions")) + 
    scale_y_continuous("Latitude", limits=latRange) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white')) + labs(x="Subsidence (m)", y="Latitude") + guides(shape=FALSE, fill=FALSE)
  
  # plot magnitude distribution
  mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault, normalizeTaper=TRUE)
  cleanMags = mags[is.finite(mags)]
  cleanMags = cleanMags[cleanMags != 0]
  tmp = data.frame(cleanMags=cleanMags)
  pl4 = ggplot(tmp) + geom_histogram(aes(cleanMags, y=..density..), binwidth=binwidth) + labs(x="Magnitudes", y="Density") + 
    geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.95), linetype=2) + 
    geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.05), linetype=2) + 
    geom_vline(col="purple", xintercept=mean(cleanMags)) + 
    ggtitle(paste0(plotNameRoot, "Histogram of magnitudes")) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='white')) + 
    scale_x_continuous(limits=magLimits)
  
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
ggplotFixedSlipTMB = function(meanSlip, medSlip=NULL, l95, u95, slipSD=NULL, plotNameRoot="full", 
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
ggCompareSubsTMB = function(params, 
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

ggplotSubsidenceGridTMB = function(allSubs, allPlotNames=NULL, savePlots=TRUE, fileNameRoot="", 
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
ggplotSlipGridTMB = function(allSlips, allPlotNames=NULL, savePlots=TRUE, 
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


ggplotSplineUncertaintyTMB = function(splinePar, covMat, nKnots=5, niter=1000, latRange=c(40,50), 
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

ggplotTaperDepthUncertaintyTMB = function(splinePar, covMat, nKnots=5, niter=1000, latRange=c(40,50), 
                                       nKnotsGPS=5, latsOnX=TRUE, depthFrac=.9, confLevel=.95, 
                                       diffGPSTaper=FALSE, main=TeX(paste0(depthFrac * 100, "% Taper Depth ", confLevel * 100, "% Confidence Band")), 
                                       uncertaintyBands=TRUE, dStar=21000) {
  
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
    lambdas = sweep(Xi %*% L %*% Zs, 1, cntr, "+")
    depths = getFracTaperDepth(lambdas, depthFrac, dStar=dStar)
    lows = apply(depths, 1, quantile, probs=(1 - confLevel) / 2)
    his = apply(depths, 1, quantile, probs=1 - (1 - confLevel) / 2)
    cntr = getFracTaperDepth(cntr, depthFrac, dStar=dStar)
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
      guides(col=FALSE) + ggtitle(main) + labs(x=TeX(paste0(depthFrac * 100, "% Taper Depth")), y="Latitude") + 
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
      guides(col=FALSE) + ggtitle(main) + labs(x=TeX(paste0(depthFrac * 100, "% Taper Depth")), y="Latitude") + 
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
ggCompareModelsTMB = function(modelFitList, 
                           nsim=100, plotNameRoot="", savePlots=TRUE, 
                           G=NULL, fileNameRoot=plotNameRoot, muVec=NULL, subDat=dr1, 
                           logScale=FALSE, fault=csz, latRange=c(40, 50), 
                           posNormalModelVec=rep(FALSE, length(modelFitList)), 
                           normalModelVec=rep(TRUE, length(modelFitList)), 
                           useMVNApprox=FALSE, noTitle=TRUE, taperedGPSDat=TRUE, 
                           magRange=NULL, lwd=.5, magTicks=NULL, plotSubPreds=FALSE, 
                           subPredMeanRange=c(-2,2), adjustedMeans=NULL, 
                           varRange=NULL, anisotropic=FALSE) {
  
  plots = list()
  for(i in 1:length(modelFitList)) {
    fit = modelFitList[[i]]
    params = fit$MLEs
    if(!is.null(adjustedMeans)) {
      params[2] = adjustedMeans[i]
    }
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
    
    # get the correlation range
    if(taperedGPSDat)
      phiZeta = params[length(params) - anisotropic]
    else
      phiZeta = NULL
    
    # generate predictions if they are left NULL by the user
    if(is.null(slipPreds))
      slipPreds = preds(params, nsim=nsim, fault=fault, muVec=c(muZetaGPS, muZetaCSZ), tvec=tvec, 
                        posNormalModel=posNormalModel, normalModel=normalModel, phiZeta=phiZeta, 
                        anisotropic=anisotropic, taperedGPSDat=taperedGPSDat)
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
                           xlim=lonRange, ylim=latRange, clab="", lwd=lwd, varRange=varRange) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
    }
    else {
      pl1 = ggPlotFaultDat(fault, plotVar=meanSlip, main=paste0(plotNameRoot, "Mean Slip (m)"), 
                           xlim=lonRange, ylim=latRange, logScale=TRUE, clab="", lwd=lwd, varRange=varRange) + geom_point(aes(x=Lon, y=Lat, col="red"), shape=3, data=subDat)
    }
    plots = c(plots, list(pl1))
    
    # simulated subsidence data from Okada model using marginal distribution
    subRange = range(c(-meanSub, -l95Noise, -u95Noise, subDat$subsidence))
    ord = order(subDat$Lat)
    ordDat = subDat[ord,]
    ordL95 = l95Noise[ord]
    ordU95 = u95Noise[ord]
    ordMeanSub = meanSub[ord]
    ordLat = subDat$Lat[ord]
    # tmp = cbind(ordDat, ordL95, ordU95)
    # pl2 = ggplot() + geom_point(aes(x=subsidence, y=Lat, col="red"), shape=3, data=tmp) + 
    #   scale_y_continuous("Latitude", limits=latRange) + 
    #   labs(x="Subsidence (m)", y="Latitude") + geom_path(aes(x=-ordL95, y=Lat), col="blue", data=tmp) +
    #   geom_path(aes(x=-ordU95, y=Lat), col="blue", data=tmp) + guides(col=FALSE) + 
    #   ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
    #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(), 
    #         panel.background = element_rect(fill='white'))
    
    if(!plotSubPreds) {
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
        scale_color_distiller("Z-Score", palette = "RdBu", direction=-1, limits=c(-qnorm(.975), qnorm(.975))) + 
        ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.background = element_rect(fill='white'))
    }
    else {
      tmp = data.frame(meanSub=ordMeanSub, Lat=ordLat)
      pl2 = ggplot() + geom_point(aes(x=subsidence, y=Lat), col="red", shape=3, size=.3, data=ordDat) + 
        geom_point(aes(x=-meanSub, y=Lat), col="blue", shape=19, size=.3, data=tmp) + 
        ggtitle(paste0(plotNameRoot, "95% Prediction Band")) + 
        scale_y_continuous("Latitude", limits=latRange) + 
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), 
              panel.background = element_rect(fill='white')) + labs(x="Subsidence (m)", y="Latitude") + guides(shape=FALSE, fill=FALSE) +
        scale_x_continuous(limits=subPredMeanRange)
    }
    # geom_point(aes(x=subsidence, y=Lat, fill=normResids), shape=17, col="black", size=.3, data=tmp[outOfBounds,], inherit.aes=FALSE) +
    # scale_shape_manual(values=c(19, 17)) + 
    plots = c(plots, list(pl2))
    
    
    # plot magnitude distribution
    mags = apply(slipPreds$slipSims, 2, getMomentFromSlip, fault=fault, dStar=dStar, normalizeTaper=normalizeTaper)
    cleanMags = mags[is.finite(mags) & mags != 0]
    if(is.null(magRange))
      magRange=range(cleanMags)
    tmp = data.frame(cleanMags=cleanMags)
    pl3 = ggplot(tmp) + geom_histogram(aes(cleanMags, y=..density..)) + labs(x="Magnitudes", y="Density") + 
      geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.975), linetype=2) + 
      geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.025), linetype=2) + 
      geom_vline(col="purple", xintercept=mean(cleanMags)) + 
      ggtitle(paste0(plotNameRoot, "Histogram of earthquake magnitudes")) + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='white'))
    # pl3 = qplot(cleanMags, xlim=magRange) + labs(x="Magnitudes", y="Frequency") + 
    #   geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.975), linetype=2) + 
    #   geom_vline(col="purple", xintercept=quantile(cleanMags, probs=.025), linetype=2) + 
    #   geom_vline(col="purple", xintercept=mean(cleanMags)) + 
    #   ggtitle(paste0(plotNameRoot, "Histogram of earthquake magnitudes")) + 
    #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
    #         panel.grid.minor = element_blank(), 
    #         panel.background = element_rect(fill='white'))
    
    if(!is.null(magTicks))
      pl3 = pl3 + scale_x_continuous(limits=magRange, breaks=magTicks)
    
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
ggplotLockingResidualsTMB = function(modelFit, tvecGPS, gpsDat, latRange=c(40,50), 
                                  main=NULL, doGammaSpline=FALSE) {
  
  params = modelFit$optPar
  muZ = params[1]
  sigmaZ = params[2]
  gamma = modelFit$gammaEst
  if(doGammaSpline)
    gamma = modelFit$gammaEst$gammaEstGps
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
ggplotSubsidenceResidualsTMB = function(modelFit, tvec, subDat, G, latRange=c(40,50), 
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

ggplotSplineUncertaintyTMB2 = function(modelInfo, parName=NULL, splineParI=NULL, latsOnX=TRUE, latRange=c(40,50), 
                                      main=TeX(paste0(parName, " 95% Confidence Band")), uncertaintyBands=TRUE, 
                                      useReport=FALSE) {
  if(!useReport) {
    allPar = modelInfo$opt$par
    if(is.null(splineParI))
      splineParI = grepl(parName, names(allPar))
    fullCovMat = solve(-modelInfo$hess)
    splinePar = allPar[splineParI]
  } else {
    allPar = modelInfo$report$par.fixed
    if(is.null(splineParI))
      splineParI = grepl(parName, names(allPar))
    fullCovMat = modelInfo$report$cov.fixed
    splinePar = allPar[splineParI]
  }
  
  modelInfo$report$cov
  lats = seq(latRange[1], latRange[2], l=500)
  Xi = getSplineBasis(csz, nKnots=length(splinePar), lats=lats, latRange=latRange)
  
  # get the conditional covariance, conditioning on the other parameters
  otherParI = (1:length(allPar))[-splineParI]
  otherPar = allPar[otherParI]
  out = conditionalNormal(otherPar, splinePar, otherPar, fullCovMat[splineParI, splineParI], fullCovMat[otherParI, otherParI], 
                          fullCovMat[splineParI, otherParI])
  muc = out$muc
  covMat = out$Sigmac
  
  # draw simulations of the taper function
  cntr = Xi %*% muc
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