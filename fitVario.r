setwd("~/git/M9")
source("loadTestData.r")
library(fields)

filterVG = function(VG, pct, minVG=-Inf, maxVG=Inf) {
  n = ceiling(length(VG$d)*pct)
  ind = (VG$vgram <= maxVG) & (VG$vgram >= minVG)
  VG$vgram = VG$vgram[ind]
  VG$d = VG$d[ind]
  if(length(VG$d) <= n) {
    return(VG)
  }
  newInd = sample(1:length(VG$d), n)
  VG$d = VG$d[newInd]
  VG$vgram = VG$vgram[newInd]
  return(VG)
}

meanVG = function(VG, minD=-Inf, maxD=Inf, statFun=mean, ...) {
  ind = (VG$d > minD) & (VG$d < maxD)
  do.call(statFun, c(list(VG$vgram[ind]), list(...)))
}

plotVGMean = function(x, N = 10, breaks = pretty(x$d, N, eps.correct = 1), 
                      add = FALSE, ...) 
{
  otherArgs = list(...)
  type = x$type
  if (is.null(otherArgs$ylab)) {
    if (type == "variogram")
      ylab = "sqrt(Variance)"
    else if (type == "covariogram" || type == "cross-covariogram") 
      ylab = "Covariance"
    else if (type == "correlogram" || type == "cross-correlogram") 
      ylab = "Correlation"
    else stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  else {
    ylab = otherArgs$ylab
    otherArgs$ylab = NULL
  }
  if (is.null(otherArgs$xlab)) 
    xlab = "Distance (Miles)"
  else {
    xlab = otherArgs$xlab
    otherArgs$xlab = NULL
  }
  if (is.null(otherArgs$main)) {
    if (type == "variogram") 
      main = "Empirical Variogram"
    else if (type == "covariogram") 
      main = "Empirical Covariogram"
    else if (type == "correlogram") 
      main = "Empirical Correlogram"
    else if (type == "cross-covariogram") 
      main = "Empirical Cross-Covariogram"
    else if (type == "cross-correlogram") 
      main = "Empirical Cross-Correlogram"
    else stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  else {
    main = otherArgs$main
    otherArgs$main = NULL
  }
  if (is.null(otherArgs$ylim)) {
    if (type == "correlogram" || type == "cross-correlogram") 
      ylim = c(-1, 1)
    else ylim = NULL
  }
  else {
    ylim = otherArgs$ylim
    otherArgs$ylim = NULL
  }
  if (is.null(otherArgs$type)) 
    type = "o"
  else {
    type = otherArgs$type
    otherArgs$type = NULL
  }
  meansFromBreak = function(breakBounds = c(-Inf, Inf)) {
    meanVG(x, minD=breakBounds[1], maxD=breakBounds[2], na.rm=TRUE)
  }
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  ys = apply(breakBounds, 1, meansFromBreak)
  if(x$type == "variogram")
    ys=sqrt(ys)
  notNas = !is.na(ys)
  centers = centers[notNas]
  ys = ys[notNas]
  if (!add) 
    do.call(plot, c(list(centers, ys, main = main, xlab = xlab, 
                         ylab = ylab, type = type, ylim = ylim), otherArgs))
  else do.call(lines, c(list(centers, ys, main = main, xlab = xlab, 
                             ylab = ylab, type = type, ylim = ylim), otherArgs))
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

fitVario = function() {
  
  #generate vgram data
  R = 3959 #in miles
  lonExtent= c(235.79781, 235.82087)
  latExtent = c(41.739671,41.762726)
  CCLon = mean(lonExtent)
  CCLat = mean(latExtent)
  lonDist = cos(CCLat)*2*pi/360*R
  latDist = 2*pi*R/360
  
  X = lon*lonDist
  Y = lat*latDist
  gridL = matrix(c(X, Y), ncol=2)
  
  nX = dim(allHMax)[2]
  nY = dim(allHMax)[3]
  maxX = max(X)
  minX = min(X)
  maxY = max(Y)
  minY = min(Y)
  distPerCellX = (maxX - minX)/nX
  distPerCellY =  (maxY - minY)/nY
  
  #subtract regressed topography data from field first
  residHMax = allHMax
  for(i in 1:dim(allHMax)[1]) {
    floodVals = allHMax[i,,]
    model = lm(c(floodVals) ~ c(topo))
    residVals = residuals(model)
    residHMax[i,,] = array(c(residVals), dim=c(1,nX,nY))
  }
  
  #set distances and how many points to compare each point with for vgram
  rectDim = 15 #size of the square around each point containing points for comparison
  maxIndexDist = floor(rectDim/2)
  #maxDist = maxIndexDist*min(c(distPerCellX, distPerCellY))
  
  #these are the indices each point is connected to in the general case, 
  #disregarding points near the edge of the grid.  Those are taken into account 
  #later with the removeMask variable.  Note that only indices beyond the current 
  #index are compared to avoid double-counting
  connectI = seq(0, nY*maxIndexDist, by=nY)
  connectI = rep(connectI, maxIndexDist+1) + rep(0:maxIndexDist, rep(maxIndexDist+1, maxIndexDist+1))
  xI = rep(0:maxIndexDist, maxIndexDist+1)
  yI = connectI %% nY
  xI = xI[connectI != 0]
  yI = yI[connectI != 0]
  connectI = connectI[connectI != 0] 
  
  #ID is the id matrix for vgram containing what points to compare to what. Generate 
  #all comparisons in this for loop
  ID = matrix(NA, nrow=nX*nY*length(connectI), ncol=2)
  removeMask = rep(FALSE, nX*nY*length(connectI))
  for(i in 1:(nX*nY)) {
    thisConnectI = connectI + i
    thisXI = ((i-1) %/% nY) + 1 #between 1 and nX
    thisYI = ((i-1) %% nY) + 1 #between 1 and nY
    
    #change from torus to R2 topology by removing certain comparisons
    thisRemoveMask = rep(FALSE, length(thisConnectI))
    if(thisXI > nX - maxIndexDist)
      thisRemoveMask = thisRemoveMask | xI > nX - thisXI
    if(thisYI > nY - maxIndexDist)
      thisRemoveMask = thisRemoveMask | yI > nY - thisYI
    #add the connections/comparisons for point i to the ID matrix
    startI = (i-1)*length(connectI) + 1
    endI = i*length(connectI)
    
    ID[startI:endI,] = cbind(i, thisConnectI)
    removeMask[startI:endI] = thisRemoveMask
  }
  
  #remove bad rows in ID matrix
  ID = ID[!removeMask,]
  # should end up having 
  #243*243*63 + (2*243)*((1+2+3+4+5+6+7)*8-7) + 7^2*(7^2-1)/2 
  #rows?
  
  #calculate vgram for each slice of HMax: each tsunami realization
  for(i in 1:dim(residHMax)[1]) {
    floodVals = c(residHMax[i,,])
    if(i == 1) {
      VG= vgram(gridL, floodVals, id=ID, type="variogram")
      coVG = vgram(gridL, floodVals, id=ID, type="covariogram")
    }
    else {
      #concetenate vgram
      VGslice = vgram(gridL, floodVals, id=ID, type="variogram")
      coVGslice = vgram(gridL, floodVals, id=ID, type="covariogram")
      VG$d = c(VG$d, VGslice$d)
      VG$vgram = c(VG$vgram, VGslice$vgram)
      coVG$d = c(coVG$d, coVGslice$d)
      coVG$vgram = c(coVG$vgram, coVGslice$vgram)
    }
  }
  
  #fit exponential variogram to data
  maxDist = maxIndexDist*min(c(distPerCellX, distPerCellY))
  s = mean(VG$vgram[VG$d > quantile(VG$d, .9)]) #sill
  n = mean(VG$vgram[VG$d < quantile(VG$d, .1)]) #nugget
  r = maxDist # range
  ys = VG$vgram
  ds = list(ds=VG$d)
  lower = list(s=0.00001, n=0.00001, r=.00001)
  #this is the variogram function. Go to http://www.seas.upenn.edu/~ese502/NOTEBOOK/Part_II/4_Variograms.pdf 
  #to get covariogram from variogram
  fit = nls(ys ~ (s - n)*(1 - exp(-(ds)/r)) + n, start=list(s=s, n=n, r=r), lower=lower, data=ds, algorithm="port")
  #lower = list(a=.00001, r=.00001)
  #a = mean(coVG$vgram[coVG$d < quantile(coVG$d, .1)]) #intercept
  #start = list(a=a, r=r)
  #covFit = nls(ys ~ a*exp(-(ds)/r) + n, start=start, lower=lower, data=ds, algorithm="port")
  summary(fit)
  
  #get variogram coefficients
  coefs = coef(fit)
  s = coefs[1]
  n = coefs[2]
  r = coefs[3]
  
  #plot variogram and covariogram fits
  expVGram = function(h) {
    (s - n)*(1 - exp(-h/r)) + n
  }
  expCoVGram = function(h) {
    a*exp(-h/r2)
  }
  expCoVGram2 = function(h) {
    (s - n)*exp(-h/r)
  }
  xs = seq(0, maxDist, length=500)
  
  filteredCoVG = filterVG(coVG, pct=.1)
  filteredVG = filterVG(VG, pct=.1)
  
  pdf("expVGramPlot.pdf", height=5, width=7)
  plotVGMean(filteredVG, main="Empirical and Exponential Variogram Fit")
  lines(xs, sqrt(expVGram(xs)), col="green")
  dev.off()
  
#   pdf("expVGramBoxplot.pdf", height=5, width=7)
#   boxplotVGram(filteredVG, main="Empirical and Exponential Variogram Fit")
#   lines(xs, sqrt(expVGram(xs)), col="green")
#   dev.off()
  
  pdf("expCoVGramPlot.pdf", height=5, width=7)
  plotVGMean(filteredCoVG, main="Empirical and Exponential Covariogram Fit")
  lines(xs, expCoVGram(xs), col="green")
  dev.off()

#   pdf("expCoVGramBoxplot.pdf", height=5, width=7)
#   boxplotVGram(filteredCoVG, main="Empirical and Exponential Covariogram Fit")
#   lines(xs, expCoVGram(xs), col="green")
#   dev.off()
  
  save(s, n, r, VG, coVG, file="fitExpVarioParams.RData")
}


