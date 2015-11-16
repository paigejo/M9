setwd("~/git/M9")
source("loadTestData.r")
library(fields)

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
  
  #set distances and how many points to compare each point with for vgram
  nX = dim(allHMax)[2]
  nY = dim(allHMax)[3]
  maxX = max(X)
  minX = min(X)
  maxY = max(Y)
  minY = min(Y)
  distPerCellX = (maxX - minX)/nX
  distPerCellY =  (maxY - minY)/nY
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
  for(i in 1:dim(allHMax)[1]) {
    floodVals = c(allHMax[i,,])
    if(i == 1)
      VG= vgram(gridL, floodVals, id=ID)
    else {
      #concetenate vgram
      VGslice = vgram(gridL, floodVals, id=ID)
      VG$d = c(VG$d, VGslice$d)
      VG$vgram = c(VG$vgram, VGslice$vgram)
    }
  }
  
  #fit exponential variogram to data
  maxDist = maxIndexDist*min(c(distPerCellX, distPerCellY))
  s = mean(VG$vgram[VG$d > maxDist*.9]) #sill
  n = mean(VG$vgram[VG$d < maxDist*.1]) #nugget
  r = maxDist # range
  ys = VG$vgram
  ds = list(ds=VG$d)
  fit = nls(ys ~ (s - n)*(1 - exp(-(ds)/r)) + n, start=list(s=s, n=n, r=r), data=ds)
  summary(fit)
  
  #get variogram coefficients
  coefs = coef(fit)
  s = coefs$s
  n = coefs$n
  r = coefs$r
  
  #plot variogram fit
  expVGram = function(h) {
    (s - n)*(1 - exp(-h/r)) + n
  }
  xs = seq(0, maxDist, length=500)
  vgramFit = expVGram(xs)
  
  pdf("expVGram.pdf", height=5, width=7)
  plot(VG, main="Empirical and Exponential Variogram Fit")
  lines(xs, expVGram(xs), col=green)
  dev.off()
  
  save(s, n, r, file="fitExpVarioParams.RData")
}
