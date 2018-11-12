##### functions for making fault

# width and height of grid cells determined by lonLatGrid resolution. Corners of fault grid 
# given by points of lonLatGrid
makeRectFault = function(lonLatGrid=NULL, maxDepth=26000) {
  # make the grid over which to make pointwise estimates (~7040 points).  These will be the top centers of the cells
  if(is.null(lonLatGrid)) {
    thisSlipDat = slipDat[slipDat$Depth < maxDepth,]
    
    latRange=range(thisSlipDat$lat)
    lonRange=range(thisSlipDat$lon)
    nx = 80
    ny = 240
    lonGrid = seq(lonRange[1], lonRange[2], l=nx)
    latGrid = seq(latRange[1], latRange[2], l=ny)
    lonLatGrid = make.surface.grid(list(lon=lonGrid, lat=latGrid))
    lonLatGrid = data.frame(list(lon=lonLatGrid[,1], lat=lonLatGrid[,2]))
  }
  
  # get predicted depth at prediction locations
  phiZeta = 232.5722
  out = fastTps(cbind(slipDat$lon, slipDat$lat), slipDat$Depth, m=3, theta=phiZeta, lon.lat=TRUE, 
                Dist.args=list(miles=FALSE, method="greatcircle"))
  depths = predict(out, lonLatGrid)
  
  # calculate depthGradient
  getDepthGrad = function(i) {
    # i is the index of lonLatGrid
    lonI = ((i-1) %% nx) + 1
    latI = floor((i-1)/nx) + 1
    
    # first assume we aren't on the edge of the grid
    upI = i + nx
    downI = i - nx
    rightI = i + 1
    leftI = i - 1
    
    # if we are, modify any nonexisting point to be original point
    if(latI == ny)
      upI = i
    if(latI == 1)
      downI = i
    if(lonI == nx)
      rightI = i
    if(lonI == 1)
      leftI = i
    
    # now get depths, distances (in meters), and gradient
    uDepth = depths[upI]
    dDepth = depths[downI]
    rDepth = depths[rightI]
    lDepth = depths[leftI]
    uDist = rdist.earth(lonLatGrid[downI,], lonLatGrid[upI,], miles=FALSE) * 1000
    rDist = rdist.earth(lonLatGrid[leftI,], lonLatGrid[rightI,], miles=FALSE) * 1000
    uGrad = (uDepth - dDepth)/udWid
    rGrad = (rDepth - lDepth)/lrWid
    
  }
  
  # make sure we only use grid cells with positive depth
  posDepth = depths > 0
  lonLatGrid = lonLatGrid[posDepth,]
  depths = depths[posDepth,]
  
  
}