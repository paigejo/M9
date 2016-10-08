#functions for loading GeoClaw Output

# lonReg = c(360-124.20219, 360-124.1791325)
# latReg = c(41.73967, 41.7627285)
# lonReg = c(360-124.203, 360-124.179)
# latReg = c(41.739, 41.763)

#format:
#a bunch of pieces in the form:
#num grid_number
#num AMR_level
#num mx
#num my
#num xlow
#num ylow
#num dx
#num dy
#line
#dat
#line
loadGeoQ = function(fname) {
  dat = scan(fname, what="character")
  
  gridStart = 0
  ngrid = 0
  allgrid_number = list()
  allAMR_level = list()
  allmx = list()
  allmy = list()
  allxlow = list()
  allylow = list()
  alldx = list()
  alldy = list()
  alldat = list()
  alllon =list()
  alllat =list()
  allwh = list() #water height
  allxvel = list()
  allyvel = list()
  allsurfDisp = list() #surface displacement (water height - height at rest)
  
  while(gridStart < length(dat)) {
    ngrid = ngrid + 1
    
    #get header information
    allgrid_number[ngrid] <- grid_number <- as.numeric(dat[gridStart+1])
    allAMR_level[ngrid] <- AMR_level <- as.numeric(dat[gridStart+3])
    allmx[ngrid] <- mx <- as.numeric(dat[gridStart+5])
    allmy[ngrid] <- my <- as.numeric(dat[gridStart+7])
    allxlow[ngrid] <- xlow <- as.numeric(dat[gridStart+9])
    allylow[ngrid] <- ylow <- as.numeric(dat[gridStart+11])
    alldx[ngrid] <- dx <- as.numeric(dat[gridStart+13])
    alldy[ngrid] <- dy <- as.numeric(dat[gridStart+15])
    
    flood <- as.numeric(dat[(gridStart+16 + 1):(gridStart+16 + mx*my*4)])
    wh = array(flood[seq(from=1, length=mx*my, by=4)], dim=c(mx, my))
    xvel = array(flood[seq(from=2, length=mx*my, by=4)], dim=c(mx, my))
    yvel = array(flood[seq(from=3, length=mx*my, by=4)], dim=c(mx, my))
    surfDisp = array(flood[seq(from=4, length=mx*my, by=4)], dim=c(mx, my))
    allwh[[ngrid]] = wh
    allxvel[[ngrid]] = xvel
    allyvel[[ngrid]] = yvel
    allsurfDisp[[ngrid]] = surfDisp
    
    alllon[[ngrid]] <- lon <- seq(xlow, xlow+mx*dx, l=mx)
    alllat[[ngrid]] <- lat <- seq(ylow, ylow+my*dy, l=my) # don't reverse order because data doesn't start from the top
    
    gridStart = gridStart+16 + mx*my*4
  }
  
  out = list(grid_number=allgrid_number, AMR_level=allAMR_level, mx=allmx, my=allmy, 
             xlow=allxlow, ylow=allylow, dx=alldx, dy=alldy, lon=alllon, lat=alllat, 
             wh=allwh, xvel=allxvel, yvel=allyvel, surfDisp=allsurfDisp)
  return(out)
}

#format:
#in the form:
#num time
#num meqn
#num ngrids
#num naux
#num ndim
#num nghost
#(but all we care about is the time)
loadGeoT = function(fname) {
  dat = scan(fname, what="character")
  return(as.numeric(dat[1]))
}

# grid_number
# AMR_level
# mx 
# my 
# xlow
# ylow
# dx 
# dy
# lon 
# lat 
# wh
# xvel
# yvel
# surfDisp
plotGeoQ = function(geoQDat, main="Water Surface Displacement (meters)", zlim=NULL, 
                    land = FALSE, maxToPlot=5000) {
  #NOTE: land == FALSE ==> only plot points with water depth > 0
  len = length(geoQDat$mx)
  
  calcXRes = function(i) {
    thismx = geoQDat$mx[[i]]
    thislon = geoQDat$lon[[i]]
    return((max(thislon)-min(thislon))/thismx)
  }
  calcYRes = function(i) {
    thismy = geoQDat$my[[i]]
    thislat = geoQDat$lat[[i]]
    return((max(thislat)-min(thislat))/thismy)
  }
  if(!land) {
    for(i in 1:length(geoQDat$surfDisp)) {
      wh = geoQDat$wh[[i]]
      geoQDat$surfDisp[[i]][wh <= 0] = NA
    }
  }
  
  #calculate range of data values (lat, lon, water height, resolution)
  rangeLon = range(sapply(geoQDat$lon, range, na.rm=TRUE))
  rangeLat = range(sapply(geoQDat$lat, range, na.rm=TRUE))
  if(is.null(zlim))
    rangeDat = range(sapply(geoQDat$surfDisp, range, na.rm=TRUE, finite=TRUE), 
                     na.rm=TRUE, finite=TRUE)
  else
    rangeDat = zlim
  xres = sapply(1:len, calcXRes)
  yres = sapply(1:len, calcYRes)
  avgres = (xres + yres)/2
  
  #plot low-resolution data first, high resolution data on top of low-res
  sortI = sort(avgres, index.return=TRUE, decreasing=TRUE)$ix
  
  #plot data for each grid in geoQDat
  first=TRUE
  maxPerGrid = floor(maxToPlot/length(sortI))
  for(i in sortI) {
    #get data for this grid
    thisLon = geoQDat$lon[[i]]
    thisLat = geoQDat$lat[[i]]
    grid = make.surface.grid(list(lon=thisLon, lat=thisLat))
    dat = geoQDat$surfDisp[[i]]
    nx = length(geoQDat$lon[[i]])
    ny = length(geoQDat$lat[[i]])
    thisRangeLon = range(thisLon)
    thisRangeLat = range(thisLat)
    
    #randomly subsample maxPerGrid points from dat and grid
    if(is.finite(maxPerGrid) && maxPerGrid < length(dat)) {
      pts = sample(1:length(dat), maxPerGrid)
      grid = grid[pts,]
      dat = dat[pts]
    }
    
    #get points for plotting rectangle around this grid
    polyX = c(thisRangeLon[1], thisRangeLon[1], thisRangeLon[2], thisRangeLon[2], thisRangeLon[1])
    polyY = c(thisRangeLat[1], thisRangeLat[2], thisRangeLat[2], thisRangeLat[1], thisRangeLat[1])
    
    #plot water height (surface displacement) and black box around data range:
    if(first) {
      quilt.plot(grid, dat, xlim=rangeLon, ylim=rangeLat, zlim=rangeDat, nx=nx, ny=ny, 
                 main=main, xlab="Longitude", ylab="Latitude")
      first = FALSE
    }
    else {
      quilt.plot(grid, dat, add=TRUE, add.legend=FALSE, zlim=rangeDat, nx=nx, ny=nx)
    }
    polygon(polyX, polyY, border="black")
  }
  
}

#x is either a directory or a list of geoQDat objects.  If list of geoQ 
#objects, must also include allGeoQTimes
geoQMovie = function(x, allGeoQTimes=NULL, land=FALSE, range=NULL, maxToPlot=10000) {
  
  #x: either list of outputs from loadGeoQ or the directory with all the fort.* files
  #allGeoQTimes: ignored if x is a directory.  Vector of times corresponding to 
  #fort.q* files
  #range: if range is NULL, read in all data to get range.
  #Else read in files one by one
  #maxToPlot: maximum number of points to plot with quilt.plot
  
  #if x is a string, then it is a directory for loading geoQDat
  if(is.character(x)) {
    #change working directory
    wd = getwd()
    setwd(x)
    
    #get data files in given data direcotry
    qFiles = system("ls fort.q*", intern=TRUE)
    tFiles = system("ls fort.t*", intern=TRUE)
    
    #load in data files
    allGeoQDat = lapply(qFiles, loadGeoQ)
    allGeoQTimes = sapply(tFiles, loadGeoT)
    
    #change working directory back
    setwd(wd)
  }
  # if x is not a string, then geoQ and geoT values were input by user
  else {
    allGeoQDat = x
    if(is.null(allGeoQTimes)) {
      stop("no times for geoQ data input")
    }
  }
  
  # if specified, ignore surfDisp values over land
  if(!land) {
    for(i in 1:length(allGeoQDat)) {
      for(j in 1:length(allGeoQDat[[i]]$surfDisp)) {
        wh = allGeoQDat[[i]]$wh[[j]]
        allGeoQDat[[i]]$surfDisp[[j]][wh <= 0] = NA
      }
    }
  }
  
  #get range of GeoQ values
  geoQRange = function(geoQDat) {
    range(geoQDat$surfDisp, na.rm=TRUE, finite=TRUE)
  }
  if(is.null(range))
    datRange = range(sapply(allGeoQDat, geoQRange), na.rm=TRUE, finite=TRUE)
  else
    datRange = range
  
  for(i in 1:length(allGeoQTimes)) {
    geoQDat = allGeoQDat[[i]]
    plotGeoQ(geoQDat, zlim=datRange, 
             main=paste0("Water Surface Displacement (meters) at Time ", allGeoQTimes[i]))
  }
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
# Region plotting and data-loading tools

#Like loadGeoQ but for a specific region with latitude and longitude ranges
# given by latReg and lonReg respectively
getGeoQRegion = function(geoQDat, lonReg, latReg) {
  #NOTE: land == FALSE ==> only plot points with water depth > 0
  len = length(geoQDat$mx)
  
  #get data for each grid in geoQDat
  first=TRUE
  for(i in 1:len) {
    #get data for this grid
    thisLon = geoQDat$lon[[i]]
    thisLat = geoQDat$lat[[i]]
    grid = make.surface.grid(list(lon=thisLon, lat=thisLat))
    dat = geoQDat$surfDisp[[i]]
    nx = length(geoQDat$lon[[i]])
    ny = length(geoQDat$lat[[i]])
    thisRangeLon = range(thisLon)
    thisRangeLat = range(thisLat)
    
    #determine if this grid is contained within region grid.  If not, skip
    lonGood = (thisRangeLon[1] >= lonReg[1]) && (thisRangeLat[2] <= lonReg[2])
    latGood = (thisRangeLat[1] >= latReg[1]) && (thisRangeLat[2] <= latReg[2])
    if(!lonGood || !latGood)
      next
    
    #plot water height (surface displacement) and black box around data range:
    if(first) {
      quilt.plot(grid, dat, xlim=lonReg, ylim=latReg, zlim=rangeDat, nx=nx, ny=ny, 
                 main=main, xlab="Longitude", ylab="Latitude")
      first = FALSE
    }
    else {
      quilt.plot(grid, dat, add=TRUE, add.legend=FALSE, zlim=rangeDat, nx=nx, ny=nx)
    }
    polygon(polyX, polyY, border="black")
  }
  
  if(first) {
    stop("no data region is within specified lon and lat ranges")
  }
  
}

#Like plotGeoQ but for a specific region with latitude and longitude ranges
# given by latReg and lonReg respectively
plotRegion = function(geoQDat, lonReg, latReg, zlim=NULL, 
                      main="Water Surface Displacement (meters)", land = FALSE) {
  #NOTE: land == FALSE ==> only plot points with water depth > 0
  len = length(geoQDat$mx)
  
  calcXRes = function(i) {
    thismx = geoQDat$mx[[i]]
    thislon = geoQDat$lon[[i]]
    return((max(thislon)-min(thislon))/thismx)
  }
  calcYRes = function(i) {
    thismy = geoQDat$my[[i]]
    thislat = geoQDat$lat[[i]]
    return((max(thislat)-min(thislat))/thismy)
  }
  if(!land) {
    for(i in 1:length(geoQDat$surfDisp)) {
      wh = geoQDat$wh[[i]]
      geoQDat$surfDisp[[i]][wh <= 0] = NA
    }
  }
  
  #calculate range of data values (lat, lon, water height, resolution)
  rangeLon = range(sapply(geoQDat$lon, range, na.rm=TRUE))
  rangeLat = range(sapply(geoQDat$lat, range, na.rm=TRUE))
  if(is.null(zlim))
    rangeDat = range(sapply(geoQDat$surfDisp, range, na.rm=TRUE, finite=TRUE), 
                     na.rm=TRUE, finite=TRUE)
  else
    rangeDat = zlim
  xres = sapply(1:len, calcXRes)
  yres = sapply(1:len, calcYRes)
  avgres = (xres + yres)/2
  
  #plot low-resolution data first, high resolution data on top of low-res
  sortI = sort(avgres, index.return=TRUE, decreasing=TRUE)$ix
  
  #plot data for each grid in geoQDat
  first=TRUE
  for(i in sortI) {
    #get data for this grid
    thisLon = geoQDat$lon[[i]]
    thisLat = geoQDat$lat[[i]]
    grid = make.surface.grid(list(lon=thisLon, lat=thisLat))
    dat = geoQDat$surfDisp[[i]]
    nx = length(geoQDat$lon[[i]])
    ny = length(geoQDat$lat[[i]])
    thisRangeLon = range(thisLon)
    thisRangeLat = range(thisLat)
    
    #determine if this grid is contained within region grid.  If not, skip
    lonGood = (thisRangeLon[1] >= lonReg[1]) && (thisRangeLat[2] <= lonReg[2])
    latGood = (thisRangeLat[1] >= latReg[1]) && (thisRangeLat[2] <= latReg[2])
    next
    
    #get points for plotting rectangle around this grid
    polyX = c(thisRangeLon[1], thisRangeLon[1], thisRangeLon[2], thisRangeLon[2], thisRangeLon[1])
    polyY = c(thisRangeLat[1], thisRangeLat[2], thisRangeLat[2], thisRangeLat[1], thisRangeLat[1])
    
    #plot water height (surface displacement) and black box around data range:
    if(first) {
      quilt.plot(grid, dat, xlim=lonReg, ylim=latReg, zlim=rangeDat, nx=nx, ny=ny, 
                 main=main, xlab="Longitude", ylab="Latitude")
      first = FALSE
    }
    else {
      quilt.plot(grid, dat, add=TRUE, add.legend=FALSE, zlim=rangeDat, nx=nx, ny=nx)
    }
    polygon(polyX, polyY, border="black")
  }
  
  if(first) {
    stop("no data region is within specified lon and lat ranges")
  }
  
}