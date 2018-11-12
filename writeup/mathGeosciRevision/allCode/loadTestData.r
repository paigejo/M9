# setwd("~/git/M9/cascadia/_output/")
# setwd("~/git/M9/chile2010/_output/")

##### Main functions for users:

# loadFloodData()
#   return(list(allHMax=allHMax, topo=topo, lon=lon, lat=lat, Mw=Mw))

# loadAllDeformations()
#   return all deformation data

# calcLonLatDist(lon, lat, R = 3959)
#   return(list(X=X, Y=Y, distPerCellX=distPerCellX, distPerCellY=distPerCellY, 
#          grid=gridL))

# test emulator on small testing dataset from M9 project SageMathCloud page
library(RcppCNPy)
library(fields)
library(abind)

#####get SageMathCloud testing data:
loadFloodData = function() {
  wd = getwd()
  setwd("~/git/M9/test_data")
  
  #allHMax
  allHMax = array(NA, dim=c(5, 250, 250)) #remove the SS realization
  for(i in 0:4) {
    fname = paste0("allHMax", i, ".npy")
    tmp = npyLoad(fname)
    allHMax[i+1, , ] = tmp
  }
  Mw = c(9.0, 9.1, 9.0, 8.9, 9.0, 8.9, 8.7, 8.8, 8.7, 8.6, 8.6, 9.1, 9.2, 9.1, 9.1, 9.2, 9.1, 8.6, 8.6)
  finI = c(15, 12, 1, 4, 7)
  Mw = Mw[finI]
  
  #topography/bathymetry, lon, lat
  topo = npyLoad("CCTopo.npy")
  lon = npyLoad("CCLon.npy")
  lat = npyLoad("CCLat.npy")
  
  setwd(wd)
  
  return(list(allHMax=allHMax, topo=topo, lon=lon, lat=lat, Mw=Mw))
}
#CSZa: CSZR_XXL1.tt3
#CSZb: CSZR_XL1.tt3
#CSZc: CSZR_L1.tt3
#CSZd: CSZR_M1.tt3
#CSZe: CSZR_SM1.tt3
#CSZf: CSZ_SS3_Defm_FINAL.tt3
#allI = c(15, 12, 1, 4, 7, 19)
#finI = c(15, 12, 1, 4, 7)

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

#utility function
catStr = function(strs, sep=" ") {
  str = strs[1]
  if(length(strs) == 1)
    return(str)
  for(i in 2:length(strs)) {
    str = paste(str, strs[i], sep=sep)
  }
  return(str)
}

#Functions for loading topography data
loadTopo = function(fname) {
  dat = scan(fname, what="character")
  
  #get header information
  ncols = as.numeric(dat[1]) #num longitude vals
  nrows = as.numeric(dat[3]) #num latitude values
  xlower = as.numeric(dat[5])
  ylower = as.numeric(dat[7])
  cellsize = as.numeric(dat[9])
  nodata_value = as.numeric(dat[11])
  lon = seq(xlower, xlower+ncols*cellsize, l=ncols)
  lat = rev(seq(ylower, ylower+nrows*cellsize, l=nrows)) # reverse order because data starts from the top
  
  #get topography/bathmetry/deformation data
  dat = array(as.numeric(dat[13:length(dat)]), dim=c(ncols, nrows))
  
  #convert missing data to NAs
  dat[dat == nodata_value] = NA
  
  out = list(dat=dat, nlon=ncols,nlat=nrows,xlower=xlower,ylower=ylower,
             cellsize=cellsize,lon=lon,lat=lat)
  return(out)
}

plotTopo = function(topoDat, ...) {
  lon = topoDat$lon
  lat = topoDat$lat
  dat = topoDat$dat
  if(length(dim(drop(dat))) == 3) {
    dat = dat[,,2]
  }
  grid = make.surface.grid(list(lon=lon, lat=lat))
  otherArgs = list(...)
  if(is.null(otherArgs$main))
    otherArgs$main = "Topography/Bathymetry Data"
  do.call("quilt.plot", c(list(grid, dat), otherArgs))
  invisible(NULL)
}

#header:
#num ncols (longitude)
#num nrows (latitude
#num xlower (longitude)
#num ylower (latitude)
#num cellsize
#num nodata_value
#data
saveTopo = function(dat, lon, lat, fname) {
  # make sure lat and lon are in vector, not matrix form from the example .npy files
  if(is.matrix(lon)) {
    lon = lon[,1]
  }
  if(is.matrix(lat)) {
    lat = lat[1,]
  }
  
  #make sure rows of data go from high to low latitudes
  tmp = sort(lat, decreasing=TRUE, index.return=TRUE)
  lat = tmp$x
  dat = dat[,tmp$ix]
  
  #make sure lon goes from -180 to 180, not 0 to 360
  if(max(lon) > 180)
    lon[lon > 180] = lon[lon > 180] - 360
  
  # calculate header variable information
  ncols = length(lon)
  nrows = length(lat)
  xlower = min(lon)
  ylower = min(lat)
  cellsize = ((max(lon) - xlower)/ncols + (max(lat) - ylower)/nrows)/2 #sizes should be same in each direction
  nodata_value = -9999
  
  #add header info to output string (and make sure it's in scientific notation with good precision)
  outStr = paste0(ncols, " ncols")
  outStr = c(outStr, paste0(nrows, " nrows"))
  outStr = c(outStr, paste0(format(xlower, digits=15, scientific=TRUE), " xlower"))
  outStr = c(outStr, paste0(format(ylower, digits=15, scientific=TRUE), " ylower"))
  outStr = c(outStr, paste0(format(cellsize, digits=15, scientific=TRUE), " cellsize"))
  outStr = c(outStr, paste0(nodata_value, " nodata_value"))
  
  #set NA's in dat to nodata_value
  dat[is.na(dat)] = nodata_value
  
  #transpose dat because it's printed in column order.  Also convert to character vector
  for(c in 1:ncol(dat)) {
    outStr = c(outStr, catStr(format(dat[,c], digits=15, scientific=TRUE)))
  }
  
  print("starting to write file...")
  
  write(outStr, file=fname, sep="\n")
  
  print("finished writing file")
}

loadDTopo = function(fname) {
  dat = scan(fname, what="character")
  
  #get header information
  mx = as.numeric(dat[1])
  my = as.numeric(dat[3])
  mt = as.numeric(dat[5])
  xlower = as.numeric(dat[7])
  ylower = as.numeric(dat[9])
  t0 = as.numeric(dat[11])
  dx = as.numeric(dat[13])
  dy = as.numeric(dat[15])
  dt = as.numeric(dat[17])
  lon = seq(xlower, xlower+mx*dx, l=mx)
  lat = rev(seq(ylower, ylower+my*dy, l=my)) # reverse order because data starts from the top
  
  #get topography/bathmetry/deformation data
  dat = array(as.numeric(dat[19:length(dat)]), dim=c(mx, my, mt))
  
  out = list(dat=dat, mx=mx,my=my,mt=mt,xlower=xlower,ylower=ylower,
             t0=t0,dx=dx,dy=dy,dt=dt,lon=lon,lat=lat)
  return(out)
}

#currently supports only instantaneous deformation data (mt=1) of type 3 (.tt3)
#header:
#num mx
#num my
#num mt
#num xlower
#num ylower
#num t0
#num dx
#num dy
#num dt
#data
saveDTopo = function(dat, lon, lat, fname) {
  #make sure dat has no extra dimensions
  dat = drop(dat)
  
  # make sure lat and lon are in vector, not matrix form from the example .npy files
  if(is.matrix(lon)) {
    lon = lon[,1]
  }
  if(is.matrix(lat)) {
    lat = lat[1,]
  }
  
  #make sure rows of data go from high to low latitudes
  tmp = sort(lat, decreasing=TRUE, index.return=TRUE)
  lat = tmp$x
  dat = dat[,tmp$ix]
  
  #make sure lon goes from -180 to 180, not 0 to 360
  if(max(lon) > 180)
    lon[lon > 180] = lon[lon > 180] - 360
  
  # calculate header variable information
  mx = length(lon)
  my = length(lat)
  mt = 1
  xlower = min(lon)
  ylower = min(lat)
  t0 = 1
  dx = (max(lon) - xlower)/mx
  dy = (max(lat) - ylower)/my
  dt = 0
  
  #add header info to output string (and make sure it's in 
  # scientific notation when necessary with good precision)
  outStr = paste0(mx, " mx")
  outStr = c(outStr, paste0(my, " my"))
  outStr = c(outStr, paste0(mt, " mt"))
  outStr = c(outStr, paste0(format(xlower, digits=15, scientific=TRUE), " xlower"))
  outStr = c(outStr, paste0(format(ylower, digits=15, scientific=TRUE), " ylower"))
  outStr = c(outStr, paste0(format(t0, scientific=TRUE), " t0"))
  outStr = c(outStr, paste0(format(dx, digits=15, scientific=TRUE), " dx"))
  outStr = c(outStr, paste0(format(dy, digits=15, scientific=TRUE), " dy"))
  outStr = c(outStr, paste0(format(dt, scientific=TRUE), " dt"))
  
  #transpose dat because it's printed in column order.  Also convert to character vector
  for(c in 1:ncol(dat)) {
    outStr = c(outStr, catStr(format(dat[,c], digits=15, scientific=TRUE)))
  }
  
  print("starting to write file...")
  
  write(outStr, file=fname, sep="\n")
  
  print("finished writing file")
}

#CSZa: CSZR_XXL1.tt3
#CSZb: CSZR_XL1.tt3
#CSZc: CSZR_L1.tt3
#CSZd: CSZR_M1.tt3
#CSZe: CSZR_SM1.tt3
#CSZf: CSZ_SS3_Defm_FINAL.tt3

loadDeformations = function() {
  wd = getwd()
  setwd("~/git/M9/CSZR")
  
  #topography data files
  files = c("CSZR_XXL1.tt3", "CSZR_XL1.tt3", "CSZR_L1.tt3", 
            "CSZR_M1.tt3", "CSZR_SM1.tt3", "CSZ_SS3_Defm_FINAL.tt3")
  
  #NOTE: CSZ_SS3_Defm_FINAL.tt3 has slightly different dx, dy, dt.
  #Instead of 0.0166, 0.0166, and 0.5, it's 0.0167, 0.0167, and 1
  #respectively
  
  for(f in 1:length(files)) {
    if(f == 1) {
      out = loadDTopo(files[f])
      out$files = files[f]
      out$dat = out$dat[,,out$mt]
    }
    else {
      #concatenate data in out list:
      tmp = loadDTopo(files[f])
      out$files = c(out$files, files[f])
      out$dat = abind(out$dat, tmp$dat[,,tmp$mt], along=3)
      out$mx = c(out$mx, tmp$mx)
      out$my = c(out$my, tmp$my)
      out$mt = c(out$mt, tmp$mt)
      out$xlower = c(out$xlower, tmp$xlower)
      out$ylower = c(out$ylower, tmp$ylower)
      out$t0 = c(out$t0, tmp$t0)
      out$dx = c(out$dx, tmp$dx)
      out$dy = c(out$dy, tmp$dy)
      out$dt = c(out$dt, tmp$dt)
      out$lon = c(out$lon, tmp$lon)
      out$lat = c(out$lat, tmp$lat)
    }
  }
  
  #set wd back to what it was before
  setwd(wd)
  
  return(out)
}

loadAllDeformations = function() {
  wd = getwd()
  setwd("~/git/M9/CSZR")
  
  #topography data files
  files = system("ls *.tt3", intern=TRUE)
  
  #NOTE: CSZ_SS3_Defm_FINAL.tt3 has slightly different dx, dy, dt.
  #Instead of 0.0166, 0.0166, and 0.5, it's 0.0167, 0.0167, and 1
  #respectively
  
  for(f in 1:length(files)) {
    if(f == 1) {
      out = loadDTopo(files[f])
      out$files = files[f]
      out$dat = out$dat[,,out$mt]
    }
    else {
      #concatenate data in out list:
      tmp = loadDTopo(files[f])
      out$files = c(out$files, files[f])
      out$dat = abind(out$dat, tmp$dat[,,tmp$mt], along=3)
      out$mx = c(out$mx, tmp$mx)
      out$my = c(out$my, tmp$my)
      out$mt = c(out$mt, tmp$mt)
      out$xlower = c(out$xlower, tmp$xlower)
      out$ylower = c(out$ylower, tmp$ylower)
      out$t0 = c(out$t0, tmp$t0)
      out$dx = c(out$dx, tmp$dx)
      out$dy = c(out$dy, tmp$dy)
      out$dt = c(out$dt, tmp$dt)
      
      #just use the lon and lat of the first realization
      #out$lon = c(out$lon, tmp$lon)
      #out$lat = c(out$lat, tmp$lat)
    }
  }
  
  #set wd back to what it was before
  setwd(wd)
  
  return(out)
}

#header:
#num ncols (longitude)
#num nrows (latitude
#num xlower (longitude)
#num ylower (latitude)
#num cellsize
#num nodata_value
#data
saveASC = function(dat, lon, lat, fname) {
  # make sure lat and lon are in vector, not matrix form from the example .npy files
  if(is.matrix(lon)) {
    lon = lon[,1]
  }
  if(is.matrix(lat)) {
    lat = lat[1,]
  }
  
  #make sure rows of data go from high to low latitudes
  tmp = sort(lat, decreasing=TRUE, index.return=TRUE)
  lat = tmp$x
  dat = dat[,tmp$ix]
  
  #make sure lon goes from -180 to 180, not 0 to 360
  if(max(lon) > 180)
    lon[lon > 180] = lon[lon > 180] - 360
  
  # calculate header variable information
  ncols = length(lon)
  nrows = length(lat)
  xlower = min(lon)
  ylower = min(lat)
  cellsize = ((max(lon) - xlower)/ncols + (max(lat) - ylower)/nrows)/2 #sizes should be same in each direction
  nodata_value = -99999
  
  #add header info to output string (and make sure it's in scientific notation with good precision)
  outStr = paste0(ncols, " ncols")
  outStr = c(outStr, paste0(nrows, " nrows"))
  outStr = c(outStr, paste0(format(xlower, digits=15, scientific=TRUE), " xlower"))
  outStr = c(outStr, paste0(format(ylower, digits=15, scientific=TRUE), " ylower"))
  outStr = c(outStr, paste0(format(cellsize, digits=15, scientific=TRUE), " cellsize"))
  outStr = c(outStr, paste0(nodata_value, " nodata_value"))
  
  #set NA's in dat to nodata_value
  dat[is.na(dat)] = nodata_value
  
  #transpose dat because it's printed in column order.  Also convert to character vector
  for(c in 1:ncol(dat)) {
    outStr = c(outStr, catStr(format(dat[,c], digits=15, scientific=TRUE)))
  }
  
  print("starting to write file...")
  
  write(outStr, file=fname, sep="\n")
  
  print("finished writing file")
}

##### convert .asc topo files from grid extract website to format desired by GeoClaw
#setwd("~/git/M9/")
#NOTE: this function is not yet used?

convertTopo = function(fname) {
  dat = scan(fname, what="character")
  
  #get header information
  mx = as.numeric(dat[1])
  my = as.numeric(dat[3])
  mt = as.numeric(dat[5])
  xlower = as.numeric(dat[7])
  ylower = as.numeric(dat[9])
  t0 = as.numeric(dat[11])
  dx = as.numeric(dat[13])
  dy = as.numeric(dat[15])
  dt = as.numeric(dat[17])
  
  #get topography/bathmetry/deformation data
  dat = array(as.numeric(dat[19:length(dat)]), dim=c(mx, my, mt))
  
  out = list(dat=dat, mx=mx,my=my,mt=mt,xlower=xlower,ylower=ylower,
             t0=t0,dx=dx,dy=dy,dt=dt)
  return(out)
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

calcLonLatDist = function(lon, lat, R = 3959) {
  #R is radius of earth.  Default is radius in miles
  lonExtent = range(lon)
  latExtent = range(lat)
  
  if(is.null(lonExtent))
    lonExtent = c(235.79781, 235.82087)
  if(is.null(latExtent))
    latExtent = c(41.739671, 41.762726)
  CCLon = mean(lonExtent)
  CCLat = mean(latExtent)
  lonDist = cos(CCLat)*2*pi/360*R
  latDist = 2*pi*R/360
  
  X = lon*lonDist
  Y = lat*latDist
  gridL = matrix(c(X, Y), ncol=2)
  
  nX = length(lon)
  nY = length(lat)
  maxX = max(X)
  minX = min(X)
  maxY = max(Y)
  minY = min(Y)
  distPerCellX = (maxX - minX)/nX
  distPerCellY =  (maxY - minY)/nY
  
  out = list(X=X, Y=Y, distPerCellX=distPerCellX, distPerCellY=distPerCellY, 
             grid=gridL)
}


