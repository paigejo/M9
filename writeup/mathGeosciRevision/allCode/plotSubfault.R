##### plot original CSZ fault geometry data
library(graphics)

# plot all subfaults in the fault (users should use this instead of plotSubfault for simplicity)
# set plotData=FALSE to only plot the polygons in the fault with no filling
# NOTE: don't worry about the warnings, they are caused by setting par to be suitable for the legend.
plotFault = function(rows, plotVar="depth", varRange=NULL, 
                     cols=tim.colors(), new=TRUE, lwd=1, plotData=TRUE, 
                     legend.mar=6, logScale=FALSE, ...) {
  if(!is.data.frame(rows))
    rows = data.frame(rows)
  
  # if log scale, modify data and range accordingly
  if(logScale) {
    plotVar = log10(plotVar)
    if(!is.null(varRange))
      varRange=log10(varRange)
  }
  
  # if the user supplies a variable to plot in plotVar:
  if(!is.character(plotVar)) {
    rows$tmp = plotVar
    plotVar = "tmp"
  }
  
  # set default plotting parameters to be reasonable
  if(is.null(varRange))
    varRange = range(rows[,plotVar], na.rm=TRUE)
  otherArgs = list(...)
  if(is.null(otherArgs$xlim))
    otherArgs$xlim=range(rows$longitude)
  if(is.null(otherArgs$ylim))
    otherArgs$ylim=range(rows$latitude)
  if(!plotData)
    cols = c(NA, NA)
  else {
    currPar = par()
    newPar = currPar
    newMar = newPar$mar
    newMar[4] = max(newMar[4], legend.mar)
    newPar$mar = newMar
    if(currPar$mar[4] != newMar[4])
      suppressWarnings({par(newPar)})
  }
  
  # plot the subfaults
  do.call("plotSubfault", c(list(rows[1,], plotVar, varRange, cols, new, polyArgs=list(lwd=lwd)), otherArgs))
  invisible(apply(rows[-1,], 1, plotSubfault, plotVar=plotVar, varRange=varRange, cols=cols, new=FALSE, polyArgs=list(lwd=lwd)))
  
  #add legend if necessary
  if(plotData) {
    if(logScale) {
      ticks = axisTicks(varRange, log=TRUE)
      image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
                 col=cols, add = TRUE, axis.args=list(at=log10(ticks), labels=ticks))
    }
    else {
      image.plot(zlim=varRange, nlevel=length(cols), legend.only=TRUE, horizontal=FALSE, 
                 col=cols, add = TRUE)
    }
    
    # suppressWarnings({par(currPar)})
  }
}

plotSubfault = function(row, plotVar="depth", varRange=NULL, 
                        cols=tim.colors(), new=FALSE, polyArgs=NULL, ...) {
  
  if(!is.list(row))
    row = as.list(row)
  if(is.null(varRange))
    varRange = range(faultGeom[,plotVar])
  plotVar = row[[plotVar]]
  
  otherArgs = list(...)
  if(is.null(otherArgs$xlim))
    otherArgs$xlim=range(faultGeom$longitude)
  if(is.null(otherArgs$ylim))
    otherArgs$ylim=range(faultGeom$latitude)
  if(is.null(otherArgs$ylab))
    otherArgs$ylab = "Latitude"
  if(is.null(otherArgs$xlab))
    otherArgs$xlab = "Longitude"
  otherArgs$type="n"
  
  geom = calcGeom(row)
  coords = geom$corners[,1:2]
  
  # generate new plot if necessary
  if(new)
    do.call("plot", c(list(1, 2), otherArgs))
  
  # get color to plot
  vals = c(varRange, plotVar)
  vals = vals-vals[1]
  vals = vals/(vals[2] - vals[1])
  col = cols[round(vals[3]*(length(cols)-1))+1]
  
  # plot subfault rectangle
  ord = c(1:4, 1)
  do.call("polygon", c(list(coords[ord,1], coords[ord,2], col=col), polyArgs))
  
}

# This function is based on the SubFault.calculate_geometry function from
# dtopotools.  Note that this is written based on coordinate specification 
# being "top center".  The Okada model requires depth at the bottom center.
calcGeom = function(subfault) {
  #   Calculate the fault geometry.
  #   
  #   Routine calculates the class attributes *corners* and 
  #   *centers* which are the corners of the fault plane and 
  #   points along the centerline respecitvely in 3D space.
  #   
  #   **Note:** *self.coordinate_specification*  specifies the location on each
  #   subfault that corresponds to the (longitude,latitude) and depth 
  #   of the subfault.
  #   Currently must be one of these strings:
  #     
  #   - "bottom center": (longitude,latitude) and depth at bottom center
  #   - "top center": (longitude,latitude) and depth at top center
  #   - "centroid": (longitude,latitude) and depth at centroid of plane
  #   - "noaa sift": (longitude,latitude) at bottom center, depth at top,  
  #   This mixed convention is used by the NOAA SIFT
  #   database and "unit sources", see:
  #     http://nctr.pmel.noaa.gov/propagation-database.html
  #   
  #   The Okada model is expressed assuming (longitude,latitude) and depth
  #   are at the bottom center of the fault plane, so values must be
  #   shifted or other specifications.
  
  row = subfault
  if(!is.list(row))
    row = data.frame(row)
  
  # Simple conversion factors
  #lat2meter = util.dist_latlong2meters(0.0, 1.0)[1]
  #LAT2METER = 110.574 #* 10^3
  LAT2METER = 111133.84012073894 #/10^3
  lat2meter = LAT2METER
  DEG2RAD = 2*pi/360
  
#   Setup coordinate arrays
#   Python format:
#   Top edge    Bottom edge
#     a ----------- b          ^ 
#     |             |          |         ^
#     |             |          |         |
#     |             |          |         | along-strike direction
#     |             |          |         |
#     0------1------2          | length  |
#     |             |          |
#     |             |          |
#     |             |          |
#     |             |          |
#     d ----------- c          v
#     <------------->
#     width
#   
#     <-- up dip direction
  
#   corners = [[None, None, None], # a 
#              [None, None, None], # b
#              [None, None, None], # c
#              [None, None, None]] # d
#   centers = [[None, None, None], # 1
#              [None, None, None], # 2 
#              [None, None, None]] # 3
  corners = matrix(nrow=4, ncol=3) # rows are a,b,c,d, and cols are lon,lat,depth
  centers = matrix(nrow=3, ncol=3) # rows are 1,2,3 (0,1,2), and cols are lon,lat,depth
  
  # Set depths
  centers[1,3] = row$depth
  centers[2,3] = row$depth + 0.5 * row$width * sin(row$dip * DEG2RAD)
  centers[3,3] = row$depth + row$width * sin(row$dip * DEG2RAD)
  
  corners[1,3] = centers[1,3]
  corners[4,3] = centers[1,3]
  corners[2,3] = centers[3,3]
  corners[3,3] = centers[3,3]
  
  # Locate fault plane in 3D space:  
  # See the class docstring for a guide to labeling of corners/centers.
  
  # Vector *up_dip* goes from bottom edge to top edge, in meters,
  # from point 2 to point 0 in the figure in the class docstring.
  
  up_dip = c(-row$width * cos(row$dip * DEG2RAD) * cos(row$strike * DEG2RAD) 
            / (LAT2METER * cos(row$latitude * DEG2RAD)),
            row$width * cos(row$dip * DEG2RAD) 
            * sin(row$strike * DEG2RAD) / LAT2METER)
  
  centers[1,1:2] = c(row$longitude, row$latitude)
  centers[2,1:2] = c(row$longitude - 0.5 * up_dip[1],
                     row$latitude - 0.5 * up_dip[2])
  centers[3,1:2] = c(row$longitude - up_dip[1],
                     row$latitude - up_dip[2])
  
  # Calculate coordinates of corners:
  # Vector *strike* goes along the top edge from point 0 to point a
  # in the figure in the class docstring.
  
  up_strike = c(0.5 * row$length * sin(row$strike * DEG2RAD) 
                / (lat2meter * cos(centers[3,2] * DEG2RAD)),
                0.5 * row$length * cos(row$strike * DEG2RAD) / lat2meter)
  
  corners[1,1:2] = c(centers[1,1] + up_strike[1],
                         centers[1,2] + up_strike[2])
  corners[2,1:2] = c(centers[3,1] + up_strike[1], 
                         centers[3,2] + up_strike[2])
  corners[3,1:2] = c(centers[3,1] - up_strike[1],
                         centers[3,2] - up_strike[2])
  corners[4,1:2] = c(centers[1,1] - up_strike[1],
                         centers[1,2] - up_strike[2])
  
  return(list(corners=corners, centers=centers))
}

plotSubfaultGeom = function(subfault, varRange=NULL, 
                            cols=tim.colors(), new=FALSE, polyArgs=NULL, ...) {
  plotVar = "depth"
  geom = calcGeom(subfault)
  coords = geom$corners
  depth = geom$centers[2,3]
  
  # set varRange if necessary
  if(is.null(varRange))
    varRange = range(faultGeom[,"depth"])
  
  # generate new plot if necessary
  if(new)
    do.call("plot.new", list(...))
  
  # get color to plot
  vals = c(varRange, depth)
  vals = vals-vals[1]
  vals = vals/(vals[2] - vals[1])
  col = cols[round(vals[3]*(length(cols)-1))+1]
  
  # plot subfault rectangle
  ord = c(1, 2, 3, 4, 1)
  do.call("polygon", c(list(x=coords[ord,1], y=coords[ord,2], col=col), polyArgs))
}

##### function for subdividing faults and subfaults.  Assumes units are in meters
divideFault = function(rows, nDown=3, nStrike=4) {
  subFaults = apply(rows, 1, divideSubfault, nDown=nDown, nStrike=nStrike)
  do.call("rbind", subFaults)
}

# This function subdivides the subfault into smaller subfaults, with nDown in the dip
# direction and nStrike in the strike direction.  Assumes the coordinates of the 
# subfaults are at the top center of each subfault (as in the CSZ fault geometry),
# the locations are in latitude and longitude, and the distance units are in meters
divideSubfault = function(row, nDown=3, nStrike=4) {
  if(!is.list(row))
    row = as.list(row)
  
  geom = calcGeom(row)
  corners = geom$corners
  
  # generate set of nDown points mid strike and going down dip with the correct depths
  centers = matrix(nrow=nDown, ncol=3) # rows are 1,2,...,nDown, and cols are lon,lat,depth
  
  # Simple conversion factors
  #lat2meter = util.dist_latlong2meters(0.0, 1.0)[1]
  #LAT2METER = 110.574 #* 10^3
  LAT2METER = 111133.84012073894 #/10^3
  lat2meter = LAT2METER
  DEG2RAD = 2*pi/360
  
  # Set depths
  centers[,3] = row$depth + (0:(nDown-1))/nDown * row$width * sin(row$dip * DEG2RAD)
  
  # Vector *up_dip* goes from bottom edge to top edge, in meters,
  # from point 2 to point 0 in the figure in the class docstring.
  # Vector *up_strike* goes along the top edge from point d to point a
  # in the figure in the class docstring. (this is different from in calcGeom)
  # up_depth is the depth difference from top to bottom of fault
  up_dip = c(-row$width * cos(row$dip * DEG2RAD) * cos(row$strike * DEG2RAD) 
             / (LAT2METER * cos(row$latitude * DEG2RAD)),
             row$width * cos(row$dip * DEG2RAD) 
             * sin(row$strike * DEG2RAD) / LAT2METER)
  up_strike = c(row$length * sin(row$strike * DEG2RAD) 
                / (lat2meter * cos(geom$centers[3,2] * DEG2RAD)),
                row$length * cos(row$strike * DEG2RAD) / lat2meter)
  
  # Set lon and lat of centers
  centers[,1] = row$longitude - (0:(nDown-1))/nDown * up_dip[1]
  centers[,2] = row$latitude - (0:(nDown-1))/nDown * up_dip[2]
  
  # get points along down strike edge of subfault by moving in down strike direction 
  # from centers
  starts = centers
  starts[,1:2] = starts[,1:2] + (- 0.5 + 1/(2*nStrike))*matrix(rep(up_strike, nrow(centers)), ncol=2, byrow = TRUE)
  
  # compute complete set of subfault coordinates by shifting starts by various vectors in
  # the up strike direction
  coords = cbind(rep(starts[,1], nStrike), rep(starts[,2], nStrike), rep(starts[,3], nStrike))
  shifts = matrix(rep(up_strike, nStrike*nDown), ncol=2, byrow = TRUE)
  mult = rep((0:(nStrike-1))/nStrike, rep(nDown, nStrike))
  shifts = sweep(shifts, MARGIN=1, STATS = mult, FUN = "*")
  coords[,1:2] = coords[,1:2] + shifts
  
  # construct subfaults
  rows = data.frame(matrix(unlist(rep(row, nStrike*nDown)), ncol=length(row), byrow = TRUE))
  names(rows) = names(row)
  rows$Fault = row$Fault
  rows$length = rows$length/nStrike
  rows$width = rows$width/nDown
  rows$depth = coords[,3]
  rows$longitude = coords[,1]
  rows$latitude = coords[,2]
  
  return(rows)
}

# get table with lon, lat, and depth coords for each subfault
# when index=2, this gives the center of the subfaults.  index = 1
# corresponds to the top-center, and index=3 corresponds to the 
# bottom center
getFaultCenters = function(rows, index=2) {
  t(apply(rows, 1, getSubfaultCenter, index=index))
}

getSubfaultCenter = function(row, index=2) {
  if(!is.list(row))
    row = as.list(row)
  
  geom = calcGeom(row)
  centers = geom$centers
  return(centers[index,])
}

getFaultPolygons = function(fault=csz) {
  if(length(unique(fault$Fault)) != length(fault$Fault)) {
    fault$Fault = 1:nrow(fault)
  }
  
  calcSubfaultPolygon = function(subfault) {
    tmp = subfault
    if(!is.list(subfault)) {
      tmp = matrix(subfault, ncol=length(subfault))
      colnames(tmp) = names(subfault)
      tmp = data.frame(tmp)
    }
    subfault = tmp
    
    subFaultPoly = calcGeom(subfault)$corners[,1:2]
    return(cbind(subfault$Fault, subFaultPoly))
  }
  
  # faultPolys = apply(fault, 1, calcSubfaultPolygon)
  subfaultList = list()
  for(i in 1:nrow(fault)) {
    subfaultList = c(subfaultList, list(fault[i,]))
  }
  faultPolys = lapply(subfaultList, calcSubfaultPolygon)
  faultPolys = do.call("rbind", faultPolys)
  faultPolys = data.frame(list(Fault=faultPolys[,1], longitude=faultPolys[,2], latitude=faultPolys[,3]))
  return(faultPolys)
}

##### now make ggplot fault plotters

# plot all subfaults in the fault using ggplot.  Don't add data, returns the geom_polygon object
ggPlotFault = function(fault=csz, color="black") {
  faultPolys = getFaultPolygons(fault)
  geom_polygon(aes(x=longitude, y=latitude, group=factor(Fault)), data=faultPolys, 
               fill=rgb(1,1,1, 0), color=color)
}

ggPlotFaultDat = function(rows, plotVar="depth", varRange=NULL, plotData=TRUE, 
                          logScale=FALSE, xlim=c(-128, -122), ylim=c(39.5, 50.5), 
                          xlab="Longitude", ylab="Latitude", main="Cascadia Subduction Zone", 
                          clab="Depth (m)", addLegend=TRUE, lwd=1) {
  
  # get relevant map data
  states <- map_data("state")
  west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
  canada = map_data("world", "Canada")
  
  if(!is.data.frame(rows))
    rows = data.frame(rows)
  
  # rename row$Fault so that each fault gets unique name
  rows$Fault=1:nrow(rows)
  
  # if the user supplies a variable to plot in plotVar:
  if(!is.character(plotVar)) {
    rows$tmp = plotVar
    plotVar = "tmp"
  }
  
  # make fault polygons
  faultDat = getFaultPolygons(rows)
  faultDat = merge(faultDat, rows[,c("Fault", plotVar)], by=c("Fault"))
  faultDat$plotVar=faultDat[,plotVar]
  
  # grey maps plot:
  bg = ggplot(faultDat, aes(x=longitude, y=latitude)) + 
    geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
    geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
    coord_fixed(xlim = xlim,  ylim = ylim, ratio = 1.3, expand=FALSE) + 
    labs(x=xlab, y=ylab) + scale_x_continuous("Longitude", c(-127, -125, -123), labels=c("-127", "", "-123"), limits=c(-360, 360)) + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill='lightblue1'))
  
  # generate fault polygons portion of plot
  
  if(!plotData) {
    # faultPoly = ggPlotFault(rows, color=rgb(.2,.2,.2))
    pl = bg + geom_polygon(aes(fill=rgb(1,1,1,0), group=factor(Fault)), color="black", size=lwd)
  }
  else {
    faultPoly = geom_polygon(aes(fill=plotVar, group=factor(Fault)), color="black", size=lwd)
    
    if(is.null(varRange)) {
      if(logScale)
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1, trans="log") + 
          ggtitle(main)
      else
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1) +
          ggtitle(main)
    }
    else {
      if(logScale)
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1, trans="log", 
                                                   limits=varRange) + 
          ggtitle(main)
      else
        pl = bg + faultPoly + scale_fill_distiller(clab, palette = "Spectral", direction=-1, 
                                                   limits=varRange) +
          ggtitle(main)
    }
    # ii <- cut(values, breaks = seq(min(values), max(values), len = 100), 
    #           include.lowest = TRUE)
    # ## Use bin indices, ii, to select color from vector of n-1 equally spaced colors
    # colors <- colorRampPalette(c("lightblue", "blue"))(99)[ii]
    # if(logScale)
    #   faultPoly = faultPoly + scale_fill_manual("", palette = "Spectral", direction=-1, trans="log")
    # else
    #   faultPoly = faultPoly + scale_fill_distiller("", palette = "Spectral", direction=-1)
  }
  
  if(addLegend)
    pl + guides(color=FALSE)
  else
    pl + guides(color=FALSE, fill=FALSE)
}





