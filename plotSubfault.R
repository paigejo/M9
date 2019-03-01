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

# same as calcGeom, except assumes the subfault has information about the coordinates of it's corners
calcGeom2 = function(subfault) {
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
  
  # up_dip = c(-row$width * cos(row$dip * DEG2RAD) * cos(row$strike * DEG2RAD) 
  #            / (LAT2METER * cos(row$latitude * DEG2RAD)),
  #            row$width * cos(row$dip * DEG2RAD) 
  #            * sin(row$strike * DEG2RAD) / LAT2METER)
  up_dip = c(row$bottomLeftX - row$bottomRightX, 0)
  
  centers[1,1:2] = c(0.5 * (row$topLeftX + row$bottomLeftX), 0.5 * (row$topLeftY + row$bottomLeftY))
  centers[2,1:2] = c(centers[1,1] - 0.5 * up_dip[1],
                     centers[1,2] - 0.5 * up_dip[2])
  centers[3,1:2] = c(centers[1,1] - up_dip[1],
                     centers[1,2] - up_dip[2])
  
  # Calculate coordinates of corners:
  # Vector *strike* goes along the top edge from point 0 to point a
  # in the figure in the class docstring.
  
  # up_strike = c(0.5 * row$length * sin(row$strike * DEG2RAD) 
  #               / (lat2meter * cos(centers[3,2] * DEG2RAD)),
  #               0.5 * row$length * cos(row$strike * DEG2RAD) / lat2meter)
  up_strike = c(0.5 * (row$topLeftX - row$bottomLeftX),
                0.5 * (row$topLeftY - row$bottomLeftY))
  
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

# same as dividFault, except takes in the polygons of the fault instead of 
# the original geometry format of, e.g., faultGeom
divideFault2 = function(rows, nDown=3, nStrike=4) {
  subFaults = apply(rows, 1, divideSubfault2, nDown=nDown, nStrike=nStrike)
  do.call("rbind", subFaults)
}

# same as dividSubfault, except takes in the polygons of the fault instead of 
# the original geometry format of, e.g., faultGeom
divideSubfault2 = function(row, nDown=3, nStrike=4) {
  if(!is.list(row))
    row = as.list(row)
  
  geom = calcGeom2(row)
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
  # up_dip = c(-row$width * cos(row$dip * DEG2RAD) * cos(row$strike * DEG2RAD) 
  #            / (LAT2METER * cos(row$latitude * DEG2RAD)),
  #            row$width * cos(row$dip * DEG2RAD) 
  #            * sin(row$strike * DEG2RAD) / LAT2METER)
  up_dip = c(row$topLeftX - row$topRightX,
             row$topLeftY - row$topRightY)
  # up_strike = c(row$length * sin(row$strike * DEG2RAD) 
  #               / (lat2meter * cos(geom$centers[3,2] * DEG2RAD)),
  #               row$length * cos(row$strike * DEG2RAD) / lat2meter)
  up_strike = c(row$topLeftX - row$bottomLeftX,
                row$topLeftY - row$bottomLeftY)
  
  # Set lon and lat of centers
  # centers[,1] = row$longitude - (0:(nDown-1))/nDown * up_dip[1]
  # centers[,2] = row$latitude - (0:(nDown-1))/nDown * up_dip[2]
  centers[,1] = row$topLeftX - (1:(nDown))/nDown * up_dip[1]
  centers[,2] = 0.5 * (row$topLeftY + row$bottomLeftY) - (0:(nDown-1))/nDown * up_dip[2]
  
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
  
  # construct subfaults (ignore latitude longitude coordinates)
  rows = data.frame(matrix(unlist(rep(row, nStrike*nDown)), ncol=length(row), byrow = TRUE))
  names(rows) = names(row)
  rows$Fault = row$Fault
  rows$length = rows$length/nStrike
  rows$width = rows$width/nDown
  rows$depth = coords[,3]
  # rows$longitude = coords[,1]
  # rows$latitude = coords[,2]
  rows$middleRightX = coords[,1]
  rows$middleRightY = coords[,2]
  
  # calculate the polygons of the sub faults (divide length by two since we start at the center of the sub fault)
  subLength = 0.5 * (row$topRightY - row$bottomRightY) / nStrike
  subWidth = (row$topRightX - row$topLeftX) / nDown
  
  rows$topRightX = rows$middleRightX
  rows$bottomRightX = rows$middleRightX
  rows$topLeftX = rows$middleRightX - subWidth
  rows$bottomLeftX = rows$middleRightX - subWidth
  rows$bottomMiddleX = rows$middleRightX - 0.5 * subWidth
  rows$topMiddleX = rows$middleRightX - 0.5 * subWidth
  
  rows$topRightY = rows$middleRightY + subLength
  rows$bottomRightY = rows$middleRightY - subLength
  rows$topLeftY = rows$middleRightY + subLength
  rows$bottomLeftY = rows$middleRightY - subLength
  rows$bottomMiddleY = rows$middleRightY - subLength
  rows$topMiddleY = rows$middleRightY + subLength
  
  return(rows)
}

# compute_subfault_distances(fault):
# modified from https://github.com/rjleveque/KLslip-paper/blob/master/KL2d_figures.py
# Estimate the distance between subfaults i and j for every pair in the data frame 
# fault. Distances are calculated in kilometers rather than in meters as in the python script.
# :Inputs:
# -  *fault* of class dtopotools.Fault or some subclass,
# 
# :Outputs:
# - *D* array of Euclidean distances based on longitudes, latitudes, and depths
# - *Dstrike* array of estimated distances along strike direction
# - *Ddip* array of estimated distances along dip direction
# with D**2 = Dstrike**2 + Ddip**2 to within roundoff.
# For each array, the [i,j] entry is distance from subfault i to j when
# ordered in the order the subfaults appear in the list fault.subfaults.
# Distance in dip direction based on differences in depth.  
compute_subfault_distances = function(fault) {
  rad = pi/180.       # conversion factor from degrees to radians
  rr = 6.378e6        # radius of earth
  lat2meter = rr*rad  # conversion factor from degrees latitude to meters
  
  nsubfaults = nrow(fault)
  D = matrix(nrow=nsubfaults, ncol=nsubfaults)
  Dstrike2 = matrix(nrow=nsubfaults, ncol=nsubfaults)
  Ddip = matrix(nrow=nsubfaults, ncol=nsubfaults)
  for(i in 1:nsubfaults) {
    xi = fault$longitude[i]
    yi = fault$latitude[i]
    zi = fault$depth[i]
    
    for(j in 1:nsubfaults) {
      xj = fault$longitude[j]
      yj = fault$latitude[j]
      zj = fault$depth[j]
      dx = abs(xi-xj)*cos(0.5*(yi+yj)*pi/180.) * lat2meter
      dy = abs(yi-yj) * lat2meter
      dz = abs(zi-zj)
      
      # Euclidean distance:
      D[i,j] = sqrt(dx**2 + dy**2 + dz**2)
      
      # estimate distance down-dip based on depths:
      dip = 0.5*(fault$dip[i] + fault$dip[j])
      ddip1 = dz / sin(dip*pi/180.)
      Ddip[i,j] = ddip1 
      
      # compute distance in strike direction to sum up properly:
      Dstrike2[i,j] = D[i,j]**2 - Ddip[i,j]**2
    }
  }
  Dstrike2[Dstrike2 < 0] = 0
  
  # return the results, converted from meters to kilometers
  list(D=D / 1000,Dstrike=sqrt(Dstrike2) / 1000,Ddip=Ddip / 1000)
}

# same as compute_subfault_distances, but for the gps data. In order to 
# estimate the distance down dip for each of these points, we project them 
# onto the fall geometry
compute_dip_strike_distance_gps = function(gpsDat, fault) {
  rad = pi/180.       # conversion factor from degrees to radians
  rr = 6.378e6        # radius of earth
  lat2meter = rr*rad  # conversion factor from degrees latitude to meters
  
  # calculate gps dips
  faultIndices = getFaultIndexFromGps(gpsDat, fault)
  gpsDips = fault$dip[faultIndices]
  
  npoints = nrow(gpsDat)
  D = matrix(nrow=npoints, ncol=npoints)
  Dstrike2 = matrix(nrow=npoints, ncol=npoints)
  Ddip = matrix(nrow=npoints, ncol=npoints)
  for(i in 1:npoints) {
    xi = gpsDat$lon[i]
    yi = gpsDat$lat[i]
    zi = gpsDat$Depth[i]
    
    for(j in 1:npoints) {
      xj = gpsDat$lon[j]
      yj = gpsDat$lat[j]
      zj = gpsDat$Depth[j]
      dx = abs(xi-xj)*cos(0.5*(yi+yj)*pi/180.) * lat2meter
      dy = abs(yi-yj) * lat2meter
      dz = abs(zi-zj)
      
      # Euclidean distance:
      D[i,j] = sqrt(dx**2 + dy**2 + dz**2)
      
      # estimate distance down-dip based on depths:
      dip = 0.5*(gpsDips[i] + gpsDips[j])
      ddip1 = dz / sin(dip*pi/180.)
      Ddip[i,j] = ddip1 
      
      # compute distance in strike direction to sum up properly:
      Dstrike2[i,j] = D[i,j]**2 - Ddip[i,j]**2
    }
  }
  Dstrike2[Dstrike2 < 0] = 0
  
  # return the results, converted from meters to kilometers
  list(D=D / 1000,Dstrike=sqrt(Dstrike2) / 1000,Ddip=Ddip / 1000)
}

# for any given point over the given fault geometry, return the index 
# of the row of default data table corresponding to the part of the fault 
# that the point is over.
getFaultIndexFromGps = function(gpsDat=slipDatCSZ, faultGeom=csz) {
  # helper function for determining if GPS observations are within a specific subfault geometry
  getSubfaultGPSDat = function(i) {
    row = faultGeom[i,]
    geom = calcGeom(row)
    corners = geom$corners[,1:2]
    in.poly(cbind(gpsDat$lon, gpsDat$lat), corners)
  }
  gpsInSubFault = sapply(1:nrow(faultGeom), getSubfaultGPSDat)
  apply(gpsInSubFault, 1, function(x) {match(TRUE, x)})
}

# get table with lon, lat, and depth coords for each subfault
# when type == "centers":
# when index=2, this gives the center of the subfaults.  index = 1
# corresponds to the top-center, and index=3 corresponds to the 
# bottom center
# when type == "corners":
# index == 1, ..., 4 corresponds to the top left corner to the bottom left corner 
# in the clockwise direction
getFaultCenters = function(rows, index=2, type=c("centers", "corners")) {
  type = match.arg(type)
  t(apply(rows, 1, getSubfaultCenter, index=index, type=type))
}

getSubfaultCenter = function(row, index=2, type=c("centers", "corners")) {
  type = match.arg(type)
  if(!is.list(row))
    row = as.list(row)
  
  geom = calcGeom(row)
  pts = geom[[type]]
  return(pts[index,])
}

getFaultPolygons = function(fault=csz, faultCornerTable=NULL) {
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
    
    if(is.null(faultCornerTable))
      subFaultPoly = calcGeom(subfault)$corners[,1:2]
    else
      subFaultPoly = rbind(c(subfault$topLeftX, subfault$topLeftY), 
                           c(subfault$topRightX, subfault$topRightY), 
                           c(subfault$bottomRightX, subfault$bottomRightY), 
                           c(subfault$bottomLeftX, subfault$bottomLeftY))
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
                          xlab=NULL, ylab=NULL, main="Cascadia Subduction Zone", 
                          clab="Depth (m)", addLegend=plotData, lwd=1, 
                          xName="longitude", yName="latitude", 
                          projection=NULL, parameters=NULL, orientation=NULL, 
                          scale=1, roundDigits=2, coordsAlreadyScaled=FALSE) {
  
  if(is.null(orientation) && !is.null(projection)) {
    warning("no orientation specified, so projection being set to the last used projection")
    projection = ""
    parameters = NULL
  }
  
  # get relevant map data
  if(!is.null(projection)) {
    states <- map_data("state", projection=projection, parameters=parameters, orientation=orientation)
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada", projection=projection, parameters=parameters, orientation=orientation)
  } else {
    states <- map_data("state")
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada")
  }
  
  # set some reasonable plotting parameters
  if(is.null(projection)) {
    xlab = "Longitude"
    ylab = "Latitude"
  } else {
    xlab = "Easting (km)"
    ylab = "Northing (km)"
  }
  
  if(!is.data.frame(rows))
    rows = data.frame(rows)
  
  # rename row$Fault so that each fault gets unique name
  rows$Fault=1:nrow(rows)
  
  # if the user supplies a variable to plot in plotVar:
  if(!is.character(plotVar)) {
    rows$tmp = plotVar
    plotVar = "tmp"
  }
  
  # 
  if("topRightX" %in% names(rows))
    faultCornerTable = rows[,9:ncol(rows)]
  else
    faultCornerTable = NULL
  
  # make fault polygons
  faultDat = getFaultPolygons(rows, faultCornerTable=faultCornerTable)
  faultDat = merge(faultDat, rows[,c("Fault", plotVar)], by=c("Fault"))
  faultDat$plotVar=faultDat[,plotVar]
  
  # rescale coordinates of the faulty geometry if necessary
  if(coordsAlreadyScaled) {
    faultDat$longitude = faultDat$longitude / scale
    faultDat$latitude = faultDat$latitude / scale
  }
  
  # set x and y limits if necessary
  if(is.null(xlim) && is.null(ylim)) {
    xlim=c(-128, -122)
    ylim=c(39.5, 50.5)
    proj<- mapproject(xlim, ylim) # if projection unspecified, last projection is used
    xlim = proj$x
    ylim = proj$y
    xScale = scale_x_continuous(xlab, c(-.02, 0, .02, .04), labels=as.character(round(scale*c(-.02, 0, .02, .04), digits=roundDigits)), limits=c(-360, 360))
    yScale = scale_y_continuous(ylab, seq(-.74, -.60, by=.02), labels=as.character(round(scale*seq(-.74, -.60, by=.02), digits=roundDigits)), limits=c(-360, 360))
  } else if(is.null(xlim)) {
    xlim = range(c(rows[[xName]], faultDat$longitude))
  } else if(is.null(ylim)) {
    ylim = range(c(rows[[yName]], faultDat$latitude))
  } else {
    xScale = scale_x_continuous(xlab, c(-127, -125, -123), labels=c("-127", "", "-123"), limits=c(-360, 360))
  }
  
  # grey maps plot:
  if(identical(projection, "")) {
    bg = ggplot(faultDat, aes(x=longitude, y=latitude)) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + yScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  } else {
    bg = ggplot(faultDat, aes(x=longitude, y=latitude)) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, ratio = 1.3, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  }
  
  # generate fault polygons portion of plot
  if(!plotData) {
    # faultPoly = ggPlotFault(rows, color=rgb(.2,.2,.2))
    pl = bg + geom_polygon(aes(group=factor(Fault)), color="black", size=lwd, fill=NA)
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

##### now do the same but for the gps data
ggplotGpsData = function(gpsDat, plotVar=gpsDat$Depth, varRange=NULL, plotData=TRUE, 
                         logScale=FALSE, xlim=NULL, ylim=NULL, 
                         xlab="Longitude", ylab="Latitude", main="Cascadia Subduction Zone", 
                         clab="", addLegend=TRUE, lwd=1, includeFault=TRUE, 
                         gpsHull=NULL, xName="lon", yName="lat", 
                         projection=NULL, parameters=NULL, orientation=NULL) {
  library(ggmap)
  library(ggplot2)
  library(scales)
  library(RColorBrewer)
  library(latex2exp)
  
  ### Get a map
  # https://mapstyle.withgoogle.com/
  # https://console.cloud.google.com/home/dashboard?consoleReturnUrl=https:%2F%2Fcloud.google.com%2Fmaps-platform%2Fterms%2F%3Fapis%3Dmaps%26project%3Dspatialstatisticalmapping&consoleUI=CLOUD&mods=metropolis_maps&project=spatialstatisticalmapping&organizationId=657476903663
  # map <- get_map(location = c(lon[1], lat[1], lon[2], lat[2]), zoom = 6,
  #                maptype = "terrain", source = "google")
  # style1=c(feature="administrative", element="labels", visibility="off")
  # style2=c("&style=", feature="road", element="geometry", visibility="off")
  # style3=c("&style=", feature="poi", element="labels", visibility="off")
  # style4=c("&style=", feature="landscape", element="labels", visibility="off")
  # style5=c("&style=", feature="administrative", element="geometry.stroke", color="black")
  # style6=c("&style=", feature="administrative", element="geometry.stroke", weight=.75)
  # api_key =  'AIzaSyAA5-1cW2j6Q0_-xpuWF0YB6alAXkCBsmM' # key to my google account
  # map <- get_googlemap(center=c(lon = mean(xlim), lat = mean(ylim)), zoom=5,
  #                      style=c(style1, style2, style3, style4, style5, style6), 
  #                      key=api_key)
  
  ##### Let's try with greyed out land:
  library(maps)
  library(mapdata)
  
  # choose color if necessary (lightblue1 or white)
  # fields.color.picker()
  
  # get relevant map data
  if(is.null(orientation) && !is.null(projection)) {
    warning("no orientation specified, so projection being set to the last used projection")
    projection = ""
    parameters = NULL
  }
  
  # get relevant map data
  if(!is.null(projection)) {
    states <- map_data("state", projection=projection, parameters=parameters, orientation=orientation)
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada", projection=projection, parameters=parameters, orientation=orientation)
  } else {
    states <- map_data("state")
    west_coast <- subset(states, region %in% c("california", "oregon", "washington"))
    canada = map_data("world", "Canada")
  }
  
  # try plotting on an interpolated grid
  library(gstat)
  library(sp)
  library(maptools)
  library(RColorBrewer)
  
  # set x and y limits if necessary
  if(is.null(xlim) && is.null(ylim)) {
    xlim=c(-128, -122)
    ylim=c(39.5, 50.5)
    proj<- mapproject(xlim, ylim) # if projection unspecified, last projection is used
    xlim = proj$x
    ylim = proj$y
    xScale = scale_x_continuous("Easting", c(-.02, 0, .02, .04), labels=c("-.02", "0", ".02", ".04"), limits=c(-360, 360))
    ylab = "Northing"
    xlab = "Easting"
  } else if(is.null(xlim)) {
    xlim = range(gpsDat[[xName]])
  } else if(is.null(ylim)) {
    ylim = range(gpsDat[[yName]])
  } else {
    xScale = scale_x_continuous("Longitude", c(-127, -125, -123), labels=c("-127", "", "-123"), limits=c(-360, 360))
  }
  
  # plot it (choose background color with fields.color.picker())
  # bg = ggplot() + 
  #   geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
  #   geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
  #   coord_fixed(xlim = lon,  ylim = lat, ratio = 1.3, expand=FALSE) + 
  #   theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #         panel.grid.minor = element_blank(), 
  #         panel.background = element_rect(fill='lightblue1')) + 
  #   labs(x="Longitude", y="Latitude")
  if(identical(projection, "")) {
    bg = ggplot() + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, ratio = 1.3, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  } else {
    bg = ggplot() + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = "grey", color = "black", data=west_coast) + 
      coord_fixed(xlim = xlim,  ylim = ylim, expand=FALSE) + 
      labs(x=xlab, y=ylab) + xScale + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_rect(fill='lightblue1'))
  }
  
  # construct bounding polygon around data:
  #make concave hull around prediction mask to get final prediction points
  locking30km = gpsDat
  if(is.null(gpsHull)) {
    source("model1/seaDefAnis.r")
    if(identical(projection, "")) {
      gpsHull = ahull(locking30km[[xName]], locking30km[[yName]], alpha=.05)
    } else {
      gpsHull = ahull(locking30km[[xName]], locking30km[[yName]], alpha=2)
    }
    
  }
  indx=gpsHull$arcs[,"end1"]  
  hullPts <- cbind(locking30km[[xName]], locking30km[[yName]])[indx,] # extract the boundary points from maskXY
  
  #plot hull and data to make sure it works
  # plotSubPoly(rbind(hullPts, hullPts[1,]), cbind(locking30km$lon, locking30km$lat))
  
  # now subset prediction data frame to only be predictions within the polygon from alphahull
  library(akima)
  
  # interpolate between observations for maximum purdyness
  predGrid = make.surface.grid(list(x=seq(xlim[1], xlim[2], l=500), lat=seq(ylim[1], ylim[2], l=500)))
  preds = interpp(locking30km[[xName]], locking30km[[yName]], plotVar, predGrid[,1], predGrid[,2])
  maskFinalPreds = in.poly(cbind(preds$x, preds$y), hullPts, convex.hull=FALSE)
  preds = data.frame(preds)
  finalPreds = preds[maskFinalPreds,]
  
  if(identical(projection, "")) {
    p1 = bg + geom_tile(data = finalPreds, aes(x = x, y = y, fill = z)) + 
      scale_fill_distiller(clab, palette = "Spectral", direction=-1) + 
      geom_point(aes(x=x, y=y), pch=20, col="black", size=.1, data=locking30km) + 
      ggtitle(main) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)
  } else {
    p1 = bg + geom_tile(data = finalPreds, aes(x = x, y = y, fill = z)) + 
      scale_fill_distiller(clab, palette = "Spectral", direction=-1) + 
      geom_point(aes(x=lon, y=lat), pch=20, col="black", size=.1, data=locking30km) + 
      ggtitle(main) + 
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=canada) +
      geom_polygon(aes(x = long, y = lat, group = group), fill = rgb(0,0,0,0), color = "black", data=west_coast)
  }
  
  p1
}

# function for rotating the entire fault counterclockwise
rotateFault = function(fault, angle, origin=c(fault$longitude[1], fault$latitude[2]), 
                       xName="longitude", yName="latitude", updateStrike=TRUE) {
  angle = angle %% 360
  angleRad = angle * (2*pi/360)
  
  # construct a rotation matrix
  A = matrix(c(cos(angleRad), sin(angleRad), -sin(angleRad), cos(angleRad)), ncol=2)
  
  # center the coordinates
  newFault = fault
  x = rbind(newFault[[xName]], newFault[[yName]])
  x = sweep(x, 1, origin, "-")
  
  # now rotate the fault and add origin back in
  newX = t(A %*% x)
  newX = sweep(newX, 2, origin, "+")
  newFault[[xName]] = newX[,1]
  newFault[[yName]] = newX[,2]
  if(updateStrike)
    newFault$strike = (newFault$strike - angle) %% 360
  
  newFault
}

# function for straightening out a planar fault geometry. Works by 
# successively rotating the fault around the bottom center point of 
# the down dip sub faults, working from bottom to top
straightenFault = function(fault, depthSplit=10000) {
  strikes = fault$strike
  strikes0 = strikes == 0
  depths = fault$depth
  
  # # fix the origin of rotation
  # # rotate around the southern midpoint of the sub fault plane
  # bottomRight = getSubfaultCenter(fault[1,], index = 3, "corners")[1:2]
  # bottomLeft = getSubfaultCenter(fault[1,], index = 4, "corners")[1:2]
  # origin = c(0.5 * (bottomRight[1] + bottomLeft[1]), 0.5 * (bottomRight[2] + bottomLeft[2]))
  
  while(!all(strikes0)) {
    # find the bottom-most sub fault with a strike that's off and the next down dip sub fault
    firstOff = match(FALSE, strikes0)
    if(is.na(firstOff))
      break
    
    if(firstOff == 1) {
      # rotate around the southern midpoint of the sub fault plane
      originI = match(TRUE, depths[firstOff:length(depths)] >= depthSplit) + firstOff - 1
      bottomRight = getSubfaultCenter(fault[originI,], index = 3, "corners")[1:2]
      bottomLeft = getSubfaultCenter(fault[originI,], index = 4, "corners")[1:2]
      origin = c(0.5 * (bottomRight[1] + bottomLeft[1]), 0.5 * (bottomRight[2] + bottomLeft[2]))
    }
    else {
      # rotate around the northern midpoint of the sub fault plane
      originI = match(TRUE, rev(depths[1:(firstOff - 1)] >= depthSplit))
      originI = firstOff - originI
      topRight = getSubfaultCenter(fault[originI,], index = 2, "corners")[1:2]
      topLeft = getSubfaultCenter(fault[originI,], index = 1, "corners")[1:2]
      origin = c(0.5 * (topRight[1] + topLeft[1]), 0.5 * (topRight[2] + topLeft[2]))
    }
    
    # rotate the rest of the fault
    rotateAngle = strikes[firstOff]
    
    if(firstOff == 1)
      fault = rotateFault(fault, rotateAngle, origin)
    else
      fault = rbind(fault[1:(firstOff - 1),], rotateFault(fault[firstOff:nrow(fault),], rotateAngle, origin))
    
    # shift the fault so the rotation leaves no gaps
    if(firstOff != 1) {
      topRightPrevious = getFaultCenters(fault[1:(firstOff - 1),], 2, "corners")[,2]
      topPrevious = max(topRightPrevious)
      bottomRightNew = getFaultCenters(fault[firstOff:nrow(fault),], 3, "corners")[,2]
      bottomNew = min(bottomRightNew)
      shift = bottomNew - topPrevious
      if(topPrevious <= bottomNew)
        fault[firstOff:nrow(fault),]$latitude = fault[firstOff:nrow(fault),]$latitude - shift
    }
    
    # update the strikes
    strikes = fault$strike
    strikes0 = strikes == 0
  }
  
  fault
}

# a new version of straightenFault that first projects the faults in euclidean space and 
# then rotates the sub faults in the same manner as above
straightenFault2 = function(fault, depthSplit=10000, projection="rectangular", parameters=NULL) {
  # first get the center and corner locations of the sub faults
  topLeft = getFaultCenters(fault, 1, "corners")[,1:2]
  topRight = getFaultCenters(fault, 2, "corners")[,1:2]
  topMiddle = getFaultCenters(fault, 1, "centers")[,1:2]
  bottomRight = getFaultCenters(fault, 3, "corners")[,1:2]
  bottomLeft = getFaultCenters(fault, 4, "corners")[,1:2]
  bottomMiddle = getFaultCenters(fault, 3, "centers")[,1:2]
  
  # set default projection parameters
  if(is.null(parameters)) {
    if(projection == "rectangular")
      parameters = mean(topMiddle[,2])
  }
  
  # now stereographically project the coordinates into euclidean space
  proj<- mapproject(topLeft[,1], topLeft[,2], projection=projection, parameters=parameters)
  topLeft = cbind(proj$x, proj$y)
  # projection = .Last.projection
  proj<- mapproject(topRight[,1], topRight[,2]) # if projection unspecified, last projection is used
  topRight = cbind(proj$x, proj$y)
  proj<- mapproject(topMiddle[,1], topMiddle[,2])
  topMiddle = cbind(proj$x, proj$y)
  proj<- mapproject(bottomRight[,1], bottomRight[,2])
  bottomRight = cbind(proj$x, proj$y)
  proj<- mapproject(bottomMiddle[,1], bottomMiddle[,2])
  bottomMiddle = cbind(proj$x, proj$y)
  proj<- mapproject(bottomLeft[,1], bottomLeft[,2])
  bottomLeft = cbind(proj$x, proj$y)
  
  # add the updated coordinates to the fault
  newFault = cbind(fault, 
                   data.frame(topLeftX=topLeft[,1], topLeftY=topLeft[,2], 
                              topRightX=topRight[,1], topRightY=topRight[,2], 
                              topMiddleX=0.5 * (topLeft[,1] + topRight[,1]), topMiddleY=0.5 * (topLeft[,2] + topRight[,2]), 
                              bottomRightX=bottomRight[,1], bottomRightY=bottomRight[,2], 
                              bottomLeftX=bottomLeft[,1], bottomLeftY=bottomLeft[,2], 
                              bottomMiddleX=0.5 * (bottomLeft[,1] + bottomRight[,1]), bottomMiddleY=0.5 * (bottomLeft[,2] + bottomRight[,2])))
  
  # now rotate and shift as in the above straightenFault function
  strikes = newFault$strike
  strikes0 = strikes == 0
  depths = newFault$depth
  
  # # fix the origin of rotation
  # # rotate around the southern midpoint of the sub fault plane
  # bottomRight = getSubfaultCenter(newFault[1,], index = 3, "corners")[1:2]
  # bottomLeft = getSubfaultCenter(newFault[1,], index = 4, "corners")[1:2]
  # origin = c(0.5 * (bottomRight[1] + bottomLeft[1]), 0.5 * (bottomRight[2] + bottomLeft[2]))
  
  while(!all(strikes0)) {
    # find the bottom-most sub fault with a strike that's off and the next down dip sub fault
    firstOff = match(FALSE, strikes0)
    if(is.na(firstOff))
      break
    
    if(firstOff == 1) {
      # rotate around the southern midpoint of the sub fault plane
      originI = match(TRUE, depths[firstOff:length(depths)] >= depthSplit) + firstOff - 1
      bottomRight = c(newFault$bottomRightX[originI], newFault$bottomRightY[originI])
      bottomLeft = c(newFault$bottomMiddleX[originI], newFault$bottomMiddleY[originI])
      origin = c(newFault$bottomMiddleX[originI], newFault$bottomMiddleY[originI])
    }
    else {
      # rotate around the northern midpoint of the sub fault plane
      originI = match(TRUE, rev(depths[1:(firstOff - 1)] >= depthSplit))
      originI = firstOff - originI
      origin = c(newFault$topMiddleX[originI], newFault$topMiddleY[originI])
    }
    
    # rotate the rest of the fault
    rotateAngle = strikes[firstOff]
    
    if(firstOff == 1) {
      newFault = rotateFault(newFault, rotateAngle, origin, xName="topLeftX", yName="topLeftY")
      newFault = rotateFault(newFault, rotateAngle, origin, xName="topMiddleX", yName="topMiddleY", updateStrike=FALSE)
      newFault = rotateFault(newFault, rotateAngle, origin, xName="topRightX", yName="topRightY", updateStrike=FALSE)
      newFault = rotateFault(newFault, rotateAngle, origin, xName="bottomLeftX", yName="bottomLeftY", updateStrike=FALSE)
      newFault = rotateFault(newFault, rotateAngle, origin, xName="bottomMiddleX", yName="bottomMiddleY", updateStrike=FALSE)
      newFault = rotateFault(newFault, rotateAngle, origin, xName="bottomRightX", yName="bottomRightY", updateStrike=FALSE)
    }
    else {
      newFault = rbind(newFault[1:(firstOff - 1),], 
                       rotateFault(newFault[firstOff:nrow(newFault),], rotateAngle, origin, xName="topLeftX", yName="topLeftY"))
      newFault = rbind(newFault[1:(firstOff - 1),], 
                       rotateFault(newFault[firstOff:nrow(newFault),], rotateAngle, origin, xName="topRightX", yName="topRightY", updateStrike=FALSE))
      newFault = rbind(newFault[1:(firstOff - 1),], 
                       rotateFault(newFault[firstOff:nrow(newFault),], rotateAngle, origin, xName="topMiddleX", yName="topMiddleY", updateStrike=FALSE))
      newFault = rbind(newFault[1:(firstOff - 1),], 
                       rotateFault(newFault[firstOff:nrow(newFault),], rotateAngle, origin, xName="bottomLeftX", yName="bottomLeftY", updateStrike=FALSE))
      newFault = rbind(newFault[1:(firstOff - 1),], 
                       rotateFault(newFault[firstOff:nrow(newFault),], rotateAngle, origin, xName="bottomMiddleX", yName="bottomMiddleY", updateStrike=FALSE))
      newFault = rbind(newFault[1:(firstOff - 1),], 
                       rotateFault(newFault[firstOff:nrow(newFault),], rotateAngle, origin, xName="bottomRightX", yName="bottomRightY", updateStrike=FALSE))
    }
    
    # shift the fault so the rotation leaves no gaps
    if(firstOff != 1) {
      topPrevious = max(newFault[1:(firstOff - 1),]$topRightY)
      bottomNew = min(newFault[firstOff:nrow(newFault),]$bottomRightY)
      shift = bottomNew - topPrevious
      if(topPrevious <= bottomNew) {
        newFault[firstOff:nrow(newFault),]$topRightY = newFault[firstOff:nrow(newFault),]$topRightY - shift
        newFault[firstOff:nrow(newFault),]$topLeftY = newFault[firstOff:nrow(newFault),]$topLeftY - shift
        newFault[firstOff:nrow(newFault),]$topMiddleY = newFault[firstOff:nrow(newFault),]$topMiddleY - shift
        newFault[firstOff:nrow(newFault),]$bottomRightY = newFault[firstOff:nrow(newFault),]$bottomRightY - shift
        newFault[firstOff:nrow(newFault),]$bottomLeftY = newFault[firstOff:nrow(newFault),]$bottomLeftY - shift
        newFault[firstOff:nrow(newFault),]$bottomMiddleY = newFault[firstOff:nrow(newFault),]$bottomMiddleY - shift
      }
    }
    
    # update the strikes
    strikes = newFault$strike
    strikes0 = strikes == 0
  }
  
  # shift the two top-left subfaults downward slightly to be 
  # aligned with the two top-right subfaults
  topLeftI = 20
  topLeftIs = 19:20
  topRightI = 18
  shift = newFault$topMiddleY[topLeftI] - newFault$topMiddleY[topRightI]
  newFault$topMiddleY[topLeftIs] = newFault$topMiddleY[topLeftIs] - shift
  newFault$bottomMiddleY[topLeftIs] = newFault$bottomMiddleY[topLeftIs] - shift
  
  # lengthen the third from top on the left subfault
  thirdTopLeftI = 16
  thirdTopRightI = 15
  newFault$topMiddleY[thirdTopLeftI] = newFault$topMiddleY[thirdTopRightI]
  newFault$bottomMiddleY[thirdTopLeftI] = newFault$bottomMiddleY[thirdTopRightI]
  
  # lengthen the fourth from top on the left subfault
  fourthTopLeftI = 14
  fourthTopRightI = 13
  newFault$topMiddleY[fourthTopLeftI] = newFault$topMiddleY[fourthTopRightI]
  
  # modify the subfaults so they have right angles
  leftMiddleX = 0.5*(newFault$topLeftX + newFault$bottomLeftX)
  rightMiddleX = 0.5*(newFault$topRightX + newFault$bottomRightX)
  
  newFault$topLeftX = leftMiddleX
  newFault$bottomLeftX = leftMiddleX
  newFault$topRightX = rightMiddleX
  newFault$bottomRightX = rightMiddleX
  
  newFault$topLeftY = newFault$topMiddleY
  newFault$topRightY = newFault$topMiddleY
  newFault$bottomLeftY = newFault$bottomMiddleY
  newFault$bottomRightY = newFault$bottomMiddleY
  
  # plot results
  ggPlotFaultDat(newFault, plotData=FALSE, projection=projection, xName="topLeftX", yName="topLeftY", parameters=parameters, 
                 xlim=NULL, ylim=NULL, main="'Straightened' Fault")
  
  newFault
}

# a new version of straightenFault2 that aligns subfaults with each other perfectly
straightenFault3 = function(projection="rectangular", parameters=NULL) {
  fault = faultGeom
  
  # first get the center and corner locations of the sub faults
  topLeft = getFaultCenters(fault, 1, "corners")[,1:2]
  topRight = getFaultCenters(fault, 2, "corners")[,1:2]
  topMiddle = getFaultCenters(fault, 1, "centers")[,1:2]
  bottomRight = getFaultCenters(fault, 3, "corners")[,1:2]
  bottomLeft = getFaultCenters(fault, 4, "corners")[,1:2]
  bottomMiddle = getFaultCenters(fault, 3, "centers")[,1:2]
  
  # set default projection parameters
  if(is.null(parameters)) {
    if(projection == "rectangular")
      parameters = mean(topMiddle[,2])
  }
  
  # now stereographically project the coordinates into euclidean space
  proj<- mapproject(topLeft[,1], topLeft[,2], projection=projection, parameters=parameters)
  topLeft = cbind(proj$x, proj$y)
  # projection = .Last.projection
  proj<- mapproject(topRight[,1], topRight[,2]) # if projection unspecified, last projection is used
  topRight = cbind(proj$x, proj$y)
  proj<- mapproject(topMiddle[,1], topMiddle[,2])
  topMiddle = cbind(proj$x, proj$y)
  proj<- mapproject(bottomRight[,1], bottomRight[,2])
  bottomRight = cbind(proj$x, proj$y)
  proj<- mapproject(bottomMiddle[,1], bottomMiddle[,2])
  bottomMiddle = cbind(proj$x, proj$y)
  proj<- mapproject(bottomLeft[,1], bottomLeft[,2])
  bottomLeft = cbind(proj$x, proj$y)
  
  # calculate the overall bottom right for plotting, shifting coordinate system at the end
  bottomRightOrigin = c(max(bottomRight[,1]), min(bottomRight[,2]))
  
  # calculate 'average' height and width of the subfaults (chose these estimates to match old distances)
  middleRight = 0.5 * (topRight + bottomRight)
  middleLeft = 0.5 * (topLeft + bottomLeft)
  l = max(rdist(topMiddle, bottomMiddle)) / 10
  w = mean(rdist.vec(middleRight, middleLeft))
  
  # modify the subfaults to have the correct length and width
  faultRow = c(1, 2, 1, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 10, 9, 10) - 1
  faultCol = 0 + (fault$depth >= 10000)
  bottomLeft = cbind(faultCol * w, faultRow * l)
  topLeft = cbind(faultCol * w, (faultRow + 1) * l)
  topRight = cbind((faultCol + 1) * w, (faultRow + 1) * l)
  bottomRight = cbind((faultCol + 1) * w, faultRow * l)
  
  # add back in the origin frame of reference
  rightShift = bottomRightOrigin[1] - 2 * w
  upShift = bottomRightOrigin[2]
  bottomLeft = sweep(bottomLeft, 2, c(rightShift, upShift), "+")
  topLeft = sweep(topLeft, 2, c(rightShift, upShift), "+")
  topRight = sweep(topRight, 2, c(rightShift, upShift), "+")
  bottomRight = sweep(bottomRight, 2, c(rightShift, upShift), "+")
  
  # add the updated coordinates to the fault
  newFault = cbind(fault, 
                   data.frame(topLeftX=topLeft[,1], topLeftY=topLeft[,2], 
                              topRightX=topRight[,1], topRightY=topRight[,2], 
                              topMiddleX=0.5 * (topLeft[,1] + topRight[,1]), topMiddleY=0.5 * (topLeft[,2] + topRight[,2]), 
                              bottomRightX=bottomRight[,1], bottomRightY=bottomRight[,2], 
                              bottomLeftX=bottomLeft[,1], bottomLeftY=bottomLeft[,2], 
                              bottomMiddleX=0.5 * (bottomLeft[,1] + bottomRight[,1]), bottomMiddleY=0.5 * (bottomLeft[,2] + bottomRight[,2])))
  
  # plot results
  ggPlotFaultDat(newFault, plotData=FALSE, projection=projection, xName="topLeftX", yName="topLeftY", parameters=parameters, 
                 xlim=NULL, ylim=NULL, main="'Straightened' Fault")
  
  newFault
}

# a new version of straightenFault2 that aligns subfaults with each other perfectly
straightenFaultLambert = function(fault=faultGeom, initParameters=range(fault$latitude)) {
  
  # get original distances between the points and the fault
  middle = getFaultCenters(fault)[,1:2]
  originalDist = rdist.earth(middle, miles = FALSE)
  originalDist = originalDist[lower.tri(originalDist)]
  
  # minimize the absolute residual distances relative to spherical distance
  fun = function(par) {
    scaleToOriginal = par[1]
    proj<- mapproject(middle[,1], middle[,2], projection="lambert", parameters=par[2:3]) # if projection unspecified, last projection is used
    newMiddle = cbind(proj$x, proj$y) * scaleToOriginal
    newDist = rdist(newMiddle)
    newDist = newDist[lower.tri(newDist)]
    mean(abs(newDist - originalDist))
  }
  
  # set initial parameters
  proj<- mapproject(middle[,1], middle[,2], projection="lambert", parameters=initParameters) # if projection unspecified, last projection is used
  newMiddle = cbind(proj$x, proj$y)
  newDist = rdist(newMiddle)
  newDist = newDist[lower.tri(newDist)]
  initScale = mean(newDist / originalDist)
  initPar = c(initScale, initParameters)
  
  # calculate optimal production parameters
  controls = list(fnscale = 1, reltol=1e-6, parscale=c(rep(1, length(initPar))))
  opt = optim(initPar, fun, control = controls)
  optPar = opt$par
  scale = optPar[1]
  LambertPar = optPar[2:3]
  
  # first get the center and corner locations of the sub faults
  topLeft = getFaultCenters(fault, 1, "corners")[,1:2]
  topRight = getFaultCenters(fault, 2, "corners")[,1:2]
  topMiddle = getFaultCenters(fault, 1, "centers")[,1:2]
  bottomRight = getFaultCenters(fault, 3, "corners")[,1:2]
  bottomLeft = getFaultCenters(fault, 4, "corners")[,1:2]
  bottomMiddle = getFaultCenters(fault, 3, "centers")[,1:2]
  
  # now project the coordinates into euclidean space
  proj<- mapproject(topLeft[,1], topLeft[,2])
  topLeft = cbind(proj$x, proj$y) * scale
  # projection = .Last.projection
  proj<- mapproject(topRight[,1], topRight[,2]) # if projection unspecified, last projection is used
  topRight = cbind(proj$x, proj$y) * scale
  proj<- mapproject(topMiddle[,1], topMiddle[,2])
  topMiddle = cbind(proj$x, proj$y) * scale
  proj<- mapproject(bottomRight[,1], bottomRight[,2])
  bottomRight = cbind(proj$x, proj$y) * scale
  proj<- mapproject(bottomMiddle[,1], bottomMiddle[,2])
  bottomMiddle = cbind(proj$x, proj$y) * scale
  proj<- mapproject(bottomLeft[,1], bottomLeft[,2])
  bottomLeft = cbind(proj$x, proj$y) * scale
  
  # Now generate strike and dip axes/coordinates using PCA
  proj<- mapproject(middle[,1], middle[,2]) # if projection unspecified, last projection is used
  newMiddle = cbind(proj$x, proj$y) * scale
  center = colMeans(newMiddle)
  X = sweep(newMiddle, 2, center)
  xtx = t(X) %*% X
  out = eigen(xtx)
  eigenvectors = out$vectors
  eigenvectors[,2] = -eigenvectors[,2] # use down dip rather than up dip direction
  eigenvalues = out$values
  
  projectedFault = cbind(fault, 
                   data.frame(topLeftX=topLeft[,1], topLeftY=topLeft[,2], 
                              topRightX=topRight[,1], topRightY=topRight[,2], 
                              topMiddleX=0.5 * (topLeft[,1] + topRight[,1]), topMiddleY=0.5 * (topLeft[,2] + topRight[,2]), 
                              bottomRightX=bottomRight[,1], bottomRightY=bottomRight[,2], 
                              bottomLeftX=bottomLeft[,1], bottomLeftY=bottomLeft[,2], 
                              bottomMiddleX=0.5 * (bottomLeft[,1] + bottomRight[,1]), bottomMiddleY=0.5 * (bottomLeft[,2] + bottomRight[,2])))
  
  # plot results
  pl = ggPlotFaultDat(projectedFault, plotData=FALSE, projection="lambert", xName="topLeftX", yName="topLeftY", parameters=parameters, 
                      xlim=NULL, ylim=NULL, main="'Straightened' Fault", scale=scale, roundDigits=-1, xlab="Easting (km)", 
                      ylab="Northing (km)", coordsAlreadyScaled=TRUE)
  
  # now add strike axis
  slope = eigenvectors[2, 1] / eigenvectors[1, 1]
  intercept = center[2] - slope * center[1]
  pl + geom_abline(col="red", lwd=1, slope = slope, intercept = intercept / scale) + ggtitle("Projected Fault with Strike Axis (Red Line)")
  
  # modify the coordinates by projecting them onto the eigenvectors (strike, dip)
  # add the center shift back in after the transformation for plotting purposes
  topLeft = sweep(sweep(topLeft, 2, center) %*% eigenvectors, 2, rev(center), "+")
  topRight = sweep(sweep(topRight, 2, center) %*% eigenvectors, 2, rev(center), "+")
  topMiddle = sweep(sweep(topMiddle, 2, center) %*% eigenvectors, 2, rev(center), "+")
  bottomRight = sweep(sweep(bottomRight, 2, center) %*% eigenvectors, 2, rev(center), "+")
  bottomMiddle = sweep(sweep(bottomMiddle, 2, center) %*% eigenvectors, 2, rev(center), "+")
  bottomLeft = sweep(sweep(bottomLeft, 2, center) %*% eigenvectors, 2, rev(center), "+")
  
  # make a function for the transformation to the new space
  transformation = function(coords) {
    proj<- mapproject(coords[,1], coords[,2]) # if projection unspecified, last projection is used
    coords = cbind(proj$x, proj$y) * scale
    coords = sweep(sweep(coords, 2, center) %*% eigenvectors, 2, rev(center), "+")
  }
  
  # add the updated coordinates to the fault
  newFault = cbind(fault, 
                   data.frame(topLeftX=topLeft[,2], topLeftY=topLeft[,1], 
                              topRightX=topRight[,2], topRightY=topRight[,1], 
                              topMiddleX=0.5 * (topLeft[,2] + topRight[,2]), topMiddleY=0.5 * (topLeft[,1] + topRight[,1]), 
                              bottomRightX=bottomRight[,2], bottomRightY=bottomRight[,1], 
                              bottomLeftX=bottomLeft[,2], bottomLeftY=bottomLeft[,1], 
                              bottomMiddleX=0.5 * (bottomLeft[,2] + bottomRight[,2]), bottomMiddleY=0.5 * (bottomLeft[,1] + bottomRight[,1])))
  
  pl = ggPlotFaultDat(newFault, plotData=FALSE, projection="lambert", xName="topLeftX", yName="topLeftY", parameters=parameters, 
                      xlim=NULL, ylim=NULL, main="'Straightened' Fault", scale=scale, roundDigits=-1, xlab="Easting (km)", 
                      ylab="Northing (km)", coordsAlreadyScaled=TRUE)
  
  # now add strike axis
  slope = eigenvectors[2, 1] / eigenvectors[1, 1]
  intercept = center[2] - slope * center[1]
  pl + geom_abline(col="red", lwd=1, slope = slope, intercept = intercept / scale) + ggtitle("Projected Fault with Strike Axis (Red Line)")
  
  list(fault=newFault, scale=scale, projPar=LambertPar, MAE=opt$value, transformation=transformation) 
}

# this function projects gps coordinates onto the new, straightened coordinate system by
# computing the gps observations' relative coordinates within each sub fault
calcStraightGpsCoords = function(newFault, oldFault=faultGeom, gpsDat=slipDatCSZ, 
                                 projection="rectangular", parameters=NULL) {
  gpsFaultIndices = getFaultIndexFromGps(gpsDat, faultGeom = oldFault)
  
  # first get the center and corner locations of the sub faults
  topLeft = getFaultCenters(oldFault, 1, "corners")[,1:2]
  topRight = getFaultCenters(oldFault, 2, "corners")[,1:2]
  topMiddle = getFaultCenters(oldFault, 1, "centers")[,1:2]
  bottomRight = getFaultCenters(oldFault, 3, "corners")[,1:2]
  bottomLeft = getFaultCenters(oldFault, 4, "corners")[,1:2]
  bottomMiddle = getFaultCenters(oldFault, 3, "centers")[,1:2]
  
  # set default projection parameters
  if(is.null(parameters)) {
    if(projection == "rectangular")
      parameters = mean(topMiddle[,2])
  }
  
  # now stereographically project the coordinates into euclidean space
  proj<- mapproject(topLeft[,1], topLeft[,2], projection=projection, parameters=parameters)
  topLeft = cbind(proj$x, proj$y)
  # projection = .Last.projection
  proj<- mapproject(topRight[,1], topRight[,2]) # if projection unspecified, last projection is used
  topRight = cbind(proj$x, proj$y)
  proj<- mapproject(topMiddle[,1], topMiddle[,2])
  topMiddle = cbind(proj$x, proj$y)
  proj<- mapproject(bottomRight[,1], bottomRight[,2])
  bottomRight = cbind(proj$x, proj$y)
  proj<- mapproject(bottomMiddle[,1], bottomMiddle[,2])
  bottomMiddle = cbind(proj$x, proj$y)
  proj<- mapproject(bottomLeft[,1], bottomLeft[,2])
  bottomLeft = cbind(proj$x, proj$y)
  
  # now do the same projection to the gps data
  proj<- mapproject(gpsDat$lon, gpsDat$lat)
  gpsCoords = cbind(proj$x, proj$y)
  
  # calculate the overall bottom right for plotting, shifting coordinate system at the end
  bottomRightOrigin = c(max(bottomRight[,1]), min(bottomRight[,2]))
  
  # calculate 'average' height and width of the subfaults (chose these estimates to match old distances)
  middleRight = 0.5 * (topRight + bottomRight)
  middleLeft = 0.5 * (topLeft + bottomLeft)
  
  # our goal is now to project every gps point on to the axes of its associated subfault. The 
  # difficult part is that the subfaults might not be perfectly rectangular, so we must account 
  # for that
  getGpsRelativeLocations = function(gpsI) {
    # get this gps observation's coordinates and associated subfault
    thisGpsCoords = gpsCoords[gpsI,]
    faultI = gpsFaultIndices[gpsI]
    
    # solution from https://stackoverflow.com/questions/28675909/relative-position-of-a-point-within-a-quadrilateral
    x1 = bottomLeft[faultI, 1]
    x2 = bottomRight[faultI, 1]
    x3 = topRight[faultI, 1]
    x4 = topLeft[faultI, 1]
    y1 = bottomLeft[faultI, 2]
    y2 = bottomRight[faultI, 2]
    y3 = topRight[faultI, 2]
    y4 = topLeft[faultI, 2]
    x = thisGpsCoords[1]
    y = thisGpsCoords[2]
    
    # shift all coordinates so they are positive
    minX = min(x1, x2, x3, x4, x) - .01
    minY = min(y1, y2, y3, y4, y) - .01
    x1 = x1 + minX
    x2 = x2 + minX
    x3 = x3 + minX
    x4 = x4 + minX
    x = x + minX
    y1 = y1 + minY
    y2 = y2 + minY
    y3 = y3 + minY
    y4 = y4 + minY
    y = y + minY
    
    # sqrtTerm = sqrt(4 * ((x3 - x4) * (y1 - y2) - (x1 - x2) * (y3 - y4)) * 
    #                   (x4 * (-y + y1) + x1 * (y - y4) + x * (-y1 + y4)) + 
    #                   (x3 * y - x4 * y - x3 * y1 + 2 * x4 * y1 - x4 * y2 + x1 * (y + y3 - 2 * y4) + 
    #                    x2 * (-y + y4) + x * (-y1 + y2 - y3 + y4))^2
    #                )
    # uDenominator = 2 * ((x3 - x4) * (y1 - y2) - (x1 - x2) * (y3 - y4))
    # vDenominator = 2 * ((x1 - x4) * (y2 - y3) - (x2 - x3) * (y1 - y4))
    # uFirstNumerator = -x2 * y + x3 * y - x4 * y - x * y1 - x3 * y1 + 2 * x4 * y1 + x * y2 - 
    #                   x4 * y2 - x * y3 + x1 * (y + y3 - 2 * y4) + x * y4 + x2 * y4
    # vFirstNumerator1 = x2 * y - x3 * y + x4 * y + x * y1 - 2 * x2 * y1 + x3 * y1 - x * y2 - 
    #                    x4 * y2 + x * y3 - x1 * (y - 2 * y2 + y3) - x * y4 + x2 * y4
    # vFirstNumerator2 = x1 * y + x3 * y - x4 * y - x * y1 - x3 * y1 + x * y2 - 2 * x1 * y2 + 
    #                    x4 * y2 - x * y3 + x1 * y3 + x * y4 - x2 * (y - 2 * y1 + y4)
    # solution1 = c(-(uFirstNumerator + sqrtTerm) / uDenominator, (vFirstNumerator1 + sqrtTerm) / vDenominator)
    # solution2 = c((-uFirstNumerator + sqrtTerm) / uDenominator, -(vFirstNumerator2 + sqrtTerm) / vDenominator)
    # 
    # # check which of the solutions are between -1 and 1
    # inRange = function(x) { (x >= -1) & (x <= 1) }
    # if(all(inRange(solution1)))
    #   solution1
    # else
    #   solution2
    
    v_sqrt = sqrt(
      4 * (
        (x3 - x4) * (y1 - y2) - (x1 - x2) * (y3 - y4)) * (x4 * (-1 * y + y1) + x1 * (y - y4) + x * (-1 * y1 + y4)) +
        (x3 * y - x4 * y - x3 * y1 + 2 * x4 * y1 - x4 * y2 + x1 * (y + y3 - 2 * y4) + x2 * (-1 * y + y4) + x * (-1 * y1 + y2 - y3 + y4))^2
    )
    
    u_sqrt = sqrt(
      4 * ((x3 - x4) * (y1 - y2) - (x1 - x2) * (y3 - y4))
      * (
        x4 * (-1 * y + y1) + x1 * (y - y4) + x * (-1 * y1 + y4)
      ) +
        (x3 * y - x4 * y - x3 * y1 + 2 * x4 * y1 - x4 * y2 + x1 * (y + y3 - 2 * y4) + x2 * (-1 * y + y4) + x * (-1 * y1 + y2 - y3 + y4))^2
    )
    
    k = 1 / (2 * ((x3 - x4) * (y1 - y2) - (x1 - x2) * (y3 - y4)))
    l = 1 / (2 * ((x1 - x4) * (y2 - y3) - (x2 - x3) * (y1 - y4)))
    
    v1 = l *
      (x2 * y - x3 * y + x4 * y + x * y1 - 2 * x2 * y1 + x3 * y1 - x * y2 - x4 * y2 + x * y3 - x1 * (y - 2 * y2 + y3) - x * y4 + x2 * y4 +
         v_sqrt)
    
    u1 = -1 * k *
      (-x2 * y + x3 * y - x * y1 - x3 * y1 + 2 * x4 * y1 + x * y2 - x4 * y2 - x * y3 + x1 * (y + y3 - 2 * y4) + x * y4 + x2 * y4 +
         u_sqrt)
    
    v2 = -1 * l *
      (x1 * y + x3 * y - x4 * y - x * y1 - 2 * x3 * y1 + x * y2 - -2 * x1 * y2 + x4 * y2 - x * y3 + x1 * y3 + x * y4 - x2 * (y - 2 * y1 + y4) +
         v_sqrt)
    
    u2 = k *
      (x2 * y - x3 * y + x4 * y + x * y1 + x3 * y1 - 2 * x4 * y1 - x * y2 + x4 * y2 + x * y3 - x1 * (y + y3 - 2 * y4) - x * y4 - x2 * y4 +
         u_sqrt)
    
    solution1 = c(u1, v1)
    solution2 = c(u2, v2)
    
    a = (x1 * y3 - x1 * y4 - x2 * y3 + x2 * y4 + (-x3) * y1 + x3 * y2 + x4 * y1 - x4 * y2)

    # check which of the solutions are between 0 and 1
    inRange = function(x) { (x >= 0) & (x <= 1) }
    if(all(inRange(solution1)))
      solution1
    else if(all(inRange(solution2)))
      solution2
    else {
      if(abs(a) >= .0001)
        warning(paste0("a is relatively large (", round(a, digits=2), "), but both quadratic solutions are bad. Using linear solution"))
      # in this case, some of the edges must've been parallel, so we had a linear equation not quadratic
      # calculate u, the relative x coordinate
      b = (-x1 * y3 + 2 * x1 * y4 - x1 * y - x2 * y4 + x2 * y + x3 * y1 - x3 * y - 2 * x4 * y1 + x4 * y2 + x4 * y + x * y1 - x * y2 + x * y3 - x * y4)
      c = (-x1 * y4 + x1 * y + x4 * y1 - x4 * y - x * y1 + x * y4)
      u = -c / b
      
      # now calculate v, the relative y coordinate
      relativeBottom = c(x1 + u * (x2 - x1), y1 + u * (y2 - y1))
      relativeTop = c(x4 + u * (x3 - x4), y4 + u * (y3 - y4))
      vDists = rdist.vec(rbind(c(x, y), relativeTop), rbind(relativeBottom, relativeBottom))
      v = vDists[1] / vDists[2]
      solution = c(u, v)
      solution
    }
  }
  relativeCoords = t(sapply(1:nrow(gpsCoords), getGpsRelativeLocations))
  
  # now that we have the relative coordinates of the gps locations, we can compute their up strike and
  # down dip distances
  getGpsStrikeDipCoords = function(gpsI) {
    # get this gps observation's relative coordinates within associated subfault
    thisGpsRelativeCoords = relativeCoords[gpsI,]
    faultI = gpsFaultIndices[gpsI]
    
    # convert from relative coordinates to strike and dip coordinates
    alphaDip = thisGpsRelativeCoords[1]
    alphaStrike = thisGpsRelativeCoords[2]
    thisHeight = newFault[faultI, ]$topRightY - newFault[faultI, ]$bottomRightY
    thisWidth = newFault[faultI, ]$topRightX - newFault[faultI, ]$topLeftX
    gpsX = newFault$bottomLeftX[faultI] + alphaDip * thisWidth
    gpsY = newFault$bottomLeftY[faultI] + alphaStrike * thisHeight
    c(gpsX, gpsY)
  }
  t(sapply(1:nrow(gpsCoords), getGpsStrikeDipCoords))
}











