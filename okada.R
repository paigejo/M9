okada = function(rows, x, y, slips=rep(0, nrow(rows)), rakes=rep(90, nrow(rows)), 
                 poisson=0.25, inds=NULL) {
  # allow user to send a constant for slips and rakes instead of vector for ease of use
  if(length(slips) == 1)
    slips = rep(slips, nrow(rows))
  if(length(rakes) == 1)
    rakes = rep(rakes, nrow(rows))
  
  fullOkada = matrix(0, nrow=length(y), ncol=length(x))
  for(i in 1:nrow(rows)) {
    fullOkada = fullOkada + okadaSubfault(rows[i,], x, y, slips[i], rakes[i], poisson, inds)
  }
  return(fullOkada)
}

# same as Okada, but return a matrix of dimension nrow(datCoords) x nrow(rows)
# with a decomposition of the seaDef at each data location induced by the slip 
# at each subfault (G matrix from model)
okadaAll = function(rows, lonGrid, latGrid, datCoords, slips=rep(0, nrow(rows)), rakes=rep(90, nrow(rows)), 
                    poisson=0.25) {
  # allow user to send a constant for slips and rakes instead of vector for ease of use
  if(length(slips) == 1)
    slips = rep(slips, nrow(rows))
  if(length(rakes) == 1)
    rakes = rep(rakes, nrow(rows))
  
  # calculate the simulated subsidence at the data locations
  # round the data locations to the lon lat grid to make this function faster
  # and to have consistent grid across different subfaults
  roundToRange = function(coordVec, coordRange) {
    inds = (coordVec - min(coordRange))/(max(coordRange) - min(coordRange))
    inds = round(inds*(length(coordRange)-1)) + 1
    roundedCoordVec = coordRange[inds]
    return(roundedCoordVec)
  }
  roundLon = roundToRange(datCoords[,1], lonGrid)
  roundLat = roundToRange(datCoords[,2], latGrid)
  roundCoords = cbind(roundLon, roundLat)
  
  # find indices of grid coords corresponding to rounded data coords
  findIndex = function(rCoords, gCoords) {
    return(match(data.frame(t(rCoords)), data.frame(t(gCoords))))
  }
  coordGrid = make.surface.grid(list(lonGrid, latGrid))
  inds = findIndex(roundCoords, coordGrid)
  
  # compute the result of the Okada model for each subfault at each data location
  fullOkada = matrix(0, nrow=nrow(datCoords), ncol=nrow(rows))
  for(i in 1:nrow(rows)) {
    # simulate the seaDef
    # dZ = okadaSubfault(rows[i,], x, y, slips[i], rakes[i], poisson)
    fullOkada[,i] = okadaSubfault(rows[i,], lonGrid, latGrid, slips[i], rakes[i], poisson, inds)
      
    # get the simulated seaDef at the data locations
    # fullOkada[i,] = c(t(dZ))[inds]
  }
  return(fullOkada)
}

# function based on dtopotools.Subfault.okada function in python
#
# x and y are 1D vectors forming the coordinates of a grid (as in 
# make.surface.grid(list(x, y))   )
# ind is a set of indices such that c(t(dZ))[inds] gives the desired seaDefs
# where dZ is the result of okadaSubfault
okadaSubfault = function(row, x, y, slip=0, rake=90, poisson=0.25, inds=NULL) {
  # Apply Okada to this subfault and return a DTopography object.
  # 
  # :Input:
  # - x,y are 1d arrays
  # :Output:
  # - DTopography object with dZ array of shape (1,len(x),len(y))
  # with single static displacement and times = [0.].
  # 
  # Currently only calculates the vertical displacement.
  # 
  # Okada model is a mapping from several fault parameters
  # to a surface deformation.
  # See Okada 1985 [Okada85]_, or Okada 1992, Bull. Seism. Soc. Am.
  # 
  # okadamap function riginally written in Python by Dave George for
  # Clawpack 4.6 okada.py routine, with some routines adapted
  # from fortran routines written by Xiaoming Wang.
  # 
  # Rewritten and made more flexible by Randy LeVeque
  # 
  # **Note:** *self.coordinate_specification* (str) specifies the location on 
  # each subfault that corresponds to the (longitude,latitude) and depth 
  # subfault.
  # 
  # See the documentation for *SubFault.calculate_geometry* for dicussion of the 
  # possible values *self.coordinate_specification* can take.
  
  if(!is.list(row))
    row = as.list(row)
  
  # Setup coordinate arrays
  #   corners = [[None, None, None], # a 
  #              [None, None, None], # b
  #              [None, None, None], # c
  #              [None, None, None]] # d
  #   centers = [[None, None, None], # 1
  #              [None, None, None], # 2 
  #              [None, None, None]] # 3
  geom = calcGeom(row)
  centers = geom$centers
  corners = geom$corners
  
  # Okada model assumes x,y are at bottom center:
  x_bottom = centers[3,1]
  y_bottom = centers[3,2]
  depth_bottom = centers[3,3]
  
  length = row$length
  width = row$width
  depth = row$depth
  
  halfL = 0.5*length
  w  =  width
  
  # convert angles to radians:
  DEG2RAD = 2*pi/360
  ang_dip = DEG2RAD * row$dip
  ang_rake = DEG2RAD * rake
  ang_strike = DEG2RAD * row$strike
  
  # this format should be the same as numpy.meshgrid
  mesh = meshgrid(x, y)   # use convention of upper case for 2d
  X = mesh$X
  Y = mesh$Y
  
  # Convert distance from (X,Y) to (x_bottom,y_bottom) from degrees to
  # meters:
  LAT2METER = 111133.84012073894 #/10^3
  xx = LAT2METER * cos(DEG2RAD * Y) * (X - x_bottom)   
  yy = LAT2METER * (Y - y_bottom)
  
  # if user only wants a subset of seaDefs given by inds, subset now to save time
  if(!is.null(inds)) {
    xx = c(t(xx))[inds]
    yy = c(t(yy))[inds]
  }
  
  # Convert to distance along strike (x1) and dip (x2):
  x1 = xx * sin(ang_strike) + yy * cos(ang_strike) 
  x2 = xx * cos(ang_strike) - yy * sin(ang_strike) 
  
  # In Okada's paper, x2 is distance up the fault plane, not down dip:
  x2 = -x2
  
  p = x2 * cos(ang_dip) + depth_bottom * sin(ang_dip)
  q = x2 * sin(ang_dip) - depth_bottom * cos(ang_dip)
  
  # to save computation time, set strike slip to 0 if rake = 90
  if(rake != 90) {
    f1 = strike_slip(x1 + halfL, p,     ang_dip, q, poisson)
    f2 = strike_slip(x1 + halfL, p - w, ang_dip, q, poisson)
    f3 = strike_slip(x1 - halfL, p,     ang_dip, q, poisson)
    f4 = strike_slip(x1 - halfL, p - w, ang_dip, q, poisson)
  }
  else {
    f1 = f2 = f3 = f4 = rep(0, length(x1))
  }
  
  g1=dip_slip(x1 + halfL, p,     ang_dip, q, poisson)
  g2=dip_slip(x1 + halfL, p - w, ang_dip, q, poisson)
  g3=dip_slip(x1 - halfL, p,     ang_dip, q, poisson)
  g4=dip_slip(x1 - halfL, p - w, ang_dip, q, poisson)
  
  # Displacement in direction of strike and dip:
  ds = slip * cos(ang_rake)
  dd = slip * sin(ang_rake)
  
  us = (f1 - f2 - f3 + f4) * ds
  ud = (g1 - g2 - g3 + g4) * dd
  
  dz = (us+ud)
  
  # I opted to only return the dZ instead of the list of objects for simplicity
  #dtopo = list()
  #dtopo$X = X
  #dtopo$Y = Y
  #dtopo$dZ = dz
  #dtopo$times = 0.
  return(dz)
}

# function to create grids of 2d arrays forming mesh based on 1d vectors
meshgrid = function(x, y, transpose=FALSE) {
  if(transpose) {
    X = matrix(rep(x, length(y)), nrow=length(x))
    Y = matrix(rep(y, length(x)), ncol=length(y), byrow = TRUE)
  }
  else { # numpy format
    X = matrix(rep(x, length(y)), nrow=length(y), byrow=TRUE)
    Y = matrix(rep(y, length(x)), ncol=length(x))
  }
  return(list(X=X, Y=Y))
}

strike_slip = function(y1, y2, ang_dip, q, poisson=0.25) {
  # Used for Okada's model
  # Methods from Yoshimitsu Okada (1985)
  
  sn = sin(ang_dip)
  cs = cos(ang_dip)
  d_bar = y2*sn - q*cs
  r = sqrt(y1**2 + y2**2 + q**2)
  xx = sqrt(y1**2 + q**2)
  a4 = 2.0*poisson/cs*(log(r+d_bar) - sn*log(r+y2))
  # rewritten to use multiplication instead of division when possible (faster)
  # f = -(d_bar*q/r/(r+y2) + q*sn/(r+y2) + a4*sn)/(2.0*pi)
  f = -(d_bar*q/(r*(r+y2)) + q*sn/(r+y2) + a4*sn)*(1/(2.0*pi))
  
  return(f)
}

dip_slip = function(y1, y2, ang_dip, q, poisson=0.25) {
  # Based on dtopotools.SubFault._strike_slip and Okada's paper (1985)
  # Added by Xiaoming Wang
  
  sn = sin(ang_dip)
  cs = cos(ang_dip)
  
  d_bar = y2*sn - q*cs;
  r = sqrt(y1**2 + y2**2 + q**2)
  xx = sqrt(y1**2 + q**2)
#   a5 = 4.*poisson/cs*atan((y2*(xx+q*cs)+xx*(r+xx)*sn)/y1/(r+xx)/cs)
#   f = -(d_bar*q/r/(r+y1) + sn*atan(y1*y2/q/r) - a5*sn*cs)/(2.0*pi)
  # rewritten to use multiplication instead of division when possible (faster)
  a5 = 4.*poisson/cs*atan((y2*(xx+q*cs)+xx*(r+xx)*sn)/(y1*(r+xx)*cs))
  f = -(d_bar*q/(r*(r+y1)) + sn*atan(y1*y2/(q*r)) - a5*sn*cs)*(1/(2.0*pi))
  
  return(f)
}








