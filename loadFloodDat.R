# loads dr1 and GPS data among other variables into global environment.  Only loads major event data 
# corresponding to full CSZ ruptures
loadMajorFloodDat = function() {
  wd = getwd()
  setwd("~/git/M9/")
  slipDat <<- read.table("Cascadia-lockingrate-hdr.txt", header=TRUE)
  names(slipDat) <<- c("lon", "lat", "slip", "Depth", "slipErr")
  slipDat$Depth <<- slipDat$Depth*10^3 # convert from km to m
  
  attach(slipDat)
  lonRange <<- range(lon)
  latRange <<- range(lat)
  
  # fault geometry data (convert from km to m and to longitude format consistent with other data)
  faultGeom <<- read.csv("CSZe01.csv")
  faultGeom$longitude <<- faultGeom$longitude - 360
  kmCols <<- c(4, 6, 7)
  faultGeom[,kmCols] <<- faultGeom[,kmCols] * 10^3
  
  csz <<- divideFault(faultGeom)
  
  # subsidence data
  load("DR1.RData", envir=globalenv())
  attach(dr1)
  
  # get unique CSZ earthquake event names
  events <<- as.character(event)
  uniqueEvents <<- unique(events)
  sortI = rev(c(1, 12, 11, 14, 2, 15, 3, 4, 18, 5, 13, 6, 19, 7, 20, 10, 21, 8, 9, 16, 17)) # T1 is most recent, so reverse (T1 is last)
  
  uniqueEvents <<- uniqueEvents[sortI]
  
  # subset uniqueEvents so it only contains major events
  eventSubset = c(1, 2, 4, 6, 8, 10, 12, 15, 17, 19, 20, 21)
  majorEvents = uniqueEvents[eventSubset]
  
  # remove all other events from dr1
  nonMajorEvents = setdiff(uniqueEvents, majorEvents)
  for(i in 1:length(nonMajorEvents)) {
    nonEventInds = events != nonMajorEvents[i]
    dr1 <<- dr1[nonEventInds,]
  }
  uniqueEvents <<- majorEvents
  
  # sort dr1 by event
  dr1$event <<- factor(dr1$event, levels=uniqueEvents)
  inds = order(dr1$event)
  dr1 <<- dr1[inds,]
  events <<- as.character(dr1$event)
  
  #set working directory back to what it was before calling this function
  setwd(wd)
  
  invisible(list(dr1=dr1, slipDat=slipDat, csz=csz))
}

# loads dr1 and GPS data among other variables into global environment
loadFloodDat = function() {
  wd = getwd()
  setwd("~/git/M9/")
  slipDat <<- read.table("Cascadia-lockingrate-hdr.txt", header=TRUE)
  names(slipDat) <<- c("lon", "lat", "slip", "Depth", "slipErr")
  slipDat$Depth <<- slipDat$Depth*10^3 # convert from km to m
  
  attach(slipDat)
  lonRange <<- range(lon)
  latRange <<- range(lat)
  
  # fault geometry data (convert from km to m and to longitude format consistent with other data)
  faultGeom <<- read.csv("CSZe01.csv")
  faultGeom$longitude <<- faultGeom$longitude - 360
  kmCols <<- c(4, 6, 7)
  faultGeom[,kmCols] <<- faultGeom[,kmCols] * 10^3
  
  csz <<- divideFault(faultGeom)
  
  # subsidence data
  load("DR1.RData", envir=globalenv())
  attach(dr1)
  event <<- dr1$event
  
  # get unique CSZ earthquake event names
  events <<- as.character(event)
  uniqueEvents <<- unique(events)
  sortI = rev(c(1, 12, 11, 14, 2, 15, 3, 4, 18, 5, 13, 6, 19, 7, 20, 10, 21, 8, 9, 16, 17)) # T1 is most recent, so reverse (T1 is last)
  uniqueEvents <<- uniqueEvents[sortI]
  
  # sort dr1 by event
  dr1$event <<- factor(dr1$event, levels=uniqueEvents)
  inds = order(dr1$event)
  dr1 <<- dr1[inds,]
  events <<- as.character(dr1$event)
  
  # divide uncertainty by 2 because Uncertainty seems to represent 1.96 sigma
  dr1$Uncertainty <<- dr1$Uncertainty/qnorm(0.975)
  
  # subset GPS data by whether it's in the CSZ fault geometry boundaries
  slipDatCSZ <<- getFaultGPSDat()
  
  # piecewise linear spline
  coastLonLat = matrix(c(-127, 50, 
                         -124.8, 49.2, 
                         -125.85, 49.1, 
                         -124.3, 48.2, 
                         -123.8, 46.3, 
                         -124.1, 43.7, 
                         -124.5, 42.85, 
                         -124.1, 40.85, 
                         -124.4, 40.5, 
                         -124, 40), byrow = TRUE, ncol=2)
  coastLon <<- coastLonLat[,1]
  coastLat <<- coastLonLat[,2]
  
  getSubPoints = function(nPoints) {
    # get some points from the line in the piecewise spline.
    
    # first get the latitudes of the points
    maxLat = max(coastLat)
    minLat = min(coastLat)
    allCoastLats = seq(minLat, maxLat, length=nPoints)
    
    # now find alpha (for the interpolation) for each section.
    # take a weighted average in between each point
    
    # for the given latitude, compute the longitude for the 
    # piecewise spline
    getLon = function(lat) {
      # first deal with base case
      if(lat == 40)
        return(coastLon[length(coastLon)])
      
      # get index of point that is the first smaller than the given lat
      last = match(TRUE, coastLat < lat)
      # now get the point right before the first point smaller
      first = last-1
      
      # get coordinates of the points to interpolate
      lastCoords = coastLonLat[last,]
      firstCoords = coastLonLat[first,]
      
      # interpolate (alpha and 1-alpha are weights for first and last)
      alpha = (lat - lastCoords[2])/(firstCoords[2] - lastCoords[2])
      # now make sure to do the interpolation
      return(alpha*firstCoords[1] + (1-alpha)*lastCoords[1])
    }
    # Now vectorize the above function
    allLon = sapply(allCoastLats, getLon)
    allCoords = cbind(allLon, allCoastLats)
    return(allCoords)
  }
  allCoastCoords <<- getSubPoints(250)
  allCoastLon <<- allCoastCoords[,1]
  allCoastLat <<- allCoastCoords[,2]
  
  # read in Goldfinger 2012 age data
  # ages <<- read.csv("goldfingerAgeEsts.csv", header=TRUE)
  
  # remove all data past T13 since it won't be used to calculate interevent times for 
  # land-based coseismic subsidence data (which only goes up to T12)
  # ages <<- ages[1:33,]
  
  # change names of T10b and T10f for consistency with dr1 dataset
#   ages[,1] <<- as.character(ages[,1])
#   ages[25, 1] <<- "T10R1"
#   ages[29, 1] <<- "T10R2"
#   ages[,1] <<- factor(ages[,1])
  
#   # remove events not in dr1 dataset:
#   # Events in ages that aren't in dr1: T2a, T5c, T6b, T8b, T9b, T10a, T10c, T10d, T10e, T10f
#   eventsToDelete = c(3, 11, 14, 19, 24, 26, 27, 28, 29)
#   ages <<- ages[-eventsToDelete,]
  
  # make sure ages variables are in correct format
  # ages$age <<- as.numeric(as.character(ages$age))
  
  # calculate interevent times
#   age1 = ages[1:(nrow(ages)-1),2]
#   age2 = ages[2:nrow(ages),2]
#   intereventTime = c(age2 - age1, NA)
#   ages$intereventTime <<- intereventTime
#   majorEvents = c(1, 2, 4, 6, 8, 12, 15, 17, 20, 23, 30, 31, 33)
#   majorAges <<- ages[majorEvents,]
#   majorAges1 = majorAges[1:(nrow(majorAges)-1),2]
#   majorAges2 = majorAges[2:nrow(majorAges),2]
#   majorIntereventTime = c(majorAges2 - majorAges1, NA)
#   majorAges$intereventTime <<- majorIntereventTime
  
  # goldMwAll <<- read.csv("goldfinger2012Table8.csv", header=TRUE)
  # # convert from dyne cm to Nm 
  # goldMwAll$seismicMoment <<- goldMwAll$seismicMoment * 10^(-7)
  # dr1Events = c(1:10, 12:17, 19:20, 22, 24, 28:29)
  # dr1Events = 1:31 # only keep events after and including T13
  # goldMw <<- goldMwAll[dr1Events,]
  
  #set working directory back to what it was before calling this function
  setwd(wd)
  
  invisible(list(dr1=dr1, csz=csz, slipDat=slipDat, slipDatCSZ=slipDatCSZ))
}

##### function for subsetting GPS data by whether it is in CSZ fault geometry
# dat is a dataframe with lon and lat columns
getFaultGPSDat = function(dat=slipDat) {
  
  # helper function for determining if GPS data is within a specific subfault geometry
  getSubfaultGPSDat = function(i) {
    row = faultGeom[i,]
    geom = calcGeom(row)
    corners = geom$corners[,1:2]
    in.poly(cbind(dat$lon, dat$lat), corners)
  }
  
  # construct logical matrix, where each column is the result of getSubfaultGPSDat(j)
  # Hence, row represents data index, column represents subfault index.  If a data 
  # observation is in any subfault, it is in the fault.
  inSubfaults = sapply(1:nrow(faultGeom), getSubfaultGPSDat)
  inFault = apply(inSubfaults, 1, any)
  
  # return the GPS data within the fault
  dat[inFault,]
}

loadFloodDat()

#where is T9a?
# [1] "T11"   "T4a"   "T10"   "T8a"   "T10R1" "T7a"   "T8"    "T5b"   "T7"    "T2"    "T6"    "T12"   "T5a"   "T5"    "T3a"   "T4"   
# [17] "T6a"   "T9"    "T3"    "T1"   
# 
# [1] "T1"  "T4"  "T5"  "T5a" "T6"  "T7"  "T8"  "T10" "T9"  "T3"  "T2"  "T6a" "T3a" "T4a" "T11" "T12" "T5b" "T7a" "T8a" "T9a"