
##### precomputations for the model fitting

library(fields)
setwd("~/git/M9/")
source('plotSubfault.R')

##### read in data
source('loadFloodDat.R')

##### split the fault geometry up and compute the Okada model for each subfault and unit slips
##### NOTE: we only need to save the seafloor deformations corresponding to each subsidence 
##### data location

# split fault into smaller subfaults
csz = divideFault(testFault)

# compute unit slip Okada seadef
nx = 300
ny=  900
lonGrid = seq(lonRange[1], lonRange[2], l=nx)
latGrid = seq(latRange[1], latRange[2], l=ny)
dZ = okadaAll(csz, lonGrid, latGrid, cbind(Lon, Lat), slip=1)

##### estimate event ages and interarrival times assuming that 
##### each event is a full length rupture (for now.  That 
##### assumption seems very problematic)
getAgeEsts = function(dr1, useBP=TRUE, useRCYBP=TRUE) {
  # NOTE: ns as returned by this function counts bp and rcybp as two separate data points
  # which double counts observations with both
  ageMeans = rep(0, length(uniqueEvents))
  ageMeds = rep(0, length(uniqueEvents))
  ns = rep(0, length(uniqueEvents))
  for(e in 1:length(uniqueEvents)) {
    eventDat = dr1[dr1$event == uniqueEvents[e],]
    # NOTE: bp is calendar years before 1950 and rcybp is radiocarbon years
    # before 1950.  Add 66 to estimated age since it's 2016
    
    # set which dataing methods should be used
    if(useBP) {
      estsAndWeights = cbind(eventDat$bpCntr, 1/eventDat$bpSE^2)
      if(useRCYBP) {
        estsAndWeights = rbind(estsAndWeights, cbind(eventDat$rcybpCntr, 1/eventDat$rcybpSE^2))
      }
    }
    else{
      if(!useRCYBP) {
        stop("at least one of useBP and useRCYBP must be set to TRUE")
      }
      estsAndWeights = cbind(eventDat$rcybpCntr, 1/eventDat$rcybpSE^2)
    }
    estsAndWeights = matrix(estsAndWeights, ncol=2)
    
    finite = is.finite(estsAndWeights[,1]) & is.finite(estsAndWeights[,2])
    estsAndWeights = matrix(estsAndWeights[finite,], ncol=2)
    ageMeans[e] = weighted.mean(estsAndWeights[,1], estsAndWeights[,2], na.rm=TRUE)
    ageMeds[e] = median(estsAndWeights[,1], na.rm=TRUE)
    ns[e] = sum(finite)
  }
  # set 1700 event age because it is known historically
  ageMeans[e] = ageMeds[e] = 250
  list(ageMeans=ageMeans, ageMeds=ageMeds, ns=ns)
}
tmp = getAgeEsts(dr1)
ageMeans = tmp$ageMeans
ageMeds = tmp$ageMeds
ns = tmp$ns
cbind(uniqueEvents, ageMeans, ageMeds, ns)
# median seems like a better estimator when comparing with
# http://activetectonics.coas.oregonstate.edu/paper_files/Kulkarni%20et%20al.,%202013.pdf

# test with only bp
tmp = getAgeEsts(dr1, useRCYBP=FALSE)
ageMeans = tmp$ageMeans
ageMeds = tmp$ageMeds
ns = tmp$ns
cbind(uniqueEvents, ageMeans, ageMeds, ns)

hist(c(dr1$bpCntr, dr1$rcybpCntr), breaks=200)
for(e in 1:length(uniqueEvents)) {
  abline(v=ageEsts[e], col="blue")
}




