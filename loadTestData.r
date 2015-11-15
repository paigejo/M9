# test emulator on small testing dataset from M9 project SageMathCloud page
setwd("~/git/M9")
library(RcppCNPy)
library(fields)

#####get SageMathCloud testing data:
setwd("~/git/M9/data")

#allHMax
allHMax = array(NA, dim=c(6, 250, 250))
for(i in 0:5) {
  fname = paste0("allHMax", i, ".npy")
  tmp = npyLoad(fname)
  allHMax[i+1, , ] = tmp
}

#topography/bathymetry, lon, lat
topo = npyLoad("CCTopo.npy")
lon = npyLoad("CCLon.npy")
lat = npyLoad("CCLat.npy")
