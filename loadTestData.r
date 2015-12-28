# test emulator on small testing dataset from M9 project SageMathCloud page
wd = getwd()
library(RcppCNPy)
library(fields)
library(abind)

#####get SageMathCloud testing data:
setwd("~/git/M9/test_data")

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

setwd(wd)

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

#Functions for loading topography data
loadTopo = function(fname) {
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
      out = loadTopo(files[f])
      out$dat = out$dat[,,out$mt]
    }
    else {
      #concatenate data in out list:
      tmp = loadTopo(files[f])
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
      out = loadTopo(files[f])
      out$dat = out$dat[,,out$mt]
    }
    else {
      #concatenate data in out list:
      tmp = loadTopo(files[f])
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
    }
  }
  
  #set wd back to what it was before
  setwd(wd)
  
  return(out)
}







