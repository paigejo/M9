##### Main functions for users:

# exploratoryAnalysis()
#   makes some test plots

# getResDat(defDat)
#   return(list(resids=resids, mod=mod2, Mw=Mw2, type=type2, probs=probs2, mask=mask))

# regress off:
#   on log scale: magnitude, intercept, grid cell
# use lasso to induce sparsity ? or use mask to specify domain of slip?
# fit spatial model to residuals

#setwd("~/Google Drive/UW/guttorp/code/")
setwd("~/git/M9/")
source("loadTestData.r")

exploratoryAnalysis = function() {
  #generate data matrix from deformation data
  defDat = loadAllDeformations()
  
  #transform to one column of data for each realization
  dat = matrix(nrow=dim(defDat$dat)[3], ncol=dim(defDat$dat)[1]*dim(defDat$dat)[2])
  for(i in 1:dim(defDat$dat)[3]) {
    dat[i,] = c(defDat$dat[,,i])
  }
  
  #grid of lat/lon values
  grid = make.surface.grid(list(x=seq(-127, -123.5, length=dim(defDat$dat)[1]), 
                                y=seq(49.5, 39.5, length=dim(defDat$dat)[2])))
  
  #get mask of grid cells deformed by CSZ quake
  datmax = apply(abs(dat), 2, max)
  m0 = datmax > .5
  m1 = (grid[,1] >= -126) & (grid[,1] <= -125) & (grid[,2] >= 48)
  m2 = (grid[,1] >= -125) & (grid[,1] <= -124) & (grid[,2] >= 46) & (grid[,2] <= 48.3)
  m3 = (grid[,1] >= -124.8) & (grid[,1] <= -124.2) & (grid[,2] >= 40.3) & (grid[,2] <= 46)
  mask = m0 | m1 | m2 | m3
  
  #plot mask:
  quilt.plot(grid, mask, main="Final CSZ mask")
  
  #####do regression
  maskdat = dat[,mask]
  #magnitudes (not sure about the 8.6's.  They correspond to the SS1 and SS3 cases)
  Mw = c(9.0, 9.1, 9.0, 8.9, 9.0, 8.9, 8.7, 8.8, 8.7, 8.6, 8.6, 9.1, 9.2, 9.1, 9.1, 9.2, 9.1, 8.6, 8.6)
  Mw10 = 10^Mw
  #type of earthquake (splay fault, deep)
  type = c("Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Deep", "Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Deep")
  type = as.factor(type)
  #probability weights of each observation (assume SS is half of equivalent SM probs)
  probs = c(.128, .016, .016, .318, .106, .106, .104, .078, .078, .052, .039, 
            .02, .0025, .0025, .02, .0025, .0025, .052, .039)
  # gridCell = as.factor(1:sum(mask))
  
  #repeat variables for each grid cell and each replication as necessary
  datPerRep = ncol(maskdat)
  nRep = nrow(maskdat)
  
  Mw10Fin = rep(Mw10, datPerRep)
  typeFin = rep(type, datPerRep)
  # gridCellFin = rep(gridCell, rep(nRep, datPerRep))
  
  #regress
  #mod1 = lm(c(maskdat) ~ Mw10Fin + typeFin + gridCellFin)
  mod1 = lm(maskdat ~ Mw10*type)
  
  #predict
  testIn = data.frame(Mw10=Mw10, type=type)
  testOut = predict(mod1, testIn)
  slips = array(data=0, dim=dim(dat))
  slips[,mask] = testOut
  
  #plot true field versus predictions
  for(i in 1:nrow(slips)) {
    rangeV = c(range(slips[i,]), range(dat[i,]))
    minV = min(rangeV)
    maxV = max(rangeV)
    rangeR = range(dat[i,] - slips[i,])
    minR = min(rangeR)
    maxR = max(rangeR)
    
    quilt.plot(grid, dat[i,], main=paste0(i, "th Realization, True"), zlim=c(minV, maxV))
    quilt.plot(grid, slips[i,], main=paste0(i, "th Realization, Predicted"), zlim=c(minV, maxV))
    quilt.plot(grid, dat[i,]-slips[i,], main=paste0(i, "th Realization, Residuals"), zlim=c(minR, maxR))
  }
  
  #####do regression without SS earthquakes
  #remove SS earthquakes
  SS=c(10, 11, 18, 19)
  dat2 = dat[-SS,]
  maskdat2 = maskdat[-SS,]
  Mw2 = Mw[-SS]
  Mw102 = Mw10[-SS]
  type2 = type[-SS]
  probs2 = probs[-SS]
  
  #repeat variables for each grid cell and each replication as necessary
  datPerRep2 = ncol(maskdat2)
  nRep2 = nrow(maskdat2)
  
  Mw10Fin2 = rep(Mw102, datPerRep2)
  typeFin2 = rep(type2, datPerRep2)
  # gridCellFin = rep(gridCell, rep(nRep, datPerRep))
  
  #regress
  #mod1 = lm(c(maskdat) ~ Mw10Fin + typeFin + gridCellFin)
  #NOTE: excluding interaction term leaves strange artifact from splay fault discontinuities
  mod2 = lm(maskdat2 ~ Mw102*type2)
  
  #predict
  testIn = data.frame(Mw102=Mw102, type2=type2)
  testOut = predict(mod2, testIn)
  slips = array(data=0, dim=dim(dat2))
  slips[,mask] = testOut
  
  #plot true field versus predictions
  for(i in 1:nrow(slips)) {
    rangeV = c(range(slips[i,]), range(dat2[i,]))
    minV = min(rangeV)
    maxV = max(rangeV)
    rangeR = range(dat2[i,] - slips[i,])
    minR = min(rangeR)
    maxR = max(rangeR)
    
    quilt.plot(grid, dat2[i,], main=paste0(i, "th Realization, True"), zlim=c(minV, maxV))
    quilt.plot(grid, slips[i,], main=paste0(i, "th Realization, Predicted"), zlim=c(minV, maxV))
    quilt.plot(grid, dat2[i,]-slips[i,], main=paste0(i, "th Realization, Residuals"), zlim=c(minR, maxR))
  }
  quilt.plot(grid, dat2[10,], main="Example Seafloor Deformation")
}

getResDat = function(defDat) {
  
  #transform to one column of data for each realization
  dat = matrix(nrow=dim(defDat$dat)[3], ncol=dim(defDat$dat)[1]*dim(defDat$dat)[2])
  for(i in 1:dim(defDat$dat)[3]) {
    dat[i,] = c(defDat$dat[,,i])
  }
  
  #grid of lat/lon values
  grid = make.surface.grid(list(x=seq(-127, -123.5, length=dim(defDat$dat)[1]), 
                                y=seq(49.5, 39.5, length=dim(defDat$dat)[2])))
  
  #get mask of grid cells deformed by CSZ quake
  datmax = apply(abs(dat), 2, max)
  m0 = datmax > .5
  m1 = (grid[,1] >= -126) & (grid[,1] <= -125) & (grid[,2] >= 48)
  m2 = (grid[,1] >= -125) & (grid[,1] <= -124) & (grid[,2] >= 46) & (grid[,2] <= 48.3)
  m3 = (grid[,1] >= -124.8) & (grid[,1] <= -124.2) & (grid[,2] >= 40.3) & (grid[,2] <= 46)
  mask = m0 | m1 | m2 | m3
  
  maskdat = dat[,mask]
  #magnitudes (not sure about the 8.6's.  They correspond to the SS1 and SS3 cases)
  Mw = c(9.0, 9.1, 9.0, 8.9, 9.0, 8.9, 8.7, 8.8, 8.7, 8.6, 8.6, 9.1, 9.2, 9.1, 9.1, 9.2, 9.1, 8.6, 8.6)
  Mw10 = 10^Mw
  #type of earthquake (splay fault, deep)
  type = c("Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Deep", "Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Deep")
  type = as.factor(type)
  #probability weights of each observation (assume SS is half of equivalent SM probs)
  probs = c(.128, .016, .016, .318, .106, .106, .104, .078, .078, .052, .039, 
            .02, .0025, .0025, .02, .0025, .0025, .052, .039)
  
  #remove SS earthquakes
  SS=c(10, 11, 18, 19)
  dat2 = dat[-SS,]
  maskdat2 = maskdat[-SS,]
  Mw2 = Mw[-SS]
  Mw102 = Mw10[-SS]
  type2 = type[-SS]
  probs2 = probs[-SS]
  
  #repeat variables for each grid cell and each replication as necessary
  datPerRep2 = ncol(maskdat2)
  nRep2 = nrow(maskdat2)
  
  Mw10Fin2 = rep(Mw102, datPerRep2)
  typeFin2 = rep(type2, datPerRep2)
  # gridCellFin = rep(gridCell, rep(nRep, datPerRep))
  
  #regress
  #mod1 = lm(c(maskdat) ~ Mw10Fin + typeFin + gridCellFin)
  #NOTE: excluding interaction term leaves strange artifact from splay fault discontinuities
  mod2 = lm(maskdat2 ~ Mw102*type2)
  
  #get residuals
  testIn = data.frame(Mw102=Mw102, type2=type2)
  testOut = predict(mod2, testIn)
  preds = array(0, dim=dim(dat2))
  preds[,mask] = testOut
  resids = preds - dat2
  resids = t(resids)
  resids = array(resids, dim=dim(defDat$dat[,,-SS]))
  
  out = list(resids=resids, mod=mod2, Mw=Mw2, type=type2, probs=probs2, mask=mask, SS=SS)
  
  #return results
  return(out)
}

getPctResDat = function(defDat) {
  
  #transform to one column of data for each realization
  dat = matrix(nrow=dim(defDat$dat)[3], ncol=dim(defDat$dat)[1]*dim(defDat$dat)[2])
  for(i in 1:dim(defDat$dat)[3]) {
    dat[i,] = c(defDat$dat[,,i])
  }
  
  #grid of lat/lon values
  grid = make.surface.grid(list(x=seq(-127, -123.5, length=dim(defDat$dat)[1]), 
                                y=seq(49.5, 39.5, length=dim(defDat$dat)[2])))
  
  #get mask of grid cells deformed by CSZ quake
  datmax = apply(abs(dat), 2, max)
  m0 = datmax > .5
  m1 = (grid[,1] >= -126) & (grid[,1] <= -125) & (grid[,2] >= 48)
  m2 = (grid[,1] >= -125) & (grid[,1] <= -124) & (grid[,2] >= 46) & (grid[,2] <= 48.3)
  m3 = (grid[,1] >= -124.8) & (grid[,1] <= -124.2) & (grid[,2] >= 40.3) & (grid[,2] <= 46)
  mask = m0 | m1 | m2 | m3
  
  maskdat = dat[,mask]
  #magnitudes (not sure about the 8.6's.  They correspond to the SS1 and SS3 cases)
  Mw = c(9.0, 9.1, 9.0, 8.9, 9.0, 8.9, 8.7, 8.8, 8.7, 8.6, 8.6, 9.1, 9.2, 9.1, 9.1, 9.2, 9.1, 8.6, 8.6)
  Mw10 = 10^Mw
  #type of earthquake (splay fault, deep)
  type = c("Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Deep", "Splay", "Shallow", "Deep", "Splay", "Shallow", "Deep", "Splay", "Deep")
  type = as.factor(type)
  #probability weights of each observation (assume SS is half of equivalent SM probs)
  probs = c(.128, .016, .016, .318, .106, .106, .104, .078, .078, .052, .039, 
            .02, .0025, .0025, .02, .0025, .0025, .052, .039)
  
  #remove SS earthquakes
  SS=c(10, 11, 18, 19)
  dat2 = dat[-SS,]
  maskdat2 = maskdat[-SS,]
  Mw2 = Mw[-SS]
  Mw102 = Mw10[-SS]
  type2 = type[-SS]
  probs2 = probs[-SS]
  
  #repeat variables for each grid cell and each replication as necessary
  datPerRep2 = ncol(maskdat2)
  nRep2 = nrow(maskdat2)
  
  Mw10Fin2 = rep(Mw102, datPerRep2)
  typeFin2 = rep(type2, datPerRep2)
  # gridCellFin = rep(gridCell, rep(nRep, datPerRep))
  
  #regress
  #mod1 = lm(c(maskdat) ~ Mw10Fin + typeFin + gridCellFin)
  #NOTE: excluding interaction term leaves strange artifact from splay fault discontinuities
  mod2 = lm(maskdat2 ~ Mw102*type2)
  
  #get residuals
  testIn = data.frame(Mw102=Mw102, type2=type2)
  testOut = predict(mod2, testIn)
  preds = array(0, dim=dim(dat2))
  preds[,mask] = testOut
  resids = preds - dat2
  
  #convert to pct residuals
  pctRes = resids/dat2
  pctRes[dat2 == 0] = 0 #both resDat and defDat are 0, so set error to 0
  pctRes[,!mask] = NA # everything outside of mask has infinite pct error
  
  #set dim=c(nx, ny, 15)
  pctRes = t(pctRes)
  pctRes = array(pctRes, dim=dim(defDat$dat[,,-SS]))
  
  out = list(pctRes=pctRes, mod=mod2, Mw=Mw2, type=type2, probs=probs2, mask=mask, SS=SS)
  
  #return results
  return(out)
}