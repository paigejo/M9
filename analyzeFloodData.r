setwd("/Users/paigejo/Google Drive/UW/guttorp/code")
library(RcppCNPy)
library(fields)

#get max flood heights
allHMax = array(NA, dim=c(6, 250, 250))
for(i in 0:5) {
  fname = paste0("allHMax", i, ".npy")
  tmp = npyLoad(fname)
  allHMax[i+1, , ] = tmp
}

alpha = .1
contVal = 5

getContourSignificances = function(contVal) {
  getContPercAtLoc = function(floodLevels) {
    fun = ecdf(floodLevels)
    fun(contVal)
  }
  
  contPercentiles = apply(allHMax, c(2, 3), getContPercAtLoc)
  
  #get significances on .5-1 scale
  contSigs = 1 - abs(.5 - contPercentiles)
  
  #get significances on 0-1 scale
  contSigs = contSigs*2 - 1
  
  return(list(contPercentiles=contPercentiles, contSigs = contSigs))
}

out = getContourSignificances(contVal)
contPercentiles = out$contPercentiles
contSigs = out$contSigs

xgrid = 1:250
ygrid = 1:250
gridList = make.surface.grid(list(x=xgrid, y=ygrid))
sigs = as.numeric(contSigs)
goodSigs = sigs > alpha
percs = as.numeric(contPercentiles)
quilt.plot(gridList[goodSigs,], sigs[goodSigs], 
           main=paste0(contVal, "m Contour Significances"), zlim=c(0,1))
quilt.plot(gridList, percs, 
           main=paste0(contVal, "m Contour Percentiles"), zlim=c(0,1))


##################
#####do excursions
##################

u = contVal
mu = apply(allHMax, 1, mean)
type="="

# tps <- fastTps(gridList, ys, theta=5)
# tpsInterp <- predictSurface(tps, nx=300, ny=300)
# surface(tpsInterp, main=paste0(contVal, "m contour significances"), 
#         xlab='longitude', ylab='latitude')

#fit exponential covariance to data




#example
n = 21
Q.x = sparseMatrix(i=c(1:n, 2:n),
                   j=c(1:n, 1:(n-1)),
                   x=c(rep(1, n), rep(-0.1, n-1)),
                   dims=c(n, n),
                   symmetric=TRUE)

## Set the mean value function
mu.x = seq(-5, 5, length=n)

## calculate the level 0 positive excursion function
res.x = excursions(alpha=1, u=0, mu=mu.x, Q=Q.x, type='>',
                   verbose=1, max.threads=2)

## Plot the excursion function and the marginal excursion probabilities
plot(res.x$F, type="l",
     main='Excursion function (black) and marginal probabilites (red)')
lines(res.x$rho, col=2)