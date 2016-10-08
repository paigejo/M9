source('~/git/M9/seaDefAnis.r')
source('~/git/M9/loadTestData.r')
source('~/git/M9/predictFields.r')
source('~/git/M9/regressMask.r')
seaDef = loadAllDeformations()
predHull = makePredHull(seaDef, 30000)
mod = fitAnisModel(seaDef, 2000)
test = krigeModel(seaDef, mod, N = 500, numPredictions = 11, predHull=predHull)

model= mod
model$phi = model$phi*2
preds = predictFields(seaDef, model, N=500, M=1, predHull=predHull)

krigMod = preds$krigMod

obj = fitFields(seaDef, model, N=500, predHull=predHull)
obj1 = obj$obj1
obj2 = obj$obj2
obj2$summary
obj2$cov.args.MLE

mu = mean(seaDef$dat)
obj = fitFields.multi(seaDef, model, N=2000, predHull=predHull, m=0)
obj$summary
obj$cov.args.MLE

preds = predictFields(seaDef, krigeMLE=obj, N=500, M=1, predHull=predHull, real=4)

# remove outliers from each realization in seaDef 
# (magnitude of residuals above 99th percentile of each realization)
source('~/git/M9/seaDefAnis.r')
source('~/git/M9/loadTestData.r')
source('~/git/M9/predictFields.r')
source('~/git/M9/regressMask.r')
source('~/git/M9/mKrig.MLE.multi.R')
source('~/git/M9/stationary.cov.multi.R')
load("~/git/M9/predHull.RData")

seaDef = loadAllDeformations()
# predHull = makePredHull(seaDef, 30000)
tmpPredHull = predHull
tmpPredHull$lon = tmpPredHull$lon - 360
plotPredHull(seaDef, predHull, toMeters=FALSE)

resids = getPctResDat(seaDef)
pctRes = resids$pctRes
outThreshes = apply(abs(pctRes), 3, quantile, probs=.99, na.rm=TRUE)
for(r in 1:dim(pctRes)[3]) {
  thresh = outThreshes[r]
  thisPctRes = pctRes[,,r]
  thisPctRes[abs(thisPctRes) > thresh] = NA
  pctRes[,,r] = thisPctRes
}
mu = mean(pctRes, na.rm=TRUE)
pctRes = pctRes - mu

mod = fitAnisModel(seaDef, pctRes=pctRes, N=3000)
model= mod
model$phi
model$sigmasq
#note: model is only used here for matern smoothness
test = fitFields.multi(seaDef=seaDef, model=model, pctRes=pctRes, 
                      N=1000, predHull=predHull)
obj1 = fitFields(seaDef=seaDef, model=model, pctRes=pctRes, 
                      N=1000, predHull=predHull, real=4)

#obj$summary
obj$cov.args.MLE
obj$lambda.MLE

preds = predictFields(seaDef, krigeMLE=obj, N=500, M=3, predHull=predHull, 
                      real=4, mu=mu, expand=3)
out = preds$krigMod
out$sigma.MLE
out$rho.MLE

#test predictions for different parameters
test = obj
test$lambda.MLE = .01
# test$cov.args.MLE$phi = 300000
test$cov.args.MLE[[1]] = 400000
preds = predictFields(seaDef, krigeMLE=test, N=500, M=3, predHull=predHull, 
                      real=4, mu=mu, expand=15)

preds = myPredictFields(seaDef, krigeMLE=test, N=500, M=3, predHull=predHull, 
                      real=1, mu=mu)

preds = myPredictFields2(seaDef, model=model, N=500, M=3, predHull=predHull, 
                        real=4, mu=mu)

# find out range of each of the realizations' data:
apply(pctRes, 3, range, na.rm=TRUE)
range(pctRes, na.rm=TRUE)

tmp = resids$pctRes
range(tmp, na.rm=TRUE)

sampleNoNA = function(dat, size=1) {
  samps = NA
  while(any(is.na(samps))) {
    samps = sample(dat, 1)
  }
  return(samps)
}
tmp = apply(pctRes, 3, sampleNoNA, size=1)
tmpMod = model
var(tmp)*(14/15) #get MLE
sqrt(var(tmp)*(14/15))
tmpMod$sigmasq = var(tmp)*(14/15) #get MLE
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=1, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=2, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=3, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=4, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=5, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=6, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=7, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=8, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=9, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=10, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=11, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=12, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=13, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=14, mu=mu)
preds = myPredictFields2(seaDef, model=tmpMod, N=500, M=3, predHull=predHull, 
                         real=15, mu=mu)

#variogram methods:
source('~/git/M9/residsVGram.R')
source('~/git/M9/fitVario.r')
source('~/git/M9/likeEval.multi.R')
source('~/git/M9/stationary.cov.multi.R')
library(fields)
library(Matrix)

#get variogram
VG = getVGram(seaDef, pctRes, N=50000)

#make plots showing model has very poor fit
plotVGram(VG, kappa=3, phi=mod$phi, sill=mod$sigmasq, nugget=mod$nugget, ylim=c(0, .1))
plotVGram(VG, kappa=3, phi=mod$phi, sill=mod$sigmasq, nugget=mod$nugget, ylim=c(0, 1))

#####fit parameters for Matern variogram function
#these would be initial guesses for Matern parameter optimization
sillInit=sd(apply(pctRes, 3, median, na.rm=TRUE))
nuggetInit = mean(VG$vgram[VG$d < 50000])
rangeInit=mod$phi
sillInit
nuggetInit
rangeInit
# nls function is having issues
# fit = fitVGram(VG, sillInit=sillInit, rangeInit = 300000, nuggetInit=nuggetInit, kappa=3)

#change Matern function from fields to variance rather than covariance 
#and include nugget and sill
f = function(d, phi=rangeInit, kap=3, s=sillInit, n=nuggetInit) {
  #if parameters don't make sense, return bad number
  if(phi <= 0 || s < n || n < 0) {
    return(rep(-10, length(d)))
  }
  d = d/phi
  s - (Matern(d, nu=kap, phi=(s - n)) + n)
}

# plot fitted functions versus empirical variogram:
xs = seq(0, max(VG$d), length=500)
plotVGMean(VG, ylim=c(0, .1), 
           main="Empirical Matern Variogram and Fit", N=200, pch=19, cex=.4)
lines(xs, sqrt(f(xs, phi=300000, s=sillInit, n=nuggetInit)), col="green")

#do grid fit for Matern parameters
nugget.grid=seq(0, .08, length=9)
sill.grid=seq(0, .1, length=21)
range.grid=10^seq(log10(50000), log10(500000), length=15)
param.grid = make.surface.grid(list(phi=range.grid, s=sill.grid, n=nugget.grid))
outLS = gridLS(VG$d, VG$vgram, f, param.grid)
OLSE = outLS$OLSE
summTable = outLS$summTable

#plot grid search MSE over nugget and sill parameters fixing
#range at the range OLSE
OLSE.range = OLSE[1]
plotGridI = param.grid[,1] == OLSE.range
quilt.plot(param.grid[plotGridI,2:3], summTable[plotGridI,1], zlim=c(0, .005), 
           xlab="sill", ylab="nugget", main="MSE")
quilt.plot(param.grid[plotGridI,2:3], summTable[plotGridI,1], zlim=c(0, .001), 
           xlab="sill", ylab="nugget", main="MSE")
xs = seq(0, max(VG$d), length=500)
plotVGMean(VG, ylim=c(0, .1), 
           main="Empirical Matern Variogram and Fit", N=200, pch=19, cex=.4)
lines(xs, sqrt(f(xs, phi=OLSE.range, s=OLSE[2], n=OLSE[3])), col="green")

#generate predictions
varioMod = list(phi=OLSE.range, sigmasq=OLSE[2], tausq=OLSE[3], kappa=3)
preds = myPredictFields2(seaDef, pctRes, model=varioMod, N=500, M=3, predHull=predHull, 
                         real=1, mu=mu)

#update pctRes so that threshold is removed for final data comparison
pctResThresh = pctRes
pctRes = resids$pctRes - mu

plotSample(seaDef=seaDef, pctRes=pctRes, model=varioMod, N=1000, predHull=predHull, mu=mu)

#show that means of each realization pctRes differ by more than you might expect
hist(apply(pctResThresh, 3, median, na.rm=TRUE), breaks=20)
varioMod$sigmasq


##### Maximum likelihood variogram fitting:

#first get sample of points
N=3000
NPerReal = round(N/15)
samplePts = c()
sampleRes = c()
for(r in 1:15) {
  tmp = samplePctRes(seaDef, pctRes, N=NPerReal, real=r, thresh=Inf)
  samplePts = rbind(samplePts, tmp$gridXY)
  sampleRes = c(sampleRes, tmp$res)
}
sampleReals = rep(1:15, rep(NPerReal, 15))

#now grid fit for Matern parameters
nugget.grid=seq(0, .08, length=9)
sill.grid=seq(.001, .05, length=15)
range.grid=10^seq(log10(50000), log10(500000), length=15)
full.grid = make.surface.grid(list(sill=sill.grid, nugget=nugget.grid, theta=range.grid))
sill.grid=full.grid[,1]
nugget.grid=full.grid[,2]
range.grid=full.grid[,3]
param.grid=cbind(theta=range.grid)
outML = gridML(samplePts, sampleRes, realizations=sampleReals, sill=sill.grid, nugget=nugget.grid, 
               paramGridList=param.grid, Covariance="Matern", smoothness=3)
MLE = outML$MLE
summTableMLE = outML$summTable

#plot grid search log likelihood over nugget and sill parameters fixing
#range at the range MLE
MLE.range = MLE[3]
plotGridI = range.grid == MLE.range
quilt.plot(full.grid[plotGridI,1:2], summTableMLE[plotGridI,1], 
           xlab="sill", ylab="nugget", main="Log Likelihood")
plotGridI = sapply(range.grid, all.equal, 219698.53) == TRUE
quilt.plot(full.grid[plotGridI,1:2], summTableMLE[plotGridI,1], 
           xlab="sill", ylab="nugget", main="Log Likelihood")
xs = seq(0, max(VG$d), length=500)
plotVGMean(VG, ylim=c(0, .1), 
           main="Empirical Matern Variogram and Fit (maximum Likelihood)", N=200, pch=19, cex=.4)
lines(xs, sqrt(f(xs, phi=MLE.range, s=MLE[1], n=MLE[2])), col="green")
