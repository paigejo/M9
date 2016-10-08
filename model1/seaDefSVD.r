source("loadTestData.r")
library(rARPACK)

#generate data matrix from deformation data
defDat = loadAllDeformations()
dat = matrix(nrow=dim(defDat$dat)[3], ncol=dim(defDat$dat)[1]*dim(defDat$dat)[2])
for(i in 1:dim(defDat$dat)[3]) {
  dat[i,] = c(defDat$dat[,,i])
}

# plot the deformation data
grid = make.surface.grid(list(x=seq(-127, -123.5, length=dim(defDat$dat)[1]), 
                              y=seq(49.5, 39.5, length=dim(defDat$dat)[2])))
for(i in 1:nrow(dat)) {
  quilt.plot(grid, c(dat[i,]), main=paste0(i, "th Realization"))
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
}

# now redo the plots on a shared z (color) scale
minDef= min(dat)
maxDef = max(dat)
for(i in 1:nrow(dat)) {
  pdf(file=paste0("realization", i, ".pdf"), width=4, height=5)
  quilt.plot(grid, c(dat[i,]), main=paste0(i, "th Realization Uplift (m)"), zlim=c(minDef, maxDef), 
             xlab="Longitude", ylab="Laittude")
  map("world", "Canada", add=TRUE, lwd=1.5)
  US(add=TRUE, lwd=1.5)
  dev.off()
}

#get center columns of dataMat by their means
mu = colMeans(dat)
datCntr = sweep(dat, 2, mu)

#take truncated SVD (only need V since t(X)%*%X = V S^2 V* for U,S,V result of svd(X))
truncation = 19
out = svds(datCntr, k=truncation)

#find total variance: square Frobenius norm of X
totalVar = norm(datCntr, type="F")^2

#####save results
save(out, totalVar, mu, file="seaDef_svd.RData")

#####plots
#plot variance explained
plot(1:length(out$d), out$d^2, pch=19, col="blue", 
     main="Deformation Variance Explained", xlab="", ylab="Singular Value")

#plot cumulative variance explained
plot(cumsum(out$d^2)/sum(out$d^2), xlab="Eigenvector", ylab="Variance Explained (Fraction)", 
     main="Cumulative Fraction Variance Explained", pch=19, col="blue")

#plot eigenmodes
for(i in 1:length(out$d)) {
  quilt.plot(grid, c(out$v[,i]), main=paste0(i, "th Eigenmode"))
}
#note that 1-9 appear to have little noise

print.matrix <- function(m){
  paste0(c("c(", paste(as.character(m), collapse=", "), ")"), collapse="")
}

##### train statistical emulator

scores = out$u %*% diag(out$d)
#these are the indices of deformation realization corresponding to test flood levels
testI = c(1, 4, 7, 12, 15, 19)
truncation = 3 #eigenmodes 1-9 appear to have little noise, first 3 capture almost 98% of variation
testScores = scores[testI, 1:truncation]
scoreMean = apply(testScores, 2, mean) #chacteristic scores
normScores = sweep(testScores, 2, scoreMean)
scoreSD = sd(normScores)
normScores = normScores*(1/scoreSD)
#assume correlation scale is inversely proportional to variance explained?
#Not yet, must explore. Not enough data to explore yet, though.

floodDat = loadFloodData()
allHMax = floodDat$allHMax
top = floodDat$topo
lon = floodDat$lon
lat = floodDat$lat
maxFlood = apply(allHMax, 1, max)

#Are the max flood levels smooth?  If so, can use only quadratic poly regressors
quilt.plot(normScores[,1], normScores[,2], maxFlood, xlab="1st Eigenvector Score", 
           ylab="2nd Eigenvector Score", nx=20, ny=20, main="Max Flood Level (m)")
quilt.plot(normScores[,1], normScores[,3], maxFlood, xlab="1st Eigenvector Score", 
           ylab="3rd Eigenvector Score", nx=20, ny=20, main="Max Flood Level (m)")
quilt.plot(normScores[,2], normScores[,3], maxFlood, xlab="2nd Eigenvector Score", 
           ylab="3rd Eigenvector Score", nx=20, ny=20, main="Max Flood Level (m)")

# regress off elevation from allHMax
nX = dim(allHMax)[2]
nY = dim(allHMax)[3]
residHMax = allHMax
for(i in 1:dim(allHMax)[1]) {
  floodVals = allHMax[i,,]
  model = lm(c(floodVals) ~ c(top))
  residVals = residuals(model)
  residHMax[i,,] = array(c(residVals), dim=c(1,nX,nY))
}
residMean = apply(residHMax, c(2, 3), mean)
residCntr = sweep(residHMax, c(2, 3), residMean)
allHMaxMean = apply(allHMax, c(2, 3), mean)
allHMaxCntr = sweep(allHMax, c(2, 3), allHMaxMean)

######
###### plot allHMax and residHMax (centered and not centered by grid box mean)
floodGrid = matrix(c(lon, lat), ncol=2)
for(i in 1:nrow(residHMax)) {
  quilt.plot(floodGrid, c(allHMax[i,,]), main=paste0(i, "th Realization"))
}
for(i in 1:nrow(residHMax)) {
  quilt.plot(floodGrid, c(allHMaxCntr[i,,]), main=paste0(i, "th Realization"))
}
for(i in 1:nrow(residHMax)) { # for non-complex models, this is the best
  quilt.plot(floodGrid, c(residHMax[i,,]), main=paste0(i, "th Realization"))
}
for(i in 1:nrow(residHMax)) { #this is the best, but is the model to complex?
  quilt.plot(floodGrid, c(residCntr[i,,]), main=paste0(i, "th Realization"))
}

