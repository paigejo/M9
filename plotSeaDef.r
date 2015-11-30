source("loadTestData.r")

#generate data matrix from deformation data
defDat = loadDeformations()
dat = matrix(nrow=dim(defDat$dat)[3], ncol=dim(defDat$dat)[1]*dim(defDat$dat)[2])
for(i in 1:dim(defDat$dat)[3]) {
  dat[i,] = c(defDat$dat[,,i])
}

#plot data
grid = make.surface.grid(list(x=1:dim(defDat$dat)[1], y=1:dim(defDat$dat)[2]))
for(i in 1:nrow(dat)) {
  quilt.plot(grid, c(dat[i,]), main=paste0(i, "th Realization"))
}

tmp = load("seaDef_svd.RData")
S = out$d
V = out$v

#plot variance explained
plot(1:length(S), S^2, pch=19, col="blue", 
     main="Deformation Variance Explained", xlab="", ylab="Singular Value")

#plot cumulative variance explained
thresh = sum(S^2)*.95
cumVarExp = cumsum(S^2)
plot(1:length(S), cumVarExp/sum(S^2), pch=19, col="blue", ylim=c(0, 1), 
     main="Percent of Total Deformation Variance Explained", xlab="", 
     ylab="Singular Value")

#plot eigenmodes
for(i in 1:length(S)) {
  quilt.plot(grid, c(V[,i]), main=paste0(i, "th Eigenmode"))
}









