source("loadTestData.r")

#generate data matrix from deformation data
defDat = loadDeformations()
dat = matrix(nrow=dim(defDat$dat)[3], ncol=dim(defDat$dat)[1]*dim(defDat$dat)[2])
for(i in 1:dim(defDat$dat)[3]) {
  dat[i,] = c(defDat$dat[,,i])
}

#take truncated SVD (only need V since t(X)%*%X = V S^2 V* for U,S,V result of svd(X))
truncation = 40
out = svd(dat, nu=0, nv=truncation)

save(out, file="seaDef_svd.RData")

