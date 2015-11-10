setwd("/Users/paigejo/Google Drive/UW/guttorp/code")
source("loadTestData.r")
library(fields)

#generate vgram data
R = 3959 #in miles
lonExtent= c(235.79781, 235.82087)
CCLon = mean(lonExtent)
CCLat = mean(latExtent)
latExtent = c(41.739671,41.762726)
lonDist = cos(lat)*2*pi/360*R
latDist = 2*pi*R/360

xgrid = lon*lonDist
ygrid = lat*LatDist
grid = make.surface.grid(list(X=xgrid, Y=ygrid))
X = grid$X
Y = grid$Y

maxX = max(X)
minX = min(X)
distPerCell = (maxX - minX)/250
rectDim = 15
maxDist = floor(rectDim/2)*distPerCell

#calculate vgram for each slice of HMax: each tsunami realization
for(i in 1:dim(allHMax)[1]) {
  if(i == 1)
    VG= vgram(grid, allHMax[i,,], dmax=maxDist)
  else {
    #concetenate vgram
    VGslice = vgram(grid, allHMax[i,,], dmax=maxDist)
    VG$d = c(VG$d, VGslice$d)
    VG$vgram = c(VG$vgram, VGslice$vgram)
  }
}

#fit exponential variogram to data
s = mean(VG$vgram[VG$d > maxDist*.9]) #sill
n = mean(VG$vgram[VG$d < maxDist*.1]) #nugget
r = maxDist # range
fit = nls(VG$vgram ~ (s - n)*(1 - exp(-VG$d/r)) + n, start=list(s=s, n=n, r=r))
summary(fit)

#get variogram coefficients
coefs = coef(fit)
s = coefs$s
n = coefs$n
r = coefs$r

#plot variogram fit
expVGram = function(h) {
  (s - n)*(1 - exp(-h/r)) + n
}
xs = seq(0, maxDist, length=500)
vgramFit = expVGram(xs)

pdf("expVGram.pdf", height=5, width=7)
plot(VG, main="Empirical and Exponential Variogram Fit")
lines(xs, expVGram(xs), col=green)
dev.off()

save(s, n, r, file="fitExpVarioParams.RData")