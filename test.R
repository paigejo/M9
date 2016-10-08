# function for testing normality assumption of G T zeta
testNormality = function(params, fault=faultGeom, nsim=10000) {
  # draw from marginal distribution of slip
  out = preds(params, fault=fault, nsim=nsim)
  
  # compute subsidences
  out = predsToSubsidence(params, out, fault=fault)
  subSims = out$subSims
  
  # make histogram and add normal approximation
  hist(subSims, freq=F, breaks=550, xlim=c(-2,1), main="Simulated Marginal Subsidences")
  subMu = mean(subSims)
  subSigma = sd(subSims)
  xs = seq(-2, 1, l=200)
  lines(xs, dnorm(xs, subMu, subSigma), col="red")
  # hist(subSims, freq=F, breaks=150, xlim=c(-3, 2))
}
# NOTE: For the 20 or 240 size faults, definitely not normal.  
#       True distribution is much more peaked in the center and heavy tailed




