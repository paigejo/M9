# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
mKrig.MLE.multi <- function(x, y, realizations=rep(1, nrow(x)), ...) {
  #NOTE: lambda.profile argument for mKrig.MLE is assumed to be FALSE here
  #along with verbose
  
  #get number of realizations
  uniqueReals = unique(realizations)
  numReals = length(uniqueReals)
  
  #call mKrig.MLE on each realization
  outList = list()
  summary = c()
  for(r in 1:length(uniqueReals)) {
    xr = x[realizations==r,]
    yr = y[realizations==r]
    
    outList[[r]] = do.call("mKrig.MLE", c(list(x=xr, y=yr, lambda.profile=FALSE, verbose=FALSE), list(...)))
    if(r == 1)
      summary = outList[[r]]$summary
    else {
      tmp = outList[[r]]$summary
      #update EffDf
      summary[,1] = summary[,1] + tmp[,1]
      #update lnProfLik
      summary[,2] = summary[,2] + tmp[,2]
      #update GCV (take average GCV over the replications)
      summary[,3] = (summary[,3]*(r-1) + tmp[,3])/r
    }
    
    outList[[r]] = c(outList[[r]], list(realization=uniqueReals[r]))
  }
  
  #combine results
  par.grid = outList[[1]]$par.grid
  maxRow = which.max(summary[,2])
  cov.args.MLE = as.list(par.grid[maxRow,])
  lambda.best <- lambda.MLE <- exp(summary[maxRow,6])
  
  #these should have the same value or will be NULL for all calls
  # to mKrig.MLE
  mKrig.args = outList[[1]]$mKrig.args
  call = outList[[1]]$call
  lnLike.eval = outList[[1]]$lnLike.eval
  
  #create final output object
  out = list(summary=summary, par.grid=par.grid, cov.args.MLE=cov.args.MLE, 
             mKrig.args=mKrig.args, lambda.best=lambda.best, lambda.MLE=lambda.MLE, 
             call=call, lnLike.eval=lnLike.eval)
  return(out)
}