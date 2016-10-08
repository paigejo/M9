#same as stationary.cov from fields, except takes into account data being 
#from multiple realizations, where data from different realizations are 
#considered independent.  Creates sparse block diagonal matrix from 
#Matrix class
stationary.cov.multi = function(x1, realizations=NULL, ...) {
  uniqueReals = unique(realizations)
  
  # if only one realization, return result of stationary.cov
  if(is.null(realizations) || length(unique(realizations)) == 1) {
    return(do.call("stationary.cov", c(list(x1=x1), list(...))))
  }
  
  #make sure we have matrices
  if(!is.matrix(x1)) {
    x1 = matrix(x1, nrow=length(x1))
  }
  
  #store covariance matrices for each realization in a list
  Sigmas = list()
  for(i in 1:length(uniqueReals)) {
    #make covariance matrix for this realization
    real = uniqueReals[i]
    
    thisReal = realizations == real
    thisx1 = x1[thisReal,]
    Sigmas[[i]] = do.call("stationary.cov", c(list(x1=thisx1), list(...)))
  }
  
  #return sparse block diagonal covariance matrix
  return(bdiag(Sigmas))
}


