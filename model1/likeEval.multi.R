#evaluates Gaussian process likelihood for data locations x1 corresponding 
#to realizations given as input using stationary.cov.multi function.  
#Additional arguments can be passed as arguments to stationary.cov.multi 
#and from that to stationary.cov.R
likeEval.multi = function(x1, y, realizations, sill, nugget, ...) {
  #get covariance matrix
  Sigma = do.call("stationary.cov.multi", c(list(x1=x1, realizations=realizations), list(...)))
  Sigma = Sigma*sill + diag(nugget, nrow=length(y))
  
  #do Cholesky decomposition
  U = chol(Sigma)
  
  #get log of determinant component of log likelihood
  lik.logDetSigma = 2*sum(log(diag(U)))
  
  #get quadratic form component of log likelihood (y^T Sigma^(-1) y)
  u = forwardsolve(U, transpose=TRUE, y, upper.tri=TRUE)
  lik.quadForm = sum(u^2)
  
  #get constant term of log likelihood
  n = length(y)
  lik.const = -(n/2)*log(2*pi)
  
  #return final log likelihood
  return(lik.const - (1/2)*lik.quadForm - (1/2)*lik.logDetSigma)
}