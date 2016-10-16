# Stan backprojection model (with lognormal prior on beta)
# see http://www.uvm.edu/~bbeckage/Teaching/DataAnalysis/Manuals/stan-reference-2.8.0.pdf
# for help
data {
  int<lower=0> n; // number observations (~500 for full dataset)
  int<lower=0> pAreal; // number of areal zeta coefficients (240)
  int<lower=0> pPoint; // number of point zeta values (631)
  
  matrix[n,pAreal] X; // design matrix (G %*% T)
  vector[n] y;   // observations
  vector[n] sigmaY; // standard deviations of y conditional on zeta
  
  // these are the lognormal canonical parameters (parameters for MVN to be exponentiated).
  // In other words, these are the parameters for the marginal distribution of log(zeta)
  vector[pPoint+pAreal] priorMu;     // mean vector
  matrix[pPoint+pAreal,pPoint+pAreal] priorSigmaL;     // covariance matrix Cholesky decomposition
}
transformed data {
  int<lower=0> pTot;
  pTot = pPoint + pAreal;
}
parameters {
  // reparameterization of the regression coefficients
  vector[pTot] alpha;
}
transformed parameters {
  // transform from uncorrelated standard normals to MVNs and MVLNs
  // get values for both CSZ grid cell averages and pointwise values
  vector[pTot] logZetaAll;
  vector[pAreal] beta;
  vector[pPoint] logZetaPoint;
  vector<lower=0>[pAreal] zeta;
  vector<lower=0>[pPoint] zetaPoint;
  logZetaAll = priorMu + priorSigmaL * alpha;
  beta = logZetaAll[(pPoint+1):pTot];
  logZetaPoint = logZetaAll[1:pPoint];
  zeta = exp(beta);
  zetaPoint = exp(logZetaPoint);
}
model {
  // vectorized prior
  alpha ~ normal(0, 1);
  // implies beta ~ multi_normal_cholesky(priorMu, priorSigmaL)
  
  // data model (likelihood)
  y ~ normal(X * zeta, sigmaY);
}