# Stan backprojection model (with lognormal prior on beta)
# see http://www.uvm.edu/~bbeckage/Teaching/DataAnalysis/Manuals/stan-reference-2.8.0.pdf
# for help
data {
  int<lower=0> n; // number observations (~500)
  int<lower=0> p; // number of coefficients (20)
  
  matrix[n,p] X; // design matrix
  vector[n] y;   // observations
  
  real<lower=0> priorMu;     // lognormal mean for prior on beta (computed based on Golfinger et al. 2012 data)
  real<lower=0> priorSigma;     // lognormal mean for prior on beta (computed based on Golfinger et al. 2012 data)
}
parameters {
  // regression coefficients
  vector<lower=0>[p] beta;
  
  // NOTE: why isn't stdev of y closer to 1?
  real<lower=0> sigma_y;
}
model {
  // vectorized prior
  beta ~ lognormal(priorMu,priorSigma);
  
  // data model
  y ~ normal(X * beta, sigma_y);
}
generated quantities {
  vector[p] logBeta;
  logBeta = log(beta);
}