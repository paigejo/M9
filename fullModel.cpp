// fit combined model allowing taper, gamma, and standard deviation to vary with latitude
#include <TMB.hpp>

template<class Type>
// define taper function
Type taper(Type depth, Type lambda, Type dStar) {
  if(sqrt(pow(lambda, Type(2))) < Type(0.0000005)) {
    // this is a limit as lambda approaches 0. Coding this explicitly avoids numerical instability
    return Type(1) - pow(depth/dStar, Type(2));
  } else if(depth > dStar) {
    return Type(0);
  } else if(depth <= Type(0)) {
    return Type(1);
  } else {
    Type scaledDepth = pow(sqrt(pow(depth/dStar, Type(2))), Type(2)) * pow(lambda, Type(2));
    Type ans = Type(1) - (Type(1) - exp(-scaledDepth))/(Type(1) - exp(-pow(lambda, Type(2))));
    return ans;
  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(x);
  DATA_VECTOR(xsd);
  
  DATA_VECTOR(y);
  DATA_VECTOR(ysd);
  DATA_VECTOR(lowI);
  
  DATA_MATRIX(G);
  DATA_MATRIX(zeroMask);
  DATA_MATRIX(DSStrikeCSZ);
  DATA_MATRIX(DSDipCSZ);
  DATA_MATRIX(DSStrikeGPS);
  DATA_MATRIX(DSDipGPS);
  DATA_MATRIX(sdBasisX);
  DATA_MATRIX(sdBasisY);
  DATA_MATRIX(lambdaBasisX);
  DATA_MATRIX(lambdaBasisY);
  DATA_MATRIX(gammaBasis);
  
  DATA_VECTOR(faultDepths);
  DATA_VECTOR(xDepths);
  
  DATA_SCALAR(dStar);
  DATA_SCALAR(dStarGPS);
  
  PARAMETER(logmu);
  PARAMETER_VECTOR(betaTaper);
  PARAMETER_VECTOR(betasd);
  PARAMETER_VECTOR(betaGamma);
  PARAMETER(logphi);
  PARAMETER(logalpha);
  PARAMETER(loglowInflate);
  PARAMETER(loghighInflate);
  
  int nx = x.size();
  int ny = y.size();
  int nTaper = betaTaper.size();
  int nsd = betasd.size();
  int nGamma = betaGamma.size();
  
  // Outline: 
  // reparameterize
  // compute taper, standard deviation, and gamma vectors from the basis matrices
  // calculate mean vectors for both data sets
  // construct covariance matrices for zeta and slips over the fault and locking rates
  // calculate covariance matrices for the subsidence and locking rate data
  // calculate negative log likelihoods
  
  // reparameterize the scalars (although non-scalars must be reparameterized later)
  Type mu = exp(logmu);
  Type phi = exp(logphi);
  Type alpha = exp(logalpha);
  Type lowInflate = exp(loglowInflate);
  Type highInflate = exp(loghighInflate);
  
  // construct lambda, standard deviation, and gamma vectors for both fault and locking rate data
  vector<Type> lambdaVecX = lambdaBasisX * betaTaper;
  vector<Type> lambdaVecY = lambdaBasisY * betaTaper;
  vector<Type> sdVecX = sdBasisX * betasd;
  vector<Type> sdVecY = sdBasisY * betasd;
  vector<Type> gammaVec = gammaBasis * betaGamma;
  
  // reparameterize standard deviation and gamma
  for(int i=0; i<sdVecX.size(); i++)
    sdVecX(i) = exp(sdVecX(i));
  for(int i=0; i<sdVecY.size(); i++)
    sdVecY(i) = exp(sdVecY(i));
  for(int i=0; i<gammaVec.size(); i++)
    gammaVec(i) = exp(gammaVec(i));
  
  // construct taper values
  vector<Type> taperX(nx);
  REPORT(mu)
  for(int i=0; i<lambdaVecX.size(); i++)
    taperX(i) = taper(xDepths(i), lambdaVecX(i), dStarGPS);
  
  vector<Type> taperFault(lambdaVecY.size());
  for(int i=0; i<lambdaVecY.size(); i++)
    taperFault(i) = taper(faultDepths(i), lambdaVecY(i), dStar);
  
  // calculate mean vectors for the subsidence and locking rate data
  // vector<Type> muY = mu * G * taperFault;
  // vector<Type> muX = mu * gammaVec * taperX;
  vector<Type> muY = G * taperFault;
  vector<Type> muX(lambdaVecX.size());
  for(int i=0; i<muY.size(); i++)
    muY(i) = muY(i) * mu;
  for(int i=0; i<muX.size(); i++)
    muX(i) = gammaVec(i) * mu * taperX(i);
  
  // calculate anisotropic covariance matrices of zeta for the fault and the locking rate data
  Type kappa = Type(3 / 2);
  Type alphaSq = pow(alpha, Type(2));
  Type alphaInvSq = Type(1) / pow(alpha, Type(2));
  
  matrix<Type> SigmaZetaCSZ(DSStrikeCSZ);
  for(int i=0; i<SigmaZetaCSZ.rows(); i++)
    for(int j=0; j<SigmaZetaCSZ.cols(); j++)
      SigmaZetaCSZ(i,j) = sdVecY(i) * sdVecY(j) * matern(sqrt(alphaSq * DSStrikeCSZ(i,j) + alphaInvSq * DSDipCSZ(i,j)), phi, kappa);
  
  matrix<Type> SigmaZetaGPS(DSStrikeGPS);
  for(int i=0; i<SigmaZetaGPS.rows(); i++)
    for(int j=0; j<SigmaZetaGPS.cols(); j++)
      SigmaZetaGPS(i,j) = sdVecX(i) * sdVecX(j) * matern(sqrt(alphaSq * DSStrikeGPS(i,j) + alphaInvSq * DSDipGPS(i,j)), phi, kappa);
  
  // calculate covariance matrix of slips for fault
  matrix<Type> SigmaSlipFault(SigmaZetaCSZ);
  for(int i=0; i<SigmaZetaCSZ.rows(); i++)
    for(int j=0; j<SigmaZetaCSZ.cols(); j++)
      SigmaSlipFault(i,j) = taperFault(i) * taperFault(j) * SigmaZetaCSZ(i,j);
  
  ///// calculate covariance matrices for the subsidence and locking rate data
  // inflate residual covariance for the subsidence data
  vector<Type> yInflateSD(ny);
  for(int i=0; i<ysd.size(); i++) {
    if(lowI(i) == Type(1))
      yInflateSD(i) = ysd(i) * lowInflate;
    else
      yInflateSD(i) = ysd(i) * highInflate;
  }
  
  // calculate covariance matrix of locking rate data
  matrix<Type> SigmaXi(SigmaZetaGPS);
  for(int i=0; i<SigmaXi.rows(); i++)
    for(int j=0; j<SigmaXi.cols(); j++)
      SigmaXi(i,j) = xsd(i) * xsd(j) * (Type(1) / sdVecX(i)) * (Type(1) / sdVecX(j)) * SigmaZetaGPS(i,j);
  
  // calculate subsidence covariance
  matrix<Type> SigmaY = G * SigmaSlipFault * G.transpose();
  for(int i=0; i<ny; i++) {
    for(int j=0; j<ny; j++) {
      if(i == j)
        SigmaY(i,i) = SigmaY(i,i) + pow(yInflateSD(i), Type(2));
      else
        SigmaY(i,j) = SigmaY(i,j) * zeroMask(i,j);
    }
  }
  
  // calculate locking rate covariance
  matrix<Type> SigmaX(SigmaZetaGPS);
  for(int i=0; i<nx; i++) {
    for(int j=0; j<nx; j++) {
      SigmaX(i,j) = gammaVec(i) * gammaVec(j) * taperX(i) * taperX(j) * SigmaZetaGPS(i,j) + SigmaXi(i,j);
    }
  }
  
  // calculate subsidence negative log likelihood
  // NOTE: density::MVNORM_t<Type> function gives the negative log likelihood, not the likelihood itself
  Type nll = density::MVNORM_t<Type>(SigmaY)(y - muY);
  nll  +=  density::MVNORM_t<Type>(SigmaX)(x - muX);
  return nll;
}



















