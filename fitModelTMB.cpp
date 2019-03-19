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
// define a smoothness penalty: numerical integral on the spline values
Type penaltyFun(vector<Type> splineValues, Type delta, Type logLambda) {
  vector<Type> integrand(splineValues.size() - 1);
  for(int i=0; i<integrand.size(); i++)
    integrand(i) = pow((splineValues(i + 1) - splineValues(i)) * (1 / delta), Type(2));
  return sum(integrand) * delta * exp(logLambda);
}

template<class Type>
// define a difference penalty: numerical integral on the splines' differences
Type penaltyDiffFun(vector<Type> splineDifferences, Type delta, Type logLambda) {
  vector<Type> integrand(splineDifferences.size());
  for(int i=0; i<integrand.size(); i++)
    integrand(i) = pow(splineDifferences(i), Type(2));
  return sum(integrand) * delta * exp(logLambda);
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
  DATA_MATRIX(DSStrikeCross);
  DATA_MATRIX(DSDipCross);
  DATA_MATRIX(sdBasisX);
  DATA_MATRIX(sdBasisY);
  DATA_MATRIX(meanBasisY);
  DATA_MATRIX(meanBasisX);
  DATA_MATRIX(meanBasisXGPS);
  DATA_MATRIX(lambdaBasisX);
  DATA_MATRIX(lambdaBasisXGPS);
  DATA_MATRIX(lambdaBasisY);
  DATA_MATRIX(gammaBasis);
  
  DATA_MATRIX(sdBasisPenalty);
  DATA_MATRIX(taperBasisPenalty);
  DATA_MATRIX(taperBasisGPSPenalty);
  DATA_MATRIX(gammaBasisPenalty);
  DATA_MATRIX(meanBasisPenalty);
  DATA_MATRIX(meanBasisGPSPenalty);
  
  DATA_VECTOR(faultDepths);
  DATA_VECTOR(xDepths);
  
  DATA_SCALAR(dStar);
  DATA_SCALAR(dStarGPS);
  DATA_SCALAR(deltaPenalty);
  DATA_SCALAR(penaltyMean);
  DATA_SCALAR(penaltySD);
  DATA_SCALAR(sharedPenalty);
  DATA_SCALAR(doTaperDiffPenalty);
  DATA_SCALAR(diffPenaltyMean);
  DATA_SCALAR(diffPenaltySD);
  DATA_SCALAR(diffMean);
  DATA_SCALAR(useHyperpriors);
  DATA_SCALAR(sharedSpatialProcess);
  DATA_SCALAR(jointShared);
  
  PARAMETER(logmu);
  PARAMETER_VECTOR(betaMean);
  PARAMETER(logMeanGPS);
  PARAMETER_VECTOR(betaMeanGPS);
  PARAMETER(betaTaperIntercept);
  PARAMETER_VECTOR(betaTaper);
  PARAMETER_VECTOR(betaTaperGPS);
  PARAMETER(betasdIntercept);
  PARAMETER_VECTOR(betasd);
  PARAMETER(betaGammaIntercept);
  PARAMETER_VECTOR(betaGamma);
  PARAMETER(logphi);
  PARAMETER(logalpha);
  PARAMETER(loglowInflate);
  PARAMETER(loghighInflate);
  PARAMETER(logitOmega);
  
  PARAMETER(betasdPenaltyLogLambda);
  PARAMETER(betaTaperPenaltyLogLambda);
  PARAMETER(betaTaperGPSPenaltyLogLambda);
  PARAMETER(betaGammaPenaltyLogLambda);
  PARAMETER(betaMeanPenaltyLogLambda);
  PARAMETER(betaMeanGPSPenaltyLogLambda);
  PARAMETER(taperDiffPenaltyLogLambda);
  PARAMETER(meanDiffPenaltyLogLambda);
  
  int nx = x.size();
  int ny = y.size();
  
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
  Type omega = 0;
  if(sharedSpatialProcess == Type(1)) {
    omega = exp(logitOmega)/(Type(1)+exp(logitOmega));
  }
  
  // construct lambda, standard deviation, and gamma vectors for both fault and locking rate data
  vector<Type> lambdaVecX = lambdaBasisX * betaTaper + lambdaBasisXGPS * betaTaperGPS;
  vector<Type> lambdaVecY = lambdaBasisY * betaTaper;
  vector<Type> sdVecX = sdBasisX * betasd;
  vector<Type> sdVecY = sdBasisY * betasd;
  vector<Type> meanVecX = meanBasisX * betaMean;
  vector<Type> meanVecY = meanBasisY * betaMean;
  vector<Type> gammaVec = gammaBasis * betaGamma;
  
  // add in intercepts
  for(int i=0; i<lambdaVecX.size(); i++)
    lambdaVecX(i) = lambdaVecX(i) + betaTaperIntercept;
  for(int i=0; i<lambdaVecY.size(); i++)
    lambdaVecY(i) = lambdaVecY(i) + betaTaperIntercept;
  for(int i=0; i<sdVecX.size(); i++)
    sdVecX(i) = sdVecX(i) + betasdIntercept;
  for(int i=0; i<sdVecY.size(); i++)
    sdVecY(i) = sdVecY(i) + betasdIntercept;
  for(int i=0; i<meanVecX.size(); i++)
    meanVecX(i) = meanVecX(i) + logmu;
  for(int i=0; i<meanVecY.size(); i++)
    meanVecY(i) = meanVecY(i) + logmu;
  for(int i=0; i<gammaVec.size(); i++)
    gammaVec(i) = gammaVec(i) + betaGammaIntercept;
  
  // add in a different mean for the gps data if necessary
  if(diffMean == Type(1)) {
    meanVecX = meanVecX + meanBasisXGPS * betaMeanGPS;
    for(int i=0; i<meanVecX.size(); i++)
      meanVecX(i) = meanVecX(i) + logMeanGPS;
  }
  
  // reparameterize standard deviation, mean, and gamma (lambda values are reparameterized later)
  for(int i=0; i<sdVecX.size(); i++)
    sdVecX(i) = exp(sdVecX(i));
  for(int i=0; i<sdVecY.size(); i++)
    sdVecY(i) = exp(sdVecY(i));
  for(int i=0; i<meanVecX.size(); i++)
    meanVecX(i) = exp(meanVecX(i));
  for(int i=0; i<meanVecY.size(); i++)
    meanVecY(i) = exp(meanVecY(i));
  for(int i=0; i<gammaVec.size(); i++)
    gammaVec(i) = exp(gammaVec(i));
  
  // if the penalty parameters are shared, set them all to the sd penalty
  if(sharedPenalty == Type(1)) {
    betaTaperPenaltyLogLambda = betasdPenaltyLogLambda;
    betaTaperGPSPenaltyLogLambda = betasdPenaltyLogLambda;
    betaGammaPenaltyLogLambda = betasdPenaltyLogLambda;
    betaMeanPenaltyLogLambda = betasdPenaltyLogLambda;
    betaMeanGPSPenaltyLogLambda = betasdPenaltyLogLambda;
  }
  
  // Now do the same for the penalty splines
  vector<Type> lambdaVecPenalty = taperBasisPenalty * betaTaper + taperBasisGPSPenalty * betaTaperGPS;
  vector<Type> lambdaVecGPSPenalty = taperBasisPenalty * betaTaper;
  vector<Type> sdVecPenalty = sdBasisPenalty * betasd;
  vector<Type> gammaVecPenalty = gammaBasisPenalty * betaGamma;
  vector<Type> meanVecPenalty = meanBasisPenalty * betaMean;
  vector<Type> meanVecGPSPenalty(meanVecPenalty);
  if(diffMean == Type(1))
    meanVecGPSPenalty = meanBasisPenalty * betaMean;
  
  // add in intercepts (since we are exponentiating the parameters, adding a constant 
  // term will affect the penalty)
  for(int i=0; i<lambdaVecPenalty.size(); i++)
    lambdaVecPenalty(i) = lambdaVecPenalty(i) + betaTaperIntercept;
  for(int i=0; i<lambdaVecGPSPenalty.size(); i++)
    lambdaVecGPSPenalty(i) = lambdaVecGPSPenalty(i) + betaTaperIntercept;
  for(int i=0; i<sdVecPenalty.size(); i++)
    sdVecPenalty(i) = sdVecPenalty(i) + betasdIntercept;
  for(int i=0; i<gammaVecPenalty.size(); i++)
    gammaVecPenalty(i) = gammaVecPenalty(i) + betaGammaIntercept;
  for(int i=0; i<meanVecPenalty.size(); i++)
    meanVecPenalty(i) = meanVecPenalty(i) + logmu;
  if(diffMean == Type(1)) {
    meanVecGPSPenalty = meanVecGPSPenalty + meanBasisGPSPenalty * betaMeanGPS;
    for(int i=0; i<meanVecGPSPenalty.size(); i++)
      meanVecGPSPenalty(i) = meanVecGPSPenalty(i) + logmu + logMeanGPS;
  }
  
  // reparameterize splines (here we penalize lambda directly, so we reparameterize it here)
  for(int i=0; i<lambdaVecGPSPenalty.size(); i++)
    lambdaVecGPSPenalty(i) = exp(lambdaVecGPSPenalty(i));
  for(int i=0; i<lambdaVecPenalty.size(); i++)
    lambdaVecPenalty(i) = exp(lambdaVecPenalty(i));
  for(int i=0; i<sdVecPenalty.size(); i++)
    sdVecPenalty(i) = exp(sdVecPenalty(i));
  for(int i=0; i<gammaVecPenalty.size(); i++)
    gammaVecPenalty(i) = exp(gammaVecPenalty(i));
  for(int i=0; i<meanVecPenalty.size(); i++)
    meanVecPenalty(i) = exp(meanVecPenalty(i));
  if(diffMean == Type(1)) {
    for(int i=0; i<meanVecGPSPenalty.size(); i++)
      meanVecGPSPenalty(i) = exp(meanVecGPSPenalty(i));
  }
  
  // now that we are done reparameterizing, print out the parameters
  ADREPORT(mu);
  ADREPORT(logmu);
  for(int i=0; i<betaMean.size(); i++)
    ADREPORT(betaMean(i));
  ADREPORT(logMeanGPS);
  for(int i=0; i<betaMeanGPS.size(); i++)
    ADREPORT(betaMeanGPS(i));
  ADREPORT(betasdIntercept);
  for(int i=0; i<betasd.size(); i++)
    ADREPORT(betasd(i));
  ADREPORT(betaTaperIntercept);
  for(int i=0; i<betaTaper.size(); i++)
    ADREPORT(betaTaper(i));
  for(int i=0; i<betaTaperGPS.size(); i++)
    ADREPORT(betaTaperGPS(i));
  ADREPORT(betaGammaIntercept);
  for(int i=0; i<betaGamma.size(); i++)
    ADREPORT(betaGamma(i));
  ADREPORT(betasdPenaltyLogLambda);
  ADREPORT(betaTaperPenaltyLogLambda);
  ADREPORT(betaTaperGPSPenaltyLogLambda);
  ADREPORT(betaGammaPenaltyLogLambda);
  ADREPORT(betaMeanPenaltyLogLambda);
  ADREPORT(betaMeanGPSPenaltyLogLambda);
  ADREPORT(taperDiffPenaltyLogLambda);
  ADREPORT(meanDiffPenaltyLogLambda);
  ADREPORT(lowInflate);
  ADREPORT(highInflate);
  ADREPORT(phi);
  ADREPORT(alpha);
  ADREPORT(omega);
  ADREPORT(logitOmega);
  
  // construct taper values
  vector<Type> taperX(nx);
  for(int i=0; i<lambdaVecX.size(); i++)
    taperX(i) = taper(xDepths(i), exp(lambdaVecX(i)), dStarGPS);
  
  vector<Type> taperFault(lambdaVecY.size());
  for(int i=0; i<lambdaVecY.size(); i++)
    taperFault(i) = taper(faultDepths(i), exp(lambdaVecY(i)), dStar);
  
  // calculate mean vectors for the subsidence and locking rate data
  vector<Type> meanSlip(taperFault.size());
  for(int i=0; i<meanSlip.size(); i++)
    meanSlip(i) = taperFault(i) * meanVecY(i);
  
  vector<Type> muX(lambdaVecX.size());
  for(int i=0; i<muX.size(); i++)
    muX(i) = gammaVec(i) * meanVecX(i) * taperX(i);
  
  vector<Type> muY = G * meanSlip;
  
  // concatenate mean vectors into a single joint vector if necessary
  vector<Type> muJoint(nx + ny);
  if(jointShared == Type(1)) {
    for(int i=0; i<muJoint.size(); i++) {
      if(i >= nx) {
        int iY = i - nx;
        muJoint(i) = muY(iY);
      }
      else {
        muJoint(i) = muX(i);
      }
    }
  }
  
  // do the same for the observations if necessary
  vector<Type> obs(nx + ny);
  if(jointShared == Type(1)) {
    for(int i=0; i<obs.size(); i++) {
      if(i >= nx) {
        int iY = i - nx;
        obs(i) = y(iY);
      }
      else {
        obs(i) = x(i);
      }
    }
  }
  
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
  
  matrix<Type> SigmaZetaCross(DSStrikeCross);
  if(jointShared == Type(1)) {
    // compute the cross covariances for zeta if necessary (from GPS to CSZ)
    for(int i=0; i<SigmaZetaCross.rows(); i++)
      for(int j=0; j<SigmaZetaCross.cols(); j++)
        SigmaZetaCross(i,j) = sdVecX(i) * sdVecY(j) * matern(sqrt(alphaSq * DSStrikeCross(i,j) + alphaInvSq * DSDipCross(i,j)), phi, kappa);
  }
  
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
  
  // calculate covariance matrix of locking rate data (assume same correlation as zeta)
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
      else {
        if(sharedSpatialProcess != Type(1)) {
          SigmaY(i,j) = SigmaY(i,j) * zeroMask(i,j);
        }
        else {
          // in this case, different earthquakes still share omega of the variance
          SigmaY(i,j) = SigmaY(i,j) * (omega + (1 - omega) * zeroMask(i,j));
        }
      }
    }
  }
  
  // calculate locking rate covariance
  matrix<Type> SigmaX(SigmaZetaGPS);
  for(int i=0; i<nx; i++) {
    for(int j=0; j<nx; j++) {
      SigmaX(i,j) = gammaVec(i) * gammaVec(j) * taperX(i) * taperX(j) * SigmaZetaGPS(i,j) + SigmaXi(i,j);
    }
  }
  
  // calculate the joint covariance between subsidence and gps data if necessary
  matrix<Type> SigmaCross(nx, ny);
  matrix<Type> SigmaJoint(nx + ny, nx + ny);
  if(jointShared == Type(1)) {
    // cross covariance from gps to csz
    for(int i=0; i<nx; i++) {
      for(int j=0; j<ny; j++) {
        SigmaCross(i,j) = gammaVec(i) * taperX(i) * omega * SigmaZetaCross(i,j) * taperFault(j);
      }
    }
    SigmaCross = SigmaCross * G.transpose();
    
    // joint covariance (gps, csz)
    for(int i=0; i<nx + ny; i++) {
      for(int j=0; j<nx + ny; j++) {
        if(i<nx && j<nx) {
          // this is the gps portion of the covariance matrix
          SigmaJoint(i,j) = SigmaX(i,j);
        }
        else if(i<nx && j >= nx) {
          // cross covariance from gps to csz
          int jY = j - nx;
          SigmaJoint(i,j) = SigmaCross(i, jY);
        }
        else if(i >= nx && j<nx) {
          // cross covariance from csz to gps
          int iY = i - nx;
          SigmaJoint(i,j) = SigmaCross(j, iY);
        }
        else {
          // this is the csz portion of the covariance matrix
          int iY = i - nx;
          int jY = j - nx;
          SigmaJoint(i,j) = SigmaY(iY,jY);
        }
      }
    }
  }
  
  // calculate subsidence negative log likelihood
  // NOTE: density::MVNORM_t<Type> function gives the negative log likelihood, not the likelihood itself
  Type nll = 0;
  if(jointShared != Type(1)) {
    nll += density::MVNORM_t<Type>(SigmaY)(y - muY);
    nll +=  density::MVNORM_t<Type>(SigmaX)(x - muX);
  }
  else {
    nll += density::MVNORM_t<Type>(SigmaJoint)(obs - muJoint);
  }
  
  // add-on the penalty functions/priors
  nll  +=  penaltyFun(lambdaVecPenalty, deltaPenalty, betaTaperPenaltyLogLambda);
  nll  +=  penaltyFun(lambdaVecGPSPenalty, deltaPenalty, betaTaperGPSPenaltyLogLambda);
  nll  +=  penaltyFun(sdVecPenalty, deltaPenalty, betasdPenaltyLogLambda);
  nll  +=  penaltyFun(gammaVecPenalty, deltaPenalty, betaGammaPenaltyLogLambda);
  nll  +=  penaltyFun(meanVecPenalty, deltaPenalty, betaMeanPenaltyLogLambda);
  if(diffMean == Type(1)) {
    for(int i=0; i<meanVecGPSPenalty.size(); i++)
      nll  +=  penaltyFun(meanVecGPSPenalty, deltaPenalty, betaMeanGPSPenaltyLogLambda);
  }
  if(doTaperDiffPenalty == Type(1)) {
    vector<Type> lambdaDiffPenalty(lambdaVecPenalty.size());
    for(int i=0; i<lambdaDiffPenalty.size(); i++)
      lambdaDiffPenalty(i) = lambdaVecPenalty(i) - lambdaVecGPSPenalty(i);
    nll  +=  penaltyDiffFun(lambdaDiffPenalty, deltaPenalty, taperDiffPenaltyLogLambda);
    
    if(diffMean == Type(1)) {
      vector<Type> meanDiffPenalty(meanVecPenalty.size());
      for(int i=0; i<meanDiffPenalty.size(); i++)
        meanDiffPenalty(i) = meanVecPenalty(i) - meanVecGPSPenalty(i);
      nll  +=  penaltyDiffFun(meanDiffPenalty, deltaPenalty, meanDiffPenaltyLogLambda);
    }
  }
  
  // add on the negative log hyperpriors (including jacobian factors)
  if(useHyperpriors == Type(1)) {
    nll  +=  -dnorm(betaTaperPenaltyLogLambda, penaltyMean, penaltySD,true) - betaTaperPenaltyLogLambda;
    nll  +=  -dnorm(betaTaperGPSPenaltyLogLambda, penaltyMean, penaltySD,true) - betaTaperGPSPenaltyLogLambda;
    nll  +=  -dnorm(betasdPenaltyLogLambda, penaltyMean, penaltySD,true) - betasdPenaltyLogLambda;
    nll  +=  -dnorm(betaGammaPenaltyLogLambda, penaltyMean, penaltySD,true) - betaGammaPenaltyLogLambda;
    nll  +=  -dnorm(betaMeanPenaltyLogLambda, penaltyMean, penaltySD,true) - betaMeanPenaltyLogLambda;
    if(doTaperDiffPenalty == Type(1)) {
      nll  +=  -dnorm(taperDiffPenaltyLogLambda, diffPenaltyMean, diffPenaltySD,true) - taperDiffPenaltyLogLambda;
      
      if(diffMean == Type(1))
        nll  +=  -dnorm(meanDiffPenaltyLogLambda, diffPenaltyMean, diffPenaltySD,true) - taperDiffPenaltyLogLambda;
    }
  }
  
  return nll;
}









