#include <TMB.hpp>

// define double version of the taper function
double taper(double depth, double lambda, double dStar) {
  if(fabs(lambda, 2.) < 0.0000005) {
    // this is a limit as lambda approaches 0. Coding this explicitly avoids numerical instability
    return 1. - pow(depth/dStar, 2);
  } else if(depth > dStar) {
    return 0;
  } else if(depth <= 0) {
    return 1;
  } else {
    Type scaledDepth = pow(depth/dStar, 2) * pow(lambda, 2);
    Type ans = 1 - (1 - exp(-scaledDepth))/(1 - exp(-pow(lambda, 2)));
    return ans;
  }
}

// define scalar version of the taper function
template<class Type>
Type taper(Type depth, Type lambda, Type dStar) {
  CppAD::vector<Type> tdepth(1);
  CppAD::vector<Type> tlambda(1);
  CppAD::vector<Type> tdStar(1);
  tdepth(0) = depth;
  tlambda(0) = lambda;
  tdStar(0) = dStar;
  return taper(tdepth, tlambda, tdStar)[0];
}

// Vectorized version
VECTORIZE1_t(taper)

// reference on reverse mode differentiation:
// https://rufflewind.com/2016-12-30/reverse-mode-automatic-differentiation
TMB_ATOMIC_VECTOR_FUNCTION(
  // ATOMIC_NAME
  taper
  ,
  // OUTPUT_DIM
  1,
  // ATOMIC_DOUBLE
  ty[0] = taper(tdepths(0), tlambda(0), tdStar(0)); // Call the 'double' version
,
// ATOMIC_REVERSE
// example from https://kaskr.github.io/adcomp/_book/AtomicFunctions.html
// Type W  = ty[0];                    // Function value from forward pass
// Type DW = 1. / (exp(W) * (1. + W)); // Derivative
// px[0] = DW * py[0];                 // Reverse mode chain rule
Type pa = -py;
Type pb = py * 1. / pow(1. - exp(-pow(tlambda(0), 2.)));

Type scaledD = dRatio(0) * tlambda(0);
Type expD = exp(-pow(scaledD, Type(2))); // expD is the unnormalized taper function

Type px[0] = pa * (pow(scaledD, 2.) * 2 * ;
Type dRatio = tdepths(0)/tdStar(0);


expLam = exp(-lambda^2)
pVec = ((1 - expD)*expLam*2*lambda)/(1-expLam)^2    -    expD*2*scaledD*dRatio/(1 - expLam)
pVec[d > dStar] = 0
)

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_SCALAR(lambda);
  DATA_SCALAR(dStar);
  DATA_VECTOR(depths);
  PARAMETER_SCALAR(lambda);
  Type f = taper(depths, lambda, dStar).sum();
  return f;
}