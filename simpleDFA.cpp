#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_ARRAY(y);
  //DATA_ARRAY(logSdO);
  DATA_ARRAY(obs);
  
  int nSp = y.dim[0];
  int nT = y.dim[1];
  
  // For one-step-ahead residuals
  //DATA_ARRAY_INDICATOR(keep, y);

  // Parameters
  PARAMETER(logSdO); // Log of st. dev. for the observation error
  
  // Loadings matrix
  PARAMETER_MATRIX(Z);
  
  // Latent trends
  PARAMETER_MATRIX(x);
  
  // Observation error standard deviation.
  Type sdo = exp(logSdO);

  
  // Optimization target: negative log-likelihood (nll)
  Type nll = 0.0;

  // Latent random walk model. x(0) = 0. 
  for(int t = 1; t < x.cols(); ++t){
    for(int f = 0; f < x.rows(); ++f){
      nll -= dnorm(x(f, t), x(f, t-1), Type(1), true);
    }

        // Simulation block for process equation
    //SIMULATE {
    //  x(i) = rnorm(x(i-1), sdp);
    //}
  }
  
  // Observation model
  for(int i = 0; i < nSp; ++i){
    for(int t = 0; t < nT; ++t){
        nll -= dnorm(y(i, t), (Z.row(i) * x.col(t+1)).sum(), obs(i, t), true);
      //*----------------------- SECTION I --------------------------*/
      // Simulation block for observation equation
      //SIMULATE {
      //  y(i,t) = rnorm(x(t+1), sdo);
      //}
    }  
  }

  // State the transformed parameters to report
  // Using ADREPORT will return the point values and the standard errors
  // Note that we only need to specify this for parameters
  // we transformed, see section D above
  // The other parameters, including the random effects (states),
  // will be returned automatically
  ADREPORT(sdo);
  
  // Report simulated values
  //SIMULATE{
  //  REPORT(x);
  //  REPORT(y);
  //}

  return nll;
}
