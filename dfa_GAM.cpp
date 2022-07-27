// Dynamic Factor Analysis for multivariate time series
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  // Data
  DATA_ARRAY(y);
  // Observation standard error
  DATA_ARRAY(obs_se);
  
  int nSp = y.dim[0];
  int nT = y.dim[1];
  
  // For one-step-ahead residuals
  DATA_ARRAY_INDICATOR(keep, y);
  
  DATA_MATRIX(Z_pred);
  DATA_UPDATE(Z_pred);
  
  DATA_MATRIX(X_GAM);
  DATA_MATRIX(S_GAM);
  
  PARAMETER_MATRIX(beta_GAM);
  MVNORM_t<Type> nllmvnorm(S_GAM);
  // Parameters
  PARAMETER_VECTOR(log_re_sp); // log of sd for random effect by species
  
  // Loadings matrix
  PARAMETER_MATRIX(Z);
  
  int nfac = Z.cols();
  
  // Latent trends
  matrix<Type> x(nfac, nT);
  
  // Cluster center
  matrix<Type> x_pred(Z_pred.rows(), nT);
  
  // Matrix to hold predicted species trends
  matrix<Type> x_sp(nSp, nT);
  
  // Vector to hold penalties
  matrix<Type> pen(nfac, 1);
  
  // Random error by species
  vector<Type> re_sp;
  re_sp = exp(log_re_sp);
  // Optimization target: negative log-likelihood (nll)
  Type nll = 0.0;
  
  // Latent GAM model

  for(int f = 0; f < nfac; ++f){
    x.row(f) = X_GAM * beta_GAM.col(f); 
  } 

  // Species trends
  
  x_sp = Z * x;
  //for(int i = 0; i < nSp; ++i) {
  //  for(int t = 0; t < nT; ++t) {
  //    x_sp(i, t) = (Z.row(i) * (x.col(t))).sum();
  //  }
  //}  
  
  // Cluster center
  for (int t=0; t < nT; ++t) {
    x_pred.col(t) = Z_pred * (x.col(t));
  } 
  
  
  // Observation model
  for(int i = 0; i < nSp; ++i){
    // Skipping t = 0 when y(i, 0) is fixed at 0. Need to change this if y(i, 0) is not 0.
    // Also had had to change the index of x from t+1 to t, so that x is fixed at zero at time t=0.
    for(int t = 0; t < nT; ++t){ 
      nll -= keep(i) * dnorm(y(i, t), x_sp(i, t), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)), true); // with random effect
      //*----------------------- SECTION I --------------------------*/
      // Simulation block for observation equation
      SIMULATE {
        y(i,t) = rnorm(x_sp(i, t), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)));
        REPORT(y);
      }
    }  
  }
  
  
  
  //Penalty for very small random effects variances, this is to avoid issues with non semidefinite variance matrices.
  
  for (int i=0; i < nSp; ++i) {
    nll += exp(-3.5 * log_re_sp(i) - 20);
    // Smoothing penalty
      //nll += lambda * ((Z.row(i) * beta_GAM.transpose()) *  (S_GAM * (beta_GAM * Z.row(i).transpose()))).sum();
      
      //for (int f = 0; f < nfac; ++f) {
     //nll += lambda * Z(i,f) * Z(i, f) * pen(f);
  }
  
  for (int f=0; f < nfac; ++f) {
    nll += nllmvnorm(beta_GAM.col(f));
  }
  
  
  // State the transformed parameters to report
  // Using ADREPORT will return the point values and the standard errors
  // Note that we only need to specify this for parameters
  // we transformed, see section D above
  // The other parameters, including the random effects (states),Q
  // will be returned automatically
  ADREPORT(re_sp);
  ADREPORT(x_sp);
  ADREPORT(x);
  ADREPORT(x_pred);
  ADREPORT(Z_pred);
  
  // Report simulated values
  //SIMULATE{
  //  REPORT(x);
  //  REPORT(y);
  //}
  
  return nll;
}
