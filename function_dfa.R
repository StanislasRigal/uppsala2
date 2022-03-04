library(TMB)

dfa_model_se <- "
// Dynamic Factor Analysis for multivariate time series
#include <TMB.hpp>

template<class Type>
  Type objective_function<Type>::operator() ()
{
  // Data
  DATA_ARRAY(y);
  // Observation standard error
  DATA_ARRAY(obs_se);
  
  int nSp = y.dim[0];
  int nT = y.dim[1];
  
  // For one-step-ahead residuals
  //DATA_ARRAY_INDICATOR(keep, y);
  
  // Parameters
  PARAMETER_VECTOR(log_re_sp); // log of random error by species
  
  // Loadings matrix
  PARAMETER_MATRIX(Z);
  
  // Latent trends
  PARAMETER_MATRIX(x);
  
  //Test adding random effects explicitly
  //PARAMETER_MATRIX(eps);
  
  // Matrix to hold predicted species trends
  matrix<Type> x_sp(nSp, nT);
  
  // Random error by species
  vector<Type> re_sp;
  re_sp = exp(log_re_sp);



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
  
  for(int i = 0; i < nSp; ++i) {
    for(int t = 0; t < nT; ++t) {
      x_sp(i, t) = (Z.row(i) * x.col(t)).sum();
    }
  }
  
  
  // Observation model
  for(int i = 0; i < nSp; ++i){
  // Skipping t = 0 when y(i, 0) is fixed at 0. Need to change this if y(i, 0) is not 0.
  // Also had had to change the index of x from t+1 to t, so that x is fixed at zero at time t=0.
    for(int t = 1; t < nT; ++t){ 
      // nll -= dnorm(y(i, t), (Z.row(i) * x.col(t+1)).sum(), obs_se(i, t), true); // without random effect
       nll -= dnorm(y(i, t), (Z.row(i) * x.col(t)).sum(), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)), true); // with random effect
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
  ADREPORT(re_sp);
  ADREPORT(x_sp);
  
  // Report simulated values
  //SIMULATE{
    //  REPORT(x);
    //  REPORT(y);
    //}
  
  return nll;
  }"

  
  # Compile .cpp file with TMB DFA model
  
  if(!exists('dfamodel')){
    write(dfa_model_se, file = "dfa_model_se.cpp")
    dfamodel<-dfa_model_se
    compile("dfa_model_se.cpp")
    dyn.load(dynlib("dfa_model_se"))
  }else if(dfa_model_se != dfamodel){
    write(dfa_model_se, file = "dfa_model_se.cpp")
    dfamodel<-dfa_model_se
    compile("dfa_model_se.cpp")
    dyn.load(dynlib("dfa_model_se"))
  }

# Function to zscore error of time series values before DFA
  
Zscore_err<-function(err_ts, val_ts){
  return(err_ts/sd(val_ts,na.rm=T))
}
  
# Function to zscore time series values before DFA
  
Zscore<-function(val_ts){
  return((val_ts-mean(val_ts,na.rm=T))/sd(val_ts,na.rm=T))
}

# Function to compute AIC of DFA models
# Only works if obj has been already optimized

AIC.tmb = function(obj, tol = 0.01) {
  # Simple convergence check
  stopifnot(max(abs(obj$gr(obj$env$last.par.best[obj$env$lfixed()]))) < tol)
  as.numeric(2 * obj$env$value.best + 2*sum(obj$env$lfixed()))
}

# Generalise DFA

make_dfa <- function(data_ts, # dataset of time series
                     data_ts_se, # dataset of standard error of time series 
                     nfac,
                     AIC=TRUE) # number of trends for the DFA
{
  
  # Save input data for plot
  data_ts_save <- as.data.frame(data_ts)
  data_ts_se_save <- as.data.frame(data_ts_se)
  data_ts_save_long <- cbind(melt(data_ts_save[,-c(2)], id.vars=names(data_ts_save)[1]),
                             se=melt(data_ts_se_save[,-c(2)], id.vars=names(data_ts_se_save)[1])[,3])
  names(data_ts_save_long)[2:4] <- c("Year","value_orig","se_orig")
  data_ts_save_long$Year <- as.numeric(as.character(data_ts_save_long$Year))
  
  
  # Remove potential column of ts names
  
  data_ts <- data_ts %>% select_if(Negate(is.character))
  data_ts_se <- data_ts_se %>% select_if(Negate(is.character))
  
  # Remove first year if all se = 0
  
  # Jonas: First columns is also needed (the model doesn't 'know' that it is zero).
  # if(mean(as.data.frame(data_ts_se)[,1], na.rm=T)==0 & sd(as.data.frame(data_ts_se)[,1], na.rm=T)==0){
  #   data_ts_se <- data_ts_se[,-c(1)]
  #   data_ts <- data_ts[,-c(1)]
  # }
  
  # Zscore ts value and put se to the same scale
  
  # # Jonas: No need to zscore since we are estimating variances for each series.
  # for(i in 1:nrow(data_ts_se)){
  #   data_ts_se[i,] <- Zscore_err(data_ts_se[i,], data_ts[i,])
  # }
  # data_ts_se <- as.matrix(data_ts_se) 
  # data_ts <- as.matrix(t(apply(data_ts,1,Zscore)))
  
  # List of data for DFA
  
  # Jonas: Data should be on log scale. The best would be if SEs are also
  # on the log scale. Here I'm using a delta method approximation.
  dataTmb <- list(y = log(as.matrix(data_ts)), obs_se = as.matrix(data_ts_se/data_ts))
  
  # Prepare parameters for DFA
  
  nfac <- nfac # Number of factors
  ny <- nrow(data_ts) # Number of time series
  nT <- ncol(data_ts) # Number of time step
  
  # Worth trying multiple starting values for the optimisation to check that the right optimum is found.
  set.seed(1) 
  log_re_sp <- runif(ny, -1, 0)
  
  Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
  
  # Set constrained elements to zero
  
  constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
  Zinit[constrInd] <- 0
  
  #List of parameters for DFA
  
  tmbPar <-  list(log_re_sp=log_re_sp, Z = Zinit,
                  x=matrix(c(rep(0, nfac), rnorm(nfac * (nT - 1))),
                           ncol = nT, nrow = nfac))
  
  # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
  
  Zmap <- matrix(ncol = nfac, nrow = ny)
  Zmap[constrInd] <- NA
  Zmap[!constrInd] <- 1:sum(!constrInd)
  xmap = matrix(ncol = nT, nrow = nfac)
  xmap[,1] <- NA
  xmap[(nfac + 1) : length(tmbPar$x)] <- 1:(length(tmbPar$x) - nfac)
  tmbMap <- list(Z = as.factor(Zmap), x = as.factor(xmap))
  
  # Make DFA
  tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= c("x"), DLL= "dfa_model_se")
  tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
  
  # Check convergence
  conv <- grepl("relative convergence",tmbOpt$message)
  if(!conv){warning("Convergence issue")}
  
  sdRep <- summary(sdreport(tmbObj))
  
  
  # Temporary plot to check results.
  cbind(y = as.numeric(dataTmb$y), sdRep[grepl( 'x_sp', rownames(sdRep)),], year = rep(1996:2020, each = 16)) %>%
    as.data.frame %>% mutate(species = factor(unlist(rep(obs_se[, 1], 25)))) %>% ggplot(aes(year, y, col = species)) + geom_point() + geom_line(aes(y = Estimate, col = species)) + 
    geom_ribbon(aes(ymin = Estimate - 2* `Std. Error`, ymax = Estimate + 2* `Std. Error`, col = NULL), alpha = .4) + facet_wrap('species', scales = 'free_y')
  
  browser()
  
  #sdRep[grepl('Z|sdo', rownames(sdRep)),]
  
  # Get point estimates
  #x_hat <- matrix(sdRep[rownames(sdRep)=="x",1], nrow = nfac)
  
  x_hat <- (tmbObj$env$parList)()$x
  
  Z_hat <- (tmbObj$env$parList)(par=tmbOpt$par)$Z
  
  Z_hat_se <- sdRep[grepl('Z', rownames(sdRep)),2]
  for(i in 1:(nfac-1)){
    index_0 <- ny*i
    Z_hat_se <- append(Z_hat_se, rep(0,i), after=index_0)
  }
  Z_hat_se <- matrix(Z_hat_se, ncol=ncol(Z_hat), nrow=nrow(Z_hat))
  x_hat_se <- matrix(sdRep[grepl('x', rownames(sdRep)),2], nrow=nfac)
  
  
  
  ## Compute AIC
  
  if(AIC){
    aic <- AIC.tmb(tmbObj)
    print(aic)
  }
  
  # Jonas: Code below is now broken due to indexing changes.
  
  # Prepare data to plot
  
  if(!is.character(data_ts_save[,1]) & !is.factor(data_ts_save[,1])){
    data_ts_save <- data.frame(code_sp=as.character(row.names(data_ts_save)),
                               data_ts)
    data_ts_se_save <- data.frame(code_sp=as.character(row.names(data_ts_se_save)),
                                  data_ts_se)
  }else{
    data_ts_save <- data.frame(code_sp=data_ts_save[,1], data_ts)
    data_ts_se_save <- data.frame(code_sp=data_ts_se_save[,1], data_ts_se)
  }
  sp_ts <- data.frame(code_sp=data_ts_save[,1], Z_hat %*% x_hat)
  sp_se_ts <- data.frame(code_sp=data_ts_save[,1], abs(Z_hat_se %*% x_hat))
  for(i in 1:nrow(sp_se_ts)){
    sp_se_ts[i,-1] <- Zscore_err(sp_se_ts[i,-1], sp_ts[i,-1])
  }
  sp_ts <- data.frame(code_sp=data_ts_save[,1], t(apply(Z_hat %*% x_hat,1,Zscore)))
  
  data_to_plot_sp <- cbind(data_ts_save_long,
                           melt(data_ts_save, id.vars=names(data_ts_save)[1])[,3],
                           se=melt(data_ts_se_save, id.vars=names(data_ts_se_save)[1])[,3],
                           pred=melt(sp_ts, id.vars=names(data_ts_se_save)[1])[,3],
                           pred_se=melt(sp_se_ts, id.vars=names(data_ts_se_save)[1])[,3])
  #data_to_plot_sp$Year <- as.numeric(gsub("X", "", as.character(data_to_plot_sp$variable)))
  data_to_plot_tr <- data.frame(t(x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr_se <- data.frame(t(x_hat_se), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr_rot <- data.frame(t(solve(varimax(Z_hat)$rotmat) %*% x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr_rot_se <- data.frame(t(solve(varimax(Z_hat)$rotmat) %*% x_hat_se), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr <- cbind(melt(data_to_plot_tr, id.vars = "Year"),
                           se=melt(data_to_plot_tr_se, id.vars = "Year")[,3],
                           rot_tr=melt(data_to_plot_tr_rot, id.vars = "Year")[,3],
                           rot_tr_se=abs(melt(data_to_plot_tr_rot_se, id.vars = "Year")[,3]))
  
  data_loadings <- cbind(melt(data.frame(code_sp=data_ts_save[,1],
                                   Z_hat %*% varimax(Z_hat)$rotmat), id.vars="code_sp"),
                         se = melt(data.frame(code_sp=data_ts_save[,1],
                                              abs(Z_hat_se %*% varimax(Z_hat)$rotmat)), id.vars="code_sp")[,3])
   
  # Plots
  
  plot_sp <- ggplot(data_to_plot_sp, aes(x=Year, y=value)) + geom_point() +
    geom_pointrange(aes(ymax = value+se.value, ymin=value-se.value)) + 
    geom_line(aes(y=pred.value)) +
    geom_ribbon(aes(y=pred.value, ymax = pred.value+pred_se.value, ymin=pred.value-pred_se.value), alpha=0.5) +
    facet_wrap(code_sp ~ ., ncol=4) +
    theme_modern()
  
  plot_tr <- ggplot(data_to_plot_tr, aes(x=Year, y=rot_tr.value)) + 
    geom_line(aes(colour=variable))+
    geom_ribbon(aes(ymax = rot_tr.value+rot_tr_se.value, ymin=rot_tr.value-rot_tr_se.value,fill=variable), alpha=0.5) +
    theme_modern()
  
  plot_ld <- ggplot(data_loadings) + 
    geom_col(aes(value, code_sp, fill=variable)) +
    #geom_col(aes(value, reorder(code_sp, -value), fill=variable))+
    geom_errorbar(aes(x=value,y=code_sp,xmax = value+se.value, xmin=value-se.value), alpha=0.5) +
    facet_wrap(variable ~ ., ncol=4) +
    theme_modern() + theme(legend.position = "none")
  
  return(list(data_to_plot_sp, data_to_plot_tr, data_loadings,
              plot_sp, plot_tr, plot_ld))
}

