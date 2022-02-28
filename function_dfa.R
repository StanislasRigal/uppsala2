library(TMB)

dfa_model_se <- "
// Dynamic Factor Analysis for multivariate time series
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
  
  # Remove potential column of ts names
  
  data_ts <- data_ts %>% select_if(Negate(is.character))
  data_ts_se <- data_ts_se %>% select_if(Negate(is.character))
  
  # Remove first year if all se = 0
  
  if(mean(as.data.frame(data_ts_se)[,1], na.rm=T)==0 & sd(as.data.frame(data_ts_se)[,1], na.rm=T)==0){
    data_ts_se <- data_ts_se[,-c(1)]
    data_ts <- data_ts[,-c(1)]
  }
  
  # Zscore ts value and put se to the same scale
  
  for(i in 1:nrow(data_ts_se)){
    data_ts_se[i,] <- Zscore_err(data_ts_se[i,], y[i,])
  }
  data_ts_se <- as.matrix(data_ts_se) 
  data_ts <- as.matrix(t(apply(data_ts,1,Zscore)))
  
  # List of data for DFA
  
  dataTmb <- list(y = data_ts, obs = data_ts_se)
  
  # Prepare parameters for DFA
  
  nfac <- nfac # Number of factors
  ny <- nrow(data_ts) # Number of time series
  nT <- ncol(data_ts) # Number of time step
  
  Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
  
  # Set constrained elements to zero
  
  constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
  Zinit[constrInd] <- 0
  
  #List of parameters for DFA
  
  tmbPar <-  list(logSdO = 0, Z = Zinit,
                  x=matrix(c(rep(0, nfac), rnorm(nfac * nT)),
                           ncol = nT+1, nrow = nfac))
  
  # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
  
  Zmap <- matrix(ncol = nfac, nrow = ny)
  Zmap[constrInd] <- NA
  Zmap[!constrInd] <- 1:sum(!constrInd)
  xmap = matrix(ncol <- nT + 1, nrow = nfac)
  xmap[,1] <- NA
  xmap[(nfac + 1) : length(tmbPar$x)] <- 1:(length(tmbPar$x) - nfac)
  tmbMap <- list(Z = as.factor(Zmap), x = as.factor(xmap))
  
  # Make DFA
  
  tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= "x", DLL= "simpleDFA")
  tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
  
  # Check convergence
  conv <- grepl("relative convergence",tmbOpt$message)
  if(!conv){warning("Convergence issue")}
  sdRep <- summary(sdreport(tmbObj))
  #sdRep[grepl('Z|sdo', rownames(sdRep)),]
  
  # Get point estimates
  #x_hat <- matrix(sdRep[rownames(sdRep)=="x",1], nrow = nfac)
  
  x_hat <- (tmbObj$env$parList)()$x
  if(mean(x_hat[,1], na.rm=T)==0 & sd(x_hat[,1], na.rm=T)==0){
    x_hat <- x_hat[,-c(1)]
  }
  Z_hat <- (tmbObj$env$parList)(par=tmbOpt$par)$Z
  
  ## Compute AIC
  
  if(AIC){
    aic <- AIC.tmb(tmbObj)
    print(aic)
  }

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
  
  data_to_plot_sp <- cbind(melt(data_ts_save, id.vars=names(data_ts_save)[1]),
                           se=melt(data_ts_se_save, id.vars=names(data_ts_se_save)[1])[,3],
                           pred=melt(sp_ts, id.vars=names(data_ts_se_save)[1])[,3])
  data_to_plot_sp$Year <- as.numeric(gsub("X", "", as.character(data_to_plot_sp$variable)))
  
  data_to_plot_tr <- data.frame(t(x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr_rot <- data.frame(t(solve(varimax(Z_hat)$rotmat) %*% x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr <- cbind(melt(data_to_plot_tr, id.vars = "Year"),
                           rot_tr=melt(data_to_plot_tr_rot, id.vars = "Year")[,3])
  
  data_loadings <- melt(data.frame(code_sp=data_ts_save[,1],
                                   Z_hat %*% varimax(Z_hat)$rotmat), id.vars="code_sp")
  
  # Plots
  
  plot_sp <- ggplot(data_to_plot_sp, aes(x=Year, y=value)) + geom_point() +
    geom_pointrange(aes(ymax = value+se.value, ymin=value-se.value)) + 
    geom_line(aes(y=pred.value)) + facet_wrap(code_sp ~ ., ncol=4) +
    theme_modern()
  
  plot_tr <- ggplot(data_to_plot_tr, aes(x=Year, y=rot_tr.value)) + 
    geom_line(aes(colour=variable))+
    theme_modern()
  
  plot_ld <- ggplot(data_loadings) + 
    geom_col(aes(value, code_sp, fill=variable)) +
    #geom_col(aes(value, reorder(code_sp, -value), fill=variable))+
    facet_wrap(variable ~ ., ncol=4) +
    theme_modern() + theme(legend.position = "none")
  
  return(list(data_to_plot_sp, data_to_plot_tr, data_loadings,
              plot_sp, plot_tr, plot_ld))
}
