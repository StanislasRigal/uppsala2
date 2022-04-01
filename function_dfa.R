library(TMB)

dfa_model_se1 <- "
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
  DATA_ARRAY_INDICATOR(keep, y);
  
  // Parameters
  PARAMETER_VECTOR(log_re_sp); // log of sd for random effect by species
  
  // Loadings matrix
  PARAMETER_MATRIX(Z);
  PARAMETER_MATRIX(Z_pred);
  
  // Latent trends
  PARAMETER_MATRIX(x);
  
  // Cluster center
  matrix<Type> x_pred(Z_pred.rows(), nT);

  // Mean of latent trends
  matrix<Type> x_mean(x.rows(), 1);
  
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
    
    // Simulation block for process equation
    SIMULATE {
        x(f, t) = rnorm(x(f, t-1), Type(1));
        REPORT(x);
     }
    }
  }

  for (int f = 0; f < x.rows(); ++f) {
    x_mean(f) = x.row(f).sum() / nT;
    SIMULATE {
      x_mean(f) = x.row(f).sum() / nT;
    }
  }

  // Species trends
  for(int i = 0; i < nSp; ++i) {
    for(int t = 0; t < nT; ++t) {
      x_sp(i, t) = (Z.row(i) * (x.col(t) - x_mean)).sum();
    }
  }  
  
  // Cluster center
  for (int t=0; t < nT; ++t) {
    x_pred.col(t) = Z_pred * (x.col(t) - x_mean);
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
          y(i,t) = rnorm((Z.row(i) * (x.col(t) - x_mean)).sum(), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)));
          REPORT(y);
        }
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
  ADREPORT(x_pred);
  
  // Report simulated values
  //SIMULATE{
    //  REPORT(x);
    //  REPORT(y);
    //}
  
  return nll;
  }"

dfa_model_se2 <- "
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
  DATA_ARRAY_INDICATOR(keep, y);
  
  // Parameters
  PARAMETER_VECTOR(log_re_sp); // log of sd for random effect by species
  
  // Loadings matrix
  PARAMETER_MATRIX(Z);
  
  // Latent trends
  PARAMETER_MATRIX(x);
  
  // Mean of latent trends
  matrix<Type> x_mean(x.rows(), 1);
  
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
    
    // Simulation block for process equation
    SIMULATE {
        x(f, t) = rnorm(x(f, t-1), Type(1));
        REPORT(x);
     }
    }
  }
  for (int f = 0; f < x.rows(); ++f) {
    x_mean(f) = x.row(f).sum() / nT;
    SIMULATE {
      x_mean(f) = x.row(f).sum() / nT;
    }
  }
  // Species trends
  for(int i = 0; i < nSp; ++i) {
    for(int t = 0; t < nT; ++t) {
      x_sp(i, t) = (Z.row(i) * (x.col(t) - x_mean)).sum();
    }
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
          y(i,t) = rnorm((Z.row(i) * (x.col(t) - x_mean)).sum(), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)));
          REPORT(y);
        }
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

  
dfa_model_se <- dfa_model_se1
  
# Compile .cpp file with TMB DFA model
  
  if(!exists('dfamodel')){
    write(dfa_model_se, file = "dfa_model_se.cpp")
    dfamodel <- dfa_model_se
    compile("dfa_model_se.cpp")
    dyn.load(dynlib("dfa_model_se"))
  }else if(dfa_model_se != dfamodel){
    write(dfa_model_se, file = "dfa_model_se.cpp")
    dfamodel <- dfa_model_se
    compile("dfa_model_se.cpp")
    dyn.load(dynlib("dfa_model_se"))
  }
  

# Function to zscore error of time series values before DFA

Zscore_err<-function(err_ts, val_ts){
  return(err_ts/sd(val_ts,na.rm=T))
}

# Function to zscore time series values before DFA

Zscore <- function(val_ts){
  return((val_ts-mean(val_ts,na.rm=T))/sd(val_ts,na.rm=T))
}

# Function to compute AIC of DFA models (called inside make_dfa)
# Only works if obj has been already optimized
# AIC is computed excluding zero variance

AIC.tmb <- function(obj, tol = 0.01, dontCount) {
  # Simple convergence check
  stopifnot(max(abs(obj$gr(obj$env$last.par.best[obj$env$lfixed()]))) < tol)
  
  # AIC
  
  as.numeric(2 * obj$env$value.best + 2*(sum(obj$env$lfixed()) - dontCount))
}

# Prevent function to stop (to be use in loops)  
AIC.tmb2 <- function(obj, tol = 0.01, dontCount){
  tryCatch(AIC.tmb(obj, tol = 0.01, dontCount),
           error=function(e) NA)}


# Generalise DFA

make_dfa <- function(data_ts, # dataset of time series
                     data_ts_se, # dataset of standard error of time series 
                     nfac, # number of trends for the DFA
                     rand_seed=1, # initial values for the sd of the random effect
                     AIC=TRUE, # display AIC
                     with_kmeans=FALSE, # if TRUE, allow for first and second steps for clustering analysis
                     first_step=FALSE, # if TRUE, first DFA.cpp is use as the first step for clustering analysis
                                       # if FALSE, second step
                     ngroup=NULL, # number of group from kmeans, required for second step of cluster analysis
                     param_first_step=NULL, # parameters from first step
                     Z_pred_from_kmeans=NULL) # cluster loadings from kmeans, required for second step of cluster analysis
{
  
  # Save input data for plot
  
  data_ts_save <- as.data.frame(data_ts)
  data_ts_se_save <- as.data.frame(data_ts_se)
  data_ts_save_long <- cbind(melt(data_ts_save, id.vars=names(data_ts_save)[1]),
                             se=melt(data_ts_se_save, id.vars=names(data_ts_se_save)[1])[,3])
  names(data_ts_save_long)[2:4] <- c("Year","value_orig","se_orig")
  data_ts_save_long$Year <- as.numeric(as.character(data_ts_save_long$Year))
  
  # Remove potential column of ts names
  
  data_ts <- data_ts %>% select_if(Negate(is.character))
  data_ts_se <- data_ts_se %>% select_if(Negate(is.character))
  
  # List of data for DFA
  
  # SE is on log-scale, so no transformation needed
  dataTmb <- list(y = log(as.matrix(data_ts)),
                  obs_se = as.matrix(data_ts_se))
  
  if(with_kmeans==T){
    if(first_step==T){
      # Prepare parameters for DFA
      
      ngroup <- 10 # Arbitrary number of groups
      nfac <- nfac # Number of factors
      ny <- nrow(data_ts) # Number of time series
      nT <- ncol(data_ts) # Number of time step
      
      # Worth trying multiple starting values for the optimisation to check that the right optimum is found.
      set.seed(rand_seed) 
      log_re_sp <- runif(ny, -1, 0)
      
      Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
      
      # Set constrained elements to zero
      
      constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
      Zinit[constrInd] <- 0
      
      Z_predinit <- matrix(rep(0, ngroup * nfac), ncol = nfac)
      
      # List of parameters for DFA
      
      tmbPar <-  list(log_re_sp=log_re_sp, Z = Zinit, Z_pred = Z_predinit,
                      x=matrix(c(rep(0, nfac), rnorm(nfac * (nT - 1))),
                               ncol = nT, nrow = nfac))
      
      # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
      
      Zmap <- matrix(ncol = nfac, nrow = ny)
      Zmap[constrInd] <- NA
      Zmap[!constrInd] <- 1:sum(!constrInd)
      xmap <- matrix(ncol = nT, nrow = nfac)
      xmap[,1] <- NA
      xmap[(nfac + 1) : length(tmbPar$x)] <- 1:(length(tmbPar$x) - nfac)
      Z_predmap <- matrix(NA, ncol = nfac, nrow = ngroup)
      
      tmbMap <- list(Z = as.factor(Zmap), Z_pred = as.factor(Z_predmap),
                     x = as.factor(xmap))
      
      # Recompile .cpp
      
      dfa_model_se <- dfa_model_se1
      
      if(!exists('dfamodel')){
        write(dfa_model_se, file = "dfa_model_se.cpp")
        dfamodel <- dfa_model_se
        compile("dfa_model_se.cpp")
        dyn.load(dynlib("dfa_model_se"))
      }else if(dfa_model_se != dfamodel){
        write(dfa_model_se, file = "dfa_model_se.cpp")
        dfamodel <- dfa_model_se
        compile("dfa_model_se.cpp")
        dyn.load(dynlib("dfa_model_se"))
      }
      
      
      # Make DFA
      
      tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= c("x"), DLL= "dfa_model_se")
      tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
      
      # Avoid infinite SE when SD are close or equal to zero
      
      oh <- optimHess(tmbOpt$par, fn=tmbObj$fn, gr=tmbObj$gr)
      singularThreshold = 1e-5
      nSingular = sum(diag(oh) < singularThreshold)
      diag(oh) <- pmax(diag(oh), singularThreshold)
      
      sdRep <- summary(sdreport(tmbObj, hessian.fixed = oh))
      
      # Save parameter for the second step
      
      param_first_step <- (tmbObj$env$parList)()
      
            
    }else{
      # remove attr from param_first_step
      
      attr(param_first_step,"check.passed") <- NULL
      
      # Prepare parameters for DFA
      
      ngroup <- ngroup # number of group from kmeans
      nfac <- nfac # Number of factors
      ny <- nrow(data_ts) # Number of time series
      nT <- ncol(data_ts) # Number of time step
      
      Z_predinit <- Z_pred_from_kmeans
      
      # List of parameters for DFA
      
      tmbPar <-  list(log_re_sp = param_first_step$log_re_sp,
                      Z = param_first_step$Z, Z_pred = Z_predinit,
                      x = param_first_step$x)
      
      # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
      constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
      Zmap <- matrix(ncol = nfac, nrow = ny)
      Zmap[constrInd] <- NA
      Zmap[!constrInd] <- 1:sum(!constrInd)
      xmap <- matrix(ncol = nT, nrow = nfac)
      xmap[,1] <- NA
      xmap[(nfac + 1) : length(tmbPar$x)] <- 1:(length(tmbPar$x) - nfac)
      
      tmbMap <- list(Z = as.factor(Zmap), # no need to fix Z_pred as it now estimated from kmeans
                     x = as.factor(xmap))
      
      # Recompile .cpp
      
      dfa_model_se <- dfa_model_se1
      
      if(!exists('dfamodel')){
        write(dfa_model_se, file = "dfa_model_se.cpp")
        dfamodel <- dfa_model_se
        compile("dfa_model_se.cpp")
        dyn.load(dynlib("dfa_model_se"))
      }else if(dfa_model_se != dfamodel){
        write(dfa_model_se, file = "dfa_model_se.cpp")
        dfamodel <- dfa_model_se
        compile("dfa_model_se.cpp")
        dyn.load(dynlib("dfa_model_se"))
      }
      
      # Make DFA
      
      tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= c("x"), DLL= "dfa_model_se")
      tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))

      # Avoid infinite SE when SD are close or equal to zero
      
      oh <- optimHess(tmbOpt$par, fn=tmbObj$fn, gr=tmbObj$gr)
      singularThreshold = 1e-5
      nSingular = sum(diag(oh) < singularThreshold)
      diag(oh) <- pmax(diag(oh), singularThreshold)
      
      # Prepare par.fixed
      tmbPar2 <- unlist(tmbPar)
      tmbPar2 <- tmbPar2[tmbPar2!=0]
      r <- tmbObj$env$random
      tmbPar3 <- tmbPar2[-r]
      
      sdRep <- summary(sdreport(tmbObj, hessian.fixed = oh, par.fixed = tmbPar3, ignore.parm.uncertainty = T))
      #sdRep <- summary(sdreport(tmbObj, par.fixed = tmbPar3 ))
      
    }
  }else{
    
    # Create NA parameter object for output when with_kmeans = F
    
    param_first_step <- NA
    
    # Prepare parameters for DFA
    
    nfac <- nfac # Number of factors
    ny <- nrow(data_ts) # Number of time series
    nT <- ncol(data_ts) # Number of time step
    
    # Worth trying multiple starting values for the optimisation to check that the right optimum is found.
    set.seed(rand_seed) 
    log_re_sp <- runif(ny, -1, 0)
    
    Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
    
    # Set constrained elements to zero
    
    constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
    Zinit[constrInd] <- 0
    
    # List of parameters for DFA
    
    tmbPar <-  list(log_re_sp=log_re_sp, Z = Zinit,
                    x=matrix(c(rep(0, nfac), rnorm(nfac * (nT - 1))),
                             ncol = nT, nrow = nfac))
    
    # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
    
    Zmap <- matrix(ncol = nfac, nrow = ny)
    Zmap[constrInd] <- NA
    Zmap[!constrInd] <- 1:sum(!constrInd)
    xmap <- matrix(ncol = nT, nrow = nfac)
    xmap[,1] <- NA
    xmap[(nfac + 1) : length(tmbPar$x)] <- 1:(length(tmbPar$x) - nfac)
    tmbMap <- list(Z = as.factor(Zmap), x = as.factor(xmap))
    
    dfa_model_se <- dfa_model_se2
    
    if(!exists('dfamodel')){
      write(dfa_model_se, file = "dfa_model_se.cpp")
      dfamodel <- dfa_model_se
      compile("dfa_model_se.cpp")
      dyn.load(dynlib("dfa_model_se"))
    }else if(dfa_model_se != dfamodel){
      write(dfa_model_se, file = "dfa_model_se.cpp")
      dfamodel <- dfa_model_se
      compile("dfa_model_se.cpp")
      dyn.load(dynlib("dfa_model_se"))
    }
    
    # Make DFA
    
    tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= c("x"), DLL= "dfa_model_se")
    tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
  
    # Avoid infinite SE when SD are close or equal to zero
    
    oh <- optimHess(tmbOpt$par, fn=tmbObj$fn, gr=tmbObj$gr)
    singularThreshold = 1e-5
    nSingular = sum(diag(oh) < singularThreshold)
    diag(oh) <- pmax(diag(oh), singularThreshold)
    sdRep <- summary(sdreport(tmbObj, hessian.fixed = oh))
  }
  
  # Check convergence
  
  conv <- grepl("relative convergence",tmbOpt$message)
  if(!conv){warning(paste0("Convergence issue:", tmbOpt$message))}
  
  # Get point estimates
  #x_hat <- matrix(sdRep[rownames(sdRep)=="x",1], nrow = nfac)
  
  x_hat <- (tmbObj$env$parList)()$x
  
  Z_hat <- (tmbObj$env$parList)(par=tmbOpt$par)$Z
  
  Z_hat_se <- sdRep[rownames(sdRep)=="Z",2]
  for(i in 1:(nfac-1)){
    index_0 <- ny*i
    Z_hat_se <- append(Z_hat_se, rep(0,i), after=index_0)
  }
  Z_hat_se <- matrix(Z_hat_se, ncol=ncol(Z_hat), nrow=nrow(Z_hat))
  x_hat_se <- matrix(c(rep(0,nfac),sdRep[rownames(sdRep)=="x",2]), nrow=nfac)
  
  # Compute AIC
  if(AIC){
    aic <- AIC.tmb2(tmbObj, dontCount = 0) 
    aic2 <- AIC.tmb2(tmbObj, dontCount = nSingular) # Not sure if this is ok, should be checked.
    writeLines(paste('AIC: ', aic))
    writeLines(paste('AIC not counting singular random effects: ', aic2))
  } else {aic <- NA}
  
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
  
  # Back-transform log values of prediction of species ts
  
  sp_ts <- data.frame(code_sp=data_ts_save[,1],
                      matrix(sdRep[rownames(sdRep)=="x_sp",1], nrow=ny))
  
  # Jonas: This does not work as it will not give the SE of exp(x_sp).
  # It is better to transform in the plot instead. I.e. plot(exp(x_sp)) and
  # use ymin = exp(x_sp - 2*x_sp_se), ymax = exp(x_sp + 2*x_sp_se) to get 
  # confidence bands. I have changed this in plot_sp below.
  sp_se_ts <- data.frame(code_sp=data_ts_save[,1],
                         matrix(sdRep[rownames(sdRep)=="x_sp",2], nrow=ny))
  
  # Data for species time-series plot
  
  data_to_plot_sp <- cbind(data_ts_save_long,
                           melt(data_ts_save, id.vars=names(data_ts_save)[1])[,3],
                           se=melt(data_ts_se_save, id.vars=names(data_ts_se_save)[1])[,3],
                           pred=melt(sp_ts, id.vars=names(data_ts_se_save)[1])[,3],
                           pred_se=melt(sp_se_ts, id.vars=names(data_ts_se_save)[1])[,3])
  
  # Data for DFA trend plot
  
  data_to_plot_tr <- data.frame(t(x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr_se <- data.frame(t(x_hat_se), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  
  # Add rotated trends
  
  data_to_plot_tr_rot <- data.frame(t(solve(varimax(Z_hat)$rotmat) %*% x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  
  #Jonas: The below does not work, i.e. you will not get the correct SE from doing the same transformation as for the prediction.
  # Getting the SEs for the varimax rotation will be somewhat complicated, better leave it for now, I don't really see that we need them.
  #data_to_plot_tr_rot_se <- data.frame(t(solve(varimax(Z_hat)$rotmat) %*% x_hat_se), Year=min(data_to_p´lot_sp$Year):max(data_to_plot_sp$Year))
  
  data_to_plot_tr <- cbind(melt(data_to_plot_tr, id.vars = "Year"),
                           se=melt(data_to_plot_tr_se, id.vars = "Year")[,3], # This SE is ok as it comes directly from TMB
                           rot_tr=melt(data_to_plot_tr_rot, id.vars = "Year")[,3])
  
  # Data for species loadings
  
  data_loadings <- cbind(melt(data.frame(code_sp=data_ts_save[,1],
                                         Z_hat %*% varimax(Z_hat)$rotmat), id.vars="code_sp"),
                         se.value = NA)
  
  # Plots
  
  plot_sp <- ggplot(data_to_plot_sp, aes(x=Year, y=value)) + geom_point() +
    geom_pointrange(aes(ymax = value*exp(1.96 * se.value), ymin=value * exp(-1.96 * se.value))) + 
    geom_line(aes(y=exp(pred.value))) +
    geom_ribbon(aes(y=exp(pred.value), ymax = exp(pred.value+1.96*pred_se.value), ymin=exp(pred.value-1.96*pred_se.value)), alpha=0.5) +
    facet_wrap(code_sp ~ ., ncol=4, scales = "free") +
    theme_modern()
  
  plot_tr <- ggplot(data_to_plot_tr, aes(x=Year, y=rot_tr.value)) + 
    geom_line(aes(colour=variable))+
    theme_modern()
  
  plot_ld <- ggplot(data_loadings) + 
    geom_col(aes(value, code_sp, fill=variable)) +
    geom_errorbar(aes(x=value,y=code_sp,xmax = value+se.value, xmin=value-se.value), alpha=0.5) +
    facet_wrap(variable ~ ., ncol=4) +
    theme_modern() + theme(legend.position = "none")
  
  return(list(data_to_plot_sp, data_to_plot_tr, data_loadings,
              plot_sp, plot_tr, plot_ld, aic, param_first_step, sdRep))
}

# Group species a posteriori

group_from_dfa <- function(dfa_res, species_sub, eco_reg=FALSE, weight=FALSE){
  
  # Get loadings from DFA
  dfa_res_val <- dcast(dfa_res[[3]], code_sp~variable, value.var = "value")
  dfa_res_se <- dcast(dfa_res[[3]], code_sp~variable, value.var = "se.value")
  names(dfa_res_se) <- c("code_sp",paste0("se_",names(dfa_res_val[,-1])))
  
  mat_loading <- as.matrix(dfa_res_val[,-1])
  nb_dim <- ncol(dfa_res_se) -1 
  
  # Calculate weight for each species as the inverse of the volume of the n dimension ellipse (one dimension by DFA trends) defined by SE of each loadings (equivallent to semi axes in the ellipse)
  if(weight==TRUE){
    weight_loading <- apply(dfa_res_se[,-1], 1, 
                            FUN= function(x){
                              vol <- 2/nb_dim * (pi^(nb_dim/2)) / gamma(nb_dim/2) * prod(x)
                              return(vol)
                            })
    weight_loading <- weight_loading/min(weight_loading)
  }else{
    weight_loading <- 1
  }
  
  
  # Calculate gap statistic to find the best number of clusters
  gap_stat <- clusGap(mat_loading,
                      FUN = kmeans,
                      nstart = 25,
                      K.max = 10,
                      B = 500)
  
  # Plot number of clusters vs. gap statistic and let the user choose the number of cluster
  print(fviz_gap_stat(gap_stat))
  print(fviz_nbclust(mat_loading, kmeans, method = "wss"))
  #print(fviz_nbclust(mat_loading, kmeans, method = "silhouette"))
  nb_group <- as.numeric(readline(prompt = "Enter number of clusters: "))
  
  # Compute kmeans
  df.kmeans <- cclust(mat_loading, nb_group, weights = 1/weight_loading, 
                      method = "hardcl")
  #plot(mat_loading, col=predict(df.kmeans))
  #points(df.kmeans@centers, pch="x", cex=2, col=3)
  
  myPCA <- prcomp(mat_loading, scale. = F, center = F)
  
  # Group all info as output
  if(eco_reg==FALSE){
    kmeans_res <- list(merge(data.frame(code_sp=dfa_res_val[,1],
                                        myPCA$x[,1:2],
                                        group=as.factor(predict(df.kmeans)),
                                        dfa_res_val[,-1],
                                        dfa_res_se[,-1]),species_sub[,c("name_long","code_sp")],by="code_sp"),
                       data.frame(group=as.factor(1:nb_group),df.kmeans@centers,
                                  df.kmeans@centers %*% myPCA$rotation[,1:2]))
  }else{
    kmeans_res <- list(merge(data.frame(code_sp=dfa_res_val[,1],
                                        myPCA$x[,1:2],
                                        group=as.factor(predict(df.kmeans)),
                                        dfa_res_val[,-1],
                                        dfa_res_se[,-1]),species_sub[,c("name_long_eco","code_sp_eco")],by.x="code_sp", by.y="code_sp_eco"),
                       data.frame(group=as.factor(1:nb_group),df.kmeans@centers,
                                  df.kmeans@centers %*% myPCA$rotation[,1:2]))
  }
  
  
  # Get area (convex hull) for each group for plotting
  find_hull <- function(x){x[chull(x$PC2, x$PC1), ]}
  hulls <- ddply(kmeans_res[[1]], "group", find_hull)
  
  # Get average trend for each group
  trend_dfa <- as.matrix(dcast(as.data.frame(dfa_res[[2]])[,c("Year","variable","value")],
                               Year~variable, value.var = "value"))
  trend_dfa_se <- as.matrix(dcast(as.data.frame(dfa_res[[2]])[,c("Year","variable","se.value")],
                                  Year~variable, value.var = "se.value"))
  mean_trend <- data.frame(trend_dfa[,1],
                           trend_dfa[,-1] %*% t(df.kmeans@centers),
                           sqrt(trend_dfa_se[,-1]^2 %*% t(df.kmeans@centers)^2))
  names(mean_trend) <- c("year",paste0("group_",1:nb_group),
                         paste0("se_group_",1:nb_group))
  
  # Prepare subgraph of these trend to add to the final graph
  centroids <- as.data.frame(hulls %>% group_by(group) %>% summarize(PC1=mean(PC1), PC2=mean(PC2))) 
  
  graph <- setNames(lapply(1:nb_group, function(i){
    test <- mean_trend[,c(1,i+1, i+nb_group+1)]
    test$Index_SE <- test[,3]
    test$Index <- test[,2]
    ggplot(test, aes(x=year, y=Index)) +
      geom_line() +
      geom_ribbon(aes(ymin=Index-Index_SE,ymax=Index+Index_SE),alpha=0.7, col="black",fill="white")+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)
  }), names(mean_trend)[2:(nb_group+1)])
  
  centroids_data <- tibble(x=centroids$PC1,
                           y=centroids$PC2,
                           width=0.03,
                           pie = graph)
  
  # Plot final output
  final_plot <- ggplot(kmeans_res[[1]], aes(PC1,PC2, col=group, fill=group)) +
    geom_point() + geom_polygon(data=hulls, alpha=.2) +
    geom_point(data=(kmeans_res[[2]]), shape=2) +
    geom_text(label=kmeans_res[[1]]$name_long, nudge_x = 0.005, nudge_y = 0.005, check_overlap = F) +
    geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroids_data) +
    theme_modern()
  
  return(list(kmeans_res,final_plot,mean_trend, centroids, hulls))
}

# Plot groups and clustering trends

plot_group <- function(nb_group, centroids, kmeans_res, dfa_second_step, hulls){
  
  data_trend_group <- data.frame(group=rep(paste0("g",1:nb_group),23),
                                 year=sort(rep(c(1998:2020), nb_group)),
                                 dfa_second_step[[9]][grepl("x_pred",row.names(dfa_second_step[[9]])),])
  
  
  graph <- setNames(lapply(1:nb_group, function(i){
    test <- data_trend_group[data_trend_group$group==paste0("g",i),]
    test$Index_SE <- test$Std..Error
    test$Index <- test$Estimate
    ggplot(test, aes(x=year, y=Index)) +
      geom_line() +
      geom_ribbon(aes(ymin=Index-Index_SE,ymax=Index+Index_SE),alpha=0.7, col="black",fill="white")+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)
  }), levels(as.factor(data_trend_group$group)))
  
  centroids_data <- tibble(x=centroids$PC1,
                           y=centroids$PC2,
                           width=0.03,
                           pie = graph)
  
  # Get polygons and area
  xys <- st_as_sf(hulls, coords=c("PC1","PC2"))
  
  polys <- xys %>% 
    dplyr::group_by(group) %>% 
    dplyr::summarise() %>%
    st_cast("POLYGON") %>%
    st_convex_hull()
  centroids$area <- st_area(polys)
  random_point <- ddply(centroids,.(group),.fun=function(x){
    x <- as.numeric(x)
    y <- data.frame(PC1=rnorm(round(10000*x[4]),x[2],2*x[4]),
                    PC2=rnorm(round(10000*x[4]),x[3],2*x[4]))
    return(y)
  })
  
  # Plot final output
  final_plot <- ggplot(kmeans_res[[1]], aes(PC1,PC2)) +
    #geom_polygon(data=hulls, alpha=.2) +
    stat_density_2d(data=random_point, aes(PC1,PC2,fill = ..level..), geom = "polygon") +
    #stat_density_2d(data=random_point, aes(PC1,PC2,fill = ..density..), geom = "raster", contour = FALSE) +
    #scale_x_continuous(expand = c(0, 0)) +
    #scale_y_continuous(expand = c(0, 0)) +
    scale_fill_gradient2(low="white", mid="yellow", high = "red",midpoint = 100) +
    geom_point() + #geom_point(data=(kmeans_res[[2]]), shape=2) +
    geom_text(label=kmeans_res[[1]]$name_long, nudge_x = 0.005, nudge_y = 0.005, check_overlap = F) +
    geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroids_data) +
    theme_modern() +
    theme(legend.position='none')
  
  return(final_plot)
}