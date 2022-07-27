# I) TMB function

library(TMB)
library(mvtnorm)


## Compile .cpp file with TMB DFA model
compile("dfa_GAM.cpp")
dyn.load(dynlib("dfa_GAM"))

#compile("dfa_GAM_fix.cpp")
#dyn.load(dynlib("dfa_GAM_fix"))


AIC.tmb <- function(obj) {
  as.numeric(2 * obj$env$value.best + 2*(sum(obj$env$lfixed())))
}




gam_dfa <- function(data_ts, # Dataset of time series
                     data_ts_se, # Dataset of standard error of time series 
                     nfac, # Number of trends for the DFA
                     rand_seed=1, # Initial values for the sd of the random effect
                     AIC=TRUE, # Display AIC
                     silent = TRUE, 
                     lambda = 0,
                     control = list()
                     ) {
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

  # Overwrite any changes to default control 
  con <- list(nstart = 3, maxit = 10000, reltol = 1e-12, factr = 1e-11, gradtol = 1e-3, nlldeltatol = 1e-4, method = c('NLMINB', 'BFGS'))
  con[names(control)] <- control
  
  # Initialise Z_pred, actual values will be provided by group_from_dfa2
  
  Z_predinit <- matrix(rep(0, 10 * nfac), ncol = nfac)
  
  # SE is on log-scale, so no transformation needed
  g = mgcv::gam(y~ -1 + s(year, k = 20, bs = 'ts'), data = data.frame(y = rnorm(ncol(data_ts)),  year = as.numeric(colnames(data_ts))))
  X_GAM = mgcv::model.matrix.gam(g)
  
  # For dfa_GAM_fix
  #tmbData <- list(y = log(as.matrix(data_ts)),
  #                 obs_se = as.matrix(data_ts_se),
  #                 Z_pred = Z_predinit,
  #                 X_GAM = X_GAM[, 1:8],
  #                 S_GAM = g$smooth[[1]]$S[[1]][1:8, 1:8],
  #                 lambda = lambda)
  
  tmbData <- list(y = log(as.matrix(data_ts)),
                   obs_se = as.matrix(data_ts_se),
                   Z_pred = Z_predinit,
                   X_GAM = X_GAM,
                   S_GAM =  1/lambda * solve(g$smooth[[1]]$S[[1]]), # Results seems largely independent of value of lambda, only affects scalings of Z and beta.
                   lambda = lambda)
  
  # Prepare parameters for DFA
  
  nfac <- nfac # Number of factors
  ny <- nrow(data_ts) # Number of time series
  nT <- ncol(data_ts) # Number of time step
  
  # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
  constrInd <- (rep(1:nfac, each = ny) > rep(1:ny,  nfac))
  Zmap <- matrix(ncol = nfac, nrow = ny)
  Zmap[constrInd] <- NA
  Zmap[!constrInd] <- 1:sum(!constrInd)
  tmbMap <- list(Z = as.factor(Zmap))
  

  optList = vector(con$nstart * length(con$method), mode = 'list')
  names(optList) = rep(con$method, con$nstart)

  
  for (i in 1:length(optList)) {
    log_re_sp <- runif(ny, -1, 0)
    
    Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
    
    # Set constrained elements to zero
    
    Zinit[constrInd] <- 0
    
    # List of parameters for DFA
    
    tmbPar <-  list(log_re_sp=log_re_sp, Z = Zinit, beta_GAM = matrix(0, ncol = nfac, nrow = ncol(tmbData$X_GAM)))
    #tmbPar <-  list(log_re_sp=log_re_sp, Z = Zinit, beta_GAM = matrix(0.01, ncol = nfac, nrow = ncol(tmbData$X_GAM)))
    #tmbPar$x_sum = matrix(rep(0, nfac), ncol = 1)
    # Make DFA
    tmbObj <- MakeADFun(data = tmbData, parameters = tmbPar, map = tmbMap, random= c('beta_GAM'), DLL= "dfa_GAM", silent = silent)
    optList[[i]] = switch(names(optList)[i],
                          NLMINB = nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = con$maxit, eval.max  =2*con$maxit, rel.tol =  con$reltol)),
                          BFGS = optim(tmbObj$par, tmbObj$fn, tmbObj$gr, method = 'BFGS', control = list(maxit = con$maxit, reltol = con$reltol)),
                          LBFGS = optim(tmbObj$par, tmbObj$fn, tmbObj$gr, method = 'L-BFGS-B', control = list(maxit = con$maxit, factr = con$factr))
    )
    if (names(optList)[i] == 'NLMINB')
      optList[[i]]$value = optList[[i]]$objective
  }
  convergence = sapply(optList, FUN = `[[`, 'convergence')
  print(convergence)
  nll = sapply(optList, FUN = `[[`, 'value')
  print(nll)
  maxgrad = sapply(optList, FUN = function(opt) {max(abs(tmbObj$gr(opt$par))) })
  print(maxgrad)
  eligible = abs(nll - min(nll)) < con$nlldeltatol & convergence == 0 & maxgrad < con$gradtol
  if (!any(eligible)) { 
    eligible = abs(nll - min(nll)) < con$nlldeltatol & convergence == 0 # Currently prioritizes optim convergence over gradient check
    if(!any(eligible))
      eligible = abs(nll - min(nll)) < con$nlldeltatol
  } 

  ind.best =  which.min((nll - 1e6 * sign(min(nll)) * !eligible)) # Return the smallest loglikelihood fit that meets other convergence criteria
  tmbOpt = optList[[ind.best]] 
  # Jonas: I suggest we return the full nll, convergence, and maxgrad vectors at the end of the function. It's useful to have for checking stability of the model.
  
  sdr_all <- sdreport(tmbObj)
  sdr <- summary(sdr_all)
  
  x_sp = matrix(sdr[rownames(sdr) == 'x_sp',1], ncol =23)
  matplot(t(x_sp), type = 'l')
#  x = matrix(sdr[rownames(sdr) == 'x',1], ncol =23)
#  matplot(t(x), type = 'l')
  
   # Check convergence
  conv <- tmbOpt$convergence
  if(tmbOpt$convergence != 0){warning(paste0("Convergence sdissue:", tmbOpt$message))}
  if (tmbOpt$convergence == 0 & maxgrad[ind.best] > con$gradtol) 
    warning(paste0('Optimization converged, but maximum gradient = ', maxgrad[ind.best]))
  # Compute AIC
  
  aic <- AIC.tmb(tmbObj) 
  writeLines(paste('AIC: ', aic))

  return(list(tmbObj = tmbObj, # TMB output
              opt = tmbOpt, # Optimisation from TMB
              data_ts, # Dataset of time series (output)
              data_ts_se, # Dataset of standard error of time series (output)
              data_ts_save, # Dataset of time series (input saved)
              data_ts_save_long, # Dataset of time series (input saved in long format)
              data_ts_se_save, # Dataset of standard error of time series (input saved)
              ny, # Number of time series
              nT, # Number of time step
              aic = aic, # AIC
              conv = conv, # Convergence check
              sdr.sum = sdr, # Summary of the TMB optimisation output
              sdr = sdr_all # Complete TMB optimisation output
              ))
}
