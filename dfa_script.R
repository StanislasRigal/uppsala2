## Dynamical Factor Analysis with TMB and uncertainty

# Load packages

source("package_used.R")

# Load functions
source("function_ts.R")

# DFA with TMB
source('function_dfa.R')

# Time-series from the SBBS

ts_bird_se_allcountry <- readRDS("output/ts_bird_se_allcountry.rds")
species_data <- readRDS("output/species_data.rds")

# Example with farmland species

ts_bird_se_allcountry_data <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance) & ts_bird_se_allcountry$CI_inf!=0),]
ts_bird_se_allcountry_data <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp != names(which(table(ts_bird_se_allcountry_data$code_sp)<=1)),]
species_sub <- droplevels(species_data[species_data$code_sp %in% c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                                                   "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                                                   "LANCOL","STUVUL","LINCAN","EMBCIT",
                                                                   "PASMON","CORFRU","ANTPRA","EMBHOR"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
y <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs <- dcast(Obs[,c("code_sp","Standard_error","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Standard_error")
obs <- obs[,-c(1:2)]
y <- y[,-c(1:2)]
for(i in 1:nrow(obs)){
  obs[i,] <- Zscore_err(obs[i,], y[i,])
}
obs <- as.matrix(obs) 
y <- as.matrix(t(apply(y,1,Zscore)))

dataTmb <- list(y =y, obs=obs)

nfac = 2 # Number of factors
ny = nrow(y) # Number of time series
nT = ncol(y) # Number of time step

Zinit = matrix(rnorm(ny * nfac), ncol = nfac)
## Set constrained elements to zero
constrInd = rep(1:nfac, each = ny) > rep(1:ny,  nfac)
Zinit[constrInd] = 0
Zinit

tmbPar =  list(logSdO = 0, Z = Zinit,
               x=matrix(c(rep(0, nfac), rnorm(nfac * nT)), ncol = nT+1, nrow = nfac))

## Set up parameter constraints. Elements set to NA will be fixed and not estimated.
Zmap = matrix(ncol = nfac, nrow = ny)
Zmap[constrInd] = NA
Zmap[!constrInd] = 1:sum(!constrInd)
xmap = matrix(ncol = nT + 1, nrow = nfac)
xmap[,1] = NA
xmap[(nfac + 1) : length(tmbPar$x)] = 1:(length(tmbPar$x) - nfac)
tmbMap = list(Z = as.factor(Zmap), x = as.factor(xmap))

tmbObj = MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= "x", DLL= "simpleDFA")
tmbOpt = nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
tmbOpt$message
sdRep <- summary(sdreport(tmbObj))
sdRep[grepl('Z|sdo', rownames(sdRep)),]

## Get point estimates
#x_hat <- matrix(sdRep[rownames(sdRep)=="x",1], nrow = nfac)

x_hat = (tmbObj$env$parList)()$x
Z_hat = (tmbObj$env$parList)(par=tmbOpt$par)$Z

matplot(t(y), pch =20)
matpoints(t(Z_hat %*% x_hat[,-1]), type = 'l', lwd = 3)
matplot(t(x_hat[,-1]), type = 'l')

## Compute AIC

AIC.tmb(tmbObj)

#aic = 2 * as.numeric(tmbObj$fn(tmbOpt$par)) + 2 * length(tmbOpt$par)
#aic

## Plot rotated trends, see https://atsa-es.github.io/atsa-labs/sec-dfa-rotating-loadings.html
matplot(t(solve(varimax(Z_hat)$rotmat) %*% x_hat[,-1]), type = 'l')
Z_hat %*% varimax(Z_hat)$rotmat

# Generalise DFA

make_dfa <- function(data_ts, # dataset of time series
                     data_ts_se, # dataset of standard error of time series 
                     nfac) # number of trends for the DFA
  {
  
  # Save input data for plot
  
  data_ts_save <- data_ts
  data_ts_se_save <- data_ts_se
  
  # Remove potential column of ts names
  
  data_ts <- data_ts %>% select_if(Negate(is.character))
  data_ts_se <- data_ts_se %>% select_if(Negate(is.character))
  
  # Remove first year if all se = 0
  
  if(mean(data_ts_se[,1], na.rm=T)==0 & sd(data_ts_se[,1], na.rm=T)==0){
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
  
  aic <- AIC.tmb(tmbObj)
  print(aic)
  
  # Plot rotated trends, see https://atsa-es.github.io/atsa-labs/sec-dfa-rotating-loadings.html
  matplot(t(y), pch =20)
  matpoints(t(Z_hat %*% x_hat[,-1]), type = 'l', lwd = 3)
  matplot(t(x_hat[,-1]), type = 'l')
  matplot(t(solve(varimax(Z_hat)$rotmat) %*% x_hat[,-1]), type = 'l')
  Z_hat %*% varimax(Z_hat)$rotmat
}
