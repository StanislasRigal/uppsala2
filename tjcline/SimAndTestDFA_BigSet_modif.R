## Dynamical Factor Analysis with TMB for ts computed for the entire country

# time-series from the SBBS

ts_bird_se_allcountry <- readRDS("output/ts_bird_se_allcountry.rds")
species_data <- readRDS("output/species_data.rds")

# load functions
source("function_ts.R")

# DFA with TMB from https://github.com/tjcline/dfaTMB/tree/master/SimulationTesting
source('tjcline/DynamicFactorAnalysis_TMB_UsingUnstructuredCorr.R')

# farmland species
ts_bird_se_allcountry_data <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance) & ts_bird_se_allcountry$CI_inf!=0),]
ts_bird_se_allcountry_data <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp != names(which(table(ts_bird_se_allcountry_data$code_sp)<=1)),]
species_sub <- droplevels(species_data[species_data$code_sp %in% c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                                                   "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                                                   "LANCOL","STUVUL","LINCAN","EMBCIT",
                                                                   "PASMON","CORFRU","ANTPRA","EMBHOR"),])
nsim <- 100
sim_agri <- get_sim(ts_bird_se_allcountry_data,
                    species_sub$code_sp,
                    nsim)

Zscore<-function(x){
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

simObs <- t(apply(sim_agri[,-1],1,FUN=Zscore))

library(MASS)
library(MARSS)

# Fit the model with TMB code
myFit<-runDFA(simObs,NumStates=2)

# results total
matplot(t(simObs), pch =20)
matpoints(t(myFit$Estimates$Z %*% myFit$Estimates$u), type = 'l', lwd = 3)
matplot(t(myFit$Estimates$u), type = 'l')

# explore result for each species
matplot(t(simObs[1:nsim,]), pch =20)
matpoints(t(myFit$Estimates$Z[1:nsim,] %*% myFit$Estimates$u), type = 'l', lwd = 3)

data_dfa_plot <- data.dfa.plot(dataset = ts_bird_se_allcountry_data,
                               species = species_sub$code_sp,
                               sim_data = sim_agri,
                               dfa_res = myFit,
                               nsim = nsim)

ggplot(data_dfa_plot[data_dfa_plot$code_sp=="ALAARV",], aes(x = year,y = relative_abundance_std)) + 
  geom_point() +
  geom_pointrange(aes(ymax = relative_abundance_std+1.96*se_std, ymin=relative_abundance_std-1.96*se_std)) +
  geom_line(aes(y = mean_ts_dfa_std)) +
  geom_ribbon(aes(ymax = mean_ts_dfa_std+1.96*sd_ts_dfa_std, ymin=mean_ts_dfa_std-1.96*sd_ts_dfa_std), alpha=0.5) +
  theme_modern()

# test with time as covariable to mimic global random effect

myFit2 <- runDFA(simObs,NumStates=2,EstCovar=T,Covars=matrix(Zscore(1:25),nrow=1))

# test with time as covariable to mimic individual random effect

myFit3 <- runDFA(simObs,NumStates=2,EstCovar=T,indivCovar=T,Covars=matrix(Zscore(1:25),nrow=1))

# test with unconstraaint observation error structure

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
Obs <- dcast(Obs[,c("code_sp","relative_abundance","year")], code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
simObs <- t(apply(Obs[,-1],1,FUN=Zscore))
myFit3 <- runDFA(simObs,NumStates=2,ErrStruc = "UNC")

# simpleDFA

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
y <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs <- dcast(Obs[,c("code_sp","Standard_error","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Standard_error")
y <- as.matrix(y[,-c(1:2)])
obs <- as.matrix(obs[,-c(1:2)])


compile("simpleDFA.cpp")
dyn.load(dynlib("simpleDFA"))

dataTmb <- list(y =y, obs=obs)

nfac = 2 # Number of factors
ny = nrow(y) # Number of time series
nT = ncol(y) # Number of time step

Zinit = matrix(rnorm(ny * nfac), ncol = nfac)
# Set constrained elements to zero
constrInd = rep(1:nfac, each = ny) > rep(1:ny,  nfac)
Zinit[constrInd] = 0
Zinit

tmbPar =  list(logSdO = 0, Z = Zinit,
               x=matrix(c(rep(0, nfac), rnorm(nfac * nT)), ncol = nT+1, nrow = nfac))

# Set up parameter constraints. Elements set to NA will be fixed and not estimated.
Zmap = matrix(ncol = nfac, nrow = ny)
Zmap[constrInd] = NA
Zmap[!constrInd] = 1:sum(!constrInd)
xmap = matrix(ncol = nT + 1, nrow = nfac)
xmap[,1] = NA
xmap[(nfac + 1) : length(tmbPar$x)] = 1:(length(tmbPar$x) - nfac)
tmbMap = list(Z = as.factor(Zmap), x = as.factor(xmap))


tmbObj = MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= "x", DLL= "simpleDFA")

names(tmbObj)

tmbOpt = nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))

tmbOpt$message

sdRep <- summary(sdreport(tmbObj))

sdRep[grepl('Z|sdo', rownames(sdRep)),]


# Get point estimates
#x_hat <- matrix(sdRep[rownames(sdRep)=="x",1], nrow = nfac)

x_hat = (tmbObj$env$parList)()$x

Z_hat = (tmbObj$env$parList)(par=tmbOpt$par)$Z

matplot(t(y), pch =20)

matpoints(t(Z_hat %*% x_hat), type = 'l', lwd = 3)

matplot(t(x_hat), type = 'l')

# Compute AIC

# Only works if obj has been already optimized
AIC.tmb = function(obj, tol = 0.01) {
  # Simple convergence check
  stopifnot(max(abs(obj$gr(obj$env$last.par.best[obj$env$lfixed()]))) < tol)
  as.numeric(2 * obj$env$value.best + 2*sum(obj$env$lfixed()))
}

AIC.tmb(tmbObj)

#aic = 2 * as.numeric(tmbObj$fn(tmbOpt$par)) + 2 * length(tmbOpt$par)
#aic


# Plot rotated trends, see https://atsa-es.github.io/atsa-labs/sec-dfa-rotating-loadings.html
matplot(t(solve(varimax(Z_hat)$rotmat) %*% x_hat), type = 'l')

Z_hat %*% varimax(Z_hat)$rotmat
