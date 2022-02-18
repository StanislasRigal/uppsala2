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

myFit2 <- runDFA(simObs,NumStates=2,EstCovar=T,indivCovar=T,Covars=matrix(Zscore(1:25),nrow=1))
