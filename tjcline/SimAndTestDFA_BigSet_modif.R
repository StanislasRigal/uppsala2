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

#RsimUNC<-as.matrix(read.csv("tjcline/BigCovarianceMatrix.csv",header=F))

Zscore<-function(x){
  return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
}

simObs <- t(apply(sim_agri[,-1],1,FUN=Zscore))

#These are the loadings for the simulation
#Zsim<-matrix(runif(nrow(RsimUNC)*2,0,1),ncol=2)#c(0.3,0.5,0.8)
#Zsim[,1]<-sort(Zsim[,1])
#Zsim[,2]<-rev(sort(Zsim[,2]))

#Simulated randomwalk common trend
#uSim<-matrix(rep(rep(NA,100),2),nrow=2)
#uSim[,1]<--3
#for(i in 2:100){
#  uSim[,i]<-uSim[,i-1]+rnorm(2,0,1)
#}
#uSim<-t(apply(uSim,1,FUN=Zscore))

#plot(uSim)

#simulate a ts for air temp. I used an autocorrelated process but not random walk
#airSim<-rep(NA,100)
#airSim[1]<-0
#for(i in 2:100){
#  airSim[i]<-0.7*airSim[i-1]+rnorm(1,0,1)
#}
#airSim<-Zscore(airSim)

#plot(airSim,type='l')

#AirTemperature loadings for simulation
#Dsim<-runif(nrow(RsimUNC),0,1)#c(0.1,0.6,0.9)

library(MASS)

#R matrices for simulation. Testing all three common forms.
#RsimDE<-diag(0.1,3)
#RsimDUE<-diag(c(0.1,0.2,0.4),3)
#RsimUNC<-matrix(c(0.004328536,0.025618901,0.023676532,0.025618901,0.1516287,0.1405050,0.023676532,0.1405050,0.38084122),nrow=3,ncol=3)

#RsimUNC<-matrix(c(0.038364445,-0.001262234,-0.079377036,-0.001262234,0.067458194,0.013427599,-0.07937704,0.01342760,0.16810893),nrow=3,ncol=3)
#Simulate data series
#simObs<- matrix(Zsim,ncol=2) %*% uSim + matrix(Dsim,ncol=1) %*% airSim + t(mvrnorm(100,rep(0,nrow(RsimUNC)),Sigma=RsimUNC))
#simObs[sample.int(300,size=50)]<-NA

#plot(simObs[1,],type='l',ylim=c(-4,4))
#points(simObs[2,],type='l',col='red')
#points(simObs[3,],type='l',col='blue')

library(MARSS)
#Fit the model using MARSS
#marssFit<-MARSS(simObs,model=list(m=1,R='unconstrained'),form='dfa',control=list(maxit=500))
#marssPred<-coef(marssFit,type='matrix')$Z %*% marssFit$states + coef(marssFit,type='matrix')$D %*% airSim

#par(mfrow=c(3,2),mar=c(2,2,1,1))
#ScaleFacM<-max(Mod(marssFit$states[1,]))/max(Mod(uSim))
#plot(uSim)
#points(marssFit$states[1,]/ScaleFacM,type='l')

#plot(marssFit$states[1,]/ScaleFacM~uSim)
#abline(0,1)

#coef(marssFit)$Z*ScaleFacM
#Zsim

#plot(coef(marssFit)$Z*ScaleFacM~Zsim,ylim=c(0,1),xlim=c(0,1))
#abline(0,1)

#coef(marssFit)$D
#Dsim

#plot(coef(marssFit)$D~Dsim,ylim=c(0,1),xlim=c(0,1))
#abline(0,1)

#plot(simObs[2,])
#points(marssPred[2,],type='l')

#plot(coef(marssFit)$R~RsimUNC[lower.tri(RsimUNC,diag=T)],ylim=c(min(RsimUNC),max(RsimUNC)),xlim=c(min(RsimUNC),max(RsimUNC)))
#abline(0,1)


# Fit the model with TMB code
myFit<-runDFA(simObs,NumStates=2)

#par(mfrow=c(2,2),mar=c(4,4,1,1))
#ScaleFac1<-max(Mod(myFit$Estimates$u[1,]))/max(Mod(uSim[2,]))
#ScaleFac2<-max(Mod(myFit$Estimates$u[2,]))/max(Mod(uSim[1,]))
#plot(uSim[1,],main='Shared Trend')
#plot(myFit$Estimates$u[2,]/ScaleFac1,type='l')
#points(uSim[2,],col='red')
#points(myFit$Estimates$u[1,]/ScaleFac2,type='l',col='red')
#legend('top',legend=c('Simulated','ModelFit'),pch=c(1,NA),lty=c(NA,1))

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

#plot(myFit$Estimates$u[1,]/ScaleFac~uSim,xlab='si')
#abline(0,1)

#myFit$Estimates$Z*ScaleFac1
#Zsim

#plot(myFit$Estimates$Z[,1]*ScaleFac2~myFit$Estimates$Z[,2],ylab='Loading 1',xlab='Loading 2',main='Loadings(Z)')
#abline(0,1)
#points(myFit$Estimates$Z[,2]*ScaleFac1~Zsim[,1])


#myFit$Estimates$D
#Dsim
#plot(myFit$Estimates$D~Dsim,ylab='Estimated',xlab='Simulated',main='Covariates(D)')
#abline(0,1)

#plot(myFit$Estimates$R~RsimUNC,ylab='Estimated',xlab='Simulated',main='CovarianceMatrix(R)')
#abline(0,1)

#myFit$Estimates$R
#coef(marssFit,type='matrix')$R
#RsimUNC



