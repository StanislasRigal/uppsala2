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

# Clean data 

ts_bird_se_allcountry_data <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance) & ts_bird_se_allcountry$CI_inf!=0),]
#ts_bird_se_allcountry_data <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp != names(which(table(ts_bird_se_allcountry_data$code_sp)<=1)),]

# Example with farmland species

species_sub <- droplevels(species_data[species_data$code_sp %in% c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                                                   "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                                                   "LANCOL","STUVUL","LINCAN","EMBCIT",
                                                                   "PASMON","CORFRU","ANTPRA","EMBHOR"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
y <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs_se <- dcast(Obs[,c("code_sp","Standard_error","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Standard_error")
#obs_se <- obs_se[,-c(1:2)]
#y <- y[,-c(1:2)]
#for(i in 1:nrow(obs_se)){
#  obs_se[i,] <- Zscore_err(obs_se[i,], y[i,])
#}
#obs_se <- as.matrix(obs_se) 
#y <- as.matrix(t(apply(y,1,Zscore)))

# Remove column of ts names

y <- y %>% select_if(Negate(is.character))
obs_se <- obs_se %>% select_if(Negate(is.character))

dataTmb <- list(y =log(as.matrix(y)),
                obs_se = as.matrix(obs_se/y))

nfac <- 2 # Number of factors
ny <- nrow(y) # Number of time series
nT <- ncol(y) # Number of time step
set.seed(1) 
log_re_sp <- runif(ny, -1, 0)
Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
## Set constrained elements to zero
constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
Zinit[constrInd] <- 0
Zinit

tmbPar <- list(log_re_sp=log_re_sp, Z = Zinit,
               x=matrix(c(rep(0, nfac), rnorm(nfac * (nT - 1))),
                        ncol = nT, nrow = nfac))
## Set up parameter constraints. Elements set to NA will be fixed and not estimated.
Zmap <- matrix(ncol = nfac, nrow = ny)
Zmap[constrInd] <- NA
Zmap[!constrInd] <- 1:sum(!constrInd)
xmap <- matrix(ncol = nT, nrow = nfac)
xmap[,1] <- NA
xmap[(nfac + 1) : length(tmbPar$x)] <- 1:(length(tmbPar$x) - nfac)
#log_re_spmap <- matrix(ncol = ny, nrow = 1)
#log_re_spmap[,1] <- NA
#log_re_spmap[nfac : ny] <- 1:(ny - nfac+1)
tmbMap <- list(#log_re_sp = as.factor(log_re_spmap),
               Z = as.factor(Zmap), x = as.factor(xmap))

tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= "x", DLL= "dfa_model_se")
tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
tmbOpt$message
sdRep <- summary(sdreport(tmbObj))
sdRep[grepl('Z|sdo', rownames(sdRep)),]

## Get point estimates
#x_hat <- matrix(sdRep[rownames(sdRep)=="x",1], nrow = nfac)

x_hat = (tmbObj$env$parList)()$x
Z_hat = (tmbObj$env$parList)(par=tmbOpt$par)$Z

matplot(t(y), pch =20)
matpoints(t(matrix(exp(sdRep[rownames(sdRep)=="x_sp",1]), nrow=ny)), type = 'l', lwd = 3)
# equivalent to
#matplot(t(y), pch =20)
#matpoints(t(exp(Z_hat %*% x_hat)), type = 'l', lwd = 3)
matplot(t(x_hat[,-1]), type = 'l')

## Compute AIC

AIC.tmb(tmbObj, dontCount = 0)

#aic = 2 * as.numeric(tmbObj$fn(tmbOpt$par)) + 2 * length(tmbOpt$par)
#aic

## Plot rotated trends, see https://atsa-es.github.io/atsa-labs/sec-dfa-rotating-loadings.html
matplot(t(solve(varimax(Z_hat)$rotmat) %*% x_hat[,-1]), type = 'l')
Z_hat %*% varimax(Z_hat)$rotmat

# Simulate the data (same code as above)
set.seed(553)
mSimTmb <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= "x", DLL= "dfa_model_se", silent = TRUE)
simTmb <- mSimTmb$simulate(complete = T)
# Fit the model to check estimation
mCheckTmb <- MakeADFun(list(y = simTmb$y, obs_se = simTmb$obs_se), # Note difference in this argument
                       tmbPar, map=tmbMap, random= "x", DLL= "dfa_model_se", silent = TRUE)
fCheckTmb <- nlminb(mCheckTmb$par, mCheckTmb$fn, mCheckTmb$gr)
fCheckTmb$message

nrep <- 200 # Number of replicates
set.seed(999) # Set the random seed to be able to reproduce the simulation
# Create matrix to keep track of replicate test statistic values
repsT <- matrix(NA, nrow=nrep, ncol=2*nfac)
repsT2 <- matrix(NA, nrow=nrep, ncol=2*ny)
name_vec <- c()
for(i in 1:nfac){
  name_vec <- c(name_vec,paste0("Mean_x",i),paste0("SD_x",i))
}
colnames(repsT) <- name_vec

name_vec <- c()
for(i in 1:ny){
  name_vec <- c(name_vec,paste0("Mean_sp",i),paste0("SD_sp",i))
}
colnames(repsT2) <- name_vec

# Parameter to use to simulate
parSp <- tmbObj$env$par
# Set to estimated values
parSp[1:length(tmbOpt$par)] <- tmbOpt$par
for(i in 1:nrep){ # For each replicate
  yrep <- mSimTmb$simulate(complete=T, par=parSp)$x # simulate observations
  for(j in 1:nfac){
    repsT[i, (2*j-1)] <- mean(yrep[j,])
    repsT[i, (2*j)] <- sd(yrep[j,])
  }
  yrep <- mSimTmb$simulate(complete=T, par=parSp)$y # simulate observations
  for(j in 1:ny){
    repsT2[i, (2*j-1)] <- mean(yrep[j,])
    repsT2[i, (2*j)] <- sd(yrep[j,])
  }
}

par(mfrow=c(2,2))
for(i in 1:nfac){
  hist(repsT[,(2*i-1)])
  abline(v=mean(x_hat[i,]), col="hotpink", lwd=4)
  hist(repsT[,(2*i)])
  abline(v=sd(x_hat[i,]), col="hotpink", lwd=4)
}

for(i in 1:ny){
  hist(repsT2[,(2*i-1)])
  abline(v=mean(Z_hat[i,]), col="hotpink", lwd=4)
  hist(repsT2[,(2*i)])
  abline(v=sd(Z_hat[i,]), col="hotpink", lwd=4)
}
par(mfrow=c(1,1))
#sim <- replicate(50, {
#  simdata <- tmbObj$simulate(par=tmbObj$env$par, complete=TRUE)
#  mSimTmb <- MakeADFun(simdata, tmbPar, map=tmbMap, random= "x", DLL= "dfa_model_se", silent = TRUE)
#  nlminb(mSimTmb$par, mSimTmb$fn, mSimTmb$gr)$par
#})


# Check model assumptions

source("function_onesteppred_modif.R")

tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= "x", DLL= "dfa_model_se")
tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
res2pTmb <- oneStepPredict2(tmbObj, observation.name ="y",
                           data.term.indicator = "keep",
                           method="oneStepGaussian",
                           trace=FALSE, discrete = F)

qqnorm(res2pTmb$residual,main="", pch=16, cex=0.8)
qqline(res2pTmb$residual, col="hotpink", lwd=2)
acf(res2pTmb$residual)

# Check Laplace approximation

source("function_checkConsistency_modif.R")
set.seed(1)
checkConsistency(tmbObj)



# Test function for farmland birds

species_sub <- species_farm <-  droplevels(species_data[species_data$code_sp %in% c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                                                   "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                                                   "LANCOL","STUVUL","LINCAN","EMBCIT",
                                                                   "PASMON","CORFRU","ANTPRA","EMBHOR"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
y <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

farm_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
farm_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3) # best AIC
farm_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4)
farm_nfac5 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 5)

# Forest bird

species_sub <- species_forest <- droplevels(species_data[species_data$code_sp %in% c("ACCNIS","TETBON","TRIOCH","COLOEN",
                                                                   "DENMAJ","DRYMAR","PICVIR","JYNTOR",
                                                                   "DRYMIN","PICTRI","NUCCAR","GARGLA",
                                                                   "PERATE","LOPCRI","POEPAL","POEMON",
                                                                   "SITEUR","CERFAM","TURVIS","PHOPHO",
                                                                   "PHYCOL","PHYSIB","REGREG","FICHYP",
                                                                   "ANTTRI","COCCOC","SPISPI","PYRPYR","EMBRUS"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

y <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

forest_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
forest_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3)
forest_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4) # best AIC
forest_nfac5 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 5)

# All common birds

species_sub <- droplevels(species_data[species_data$code_sp %in%
                                         levels(as.factor(species_data$code_sp))[c(1:24,26:48,50:65,67:93,95:116,118:137,139:144,147:169,171:187,189:253)],])
ab_sp2 <- readRDS("output/ab_sp2.rds")

# select species representing more than 99 % of the total abundance
Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
Obs <- droplevels(Obs[Obs$code_sp %in% levels(as.factor(droplevels(ab_sp2[ab_sp2$perc_cum<0.99,])$code_sp)),])
#Obs <- droplevels(Obs[Obs$uncertanity_reason!="too rare species",])

species_all <- droplevels(species_data[species_data$code_sp %in%
                                         levels(as.factor(Obs$code_sp)),])

y <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

all_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
all_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3)
all_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4) 
all_nfac5 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 5)

# Distinguish north populations from south populations
ts_bird_se_byecoreg <- readRDS("output/ts_bird_se_byecoreg.rds")
ts_bird_se_byecoreg_data <- ts_bird_se_byecoreg[which(!is.na(ts_bird_se_byecoreg$relative_abundance)),]

ts_bird_se_byecoreg_data$code_sp <- paste0(ts_bird_se_byecoreg_data$code_sp,
                                           sep="_", ts_bird_se_byecoreg_data$ecoreg)

species_sub <- expand.grid(code_sp = c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                       "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                       "LANCOL","STUVUL","LINCAN","EMBCIT",
                                       "PASMON","CORFRU","ANTPRA","EMBHOR"),
                           eco_reg = c("north","south"))
species_sub <- merge(species_sub,species_data,by="code_sp",all.x=T)
species_sub$code_sp_eco <- paste0(species_sub$code_sp,
                                  sep="_", species_sub$eco_reg)
species_sub$name_long_eco <- paste0(species_sub$name_long,
                                    sep="_", species_sub$eco_reg)
species_sub <- species_farm_eco <- merge(species_sub,
                     data.frame(code_sp_eco=levels(as.factor(ts_bird_se_byecoreg_data[,"code_sp"]))),
                     by="code_sp_eco")

Obs <- ts_bird_se_byecoreg_data[ts_bird_se_byecoreg_data$code_sp %in% species_sub$code_sp_eco,]

Obs <- droplevels(Obs[Obs$uncertanity_reason!="too rare species",])

y <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs_se <- dcast(Obs[,c("code_sp","Standard_error","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Standard_error")

farm_eco_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
farm_eco_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3)
farm_eco_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4) # best without convergence issue
farm_eco_nfac5 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 5)

# forest by eco region

species_sub <- expand.grid(code_sp = c("TRIOCH","DENMAJ","DRYMAR","JYNTOR",
                                       "GARGLA","PERATE","LOPCRI","POEMON",
                                       "CERFAM","TURVIS","PHOPHO","PHYCOL",
                                       "PHYSIB","REGREG","FICHYP","ANTTRI",
                                       "SPISPI","PYRPYR"),
                           eco_reg = c("north","south"))
species_sub <- merge(species_sub,species_data,by="code_sp",all.x=T)
species_sub$code_sp_eco <- paste0(species_sub$code_sp,
                                  sep="_", species_sub$eco_reg)
species_sub$name_long_eco <- paste0(species_sub$name_long,
                                    sep="_", species_sub$eco_reg)
species_sub <- species_forest_eco <- merge(species_sub,
                     data.frame(code_sp_eco=levels(as.factor(ts_bird_se_byecoreg_data[,"code_sp"]))),
                     by="code_sp_eco")

Obs <- ts_bird_se_byecoreg_data[ts_bird_se_byecoreg_data$code_sp %in% species_sub$code_sp_eco,]

Obs <- droplevels(Obs[Obs$uncertanity_reason!="too rare species",])

y <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs_se <- dcast(Obs[,c("code_sp","Standard_error","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Standard_error")

forest_eco_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
forest_eco_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3) # best without convergence issue
forest_eco_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4) 
forest_eco_nfac5 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 5)

# Simulation to analyse parameter influence on dfa fit

n_y <- 25 # number of year
y <- data.frame(t(rep(NA,(n_y+1))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 20 # number of species
#sd_rand <- 0.05# sd of normal draw to simulate observation standard error
sd_rand <-0.15
for(i in 1:n_sp){
  set.seed(i)
  a1 <- rnorm(1,0,0.1) # linear regression coef
  a2 <- rnorm(1,0,0.0005) # second order coef
  a3 <- rnorm(1,0,0.0005) # third order coef
  b <- runif(1, -1, 1) # intercept
  d <- rnorm(n_y,0,2) # noise
  set.seed(100*i)
  #f <- abs(rnorm(n_y,0,sd_rand)) # standard deviation
  f <- runif(n_y,0.1,sd_rand)
  y[i,1] <- obs_se[i,1] <- paste0("SP",i)
  y[i,2:(n_y+1)] <- a3*(c(1:n_y))^3+a2*c(1:n_y)^2+a1*c(1:n_y)+b+d
  y[i,2] <- 0
  y[i,2:(n_y+1)] <- abs((y[i,2:(n_y+1)]+10)/10)
  obs_se[i,2:(n_y+1)] <- f
  obs_se[i,2] <- 0
}

y <- data.table(y)
obs_se <- data.table(obs_se)
names(y) <- names(obs_se) <- c("code_sp",1:n_y)

test_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
test_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3)
test_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4)
test_nfac5 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 5)

best_aic <- matrix(nrow=4,ncol=100)
for(i in 1:100){
  test_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2, rand_seed = i)
  test_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3, rand_seed = i)
  test_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4, rand_seed = i)
  test_nfac5 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 5, rand_seed = i)
  best_aic[1,i] <- test_nfac2[[7]]
  best_aic[2,i] <- test_nfac3[[7]]
  best_aic[3,i] <- test_nfac4[[7]]
  best_aic[4,i] <- test_nfac5[[7]]
}
seed_2 <- which.min(best_aic[1,])
seed_3 <- which.min(best_aic[2,])
seed_4 <- which.min(best_aic[3,])
seed_5 <- which.min(best_aic[4,])

# Analyse DFA results

farm_nfac3_val <- dcast(farm_nfac3[[3]], code_sp~variable, value.var = "value")
farm_nfac3_se <- dcast(farm_nfac3[[3]], code_sp~variable, value.var = "se.value")
names(farm_nfac3_se)[2:4] <- c("se_X1","se_X2","se_X3")
farm_nfac3_load <- merge(farm_nfac3_val, farm_nfac3_se, by="code_sp")

ggplot(farm_nfac3_load, aes(x=X1, y=X2)) + geom_point() +
  geom_pointrange(aes(ymax = X2+se_X2, ymin=X2-se_X2)) +
  geom_errorbarh(aes(xmax = X1+se_X1, xmin=X1-se_X1)) +
  theme_modern()

forest_nfac4_val <- dcast(forest_nfac4[[3]], code_sp~variable, value.var = "value")
forest_nfac4_se <- dcast(forest_nfac4[[3]], code_sp~variable, value.var = "se.value")
names(forest_nfac4_se)[2:5] <- c("se_X1","se_X2","se_X3","se_X4")
forest_nfac4_load <- merge(forest_nfac4_val, forest_nfac4_se, by="code_sp")

ggplot(forest_nfac4_load, aes(x=X1, y=X2)) + geom_point() +
  geom_pointrange(aes(ymax = X2+se_X2, ymin=X2-se_X2)) +
  geom_errorbarh(aes(xmax = X1+se_X1, xmin=X1-se_X1)) +
  theme_modern()

farm_eco_nfac3_val <- dcast(farm_eco_nfac3[[3]], code_sp~variable, value.var = "value")
farm_eco_nfac3_se <- dcast(farm_eco_nfac3[[3]], code_sp~variable, value.var = "se.value")
names(farm_eco_nfac3_se)[2:4] <- c("se_X1","se_X2","se_X3")
farm_eco_nfac3_load <- merge(farm_eco_nfac3_val, farm_eco_nfac3_se, by="code_sp")
farm_eco_nfac3_load$eco_reg <- sub(".*_", "", farm_eco_nfac3_load$code_sp)

ggplot(farm_eco_nfac3_load, aes(x=X1, y=X2, col=eco_reg)) + geom_point() +
  geom_pointrange(aes(ymax = X2+se_X2, ymin=X2-se_X2)) +
  geom_errorbarh(aes(xmax = X1+se_X1, xmin=X1-se_X1)) +
  theme_modern()

forest_eco_nfac3_val <- dcast(forest_eco_nfac3[[3]], code_sp~variable, value.var = "value")
forest_eco_nfac3_se <- dcast(forest_eco_nfac3[[3]], code_sp~variable, value.var = "se.value")
names(forest_eco_nfac3_se)[2:4] <- c("se_X1","se_X2","se_X3")
forest_eco_nfac3_load <- merge(forest_eco_nfac3_val, forest_eco_nfac3_se, by="code_sp")
forest_eco_nfac3_load$eco_reg <- sub(".*_", "", forest_eco_nfac3_load$code_sp)

ggplot(forest_eco_nfac3_load, aes(x=X1, y=X2, col=eco_reg)) + geom_point() +
  geom_pointrange(aes(ymax = X2+se_X2, ymin=X2-se_X2)) +
  geom_errorbarh(aes(xmax = X1+se_X1, xmin=X1-se_X1)) +
  theme_modern()

# reduce dimension with NMDS

farm_NMDS <- monoMDS(vegdist(farm_nfac3_val[,-1], method = "euclidean"), k=2)
stressplot(farm_NMDS)
ggplot(data.frame(farm_NMDS$points), aes(x=MDS1, y=MDS2)) +
  geom_point() + theme_modern()

farm_eco_NMDS <- monoMDS(vegdist(farm_eco_nfac3_val[,-1], method = "euclidean"), k=2)
stressplot(farm_eco_NMDS)
ggplot(data.frame(farm_eco_NMDS$points,eco_reg=sub(".*_", "", farm_eco_nfac3_val$code_sp)),
       aes(x=MDS1, y=MDS2, col=eco_reg)) +
  geom_point() + theme_modern()

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
    test<-mean_trend[,c(1,i+1, i+nb_group+1)]
    test$Index_SE<-test[,3]
    test$Index<-test[,2]
    ggplot(test, aes(x=year, y=Index)) +
      geom_line() +
      geom_ribbon(aes(ymin=Index-Index_SE,ymax=Index+Index_SE),alpha=0.7, col="black",fill="white")+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)
  }), names(mean_trend)[2:(nb_group+1)])
  
  centroids_data<-tibble(x=centroids$PC1,
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
  
  return(list(kmeans_res,final_plot,mean_trend))
}

group_test <- group_from_dfa(test_nfac3, data.frame(code_sp=paste0("SP",1:n_sp), name_long=paste0("Species",1:n_sp)))
group_farm <- group_from_dfa(farm_nfac3,species_farm)
group_forest <- group_from_dfa(forest_nfac4,species_forest)
group_all <- group_from_dfa(all_nfac2,species_all)
group_farm_eco <- group_from_dfa(farm_eco_nfac4,species_farm_eco, eco_reg = T)
group_forest_eco <- group_from_dfa(forest_eco_nfac3,species_forest_eco, eco_reg = T)
