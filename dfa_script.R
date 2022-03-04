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
obs_se <- dcast(Obs[,c("code_sp","Standard_error","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Standard_error")
obs_se <- obs_se[,-c(1:2)]
y <- y[,-c(1:2)]
for(i in 1:nrow(obs_se)){
  obs_se[i,] <- Zscore_err(obs_se[i,], y[i,])
}
obs_se <- as.matrix(obs_se) 
y <- as.matrix(t(apply(y,1,Zscore)))

dataTmb <- list(y =y, obs_se=obs_se)

nfac = 2 # Number of factors
ny = nrow(y) # Number of time series
nT = ncol(y) # Number of time step

Zinit = matrix(rnorm(ny * nfac), ncol = nfac)
## Set constrained elements to zero
constrInd = rep(1:nfac, each = ny) > rep(1:ny,  nfac)
Zinit[constrInd] = 0
Zinit

tmbPar =  list(re_sp_para = 1, Z = Zinit,
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


# Test function for farmland birds
species_sub <- droplevels(species_data[species_data$code_sp %in% c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                                                   "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                                                   "LANCOL","STUVUL","LINCAN","EMBCIT",
                                                                   "PASMON","CORFRU","ANTPRA","EMBHOR"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
y <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs_se <- dcast(Obs[,c("code_sp","Standard_error","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Standard_error")

farm_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
farm_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3) # best AIC
farm_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4)

# Forest bird

species_sub <- droplevels(species_data[species_data$code_sp %in% c("ACCNIS","TETBON","TRIOCH","COLOEN",
                                                                   "DENMAJ","DRYMAR","PICVIR","JYNTOR",
                                                                   "DRYMIN","PICTRI","NUCCAR","GARGLA",
                                                                   "PERATE","LOPCRI","POEPAL","POEMON",
                                                                   "SITEUR","CERFAM","TURVIS","PHOPHO",
                                                                   "PHYCOL","PHYSIB","REGREG","FICHYP",
                                                                   "ANTTRI","COCCOC","SPISPI","PYRPYR","EMBRUS"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

y <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs_se <- dcast(Obs[,c("code_sp","Standard_error","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Standard_error")

forest_nfac2 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 2)
forest_nfac3 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 3) # best AIC
forest_nfac4 <- make_dfa(data_ts = y, data_ts_se = obs_se, nfac = 4)

# Simulation to analyse parameter influence on dfa fit

n_y <- 25 # number of year
y <- data.frame(t(rep(NA,(n_y+1))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 15 # number of species
sd_rand <- 0.1# sd of normal draw to simulate observation standard error
for(i in 1:n_sp){
  set.seed(i)
  a <- rnorm(1,0,0.1) # linear regression coef
  b <- runif(1, -1, 1) # intercept
  d <- rnorm(n_y) # noise
  set.seed(10*i)
  f <- abs(rnorm(n_y,0,sd_rand)) # standard deviation
  y[i,1] <- obs_se[i,1] <- paste0("SP",i)
  y[i,2:(n_y+1)] <- a*c(1:n_y)+b+d
  y[i,2] <- 0
  y[i,2:(n_y+1)] <- (y[i,2:(n_y+1)]+10)/10
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

