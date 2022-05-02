## Dynamical Factor Analysis with TMB and uncertainty

# Load packages

source("package_used.R")

# Load functions
source("function_ts.R")

# DFA with TMB
source('function_dfa_test.R')

# Time-series from the SBBS

ts_bird_se_allcountry <- readRDS("output/ts_bird_se_allcountry.rds")
species_data <- readRDS("output/species_data.rds")

# Clean data 

ts_bird_se_allcountry_data <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance) & ts_bird_se_allcountry$CI_inf!=0),]


# Farmland birds

species_sub <- species_farm <-  droplevels(species_data[species_data$code_sp %in% c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                                                   "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                                                   "LANCOL","STUVUL","LINCAN","EMBCIT",
                                                                   "PASMON","CORFRU","ANTPRA","EMBHOR"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
y_farm <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se_farm <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")


# Forest bird

species_sub <- species_forest <- droplevels(species_data[species_data$code_sp %in% c("ACCNIS","TETBON","TRIOCH","COLOEN",
                                                                   "DENMAJ","DRYMAR","PICVIR","JYNTOR",
                                                                   "DRYMIN","PICTRI","NUCCAR","GARGLA",
                                                                   "PERATE","LOPCRI","POEPAL","POEMON",
                                                                   "SITEUR","CERFAM","TURVIS","PHOPHO",
                                                                   "PHYCOL","PHYSIB","REGREG","FICHYP",
                                                                   "ANTTRI","COCCOC","SPISPI","PYRPYR","EMBRUS"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

y_forest <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se_forest <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

# All common birds

## Selection by abundance (species representing 99 % of the total abundance)
species_sub <- droplevels(species_data[species_data$code_sp %in%
                                         levels(as.factor(species_data$code_sp))[c(1:24,26:48,50:65,67:93,95:116,118:137,139:144,147:169,171:187,189:253)],])
ab_sp2 <- readRDS("output/ab_sp2.rds")

## Or selecting the 107 species in the SBS report

species_sub <- species_all <- droplevels(species_data[species_data$code_sp %in% c("ACCNIS","TETBON","TRIOCH","COLOEN", # forest
                                                                                  "DENMAJ","DRYMAR","PICVIR","JYNTOR",
                                                                                  "DRYMIN","PICTRI","NUCCAR","GARGLA",
                                                                                  "PERATE","LOPCRI","POEPAL","POEMON",
                                                                                  "SITEUR","CERFAM","TURVIS","PHOPHO",
                                                                                  "PHYCOL","PHYSIB","REGREG","FICHYP",
                                                                                  "ANTTRI","COCCOC","SPISPI","PYRPYR",
                                                                                  "EMBRUS",
                                                                                  "VANVAN","NUMARQ","ALAARV","HIRRUS", # farmland
                                                                                  "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                                                                  "LANCOL","STUVUL","LINCAN","EMBCIT",
                                                                                  "PASMON","CORFRU","ANTPRA","EMBHOR",
                                                                                  "PODCRI","ARDCIN","ANAPLA","TADTAD",# others
                                                                                  "CYGOLO","BUTBUT","CIRAER","LYRTET",
                                                                                  "PHACOL","GRUGRU","GALCHL","FULATR",
                                                                                  "HAEOST","PLUAPR","GALGAL","NUMARQ",
                                                                                  "TRIGLA","ACTHYP","TRITOT","TRINEB",
                                                                                  "CHRRID","COLPAL","STRDEC","CUCCAN",
                                                                                  "APUAPU","JYNTOR","LULARB","DELURB",
                                                                                  "CORCOX","CORCOR","COLMON","PICPIC",
                                                                                  "AEGCAU","PARMAJ","CYACAE","TROTRO",
                                                                                  "TURPIL","TURPHI","TURILI","TURMER",
                                                                                  "OENOEN","LUSLUS","ERIRUB","LOCNAE",
                                                                                  "ACRSCI","ACRPAL","ACRSCH","HIPICT",
                                                                                  "SYLATR","SYLBOR","SYLCUR","PHYTRO",
                                                                                  "MUSSTR","PRUMOD","MOTALB","MOTCIN",
                                                                                  "CHLCHL","CARCAR","ACAFLA","CARERY",
                                                                                  "FRICOE","FRIMON","EMBSCH","PASDOM"),])


Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

# if selection by abundance:
#Obs <- droplevels(Obs[Obs$code_sp %in% levels(as.factor(droplevels(ab_sp2[ab_sp2$perc_cum<0.99,])$code_sp)),])

species_all <- droplevels(species_data[species_data$code_sp %in%
                                         levels(as.factor(Obs$code_sp)),])



y_all <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se_all <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

# Farmland north populations and south populations
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

y_farm_eco <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs_se_farm_eco <- dcast(Obs[,c("code_sp","Standard_error","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Standard_error")

# Farmland north populations and south populations

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

y_forest_eco <- dcast(Obs[,c("code_sp","relative_abundance","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance")
obs_se_forest_eco <- dcast(Obs[,c("code_sp","Standard_error","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Standard_error")


# Simulation to analyse parameter influence on dfa fit
# General framework with linear process
n_y <- 25 # number of year
y <- data.frame(t(rep(NA,(n_y+1))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 20 # number of species
sd_rand <-0.01
for(i in 1:n_sp){
  set.seed(i)
  y[i,1] <- obs_se[i,1] <- sprintf("SP%03d",i)
  y[i,2:(n_y+1)] <- SimTs(M=1, Tslope=runif(1,-0.002,0.002), Isd=0.015, Irange=0.03, Srange=0.5, Rsd=0.05, 
                          Rrange=0.1, breaks=NULL, abrupt=F, n=n_y, start=c(1998, 1), freq=1)[,1]
  #y[i,-1] <- y[i,-1]/y[i,2]
  obs_se[i,2:(n_y+1)] <- runif(n_y,sd_rand,2*sd_rand)
  #obs_se[i,2] <- 0
}

y_rand <- data.table(y)
obs_se_rand <- data.table(obs_se)
names(y_rand) <- names(obs_se_rand) <- c("code_sp",1:n_y)
species_rand <- data.frame(name_long=sprintf("species %03d",1:nrow(y_rand)), code_sp=y_rand$code_sp)


# General framework with random walk process
n_y <- 25 # number of year
y <- data.frame(t(rep(NA,(n_y+1))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 20 # number of species
sd_rand <-0.01
for(i in 1:n_sp){
  set.seed(i)
  max_new <- runif(1, 1.2, 1.8)
  min_new <- runif(1, 0.2, 0.8)
  y[i,1] <- obs_se[i,1] <- sprintf("SP%03d",i)
  #y[i,2:(n_y+1)] <- scales::rescale(c(arima.sim(model = list(order = c(0, 1, 0)), n = (n_y-1))), to=c(min_new, max_new))
  y_ts <- numeric(n_y)
  y_ts[1] <- rnorm(n = 1, mean = 0, sd = 1)
  for (t in 2:n_y) {
    r.w <- rnorm(n = 1, mean = 0, sd = 1)
    y_ts[t] <- y_ts[t - 1] + r.w
  }
  y[i,2:(n_y+1)] <- scales::rescale(y_ts, to=c(min_new, max_new))
  obs_se[i,2:(n_y+1)] <- runif(n_y,sd_rand,2*sd_rand)
}

y_rand <- data.table(y)
obs_se_rand <- data.table(obs_se)
names(y_rand) <- names(obs_se_rand) <- c("code_sp",1:n_y)
species_rand <- data.frame(name_long=sprintf("species %03d",1:nrow(y_rand)), code_sp=y_rand$code_sp)


# Testing  influence of parameters
# selection of different groups a posteriori
n_y <- 25 # number of year
y <- data.frame(t(rep(NA,(n_y+2))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 200 # number of simulations before selection
sd_rand <-0.01
for(i in 1:n_sp){
  set.seed(i)
  y[i,1] <- obs_se[i,1] <- sprintf("SP%03d",i)
  #y_ts[1] <- rnorm(n = 1, mean = 0, sd = 1)
  #for (t in 2:n_y) {
  #  r.w <- rnorm(n = 1, mean = 0, sd = 1)
  #  y_ts[t] <- y_ts[t - 1] + r.w
  #}
  y_ts <- c(arima.sim(model = list(order = c(0, 1, 0)), n = (n_y-1)))
  y_ts <- y_ts+abs(min(y_ts))+1
  y_ts <- exp(scale(log(y_ts)))
  max_new <- max(y_ts)-mean(y_ts)/4
  min_new <- min(y_ts)+mean(y_ts)/4
  y_ts <- scales::rescale(y_ts, to=c(min_new, max_new))
  y_ts_slope <- summary(lm(y_ts~c(1:n_y)))$coef[2,1]
  y[i,2:(n_y+1)] <- y_ts
  y[i,(n_y+2)] <- y_ts_slope
  obs_se[i,2:(n_y+1)] <- abs(rnorm(n_y,0.1*1/y_ts,sd_rand))
}
n_group <- 5# number of group
y$group <- cut(y[,(n_y+2)], n_group, labels=LETTERS[1:n_group])
n_sp <- 30 # number of species
y_red <- Reduce(rbind, by(y,y["group"],head,n = round(n_sp/n_group)))
obs_se_red <- obs_se[obs_se$X1 %in% y_red$X1,]
y_red <- y_red[order(y_red$X1),]

y_rand <- data.table(y_red[,1:(n_y+1)])
obs_se_rand <- data.table(obs_se_red)
names(y_rand) <- names(obs_se_rand) <- c("code_sp",1:n_y)
y_rand$code_sp <- obs_se_rand$code_sp <- sprintf("SP%03d",1:nrow(y_rand))
species_rand <- data.frame(name_long=sprintf("species %03d",1:nrow(y_rand)), code_sp=y_rand$code_sp)



# DFA

# if number of trends unknown
farm_nfac <- make_dfa2(data_ts = y_farm, data_ts_se = obs_se_farm,
                       species_sub = species_farm)
forest_nfac <- make_dfa2(data_ts = y_forest, data_ts_se = obs_se_forest,
                         species_sub = species_forest)
all_nfac <- make_dfa2(data_ts = y_all, data_ts_se = obs_se_all,
                       species_sub = species_all)
rand_nfac <- make_dfa2(data_ts = y_rand, data_ts_se = obs_se_rand,
                      species_sub = species_rand)

# examples if number of trends known 
farm_nfac3 <- make_dfa2(data_ts = y_farm, data_ts_se = obs_se_farm,
                        nfac = 3, species_sub = species_farm)
forest_nfac4 <- make_dfa2(data_ts = y_forest, data_ts_se = obs_se_forest,
                          nfac = 4, species_sub = species_forest)
all_nfac2 <- make_dfa2(data_ts = y_all, data_ts_se = obs_se_all,
                       nfac = 2, species_sub = species_all)
