## Dynamical Factor Analysis with TMB and uncertainty

# Load packages

source("package_used.R")

# Load functions
source("function_ts.R")

# DFA with TMB
source('function_dfa_clean.R')

# Time-series from the SBBS

ts_bird_se_allcountry <- readRDS("output/ts_bird_se_allcountry.rds")
species_data <- readRDS("output/species_data.rds")
species_data_en_se <- read.csv("output/species_data_en_se.csv", header = T)

# Clean data 

ts_bird_se_allcountry_data <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance) & ts_bird_se_allcountry$CI_inf!=0),]


# Farmland birds

species_sub <- species_farm <-  droplevels(species_data[species_data$code_sp %in% c(
  "FALTIN","VANVAN","ALAARV","HIRRUS","CORFRU",
  "SAXRUB","SYLCOM","ANTPRA","MOTFLA","LANCOL",
  "STUVUL","LINCAN","EMBCIT","EMBHOR","PASMON"),])
  #"VANVAN","NUMARQ","ALAARV","HIRRUS","MOTFLA","OENOEN","SAXRUB","SYLCOM",
  #"LANCOL","STUVUL","LINCAN","EMBCIT","PASMON","CORFRU","ANTPRA","EMBHOR"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]
y_farm <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se_farm <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")


# Forest bird

species_sub <- species_forest <- droplevels(species_data[species_data$code_sp %in% c(
  "ACCNIS","TETBON","TRIOCH","COLOEN","DENMAJ","DRYMAR","PICVIR","JYNTOR",
  "DRYMIN","PICTRI","NUCCAR","GARGLA","PERATE","LOPCRI","POEPAL","POEMON",
  "SITEUR","CERFAM","TURVIS","PHOPHO","PHYCOL","PHYSIB","REGREG","FICHYP",
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
                                                                                  "VANVAN","FALTIN","ALAARV","HIRRUS", # farmland
                                                                                  "MOTFLA","SAXRUB","SYLCOM",
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
y <- data.frame(t(rep(NA,(n_y+3))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 1000 # number of simulations before selection
sd_rand <- 0.01

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
  classification <- class.trajectory(c(y_ts), c(1:n_y))
  if(i==1){
    y_class <- classification$shape_class
  }else{
    y_class <- c(y_class,classification$shape_class)
  }
  y[i,2:(n_y+1)] <- y_ts
  y[i,(n_y+2)] <- classification$linear_slope
  y[i,(n_y+3)] <- classification$linear_slope_pvalue
  obs_se[i,2:(n_y+1)] <- abs(rnorm(n_y,0.1*1/y_ts,sd_rand))
}
#n_group <- 5# number of group
#y$group <- cut(y[,(n_y+2)], n_group, labels=LETTERS[1:n_group])
n_sp <- 30 # number of species
#y <- Reduce(rbind, by(y,y["group"],head,n = round(n_sp/n_group)))
set.seed(1)
y_num_red <- c(sample(which(y_class=="decrease_constant"),10),
            sample(which(y_class=="increase_constant"),10),
            sample(which(y_class=="stable_constant"),10))
y_num_red <- c(sample(which(y_class=="decrease_constant"),3),
               sample(which(y_class=="increase_constant"),3),
               sample(which(y_class=="stable_constant"),3),
               sample(which(y_class=="decrease_accelerated"),3),
               sample(which(y_class=="increase_accelerated"),3),
               sample(which(y_class=="stable_concave"),3),
               sample(which(y_class=="decrease_decelerated"),3),
               sample(which(y_class=="increase_decelerated"),3),
               sample(which(y_class=="stable_convex"),3))
y <- data.frame(y[y_num_red,], class=y_class[y_num_red])
obs_se_red <- obs_se[obs_se$X1 %in% y$X1,]
y <- y[order(y$X1),]

y_rand <- data.table(y[,1:(n_y+1)])
obs_se_rand <- data.table(obs_se_red)
names(y_rand) <- names(obs_se_rand) <- c("code_sp",1:n_y)
y_rand$code_sp <- obs_se_rand$code_sp <- sprintf("SP%03d",1:nrow(y_rand))
species_rand <- data.frame(name_long=sprintf("species %03d",1:nrow(y_rand)), code_sp=y_rand$code_sp)

# swedish indicator birds

species_sub <- species_forest_sw <- droplevels(species_data[species_data$code_sp %in% c("TETURO","TETBON","COLOEN","PICVIR","DRYMIN","PICTRI","NUCCAR","PERINF",
                                                                                        "AEGCAU","PERATE","LOPCRI","POECIN","POEPAL","POEMON","CERFAM","PYRPYR"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

y_forest_sw <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
                     code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se_forest_sw <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                          code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

# mountrain birds
species_sub <- species_rock <- droplevels(species_data[species_data$code_sp %in% c("LAGMUT","PLUAPR","STELON","OENOEN",
                                                                                   "ANTPRA","CALLAP","PLENIV","LAGLAG",
                                                                                   "TURILI","PHOPHO","LUSSVE","PHYTRO",
                                                                                   "ACAFLA","FRIMON"),])

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

y_rock <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
                code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
obs_se_rock <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                     code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")


# DFA

# if number of trends unknown
farm_nfac <- make_dfa(data_ts = y_farm, data_ts_se = obs_se_farm,
                       species_sub = species_farm)
forest_nfac <- make_dfa(data_ts = y_forest, data_ts_se = obs_se_forest,
                         species_sub = species_forest)
all_nfac <- make_dfa(data_ts = y_all, data_ts_se = obs_se_all,
                       species_sub = species_all)
rand_nfac <- make_dfa(data_ts = y_rand, data_ts_se = obs_se_rand,
                      species_sub = species_rand)
forest_sw_nfac <- make_dfa(data_ts = y_forest_sw, data_ts_se = obs_se_forest_sw,
                            species_sub = species_forest_sw)
rock_nfac <- make_dfa(data_ts = y_rock, data_ts_se = obs_se_rock,
                       species_sub = species_rock)


ggsave("output/farm_nfac.png",
       dpi=300,
       width = 8,
       height = 8
       )

# examples if number of trends known 
farm_nfac3 <- make_dfa2(data_ts = y_farm, data_ts_se = obs_se_farm,
                        nfac = 3, species_sub = species_farm)
forest_nfac4 <- make_dfa2(data_ts = y_forest, data_ts_se = obs_se_forest,
                          nfac = 4, species_sub = species_forest)
all_nfac2 <- make_dfa2(data_ts = y_all, data_ts_se = obs_se_all,
                       nfac = 2, species_sub = species_all)

# random simulation then classification

cum_perc <- expand.grid(c(0,20,40,60,80,100),
                        c(0,20,40,60,80,100),
                        c(0,20,40,60,80,100))
cum_perc[,4] <- apply(cum_perc,1,sum)
cum_perc <- cum_perc[cum_perc$V4==100,1:3]

n_y <- 25 # number of year
y <- data.frame(t(rep(NA,(n_y+3))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 1000 # number of simulations before selection
sd_rand <- 0.01

for(i in 1:n_sp){
  set.seed(i)
  y[i,1] <- obs_se[i,1] <- sprintf("SP%03d",i)
  y_ts <- c(arima.sim(model = list(order = c(0, 1, 0)), n = (n_y-1)))
  y_ts <- y_ts+abs(min(y_ts))+1
  y_ts <- exp(scale(log(y_ts)))
  max_new <- max(y_ts)-mean(y_ts)/4
  min_new <- min(y_ts)+mean(y_ts)/4
  y_ts <- scales::rescale(y_ts, to=c(min_new, max_new))
  classification <- class.trajectory(c(y_ts), c(1:n_y))
  if(i==1){
    y_class <- classification$shape_class
  }else{
    y_class <- c(y_class,classification$shape_class)
  }
  y[i,2:(n_y+1)] <- y_ts
  y[i,(n_y+2)] <- classification$linear_slope
  y[i,(n_y+3)] <- classification$linear_slope_pvalue
  obs_se[i,2:(n_y+1)] <- abs(rnorm(n_y,0.1*1/y_ts,sd_rand))
}

n_sp <- 30

rand_nfac_list <- list()
for(a in 1:nrow(cum_perc)){
  ud <- round(cum_perc[a,1]*n_sp/100)
  ui <- round(cum_perc[a,2]*n_sp/100)
  uc <- round(cum_perc[a,3]*n_sp/100)
  
  set.seed(a)
  
  y_num_red <- c(sample(which(y_class=="decrease_constant"),ud),
                 sample(which(y_class=="increase_constant"),ui),
                 sample(which(y_class=="stable_constant"),uc))
  y <- data.frame(y[y_num_red,], class=y_class[y_num_red])
  obs_se_red <- obs_se[obs_se$X1 %in% y$X1,]
  y <- y[order(y$X1),]
  
  y_rand <- data.table(y[,1:(n_y+1)])
  obs_se_rand <- data.table(obs_se_red)
  names(y_rand) <- names(obs_se_rand) <- c("code_sp",1:n_y)
  y_rand$code_sp <- obs_se_rand$code_sp <- sprintf("SP%03d",1:nrow(y_rand))
  species_rand <- data.frame(name_long=sprintf("species %03d",1:nrow(y_rand)), code_sp=y_rand$code_sp)
  
  # DFA
  
  rand_nfac <- make_dfa2(data_ts = y_rand, data_ts_se = obs_se_rand,
                         species_sub = species_rand)
  
  rand_nfac_list[[a]] <- rand_nfac
}


# do it many time (but output to big)

n_y <- 25 # number of year
y <- data.frame(t(rep(NA,(n_y+3))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))
n_sp <- 10000 # number of simulations before selection
sd_rand <- 0.01

for(i in 1:n_sp){
  set.seed(i)
  y[i,1] <- obs_se[i,1] <- sprintf("SP%03d",i)
  y_ts <- c(arima.sim(model = list(order = c(0, 1, 0)), n = (n_y-1)))
  y_ts <- y_ts+abs(min(y_ts))+1
  y_ts <- exp(scale(log(y_ts)))
  max_new <- max(y_ts)-mean(y_ts)/4
  min_new <- min(y_ts)+mean(y_ts)/4
  y_ts <- scales::rescale(y_ts, to=c(min_new, max_new))
  classification <- class.trajectory(c(y_ts), c(1:n_y))
  if(i==1){
    y_class <- classification$shape_class
  }else{
    y_class <- c(y_class,classification$shape_class)
  }
  y[i,2:(n_y+1)] <- y_ts
  y[i,(n_y+2)] <- classification$linear_slope
  y[i,(n_y+3)] <- classification$linear_slope_pvalue
  obs_se[i,2:(n_y+1)] <- abs(rnorm(n_y,0.1*1/y_ts,sd_rand))
}

n_sp <- 30
n_simul <- 100
rand_nfac_list_all <- list()
for(g in 1:n_simul){
  print(g)
  rand_nfac_list <- list()
  for(a in 1:nrow(cum_perc)){
    ud <- round(cum_perc[a,1]*n_sp/100)
    ui <- round(cum_perc[a,2]*n_sp/100)
    uc <- round(cum_perc[a,3]*n_sp/100)
    
    #set.seed(a)
    
    y_num_red <- c(sample(which(y_class=="decrease_constant"),ud),
                   sample(which(y_class=="increase_constant"),ui),
                   sample(which(y_class=="stable_constant"),uc))
    y <- data.frame(y[y_num_red,], class=y_class[y_num_red])
    obs_se_red <- obs_se[obs_se$X1 %in% y$X1,]
    y <- y[order(y$X1),]
    
    y_rand <- data.table(y[,1:(n_y+1)])
    obs_se_rand <- data.table(obs_se_red)
    names(y_rand) <- names(obs_se_rand) <- c("code_sp",1:n_y)
    y_rand$code_sp <- obs_se_rand$code_sp <- sprintf("SP%03d",1:nrow(y_rand))
    species_rand <- data.frame(name_long=sprintf("species %03d",1:nrow(y_rand)), code_sp=y_rand$code_sp)
    
    # DFA
    
    rand_nfac <- make_dfa2(data_ts = y_rand, data_ts_se = obs_se_rand,
                           species_sub = species_rand)
    
    rand_nfac_list[[a]] <- rand_nfac
  }
  rand_nfac_list_all[[g]] <- rand_nfac_list
}

# random simulation with classification a priori

simul_rand_dfa_intern <- function(a,
                                  cum_perc,
                                  n_sp_init,
                                  nb_group_exp,
                                  thres,
                                  n_y,
                                  n_sp,
                                  y_init,
                                  sd_rand,
                                  sd_rand2,
                                  sd_ci,
                                  nboot){
  ## Simulate latent trends
  
  y_init <- data.frame(t(rep(NA,(n_y)))) # latent trends
  test_cor <- 1
  while(abs(test_cor)>0.5){
    for(i in 1:n_sp_init){
      #set.seed(i+10)
      y_ts <- c()
      y_ts[1] <- rnorm(n = 1, mean = 0, sd = 1)
      for (t in 2:n_y) {
        r.w <- rnorm(n = 1, mean = 0, sd = 1)
        y_ts[t] <- y_ts[t - 1] + r.w
      }
      y_ts <- y_ts + abs(min(y_ts))+1
      y_ts <- exp(scale(log(y_ts)))
      #max_new <- max(y_ts) - mean(y_ts)/4
      #min_new <- min(y_ts) + mean(y_ts)/4
      #y_ts <- scales::rescale(y_ts, to=c(min_new, max_new))
      y_init[i,] <- y_ts
    }
    test_cor <- cor.test(as.numeric(y_init[1,]),as.numeric(y_init[2,]),method = "spearman")$estimate
  }
  
  ## from these n_sp_init latent trend, simulate n_sp ts from loading factors
  
    seed_id <- 0
    id_vec <- c()
    dist_bary <- 0
    mat_dist <- matrix(NA, ncol=n_sp_init, nrow=nb_group_exp)
    while(dist_bary<thres){ # check if enough distance between the groups
      for(g in 1:nb_group_exp){
        nb_sp_g <- round(cum_perc[a,g]*n_sp/100)
        assign(paste0("nb_sp_g",g),nb_sp_g)
        id_vec <- c(id_vec,rep(g,nb_sp_g))
        for(lt in 1:n_sp_init){
          seed_id <- seed_id + 1
          #set.seed(seed_id)
          mean_u_g <- runif(1, -1, 1)
          lf_u_g <- rnorm(nb_sp_g, mean_u_g, sd_ci)
          assign(paste0("mean_u",lt,"_g",g),mean_u_g) # mean of loading factors in group g for latend trend lt
          assign(paste0("lf_u",lt,"_g",g),lf_u_g) # loading factors for each ts of group g for latend trend lt
          mat_dist[g,lt] <- mean_u_g
        }
      }
      id_vec <- id_vec[1:n_sp]
      dist_bary <- min(dist(mat_dist))
    }
    
    y <- data.frame(t(rep(NA,(n_y+2))))
    obs_se <- data.frame(t(rep(NA,(n_y+1))))
    
    for(i in 1:n_sp){ # get simulated ts from loadings
      #set.seed(i)
      noise <- rnorm(n_y,0,sd_rand2)
      y[i,1] <- obs_se[i,1] <- sprintf("SP%03d",i)
      y_ts <- rep(0,n_y)
      g <- id_vec[i]
      i_g <- which(which(id_vec==g)==i) # new index for i in group g
      for(lt in 1:n_sp_init){
        lf_u_g <- get(paste0("lf_u",lt,"_g",g))
        y_ts <- y_ts + as.numeric(y_init[lt,])*lf_u_g[i_g]
      }
      y_ts <- y_ts + noise
      y_ts <- y_ts + abs(min(y_ts)) + 1
      y_ts <- exp(scale(log(y_ts)))
      y[i,2:(n_y+1)] <- y_ts
      y[i,(n_y+2)] <- id_vec[i]
      obs_se[i,2:(n_y+1)] <- abs(rnorm(n_y,0.1*1/y_ts,sd_rand))
      obs_se[obs_se>1] <- 1
    }  
    
    
    y_rand <- data.table(y[,1:(n_y+1)])
    obs_se_rand <- data.table(obs_se)
    names(y_rand) <- names(obs_se_rand) <- c("code_sp",1:n_y)
    y_rand$code_sp <- obs_se_rand$code_sp <- sprintf("SP%03d",1:nrow(y_rand))
    species_rand <- data.frame(name_long=sprintf("species %03d",1:nrow(y_rand)), code_sp=y_rand$code_sp)
    
    # DFA
    
    rand_nfac <- make_dfa(data_ts = y_rand, data_ts_se = obs_se_rand,
                           species_sub = species_rand, nboot=nboot)
    
    # compare DFA results to expected
    
    obs_group <- rand_nfac[[10]][[1]][[1]]
    y[,ncol(y)] <- as.numeric(as.factor(y[,ncol(y)]))
    
    if(length(obs_group)==1){
      obs_group_new <- rep(1,nrow(y_rand))
      clust_nb <- 1
      clust_stab <- 1
    }else{
      clust_stab <- gsub(", ","-",toString(paste0(round(rand_nfac[[10]][[3]],2))))
      
      clust_nb <- length(unique(obs_group$group))
      jac_sim_res <- matrix(NA, ncol=length(unique(y[,ncol(y)])),
                            nrow=length(unique(obs_group$group)))
      for(k in sort(unique(y[,ncol(y)]))){
        for(l in sort(unique(obs_group$group))){
          jac_sim_mat <- rbind(y[,ncol(y)],obs_group$group)
          jac_sim_mat[1,][which(jac_sim_mat[1,]!=k)] <- 0
          jac_sim_mat[2,][which(jac_sim_mat[2,]!=l)] <- 0
          jac_sim_mat[jac_sim_mat>0] <- 1
          jac_sim <- c(1 - vegdist(jac_sim_mat, method="jaccard"))
          jac_sim_res[l,k] <- jac_sim
        }
      }
      
      obs_group_new <- rep(NA,length(obs_group$group))
      
      # If same number of clusters
      
      if(length(unique(y[,ncol(y)]))==length(unique(obs_group$group))){
        for(l in sort(unique(obs_group$group))){
          obs_group_new[which(obs_group$group==l)] <- which.max(jac_sim_res[l,])
        }
      }
      
      # If more clusters in the observed clustering
      
      if(length(unique(y[,ncol(y)]))<length(unique(obs_group$group))){
        l_data <- c()
        for(k in sort(unique(y[,ncol(y)]))){
          l_data <- c(l_data,which.max(jac_sim_res[,k]))
        }
        k <- 0
        for(l in l_data){
          k <- k+1
          obs_group_new[which(obs_group$group==l)] <- k
        }
        extra_clus <- sort(unique(obs_group$group))[which(!(sort(unique(obs_group$group)) %in% l_data))]
        for(g_sup in 1:length(extra_clus)){
          k <- k +1
          obs_group_new[which(obs_group$group==extra_clus[g_sup])] <- k
        }
      }
      
      # If less clusters in the bootstrap clustering
      
      if(length(unique(y[,ncol(y)]))>length(unique(obs_group$group))){
        k_data <- c()
        for(l in sort(unique(obs_group$group))){
          k_data <- c(k_data,which.max(jac_sim_res[l,]))
        }
        l <- 0
        for(k in k_data){
          l <- l+1
          obs_group_new[which(obs_group$group==l)] <- k
        }
      }
    }
    
    res_rand <- c(1 - vegdist(rbind(y[,ncol(y)],obs_group_new), method="jaccard"))
    return(list(res_rand, clust_nb, clust_stab))
}

simul_rand_dfa_intern2 <- function(a,cum_perc,n_sp_init,
                     nb_group_exp,thres,n_y,
                     n_sp,y_init,sd_rand,
                     sd_rand2,sd_ci,nboot){
  tryCatch(simul_rand_dfa_intern(a,cum_perc,n_sp_init,
                                 nb_group_exp,thres,n_y,
                                 n_sp,y_init,sd_rand,
                                 sd_rand2,sd_ci,nboot),
           error=function(e) list(NA,NA,NA))}

simul_rand_dfa <- function(n_y = 25, # number of year
                           n_sp = 30, # number of species ts
                           n_sp_init = 2, # number of latent trends
                           nb_group_exp = 2, # number of expected clusters
                           thres = 1, # min distance between barycenters of clusters
                           sd_rand = 0.01, # observation error on data
                           sd_rand2 = 0.5, # random noise on ts
                           sd_ci = 0.1, # standard deviation of the loading factors
                           nboot = 100 # number of bootstrap for clustering
                           ){
  
  list_to_expend <- list()
  for(g in 1:nb_group_exp){
    list_to_expend[[g]] <- c(0,20,40,60,80,100)
  }
  cum_perc <- expand.grid(list_to_expend)
  cum_perc[,(nb_group_exp+1)] <- apply(cum_perc,1,sum)
  cum_perc <- cum_perc[which(cum_perc[,ncol(cum_perc)]==100),1:nb_group_exp]
  cum_perc$nb_sup_0 <- apply(cum_perc, 1, function(x){sum(sign(x))})
  
  if(nb_group_exp>2){
    cum_perc <- cum_perc[which(cum_perc[,ncol(cum_perc)]==nb_group_exp),1:nb_group_exp]
  }else{
    cum_perc <- cum_perc[,1:nb_group_exp]
  }
  
  
  library(parallel)
  
  # Calculate the number of cores
  no_cores <- detectCores() - 1
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type="FORK")
  
  rand_nfac_list <- parSapply(cl, c(1:nrow(cum_perc)),
                               FUN=function(x){unlist(simul_rand_dfa_intern2(a=x,
                                                                     cum_perc,
                                                                     n_sp_init,
                                                                     nb_group_exp,
                                                                     thres,
                                                                     n_y,
                                                                     n_sp,
                                                                     y_init,
                                                                     sd_rand,
                                                                     sd_rand2,
                                                                     sd_ci,
                                                                     nboot))})
  stopCluster(cl)
  
  res_sim <- data.frame(perc_group=NA, value=as.numeric(rand_nfac_list[1,]),
                        nb_group=rand_nfac_list[2,], clust_stab=rand_nfac_list[3,])
  for(a in 1:nrow(cum_perc)){
    res_sim[a,1] <- gsub(", ","-",toString(paste0(cum_perc[a,])))
  }
  
  return(res_sim)
}

# test function before simulation
rand_nfac_test <- simul_rand_dfa(nb_group_exp = 2, thres = 1)
rand_nfac_test2 <- lapply(c(1,2,3,4),
                          FUN=function(x){simul_rand_dfa(nb_group_exp = x,
                                                         thres = 1,
                                                         sd_ci = 0.01)})
rand_nfac_test3 <- lapply(c(0.25,0.5,0.75),
                          FUN=function(x){simul_rand_dfa(nb_group_exp = 2,
                                                                thres = x)})
rand_nfac_test4 <- lapply(c(0.01,1),
                          FUN=function(x){simul_rand_dfa(nb_group_exp = 2,
                                                         thres = 1,
                                                         sd_ci = x)})


# full simulations
rep_sim <- 100

sim_data <- data.frame(num_sim=1:(rep_sim*4),nb_group_exp=sort(rep(c(1,2,3,4),rep_sim)))

library(doParallel)
library(parallel)
no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
rand_nfac_test_g <- dlply(sim_data,.(num_sim),
                              .fun=function(x){simul_rand_dfa(nb_group_exp = x$nb_group_exp,
                                                             thres = 1)},
                          .parallel = T)



library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
# Initiate cluster
cl <- makeCluster(no_cores, type="FORK")
rand_nfac_test_g <- parLapply(cl,sort(rep(c(1,2,3,4),rep_sim)),
                          fun=function(x){simul_rand_dfa(nb_group_exp = x,
                                                                thres = 1)})
stopCluster(cl)


no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
rand_nfac_test_d <- parLapply(cl,sort(rep(c(0.25,0.5,0.75,1.5,2),rep_sim)),
                           fun=function(x){simul_rand_dfa(nb_group_exp = 2,
                                                                 thres = x)})
stopCluster(cl)

no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")
rand_nfac_test_l <- parLapply(cl,sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)),
                           fun=function(x){simul_rand_dfa(nb_group_exp = 2,
                                                          thres = 1,
                                                          sd_ci = x)})

stopCluster(cl)



# result simulation for number of groups
res_tot_g <- array(NA,dim=c(6,4,length(sort(rep(c(2,3,4),rep_sim)))))
for(i in 1:length(sort(rep(c(2,3,4),rep_sim)))){
  exp_group <- sort(rep(c(2,3,4),rep_sim))[i]
  
  res_i <- rand_nfac_test_g[[i]]
  res_i$clust_stab <- apply(str_split_fixed(res_i$clust_stab, '-', max(as.numeric(res_i$nb_group))),1,function(y){mean(as.numeric(y), na.rm=T)})
  res_i$nb_group <- as.numeric(res_i$nb_group)
  res_i2 <- as.matrix(res_i[,2:4])
  res_tot_g[1:nrow(res_i2),1:ncol(res_i2),i] <- res_i2
  res_tot_g[,4,i] <- exp_group
}

for(j in 2:4){
  res_mean_g <- apply(res_tot_g[1:nrow(rand_nfac_test_g[[((j-1)*rep_sim)]]),,(((j-2)*rep_sim+1):((j-1)*rep_sim))], c(1,2), function(x){mean(x, na.rm=T)})
  res_sd_g <- apply(res_tot_g[1:nrow(rand_nfac_test_g[[((j-1)*rep_sim)]]),,(((j-2)*rep_sim+1):((j-1)*rep_sim))], c(1,2), function(x){sd(x, na.rm=T)})
  res_mean_g <- data.frame(perc_group=rand_nfac_test_g[[((j-1)*rep_sim)]]$perc_group,
                            nb_group_exp=j,
                            simil=res_mean_g[,1],
                            nb_group_obs=res_mean_g[,2],
                            clust_stab=res_mean_g[,3])
  res_sd_g <- data.frame(perc_group=rand_nfac_test_g[[((j-1)*rep_sim)]]$perc_group,
                          nb_group_exp=j,
                          simil=res_sd_g[,1],
                          nb_group_obs=res_sd_g[,2],
                          clust_stab=res_sd_g[,3])
  assign(paste0("res_mean_g",j),res_mean_g)
  assign(paste0("res_sd_g",j),res_sd_g)
}

res_mean_g <- rbind(res_mean_g2,res_mean_g3,res_mean_g4)
res_mean_g <- melt(res_mean_g,id.vars = c("perc_group", "nb_group_exp"))

res_sd_g <- rbind(res_sd_g2,res_sd_g3,res_sd_g4)
res_sd_g <- melt(res_sd_g,id.vars = c("perc_group", "nb_group_exp"))
names(res_sd_g)[4] <- "sd"

# global results
res_g <- merge(res_mean_g, res_sd_g, by=c("perc_group", "nb_group_exp", "variable"))

ggplot(res_g, aes(nb_group_exp, value)) +
  geom_point(aes(col=variable))

# detailed results
res_det_g <- data.frame(perc_group=rand_nfac_test_g[[1]]$perc_group,res_tot_g[,,1])
for(i in 2:dim(res_tot_g)[3]){
  res_det_g <- rbind(res_det_g, data.frame(perc_group=rand_nfac_test_g[[i]]$perc_group,res_tot_g[1:nrow(rand_nfac_test_g[[i]]),,i]))
}
names(res_det_g) <- c("perc_group","simil","nb_group_obs","clust_stab","nb_group_exp")

res_det_g$nb_group_exp[res_det_g$perc_group=="100-0" | res_det_g$perc_group=="0-100"] <- 1

ggplot(res_det_g, aes(nb_group_exp, simil)) +
  geom_boxplot(aes(group=nb_group_exp))
ggplot(res_det_g, aes(nb_group_exp, clust_stab)) +
  geom_boxplot(aes(group=nb_group_exp))
ggplot(res_det_g, aes(nb_group_exp, nb_group_obs)) +
  geom_point(aes(group=nb_group_exp),position="jitter")


# result simulation for distance between group
res_tot_d <- array(NA,dim=c(6,4,length(sort(rep(c(0.25,0.5,0.75,1.5,2),rep_sim)))))
for(i in 1:length(sort(rep(c(0.25,0.5,0.75,1.5,2),rep_sim)))){
  exp_dist_group <- sort(rep(c(0.25,0.5,0.75,1.5,2),rep_sim))[i]
  
  res_i <- rand_nfac_test_d[[i]]
  res_i$clust_stab <- apply(str_split_fixed(res_i$clust_stab, '-', max(as.numeric(res_i$nb_group))),1,function(y){mean(as.numeric(y), na.rm=T)})
  res_i$nb_group <- as.numeric(res_i$nb_group)
  res_i2 <- as.matrix(res_i[,2:4])
  res_tot_d[1:nrow(res_i2),1:ncol(res_i2),i] <- res_i2
  res_tot_d[,4,i] <- exp_dist_group
}

for(j in 1:5){
  seq_dist_group <- c(0.25,0.5,0.75,1.5,2)
  res_mean_d <- apply(res_tot_d[1:nrow(rand_nfac_test_d[[(j*rep_sim)]]),,(((j-1)*rep_sim+1):(j*rep_sim))], c(1,2), function(x){mean(x, na.rm=T)})
  res_sd_d <- apply(res_tot_d[1:nrow(rand_nfac_test_d[[(j*rep_sim)]]),,(((j-1)*rep_sim+1):(j*rep_sim))], c(1,2), function(x){sd(x, na.rm=T)})
  res_mean_d <- data.frame(perc_group=rand_nfac_test_d[[(j*rep_sim)]]$perc_group,
                           dist_group_exp=seq_dist_group[j],
                           simil=res_mean_d[,1],
                           nb_group_obs=res_mean_d[,2],
                           clust_stab=res_mean_d[,3])
  res_sd_d <- data.frame(perc_group=rand_nfac_test_d[[(j*rep_sim)]]$perc_group,
                         dist_group_exp=seq_dist_group[j],
                         simil=res_sd_d[,1],
                         nb_group_obs=res_sd_d[,2],
                         clust_stab=res_sd_d[,3])
  assign(paste0("res_mean_d",j),res_mean_d)
  assign(paste0("res_sd_d",j),res_sd_d)
}

res_mean_d <- rbind(res_mean_d1,res_mean_d2,res_mean_d3,res_mean_d4,res_mean_d5)
res_mean_d <- melt(res_mean_d,id.vars = c("perc_group", "dist_group_exp"))

res_sd_d <- rbind(res_sd_d1,res_sd_d2,res_sd_d3,res_sd_d4,res_sd_d5)
res_sd_d <- melt(res_sd_d,id.vars = c("perc_group", "dist_group_exp"))
names(res_sd_d)[4] <- "sd"

# global results
res_d <- merge(res_mean_d, res_sd_d, by=c("perc_group", "dist_group_exp", "variable"))

ggplot(res_d, aes(dist_group_exp, value)) +
  geom_point(aes(col=variable))

# detailed results
res_det_d <- data.frame(perc_group=rand_nfac_test_d[[1]]$perc_group,res_tot_d[,,1])
for(i in 2:dim(res_tot_d)[3]){
  res_det_d <- rbind(res_det_d, data.frame(perc_group=rand_nfac_test_d[[i]]$perc_group,res_tot_d[1:nrow(rand_nfac_test_d[[i]]),,i]))
}
for(i in 1:rep_sim){
  to_add <- data.frame(perc_group=rand_nfac_test_g[[i]]$perc_group,res_tot_g[1:nrow(rand_nfac_test_g[[i]]),,i])
  to_add[,5] <- 1
  res_det_d <- rbind(res_det_d, to_add)
}
names(res_det_d) <- c("perc_group","simil","nb_group_obs","clust_stab","dist_group")

res_det_d2 <- res_det_d[res_det_d$perc_group!="100-0" & res_det_d$perc_group!="0-100",]

ggplot(res_det_d2, aes(dist_group, simil)) +
  geom_boxplot(aes(group=dist_group))
ggplot(res_det_d2, aes(dist_group, clust_stab)) +
  geom_boxplot(aes(group=dist_group))


# result simulation for standard error of factor loadings
res_tot_l <- array(NA,dim=c(6,4,length(sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)))))
for(i in 1:length(sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)))){
  exp_loadfact_sd <- sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim))[i]
  
  res_i <- rand_nfac_test_l[[i]]
  res_i$clust_stab <- apply(str_split_fixed(res_i$clust_stab, '-', max(as.numeric(res_i$nb_group))),1,function(y){mean(as.numeric(y), na.rm=T)})
  res_i$nb_group <- as.numeric(res_i$nb_group)
  res_i2 <- as.matrix(res_i[,2:4])
  res_tot_l[1:nrow(res_i2),1:ncol(res_i2),i] <- res_i2
  res_tot_l[,4,i] <- exp_loadfact_sd
}

for(j in 1:5){
  seq_loadfact_sd <- c(0.01,0.05,0.2,0.5,1)
  res_mean_l <- apply(res_tot_l[1:nrow(rand_nfac_test_l[[(j*rep_sim)]]),,(((j-1)*rep_sim+1):(j*rep_sim))], c(1,2), function(x){mean(x, na.rm=T)})
  res_sd_l <- apply(res_tot_l[1:nrow(rand_nfac_test_l[[(j*rep_sim)]]),,(((j-1)*rep_sim+1):(j*rep_sim))], c(1,2), function(x){sd(x, na.rm=T)})
  res_mean_l <- data.frame(perc_group=rand_nfac_test_l[[(j*rep_sim)]]$perc_group,
                           loadfact_sd_exp=seq_loadfact_sd[j],
                           simil=res_mean_l[,1],
                           nb_group_obs=res_mean_l[,2],
                           clust_stab=res_mean_l[,3])
  res_sd_l <- data.frame(perc_group=rand_nfac_test_l[[(j*rep_sim)]]$perc_group,
                         loadfact_sd_exp=seq_loadfact_sd[j],
                         simil=res_sd_l[,1],
                         nb_group_obs=res_sd_l[,2],
                         clust_stab=res_sd_l[,3])
  assign(paste0("res_mean_l",j),res_mean_l)
  assign(paste0("res_sd_l",j),res_sd_l)
}

res_mean_l <- rbind(res_mean_l1,res_mean_l2,res_mean_l3,res_mean_l4,res_mean_l5)
res_mean_l <- melt(res_mean_l,id.vars = c("perc_group", "loadfact_sd_exp"))

res_sd_l <- rbind(res_sd_l1,res_sd_l2,res_sd_l3,res_sd_l4,res_sd_l5)
res_sd_l <- melt(res_sd_l,id.vars = c("perc_group", "loadfact_sd_exp"))
names(res_sd_l)[4] <- "sd"

# global results
res_l <- merge(res_mean_l, res_sd_l, by=c("perc_group", "loadfact_sd_exp", "variable"))

ggplot(res_l, aes(loadfact_sd_exp, value)) +
  geom_point(aes(col=variable))

# detailed results
res_det_l <- data.frame(perc_group=rand_nfac_test_l[[1]]$perc_group,res_tot_l[,,1])
for(i in 2:dim(res_tot_l)[3]){
  res_det_l <- rbind(res_det_l, data.frame(perc_group=rand_nfac_test_l[[i]]$perc_group,res_tot_l[1:nrow(rand_nfac_test_l[[i]]),,i]))
}
for(i in 1:rep_sim){
  to_add <- data.frame(perc_group=rand_nfac_test_g[[i]]$perc_group,res_tot_g[1:nrow(rand_nfac_test_g[[i]]),,i])
  to_add[,5] <- 0.1
  res_det_l <- rbind(res_det_l, to_add)
}
names(res_det_l) <- c("perc_group","simil","nb_group_obs","clust_stab","loadfact_sd")

res_det_l2 <- res_det_l[res_det_l$perc_group!="100-0" & res_det_l$perc_group!="0-100",]

ggplot(res_det_l2, aes(loadfact_sd, simil)) +
  geom_boxplot(aes(group=loadfact_sd))
ggplot(res_det_l2, aes(loadfact_sd, clust_stab)) +
  geom_boxplot(aes(group=loadfact_sd))
