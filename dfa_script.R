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


# DFA

# if number of trends unknown
farm_nfac <- make_dfa(data_ts = y_farm, data_ts_se = obs_se_farm,
                       species_sub = species_farm)
forest_nfac <- make_dfa(data_ts = y_forest, data_ts_se = obs_se_forest,
                         species_sub = species_forest)
all_nfac <- make_dfa(data_ts = y_all, data_ts_se = obs_se_all,
                       species_sub = species_all)

ggsave("output/farm_nfac.png",
       dpi=300,
       width = 8,
       height = 8
       )


# test simulation function

rand_nfac_test <- simul_rand_dfa(nb_group_exp = 2, thres = 1)

# full simulations
rep_sim <- 100

sim_data <- data.frame(num_sim=1:(rep_sim*5),nb_group_exp=sort(rep(c(1,2,3,4,5),rep_sim)))

library(doParallel)
library(parallel)
no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
rand_nfac_test_g <- dlply(sim_data,.(num_sim),
                              .fun=function(x){simul_rand_dfa(nb_group_exp = x$nb_group_exp,
                                                             thres = 1,
                                                             seed= x$num_sim)},
                          .parallel = T)

sim_data <- data.frame(num_sim=1:(rep_sim*5),nb_group_exp=sort(rep(c(1,2,3,4,5),rep_sim)), equal=FALSE)
registerDoParallel(cores=no_cores)
rand_nfac_test_g2 <- dlply(sim_data,.(num_sim),
                          .fun=function(x){simul_rand_dfa(nb_group_exp = x$nb_group_exp,
                                                          thres = 1,
                                                          seed= x$num_sim)},
                          .parallel = T)


sim_data <- data.frame(num_sim=1:(rep_sim*5),sd_ci=sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)))
registerDoParallel(cores=no_cores)
rand_nfac_test_l <- dlply(sim_data,.(num_sim),
                           .fun=function(x){simul_rand_dfa(sd_ci = x$sd_ci,
                                                           seed= x$num_sim)},
                           .parallel = T)

sim_data <- data.frame(num_sim=1:(rep_sim*5),sd_ci=sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)), equal=FALSE)
registerDoParallel(cores=no_cores)
rand_nfac_test_l2 <- dlply(sim_data,.(num_sim),
                          .fun=function(x){simul_rand_dfa(sd_ci = x$sd_ci,
                                                          seed= x$num_sim)},
                          .parallel = T)

saveRDS(rand_nfac_test_l2,"output/rand_nfac_test_l2_new.rds")




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
