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
       width = 8, # 10 # 6
       height = 8 # 8 # 3
       )


# test simulation function

rand_nfac_test <- simul_rand_dfa(nb_group_exp = 2, thres = 1)

# full simulations
rep_sim <- 200

sim_data <- data.frame(num_sim=1:(rep_sim*4),nb_group_exp=sort(rep(c(1,2,3,4),rep_sim)))

library(doParallel)
library(parallel)
no_cores <- detectCores() - 1
registerDoParallel(cores=no_cores)
rand_nfac_test_g <- dlply(sim_data,.(num_sim),
                              .fun=function(x){simul_rand_dfa(nb_group_exp = x$nb_group_exp,
                                                             rep_rand_seed = 2,
                                                             seed= x$num_sim)},
                          .parallel = T)
saveRDS(rand_nfac_test_g,"output/rand_nfac_test_g_new.rds")

sim_data <- data.frame(num_sim=1:(rep_sim*4),nb_group_exp=sort(rep(c(1,2,3,4),rep_sim)), equal=FALSE)
registerDoParallel(cores=no_cores)
rand_nfac_test_g2 <- dlply(sim_data,.(num_sim),
                          .fun=function(x){simul_rand_dfa(nb_group_exp = x$nb_group_exp,
                                                          rep_rand_seed = 2,
                                                          seed= x$num_sim,
                                                          equi = FALSE)},
                          .parallel = T)
saveRDS(rand_nfac_test_g2,"output/rand_nfac_test_g2_new.rds")

sim_data <- data.frame(num_sim=1:(rep_sim*5),sd_ci=sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)))
registerDoParallel(cores=no_cores)
rand_nfac_test_l <- dlply(sim_data,.(num_sim),
                           .fun=function(x){simul_rand_dfa(sd_ci = x$sd_ci,
                                                           rep_rand_seed = 2,
                                                           seed= x$num_sim)},
                           .parallel = T)
saveRDS(rand_nfac_test_l,"output/rand_nfac_test_l_new.rds")

sim_data <- data.frame(num_sim=1:(rep_sim*5),sd_ci=sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)), equal=FALSE)
registerDoParallel(cores=no_cores)
rand_nfac_test_l2 <- dlply(sim_data,.(num_sim),
                          .fun=function(x){simul_rand_dfa(sd_ci = x$sd_ci,
                                                          rep_rand_seed = 2,
                                                          seed= x$num_sim,
                                                          equi = FALSE)},
                          .parallel = T)

saveRDS(rand_nfac_test_l2,"output/rand_nfac_test_l2_new.rds")




# result simulation for number of groups
res_tot_g <- matrix(NA,nrow=length(sort(rep(c(1,2,3,4),rep_sim))),ncol=5)
for(i in 1:length(sort(rep(c(1,2,3,4),rep_sim)))){
  exp_group <- sort(rep(c(1,2,3,4),rep_sim))[i]
  res_i <- rand_nfac_test_g[[i]]
  res_i$clust_stab <- apply(str_split_fixed(res_i$clust_stab, '-', max(as.numeric(res_i$nb_group))),1,function(y){mean(as.numeric(y), na.rm=T)})
  res_i$nb_group <- as.numeric(res_i$nb_group)
  res_i$nb_group2 <- as.numeric(res_i$nb_group2)
  res_i2 <- as.matrix(res_i[,2:5])
  res_tot_g[i,1:ncol(res_i2)] <- res_i2
  res_tot_g[i,5] <- exp_group
}

res_tot_g <- data.frame(res_tot_g)
names(res_tot_g) <- c("similarity","nb_group_obs","stability","nb_group_obs2","nb_group_exp")

table(res_tot_g$nb_group_obs,res_tot_g$nb_group_exp)
table(res_tot_g$nb_group_obs2,res_tot_g$nb_group_exp)

res_tot_g_sum <- data.frame(res_tot_g %>% group_by(nb_group_exp) %>% summarize(mean_similarity = mean(similarity, na.rm=T),
                                                                                 sd_similarity = sd(similarity, na.rm=T),
                                                                                 mean_stability = mean(stability, na.rm=T),
                                                                                 sd_stability = sd(stability, na.rm=T)))

res_tot_g2 <- matrix(NA,nrow=length(sort(rep(c(1,2,3,4),rep_sim))),ncol=5)
for(i in 1:length(sort(rep(c(1,2,3,4),rep_sim)))){
  exp_group <- sort(rep(c(1,2,3,4),rep_sim))[i]
  res_i <- rand_nfac_test_g2[[i]]
  res_i$clust_stab <- apply(str_split_fixed(res_i$clust_stab, '-', max(as.numeric(res_i$nb_group))),1,function(y){mean(as.numeric(y), na.rm=T)})
  res_i$nb_group <- as.numeric(res_i$nb_group)
  res_i2 <- as.matrix(res_i[,2:5])
  res_tot_g2[i,1:ncol(res_i2)] <- res_i2
  res_tot_g2[i,5] <- exp_group
}

res_tot_g2 <- data.frame(res_tot_g2)
names(res_tot_g2) <- c("similarity","nb_group_obs","stability","nb_group_obs2","nb_group_exp")

table(res_tot_g2$nb_group_obs,res_tot_g2$nb_group_exp)
table(res_tot_g2$nb_group_obs2,res_tot_g2$nb_group_exp)


res_tot_g2_sum <- data.frame(res_tot_g2 %>% group_by(nb_group_exp) %>% summarize(mean_similarity = mean(similarity, na.rm=T),
                                                                               sd_similarity = sd(similarity, na.rm=T),
                                                                               mean_stability = mean(stability, na.rm=T),
                                                                               sd_stability = sd(stability, na.rm=T)))


res_tot_l <- matrix(NA,nrow=length(sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim))),ncol=5)
for(i in 1:length(sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)))){
  exp_group <- sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim))[i]
  res_i <- rand_nfac_test_l[[i]]
  res_i$clust_stab <- apply(str_split_fixed(res_i$clust_stab, '-', max(as.numeric(res_i$nb_group))),1,function(y){mean(as.numeric(y), na.rm=T)})
  res_i$nb_group <- as.numeric(res_i$nb_group)
  res_i2 <- as.matrix(res_i[,2:5])
  res_tot_l[i,1:ncol(res_i2)] <- res_i2
  res_tot_l[i,5] <- exp_group
}

res_tot_l <- data.frame(res_tot_l)
names(res_tot_l) <- c("similarity","nb_group_obs","stability","nb_group_obs2","proximity")

table(res_tot_l$nb_group_obs,res_tot_l$proximity)
table(res_tot_l$nb_group_obs2,res_tot_l$proximity)


res_tot_l_sum <- data.frame(res_tot_l %>% group_by(proximity) %>% summarize(mean_similarity = mean(similarity, na.rm=T),
                                                                               sd_similarity = sd(similarity, na.rm=T),
                                                                               mean_stability = mean(stability, na.rm=T),
                                                                               sd_stability = sd(stability, na.rm=T)))

res_tot_l2 <- matrix(NA,nrow=length(sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim))),ncol=5)
for(i in 1:length(sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)))){
  exp_group <- sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim))[i]
  res_i <- rand_nfac_test_l2[[i]]
  res_i$clust_stab <- apply(str_split_fixed(res_i$clust_stab, '-', max(as.numeric(res_i$nb_group))),1,function(y){mean(as.numeric(y), na.rm=T)})
  res_i$nb_group <- as.numeric(res_i$nb_group)
  res_i2 <- as.matrix(res_i[,2:5])
  res_tot_l2[i,1:ncol(res_i2)] <- res_i2
  res_tot_l2[i,5] <- exp_group
}

res_tot_l2 <- data.frame(res_tot_l2)
names(res_tot_l2) <- c("similarity","nb_group_obs","stability","nb_group_obs2","proximity")

table(res_tot_l2$nb_group_obs,res_tot_l2$proximity)
table(res_tot_l2$nb_group_obs2,res_tot_l2$proximity)


res_tot_l2_sum <- data.frame(res_tot_l2 %>% group_by(proximity) %>% summarize(mean_similarity = mean(similarity, na.rm=T),
                                                                            sd_similarity = sd(similarity, na.rm=T),
                                                                            mean_stability = mean(stability, na.rm=T),
                                                                            sd_stability = sd(stability, na.rm=T)))

