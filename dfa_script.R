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



# mountain birds (with 0)

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
farm_nfac <- make_dfa(data_ts = y_farm, data_ts_se = obs_se_farm, #nfac=3,
                       species_sub = species_farm)
forest_nfac <- make_dfa(data_ts = y_forest, data_ts_se = obs_se_forest, # nfac=4,
                         species_sub = species_forest)
all_nfac <- make_dfa(data_ts = y_all, data_ts_se = obs_se_all,
                       species_sub = species_all)

ggsave("output/farm_nfac.png",
       dpi=300,
       width = 8, # 10 # 6
       height = 8 # 8 # 3
       )

table_res_farma <- dcast(data=farm_nfac$data_loadings, name_long~variable, id.vars="value")
table_res_farmb <- farm_nfac$group[[1]][[1]]
table_res_farm <- merge(table_res_farmb[,c("name_long", "group", "uncert")], table_res_farma, by="name_long", all.x=T )
table_res_farm$pvalue <- table_res_farm$st_error <- table_res_farm$slope <-  NA
for(i in 1:nrow(table_res_farm)){
  res_lm <- summary(lm(pred.value~Year,
                       data=farm_nfac$data_to_plot_sp[farm_nfac$data_to_plot_sp$code_sp==unique(farm_nfac$data_to_plot_sp$code_sp)[i],]))$coef
  table_res_farm$slope[i] <- res_lm[2,1]
  table_res_farm$st_error[i] <- res_lm[2,2]
  table_res_farm$pvalue[i] <- res_lm[2,4]
}
table_res_farmc <- farm_nfac$exp_var_lt[,1:5]
names(table_res_farmc) <- c("name_long","% Latent trend 1","% Latent trend 2",
"% Latent trend 3","% eta")
table_res_farm <- merge(table_res_farm, table_res_farmc, by="name_long", all.x=T )
table_res_farm[,3:ncol(table_res_farm)] <- round(table_res_farm[,3:ncol(table_res_farm)],4)
table_res_farm[,10:ncol(table_res_farm)] <- abs(table_res_farm[,10:ncol(table_res_farm)])


summary(lm(Estimate~year, data=farm_nfac$trend_group2[farm_nfac$trend_group2$group=="g1",]))$coef
summary(lm(MSI~year, data=RES_farmland))$coef
cor.test(farm_nfac$trend_group2[farm_nfac$trend_group2$group=="g1","Estimate"],RES_farmland$MSI)
cor.test(farm_nfac$trend_group2[farm_nfac$trend_group2$group=="g1","Estimate"],farm_nfac$trend_group2[farm_nfac$trend_group2$group=="all","Estimate"])
cor.test(farm_nfac$data_msi[farm_nfac$data_msi$group=="g1","Index_c"],farm_nfac$data_msi[farm_nfac$data_msi$group=="all","Index_c"])


table_res_foresta <- dcast(data=forest_nfac$data_loadings, name_long~variable, id.vars="value")
table_res_forestb <- forest_nfac$group[[1]][[1]]
table_res_forest <- merge(table_res_forestb[,c("name_long", "group", "uncert")], table_res_foresta, by="name_long", all.x=T )
table_res_forest$pvalue <- table_res_forest$st_error <- table_res_forest$slope <-  NA
for(i in 1:nrow(table_res_forest)){
  res_lm <- summary(lm(pred.value~Year,
                       data=forest_nfac$data_to_plot_sp[forest_nfac$data_to_plot_sp$code_sp==unique(forest_nfac$data_to_plot_sp$code_sp)[i],]))$coef
  table_res_forest$slope[i] <- res_lm[2,1]
  table_res_forest$st_error[i] <- res_lm[2,2]
  table_res_forest$pvalue[i] <- res_lm[2,4]
}
table_res_forestc <- forest_nfac$exp_var_lt[,1:6]
names(table_res_forestc) <- c("name_long","% Latent trend 1","% Latent trend 2",
                            "% Latent trend 3","% Latent trend 4","% eta")
table_res_forest <- merge(table_res_forest, table_res_forestc, by="name_long", all.x=T )
table_res_forest[,3:ncol(table_res_forest)] <- round(table_res_forest[,3:ncol(table_res_forest)],4)
table_res_forest[,11:ncol(table_res_forest)] <- abs(table_res_forest[,11:ncol(table_res_forest)])

summary(lm(Estimate~year, data=forest_nfac$trend_group2[forest_nfac$trend_group2 $group=="g1",]))$coef
summary(lm(MSI~year, data=RES_forest))$coef
cor.test(forest_nfac$trend_group2[forest_nfac$trend_group2$group=="g1","Estimate"],RES_forest$MSI)
cor.test(forest_nfac$trend_group2[forest_nfac$trend_group2$group=="g1","Estimate"],forest_nfac$trend_group2[forest_nfac$trend_group2$group=="all","Estimate"])
cor.test(forest_nfac$data_msi[forest_nfac$data_msi$group=="g1","Index_c"],forest_nfac$data_msi[forest_nfac$data_msi$group=="all","Index_c"])
summary(lm(Index_c~year, data=forest_nfac$data_msi[forest_nfac$data_msi$group=="g1",]))$coef


table_res_alla <- dcast(data=all_nfac[[3]], name_long~variable, id.vars="value")
table_res_allb <- all_nfac[[12]][[1]][[1]]
table_res_all <- merge(table_res_allb[,c("name_long", "group", "uncert")], table_res_alla, by="name_long", all.x=T )
table_res_all$pvalue <- table_res_all$st_error <- table_res_all$slope <-  NA
for(i in 1:nrow(table_res_all)){
  res_lm <- summary(lm(pred.value~Year,
                       data=all_nfac[[1]][all_nfac[[1]]$code_sp==unique(all_nfac[[1]]$code_sp)[i],]))$coef
  table_res_all$slope[i] <- res_lm[2,1]
  table_res_all$st_error[i] <- res_lm[2,2]
  table_res_all$pvalue[i] <- res_lm[2,4]
}
table_res_all[,3:25] <- round(table_res_all[,3:25],4)

summary(lm(Estimate~year, data=all_nfac[[13]][all_nfac[[13]]$group=="g1",]))$coef
summary(lm(MSI~year, data=RES_all))$coef
cor.test(all_nfac[[13]][all_nfac[[13]]$group=="g1","Estimate"],RES_all$MSI)
cor.test(all_nfac[[13]][all_nfac[[13]]$group=="g1","Estimate"],all_nfac[[14]][all_nfac[[14]]$group=="all","Estimate"])


all_farm_for <- rbind(data.frame(all_nfac[[2]],group="all"),
                      data.frame(farm_nfac[[2]],group="farm"),
                      data.frame(forest_nfac[[2]],group="forest"))

lt1 <- ggplot(all_farm_for,aes(Year,value)) +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 1" & all_farm_for$group=="farm",], col="red") +
  geom_line(aes(y=-scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 8" & all_farm_for$group=="all",], col="black") +
  annotate("text", label = paste0("rho = ", abs(round(cor(farm_nfac[[2]][farm_nfac[[2]]$variable=="Latent trend 1",]$value,
                                                      all_nfac[[2]][all_nfac[[2]]$variable=="Latent trend 8",]$value),3))), x=2010, y=-1) +
  ggtitle("FB latent trend 1 and CB latent trend 8") +
  theme_modern() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

lt2 <- ggplot(all_farm_for,aes(Year,value)) +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 2" & all_farm_for$group=="farm",], col="red") +
  geom_line(aes(y=-scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 10" & all_farm_for$group=="all",], col="black") +
  annotate("text", label = paste0("rho = ", abs(round(cor(farm_nfac[[2]][farm_nfac[[2]]$variable=="Latent trend 2",]$value,
                                                          all_nfac[[2]][all_nfac[[2]]$variable=="Latent trend 10",]$value),3))), x=2010, y=-1) +
  ggtitle("FB latent trend 2 and CB latent trend 10") +
  theme_modern() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

lt3 <- ggplot(all_farm_for,aes(Year,value)) +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 3" & all_farm_for$group=="farm",], col="red") +
  geom_line(aes(y=-scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 17" & all_farm_for$group=="all",], col="black") +
  annotate("text", label = paste0("rho = ", abs(round(cor(farm_nfac[[2]][farm_nfac[[2]]$variable=="Latent trend 3",]$value,
                                                          all_nfac[[2]][all_nfac[[2]]$variable=="Latent trend 17",]$value),3))), x=2010, y=-1) +
  ggtitle("FB latent trend 3 and CB latent trend 17") +
  theme_modern() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

lt4 <- ggplot(all_farm_for,aes(Year,value)) +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 1" & all_farm_for$group=="forest",], col="green") +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 11" & all_farm_for$group=="all",], col="black") +
  annotate("text", label = paste0("rho = ", abs(round(cor(forest_nfac[[2]][forest_nfac[[2]]$variable=="Latent trend 1",]$value,
                                                          all_nfac[[2]][all_nfac[[2]]$variable=="Latent trend 11",]$value),3))), x=2010, y=-1) +
  ggtitle("WB latent trend 1 and CB latent trend 11") +
  theme_modern() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

lt5 <- ggplot(all_farm_for,aes(Year,value)) +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 2" & all_farm_for$group=="forest",], col="green") +
  geom_line(aes(y=-scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 7" & all_farm_for$group=="all",], col="black") +
  annotate("text", label = paste0("rho = ", abs(round(cor(forest_nfac[[2]][forest_nfac[[2]]$variable=="Latent trend 2",]$value,
                                                          all_nfac[[2]][all_nfac[[2]]$variable=="Latent trend 7",]$value),3))), x=2010, y=1) +
  ggtitle("WB latent trend 2 and CB latent trend 7") +
  theme_modern() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

lt6 <- ggplot(all_farm_for,aes(Year,value)) +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 3" & all_farm_for$group=="forest",], col="green") +
  geom_line(aes(y=-scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 18" & all_farm_for$group=="all",], col="black") +
  annotate("text", label = paste0("rho = ", abs(round(cor(forest_nfac[[2]][forest_nfac[[2]]$variable=="Latent trend 3",]$value,
                                                          all_nfac[[2]][all_nfac[[2]]$variable=="Latent trend 18",]$value),3))), x=2010, y=1) +
  ggtitle("WB latent trend 3 and CB latent trend 18") +
  theme_modern() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))

lt7 <- ggplot(all_farm_for,aes(Year,value)) +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 4" & all_farm_for$group=="forest",], col="green") +
  geom_line(aes(y=scale(value)),data=all_farm_for[all_farm_for$variable=="Latent trend 9" & all_farm_for$group=="all",], col="black") +
  annotate("text", label = paste0("rho = ", abs(round(cor(forest_nfac[[2]][forest_nfac[[2]]$variable=="Latent trend 4",]$value,
                                                          all_nfac[[2]][all_nfac[[2]]$variable=="Latent trend 9",]$value),3))), x=2010, y=1) +
  ggtitle("WB latent trend 4 and CB latent trend 9") +
  theme_modern() + theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))


png(file="output/saving_plot.png", width=3000, height=2200,res = 300)
multiplot(lt1, lt2, lt3, lt4,
          lt5, lt6, lt7, cols=2)
dev.off()



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
                                                              seed= x$num_sim)},
                          .parallel = T)
saveRDS(rand_nfac_test_g,"output/rand_nfac_test_g_new.rds")

sim_data <- data.frame(num_sim=1:(rep_sim*4),nb_group_exp=sort(rep(c(1,2,3,4),rep_sim)), equal=FALSE)
registerDoParallel(cores=no_cores)
rand_nfac_test_g2 <- dlply(sim_data,.(num_sim),
                          .fun=function(x){simul_rand_dfa(nb_group_exp = x$nb_group_exp,
                                                          seed= x$num_sim,
                                                          equi = FALSE)},
                          .parallel = T)
saveRDS(rand_nfac_test_g2,"output/rand_nfac_test_g2_new.rds")

sim_data <- data.frame(num_sim=1:(rep_sim*5),sd_ci=sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)))
registerDoParallel(cores=no_cores)
rand_nfac_test_l <- dlply(sim_data,.(num_sim),
                           .fun=function(x){simul_rand_dfa(sd_ci = x$sd_ci,
                                                           seed= x$num_sim)},
                           .parallel = T)
saveRDS(rand_nfac_test_l,"output/rand_nfac_test_l_new.rds")

sim_data <- data.frame(num_sim=1:(rep_sim*5),sd_ci=sort(rep(c(0.01,0.05,0.2,0.5,1),rep_sim)), equal=FALSE)
registerDoParallel(cores=no_cores)
rand_nfac_test_l2 <- dlply(sim_data,.(num_sim),
                          .fun=function(x){simul_rand_dfa(sd_ci = x$sd_ci,
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
### Complementary test with life history traits

trait <- read.table("raw_data/Life-history characteristics of European birds.txt", header = T, sep = "\t")

data_for_trait <- merge(farm_nfac$group[[1]][[1]],
                        trait[,c("Species","LengthU_MEAN","WeightU_MEAN","Clutch_MEAN",
                                 "Nest.type","Life.span","Long.distance.migrant",
                                 "Broods.per.year","Incubation.period",
                                 "Age.of.first.breeding","Sedentary",
                                 "Folivore_Y", "Frugivore_Y", "Granivore_Y", "Arthropods_Y",
                                 "Folivore_B", "Frugivore_B", "Granivore_B", "Arthropods_B",
                                 "Deciduous.forest", "Coniferous.forest", "Woodland",
                                 "Shrub","Savanna","Tundra","Grassland",
                                 "Mountain.meadows","Reed","Human.settlements")],
                        by.x="name_long", by.y="Species", all.x=T)
#data_for_trait <- data_for_trait[data_for_trait$group!=3,]
data_for_trait$group <- as.character(data_for_trait$group)
data_for_trait$Long.distance.migrant <- as.character(data_for_trait$Long.distance.migrant)
summary(lm(WeightU_MEAN~group,data_for_trait))
summary(lm(Clutch_MEAN~group,data_for_trait))
summary(lm(Life.span~group,data_for_trait))
summary(lm(Broods.per.year~group,data_for_trait))
summary(lm(Incubation.period~group,data_for_trait))
summary(lm(Age.of.first.breeding~group,data_for_trait))
table(data_for_trait$Nest.type,data_for_trait$group)
table(data_for_trait$Long.distance.migrant,data_for_trait$group)
table(data_for_trait$Sedentary,data_for_trait$group)
table(data_for_trait$Frugivore_Y,data_for_trait$group)
table(data_for_trait$Granivore_Y,data_for_trait$group)
table(data_for_trait$Arthropods_Y,data_for_trait$group)
table(data_for_trait$Granivore_B,data_for_trait$group)
table(data_for_trait$Arthropods_B,data_for_trait$group)

data_for_trait[,20:28] <- apply(data_for_trait[,20:28],2,function(x){return(as.character(x))})

summary(lm(PC2~WeightU_MEAN+Clutch_MEAN+Life.span+Long.distance.migrant+
  Broods.per.year+Incubation.period+Age.of.first.breeding+Sedentary+
  Frugivore_Y+Granivore_Y+Arthropods_Y,data_for_trait))

library(RVAideMemoire)

m1 <- lm(PC2~LengthU_MEAN+Clutch_MEAN+Life.span+Long.distance.migrant+
     Broods.per.year+Incubation.period+Age.of.first.breeding+Nest.type+
     Granivore_B+Arthropods_B+Shrub+Grassland+Human.settlements,data_for_trait)

Anova(m1)


m1 <- lm(PC2~LengthU_MEAN+Long.distance.migrant+
           Incubation.period+Nest.type+
           Granivore_B+Arthropods_B,data_for_trait)
summary(m1)

data_for_trait <- merge(forest_nfac$group[[1]][[1]],
                        trait[,c("Species","LengthU_MEAN","WeightU_MEAN","Clutch_MEAN",
                                 "Nest.type","Life.span","Long.distance.migrant",
                                 "Broods.per.year","Incubation.period",
                                 "Age.of.first.breeding","Sedentary",
                                 "Folivore_Y", "Frugivore_Y", "Granivore_Y", "Arthropods_Y",
                                 "Folivore_B", "Frugivore_B", "Granivore_B", "Arthropods_B",
                                 "Deciduous.forest", "Coniferous.forest", "Woodland",
                                 "Shrub","Savanna","Tundra","Grassland",
                                 "Mountain.meadows","Reed","Human.settlements")],
                        by.x="name_long", by.y="Species", all.x=T)
data_for_trait$group <- as.character(data_for_trait$group)
data_for_trait$Long.distance.migrant <- as.character(data_for_trait$Long.distance.migrant)
data_for_trait[,22:30] <- apply(data_for_trait[,22:30],2,function(x){return(as.character(x))})


m1 <- lm(PC1~LengthU_MEAN+Clutch_MEAN+Life.span+Long.distance.migrant+
           Broods.per.year+Age.of.first.breeding+Nest.type+
           Granivore_B+Arthropods_B+Deciduous.forest+
           Coniferous.forest+Woodland,data_for_trait)

Anova(m1)

m1 <- lm(PC1~Clutch_MEAN+
           Broods.per.year+Age.of.first.breeding+Nest.type+
           Granivore_B+Arthropods_B,data_for_trait)
summary(m1)


# extract % variance of sp ts explained by latent trends

exp_var_lt <- farm_nfac$data_loadings[,c("variable","value","name_long")]
exp_var_lt <- dcast(exp_var_lt, name_long~variable, value.var = "value", fun.aggregate = sum)
eta_sp <- data.frame(name_long=species_farm$name_long, eta=farm_nfac$sdRep[!grepl("log_re_sp", row.names(farm_nfac$sdRep)) & grepl("re_sp", row.names(farm_nfac$sdRep)) ,1])
exp_var_lt <- merge(exp_var_lt,eta_sp, by="name_long", all.x=T)

exp_var_lt$all <- apply(exp_var_lt[,-1],1,function(x){return(sum(abs(x)))})
exp_var_lt[,2:(ncol(exp_var_lt)-1)] <- exp_var_lt[,2:(ncol(exp_var_lt)-1)]/exp_var_lt$all
exp_var_lt$name_long <- fct_reorder(exp_var_lt$name_long,exp_var_lt$eta)
exp_var_lt_long <- melt(exp_var_lt[,1:(ncol(exp_var_lt)-1)])

ggplot(exp_var_lt_long) + 
  geom_col(aes(value, name_long, fill=variable)) +
  facet_wrap(variable ~ ., ncol=length(unique(exp_var_lt_long$variable))) +
  theme_modern() + 
  theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face="italic"))
