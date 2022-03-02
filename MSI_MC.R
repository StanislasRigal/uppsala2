# MultiSpecies Index with uncertainty using Monte-Carlo approach
require(ggplot2)

# time-series from the SBBS
ts_bird_se_allcountry <- readRDS("output/ts_bird_se_allcountry.rds")
ts_bird_se_allcountry_clean <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance)),]
species_data <- readRDS("output/species_data.rds")

ts_farmland <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                             c("VANVAN","NUMARQ","ALAARV","HIRRUS",
                                               "MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                               "LANCOL","STUVUL","LINCAN","EMBCIT",
                                               "PASMON","CORFRU","ANTPRA","EMBHOR"),])

ts_forest <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                                      c("ACCNIS","TETBON","TRIOCH","COLOEN",
                                                        "DENMAJ","DRYMAR","PICVIR","JYNTOR",
                                                        "DRYMIN","PICTRI","NUCCAR","GARGLA",
                                                        "PERATE","LOPCRI","POEPAL","POEMON",
                                                        "SITEUR","CERFAM","TURVIS","PHOPHO",
                                                        "PHYCOL","PHYSIB","REGREG","FICHYP",
                                                        "ANTTRI","COCCOC","SPISPI","PYRPYR","EMBRUS"),])



# load main function
source("function_ts.R")

RES_farmland <- MSI_MC_func(ts_farmland)

RES_forest <- MSI_MC_func(ts_forest)

ggplot(RES_farmland, aes(x=year, y=MSI))+
  geom_point(colour = "dark green",size=3)+
  ylim(50, NA)+ ylab("MSI (1996= 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="dark green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="dark green")

# MSI with DFA results
farm_nfac3 <- readRDS("output/farm_nfac3.rds")

prep_data_msi <- function(x){
  x <- droplevels(x)
  x_reshape <- data.frame(code_sp=x$code_sp[1],
                          year=c(1996,x$Year),
                          relative_abundance=c(1,x$pred.value*sd(c(1,x$value_orig), na.rm=T)+mean(c(1,x$value_orig), na.rm=T)),
                          CI_inf=c(1,x$pred.value*sd(c(1,x$value_orig), na.rm=T)+mean(c(1,x$value_orig), na.rm=T)-x$pred_se.value*1.96*sd(c(1,x$value_orig), na.rm=T)))
  return(x_reshape)
}

farm_nfac3_reshape <- ddply(as.data.frame(farm_nfac3[[1]]), .(code_sp), .fun = prep_data_msi, .progress="text")

RES_DFA_farmland <- MSI_MC_func(farm_nfac3_reshape)

ggplot(RES_DFA_farmland, aes(x=year, y=MSI))+
  geom_point(colour = "dark green",size=3)+
  ylim(50, NA)+ ylab("MSI (1996 = 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="dark green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="dark green")

ggplot(RES_DFA_farmland, aes(x=year, y=MSI))+
  geom_point(colour = "green",size=3)+
  ylim(50, NA)+ ylab("MSI (1996 = 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="green")+
  geom_point(data=RES_farmland,colour = "blue",size=3)+
  geom_ribbon(data=RES_farmland,aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(data=RES_farmland,aes(y=Trend), colour="blue", size=1)+ 
  geom_pointrange(data=RES_farmland,aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="blue")

# Influence of each species on MSI

mean_load <- as.data.frame(as.data.frame(farm_nfac3[[3]]) %>% 
                             group_by(variable) %>% summarize(mean_load=mean(value)))

trend_dfa <- as.matrix(dcast(as.data.frame(farm_nfac3[[2]])[,c("Year","variable","rot_tr.value")],
                             Year~variable, value.var = "rot_tr.value"))

msi_from_trend <- trend_dfa[,-1] %*% (matrix(mean_load$mean_load, ncol=1)/sum(matrix(mean_load$mean_load, ncol=1)))
