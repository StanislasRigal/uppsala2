# MultiSpecies Index with uncertainty using Monte-Carlo approach
require(ggplot2)

# time-series from the SBBS
ts_bird_se_allcountry <- readRDS("output/ts_bird_se_allcountry.rds")
ts_bird_se_allcountry_clean <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance)),]
species_data <- readRDS("output/species_data.rds")

ts_farmland <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                             c( "FALTIN","VANVAN","ALAARV","HIRRUS","CORFRU",
                                                "SAXRUB","SYLCOM","ANTPRA","MOTFLA","LANCOL",
                                                "STUVUL","LINCAN","EMBCIT","EMBHOR","PASMON"),])
                                               #"VANVAN","NUMARQ","ALAARV","HIRRUS",
                                               #"MOTFLA","OENOEN","SAXRUB","SYLCOM",
                                               #"LANCOL","STUVUL","LINCAN","EMBCIT",
                                               #"PASMON","CORFRU","ANTPRA","EMBHOR"),])

ts_forest <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                                      c("ACCNIS","TETBON","TRIOCH","COLOEN",
                                                        "DENMAJ","DRYMAR","PICVIR","JYNTOR",
                                                        "DRYMIN","PICTRI","NUCCAR","GARGLA",
                                                        "PERATE","LOPCRI","POEPAL","POEMON",
                                                        "SITEUR","CERFAM","TURVIS","PHOPHO",
                                                        "PHYCOL","PHYSIB","REGREG","FICHYP",
                                                        "ANTTRI","COCCOC","SPISPI","PYRPYR","EMBRUS"),])

ts_all <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                                   species_all$code_sp,])



# load main function
source("function_ts.R")

RES_farmland <- MSI_MC_func(ts_farmland,SEbaseyear=1998, plotbaseyear=1998)

RES_forest <- MSI_MC_func(ts_forest,SEbaseyear=1998, plotbaseyear=1998)

RES_all <- MSI_MC_func(ts_all,SEbaseyear=1998, plotbaseyear=1998)

ggplot(RES_farmland, aes(x=year, y=MSI))+
  geom_point(colour = "dark green",size=3)+
  ylim(70, NA)+ ylab("MSI (1998 = 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="dark green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="dark green")


ggplot(RES_farmland, aes(x=year, y=MSI))+
  geom_point(size=3)+
  ylim(75, NA)+ ylab("MSI (1998 = 100 )")+ xlab("Year")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI))

ggsave("output/farm_msi.png",
       dpi=300,
       width = 6,
       height = 4
)

ggplot(RES_forest, aes(x=year, y=MSI))+
  geom_point(size=3)+
  ylim(90, NA)+ ylab("MSI (1998 = 100 )")+ xlab("Year")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI))

ggsave("output/forest_msi.png",
       dpi=300,
       width = 6,
       height = 4
)

ggplot(RES_all, aes(x=year, y=MSI))+
  geom_point(size=3)+
  ylim(90, NA)+ ylab("MSI (1998 = 100 )")+ xlab("Year")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI))

ggsave("output/all_msi.png",
       dpi=300,
       width = 6,
       height = 4
)

ts_farmland <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                                        c( "FALTIN","ALAARV","CORFRU",
                                                           "SAXRUB","ANTPRA","MOTFLA",
                                                           "LINCAN","EMBHOR","PASMON"),])
RES_farmland1 <- MSI_MC_func(ts_farmland,SEbaseyear=1998, plotbaseyear=1998)

ts_farmland <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                                        c( "FALTIN","VANVAN","HIRRUS","CORFRU",
                                                           "SYLCOM","LANCOL",
                                                           "STUVUL","LINCAN","EMBCIT","EMBHOR"),])
RES_farmland2 <- MSI_MC_func(ts_farmland,SEbaseyear=1998, plotbaseyear=1998)

ts_farmland <- droplevels(ts_bird_se_allcountry_clean[ts_bird_se_allcountry_clean$code_sp %in%
                                                        c( "VANVAN","ALAARV","HIRRUS",
                                                           "SAXRUB","SYLCOM","ANTPRA","MOTFLA","LANCOL",
                                                           "STUVUL","EMBCIT","EMBHOR","PASMON"),])
RES_farmland3 <- MSI_MC_func(ts_farmland,SEbaseyear=1998, plotbaseyear=1998)

RES_farmland_all <- rbind(RES_farmland,RES_farmland1,RES_farmland2,RES_farmland3)
RES_farmland_all$group <- c(rep("All",23),rep("Without cluster 1",23),
                            rep("Without cluster 2",23),rep("Without cluster 3",23))

ggplot(RES_farmland_all, aes(x=year, y=MSI, fill=group))+
  ylim(75, NA)+ ylab("MSI (1998 = 100 )")+ xlab("Year")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend, col=group), size=1)+ theme_modern() +
  theme(legend.title = element_blank())

ggsave("output/farm_msi_compare.png",
       dpi=300,
       width = 8,
       height = 4
)

# MSI with DFA results
farm_nfac3 <- readRDS("output/farm_nfac3.rds")

prep_data_msi <- function(x){
  x <- droplevels(x)
  x_reshape <- data.frame(code_sp=x$code_sp[1],
                          year=x$Year,
                          #relative_abundance=c(x$pred.value*sd(x$value_orig, na.rm=T)+mean(x$value_orig, na.rm=T)),
                          #CI_inf=c(x$pred.value*sd(x$value_orig, na.rm=T)+mean(x$value_orig, na.rm=T)-x$pred_se.value*1.96*sd(x$value_orig, na.rm=T)))
                          relative_abundance=x$pred.value,
                          CI_inf=c(x$pred.value-1.96*x$pred_se.value))
  return(x_reshape)
}

farm_nfac3_reshape <- ddply(as.data.frame(farm_nfac3[[1]]), .(code_sp), .fun = prep_data_msi, .progress="text")

RES_DFA_farmland <- MSI_MC_func(farm_nfac3_reshape, SEbaseyear=1998, plotbaseyear=1998)

ggplot(RES_DFA_farmland, aes(x=year, y=MSI))+
  geom_point(colour = "dark green",size=3)+
  ylim(50, NA)+ ylab("MSI (1998 = 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="dark green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="dark green")

ggplot(RES_DFA_farmland, aes(x=year, y=MSI))+
  geom_point(colour = "green",size=3)+
  ylim(50, NA)+ ylab("MSI (1998 = 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="green")+
  geom_point(data=RES_farmland,colour = "blue",size=3)+
  geom_ribbon(data=RES_farmland,aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(data=RES_farmland,aes(y=Trend), colour="blue", size=1)+ 
  geom_pointrange(data=RES_farmland,aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="blue")

forest_nfac4_reshape <- ddply(as.data.frame(forest_nfac4[[1]]), .(code_sp), .fun = prep_data_msi, .progress="text")

RES_DFA_forest <- MSI_MC_func(forest_nfac4_reshape, SEbaseyear=1998, plotbaseyear=1998)

ggplot(RES_DFA_forest, aes(x=year, y=MSI))+
  geom_point(colour = "green",size=3)+
  ylim(50, NA)+ ylab("MSI (1998 = 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="green")+
  geom_point(data=RES_forest,colour = "blue",size=3)+
  geom_ribbon(data=RES_forest,aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(data=RES_forest,aes(y=Trend), colour="blue", size=1)+ 
  geom_pointrange(data=RES_forest,aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="blue")

# Influence of each species on MSI and groups of species

mean_load <- as.data.frame(as.data.frame(farm_nfac3[[3]]) %>% 
                             group_by(variable) %>% summarize(mean_load=mean(value)))

trend_dfa <- as.matrix(dcast(as.data.frame(farm_nfac3[[2]])[,c("Year","variable","rot_tr.value")],
                             Year~variable, value.var = "rot_tr.value"))

msi_from_trend <- trend_dfa[,-1] %*% (matrix(mean_load$mean_load, ncol=1)/sum(matrix(mean_load$mean_load, ncol=1)))

ggplot(RES_DFA_farmland[-1,], aes(x=year, y=Zscore(MSI))) +
  geom_point(colour = "green",size=3) +
  geom_point(data=data.frame(x=1998:2020,y=Zscore(msi_from_trend)),
             aes(x,y), colour = "blue",size=3) + theme_modern()

dev_from_load <- ddply(as.data.frame(farm_nfac3[[3]]), .(code_sp),
                       .fun = function(x){dev1 <- x$value[x$variable=="X1"]-mean_load$mean_load[mean_load$variable=="X1"]
                       dev2 <- x$value[x$variable=="X2"]-mean_load$mean_load[mean_load$variable=="X2"]
                       dev3 <- x$value[x$variable=="X3"]-mean_load$mean_load[mean_load$variable=="X3"]
                       return(data.frame(dev1,dev2,dev3))})
dev_from_load$sig <- paste0(sign(dev_from_load$dev1),sign(dev_from_load$dev2),
                            sign(dev_from_load$dev3))
  