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
  ylim(70, NA)+ ylab("MSI (1996= 100 )")+
  geom_ribbon(aes(ymin=lower_CL_trend, ymax=upper_CL_trend), alpha=0.2)+
  geom_line(aes(y=Trend), colour="dark green", size=1)+ theme_modern() +
  geom_pointrange(aes(ymax = MSI+sd_MSI, ymin=MSI-sd_MSI), colour="dark green")

# MSI with DFA results
farm_nfac3 <- readRDS("output/farm_nfac3.rds")
farm_nfac3_reshape <- data.frame(code_sp=farm_nfac3[[1]]$code_sp,
                                 year=farm_nfac3[[1]]$Year,
                                 relative_abundance=farm_nfac3[[1]]$pred.value, # to multiply by sd and add mean orgi values
                                 CI_inf=farm_nfac3[[1]]$pred.value-farm_nfac3[[1]]$pred_se.value*1.96)

RES_DFA_farmland <- MSI_MC_func(farm_nfac3_reshape)
