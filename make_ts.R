# look at abundance for each species

bird_se <- readRDS("output/bird_se.rds")

# Use data from 1998 because low number of routes monitored in 1996 and 1997: https://doi.org/10.34080/os.v17.22684 )

bird_se_1998 <- droplevels(bird_se[bird_se$year>1997,])
saveRDS(bird_se_1998,"output/bird_se_1998.rds")

ab_sp <- as.data.frame(bird_se_1998 %>% group_by(code_sp) %>% summarize(ab_tot=sum(abund)))

ggplot(ab_sp, aes(x=reorder(code_sp, -ab_tot, sum), y=ab_tot)) + 
  geom_bar(stat = "identity")+ theme_classic() + scale_y_continuous(trans = "log")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Species")+
  ylab("Abundance (log)")

# time-series from the SBBS

## example for one species

data_ex <- droplevels(bird_se_1998[bird_se_1998$code_sp=="CORFRU",])
timestep <- length(levels(as.factor(data_ex$year)))-1

glm1 <- glm(abund~as.factor(code_route)+as.factor(year), data=data_ex, family = quasipoisson)
sglm1 <- summary(glm1)

## as link function is log, estimates need to be back transformed

coef_yr <- tail(matrix(sglm1$coefficients[,1]), timestep)
coef_yr <- rbind(1, exp(coef_yr))

error_yr <- tail(matrix(sglm1$coefficients[,2]), timestep)
error_yr <- rbind(0, error_yr)

## CIs

require(arm)
glm1.sim <- sim(glm1)
lc_inf_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .025), timestep)))
lc_sup_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .975), timestep)))

## plot

require(see)

tab_plot <- data.frame(year=as.numeric(as.character(levels(as.factor(data_ex$year)))),
                       value=coef_yr, error=error_yr, LL=lc_inf_sim, UL=lc_sup_sim)
ggplot(tab_plot, aes(year, value)) + geom_line() +
  geom_line(aes(y=LL), linetype="dashed") +
  geom_line(aes(y=UL), linetype="dashed") +
  ylab("Relative abundance") + xlab("Years") +
  theme_modern()

## for all species

### for loop from French birds analysis

listSp<-levels(as.factor(bird_se_1998$code_sp))
  
i<-0
for (sp in listSp) {
    print(i)
  
    i <- i + 1
    
    ## d data for species i
    sp<-listSp[i]
    d <- droplevels(bird_se_1998[bird_se_1998$code_sp==sp,])
    
    ## Occurrence
    ## number of route followed by year
    nb_route <- tapply(rep(1,nrow(d)),d$year,sum)
    ## number of route with species i by year
    nb_route_presence <- tapply(ifelse(d$abund>0,1,0),d$year,sum)
    year<-as.numeric(as.character(levels(as.factor(d$year))))
    firstY <- min(year)
    lastY <- max(year)
    timestep <- length(year)-1
    
    ## table for analysis result
    threshold_occurrence<-3
    tab_ana <- data.frame(year=rep(year,2),val=c(nb_route,nb_route_presence),LL = NA,UL=NA,
                       catPoint=NA,pval=NA,
                       curve=rep(c("route","presence"),each=length(year)))
    tab_ana$catPoint <- ifelse(tab_ana$val == 0,"0", ifelse(tab_ana$val < threshold_occurrence,
                                                           "inf_threshold",NA))
    
    ## raw abundance
    ## abundance by year
    abund <- tapply(d$abund,d$year,sum)
    ## table for figure
    threshold_abundance <- 5
    tab_fig <- data.frame(year=year,val=abund,LL = NA,UL=NA,catPoint=NA,pval=NA)
    tab_fig$catPoint <- ifelse(tab_fig$val == 0,"0",ifelse(tab_fig$val < threshold_abundance,
                                                           "inf_threshold",NA))
    
    # remove criteria
    remove_sp<-FALSE
    
    ## if first year empty
    if(tab_fig$val[1]==0){remove_sp<-TRUE}
    
    ## if four consecutive years empty
    ab_vec<-paste(tab_fig$val,collapse="")
    if(str_detect(ab_vec, "0000")){remove_sp<-TRUE}
    
    ## if less than consecutive years
    ab_vec2<-paste(sign(tab_fig$val),collapse="")
    if(!str_detect(ab_vec2, "111")){remove_sp<-TRUE}
    
    if(anyNA(tab_fig$catPoint) & anyNA(tab_ana$catPoint[tab_ana$curve=="presence"]) & remove_sp==F){
      
      ## GLM abundance variation
      glm1 <- glm(abund~as.factor(code_route)+as.factor(year),data=d,family=quasipoisson)
      sglm1 <- summary(glm1)
      
      ## as link function is log, estimates need to be back transformed
      coef_yr <- tail(matrix(sglm1$coefficients[,1]), timestep)
      coef_yr <- rbind(1, exp(coef_yr))
      error_yr <- tail(matrix(sglm1$coefficients[,2]), timestep)
      error_yr <- rbind(0, error_yr)
      pval <- c(1,tail(matrix(coefficients(sglm1)[,4]),timestep))
      
      ## CIs
      glm1.sim <- sim(glm1)
      ci_inf_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .025), timestep)))
      ci_sup_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .975), timestep)))

      ## table for result and figures
      thresold_signif <- 0.05
      tab_res <- data.frame(year,val=coef_yr,
                         LL=ci_inf_sim,UL=ci_sup_sim,
                         catPoint=ifelse(pval<thresold_signif,"significatif",NA),pval)
      ## cleaning out of range CIs			   
      tab_res$UL <- ifelse(nb_route_presence==0,NA,tab_res$UL)
      tab_res$UL <-  ifelse(tab_res$UL == Inf, NA,tab_res$UL)
      tab_res$UL <-  ifelse(tab_res$UL > 1.000000e+20, NA,tab_res$UL)
      tab_res$UL[1] <- 1
      tab_res$val <-  ifelse(tab_res$val > 1.000000e+20,1.000000e+20,tab_res$val)
      
      ## overdispersion index
      dispAn <- sglm1$deviance/sglm1$null.deviance
      
      ## class uncertainity
      if(dispAn > 2 | (median(nb_route_presence)<threshold_occurrence & median(abund)<threshold_abundance)) catIncert <- "Uncertain" else catIncert <-"Good"
      vecLib <-  NULL
      if(dispAn > 2 | median(nb_route_presence)<threshold_occurrence) {
        if(median(nb_route_presence)<threshold_occurrence) {
          vecLib <- c(vecLib,"too rare species")
        }
        if(dispAn > 2) {
          vecLib <- c(vecLib,"deviance")
        }
      }
      reason_uncert <-  paste(vecLib,collapse=" and ")
      
      ## table for saving results      
      tab_tot <- data.frame(code_sp=sp, year = tab_res$year, nb_year=timestep,
                          firstY=firstY, lastY=lastY,
                          relative_abundance=round(tab_res$val,3),
                          CI_inf = round(tab_res$LL,3), CI_sup = round(tab_res$UL,3),
                          Standard_error = round(error_yr,4),
                          p_value = round(tab_res$pval,3), signif = !is.na(tab_res$catPoint),
                          nb_route,nb_route_presence,abundance=abund,
                          mediane_occurrence=median(nb_route_presence), mediane_ab=median(abund) ,
                          valid = catIncert, uncertanity_reason = reason_uncert)
      
    }
    else{
      tab_tot <- data.frame(code_sp=sp, year = tab_res$year, nb_year=timestep,
                            firstY=firstY, lastY=lastY,
                            relative_abundance=NA,
                            CI_inf = NA, CI_sup = NA,
                            Standard_error = NA,
                            p_value = NA, signif = NA,
                            nb_route,nb_route_presence,abundance=abund,
                            mediane_occurrence=median(nb_route_presence), mediane_ab=median(abund) ,
                            valid = NA, uncertanity_reason = NA)
    }
    
    
    if(sp==listSp[1]) {
      glm_res <- tab_tot
    } else  {
      glm_res <- rbind(glm_res,tab_tot)
    }
} 

# compare with https://www.fageltaxering.lu.se/resultat/trender/allatrendertillsammans
ggplot(glm_res[glm_res$code_sp==levels(as.factor(glm_res$code_sp))[1],], aes(year, relative_abundance)) +
  geom_line() + geom_line(aes(y=CI_inf), linetype="dashed") + geom_line(aes(y=CI_sup), linetype="dashed") +
  ylab("Relative abundance") + xlab("Years") + theme_modern()


# function for Swedish data

get_ts <- function(data_bird_input){
  
    ## d data for species i
    d <- droplevels(data_bird_input)
    sp <- levels(as.factor(d$code_sp))
    
    ## Occurrence
    ## number of route followed by year
    nb_route <- tapply(rep(1,nrow(d)),d$year,sum)
    ## number of route with species i by year
    nb_route_presence <- tapply(ifelse(d$abund>0,1,0),d$year,sum)
    year<-as.numeric(as.character(levels(as.factor(d$year))))
    firstY <- min(year)
    lastY <- max(year)
    timestep <- length(year)-1
    
    ## table for analysis result
    threshold_occurrence<-3
    tab_ana <- data.frame(year=rep(year,2),val=c(nb_route,nb_route_presence),LL = NA,UL=NA,
                          catPoint=NA,pval=NA,
                          curve=rep(c("route","presence"),each=length(year)))
    tab_ana$catPoint <- ifelse(tab_ana$val == 0,"0", ifelse(tab_ana$val < threshold_occurrence,
                                                            "inf_threshold",NA))
    
    ## raw abundance
    ## abundance by year
    abund <- tapply(d$abund,d$year,sum)
    ## table for figure
    threshold_abundance <- 5
    tab_fig <- data.frame(year=year,val=abund,LL = NA,UL=NA,catPoint=NA,pval=NA)
    tab_fig$catPoint <- ifelse(tab_fig$val == 0,"0",ifelse(tab_fig$val < threshold_abundance,
                                                           "inf_threshold",NA))
    
    # remove criteria
    remove_sp<-FALSE
    
    ## if first year empty
    if(tab_fig$val[1]==0){remove_sp<-TRUE}
    
    ## if four consecutive years empty
    ab_vec<-paste(tab_fig$val,collapse="")
    if(str_detect(ab_vec, "0000")){remove_sp<-TRUE}
    
    ## if less than consecutive years
    ab_vec2<-paste(sign(tab_fig$val),collapse="")
    if(!str_detect(ab_vec2, "111")){remove_sp<-TRUE}
    
    if(anyNA(tab_fig$catPoint) & anyNA(tab_ana$catPoint[tab_ana$curve=="presence"]) & remove_sp==F){
      
      ## GLM abundance variation
      glm1 <- glm(abund~as.factor(code_route)+as.factor(year),data=d,family=quasipoisson)
      sglm1 <- summary(glm1)
      
      ## as link function is log, estimates need to be back transformed
      coef_yr <- tail(matrix(sglm1$coefficients[,1]), timestep)
      coef_yr <- rbind(1, exp(coef_yr))
      error_yr <- tail(matrix(sglm1$coefficients[,2]), timestep)
      error_yr <- rbind(0, error_yr)
      pval <- c(1,tail(matrix(coefficients(sglm1)[,4]),timestep))
      
      ## CIs
      glm1.sim <- sim(glm1)
      ci_inf_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .025), timestep)))
      ci_sup_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .975), timestep)))
      
      ## table for result and figures
      thresold_signif <- 0.05
      tab_res <- data.frame(year,val=coef_yr,
                            LL=ci_inf_sim,UL=ci_sup_sim,
                            catPoint=ifelse(pval<thresold_signif,"significatif",NA),pval)
      ## cleaning out of range CIs			   
      tab_res$UL <- ifelse(nb_route_presence==0,NA,tab_res$UL)
      tab_res$UL <-  ifelse(tab_res$UL == Inf, NA,tab_res$UL)
      tab_res$UL <-  ifelse(tab_res$UL > 1.000000e+20, NA,tab_res$UL)
      tab_res$UL[1] <- 1
      tab_res$val <-  ifelse(tab_res$val > 1.000000e+20,1.000000e+20,tab_res$val)
      
      ## overdispersion index
      dispAn <- sglm1$deviance/sglm1$null.deviance
      
      ## class uncertainity
      if(dispAn > 2 | (median(nb_route_presence)<threshold_occurrence & median(abund)<threshold_abundance)) catIncert <- "Uncertain" else catIncert <-"Good"
      vecLib <-  NULL
      if(dispAn > 2 | median(nb_route_presence)<threshold_occurrence) {
        
        if(median(nb_route_presence)<threshold_occurrence) {
          vecLib <- c(vecLib,"too rare species")
        }
        if(dispAn > 2) {
          vecLib <- c(vecLib,"deviance")
        }
      }
      reason_uncert <-  paste(vecLib,collapse=" and ")
      
      ## table for saving results      
      tab_tot <- data.frame(code_sp=sp, year = tab_res$year, nb_year=timestep,
                            firstY=firstY, lastY=lastY,
                            relative_abundance=round(tab_res$val,3),
                            CI_inf = round(tab_res$LL,3), CI_sup = round(tab_res$UL,3),
                            Standard_error = round(error_yr,4),
                            p_value = round(tab_res$pval,3), signif = !is.na(tab_res$catPoint),
                            nb_route,nb_route_presence,abundance=abund,
                            mediane_occurrence=median(nb_route_presence), mediane_ab=median(abund) ,
                            valid = catIncert, uncertanity_reason = reason_uncert)
      
    }
    else{
      tab_tot <- data.frame(code_sp=sp, year = year, nb_year=timestep,
                            firstY=firstY, lastY=lastY,
                            relative_abundance=NA,
                            CI_inf = NA, CI_sup = NA,
                            Standard_error = NA,
                            p_value = NA, signif = NA,
                            nb_route,nb_route_presence,abundance=abund,
                            mediane_occurrence=median(nb_route_presence), mediane_ab=median(abund) ,
                            valid = NA, uncertanity_reason = NA)
    }
    
    return(tab_tot)
} 

# require(plyr) but need to be done before loading dplyr, see package script

ts_test <- ddply(droplevels(bird_se_1998[bird_se_1998$code_sp %in% levels(as.factor(bird_se$code_sp))[1:5]]),
                 .(code_sp), .fun=get_ts, .progress="text")

ts_bird_se_allcountry <- ddply(bird_se_1998, .(code_sp), .fun=get_ts, .progress="text")

ts_bird_se_byreg <- ddply(bird_se_1998, .(code_sp, ln_kod), .fun=get_ts, .progress="text")
ts_bird_se_byreg_clean <- ts_bird_se_byreg[which(!is.na(ts_bird_se_byreg$relative_abundance)),]

ts_bird_se_byecoreg <- ddply(bird_se_1998, .(code_sp, ecoreg), .fun=get_ts, .progress="text")
ts_bird_se_byecoreg_clean <- ts_bird_se_byecoreg[which(!is.na(ts_bird_se_byecoreg$relative_abundance)),]

ts_bird_se_bysubecoreg <- ddply(bird_se_1998, .(code_sp, subecoreg), .fun=get_ts, .progress="text")
ts_bird_se_bysubecoreg_clean <- ts_bird_se_bysubecoreg[which(!is.na(ts_bird_se_bysubecoreg$relative_abundance)),]


# save outputs
saveRDS(ts_bird_se_allcountry, file = "output/ts_bird_se_allcountry.rds")
saveRDS(ts_bird_se_byreg, file = "output/ts_bird_se_byreg.rds")
saveRDS(ts_bird_se_byecoreg, file = "output/ts_bird_se_byecoreg.rds")
saveRDS(ts_bird_se_bysubecoreg, file = "output/ts_bird_se_bysubecoreg.rds")
saveRDS(glm_res, file = "output/glm_res.rds")

# plot
for(i in 1:length(levels(as.factor(bird_se$code_sp)))){
  sp <- levels(as.factor(ts_bird_se_allcountry$code_sp))[i]
  gp <- ggplot(ts_bird_se_allcountry[ts_bird_se_allcountry$code_sp==sp,],
         aes(year, relative_abundance)) + geom_line() + geom_text(x=2010, y=1, label=sp) +
    geom_line(aes(y=CI_inf), linetype="dashed") +
    geom_line(aes(y=CI_sup), linetype="dashed") +
    ylab("Relative abundance") + xlab("Years") +
    theme_modern()
  print(gp)
}

