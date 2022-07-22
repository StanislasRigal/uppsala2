# function for computing bird time series

get_ts <- function(data_bird_input){
  
  # d data for species i
  
  d <- droplevels(data_bird_input)
  
  sp <- levels(as.factor(d$code_sp))
  
  # number of route followed by year
  
  nb_route <- tapply(rep(1,nrow(d)),d$year,sum)
  
  # number of route with species i by year
  
  nb_route_presence <- tapply(ifelse(d$abund>0,1,0),d$year,sum)
  
  year<-as.numeric(as.character(levels(as.factor(d$year))))
  
  firstY <- min(year)
  
  lastY <- max(year)
  
  timestep <- length(year)-1
  
  # table for analysis result
  
  threshold_occurrence <- 3
  
  tab_ana <- data.frame(year=rep(year,2),val=c(nb_route,nb_route_presence),LL = NA,UL=NA,
                        catPoint=NA,pval=NA,
                        curve=rep(c("route","presence"),each=length(year)))
  
  tab_ana$catPoint <- ifelse(tab_ana$val == 0,"0", ifelse(tab_ana$val < threshold_occurrence,
                                                          "inf_threshold",NA))
  
  
  # abundance by year
  
  abund <- tapply(d$abund,d$year,sum)
  
  threshold_abundance <- 5
  
  tab_fig <- data.frame(year=year,val=abund,LL = NA,UL=NA,catPoint=NA,pval=NA)
  
  tab_fig$catPoint <- ifelse(tab_fig$val == 0,"0",ifelse(tab_fig$val < threshold_abundance,
                                                         "inf_threshold",NA))
  
  # remove criteria
  
  remove_sp <- FALSE
  
  # if first year empty
  
  if(tab_fig$val[1]==0){remove_sp <- TRUE}
  
  # if four consecutive years empty
  
  ab_vec <- paste(tab_fig$val,collapse="")
  
  if(str_detect(ab_vec, "0000")){remove_sp <- TRUE}
  
  # if less than three consecutive years
  
  ab_vec2 <- paste(sign(tab_fig$val),collapse="")
  
  if(!str_detect(ab_vec2, "111")){remove_sp <- TRUE}
  
  if(anyNA(tab_fig$catPoint) & anyNA(tab_ana$catPoint[tab_ana$curve=="presence"]) & remove_sp==F){
    
    # GLM abundance variation
    glm1 <- glm(abund~as.factor(code_route)+as.factor(year),data=d,family=quasipoisson)
    
    sglm1 <- summary(glm1)
    
    # mean-centered values
    
    con.mat <- diag(length(year)) - 1/length(year)
    
    colnames(con.mat) <- year # firstY:lastY
    
    rg <- ref_grid(glm1, nuisance = 'code_route')
    
    sglm2 <- summary(contrast(rg, as.data.frame(con.mat)))
    
    # as link function is log, estimates need to be back transformed from sglm1 (first year set to 1 and se to 0)
    
    coef_yr <- tail(matrix(sglm1$coefficients[,1]), timestep)
    
    coef_yr <- rbind(1, exp(coef_yr))
    
    error_yr <- tail(matrix(sglm1$coefficients[,2]), timestep)
    
    error_yr <- rbind(0, error_yr)*coef_yr # approximated se values
    
    log_error_yr <- tail(matrix(sglm1$coefficients[,2]), timestep)
    
    log_error_yr <- rbind(0, log_error_yr)
    
    pval <- c(1,tail(matrix(coefficients(sglm1)[,4]),timestep))
    
    # from sglm2 (mean value to 0)
    
    coef_yr_m0 <- exp(sglm2$estimate)
    
    error_yr_m0 <- sglm2$SE*coef_yr_m0 # approximated se values
    
    log_error_yr_m0 <- sglm2$SE
    
    pval_m0 <- sglm2$p.value
    
    # CIs
    
    glm1.sim <- sim(glm1)
    
    ci_inf_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .025), timestep)))
    
    ci_sup_sim <- c(1, exp(tail(apply(coef(glm1.sim),2, quantile, .975), timestep)))
    
    thresold_signif <- 0.05
    
    tab_res <- data.frame(year,val=coef_yr,val_m0=coef_yr_m0,
                          LL=ci_inf_sim,UL=ci_sup_sim,
                          catPoint=ifelse(pval<thresold_signif,"significatif",NA),pval)
    
    # cleaning out of range CIs			   
    
    tab_res$UL <- ifelse(nb_route_presence==0,NA,tab_res$UL)
    
    tab_res$UL <-  ifelse(tab_res$UL == Inf, NA,tab_res$UL)
    
    tab_res$UL <-  ifelse(tab_res$UL > 1.000000e+20, NA,tab_res$UL)
    
    tab_res$UL[1] <- 1
    
    tab_res$val <-  ifelse(tab_res$val > 1.000000e+20,1.000000e+20,tab_res$val)
    
    tab_res$val_m0 <-  ifelse(tab_res$val_m0 > 1.000000e+20,1.000000e+20,tab_res$val_m0)
    
    # overdispersion index
    dispAn <- sglm1$deviance/sglm1$null.deviance
    
    # class uncertainity
    
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
    
    # Store results      
    
    tab_tot <- data.frame(code_sp=sp, year = tab_res$year, nb_year=timestep,
                          firstY = firstY, lastY = lastY,
                          relative_abundance = tab_res$val,
                          CI_inf = tab_res$LL, CI_sup = tab_res$UL,
                          Standard_error = error_yr,
                          Log_SE = log_error_yr,
                          p_value = tab_res$pval,
                          relative_abundance_m0 = tab_res$val_m0,
                          Standard_error_m0 = error_yr_m0,
                          Log_SE_m0 = log_error_yr_m0,
                          p_value_m0 = pval_m0, signif = !is.na(tab_res$catPoint),
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
                          p_value = NA, 
                          relative_abundance_m0 = NA,
                          Standard_error_m0 = NA,
                          Log_SE_m0 = NA,
                          p_value_m0 = NA,signif = NA,
                          nb_route,nb_route_presence,abundance=abund,
                          mediane_occurrence=median(nb_route_presence), mediane_ab=median(abund) ,
                          valid = NA, uncertanity_reason = NA)
  }
  
  return(tab_tot)
}

# simulate ts from mean and se

mc_sim <- function(x, niter){
  
  sim_ts <- data.frame(year=x$year)
  
  for(i in 1:niter){
  
      a <- rnorm(nrow(x), mean=x$relative_abundance, sd=x$se)
    
      sim_ts <- as.data.frame(cbind(sim_ts,a))
    
      names(sim_ts)[i+1] <- paste0("sim",i)
  
    }
  
  return(sim_ts)
}

# simulate ts for a set of species

get_sim <- function(dataset, species, niter){
  
  d <- droplevels(na.omit(dataset[dataset$code_sp %in% species,]))
  
  sp <- levels(as.factor(d$code_sp))
  
  d2 <- c()
  
  for(j in 1: length(sp)){
  
      d2 <- rbind(d2,d[d$code_sp==sp[j],c("code_sp","relative_abundance","year","CI_inf")])
    
  
    }
  
  d2$se <- (d2$relative_abundance-d2$CI_inf)/1.96
  
  d_sim <- ddply(d2, .(code_sp), .fun=mc_sim, niter=niter)
  
  d_sim_long <- melt(d_sim, id.vars = c("code_sp", "year"))
  
  d_sim_long$id <- paste0(d_sim_long$code_sp, sep="_", d_sim_long$variable)
  
  d_sim_wide <- dcast(d_sim_long, id~year, value.var = "value", fun.aggregate = sum)
  
  return(d_sim_wide)
}

# format data for plotting after dfa

data.dfa.plot <- function(dataset, species, sim_data, dfa_res, nsim){
  
  d <- droplevels(na.omit(dataset[dataset$code_sp %in% species,]))
  
  sp <- levels(as.factor(d$code_sp))
  
  d2 <- c()
  
  for(j in 1: length(sp)){
  
      d2 <- rbind(d2,d[d$code_sp==sp[j],c("code_sp","relative_abundance","year","CI_inf")])
  
    }
  
  d2$se <- (d2$relative_abundance-d2$CI_inf)/1.96
  
  d_sim <- ddply(d2, .(code_sp),
                 .fun=function(dataset2,species2, sim_data, dfa_res, nsim){
                   mean_rel_ab <- mean(dataset2$relative_abundance,na.rm=T)
                   sd_rel_ab <- sd(dataset2$relative_abundance,na.rm=T)
                   dataset2$relative_abundance_std <- (dataset2$relative_abundance-mean_rel_ab)/sd_rel_ab
                   dataset2$se_std <- (dataset2$se)/sd_rel_ab
                   id_sp <- levels(as.factor(dataset2$code_sp))
                   nb_sp <-  which(species2 == levels(as.factor(dataset2$code_sp)))
                   dataset2 <- cbind(dataset2,t(sim_data[grepl(id_sp,sim_data$id),-1]))
                   mean_ts_dfa <- apply((dfa_res$Estimates$Z[((nb_sp-1)*nsim+1):(nb_sp*nsim),] %*% dfa_res$Estimates$u),2,mean)
                   sd_ts_dfa <- apply((dfa_res$Estimates$Z[((nb_sp-1)*nsim+1):(nb_sp*nsim),] %*% dfa_res$Estimates$u),2,sd)
                   mean_mean_ts_dfa <- mean(mean_ts_dfa,na.rm=T)
                   sd_mean_ts_dfa <- sd(mean_ts_dfa,na.rm=T)
                   dataset2$mean_ts_dfa_std <- (mean_ts_dfa-mean_mean_ts_dfa)/sd_mean_ts_dfa
                   dataset2$sd_ts_dfa_std <- (sd_ts_dfa)/sd_mean_ts_dfa
                   return(dataset2)
                 }, species2=sp, sim_data=sim_data, dfa_res=dfa_res, nsim=nsim)
  
  return(d_sim)
}

# calculate MSI (from https://doi.org/10.1016/j.ecolind.2017.05.033)

MSI_MC_func <- function(dataset,
                        nsim=1000, # number of Monte Carlo simulations
                        SEbaseyear=1996, # year to set SE to 0. SE of other years is relative to the SE of SEbaseyear
                        plotbaseyear=1996, # year to set MSI or smoothed trend (dependend on 'index_smooth') to 100 in final MSI-plot
                        index_smooth="SMOOTH", # INDEX / SMOOTH: "INDEX" will cause the MSI in plotbaseyear set to 100; "SMOOTH" will set the smoothed trend value in the plotbaseyear to 100
                        maxCV=3,       # maximum allowed mean Coefficient of Variation (CV) of species indices (0.5 = 50%). Species with higher mean CV are excluded
                        truncfac=10,    # truncation factor (=maximum allowed year-to-year index ratio). Default for Living Planet Index = 10
                        TRUNC=1 # set all indices below TRUNC to this value and their SE to 0. TRUNC <- 0 means no truncation
){
  
  # Define and calculate model parameters
  
  dataset <- dataset[!is.na(dataset$code_sp) & !is.na(dataset$year) &
                       !is.na(dataset$relative_abundance) & !is.na(dataset$CI_inf),
  ]
  
  species <- dataset$code_sp
  
  year <- dataset$year
  
  index <- dataset$relative_abundance*100
  
  se <- 100*(dataset$relative_abundance-dataset$CI_inf)/1.96
  
  uspecies <- sort(unique(species))
  
  nspecies <- length(uspecies)
  
  uyear <- sort(unique(year))
  
  nyear <- length(unique(year))
  
  meanindex <- tapply(index, species, mean, na.rm=TRUE)/100
  
  mnindex1 <- as.data.frame(rep(meanindex,each=nyear))
  
  minyear <- min(year)
  
  maxyear <- max(year)
  
  plotbaseyear <- plotbaseyear-minyear+1
  
  baseyear <- max(1,SEbaseyear-minyear+1)
  
  INP <- data.frame(cbind(species, year, index, se))
  
  INP$index <- as.numeric(INP$index)
  
  INP$se <- as.numeric(INP$se)
  
  SPEC <- as.matrix(species)
  
  species <- sort(rep(uspecies,each=nyear))
  
  year <- rep(uyear,nspecies)
  
  INP1 <- data.frame(cbind(species,year))
  
  INP2 <- merge(INP, INP1, by=c("species","year"), sort=TRUE, all=TRUE)
  
  # Calculate and plot mean CV for indices per species
  
  CVtemp <- INP2[INP2$index >= 10, ] # select records with index >= 10 for CV calculation
  
  CVtemp <- CVtemp[!is.na(CVtemp$index), ] # reject missing values
  
  CV1 <- CVtemp$se/CVtemp$index
  
  CV1[CV1== 0] <- NA
  
  CV1[CV1== Inf] <- NA
  
  species2 <- CVtemp$species
  
  mnCV <- tapply(CV1, species2, mean, na.rm=TRUE)
  
  CV <- as.data.frame(rep(mnCV, each=nyear))
  
  # replace small indices by TRUNC and their SE by zero
  
  INP3 <- data.frame(cbind(species, year, CV, mnindex1))
  
  colnames(INP3) <- c("species","year", "CV", "mnindex1")
  
  INP4 <- merge(INP2, INP3, by=c("species","year"), sort=FALSE, all=TRUE)
  
  INP5 <- subset(INP4, CV < maxCV, select = c(species, year, index, se, mnindex1))
  
  INP5$index[INP5$index < TRUNC & !is.na(INP5$index)] <- TRUNC 
  
  INP5$se[INP5$index == TRUNC & !is.na(INP5$index)] <- 0
  
  # reset parameters
  
  nobs <- NROW(INP5)
  
  uspecies <- sort(unique(INP5$species))
  
  nspecies <- length(uspecies)
  
  year <- rep(uyear, nspecies)
  
  # Transform indices and standard deviations to log scale (Delta method)
  
  LNindex <- as.vector(log(index))
  
  LNse <- as.vector(se/index)
  
  # Monte Carlo simulations of species indices
  
  MC <- matrix(NA, nobs, nsim)
  
  for (s in 1:nsim) { 
  
      MC[,s] <- rnorm(nobs, LNindex, LNse)
  
  } 
  
  MC[MC < log(TRUNC)] <- log(TRUNC)
  
  # impute missing values using chain method
  
  CHAIN <- matrix(NA, nobs, nsim)
  
  for (s in 1:nsim ){
  
      for (o in 1:nobs-1) {
    
          CHAIN[o,s] <- MC[o+1,s]-MC[o,s]
    
          }
    
    for (sp in 1:nspecies) {
    
        CHAIN[sp*nyear,s] <- NA
    
    }
  }
  
  CHAIN[CHAIN > log(truncfac)] <- log(truncfac)
  
  CHAIN[CHAIN < log(1/truncfac)] <- log(1/truncfac)
  
  mnCHAIN <- matrix(NA, nyear, nsim)
  
  for (s in 1:nsim) {
  
      mnCHAIN[,s] <- tapply(CHAIN[,s], year, mean, na.rm=TRUE)
  
    }
  
  simMSI <- matrix(NA, nyear, nsim)
  
  bmin <- baseyear-1
  
  bplus <- baseyear+1
  
  for (s in 1:nsim){
  
      simMSI[baseyear,s] <- log(100)
  
        for (y in bmin:min(bmin,1)){
      
          simMSI[y,s] <- simMSI[y+1,s] - mnCHAIN[y,s]
    
          }
    
      for (y in min(bplus,nyear):nyear){
      
        simMSI[y,s] <- simMSI[y-1,s] + mnCHAIN[y-1,s]
    
        }
  }
  
  # calculate MSI and SE
  
  meanMSI <- array(NA, dim=c(1, nyear))
  
  stdevMSI <- array(NA, dim=c(1, nyear))   
  
  for (y in 1:nyear) {
  
      meanMSI[y] <- mean(simMSI[y,])
    
      stdevMSI[y] <- sd(simMSI[y,])
  
  }
  
  # Monte Carlo simulations of MSI (for trend calculation)
  
  simMSI <- matrix(NA, nyear, nsim)
  
  for (s in 1:nsim) { 
  
      for (y in 1:nyear) {
    
          simMSI[y,s] <- rnorm(1, meanMSI[y], stdevMSI[y])
    
          } 
  }
  
  # Back-transformation of MSI to index scale (Delta method)
  
  meanMSI <- array(NA, dim=c(1, nyear))
  
  stdevMSI <- array(NA, dim=c(1, nyear))   
  
  for (y in 1:nyear) {
  
      meanMSI[y] <- round(exp(mean(simMSI[y,])), digits=2)
    
      stdevMSI[y] <- round(sd(simMSI[y,])*meanMSI[y], digits=2)
  
    }
  
  # Confidence interval of MSI on index-scale
  
  CI <- array(NA, dim=c(nyear, 2))
  
  for (y in 1:nyear){
  
      CI[y,1] <- exp(mean(simMSI[y,])-1.96*sd(simMSI[y,]))
    
      CI[y,2] <- exp(mean(simMSI[y,])+1.96*sd(simMSI[y,]))
  }
  
  # Smoothing

  loessMSI <- array(NA, dim=c(nyear, nsim))
  
  Diff <- array(NA, dim=c(nyear, nsim))
  
  for (s in 1:nsim) {
  
      smooth <- predict(loess(simMSI[,s]~uyear, span=0.75, degree=2, na.action=na.exclude),
                      data.frame(uyear), se=TRUE)
    
      loessMSI[,s] <- round(smooth$fit, digits=4)
    
    
      for (y in 1:nyear) {
      
        Diff[y,s] <- loessMSI[nyear,s] - loessMSI[y,s]
        
    } 
  }
  
  
  # create output for flexible trend estimates
  smoothMSI <- array(NA, dim=c(nyear, 13))
  for (y in 1:nyear) {
    smoothMSI[y,1] <- round(mean(loessMSI[y,]), digits=4) # smooth MSI on log scale
    smoothMSI[y,2] <- round(sd(loessMSI[y,]), digits=4) # SE of smooth MSI on log scale
    smoothMSI[y,3] <- round(mean(Diff[y,]), digits=4) # Difference smooth MSI with last year on log scale
    smoothMSI[y,4] <- round(sd(Diff[y,]), digits=4) # SE of difference with ast year
    smoothMSI[y,12] <- round(smoothMSI[y,1]-1.96*smoothMSI[y,2], digits=2) # lower CI of smooth MSI on log scale
    smoothMSI[y,13] <- round(smoothMSI[y,1]+1.96*smoothMSI[y,2], digits=2) # upper CI of smooth MSI on log scale
  }
  
  # Trendclassification, based on smoothed trends (Soldaat et al. 2007 (J. Ornithol. DOI 10.1007/s10336-007-0176-7))
  for (y in 1:nyear) {
    smoothMSI[y,5] <- round(smoothMSI[nyear,1]/smoothMSI[y,1], digits=4)		# = TCR
    smoothMSI[y,6] <- round(smoothMSI[y,5]-1.96*(smoothMSI[y,4]/smoothMSI[y,1]), digits=4)	# = CI- TCR
    smoothMSI[y,7] <- round(smoothMSI[y,5]+1.96*(smoothMSI[y,4]/smoothMSI[y,1]), digits=4)	# = CI+ TCR
    smoothMSI[y,8] <- round(exp(log(smoothMSI[y,5])/(nyear-y)), digits=4)		# = YCR
    smoothMSI[y,9] <- round(exp(log(smoothMSI[y,6])/(nyear-y)), digits=4)		# = CI- YCR
    smoothMSI[y,10] <- round(exp(log(smoothMSI[y,7])/(nyear-y)), digits=4)		# = CI+ YCR
  }
  for (y in 1:(nyear-1)) {
    if (smoothMSI[y,9] > 1.05) {smoothMSI[y,11] <- 1} else
      if (smoothMSI[y,10] < 0.95) {smoothMSI[y,11] <- 6} else
        if (smoothMSI[y,9] > 1.00) {smoothMSI[y,11] <- 2} else
          if (smoothMSI[y,10] < 1.00) {smoothMSI[y,11] <- 5} else
            if ((smoothMSI[y,9] - 0.95)*(1.05 - smoothMSI[y,10]) < 0.00) {smoothMSI[y,11] <- 3} else
              if ((smoothMSI[y,10]) - (smoothMSI[y,9]) > 0.10) {smoothMSI[y,11] <- 3} else
              {smoothMSI[y,11] <- 4}
  } 
  TrendClass_flex <- matrix(NA, nrow= nyear, ncol=1)
  for (y in 1:(nyear-1)) {
    if (smoothMSI[y,9] > 1.05) {TrendClass_flex[y] <- "strong_increase" } else
      if (smoothMSI[y,10] < 0.95) { TrendClass_flex[y] <- "steep_decline" } else
        if (smoothMSI[y,9] > 1.00) { TrendClass_flex[y] <- "moderate_increase"} else
          if (smoothMSI[y,10] < 1.00) { TrendClass_flex[y] <- "moderate_decline" } else
            if ((smoothMSI[y,9] - 0.95)*(1.05 - smoothMSI[y,10]) < 0.00) { TrendClass_flex[y] <- "uncertain"} else
              if ((smoothMSI[y,10]) - (smoothMSI[y,9]) > 0.10) { TrendClass_flex[y] <- "uncertain" } else
              {TrendClass_flex[y] <- "stable"}
  } 
  
  # rescaling (for presentation of plot)
  
  rescale <- NA
  
  if (index_smooth =="INDEX") {rescale <- 100/meanMSI[plotbaseyear]} else
  
      if (index_smooth =="SMOOTH") {rescale <- 100/exp(smoothMSI[plotbaseyear,1])} else
    
        {rescale <- NA}
  
  simMSImean <- round(as.vector(rescale*meanMSI), digits=2)
  
  simMSIsd <- round(as.vector(rescale*stdevMSI), digits=2)
  
  uppCI_MSI <- round(rescale*CI[,2], digits=2)
  
  lowCI_MSI <- round(rescale*CI[,1], digits=2)
  
  trend_flex <- round(rescale*exp(smoothMSI[,1]), digits=2)
  
  lowCI_trend_flex <- round(rescale*exp(smoothMSI[,12]), digits=2)
  
  uppCI_trend_flex <- round(rescale* exp(smoothMSI[,13]), digits=2)
  

  # create output file for MSI + smoothed trend
  
  RES <- as.data.frame(cbind(uyear, simMSImean, simMSIsd, lowCI_MSI, uppCI_MSI, trend_flex, lowCI_trend_flex, uppCI_trend_flex))
  
  RES$trend_class <- TrendClass_flex
  
  names(RES) <- c("year", "MSI", "sd_MSI", "lower_CL_MSI", "upper_CL_MSI", "Trend", "lower_CL_trend", "upper_CL_trend", "trend_class")
  
  return(RES)
}

