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

# calculate MSI
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
  # span = 0.75 is default in R and usually fits well. But trying other values may give better results, depending on your data
  loessMSI <- array(NA, dim=c(nyear, nsim))
  Diff <- array(NA, dim=c(nyear, nsim))
  # Difftest <- array(NA, dim=c(nyear, nsim))
  for (s in 1:nsim) {
    smooth <- predict(loess(simMSI[,s]~uyear, span=0.75, degree=2, na.action=na.exclude),
                      data.frame(uyear), se=TRUE)
    loessMSI[,s] <- round(smooth$fit, digits=4)
    
    for (y in 1:nyear) {
      Diff[y,s] <- loessMSI[nyear,s] - loessMSI[y,s]
      #    Difftest[y,s] <- loessMSI[testyear,s] - loessMSI[y,s]
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





# cla trajectory
class.trajectory <- function (Y = NULL, X = NULL, dataset = NULL, interval_size = 0.5)
{
  if (is.null(Y) == TRUE & is.null(Y) == TRUE & is.null(dataset) == TRUE){
    stop("either 'dataset' or at least 'Y' and 'X' must be specified")
  }
  if (is.null(Y) == TRUE & is.null(Y) == TRUE) {
    Y <- dataset[,1]
    X <- dataset[,2]
  }else{
    if (class(Y) == "character" & class(X) == "character") {
      if (is.null(dataset) == TRUE) {
        stop("if 'Y' and 'X' are character, 'dataset' must exist")
      }else{
        Y <- dataset[, Y]
        X <- dataset[, X]
      }
    }else{
      if (!(class(Y) %in% c("numeric","integer")) == TRUE & !(class(X) %in% c("numeric","integer")) == TRUE) {stop("'Y' and 'X' must be either characters or vector but 'class' must be similar")}
    }
  }
  
  data <- data.frame(cbind(Y, X))
  data <- data[order(data$X),]                                                                      # ordering the X values
  
  if (length(X)<4){
    stop("time series length must be at least 4")
  }
  
  Y <- data$Y
  X <- data$X
  
  linear.model <- lm(Y~X)
  
  orthogonal_polynomial <- lm(Y~poly(X,2, raw=F))                                                   # After getting Y = gamma*chi + delta*X' + epsilon with orthogonal polynomial
  # we have to perform a variable change to obtain relevant values in the X interval 
  # for first_order_coefficient, second_order_coefficient and intercept,
  # knowing that X'= alpha*X + beta 
  # and chi = eta*X'^2 + theta
  
  gammab  <-  orthogonal_polynomial$coefficients[3]
  delta  <-  orthogonal_polynomial$coefficients[2]
  epsilon  <-  orthogonal_polynomial$coefficients[1]
  
  alpha  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[2]
  beta  <-  lm(orthogonal_polynomial$model[, 2][, 1]~X)$coef[1]
  
  eta  <-  1/lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[2]
  theta  <-  (-lm((orthogonal_polynomial$model[, 2][, 1])^2~orthogonal_polynomial$model[, 2][, 2])$coef[1])*eta
  
  Y2<-Y*(max(X)-min(X))/(max(Y)-min(Y))                                                             # p2 and p3 are relevant when Y and X amplitudes are equivalent,
  # in particular when studying scaled-to-1 indices, Y and X amplitudes
  # may be very different, so we scaled the amplitudes to calculate p2 and p3 
  polynomial_orthonormal_basis<-lm(Y2~poly(X,2, raw=T))$coefficients
  
  if(summary(orthogonal_polynomial)$coefficients[3, 4] <= 0.05){                                     # non linear case
    classification <- data.frame(first_order_coefficient = (delta+2*beta*gammab*eta)*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = (alpha^2)*gammab*eta,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+beta*delta+(beta^2)*gammab*eta+gammab*theta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = -(delta+2*beta*gammab*eta)/(2*alpha*gammab*eta),                    # points of interest
                                 p2 = (-polynomial_orthonormal_basis[2]+1)/(2*polynomial_orthonormal_basis[3]),
                                 p3 = (-polynomial_orthonormal_basis[2]-1)/(2*polynomial_orthonormal_basis[3]))
  }else{                                                                                            # linear case
    classification <- data.frame(first_order_coefficient = delta*alpha,
                                 first_order_pvalue = summary(orthogonal_polynomial)$coefficients[2, 4],
                                 second_order_coefficient = 0,
                                 second_order_pvalue = summary(orthogonal_polynomial)$coefficients[3, 4],
                                 strd_error=summary(orthogonal_polynomial)$coefficients[2, 2],
                                 intercept = epsilon+delta*beta,
                                 x_m = (X[length(X)]-X[1])/2+X[1],
                                 p1 = NA,
                                 p2 = NA,
                                 p3 = NA)
  }
  
  classification$r.sq <- summary(orthogonal_polynomial)$adj.r.squared                                # retrieve the adjusted coefficient of determination
  
  # compute the derivaive at xm-delta and at xm + delta with delta being half of the input interval size
  derivative  <-  2*(classification$x_m-(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  derivative2  <-  2*(classification$x_m+(X[length(X)]-X[1])*(interval_size/2))*classification$second_order_coefficient+classification$first_order_coefficient
  
  
  if(sign(derivative) != sign(derivative2)){                                                        # non consistent direction around x_m
    classification$derivative  <-  NA
    classification$intercept_derivative  <-  NA
  }else{                                                                                            # consistent direction around x_m
    classification$derivative  <-  mean(c(derivative, derivative2))
    classification$intercept_derivative  <-  (classification$second_order_coefficient*classification$x_m^2+classification$first_order_coefficient*classification$x_m+classification$intercept)-classification$x_m*classification$derivative
  }
  
  # compute the derivative of the curvature function
  classification$derivated_curvature  <-  -12*(classification$second_order_coefficient^2)*(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)*(classification$second_order_coefficient/abs(classification$second_order_coefficient))/
    ((1+(2*classification$second_order_coefficient*classification$x_m+classification$first_order_coefficient)^2)^(2.5))
  
  if(classification$second_order_pvalue>0.05){classification$derivated_curvature <- NA}
  
  classification$direction <- NA                                                                    # classify the direction
  classification$direction[which(classification$derivative > 0)] <- "increase"
  classification$direction[which(classification$derivative < 0)] <- "decrease"
  classification$direction[which(is.na(classification$derivative))] <- "stable"
  classification$direction[which(as.numeric(classification$first_order_pvalue)>0.05 & as.numeric(classification$second_order_pvalue)>0.05)] <- "stable"
  
  classification$acceleration <- NA                                                                 # classify the acceleration
  classification$acceleration[which(classification$derivated_curvature < 0)] <- "accelerated"
  classification$acceleration[which(classification$derivated_curvature > 0)] <- "decelerated"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient < 0)] <- "concave"
  classification$acceleration[which(classification$direction == "stable" &
                                      classification$second_order_coefficient > 0)] <- "convex"
  classification$acceleration[which(is.na(classification$derivated_curvature))] <- "constant"
  
  classification$shape_class <- paste(classification$direction,                                       # give the final classification combining direction and acceleration
                                      classification$acceleration,
                                      sep="_")
  
  linear.model.summary <- summary(linear.model)                                                       # provide the linear approach results for comparison
  
  classification$linear_slope <- linear.model.summary$coefficients[2, 1]
  classification$linear_slope_pvalue <- linear.model.summary$coefficients[2, 4]
  classification$linear_intercept <- linear.model.summary$coefficients[1, 1]
  
  classification$first_X_value <- X[1]
  classification$last_X_value <- X[length(X)]
  
  row.names(classification) <- "Y"
  
  return(classification)
  
}



















