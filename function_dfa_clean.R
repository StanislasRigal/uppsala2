# I) TMB function

dfa_model_se <- "

// Dynamic Factor Analysis for multivariate time series
#include <TMB.hpp>
template<class Type>
  Type objective_function<Type>::operator() ()
{
  // Data
  DATA_ARRAY(y);
  // Observation standard error
   DATA_ARRAY(obs_se);
  
  int nSp = y.dim[0];
  int nT = y.dim[1];
  
  // For one-step-ahead residuals
  DATA_ARRAY_INDICATOR(keep, y);
  
  DATA_MATRIX(Z_pred);
  DATA_UPDATE(Z_pred);
  
  DATA_MATRIX(W); // Weighting matrix to compute mean trends over species.
  DATA_UPDATE(W);
  
  // Parameters
  PARAMETER_VECTOR(log_re_sp); // log of sd for random effect by species
  
  // Loadings matrix
  PARAMETER_MATRIX(Z);
  
  // Latent trends
  PARAMETER_MATRIX(x);
  
  // Cluster center
  matrix<Type> x_pred(Z_pred.rows(), nT);
  matrix<Type> x_pred2(W.rows(), nT);
  matrix<Type> WZ(W.rows(), Z.cols());

  // Mean of latent trends
  matrix<Type> x_sum(x.rows(), 1);
  
  // Matrix to hold predicted species trends
  matrix<Type> x_sp(nSp, nT);
  
  // Random error by species
  vector<Type> re_sp;
  re_sp = exp(log_re_sp);
  // Optimization target: negative log-likelihood (nll)
  Type nll = 0.0;
  
  
    // Latent random walk model. x(0) = 0. 
    for(int t = 1; t < x.cols(); ++t){
      for(int f = 0; f < x.rows(); ++f){
        nll -= dnorm(x(f, t), x(f, t-1), Type(1), true);
    
      // Simulation block for process equation
      SIMULATE {
          x(f, t) = rnorm(x(f, t-1), Type(1));
          REPORT(x);
      }
      }
    }
  
  for (int f = 0; f < x.rows(); ++f) {
    x_sum(f) = x.row(f).sum();
    SIMULATE {
      x_sum(f) = x.row(f).sum();
    }
  }

  // Species trends
  for(int i = 0; i < nSp; ++i) {
    x_sp(i, 0) = (Z.row(i) * (-x_sum)).sum();
    for(int t = 1; t < nT; ++t) {
      x_sp(i, t) = (Z.row(i) * (x.col(t))).sum();
    }
  }  
  
  // Cluster center
  WZ = W * Z;
  x_pred.col(0) = Z_pred * (-x_sum);
  x_pred2.col(0) = WZ * (-x_sum);
  for (int t=1; t < nT; ++t) {
    x_pred.col(t) = Z_pred * (x.col(t));
    x_pred2.col(t) = WZ * (x.col(t));
  }  
  
  
  // Observation model
  for(int i = 0; i < nSp; ++i){
  // Skipping t = 0 when y(i, 0) is fixed at 0. Need to change this if y(i, 0) is not 0.
  // Also had had to change the index of x from t+1 to t, so that x is fixed at zero at time t=0.
    for(int t = 0; t < nT; ++t){
    if(!R_IsNA(asDouble(y(i,t)))){
    nll -= keep(i,t) * dnorm(y(i, t), x_sp(i, t), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)), true); // with random effect
    }
      
      
      //*----------------------- SECTION I --------------------------*/
        // Simulation block for observation equation
      SIMULATE {
          y(i,t) = rnorm(x_sp(i, t), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)));
          REPORT(y);
        }
    }  
  }
  
  //Penalty for very small random effects variances, this is to avoid issues with non semidefinite variance matrices.
  for (int i=0; i < nSp; ++i) {
    nll += exp(-3.5 * log_re_sp(i) - 20);
  }
  
  // State the transformed parameters to report
  // Using ADREPORT will return the point values and the standard errors
  // Note that we only need to specify this for parameters
  // we transformed, see section D above
  // The other parameters, including the random effects (states),Q
  // will be returned automatically
  ADREPORT(re_sp);
  ADREPORT(x_sp);
  ADREPORT(x_pred);
  ADREPORT(x_pred2);
  ADREPORT(Z_pred);
  
  // Report simulated values
  //SIMULATE{
    //  REPORT(x);
    //  REPORT(y);
    //}
  
  return nll;
}"
  
## Compile .cpp file with TMB DFA model
  
  if(!exists('dfamodel')){
    write(dfa_model_se, file = "dfa_model_se.cpp")
    dfamodel <- dfa_model_se
    compile("dfa_model_se.cpp")
    dyn.load(dynlib("dfa_model_se"))
  }else if(dfa_model_se != dfamodel){
    write(dfa_model_se, file = "dfa_model_se.cpp")
    dfamodel <- dfa_model_se
    compile("dfa_model_se.cpp")
    dyn.load(dynlib("dfa_model_se"))
  }
  


# II) Functions used in the main function make_dfa  

## 1) Function to compute AIC of DFA models (called inside make_dfa)
## Only works if obj has been already optimized
## AIC is computed excluding zero variance

AIC.tmb <- function(obj) {
  
  # AIC
  
  as.numeric(2 * obj$env$value.best + 2*(sum(obj$env$lfixed())))
}


## 2) Group species

group_from_dfa_boot1 <- function(data_loadings, # Species initial factor loadings
                                cov_mat_Z, # Covariance matrix of species factor loadings
                                species_sub, # Species names
                                nboot=100, # Number of bootstrap iteration
                                ny, # Number of time series
                                nfac # Number of latent trends
){
  
  # Indices of fixed loadings
  
  row_col_0 <- which(data_loadings$value==0)
  
  # Loadings in matrix
  
  dfa_res_val <- dcast(data_loadings, code_sp~variable, value.var = "value")
  mat_loading <- as.matrix(dfa_res_val[,-1])
  
  # Find the best number of clusters in the original data
  
  NbClust2 <- function(data, diss=NULL, distance = "euclidean",
                       method = "kmeans", min.nc=2, max.nc=min(max(c(2,round(ny/3))),10), 
                       index = "alllong", alphaBeale = 0.1){
    tryCatch(
      #try to do this
      {
        NbClust(data, diss=NULL, distance = "euclidean",
                method = "kmeans", min.nc=2, max.nc=min(max(c(2,round(ny/3))),10), 
                index = "alllong", alphaBeale = 0.1)
      },
      #if an error occurs
      error=function(e) {
        data.frame(Best.partition=rep(NA,nrow(mat_loading)))
      }
    )
  }
  
  if(min(max(c(2,round(ny/3))),10)>4){
    nb <- NbClust2(mat_loading, diss=NULL, distance = "euclidean",
                   method = "kmeans", min.nc=2, max.nc=min(max(c(2,round(ny/3))),10), 
                   index = "alllong", alphaBeale = 0.1)
    nb_group_best <- max(nb$Best.partition)
  }else{
    idx <- 0
    nb_group_best_part <- c()
    for(index in c("kl", "ch", "hartigan", "ccc", "scott", "marriot", "trcovw", "tracew",
                   "friedman", "rubin", "cindex", "db", "silhouette", "duda",
                   "ratkowsky", "ball", "ptbiserial", "gap", "mcclain", "gamma",
                   "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")){
      nb_part <- NbClust(mat_loading, diss=NULL, distance = "euclidean",
                    method = "kmeans", min.nc=2, max.nc=min(max(c(2,round(ny/3))),10), 
                    index = index, alphaBeale = 0.1)
      idx <- idx + 1
      nb_group_best_part[idx] <- max(nb_part$Best.partition)
    }
    which.max(table(nb_group_best_part))
    nb_group_best <- as.numeric(which.max(table(nb_group_best_part)))
  }
  


  # Check cluster stability

  stability_cluster_gr <- list()
  all_partition_gr <- list()
  
  # Bootstrap for cluster stability
  
  for(i in 1:nboot){
    
    # Draw factor loadings using covariance matrix
    set.seed(i)
    rand_load <- mvtnorm::rmvnorm(1, mean=data_loadings[!row_col_0,]$value, sigma=cov_mat_Z) 
    
    # Complete loading vector with fixed values
    for(j in 1:(nfac-1)){
      index_0 <- ny*j
      rand_load <- append(rand_load, rep(0,j), after=index_0)
    }
    
    rand_load <- matrix(rand_load, ncol=nfac, nrow=ny)
    
    # Find the best number of clusters in the bootstrap loadings
    
    nb <- NbClust2(rand_load, diss=NULL, distance = "euclidean",
                   method = "kmeans", min.nc=2, max.nc=min(max(c(2,round(ny/3))),10), 
                   index = "alllong", alphaBeale = 0.1)
    
    # Calculate reference partitions as a function of the number of cluster
    
    for(gr in nb_group_best:1){
      
      nb_group <- gr
      
      all_partition2_all <- kmeans(mat_loading, nb_group, iter.max = 1000, nstart = 10)
      all_partition2 <- all_partition2_all$cluster
      all_partition2_centers <- all_partition2_all$centers
      
      
      if(!anyNA(nb$Best.partition)){
        all_partition2 <- rbind(all_partition2,nb$Best.partition)
        
        # Compute Jaccard similarity to relabel clusters as in the original clustering
        
        jac_sim_res <- matrix(NA, ncol=length(unique(all_partition2[1,])),
                              nrow=length(unique(nb$Best.partition)))
        for(k in sort(unique(all_partition2[1,]))){
          for(l in sort(unique(nb$Best.partition))){
            jac_sim_mat <- all_partition2[c(1,(2)),]
            jac_sim_mat[1,][which(jac_sim_mat[1,]!=k)] <- 0
            jac_sim_mat[2,][which(jac_sim_mat[2,]!=l)] <- 0
            jac_sim_mat[jac_sim_mat>0] <- 1
            jac_sim <- c(1 - vegdist(jac_sim_mat, method="jaccard"))
            jac_sim_res[l,k] <- jac_sim
          }
        }
        
        # If same number of clusters
        
        if(length(unique(all_partition2[1,]))==length(unique(nb$Best.partition))){
          for(l in sort(unique(nb$Best.partition))){
            all_partition2[2,][which(nb$Best.partition==l)] <- which.max(jac_sim_res[l,])
          }
        }
        
        # If more clusters in the bootstrap clustering
        
        if(length(unique(all_partition2[1,]))<length(unique(nb$Best.partition))){
          l_data <- c()
          for(k in sort(unique(all_partition2[1,]))){
            l_data <- c(l_data,which.max(jac_sim_res[,k]))
          }
          k <- 0
          for(l in l_data){
            k <- k+1
            all_partition2[2,][which(nb$Best.partition==l)] <- k
          }
          extra_clus <- sort(unique(nb$Best.partition))[which(!(sort(unique(nb$Best.partition)) %in% l_data))]
          for(g_sup in 1:length(extra_clus)){
            k <- k +1
            all_partition2[2,][which(nb$Best.partition==extra_clus[g_sup])] <- k
          }
          
        }
        
        # If less clusters in the bootstrap clustering
        
        if(length(unique(all_partition2[1,]))>length(unique(nb$Best.partition))){
          k_data <- c()
          for(l in sort(unique(nb$Best.partition))){
            k_data <- c(k_data,which.max(jac_sim_res[l,]))
          }
          l <- 0
          for(k in k_data){
            l <- l+1
            all_partition2[2,][which(nb$Best.partition==l)] <- k
          }
        }
        
        if(i == 1){
          stability_cluster <- apply(jac_sim_res,2,max)
          stability_cluster_gr[[gr]] <- stability_cluster
          all_partition_gr[[gr]] <- all_partition2
        }else{
          stability_cluster <- stability_cluster_gr[[gr]]
          stability_cluster <- rbind(stability_cluster,apply(jac_sim_res,2,max))
          stability_cluster_gr[[gr]] <- stability_cluster
          all_partition <- all_partition_gr[[gr]]
          all_partition <- rbind(all_partition,all_partition2[2,])
          all_partition_gr[[gr]] <- all_partition
        }
      }else{if(i == 1){
        stability_cluster_gr[[gr]]<-1
        all_partition_gr[[gr]] <- rep(1,nrow(species_sub))
      }}
    }
  }
  
  stability_cluster_final <- 0
  
  gr <- nb_group_best + 1
  
  while(min(stability_cluster_final)<0.5 & gr > 1){  # dilution cluster if stability < 0.5
    
    gr <- gr - 1
    
    stability_cluster_final <- apply(stability_cluster_gr[[gr]],2,mean)
    
    }
  
  nb_group <- gr
  
  # Frequency of classification of each species in its reference cluster
  
  all_partition2 <- all_partition_gr[[nb_group]]
  all_partition2 <- na.omit(all_partition2)
  
  all_partition_uncertainty <- apply(all_partition2, 2,
                                     FUN = function(x){
                                       y <- max(table(x))/length(x)
                                       return(y)
                                     })

  all_partition_group <- all_partition2[1,]
  num_row <- 1
  while(length(unique(all_partition_group))!=nb_group){
    num_row <- num_row + 1
    all_partition_group <- all_partition2[num_row,]
    
  }
  
  
  # Compute PCA to get axes of the graph
  
  myPCA <- prcomp(mat_loading, scale. = F, center = T)
  
  # Group all info as output
  
  kmeans_1 <- merge(data.frame(code_sp = dfa_res_val[,1],
                               myPCA$x,
                               group = all_partition_group,
                               dfa_res_val[,-1],
                               uncert = all_partition_uncertainty),species_sub[,c("name_long","code_sp")],by="code_sp")
  kmeans_center <- rep(NA,nfac)
  for(i in 1:nb_group){
    kmeans_center_row <- c()
    for(j in 1:nfac){
      kmeans_center_row <- c(kmeans_center_row,weighted.mean(kmeans_1[kmeans_1$group==i,paste0("X",j)],
                                                             kmeans_1[kmeans_1$group==i,"uncert"]))
    }
    kmeans_center <- rbind(kmeans_center,kmeans_center_row)
  }
  kmeans_center <- kmeans_center[-1,]
  if(is.null(nrow(kmeans_center))){
    kmeans_2 <- data.frame(group=as.factor(1:nb_group),t(kmeans_center),
                           ((kmeans_center - myPCA$center) %*% myPCA$rotation))
    
  }else{
    kmeans_2 <- data.frame(group=as.factor(1:nb_group),kmeans_center,
                           (t(apply(kmeans_center, 1, function(x){x - myPCA$center})) %*% myPCA$rotation))
  }
  
  kmeans_3 <- myPCA$sdev/sum(myPCA$sdev)
  
  kmeans_res <- list(kmeans_1,kmeans_2,kmeans_3)
  
  # PCA centres
  # PC1 and PC2
  pca_centre <- myPCA$rotation[,1:2] %*% matrix(data = c(0,min(myPCA$x[,2]),
                                                         0,max(myPCA$x[,2]),
                                                         min(myPCA$x[,1]),0,
                                                         max(myPCA$x[,1]),0),nrow=2)
  
  pca_centre <- apply(pca_centre,2,function(x){x + myPCA$center})
  
  if(length(kmeans_3)>2){
    # PC1 and PC3
    pca_centre2 <- myPCA$rotation[,c(1,3)] %*% matrix(data = c(0,min(myPCA$x[,3]),
                                                               0,max(myPCA$x[,3]),
                                                               min(myPCA$x[,1]),0,
                                                               max(myPCA$x[,1]),0),nrow=2)
    
    pca_centre2 <- apply(pca_centre2,2,function(x){x + myPCA$center})
    
    # PC2 and PC3
    pca_centre3 <- myPCA$rotation[,2:3] %*% matrix(data = c(0,min(myPCA$x[,3]),
                                                            0,max(myPCA$x[,3]),
                                                            min(myPCA$x[,2]),0,
                                                            max(myPCA$x[,2]),0),nrow=2)
    
    pca_centre3 <- apply(pca_centre3,2,function(x){x + myPCA$center})
  }else{
    pca_centre2 <- pca_centre3 <- NA
  }
  
  pca_centre_list <- list(pca_centre, pca_centre2, pca_centre3)
  
  # Get weigthed centroid of groups
  
  centroids <- as.data.frame(kmeans_2[,c("group",names(kmeans_2)[grepl("PC",names(kmeans_2))])]) 
  
  # Average distance between species and cluster centres
  
  mean_dist_clust <- data.frame(mean_dist=rep(NA,length(unique(kmeans_1$group))))
  if(length(unique(kmeans_1$group))>1){
    for(g in 1:length(unique(kmeans_1$group))){
      kmeans_scale <- as.data.frame(scale(rbind(kmeans_1[kmeans_1$group==g, grepl("X",names(kmeans_1))],kmeans_2[kmeans_2$group==g, grepl("X",names(kmeans_2))])))
      sp_coord <- kmeans_scale[1:(nrow(kmeans_scale)-1),]
      cluster_coord <- kmeans_scale[nrow(kmeans_scale),]
      dist_clust <- c()
      for(i in 1:nrow(sp_coord)){
        mat_dist_clust <- as.matrix(rbind(cluster_coord,sp_coord[i,]))
        dist_clust <- c(dist_clust, dist(mat_dist_clust))
      }
      data_weight_mean <- data.frame(all_dist = dist_clust,
                                     uncert = kmeans_1[kmeans_1$group==g,"uncert"])
      mean_dist_clust[g,1] <- weighted.mean(data_weight_mean$all_dist,data_weight_mean$uncert)
      row.names(mean_dist_clust)[g] <- paste0("cluster_",g)
    }
  }else{
    kmeans_scale <- as.data.frame(scale(rbind(kmeans_1[, grepl("X",names(kmeans_1))],kmeans_center)))
    sp_coord <- kmeans_scale[1:(nrow(kmeans_scale)-1),]
    cluster_coord <- kmeans_scale[nrow(kmeans_scale),]
    dist_clust <- c()
    for(i in 1:nrow(sp_coord)){
      mat_dist_clust <- as.matrix(rbind(cluster_coord,sp_coord[i,]))
      dist_clust <- c(dist_clust, dist(mat_dist_clust))
    }
    mean_dist_clust[1,1] <- mean(dist_clust)
    row.names(mean_dist_clust) <- paste0("cluster_",1)
  }
  
  
  return(list(kmeans_res, # Results of clustering
              centroids, # Position of cluster barycentres
              stability_cluster_final, # Stability of clusters
              mean_dist_clust, # Average distance between species and barycentre
              pca_centre_list # Coordinates to plot trends of PCA axes
  ))
}



## 3) Plot groups and clustering trends

plot_group_boot <- function(nb_group, # Number of clusters
                            centroids, # Position of cluster barycentres
                            kmeans_res, # Results of clustering
                            sdrep, # Optimisation output from DFA
                            nT, # Number of time step
                            min_year, # Oldest year in time-series
                            stability_cluster_final,  # Stability of clusters
                            mean_dist_clust, # Average distance between species and barycentre
                            pca_centre, # Coordinates to plot trends of PCA axes
                            Z_hat, # Factor loadings
                            x_hat, # Latent trends
                            data_ts,
                            data_ts_se
                            ){
  
  # Combine all data for x_pred
  
  data_trend_group <- data.frame(group=rep(paste0("g",1:nb_group),nT),
                                 year=sort(rep(c(min_year:(nT+min_year-1)), nb_group)),
                                 sdrep[grepl("x_pred",row.names(sdrep)) & !grepl("x_pred2",row.names(sdrep)),])
  
  # Set colour code by group
  
  n1 <- length(unique(data_trend_group$group))                                    
  hex_codes1 <- hue_pal()(n1)
  
  # Plot time-series of barycentres
  
  graph <- setNames(lapply(1:nb_group, function(i){
    test <- data_trend_group[data_trend_group$group==paste0("g",i),]
    test$Index_SE <- test$Std..Error
    test$Index <- test$Estimate
    min1 <- min(test$Index-1.96*test$Index_SE) + (max(test$Index+1.96*test$Index_SE)-min(test$Index-1.96*test$Index_SE))/10
    min2 <- min(test$Index-1.96*test$Index_SE)
    ggplot(test, aes(x=year, y=Index)) +
      geom_line(col=hex_codes1[i], size=2) +
      geom_ribbon(aes(ymin=Index-1.96*Index_SE,ymax=Index+1.96*Index_SE),alpha=0.2,fill=hex_codes1[i])+
      xlab(NULL) + 
      ylab(NULL) + 
      annotate("text", x=mean(test$year), y=min1, label= paste0("Stability = ", round(stability_cluster_final[i],3))) +
      annotate("text", x=mean(test$year), y=min2, label= paste0("Mean distance = ", round(mean_dist_clust[i,1],3))) +
      theme_modern() + 
      theme(plot.margin=unit(c(0,0,0,0),"mm"),aspect.ratio = 2/(nb_group+1))
  
    }), levels(as.factor(data_trend_group$group)))
  
  
  # Combine all data for x_pred2
  
  data_trend_group2 <- data.frame(group=rep(c("all",paste0("g",1:nb_group)),nT),
                                 year=sort(rep(c(min_year:(nT+min_year-1)), (nb_group+1))),
                                 sdrep[grepl("x_pred2",row.names(sdrep)),])
  
  geom_mean <- data.frame(Index = apply(data_ts, 2, function(x){ exp(sum(log(x), na.rm=T)/length(x)) }),
                          Index_SE = log(apply(data_ts_se, 2, function(x){ sqrt(sum(exp(2*x))) })),
                          year = as.numeric(names(data_ts)))

  # Plot time-series of barycentres
  
  graph2 <- setNames(lapply(1:(nb_group+1), function(i){
    if(i==1){
      test <- data_trend_group2[data_trend_group2$group=="all",]
    }else{
      if(nb_group==1){
        test <- data_trend_group[data_trend_group$group==paste0("g",(i-1)),]
        min_scale <- min(data_trend_group$Estimate)-2*data_trend_group$Std..Error[which.min(data_trend_group$Estimate)]
        max_scale <- max(data_trend_group$Estimate)+2*data_trend_group$Std..Error[which.max(data_trend_group$Estimate)]
        min1 <- min(data_trend_group$Estimate-1.96*data_trend_group$Std..Error) + (max(data_trend_group$Estimate+1.96*data_trend_group$Std..Error)-min(data_trend_group$Estimate-1.96*data_trend_group$Std..Error))/10
        min2 <- min(data_trend_group$Estimate-1.96*data_trend_group$Std..Error)
      }else{
        test <- data_trend_group2[data_trend_group2$group==paste0("g",(i-1)),]
        min_scale <- min(data_trend_group2$Estimate)-2*data_trend_group2$Std..Error[which.min(data_trend_group2$Estimate)]
        max_scale <- max(data_trend_group2$Estimate)+2*data_trend_group2$Std..Error[which.max(data_trend_group2$Estimate)]
        min1 <- min(data_trend_group2$Estimate-1.96*data_trend_group2$Std..Error) + (max(data_trend_group2$Estimate+1.96*data_trend_group2$Std..Error)-min(data_trend_group2$Estimate-1.96*data_trend_group2$Std..Error))/10
        min2 <- min(data_trend_group2$Estimate-1.96*data_trend_group2$Std..Error)
      }
    }
    test$Index_SE <- test$Std..Error
    test$Index <- test$Estimate
    test$Index_scal <- (test$Estimate-mean(test$Estimate))/sd(test$Estimate)
    test$Index_SE_scal <- test$Std..Error/sd(test$Estimate)
    

    if(i==1){
      
      geom_mean$Index_scal <- scale(geom_mean$Index)
      
      ggplot(test, aes(x=year, y=Index_scal)) +
        geom_line(col="black", size=2) +
        geom_ribbon(aes(ymin=Index_scal-1.96*Index_SE_scal,ymax=Index_scal+1.96*Index_SE_scal),alpha=0.2,fill="black")+
        geom_point(data=geom_mean, col="black", size=2) +
        #geom_pointrange(data=geom_mean, aes(min=Index-1.96*Index_SE,ymax=Index+1.96*Index_SE)) +
        xlab(NULL) + 
        ylab(NULL) + 
        theme_modern() + 
        theme(plot.margin=unit(c(0,0,0,0),"mm"),aspect.ratio = 2/(nb_group+1))
    }else{
      ggplot(test, aes(x=year, y=Index)) +
        geom_line(col=hex_codes1[(i-1)], size=2) +
        geom_ribbon(aes(ymin=Index-1.96*Index_SE,ymax=Index+1.96*Index_SE),alpha=0.2,fill=hex_codes1[(i-1)])+
        xlab(NULL) + 
        ylab(NULL) + 
        ylim(c(min_scale,max_scale)) +
        annotate("text", x=mean(test$year), y=min1, label= paste0("Stability = ", round(stability_cluster_final[(i-1)],3))) +
        annotate("text", x=mean(test$year), y=min2, label= paste0("Mean distance = ", round(mean_dist_clust[(i-1),1],3))) +
        theme_modern() + 
        theme(plot.margin=unit(c(0,0,0,0),"mm"),aspect.ratio = 2/(nb_group+1))
    }
    
    
  }), levels(as.factor(data_trend_group2$group)))
  
  
  # Combine data for PCA centres
  
  mat_tr_rot <- t(solve(varimax(Z_hat)$rotmat) %*% x_hat)
    
  ts_pca <- apply(pca_centre[[1]], 2, function(x){mat_tr_rot %*% matrix(x)})
  min_y_graph3 <- min(apply(ts_pca, 2, min))
  max_y_graph3 <- max(apply(ts_pca, 2, max))
  
  graph3 <- setNames(lapply(1:4, function(i){
    test <- data.frame(Index=ts_pca[,i], year=1:length(ts_pca[,i]))
    ggplot(test, aes(x=year, y=Index)) +
      geom_line(size=1.5, alpha=0.4) + xlab(NULL) + ylab(NULL) +
      theme_modern() + theme_transparent() +
      ylim(c(min_y_graph3,max_y_graph3)) +
      theme(plot.margin=unit(c(0,0,0,0),"mm"),aspect.ratio = 2/3,
            axis.title = element_blank(),
            axis.text = element_blank(),
            axis.text.x = element_blank())
    }), paste0("pca centre ", 1:4))
  
  if(length(kmeans_res[[3]])>2){
    ts_pca2 <- apply(pca_centre[[2]], 2, function(x){mat_tr_rot %*% matrix(x)})
    min_y_graph32 <- min(apply(ts_pca2, 2, min))
    max_y_graph32 <- max(apply(ts_pca2, 2, max))
    
    graph32 <- setNames(lapply(1:4, function(i){
      test <- data.frame(Index=ts_pca2[,i], year=1:length(ts_pca2[,i]))
      ggplot(test, aes(x=year, y=Index)) +
        geom_line(size=1.5, alpha=0.4) + xlab(NULL) + ylab(NULL) +
        theme_modern() + theme_transparent() +
        ylim(c(min_y_graph32,max_y_graph32)) +
        theme(plot.margin=unit(c(0,0,0,0),"mm"),aspect.ratio = 2/3,
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.text.x = element_blank())
    }), paste0("pca centre ", 1:4))
    
    ts_pca3 <- apply(pca_centre[[3]], 2, function(x){mat_tr_rot %*% matrix(x)})
    min_y_graph33 <- min(apply(ts_pca3, 2, min))
    max_y_graph33 <- max(apply(ts_pca3, 2, max))
    
    graph33 <- setNames(lapply(1:4, function(i){
      test <- data.frame(Index=ts_pca3[,i], year=1:length(ts_pca3[,i]))
      ggplot(test, aes(x=year, y=Index)) +
        geom_line(size=1.5, alpha=0.4) + xlab(NULL) + ylab(NULL) +
        theme_modern() + theme_transparent() +
        ylim(c(min_y_graph33,max_y_graph33)) +
        theme(plot.margin=unit(c(0,0,0,0),"mm"),aspect.ratio = 2/3,
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.text.x = element_blank())
    }), paste0("pca centre ", 1:4))
  }
  
  
  # Plot species cluster in the first factorial plane
  
  res_to_plot <- kmeans_res[[1]]
  res_to_plot$group2 <- res_to_plot$group
  res_to_plot$group2 <- as.factor(res_to_plot$group2)
  
  # Add some space between species name and dot
  width_nudge <- (max(res_to_plot$PC1)-min(res_to_plot$PC1))/50
  
  # Species latin name in italic
  res_to_plot$name_long2 <- paste0("italic('",res_to_plot$name_long,"')")
  
  pca_centre_data <- t(matrix(data = c(0,min(kmeans_res[[1]]$PC2)-(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/8,
                                       0,max(kmeans_res[[1]]$PC2)+(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/8,
                                       min(kmeans_res[[1]]$PC1)-(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/8,0,
                                       max(kmeans_res[[1]]$PC1)+(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/8,0),nrow=2))
  
  pca_centre_data2 <- tibble(x=pca_centre_data[,1],
                           y=pca_centre_data[,2],
                           width=0.08,
                           pie = graph3)
  
  if(length(kmeans_res[[3]])>2){
    
    pca_centre_datab <- t(matrix(data = c(0,min(kmeans_res[[1]]$PC3)-(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/8,
                                         0,max(kmeans_res[[1]]$PC3)+(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/8,
                                         min(kmeans_res[[1]]$PC1)-(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/8,0,
                                         max(kmeans_res[[1]]$PC1)+(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/8,0),nrow=2))
    
    pca_centre_data2b <- tibble(x=pca_centre_datab[,1],
                               y=pca_centre_datab[,2],
                               width=0.08,
                               pie = graph32)
    
    
    pca_centre_datac <- t(matrix(data = c(0,min(kmeans_res[[1]]$PC3)-(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/8,
                                          0,max(kmeans_res[[1]]$PC3)+(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/8,
                                          min(kmeans_res[[1]]$PC2)-(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/8,0,
                                          max(kmeans_res[[1]]$PC2)+(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/8,0),nrow=2))
    
    pca_centre_data2c <- tibble(x=pca_centre_datac[,1],
                                y=pca_centre_datac[,2],
                                width=0.08,
                                pie = graph33)
  }
  
  
  final_plot <- ggplot(res_to_plot, aes(PC1,PC2)) +
    geom_point(aes(colour=group2, size=(1-uncert),alpha=uncert)) + 
    geom_text_repel(label=res_to_plot$name_long2, nudge_x = width_nudge, nudge_y = width_nudge, parse = TRUE, max.overlaps = 30) +
    geom_point(data=centroids,aes(x=PC1,y=PC2), shape=15) +
    geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=pca_centre_data2) +
    theme_modern() + xlab(paste0("PC1 (",round(kmeans_res[[3]][1]*100,1)," %)")) +
    ylab(paste0("PC2 (",round(kmeans_res[[3]][2]*100,1)," %)")) +
    xlim(c(min(kmeans_res[[1]]$PC1)-(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/5,
           max(kmeans_res[[1]]$PC1)+(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/5))+
    ylim(c(min(kmeans_res[[1]]$PC2)-(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/5,
           max(kmeans_res[[1]]$PC2)+(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/5))+
    theme(legend.position='none')
  
  if(length(kmeans_res[[3]])>2){
    
    final_plot2 <- ggplot(res_to_plot, aes(PC1,PC3)) +
      geom_point(aes(colour=group2, size=(1-uncert),alpha=uncert)) + 
      geom_text_repel(label=res_to_plot$name_long2, nudge_x = width_nudge, nudge_y = width_nudge, parse = TRUE, max.overlaps = 30) +
      geom_point(data=centroids,aes(x=PC1,y=PC3), shape=15) +
      geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=pca_centre_data2b) +
      theme_modern() + xlab(paste0("PC1 (",round(kmeans_res[[3]][1]*100,1)," %)")) +
      ylab(paste0("PC3 (",round(kmeans_res[[3]][3]*100,1)," %)")) +
      xlim(c(min(kmeans_res[[1]]$PC1)-(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/5,
             max(kmeans_res[[1]]$PC1)+(max(kmeans_res[[1]]$PC1)-min(kmeans_res[[1]]$PC1))/5))+
      ylim(c(min(kmeans_res[[1]]$PC3)-(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/5,
             max(kmeans_res[[1]]$PC3)+(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/5))+
      theme(legend.position='none')
    
    final_plot3 <- ggplot(res_to_plot, aes(PC2,PC3)) +
      geom_point(aes(colour=group2, size=(1-uncert),alpha=uncert)) + 
      geom_text_repel(label=res_to_plot$name_long2, nudge_x = width_nudge, nudge_y = width_nudge, parse = TRUE, max.overlaps = 30) +
      geom_point(data=centroids,aes(x=PC2,y=PC3), shape=15) +
      geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=pca_centre_data2c) +
      theme_modern() + xlab(paste0("PC2 (",round(kmeans_res[[3]][2]*100,1)," %)")) +
      ylab(paste0("PC3 (",round(kmeans_res[[3]][3]*100,1)," %)")) +
      xlim(c(min(kmeans_res[[1]]$PC2)-(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/5,
             max(kmeans_res[[1]]$PC2)+(max(kmeans_res[[1]]$PC2)-min(kmeans_res[[1]]$PC2))/5))+
      ylim(c(min(kmeans_res[[1]]$PC3)-(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/5,
             max(kmeans_res[[1]]$PC3)+(max(kmeans_res[[1]]$PC3)-min(kmeans_res[[1]]$PC3))/5))+
      theme(legend.position='none')
    
  }else{
    final_plot2 <- final_plot3 <- NA
  }
  
  final_plot_list <- list(final_plot, final_plot2, final_plot3)
  
  return(list(final_plot_list, # Plot of species clusters
              graph, # Plot of time-series of cluster barycentres 
              data_trend_group, # Data of time-series of cluster barycentres 
              graph2, # Plot of time-series of cluster barycentres from sdRep
              data_trend_group2 # Data of time-series of cluster barycentres from sdRep
              ))
}

## 4) Core function running the DFA and find the best number of latent trends

core_dfa <- function(data_ts, # Dataset of time series
                     data_ts_se, # Dataset of standard error of time series 
                     nfac, # Number of trends for the DFA
                     AIC=TRUE, # Display AIC
                     silent = TRUE, 
                     control = list()
) {
  
  # Save input data for plot
  
  data_ts_save <- as.data.frame(data_ts)
  data_ts_se_save <- as.data.frame(data_ts_se)
  data_ts_save_long <- cbind(melt(data_ts_save, id.vars=names(data_ts_save)[1]),
                             se=melt(data_ts_se_save, id.vars=names(data_ts_se_save)[1])[,3])
  names(data_ts_save_long)[2:4] <- c("Year","value_orig","se_orig")
  data_ts_save_long$Year <- as.numeric(as.character(data_ts_save_long$Year))
  
  # Remove potential column of ts names
  
  data_ts <- data_ts %>% select_if(Negate(is.character))
  data_ts_se <- data_ts_se %>% select_if(Negate(is.character))
  
  # Overwrite any changes to default control 
  con <- list(nstart = 3, maxit = 10000, reltol = 1e-12, factr = 1e-11, gradtol = 1e-3, nlldeltatol = 1e-4, method = c('NLMINB', 'BFGS'))
  con[names(control)] <- control
  
  # Initialise Z_pred, actual values will be provided by group_from_dfa2
  
  Z_predinit <- matrix(rep(0, 10 * nfac), ncol = nfac)
  W_init <- rbind(1/nrow(data_ts),       # Overall mean, useful to get geometric mean trend across all species
                  rep(0, nrow(data_ts)),rep(0, nrow(data_ts)),
                  rep(0, nrow(data_ts)),rep(0, nrow(data_ts)),
                  rep(0, nrow(data_ts)),rep(0, nrow(data_ts)),
                  rep(0, nrow(data_ts)),rep(0, nrow(data_ts)),
                  rep(0, nrow(data_ts)),rep(0, nrow(data_ts))) # Two mean trends can currently computed, add more rows to W_init to compute more mean trends.
  
  # SE is on log-scale, so no transformation needed
  
  dataTmb <- list(y = log(as.matrix(data_ts)),
                  obs_se = as.matrix(data_ts_se),
                  Z_pred = Z_predinit,
                  W = W_init)
  
  # Prepare parameters for DFA
  
  nfac <- nfac # Number of factors
  ny <- nrow(data_ts) # Number of time series
  nT <- ncol(data_ts) # Number of time step
  
  # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
  constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
  Zmap <- matrix(ncol = nfac, nrow = ny)
  Zmap[constrInd] <- NA
  Zmap[!constrInd] <- 1:sum(!constrInd)
  xmap <- matrix(ncol = nT, nrow = nfac)
  xmap[,1] <- NA
  xmap[(nfac + 1) : (nT * nfac)] <- 1:((nT - 1)* nfac)
  tmbMap <- list(Z = as.factor(Zmap),
                 x = as.factor(xmap))
  
  # Starting values for the optimisation
  
  optList <- vector(con$nstart * length(con$method), mode = 'list')
  names(optList) <- rep(con$method, con$nstart)
  tmbList <- vector(con$nstart * length(con$method), mode = 'list')
  names(tmbList) <- rep(con$method, con$nstart)
  
  for (i in 1:length(optList)) {
    log_re_sp <- runif(ny, -1, 0)
    
    Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
    
    # Set constrained elements to zero
    
    Zinit[constrInd] <- 0
    
    # List of parameters for DFA
    
    tmbPar <-  list(log_re_sp=log_re_sp, Z = Zinit,
                    x=matrix(c(rep(0, nfac), rnorm(nfac * (nT - 1))),
                             ncol = nT, nrow = nfac))
   
    
    # Make DFA
    tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= c("x"), DLL= "dfa_model_se", silent = silent)
    optList[[i]] = switch(names(optList)[i],
                          NLMINB = nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = con$maxit, eval.max  =2*con$maxit, rel.tol =  con$reltol)),
                          BFGS = optim(tmbObj$par, tmbObj$fn, tmbObj$gr, method = 'BFGS', control = list(maxit = con$maxit, reltol = con$reltol)),
                          LBFGS = optim(tmbObj$par, tmbObj$fn, tmbObj$gr, method = 'L-BFGS-B', control = list(maxit = con$maxit, factr = con$factr))
    )
    if (names(optList)[i] == 'NLMINB')
      optList[[i]]$value <- optList[[i]]$objective
      tmbList[[i]] <- tmbObj
  }
  convergence = sapply(optList, FUN = `[[`, 'convergence')
  print(convergence)
  nll = sapply(optList, FUN = `[[`, 'value')
  print(nll)
  maxgrad = sapply(optList, FUN = function(opt) {max(abs(tmbObj$gr(opt$par))) })
  print(maxgrad)
  eligible = abs(nll - min(nll)) < con$nlldeltatol & convergence == 0 & maxgrad < con$gradtol
  if (!any(eligible)) { 
    eligible = abs(nll - min(nll)) < con$nlldeltatol & convergence == 0 # Currently prioritizes optim convergence over gradient check
    if(!any(eligible))
      eligible = abs(nll - min(nll)) < con$nlldeltatol
  } 
  
  ind.best <-  which.min((nll - 1e6 * sign(min(nll)) * !eligible)) # Return the smallest loglikelihood fit that meets other convergence criteria
  tmbOpt <- optList[[ind.best]] 
  tmbObj <- tmbList[[ind.best]] 
  
   sdRep_test_all <- sdreport(tmbObj)
  sdRep_test <- summary(sdRep_test_all)
  
  # Check convergence
  conv <- tmbOpt$convergence
  if(tmbOpt$convergence != 0){warning(paste0("Convergence issue:", tmbOpt$message))}
  if (tmbOpt$convergence == 0 & maxgrad[ind.best] > con$gradtol) 
    warning(paste0('Optimization converged, but maximum gradient = ', maxgrad[ind.best]))
  
  # Compute AIC
  
  if(AIC){
    aic <- AIC.tmb(tmbObj) 
    writeLines(paste('AIC: ', aic))
  } else {aic <- bic <- NA}
  
  return(list(tmbObj, # TMB output
              tmbOpt, # Optimisation from TMB
              data_ts, # Dataset of time series (output)
              data_ts_se, # Dataset of standard error of time series (output)
              data_ts_save, # Dataset of time series (input saved)
              data_ts_save_long, # Dataset of time series (input saved in long format)
              data_ts_se_save, # Dataset of standard error of time series (input saved)
              ny, # Number of time series
              nT, # Number of time step
              aic, # AIC
              conv, # Convergence check
              sdRep_test, # Summary of the TMB optimisation output
              sdRep_test_all # Complete TMB optimisation output
              ))
}





# III) Main function for the DFA-clust analysis

make_dfa <- function(data_ts, # Dataset of time series (species in row, year in column, first column with species name)
                     data_ts_se, # Dataset of log observation error of time series, same dimensions as data_ts
                     nfac = 0, # Number of trends for the DFA, 0 to estimate the best number of trends
                     mintrend = 1, # Minimum number of trends to test
                     maxtrend = 5, # Maximum number of trends to test
                     AIC = TRUE, # Display AIC
                     species_sub,  # Species names
                     nboot = 500, # Number of bootstrap for clustering
                     silent = TRUE, # Silence optimisation
                     control = list(), # Specify changes for DFA control options
                     se_log = TRUE,
                     is_mean_centred = TRUE
                     )
{
  #data_ts=y_farm;data_ts_se=obs_se_farm;nfac=3;mintrend=1;maxtrend=5;AIC=TRUE;species_sub=species_farm;nboot=500;silent = TRUE;control = list();se_log = TRUE
  
  # Save first year for plot
  
  min_year <- as.numeric(colnames(data_ts)[2])
  
  # Log transformed standard errors if they are not (from Taylor expansion)
  # and check missing values in log transformed standard errors
  
  if(anyNA(data_ts_se)){
    if(sum(complete.cases(as.data.frame(t(data_ts_se))))<1){
      stop("Species names must be provided in data_ts_se.")
    }
    if(sum(complete.cases(as.data.frame(t(data_ts_se))))==1){
      warning("Only NAs in data_ts_se, standard errors set to 0.")
      data_ts_se <- as.data.frame(data_ts_se)
      data_ts_se[,-1] <- 0
      data_ts_se <- as.data.table(data_ts_se)
    }
    if(sum(complete.cases(as.data.frame(t(data_ts_se))))>1){
      warning("NAs in data_ts_se, missing values replaced by mean of standard errors.")
      data_ts_se <- as.data.frame(data_ts_se)
      na_coord <- which(is.na(data_ts_se),arr.ind = T)
      data_ts_se[na_coord] <- apply(data_ts_se[na_coord[,1],-1], 1, function(x){return(mean(x, na.rm=T))})
      data_ts_se <- as.data.table(data_ts_se)
    }
  }
  if(se_log == FALSE){
    data_ts_se <- as.data.frame(data_ts_se)
    for(i in 1:nrow(data_ts_se)){
      data_ts_se[i,-1] <- 1/as.numeric(data_ts[i,-1])*as.numeric(data_ts_se[1,-1])
    }
    data_ts_se <- as.data.table(data_ts_se)
  }
  
  # Handle 0 values in time-series (replacing 0 by 1 % of the reference year value)
  
  if(length(which(data_ts==0))>0){
    warning("At least one zero in time-series")
    data_ts <- as.data.frame(data_ts)
    zero_index <- which(data_ts==0, arr.ind = T)
    data_ts[,-1] <- t(apply(data_ts[,-1], 1, function(x){if(length(which(x==0))>0){x[which(x==0)] <- mean(x,na.rm=T)/100}; return(x)}))
    data_ts <- as.data.table(data_ts)
    
    data_ts_se <- as.data.frame(data_ts_se)
    if(anyNA(data_ts_se[zero_index[,1],zero_index[,2]])){
      data_ts_se[zero_index[,1],zero_index[,2]] <- 0
    }
  }
  
  # Mean-centre values if they are not
  
  if(is_mean_centred == FALSE){
    data_ts_prov <- t(apply(data_ts,1,function(x){return(exp(scale(log(as.numeric(x[-1])), scale = F)))}))
    data_ts_prov <- data.table(data_ts[,1],data_ts_prov)
    colnames(data_ts_prov) <- colnames(data_ts)
    
    data_ts_se_prov <- as.matrix(data_ts_se[,-1]) # from Taylor expansion
    for(i in 1:nrow(data_ts_se)){
      data_ts_se_prov[i,] <- 1/as.numeric(data_ts[i,-1])*as.numeric(data_ts_prov[i,-1])*as.numeric(data_ts_se[1,-1])
    }
    data_ts_se_prov <- data.table(data_ts_se[,1],data_ts_se_prov)
    
    data_ts <- data_ts_prov
    data_ts_se <- data_ts_se_prov
    
  }
  
  
  # Run the core_dfa function to find the best number of latent trend if not specified

  if(nfac==0){
    aic_best <- c()
    bic_best <- c()
    for(i in mintrend:maxtrend){
      core_dfa_res <- assign(paste0("core_dfa",i), core_dfa(data_ts=data_ts, data_ts_se=data_ts_se, nfac=i, silent = silent, control = control))
      
      if(core_dfa_res[[11]]==0){
        aic_best <- c(aic_best, core_dfa_res[[10]])
      } else {
        warning(paste0("Convergence issue:", core_dfa_res[[2]]$message))
        aic_best <- c(aic_best, core_dfa_res[[10]])
      }
    }
    if(length(which.min(aic_best))==0){stop("Convergence issues")}
    aic_min <- min(aic_best)
    delta_aic <- c()
    for(aic_ind in 1:length(aic_best)){
      delta_aic[aic_ind] <- aic_best[aic_ind] - aic_min
    }
    nfac <- min(which(delta_aic<2)) # which.min(aic_best)
    
    core_dfa_res <- get(paste0("core_dfa",nfac))
  }else{
    core_dfa_res <- core_dfa(data_ts=data_ts, data_ts_se=data_ts_se, nfac=nfac, silent = silent, control = control)
  }
  tmbObj <- core_dfa_res[[1]]
  tmbOpt <- core_dfa_res[[2]]
  data_ts <- core_dfa_res[[3]]
  data_ts_se <- core_dfa_res[[4]]
  data_ts_save <- core_dfa_res[[5]]
  data_ts_save_long <- core_dfa_res[[6]]
  data_ts_se_save <- core_dfa_res[[7]]
  ny <- core_dfa_res[[8]]
  nT <- core_dfa_res[[9]]
  aic <- core_dfa_res[[10]]
  conv <- core_dfa_res[[11]]
  sdRep_test <- core_dfa_res[[12]]
  sdRep_test_all <- core_dfa_res[[13]]

  if(is.na(aic)){
    stop("Convergence issues")
  }
  
  # Get point estimates

  x_hat <- (tmbObj$env$parList)()$x
  
  Z_hat <- (tmbObj$env$parList)(par=tmbOpt$par)$Z
  
  Z_hat_se <- sdRep_test[rownames(sdRep_test)=="Z",2]
  if(nfac>1){
    for(i in 1:(nfac-1)){
      index_0 <- ny*i
      Z_hat_se <- append(Z_hat_se, rep(0,i), after=index_0)
    }
  }
  Z_hat_se <- matrix(Z_hat_se, ncol=ncol(Z_hat), nrow=nrow(Z_hat))
  
  Z_hat_orig <- sdRep_test[rownames(sdRep_test)=="Z",1]
  if(nfac>1){
    for(i in 1:(nfac-1)){
      index_0 <- ny*i
      Z_hat_orig <- append(Z_hat_orig, rep(0,i), after=index_0)
    }
  }
  Z_hat_orig <- matrix(Z_hat_orig, ncol=ncol(Z_hat), nrow=nrow(Z_hat))
  
  x_hat_se <- matrix(c(rep(0,nfac),sdRep_test[rownames(sdRep_test)=="x",2]), nrow=nfac)
  
  # Covariance matrix of species loadings
  
  cov_mat_Z <- sdRep_test_all$cov.fixed[which(rownames(sdRep_test)=="Z"),which(rownames(sdRep_test)=="Z")]
  
  # Initial data for species loadings
  
  data_loadings <- melt(data.frame(code_sp=data_ts_save[,1],Z_hat_orig),
                        id.vars="code_sp")
  
  # Run group_from_dfa_boot to obtain species clusters
  
  if(nfac>1){
    group_dfa <- group_from_dfa_boot1(data_loadings, cov_mat_Z, species_sub, nboot=nboot, ny, nfac)
    
    if(length(group_dfa[[3]])>1){
      Z_pred_from_kmeans <- as.matrix(group_dfa[[1]][[2]][grepl("X",names(group_dfa[[1]][[2]]))])
      W_from_kmeans <- t(matrix(rep(t(as.matrix(group_dfa[[1]][[1]][grepl("uncert",names(group_dfa[[1]][[1]]))])),(1+length(group_dfa[[3]]))), nrow=nrow(data_ts)))
      for(wg in 1:length(group_dfa[[3]])){
        W_from_kmeans[(wg+1),][group_dfa[[1]][[1]]$group!=wg] <- 0
      }                       
      
    }else{
      
      Z_pred_from_kmeans <- as.matrix(group_dfa[[1]][[2]][grepl("X",names(group_dfa[[1]][[2]]))])
      W_from_kmeans <- rbind(1/nrow(data_ts), rep(0, nrow(data_ts)))
      
    }
    
  }else{
    group_dfa <- 1
    
    Z_pred_from_kmeans <- matrix(rep(0, 10 * nfac), ncol = nfac)
    W_from_kmeans <- rbind(1/nrow(data_ts), rep(0, nrow(data_ts)))
  }
  
  # Update z_pred with actual values after clustering
  
  tmbObj$env$data$Z_pred <- Z_pred_from_kmeans
  tmbObj$env$data$W <- W_from_kmeans
  
  # Recalcul sdreport
  
  sdRep <- summary(sdreport(tmbObj))
  
  # Prepare data to plot
  
  if(!is.character(data_ts_save[,1]) & !is.factor(data_ts_save[,1])){
    data_ts_save <- data.frame(code_sp=as.character(row.names(data_ts_save)),
                               data_ts)
    data_ts_se_save <- data.frame(code_sp=as.character(row.names(data_ts_se_save)),
                                  data_ts_se)
  }else{
    data_ts_save <- data.frame(code_sp=data_ts_save[,1], data_ts)
    data_ts_se_save <- data.frame(code_sp=data_ts_se_save[,1], data_ts_se)
  }
  
  # Back-transform log values of prediction of species ts
  
  sp_ts <- data.frame(code_sp=data_ts_save[,1],
                      matrix(sdRep[rownames(sdRep)=="x_sp",1], nrow=ny))
  
  sp_se_ts <- data.frame(code_sp=data_ts_save[,1],
                         matrix(sdRep[rownames(sdRep)=="x_sp",2], nrow=ny))
  
  # Data for species time-series plot
  
  data_to_plot_sp <- cbind(data_ts_save_long,
                           melt(data_ts_save, id.vars=names(data_ts_save)[1])[,3],
                           se=melt(data_ts_se_save, id.vars=names(data_ts_se_save)[1])[,3],
                           pred=melt(sp_ts, id.vars=names(data_ts_se_save)[1])[,3],
                           pred_se=melt(sp_se_ts, id.vars=names(data_ts_se_save)[1])[,3])
  
  data_to_plot_sp <- merge(data_to_plot_sp, species_sub[,c("name_long","code_sp")],by="code_sp")
  
  
  # Data for DFA trend plot
  
  data_to_plot_tr <- data.frame(t(x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr_se <- data.frame(t(x_hat_se), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  
  
  if(nfac > 1){
    
    # Add rotated trends
    
    data_to_plot_tr_rot <- data.frame(t(solve(varimax(Z_hat)$rotmat) %*% x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
    
    data_to_plot_tr <- cbind(melt(data_to_plot_tr, id.vars = "Year"),
                             se=melt(data_to_plot_tr_se, id.vars = "Year")[,3], # This SE is ok as it comes directly from TMB
                             rot_tr=melt(data_to_plot_tr_rot, id.vars = "Year")[,3])
    
    data_to_plot_tr$variable <- as.character(data_to_plot_tr$variable) %>% gsub(pattern="X", replacement = "Latent trend ") %>% as.factor()
    
    # Data for species loadings
    
    data_loadings <- cbind(melt(data.frame(code_sp=data_ts_save[,1],
                                           Z_hat %*% varimax(Z_hat)$rotmat), id.vars="code_sp"),
                           se.value = NA)
    
    data_loadings <- merge(data_loadings, species_sub[,c("name_long","code_sp")],by="code_sp")
    data_loadings <- merge(data_loadings, group_dfa[[1]][[1]][,c("code_sp","PC1")],by="code_sp")
    data_loadings$name_long <- fct_reorder(data_loadings$name_long,data_loadings$PC1)

    data_loadings$variable <- as.character(data_loadings$variable) %>% gsub(pattern="X", replacement = "Latent trend ") %>% as.factor()
    
    
    # Plots
    
    plot_sp <- ggplot(data_to_plot_sp, aes(x=Year, y=value)) + geom_point() +
      geom_pointrange(aes(ymax = value*exp(1.96 * se.value), ymin=value * exp(-1.96 * se.value))) + 
      geom_line(aes(y=exp(pred.value))) +
      geom_ribbon(aes(y=exp(pred.value), ymax = exp(pred.value+1.96*pred_se.value), ymin=exp(pred.value-1.96*pred_se.value)), alpha=0.5) +
      facet_wrap(name_long ~ ., ncol=round(sqrt(length(unique(data_to_plot_sp$code_sp)))), scales = "free", labeller = label_bquote(col = italic(.(name_long)))) +
      theme_modern() + theme(axis.title.x = element_blank(), axis.title.y = element_blank())
    
    plot_tr <- ggplot(data_to_plot_tr, aes(x=Year, y=rot_tr.value)) + 
      geom_line(aes(colour=variable))+ylab("Rotated value") +
      geom_ribbon(aes(ymax = (rot_tr.value+1.96*se.value), ymin=(rot_tr.value-1.96*se.value), fill=variable), alpha=0.1) +
      facet_wrap(variable ~ ., ncol=min(3,length(unique(data_to_plot_tr$variable)))) +
      theme_modern() + theme(legend.position = "none")
    
    plot_ld <- ggplot(data_loadings) + 
      geom_col(aes(value, name_long, fill=variable)) +
      geom_errorbar(aes(x=value,y=name_long,xmax = value+se.value, xmin=value-se.value), alpha=0.5) +
      facet_wrap(variable ~ ., ncol=length(unique(data_loadings$variable))) +
      theme_modern() + 
      theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(face="italic"))
    
  }else{
    
    data_to_plot_tr <- cbind(melt(data_to_plot_tr, id.vars = "Year"),
                             se=melt(data_to_plot_tr_se, id.vars = "Year")[,3])
    
    # Data for species loadings
    
    data_loadings <- cbind(melt(data.frame(code_sp=data_ts_save[,1],
                                           sdRep[rownames(sdRep)=="Z",1]), id.vars="code_sp"),
                           se.value = NA)
    
    # Plots
    
    plot_sp <- ggplot(data_to_plot_sp, aes(x=Year, y=value)) + geom_point() +
      geom_pointrange(aes(ymax = value*exp(1.96 * se.value), ymin=value * exp(-1.96 * se.value))) + 
      facet_wrap(code_sp ~ ., ncol=4, scales = "free") +
      theme_modern()
    
    plot_tr <- ggplot(data_to_plot_tr, aes(x=Year, y=value)) + 
      geom_line(aes(colour=variable))+
      theme_modern()
    
    plot_ld <- ggplot(data_loadings) + 
      geom_col(aes(value, code_sp, fill=variable)) +
      geom_errorbar(aes(x=value,y=code_sp,xmax = value+se.value, xmin=value-se.value), alpha=0.5) +
      facet_wrap(variable ~ ., ncol=4) +
      theme_modern() + theme(legend.position = "none")
    
  }
  
  if(is.list(group_dfa)){
    if(length(group_dfa[[3]])>1){
      plot_sp_group_all <- plot_group_boot(nb_group = nrow(group_dfa[[1]][[2]]),
                                           centroids = group_dfa[[2]],
                                           kmeans_res = group_dfa[[1]],
                                           sdrep = sdRep, nT = nT,
                                           min_year = min_year,
                                           stability_cluster_final = group_dfa[[3]], 
                                           mean_dist_clust = group_dfa[[4]],
                                           pca_centre = group_dfa[[5]],
                                           Z_hat = Z_hat,
                                           x_hat = x_hat,
                                           data_ts = data_ts,
                                           data_ts_se = data_ts_se)
      plot_sp_group <- plot_sp_group_all[[1]]
      plot_group_ts <- plot_sp_group_all[[2]]
      trend_group <- plot_sp_group_all[[3]]
      plot_group_ts2 <- plot_sp_group_all[[4]]
      trend_group2 <- plot_sp_group_all[[5]]
    }
    if(length(group_dfa[[3]])==1){
      plot_sp_group_all <- plot_group_boot(nb_group = 1,
                                           centroids = group_dfa[[2]],
                                           kmeans_res = group_dfa[[1]],
                                           sdrep = sdRep, nT = nT,
                                           min_year = min_year,
                                           stability_cluster_final = group_dfa[[3]], 
                                           mean_dist_clust = group_dfa[[4]],
                                           pca_centre = group_dfa[[5]],
                                           Z_hat = Z_hat,
                                           x_hat = x_hat,
                                           data_ts = data_ts,
                                           data_ts_se = data_ts_se)
      plot_sp_group <- plot_sp_group_all[[1]]
      plot_group_ts <- plot_sp_group_all[[2]]
      trend_group <- plot_sp_group_all[[3]]
      plot_group_ts2 <- plot_sp_group_all[[4]]
      trend_group2 <- plot_sp_group_all[[5]]
    }
  }
  
  if(!is.list(group_dfa)){
    plot_sp_group_all <- NA
    plot_sp_group <- NA
    plot_group_ts <- NA
    trend_group <- NA
    plot_group_ts2 <- NA
    trend_group2 <- NA
  }
  
  return(list(data_to_plot_sp, # Data on species time-series and fit
              data_to_plot_tr, # Data on latent trends
              data_loadings, # Data on factor loadings
              plot_sp, # Plot of species time-series and fit
              plot_tr, # Plot of latent trends
              plot_ld, # Plot of factor loadings
              plot_sp_group, # Plot clusters in factorial plan
              plot_group_ts, # Plot clustertime-series
              plot_group_ts2, # Plot clustertime-series from sdRep
              aic = aic, # Best AIC
              sdRep = sdRep, # Optimisation output
              group = group_dfa, # Cluster results
              trend_group, # Cluster barycentre times-series
              trend_group2 # Cluster barycentre times-series from sdRep
              ))
}


# Simulation

simul_rand_dfa_intern <- function(cum_perc,
                                  n_sp_init,
                                  nb_group_exp,
                                  thres,
                                  n_y,
                                  n_sp,
                                  sd_rand,
                                  sd_rand2,
                                  sd_ci,
                                  nboot,
                                  seed){
  ## Simulate latent trends
  
  y_init <- data.frame(t(rep(NA,(n_y)))) # latent trends
  test_cor <- 1
  seednum <- 0
  while(abs(test_cor)>0.8){ # check difference between latent trend
    for(i in 1:n_sp_init){
      seednum <- seednum + 1
      y_ts <- c()
      set.seed((seed+seednum))
      y_ts[1] <- rnorm(n = 1, mean = 0, sd = 1)
      for (t in 2:n_y) {
        r.w <- rnorm(n = 1, mean = 0, sd = 1)
        y_ts[t] <- y_ts[t - 1] + r.w
      }
      y_ts <- y_ts + abs(min(y_ts))+1
      y_ts <- exp(scale(log(y_ts)))
      y_init[i,] <- y_ts
    }
    test_cor <- cor.test(as.numeric(y_init[1,]),as.numeric(y_init[2,]),method = "spearman")$estimate
  }
  
  ## From these n_sp_init latent trend, simulate n_sp ts from loading factors
  
  if(nb_group_exp==1){
    id_vec <- c()
    nb_sp_g <- n_sp
    g <- 1
    assign(paste0("nb_sp_g",g),nb_sp_g)
    id_vec <- c(id_vec,rep(g,nb_sp_g))
    for(lt in 1:n_sp_init){
      mean_u_g <- 0 
      lf_u_g <- rnorm(nb_sp_g, mean_u_g, sd_ci)
      assign(paste0("mean_u",lt,"_g",g),mean_u_g) # mean of loading factors in group g for latend trend lt
      assign(paste0("lf_u",lt,"_g",g),lf_u_g) # loading factors for each ts of group g for latend trend lt
    }
  }else{
    id_vec <- c()
    seednum <- 0
    dist_clust <- thres
    coord_clust <- t(matrix(c(dist_clust/sqrt(2),0,0, # distance between summit of a tetrahedron
                              dist_clust/sqrt(2),dist_clust/sqrt(2),dist_clust/sqrt(2),
                              0,dist_clust/sqrt(2),0,
                              0,0,dist_clust/sqrt(2)), ncol=4))
    # t(matrix(c(0,sqrt(3)/(2*sqrt(2))*dist_clust,0, # also distance between summit of a tetrahedron
    # -dist_clust/2,-dist_clust/(2*sqrt(6)),-sqrt(3)/6*dist_clust,
    # dist_clust/2,-dist_clust/(2*sqrt(6)),-sqrt(3)/6*dist_clust,
    # 0,-dist_clust/(2*sqrt(6)),sqrt(3)/3*dist_clust), ncol=4))
    mat_dist <- matrix(NA, ncol=n_sp_init, nrow=nb_group_exp)
    for(g in 1:nb_group_exp){
      nb_sp_g <- round(cum_perc[g]*n_sp/100)
      assign(paste0("nb_sp_g",g),nb_sp_g)
      id_vec <- c(id_vec,rep(g,nb_sp_g))
      for(lt in 1:n_sp_init){
        seednum <- seednum + 1
        set.seed((seed+seednum))
        mean_u_g <- coord_clust[g,lt]
        lf_u_g <- rnorm(nb_sp_g, mean_u_g, sd_ci)
        assign(paste0("mean_u",lt,"_g",g),mean_u_g) # mean of loading factors in group g for latend trend lt
        assign(paste0("lf_u",lt,"_g",g),lf_u_g) # loading factors for each ts of group g for latend trend lt
        mat_dist[g,lt] <- mean_u_g
      }
    }
    if(length(id_vec)>n_sp){ # it may happen because of rounding
      id_vec <- id_vec[1:n_sp]
    }
    if(length(id_vec)<n_sp){ # it may happen because of rounding
      n_sp <- length(id_vec)
    }
  }
  
  y <- data.frame(t(rep(NA,(n_y+2))))
  obs_se <- data.frame(t(rep(NA,(n_y+1))))
  seednum <- 0
  
  for(i in 1:n_sp){ # get simulated ts from loadings
    seednum <- seednum + 1
    set.seed((seed+seednum))
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
    clust_nb <- clust_nb2 <- 1
    clust_stab <- 1
  }else{
    clust_stab <- gsub(", ","-",toString(paste0(round(rand_nfac[[10]][[3]],2))))
    
    clust_nb <- length(unique(obs_group$group))
    clust_nb2 <- data.frame(obs_group %>% group_by(group) %>% summarize(count=n()))
    clust_nb2 <- length(unique(clust_nb2$group[clust_nb2$count>1]))
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
  return(list(res_rand, clust_nb, clust_stab, clust_nb2))
}

simul_rand_dfa_intern2 <- function(cum_perc,n_sp_init,
                                   nb_group_exp,thres,n_y,
                                   n_sp,
                                   sd_rand,
                                   sd_rand2,sd_ci,nboot,seed){
  tryCatch(simul_rand_dfa_intern(cum_perc,n_sp_init,
                                 nb_group_exp,thres,n_y,
                                 n_sp,
                                 sd_rand,
                                 sd_rand2,sd_ci,nboot,seed),
           error=function(e) list(NA,NA,NA,NA))}

simul_rand_dfa <- function(n_y = 20, # number of year
                           n_sp = 15, # number of species ts
                           n_sp_init = 3, # number of latent trends
                           nb_group_exp = 2, # number of expected clusters
                           thres = 1, # distance between barycenters of clusters
                           sd_rand = 0.01, # observation error on data
                           sd_rand2 = 0.5, # random noise on ts
                           sd_ci = 0.1, # standard deviation of the loading factors
                           nboot = 500, # number of bootstrap for clustering
                           equi = TRUE, # equal size of cluster
                           seed = 1 # set seed for simulation data
){
  if(nb_group_exp>1){
    cum_perc <- rep(round(100/nb_group_exp),nb_group_exp)
    skew_val <- 100 - nb_group_exp*10
    cum_perc <- rbind(cum_perc,
                      c(skew_val,rep(round((100-skew_val)/(nb_group_exp-1)),(nb_group_exp-1))))
    if(equi==TRUE){
      cum_perc <- cum_perc[1,]
    }else{
      cum_perc <- cum_perc[2,]
    }
  }else{
    cum_perc <- 100
  }

  
  rand_nfac_list <- simul_rand_dfa_intern2(cum_perc,n_sp_init,nb_group_exp,thres,
                                           n_y,n_sp,sd_rand,sd_rand2,sd_ci,nboot,seed)
  
  res_sim <- data.frame(perc_group=gsub(", ","-",toString(paste0(cum_perc))), value=as.numeric(rand_nfac_list[[1]]),
                        nb_group=rand_nfac_list[[2]], clust_stab=rand_nfac_list[[3]],
                        nb_group2=rand_nfac_list[[4]])
  
  
  return(res_sim)
}

# Multiple plot function

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

