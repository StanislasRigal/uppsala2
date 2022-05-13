library(TMB)

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
  
  // Parameters
  PARAMETER_VECTOR(log_re_sp); // log of sd for random effect by species
  
  // Loadings matrix
  PARAMETER_MATRIX(Z);

  // Latent trends
  PARAMETER_MATRIX(x);
  
  // Cluster center
  matrix<Type> x_pred(Z_pred.rows(), nT);

  // Mean of latent trends
  matrix<Type> x_mean(x.rows(), 1);
  
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
    x_mean(f) = x.row(f).sum() / nT;
    SIMULATE {
      x_mean(f) = x.row(f).sum() / nT;
    }
  }

  // Species trends
  for(int i = 0; i < nSp; ++i) {
    for(int t = 0; t < nT; ++t) {
      x_sp(i, t) = (Z.row(i) * (x.col(t) - x_mean)).sum();
    }
  }  
  
  // Cluster center
  for (int t=0; t < nT; ++t) {
    x_pred.col(t) = Z_pred * (x.col(t) - x_mean);
  } 
  
  
  // Observation model
  for(int i = 0; i < nSp; ++i){
  // Skipping t = 0 when y(i, 0) is fixed at 0. Need to change this if y(i, 0) is not 0.
  // Also had had to change the index of x from t+1 to t, so that x is fixed at zero at time t=0.
    for(int t = 0; t < nT; ++t){ 
       nll -= keep(i) * dnorm(y(i, t), x_sp(i, t), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)), true); // with random effect
      //*----------------------- SECTION I --------------------------*/
        // Simulation block for observation equation
      SIMULATE {
          y(i,t) = rnorm((Z.row(i) * (x.col(t) - x_mean)).sum(), sqrt(obs_se(i, t)*obs_se(i, t)+re_sp(i)*re_sp(i)));
          REPORT(y);
        }
    }  
  }
  
  // State the transformed parameters to report
  // Using ADREPORT will return the point values and the standard errors
  // Note that we only need to specify this for parameters
  // we transformed, see section D above
  // The other parameters, including the random effects (states),
  // will be returned automatically
  ADREPORT(re_sp);
  ADREPORT(x_sp);
  ADREPORT(x_pred);
  ADREPORT(Z_pred);
  
  // Report simulated values
  //SIMULATE{
    //  REPORT(x);
    //  REPORT(y);
    //}
  
  return nll;
  }"

# Compile .cpp file with TMB DFA model
  
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
  

# Function to compute AIC of DFA models (called inside make_dfa)
# Only works if obj has been already optimized
# AIC is computed excluding zero variance

AIC.tmb <- function(obj, tol = 0.01, dontCount) {
  # Simple convergence check
  stopifnot(max(abs(obj$gr(obj$env$last.par.best[obj$env$lfixed()]))) < tol)
  
  # AIC
  
  as.numeric(2 * obj$env$value.best + 2*(sum(obj$env$lfixed()) - dontCount))
}

# Prevent function to stop (to be use in loops)  
AIC.tmb2 <- function(obj, tol = 0.01, dontCount){
  tryCatch(AIC.tmb(obj, tol = 0.01, dontCount),
           error=function(e) NA)}

# Group species

group_from_dfa2 <- function(dfa_res, species_sub, eco_reg=FALSE, weight=FALSE){
  
  # Get loadings from DFA
  dfa_res_val <- dcast(dfa_res, code_sp~variable, value.var = "value")
  dfa_res_se <- dcast(dfa_res, code_sp~variable, value.var = "se.value")
  names(dfa_res_se) <- c("code_sp",paste0("se_",names(dfa_res_val[,-1])))
  
  mat_loading <- as.matrix(dfa_res_val[,-1])
  nb_dim <- ncol(dfa_res_se) -1 
  
  # Calculate weight for each species as the inverse of the volume of the n dimension ellipse (one dimension by DFA trends) defined by SE of each loadings (equivallent to semi axes in the ellipse)
  if(weight==TRUE){
    weight_loading <- apply(dfa_res_se[,-1], 1, 
                            FUN= function(x){
                              vol <- 2/nb_dim * (pi^(nb_dim/2)) / gamma(nb_dim/2) * prod(x)
                              return(vol)
                            })
    weight_loading <- weight_loading/min(weight_loading)
  }else{
    weight_loading <- 1
  }
  
  
  # Calculate gap statistic to find the best number of clusters
  nb <- NbClust(mat_loading, diss=NULL, distance = "euclidean",
                method = "kmeans", min.nc=2, max.nc=10, 
                index = "alllong", alphaBeale = 0.1)
  
  # Plot number of clusters and let the user choose the number of cluster
  print(hist(nb$Best.nc[1,], breaks = max(na.omit(nb$Best.nc[1,]))))
  
  nb_group <- as.numeric(readline(prompt = "Enter number of clusters: "))
  
  # Compute kmeans
  df.kmeans <- cclust(mat_loading, nb_group, weights = 1/weight_loading, 
                      method = "hardcl")

  
  myPCA <- prcomp(mat_loading, scale. = F, center = F)
  
  # Group all info as output
  
  if(eco_reg==FALSE){
    kmeans_res <- list(merge(data.frame(code_sp=dfa_res_val[,1],
                                        myPCA$x[,1:2],
                                        group=as.factor(predict(df.kmeans)),
                                        dfa_res_val[,-1],
                                        dfa_res_se[,-1]),species_sub[,c("name_long","code_sp")],by="code_sp"),
                       data.frame(group=as.factor(1:nb_group),df.kmeans@centers,
                                  df.kmeans@centers %*% myPCA$rotation[,1:2]))
  }else{
    kmeans_res <- list(merge(data.frame(code_sp=dfa_res_val[,1],
                                        myPCA$x[,1:2],
                                        group=as.factor(predict(df.kmeans)),
                                        dfa_res_val[,-1],
                                        dfa_res_se[,-1]),species_sub[,c("name_long_eco","code_sp_eco")],by.x="code_sp", by.y="code_sp_eco"),
                       data.frame(group=as.factor(1:nb_group),df.kmeans@centers,
                                  df.kmeans@centers %*% myPCA$rotation[,1:2]))
  }
  
  
  # Get area (convex hull) for each group for plotting
  
  find_hull <- function(x){x[chull(x$PC2, x$PC1), ]}
  hulls <- ddply(kmeans_res[[1]], "group", find_hull)
  
  # Get centroid of groups
  
  centroids <- as.data.frame(hulls %>% group_by(group) %>% summarize(PC1=mean(PC1), PC2=mean(PC2))) 
  
  return(list(kmeans_res, centroids, hulls))
}


group_from_dfa_boot <- function(data_loadings, cov_mat_Z, species_sub, nboot=100, eco_reg=FALSE, weight=FALSE, ny, nfac){
  
  # Indices of fixed loadings
  
  row_col_0 <- which(data_loadings$value==0)
  
  # Loadings in matrix
  
  dfa_res_val <- dcast(data_loadings, code_sp~variable, value.var = "value")
  mat_loading <- as.matrix(dfa_res_val[,-1])
  
  # Find the best number of clusters in the original data
  
  NbClust2 <- function(data, diss=NULL, distance = "euclidean",
                       method = "kmeans", min.nc=2, max.nc=max(c(2,round(ny/3))), 
                       index = "alllong", alphaBeale = 0.1){
    tryCatch(
      #try to do this
      {
        NbClust(data, diss=NULL, distance = "euclidean",
                method = "kmeans", min.nc=2, max.nc=max(c(2,round(ny/3))), 
                index = "alllong", alphaBeale = 0.1)
      },
      #if an error occurs
      error=function(e) {
        data.frame(Best.partition=rep(NA,nrow(rand_load)))
      }
    )
  }
  
  nb <- NbClust2(mat_loading, diss=NULL, distance = "euclidean",
                 method = "kmeans", min.nc=2, max.nc=max(c(2,round(ny/3))), 
                 index = "alllong", alphaBeale = 0.1)
  
  nb_group_best <- max(nb$Best.partition)
  all_partition <- nb$Best.partition
  
  # Bootstrap to get the optimal number of clusters
  
  for(i in 1:nboot){
    
    # Draw factor loadings using covariance matrix
    set.seed(i)
    rand_load <- rmvnorm(1, mean=data_loadings[!row_col_0,]$value, cov=cov_mat_Z) 
    
    # Complete loading vector with fixed values
    for(j in 1:(nfac-1)){
      index_0 <- ny*j
      rand_load <- append(rand_load, rep(0,j), after=index_0)
    }
    
    rand_load <- matrix(rand_load, ncol=nfac, nrow=ny)

    # Find the best number of clusters in the bootstrap loadings

    nb <- NbClust2(rand_load, diss=NULL, distance = "euclidean",
                  method = "kmeans", min.nc=2, max.nc=max(c(2,round(ny/3))), 
                  index = "alllong", alphaBeale = 0.1)
    
    nb_group_best <- c(nb_group_best, max(nb$Best.partition))
    all_partition <- rbind(all_partition,nb$Best.partition)
    
    # Compute Jaccard similarity to relabel clusters as in the original clustering
    
    jac_sim_res <- matrix(NA, ncol=length(unique(all_partition[1,])),
                          nrow=length(unique(nb$Best.partition)))
    for(k in sort(unique(all_partition[1,]))){
      for(l in sort(unique(nb$Best.partition))){
        jac_sim_mat <- all_partition[c(1,(i+1)),]
        jac_sim_mat[1,][which(jac_sim_mat[1,]!=k)] <- 0
        jac_sim_mat[2,][which(jac_sim_mat[2,]!=l)] <- 0
        jac_sim_mat[jac_sim_mat>0] <- 1
        jac_sim <- c(1 - vegdist(jac_sim_mat, method="jaccard"))
        jac_sim_res[l,k] <- jac_sim
      }
    }
    
    # If same number of clusters
    
    if(length(unique(all_partition[1,]))==length(unique(nb$Best.partition))){
      for(l in sort(unique(nb$Best.partition))){
        all_partition[(i+1),][which(nb$Best.partition==l)] <- which.max(jac_sim_res[l,])
      }
    }
    
    # If more clusters in the bootstrap clustering
    
    if(length(unique(all_partition[1,]))<length(unique(nb$Best.partition))){
      l_data <- c()
      for(k in sort(unique(all_partition[1,]))){
        l_data <- c(l_data,which.max(jac_sim_res[,k]))
      }
      k <- 0
      for(l in l_data){
        k <- k+1
        all_partition[(i+1),][which(nb$Best.partition==l)] <- k
      }
      extra_clus <- sort(unique(nb$Best.partition))[which(!(sort(unique(nb$Best.partition)) %in% l_data))]
      all_partition[(i+1),][which(nb$Best.partition==extra_clus)] <- k +1
    }
    
    # If less clusters in the bootstrap clustering
    
    if(length(unique(all_partition[1,]))>length(unique(nb$Best.partition))){
      k_data <- c()
      for(l in sort(unique(nb$Best.partition))){
        k_data <- c(k_data,which.max(jac_sim_res[l,]))
      }
      l <- 0
      for(k in k_data){
        l <- l+1
        all_partition[(i+1),][which(nb$Best.partition==l)] <- k
      }
    }
  }
  
  all_partition <- na.omit(all_partition)
  nb_group_best <- na.omit(nb_group_best)
  
  # Best number of clusters
  
  nb_group <- as.numeric(names(which.max(table(nb_group_best))))
  
  #all_partition_nb_group <- all_partition[nb_group_best==nb_group,]
  #all_partition_nb_group_df <- as.data.frame(all_partition_nb_group)
  #all_partition_nb_group_df$part <- apply(all_partition_nb_group, 1, function(x){str_c(x, collapse = "-")})
  
  #all_partition2 <- all_partition_nb_group[all_partition_nb_group_df$part==names(which.max(table(all_partition_nb_group_df$part))),][1,]
  #group_num <- sort(unique(all_partition2))
  #for(t in 1:length(group_num)){
  #  all_partition2[which(all_partition2==group_num[t])] <- t
  #}
  
  all_partition2 <- kmeans(mat_loading, nb_group, iter.max = 100)$cluster

  # Bootstrap for cluster stability
  
  for(i in 1:nboot){
    
    # Draw factor loadings using covariance matrix
    set.seed(i)
    rand_load <- rmvnorm(1, mean=data_loadings[!row_col_0,]$value, cov=cov_mat_Z) 
    
    # Complete loading vector with fixed values
    for(j in 1:(nfac-1)){
      index_0 <- ny*j
      rand_load <- append(rand_load, rep(0,j), after=index_0)
    }
    
    rand_load <- matrix(rand_load, ncol=nfac, nrow=ny)
    
    # Find the best number of clusters in the bootstrap loadings
    
    nb <- kmeans(rand_load, nb_group, iter.max = 100)
    
    all_partition2 <- rbind(all_partition2,nb$cluster)
    
    # Compute Jaccard similarity to relabel clusters as in the original clustering
    
    jac_sim_res <- matrix(NA, ncol=length(unique(all_partition2[1,])),
                          nrow=length(unique(nb$cluster)))
    for(k in sort(unique(all_partition2[1,]))){
      for(l in sort(unique(nb$cluster))){
        jac_sim_mat <- all_partition2[c(1,(i+1)),]
        jac_sim_mat[1,][which(jac_sim_mat[1,]!=k)] <- 0
        jac_sim_mat[2,][which(jac_sim_mat[2,]!=l)] <- 0
        jac_sim_mat[jac_sim_mat>0] <- 1
        jac_sim <- c(1 - vegdist(jac_sim_mat, method="jaccard"))
        jac_sim_res[l,k] <- jac_sim
      }
    }
    
    if(i == 1){
      stability_cluster <- apply(jac_sim_res,2,max)
    }else{
      stability_cluster <- rbind(stability_cluster,apply(jac_sim_res,2,max))
    }
    
    if(length(unique(all_partition2[1,]))==length(unique(nb$cluster))){
      for(l in sort(unique(nb$cluster))){
        all_partition2[(i+1),][which(nb$cluster==l)] <- which.max(jac_sim_res[l,])
      }
    }
    
  }
  
  stability_cluster_final <- apply(stability_cluster,2, mean)
  
  all_partition_uncertainty <- apply(all_partition2, 2,
                                     FUN = function(x){
                                       y <- max(table(x))/length(x)
                                       return(y)
                                     })
  all_partition_group <- apply(all_partition2, 2,
                                     FUN = function(x){
                                       xmax <- as.numeric(names(which.max(table(x))))
                                       return(xmax)
                                     })
  if(length(unique(all_partition_group))<nb_group){
    all_partition_group <- all_partition2[1,]
  }
  # Compute PCA to get axes of the graph
  
  myPCA <- prcomp(mat_loading, scale. = F, center = F)
  
  # Group all info as output
  
  if(eco_reg==FALSE){
    kmeans_1 <- merge(data.frame(code_sp = dfa_res_val[,1],
                                 myPCA$x[,1:2],
                                 group = all_partition_group,#as.factor(predict(df.kmeans)),
                                 dfa_res_val[,-1],
                                 #dfa_res_se[,-1],
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
    kmeans_2 <- data.frame(group=as.factor(1:nb_group),kmeans_center,
                           kmeans_center %*% myPCA$rotation[,1:2])
    kmeans_res <- list(kmeans_1,kmeans_2)
  }else{
    kmeans_1 <- merge(data.frame(code_sp = dfa_res_val[,1],
                                 myPCA$x[,1:2],
                                 group = as.factor(predict(df.kmeans)),
                                 dfa_res_val[,-1],
                                 #dfa_res_se[,-1],
                                 uncert = all_partition_uncertainty),species_sub[,c("name_long_eco","code_sp_eco")],by.x="code_sp", by.y="code_sp_eco")
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
    kmeans_2 <- data.frame(group=as.factor(1:nb_group),kmeans_center,
                           kmeans_center %*% myPCA$rotation[,1:2])
    kmeans_res <- list(kmeans_1,kmeans_2)
  }
  
  
  # Get area (convex hull) for each group for plotting
  
  find_hull <- function(x){x[chull(x$PC2, x$PC1), ]}
  hulls <- ddply(kmeans_res[[1]], "group", find_hull)
  
  # Get centroid of groups
  
  centroids <- as.data.frame(hulls %>% group_by(group) %>% summarize(PC1=mean(PC1), PC2=mean(PC2))) 
  
  return(list(kmeans_res, centroids, hulls, stability_cluster_final))
}

# Plot groups and clustering trends

plot_group2 <- function(nb_group, centroids, kmeans_res, hulls, sdrep){
  
  data_trend_group <- data.frame(group=rep(paste0("g",1:nb_group),23),
                                 year=sort(rep(c(1998:2020), nb_group)),
                                 sdrep[grepl("x_pred",row.names(sdrep)),])
  
  
  graph <- setNames(lapply(1:nb_group, function(i){
    test <- data_trend_group[data_trend_group$group==paste0("g",i),]
    test$Index_SE <- test$Std..Error
    test$Index <- test$Estimate
    ggplot(test, aes(x=year, y=Index)) +
      geom_line() +
      geom_ribbon(aes(ymin=Index-Index_SE,ymax=Index+Index_SE),alpha=0.7, col="black",fill="white")+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)
  }), levels(as.factor(data_trend_group$group)))
  
  centroids_data <- tibble(x=centroids$PC1,
                           y=centroids$PC2,
                           width=0.03,
                           pie = graph)
  
  # Get polygons and area
  to_remove <- data.frame(hulls %>% group_by(group) %>% summarize(count=n()))
  
  if(length(levels(droplevels(to_remove$group[to_remove$count>3])))>0){
    hulls <- droplevels(hulls[hulls$group %in% levels(droplevels(to_remove$group[to_remove$count>3])),])
    
    xys <- st_as_sf(hulls, coords=c("PC1","PC2"))
    
    polys <- xys %>% 
      dplyr::group_by(group) %>% 
      dplyr::summarise() %>%
      st_cast("POLYGON") %>%
      st_convex_hull()
    centroids$area <- 0
    centroids$area[centroids$group %in% levels(droplevels(to_remove$group[to_remove$count>3]))] <- st_area(polys)
    
    random_point <- ddply(centroids,.(group),.fun=function(x){
      x <- as.numeric(x)
      y <- data.frame(PC1=rnorm(max(50,round(10000*x[4])),x[2],max(0.01,x[4])),
                      PC2=rnorm(max(50,round(10000*x[4])),x[3],max(0.01,x[4])))
      return(y)
    })
  }else{
    centroids$area <- 0.01
    
    random_point <- ddply(centroids,.(group),.fun=function(x){
      x <- as.numeric(x)
      y <- data.frame(PC1=rnorm(max(50,round(10000*x[4])),x[2],max(0.01,2*x[4])),
                      PC2=rnorm(max(50,round(10000*x[4])),x[3],max(0.01,2*x[4])))
      return(y)
    })
  }
  
  # Plot final output
  final_plot <- ggplot(kmeans_res[[1]], aes(PC1,PC2)) +
    stat_density_2d(data=random_point, aes(PC1,PC2,fill = after_stat(level)), geom = "polygon", alpha=0.2, contour = TRUE,bins = 5) +
    scale_fill_gradient2(low="white", mid="yellow", high = "red", midpoint = nrow(random_point)/10) +
    geom_point() + 
    geom_text(label=kmeans_res[[1]]$name_long, nudge_x = 0.005, nudge_y = 0.005, check_overlap = F) +
    geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroids_data) +
    theme_modern() +
    scale_x_continuous(limits = c(round(min(min(random_point$PC1), min(hulls$PC1)),2)-0.02,round(max(max(random_point$PC1), max(hulls$PC1)),2)+0.02)) +
    scale_y_continuous(limits = c(round(min(min(random_point$PC2), min(hulls$PC2)),2)-0.02,round(max(max(random_point$PC2), max(hulls$PC2)),2)+0.02)) +
    theme(legend.position='none')
  
  return(final_plot)
}


plot_group_boot <- function(nb_group, centroids, kmeans_res, hulls, sdrep, nT){
  
  data_trend_group <- data.frame(group=rep(paste0("g",1:nb_group),nT),
                                 year=sort(rep(c(1998:(nT+1997)), nb_group)),
                                 sdrep[grepl("x_pred",row.names(sdrep)),])
  
  
  graph <- setNames(lapply(1:nb_group, function(i){
    test <- data_trend_group[data_trend_group$group==paste0("g",i),]
    test$Index_SE <- test$Std..Error
    test$Index <- test$Estimate
    ggplot(test, aes(x=year, y=Index)) +
      geom_line() +
      geom_ribbon(aes(ymin=Index-Index_SE,ymax=Index+Index_SE),alpha=0.7, col="black",fill="white")+
      xlab(NULL) + 
      ylab(NULL) + 
      theme_modern() + theme_transparent()+
      theme(plot.margin=unit(c(0,0,0,0),"mm"),axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank(),aspect.ratio = 2/3)
  
    }), levels(as.factor(data_trend_group$group)))
  
  # Plot final output
  res_to_plot <- kmeans_res[[1]]
  res_to_plot$group2 <- res_to_plot$group
  res_to_plot$group2 <- as.factor(res_to_plot$group2)
  
  width_graph <- (max(res_to_plot$PC1)-min(res_to_plot$PC1))/10
  
  centroids_data <- tibble(x=centroids$PC1,
                           y=centroids$PC2,
                           width=width_graph,
                           pie = graph)
  
  width_nudge <- (max(res_to_plot$PC1)-min(res_to_plot$PC1))/50
  
  final_plot <- ggplot(res_to_plot, aes(PC1,PC2)) +
    geom_point(aes(colour=group2, size=(1-uncert),alpha=uncert)) + 
    geom_text(label=res_to_plot$name_long, nudge_x = width_nudge, nudge_y = width_nudge, check_overlap = F) +
    geom_subview(aes(x=x, y=y, subview=pie, width=width, height=width), data=centroids_data) +
    theme_modern() +
    theme(legend.position='none')
  
  return(final_plot)
}

# Generalise DFA

core_dfa <- function(data_ts, # dataset of time series
                     data_ts_se, # dataset of standard error of time series 
                     nfac, # number of trends for the DFA
                     rand_seed=1, # initial values for the sd of the random effect
                     AIC=TRUE # display AIC
                     ){
  
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
  
  # List of data for DFA
  
  # initialise Z_pred, actual values will be provided by group_from_dfa2
  Z_predinit <- matrix(rep(0, 10 * nfac), ncol = nfac)
  
  # SE is on log-scale, so no transformation needed
  dataTmb <- list(y = log(as.matrix(data_ts)),
                  obs_se = as.matrix(data_ts_se),
                  Z_pred = Z_predinit)
  # Prepare parameters for DFA
  
  nfac <- nfac # Number of factors
  ny <- nrow(data_ts) # Number of time series
  nT <- ncol(data_ts) # Number of time step
  
  # Worth trying multiple starting values for the optimisation to check that the right optimum is found.
  set.seed(rand_seed) 
  log_re_sp <- runif(ny, -1, 0)
  
  Zinit <- matrix(rnorm(ny * nfac), ncol = nfac)
  
  # Set constrained elements to zero
  
  constrInd <- rep(1:nfac, each = ny) > rep(1:ny,  nfac)
  Zinit[constrInd] <- 0
  
  # List of parameters for DFA
  
  tmbPar <-  list(log_re_sp=log_re_sp, Z = Zinit,
                  x=matrix(c(rep(0, nfac), rnorm(nfac * (nT - 1))),
                           ncol = nT, nrow = nfac))
  
  # Set up parameter constraints. Elements set to NA will be fixed and not estimated.
  
  Zmap <- matrix(ncol = nfac, nrow = ny)
  Zmap[constrInd] <- NA
  Zmap[!constrInd] <- 1:sum(!constrInd)
  xmap <- matrix(ncol = nT, nrow = nfac)
  xmap[,1] <- NA
  xmap[(nfac + 1) : length(tmbPar$x)] <- 1:(length(tmbPar$x) - nfac)
  tmbMap <- list(Z = as.factor(Zmap),
                 x = as.factor(xmap))
  
  # Make DFA
  
  tmbObj <- MakeADFun(data = dataTmb, parameters = tmbPar, map = tmbMap, random= c("x"), DLL= "dfa_model_se")
  tmbOpt <- nlminb(tmbObj$par, tmbObj$fn, tmbObj$gr, control = list(iter.max = 2000, eval.max  =3000))
  
  # Avoid infinite SE when SD are close or equal to zero
  
  oh <- optimHess(tmbOpt$par, fn=tmbObj$fn, gr=tmbObj$gr)
  singularThreshold = 1e-5
  nSingular = sum(diag(oh) < singularThreshold)
  diag(oh) <- pmax(diag(oh), singularThreshold)
  sdRep_test <- summary(sdreport(tmbObj, hessian.fixed = oh))
  sdRep_test_all <- sdreport(tmbObj, hessian.fixed = oh)
  
  # Check convergence
  
  conv <- grepl("relative convergence",tmbOpt$message)
  if(!conv){warning(paste0("Convergence issue:", tmbOpt$message))}
  
  # Compute AIC
  
  if(AIC){
    aic <- AIC.tmb2(tmbObj, dontCount = 0) 
    aic2 <- AIC.tmb2(tmbObj, dontCount = nSingular) 
    writeLines(paste('AIC: ', aic))
    writeLines(paste('AIC not counting singular random effects: ', aic2))
  } else {aic <- aic2 <- NA}
  
  return(list(tmbObj, tmbOpt, data_ts, data_ts_se, data_ts_save,
              data_ts_save_long, data_ts_se_save, ny, nT,
              aic, aic2, conv, sdRep_test, sdRep_test_all))
}

make_dfa2 <- function(data_ts, # dataset of time series
                     data_ts_se, # dataset of standard error of time series 
                     nfac=0, # number of trends for the DFA, 0 to estimate the best number of trends
                     mintrend=1, # minimum number of trends to test
                     maxtrend=5, # maximum number of trends to test
                     rand_seed=1, # initial values for the sd of the random effect
                     AIC=TRUE, # display AIC
                     species_sub,  # list of species
                     nboot=100,
                     eco_reg=FALSE, weight=FALSE)  # option for group_from_dfa2
{
  #data_ts=y_rand;data_ts_se=obs_se_rand;nfac=0;mintrend=1;maxtrend=5;rand_seed=1;AIC=TRUE;species_sub=species_rand;nboot=100;eco_reg=FALSE;weight=FALSE
  
  if(nfac==0){
    aic2_best <- aic_best <- c()
    for(i in mintrend:maxtrend){
      core_dfa_res <- assign(paste0("core_dfa",i), core_dfa(data_ts=data_ts, data_ts_se=data_ts_se, nfac=i))
      
      if(core_dfa_res[[12]]==T){
        aic_best <- c(aic_best, core_dfa_res[[10]])
        aic2_best <- c(aic2_best, core_dfa_res[[11]])
      } else {
        aic_best <- c(aic_best, NA)
        aic2_best <- c(aic2_best, NA)
      }
    }
    if(length(which.min(aic2_best))==0){stop("Convergence issues")}
    nfac <- which.min(aic2_best)
    core_dfa_res <- get(paste0("core_dfa",nfac))
  }else{
    core_dfa_res <- core_dfa(data_ts=data_ts, data_ts_se=data_ts_se, nfac=nfac)
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
  aic2 <- core_dfa_res[[11]]
  conv <- core_dfa_res[[12]]
  sdRep_test <- core_dfa_res[[13]]
  sdRep_test_all <- core_dfa_res[[14]]
  
  if(is.na(aic2)){
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
  
  # Data for species loadings
  
  cov_mat_Z <- sdRep_test_all$cov.fixed[which(rownames(sdRep_test)=="Z"),which(rownames(sdRep_test)=="Z")]
  
  data_loadings <- melt(data.frame(code_sp=data_ts_save[,1],Z_hat_orig),
                        id.vars="code_sp")
  
  # Get groups from clustering
  if(nfac>1){
    group_dfa <- group_from_dfa_boot(data_loadings, cov_mat_Z, species_sub, nboot=nboot, eco_reg, weight, ny, nfac)
    
    Z_pred_from_kmeans <- as.matrix(group_dfa[[1]][[2]][grepl("X",names(group_dfa[[1]][[2]]))])
    
  }else{
    Z_pred_from_kmeans <- matrix(rep(0, 10 * nfac), ncol = nfac)
  }
  
  # update z_pred with actual values
  
  tmbObj$env$data$Z_pred <- Z_pred_from_kmeans
  
  # Recalcul sdreport
  
  oh <- optimHess(tmbOpt$par, fn=tmbObj$fn, gr=tmbObj$gr)
  singularThreshold = 1e-5
  nSingular = sum(diag(oh) < singularThreshold)
  diag(oh) <- pmax(diag(oh), singularThreshold)
  sdRep <- summary(sdreport(tmbObj, hessian.fixed = oh))
  
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
  
  # Data for DFA trend plot
  
  data_to_plot_tr <- data.frame(t(x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  data_to_plot_tr_se <- data.frame(t(x_hat_se), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
  
  
  if(nfac>1){
    # Add rotated trends
    data_to_plot_tr_rot <- data.frame(t(solve(varimax(Z_hat)$rotmat) %*% x_hat), Year=min(data_to_plot_sp$Year):max(data_to_plot_sp$Year))
    
    data_to_plot_tr <- cbind(melt(data_to_plot_tr, id.vars = "Year"),
                             se=melt(data_to_plot_tr_se, id.vars = "Year")[,3], # This SE is ok as it comes directly from TMB
                             rot_tr=melt(data_to_plot_tr_rot, id.vars = "Year")[,3])
    
    # Data for species loadings
    
    data_loadings <- cbind(melt(data.frame(code_sp=data_ts_save[,1],
                                           Z_hat %*% varimax(Z_hat)$rotmat), id.vars="code_sp"),
                           se.value = NA)
    
    # Plots
    
    plot_sp <- ggplot(data_to_plot_sp, aes(x=Year, y=value)) + geom_point() +
      geom_pointrange(aes(ymax = value*exp(1.96 * se.value), ymin=value * exp(-1.96 * se.value))) + 
      geom_line(aes(y=exp(pred.value))) +
      geom_ribbon(aes(y=exp(pred.value), ymax = exp(pred.value+1.96*pred_se.value), ymin=exp(pred.value-1.96*pred_se.value)), alpha=0.5) +
      facet_wrap(code_sp ~ ., ncol=4, scales = "free") +
      theme_modern()
    
    plot_tr <- ggplot(data_to_plot_tr, aes(x=Year, y=rot_tr.value)) + 
      geom_line(aes(colour=variable))+
      theme_modern()
    
    plot_ld <- ggplot(data_loadings) + 
      geom_col(aes(value, code_sp, fill=variable)) +
      geom_errorbar(aes(x=value,y=code_sp,xmax = value+se.value, xmin=value-se.value), alpha=0.5) +
      facet_wrap(variable ~ ., ncol=4) +
      theme_modern() + theme(legend.position = "none")
    
    plot_sp_group <- plot_group_boot(nb_group = nrow(group_dfa[[1]][[2]]),
                                     centroids = group_dfa[[2]],
                                     kmeans_res = group_dfa[[1]],
                                     hulls = group_dfa[[3]],
                                     sdrep = sdRep, nT = nT)
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
    
    plot_sp_group <- plot_tr
  }
  
  return(list(data_to_plot_sp, data_to_plot_tr, data_loadings,
              plot_sp, plot_tr, plot_ld, plot_sp_group, aic, sdRep))
}

