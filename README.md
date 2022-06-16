# Dynamical Factor Analysis and Multispecies clustering

## Reproducible example

### Load packages and functions

```{r}

# Load packages

source("package_used.R")

# Load functions

source('function_dfa_CLEAN.R')

```

### Reproducible data

```{r}
seed_id <- 0
id_vec <- c()
n_sp_init <- 15
nb_group_exp <-2
cum_perc <- c(10,5)
n_y <- 20
n_lt <- 4

y_init <- data.frame(t(rep(NA,n_y))) # latent trends
for(i in 1:n_lt){
    set.seed(i+10)
    y_ts <- c()
    y_ts[1] <- rnorm(n = 1, mean = 0, sd = 1)
    for (t in 2:n_y) {
        r.w <- rnorm(n = 1, mean = 0, sd = 1)
        y_ts[t] <- y_ts[t - 1] + r.w
    }
    y_ts <- y_ts + abs(min(y_ts))+1
    y_ts <- exp(scale(log(y_ts)))
    #max_new <- max(y_ts) - mean(y_ts)/4
    #min_new <- min(y_ts) + mean(y_ts)/4
    #y_ts <- scales::rescale(y_ts, to=c(min_new, max_new))
    y_init[i,] <- y_ts
}

for(g in 1:nb_group_exp){
    nb_sp_g <- cum_perc[g]
    assign(paste0("nb_sp_g",g),nb_sp_g)
    id_vec <- c(id_vec,rep(g,nb_sp_g))
    for(lt in 1:n_sp_init){
        seed_id <- seed_id + 1
        set.seed(seed_id)
        mean_u_g <- runif(1, -1, 1)
        lf_u_g <- rnorm(nb_sp_g, mean_u_g, 0.1)
        assign(paste0("mean_u",lt,"_g",g),mean_u_g) # mean of loading factors in group g for latend trend lt
        assign(paste0("lf_u",lt,"_g",g),lf_u_g) # loading factors for each ts of group g for latend trend lt

    }
}
id_vec <- id_vec[1:n_sp_init]

y <- data.frame(t(rep(NA,(n_y+2))))
obs_se <- data.frame(t(rep(NA,(n_y+1))))

for(i in 1:n_sp_init){ # get simulated ts from loadings
    set.seed(i)
    noise <- rnorm(n_y,0,0.01)
    y[i,1] <- obs_se[i,1] <- sprintf("SP%03d",i)
    y_ts <- rep(0,n_y)
    g <- id_vec[i]
    i_g <- which(which(id_vec==g)==i) # new index for i in group g
    for(lt in 1:n_lt){
        lf_u_g <- get(paste0("lf_u",lt,"_g",g))
        y_ts <- y_ts + as.numeric(y_init[lt,])*lf_u_g[i_g]
    }
    y_ts <- y_ts + noise
    y_ts <- y_ts + abs(min(y_ts)) + 1
    y_ts <- exp(scale(log(y_ts)))
    y[i,2:(n_y+1)] <- y_ts
    y[i,(n_y+2)] <- id_vec[i]
    obs_se[i,2:(n_y+1)] <- abs(rnorm(n_y,0,0.01))
    obs_se[obs_se>1] <- 1
}  


y_ex <- data.table(y[,1:(n_y+1)])
obs_se_ex <- data.table(obs_se)
names(y_ex) <- names(obs_se_ex) <- c("code_sp",1995:(1995+n_y-1))
y_ex$code_sp <- obs_se_ex$code_sp <- sprintf("SP%03d",1:nrow(y_ex))
species_ex <- data.frame(name_long=sprintf("species %03d",1:nrow(y_ex)), code_sp=y_ex$code_sp)

```

### Data specification

Three dataset should be given as input for the main function 'make_dfa. 'data_ts' is the dataset of species time-series and should be provided as a 'data.table' with species time-series in row and the first colum for species names'codes and column names given the year:

```{r}

data_ts <- y_ex

head(data_ts)

```

'data_ts_se' is the dataset of the standard error of species time-series. It should be provided as a 'data.table' with species time-series in row and the first colum for species names'codes and column names given the year:

```{r}

data_ts_se <- obs_se_ex

head(data_ts_se)

```

'species_sub' is the dataset of species names. It should be provided as a 'data.frame' with species in row, the first colum for species latin names and the second column for species names'codes:

```{r}

species_sub <- species_ex

head(species_sub)

```

### Run the analysis

```{r}

ex_dfa_clust <- make_dfa(data_ts = data_ts,
data_ts_se = data_ts_se,
species_sub = species_sub)

```


## Empirical data


