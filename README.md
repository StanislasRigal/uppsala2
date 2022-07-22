# Dynamical Factor Analysis and Multispecies clustering

## Reproducible example

### Load packages and functions

```{r}

# Load packages

source("package_used.R")

# Load functions for time-series

source("function_ts.R")

# Load functions for Dynamic Factor Analysis (DFA) with TMB

source('function_dfa_clean.R')

```

### Reproducible time series

#### Set parameters

```{r}

seed_id <- 0 # starting seed

id_vec <- c() # vactor of species id

n_sp_init <- 15 # number of species time series

nb_group_exp <- 2 # number of expected clusters

cum_perc <- c(10,5) # distribution of species among clusters

n_y <- 20 # number of time steps

n_lt <- 4 # number of latent trends

```

#### Simulate latent trends

```{r}

y_init <- data.frame(t(rep(NA,n_y))) 

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

    y_init[i,] <- y_ts

}

```

#### Simulate cluster barycentres

```{r}

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

```

#### Simulate species time series

```{r}

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
    
    obs_se[i,2:(n_y+1)] <- abs(rnorm(n_y,0,0.1))
    
    obs_se[obs_se>1] <- 1
}  


```

#### Specify data in the riht format

```{r}

y_ex <- data.table(y[,1:(n_y+1)]) # species time series

obs_se_ex <- data.table(obs_se) # observation error on time series

names(y_ex) <- names(obs_se_ex) <- c("code_sp",1995:(1995+n_y-1)) # add years as column name

y_ex$code_sp <- obs_se_ex$code_sp <- sprintf("SP%03d",1:nrow(y_ex)) # add species code

species_ex <- data.frame(name_long=sprintf("species %03d",1:nrow(y_ex)), code_sp=y_ex$code_sp) # species names and code

```

### Data specification

Three datasets should be given as input for the main function `make_dfa`:

- `data_t` contains species time-series and should be provided as a `data.table` with species time-series in row and the first column for species names' codes and years as column names :

```{r}

data_ts <- y_ex

head(data_ts)

```

- `data_ts_se` contains observation errors of species time-series. It should be provided as a `data.table` with species time-series in row and the first column for species names' codes and years as column names:

```{r}

data_ts_se <- obs_se_ex

head(data_ts_se)

```

- `species_sub` contains species names. It should be provided as a `data.frame` with species in row, the first column for complete species names and the second column for species names' codes:

```{r}

species_sub <- species_ex

head(species_sub)

```

### Run the analysis

```{r}

ex_dfa_clust <- make_dfa(data_ts = data_ts, # Dataset of time series

data_ts_se = data_ts_se, # Dataset of observation error of time series 

species_sub = species_sub, # Species names

# optional variables

nfac=0, # Number of trends for the DFA, 0 to estimate the best number of trends

mintrend=1, # Minimum number of trends to test

maxtrend=5, # Maximum number of trends to test

AIC=TRUE, # Display AIC

nboot=500, # Number of bootstrap for clustering

silent = TRUE, # Silence optimisation

control = list()) # Specify changes for DFA control options

```

### Display results

#### Species time-series

The columns are set as follows:

- `code_sp`: code for species names
- `Year`: year
- `value_orig`: input values for species time-series
- `se_orig`: input values for observation error of species time-series
- `value`: back transformed values for species time-series (should be identical to `value_orig`)
- `se.value`:  back transformed values for species time-series (should be identical to `se_orig`)
- `pred.value`: predicted values for species time-series from DFA
- `pred_se.value`: predicted values for standard error of species time-series from DFA
- `name_long`: species names

```{r}

head(ex_dfa_clust[[1]])

```

#### DFA latent trends

The columns are set as follows:

- `Year`: year
- `variable`: latent trend id
- `value`: latent trend values
- `se.value`: standard error of latent trends
- `rot_tr.value`: rotated values for latent trends


```{r}

head(ex_dfa_clust[[2]])

```

#### DFA loading factors

The columns are set as follows:

- `code_sp`: code for species names
- `variable`: latent trend id
- `value`: loading factors
- `se.value`: standard error of loading factors
- `name_long`: species names

```{r}

head(ex_dfa_clust[[3]])

```

#### Plots

```{r}

# Plot species time-series

ex_dfa_clust[[4]]

# Plot latent trends

ex_dfa_clust[[5]]

# Plot loading factors

ex_dfa_clust[[6]]

# Plot species clusters

ex_dfa_clust[[7]][[1]]

# Plot time-series of cluster barycentres

ex_dfa_clust[[7]][[2]]$g1

ex_dfa_clust[[7]][[2]]$g2

```

#### Detailed information on DFA

```{r}

# Akaike information criterion (AIC) of DFA

ex_dfa_clust$aic # or ex_dfa_clust[[8]]

# Optimisation output of DFA

head(ex_dfa_clust$sdRep) # or head(ex_dfa_clust[[9]])

```

#### Main results of species clustering

The columns are set as follows:

- `code_sp`: code for species names
- `PC1`: coordinate on PCA first axis
- `PC2`: coordinate on PCA second axis
- `group`: cluster id
- `Xn`: rotated loading factors for each latent trend *n*
- `uncert`: species stability into its cluster
- `name_long`: species names

```{r}

head(ex_dfa_clust$group[[1]][[1]])

```

#### Detail information on clustering

```{r}

# Cluster barycentres

ex_dfa_clust$group[[1]][[2]]

# Variance captured by PCA first two axes

ex_dfa_clust$group[[1]][[3]]

# Cluster position in the first factorial plan

ex_dfa_clust$group[[2]]

# Cluster stability

ex_dfa_clust$group[[3]]

# Cluster dispersion

ex_dfa_clust$group[[4]]

# Time-series of cluster barycentres

head(ex_dfa_clust[[11]])

```


## Empirical data

### Load and prepare data

#### Bird data

```{r}

# Download and extract data from https://www.gbif.org/occurrence/download?dataset_key=91fa1a0d-a208-40aa-8a6e-f2c0beb9b253 (an free account is necessary) (Type: Darwin Core Archive)

bird_se_raw <- read.csv("raw_data/occurrence.txt", header = T, sep="\t")

# Cleaning data

bird_se_clean <- bird_se_raw[bird_se_raw$class=="Aves",c("class","order","family", "genus","species",
                                                            "specificEpithet","infraspecificEpithet","taxonRank",
                                                            "organismQuantity","decimalLatitude","decimalLongitude",
                                                            "day","month","year","taxonKey","speciesKey","countryCode",
                                                            "level1Gid","level2Gid","iucnRedListCategory")]

# Add a code by species from species name

species_data <- data.frame(name_long = unique(bird_se_clean$species[bird_se_clean$taxonRank != "GENUS"]))

species_data$code_sp <- paste0(toupper(substr(species_data$name_long, 1, 3)),
                               toupper(substr(sub(".* ", "", species_data$name_long), 1, 3)))
                               
species_data$code_sp[species_data$name_long=="Corvus corax"] <- "CORCOX"

species_data$code_sp[species_data$name_long=="Phylloscopus trochilus"] <- "PHYTRU"

species_data$code_sp[species_data$name_long=="Sterna paradisaea"] <- "STEPAD"

species_data$genus <- sub(" .*", "", species_data$name_long)

species_data$species <- sub(".* ", "", species_data$name_long)

bird_se_clean <- merge(bird_se_clean, species_data[,c("code_sp", "name_long")],
                       by.x = c("species"), by.y = c("name_long"), all.x = T)

species_data$class <- bird_se_clean$class[match(species_data$code_sp, bird_se_clean$code_sp)]

species_data$order <- bird_se_clean$order[match(species_data$code_sp, bird_se_clean$code_sp)]

species_data$family <- bird_se_clean$family[match(species_data$code_sp, bird_se_clean$code_sp)]

species_data$iucnRedListCategory <- bird_se_clean$iucnRedListCategory[match(species_data$code_sp, bird_se_clean$code_sp)]

```

#### Geographical coordinates and routes

```{r}

# Create route numbers from coordinates

route_data <- paste0(bird_se_clean$decimalLatitude, sep="_", bird_se_clean$decimalLongitude)

route_data <- data.frame(code_route = paste0("R",str_pad(1:length(unique(route_data)), 3, pad = "0")),
                         coordinate_chr = unique(route_data))
                         
route_data$lat <- as.numeric(sub("_.*", "", route_data$coordinate_chr))

route_data$lon <- as.numeric(sub(".*_", "", route_data$coordinate_chr))

# View routes on the map

worldmap <- ne_countries(scale = 'medium', type = 'countries',returnclass = 'sf')

sweden_map_wgs84 <- worldmap[worldmap$sovereign=="Sweden",]

sweden_map_moll <- sf::st_transform(sweden_map_wgs84, "+proj=moll")

sweden_map_swe <- sf::st_transform(sweden_map_wgs84, "+init=epsg:3006")

route_data_coord <- route_data

coordinates(route_data_coord) <- ~lon+lat

proj4string(route_data_coord) <- CRS("+proj=longlat +datum=WGS84")

route_data_coord <- spTransform(route_data_coord, CRSobj = "+proj=moll")

route_data_moll <- as.data.frame(coordinates(route_data_coord))

route_data_moll$code_route <- route_data_coord$code_route

route_data_coord <- spTransform(route_data_coord, CRSobj = "+init=epsg:3006")

route_data_swe <- as.data.frame(coordinates(route_data_coord))

route_data_swe$code_route <- route_data_coord$code_route

ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data_swe, aes(x=lon, y=lat), width=25000, height=25000, alpha=0.5) +
  theme_void() + coord_sf(datum=NA)

route_data <- merge(route_data, route_data_moll, by="code_route", all=T)

route_data <- merge(route_data, route_data_swe, by="code_route", all=T)

names(route_data)[3:8] <- c("lat_wgs", "lon_wgs", "lon_moll", "lat_moll", "lon_swe", "lat_swe")

route_data$coordinate_chr <- NULL

# Associate routes with the main dataset

bird_se_clean <- merge(bird_se_clean, route_data[,c("code_route", "lat_wgs", "lon_wgs")],
                       by.x = c("decimalLatitude", "decimalLongitude"), by.y = c("lat_wgs", "lon_wgs"), all.x = T)

route_data$level1Gid <- bird_se_clean$level1Gid[match(route_data$code_route, bird_se_clean$code_route)]

route_data$level2Gid <- bird_se_clean$level2Gid[match(route_data$code_route, bird_se_clean$code_route)]

```

#### Ecoregions

````{r}

# Associate routes and regions

municip_data <- municipality # from package swemaps2

municip_data <- sf::st_transform(municip_data, "+init=epsg:3006")

municip_data <- as(municip_data, 'Spatial')

municip_data$area_sqkm <- gArea(municip_data, byid = T)/1000000

route_data_sp <- SpatialPoints(route_data[,c("lon_swe","lat_swe")],
                               proj4string = CRS("+init=epsg:3006"))
route_data <- data.frame(route_data, over(route_data_sp, municip_data))

equival <- data.frame(ln_kod2 = as.character(levels(as.factor(route_data$ln_kod))),
                      level1Gid = as.character(apply(table(route_data$level1Gid,route_data$ln_kod),
                                                     2,function(x){names(which.max(x))})))
                                                     
route_data <- merge(route_data,equival, by="level1Gid", all.x=T)

route_data$ln_kod[is.na(route_data$ln_kod)] <- route_data$ln_kod2[is.na(route_data$ln_kod)]

# Add ecoregions from https://doi.org/10.1016/j.landurbplan.2020.103838

equival2 <- data.frame(ln_kod = as.character(levels(as.factor(route_data$ln_kod))),
                       ecoreg = c(rep("south",12),"north",rep("south",2),rep("north",6)),
                       subecoreg = c(rep("hemiboreal",8),rep("nemoral",3),"hemiboreal",
                                     "s_boreal",rep("hemiboreal",2),rep("s_boreal",4),
                                     rep("n_boreal",2)))

route_data <- merge(route_data,equival2, by="ln_kod", all.x=T)

# Plot map of regions and ecoregions

## Regions

ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data, aes(x=lon_swe, y=lat_swe, fill=ln_kod),
            width=25000, height=25000, alpha=0.5) + scale_fill_viridis(discrete = T) +
  theme_void() + coord_sf(datum=NA)
  
## Ecoregions
  
ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data, aes(x=lon_swe, y=lat_swe, fill=ecoreg),
            width=25000, height=25000, alpha=0.5) + scale_fill_viridis(discrete = T) +
  theme_void() + coord_sf(datum=NA)
  
# Subecoregions

ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data, aes(x=lon_swe, y=lat_swe, fill=subecoreg),
            width=25000, height=25000, alpha=0.5) + scale_fill_viridis(discrete = T) +
  theme_void() + coord_sf(datum=NA)

```

#### Finalise bird dataset

```{r}

# Incorporate 0 in the dataset (species not present while the route was sampled a given year)

bird_se_clean_tot <- dcast(bird_se_clean, countryCode+code_route+year~code_sp,
                           fun.aggregate = sum, value.var="organismQuantity")

bird_se <- melt(bird_se_clean_tot, id.vars = c("countryCode", "code_route", "year"))

names(bird_se)[4:5] <- c("code_sp", "abund")

bird_se$code_sp <- as.character(bird_se$code_sp)

bird_se <- bird_se[bird_se$code_sp!="NA",]

bird_se <- merge(bird_se, route_data, by="code_route", all.x = T)

bird_se$order <- species_data$order[match(bird_se$code_sp, species_data$code_sp)]

bird_se$family <- species_data$family[match(bird_se$code_sp, species_data$code_sp)]

bird_se$genus <- species_data$genus[match(bird_se$code_sp, species_data$code_sp)]

bird_se$species <- species_data$species[match(bird_se$code_sp, species_data$code_sp)]

bird_se$name_long <- species_data$name_long[match(bird_se$code_sp, species_data$code_sp)]

bird_se$iucnRedListCategory <- species_data$iucnRedListCategory[match(bird_se$code_sp, species_data$code_sp)]

```


### Compute species time-series

#### Check species abundance

```{r}

# Look at abundance for each species

# Use data from 1998 because low number of routes monitored in 1996 and 1997: https://doi.org/10.34080/os.v17.22684 )

bird_se_1998 <- droplevels(bird_se[bird_se$year>1997,])

ab_sp <- as.data.frame(bird_se_1998 %>% group_by(code_sp) %>% summarize(ab_tot=sum(abund)))

# Plot by abundance

ggplot(ab_sp, aes(x=reorder(code_sp, -ab_tot, sum), y=ab_tot)) + 
  geom_bar(stat = "identity")+ theme_classic() + scale_y_continuous(trans = "log")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + xlab("Species")+
  ylab("Abundance (log)")

ab_sp2 <- ab_sp[order(ab_sp$ab_tot, decreasing = T),]

ab_sp2$perc_cum <- cumsum(ab_sp2$ab_tot)/sum(ab_sp2$ab_tot)

```

#### Compute species time series

```{r}

# Compute time-series from the SBBS

ts_bird_se_allcountry <- ddply(bird_se_1998, .(code_sp), .fun=get_ts, .progress="text")


# Plot time-series for each species

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
```

### DFA cluster analysis for Swedish birds

#### Complete data

```{r}

# Species names with Swedish and English names

species_data_en_se <- read.csv("output/species_data_en_se.csv", header = T)

# Clean data 

ts_bird_se_allcountry_data <- ts_bird_se_allcountry[which(!is.na(ts_bird_se_allcountry$relative_abundance) & ts_bird_se_allcountry$CI_inf!=0),]

```

#### Farmland birds

```{r}

# Species names and selection

species_sub <- species_farm <-  droplevels(species_data[species_data$code_sp %in% c(
  "FALTIN","VANVAN","ALAARV","HIRRUS","CORFRU",
  "SAXRUB","SYLCOM","ANTPRA","MOTFLA","LANCOL",
  "STUVUL","LINCAN","EMBCIT","EMBHOR","PASMON"),])

# Observation data for farmland birds

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

# Species time series

y_farm <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
           
# Observation error on species time series

obs_se_farm <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
             code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

```


#### Woodland birds

```{r}

# Species names and selection

species_sub <- species_forest <- droplevels(species_data[species_data$code_sp %in% c(
  "ACCNIS","TETBON","TRIOCH","COLOEN","DENMAJ","DRYMAR","PICVIR","JYNTOR",
  "DRYMIN","PICTRI","NUCCAR","GARGLA","PERATE","LOPCRI","POEPAL","POEMON",
  "SITEUR","CERFAM","TURVIS","PHOPHO","PHYCOL","PHYSIB","REGREG","FICHYP",
  "ANTTRI","COCCOC","SPISPI","PYRPYR","EMBRUS"),])

# Observation data for woodland birds

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

# Species time series

y_forest <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
           
# Observation error on species time series

obs_se_forest <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")
                
```

#### All common birds

```{r}

# Species names and selection (107 species in the SBS report)

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
                                                                                  "PODCRI","ARDCIN","ANAPLA","TADTAD", # others
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

# Observation data for common birds

Obs <- ts_bird_se_allcountry_data[ts_bird_se_allcountry_data$code_sp %in% species_sub$code_sp,]

# Reselect species (104) as some are not in the estimated time series

species_all <- droplevels(species_data[species_data$code_sp %in%
                                         levels(as.factor(Obs$code_sp)),])
                                         
# Species time series

y_all <- dcast(Obs[,c("code_sp","relative_abundance_m0","year")],
           code_sp~year, fun.aggregate = sum, value.var = "relative_abundance_m0")
           
# Observation error on species time series
           
obs_se_all <- dcast(Obs[,c("code_sp","Log_SE_m0","year")],
                code_sp~year, fun.aggregate = sum, value.var = "Log_SE_m0")

```

#### Run the DFA-clustering analysis

```{r}

farm_nfac <- make_dfa(data_ts = y_farm, data_ts_se = obs_se_farm,
                       species_sub = species_farm)
                       
forest_nfac <- make_dfa(data_ts = y_forest, data_ts_se = obs_se_forest,
                         species_sub = species_forest)
                         
all_nfac <- make_dfa(data_ts = y_all, data_ts_se = obs_se_all,
                       species_sub = species_all)
                       
```

#### Display results

```{r}

farm_nfac[[7]]

forest_nfac[[7]]

all_nfac[[7]]

```


