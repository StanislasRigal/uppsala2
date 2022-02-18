setwd("/home/rigal/Documents/R/uppsala2")

bird_se_raw2 <- read.csv("raw_data/occurrence.txt", header = T, sep="\t")

# cleaning up data

bird_se_clean <- bird_se_raw2[bird_se_raw2$class=="Aves",c("class","order","family", "genus","species",
                                                            "specificEpithet","infraspecificEpithet","taxonRank",
                                                            "organismQuantity","decimalLatitude","decimalLongitude",
                                                            "day","month","year","taxonKey","speciesKey","countryCode",
                                                            "level1Gid","level2Gid","iucnRedListCategory")]

# add a code by species from species name

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



# create route number from coordinate

route_data <- paste0(bird_se_clean$decimalLatitude, sep="_", bird_se_clean$decimalLongitude)

require(stringr)
route_data <- data.frame(code_route = paste0("R",str_pad(1:length(unique(route_data)), 3, pad = "0")),
                         coordinate_chr = unique(route_data))
route_data$lat <- as.numeric(sub("_.*", "", route_data$coordinate_chr))
route_data$lon <- as.numeric(sub(".*_", "", route_data$coordinate_chr))

# view route on map, need to know the projection

require(sf)
require(rnaturalearth)
require(ggplot2)
require(sp)
require(rgdal)

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

ggplot() + geom_sf(data=sweden_map_moll) +
  geom_tile(data = route_data_moll, aes(x=lon, y=lat), width=25000, height=25000, alpha=0.5) +
  theme_void() + coord_sf(datum=NA)

ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data_swe, aes(x=lon, y=lat), width=25000, height=25000, alpha=0.5) +
  theme_void() + coord_sf(datum=NA)

route_data <- merge(route_data, route_data_moll, by="code_route", all=T)
route_data <- merge(route_data, route_data_swe, by="code_route", all=T)
names(route_data)[3:8] <- c("lat_wgs", "lon_wgs", "lon_moll", "lat_moll", "lon_swe", "lat_swe")
route_data$coordinate_chr <- NULL

# associate routes with the main dataset

bird_se_clean <- merge(bird_se_clean, route_data[,c("code_route", "lat_wgs", "lon_wgs")],
                       by.x = c("decimalLatitude", "decimalLongitude"), by.y = c("lat_wgs", "lon_wgs"), all.x = T)

route_data$level1Gid <- bird_se_clean$level1Gid[match(route_data$code_route, bird_se_clean$code_route)]
route_data$level2Gid <- bird_se_clean$level2Gid[match(route_data$code_route, bird_se_clean$code_route)]


# associate route and region

# require(remotes)
# remotes::install_github("filipwastberg/swemaps2")
require(swemaps2)

municip_data <- municipality
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

# ecoregion from https://doi.org/10.1016/j.landurbplan.2020.103838

equival2 <- data.frame(ln_kod = as.character(levels(as.factor(route_data$ln_kod))),
                       ecoreg = c(rep("south",12),"north",rep("south",2),rep("north",6)),
                       subecoreg = c(rep("hemiboreal",8),rep("nemoral",3),"hemiboreal",
                                     "s_boreal",rep("hemiboreal",2),rep("s_boreal",4),
                                     rep("n_boreal",2)))

route_data <- merge(route_data,equival2, by="ln_kod", all.x=T)

# check
require(viridis)
ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data, aes(x=lon_swe, y=lat_swe, fill=ln_kod),
            width=25000, height=25000, alpha=0.5) + scale_fill_viridis(discrete = T) +
  theme_void() + coord_sf(datum=NA)
ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data, aes(x=lon_swe, y=lat_swe, fill=ecoreg),
            width=25000, height=25000, alpha=0.5) + scale_fill_viridis(discrete = T) +
  theme_void() + coord_sf(datum=NA)
ggplot() + geom_sf(data=sweden_map_swe) +
  geom_tile(data = route_data, aes(x=lon_swe, y=lat_swe, fill=subecoreg),
            width=25000, height=25000, alpha=0.5) + scale_fill_viridis(discrete = T) +
  theme_void() + coord_sf(datum=NA)

# incorporate 0 in the data (species not present while the route is sample this given year)

require(maditr)
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


# save output of this script

saveRDS(bird_se_clean_tot, file = "output/bird_se_clean_tot.rds")
saveRDS(route_data, file = "output/route_data.rds")
saveRDS(species_data, file = "output/species_data.rds")
saveRDS(bird_se, file = "output/bird_se.rds")

