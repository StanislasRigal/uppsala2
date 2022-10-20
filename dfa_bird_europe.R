# Retrieve data for time series of all European birds by country

# Download Species national indices from Brlik et al. (2021)

df <- fread('https://zenodo.org/record/4590199/files/national_indices2017.csv?download=1')

# Pass from wide to long format

df <- melt(df, id.vars=c("species","euring_code","scheme","type"))
df <- dcast(df, species+euring_code+scheme+variable~type, fun.aggregate = sum)
df <- df[,c("euring_code","species","scheme","variable","index","se")]
names(df) <- c("Code","Species","CountryGroup","Year","Index","Index_SE")
df <- droplevels(na.omit(df))
df$Species <- as.factor(df$Species)
df$CountryGroup <- as.factor(df$CountryGroup)
df$Year <- as.numeric(as.character(df$Year))


# Get species and country names

S <- levels(df$Species)
C <- levels(df$CountryGroup)

# Update species names

df$Species <- as.character(df$Species)
df$Species[df$Species=="Bonasa bonasia"] <- "Tetrastes bonasia"
df$Species[df$Species=="Carduelis cannabina"] <- "Linaria cannabina"
df$Species[df$Species=="Carduelis chloris"] <- "Chloris chloris"
df$Species[df$Species=="Carduelis flammea"] <- "Acanthis flammea"
df$Species[df$Species=="Carduelis spinus"] <- "Spinus spinus"
df_cornix <- df[df$Species=="Corvus corone+cornix",]
df$Species[df$Species=="Corvus corone+cornix"] <- "Corvus corone"
df_cornix$Species[df_cornix$Species=="Corvus corone+cornix"] <- "Corvus cornix"
df <- rbind(df,df_cornix)
df$Species[df$Species=="Corvus monedula"] <- "Coloeus monedula"
df$Species[df$Species=="Delichon urbica"] <- "Delichon urbicum"
df$Species[df$Species=="Dendrocopos medius"] <- "Dendrocoptes medius"
df$Species[df$Species=="Dendrocopos minor"] <- "Dryobates minor"
df$Species[df$Species=="Hippolais pallida"] <- "Iduna pallida"
df$Species[df$Species=="Hirundo daurica"] <- "Cecropis daurica"
df$Species[df$Species=="Hirundo rupestris"] <- "Ptyonoprogne rupestris"
df$Species[df$Species=="Larus ridibundus"] <- "Chroicocephalus ridibundus"
df$Species[df$Species=="Miliaria calandra"] <- "Emberiza calandra"
df$Species[df$Species=="Parus ater"] <- "Periparus ater"
df$Species[df$Species=="Parus caeruleus"] <- "Cyanistes caeruleus"
df$Species[df$Species=="Parus cristatus"] <- "Lophophanes cristatus"
df$Species[df$Species=="Parus montanus"] <- "Poecile montanus"
df$Species[df$Species=="Parus palustris"] <- "Poecile palustris"
df$Species[df$Species=="Saxicola torquata"] <- "Saxicola torquatus"
df$Species[df$Species=="Serinus citrinella"] <- "Carduelis citrinella"
df$Species[df$Species=="Sylvia cantillans"] <- "Curruca cantillans"
df$Species[df$Species=="Sylvia communis"] <- "Curruca communis"
df$Species[df$Species=="Sylvia curruca"] <- "Curruca curruca"
df$Species[df$Species=="Sylvia hortensis"] <- "Curruca hortensis"
df$Species[df$Species=="Sylvia melanocephala"] <- "Curruca melanocephala"
df$Species[df$Species=="Sylvia melanothorax"] <- "Curruca melanothorax"
df$Species[df$Species=="Sylvia nisoria"] <- "Curruca nisoria"
df$Species[df$Species=="Sylvia undata"] <- "Curruca undata"
df$Species[df$Species=="Tetrao tetrix"] <- "Lyrurus tetrix"
df$Species <- as.factor(df$Species)

# Special case of Belgium and Germany

subd <- data.table(reg = c("Belgium-Brussels", "Belgium-Wallonia", "Germany East", "Germany West"), 
                   frac = c(161/30528, 16901/30528, 108333/357022, 248577/357022))

# Belgium

df_bel <- df[df$CountryGroup %in% c("Belgium-Brussels", "Belgium-Wallonia"),]

df_bel <- dcast(df_bel, Species + Year ~ CountryGroup, value.var = c("Index","Index_SE"))

df_bel$Index <- round((df_bel$`Index_Belgium-Brussels`*subd$frac[subd$reg=="Belgium-Brussels"] +
                df_bel$`Index_Belgium-Wallonia`*subd$frac[subd$reg=="Belgium-Wallonia"]) /
                (subd$frac[subd$reg=="Belgium-Brussels"] + subd$frac[subd$reg=="Belgium-Wallonia"]))

df_bel$Index[is.na(df_bel$Index)] <- df_bel$`Index_Belgium-Wallonia`[is.na(df_bel$Index)]

df_bel$Index_SE <- round(sqrt(((df_bel$`Index_SE_Belgium-Brussels`)^2*subd$frac[subd$reg=="Belgium-Brussels"] +
                    (df_bel$`Index_SE_Belgium-Wallonia`)^2*subd$frac[subd$reg=="Belgium-Wallonia"]) /
                        (subd$frac[subd$reg=="Belgium-Brussels"] + subd$frac[subd$reg=="Belgium-Wallonia"])))

df_bel$Index_SE[is.na(df_bel$Index_SE)] <- df_bel$`Index_SE_Belgium-Wallonia`[is.na(df_bel$Index_SE)]

# Germany

df_ger <- df[df$CountryGroup %in% c("Germany East", "Germany West"),]

df_ger <- dcast(df_ger, Species + Year ~ CountryGroup, value.var = c("Index","Index_SE"))

df_ger$Index <- round((df_ger$`Index_Germany East`*subd$frac[subd$reg=="Germany East"] +
                         df_ger$`Index_Germany West` *subd$frac[subd$reg=="Germany West"]) /
                        (subd$frac[subd$reg=="Germany East"] + subd$frac[subd$reg=="Germany West"]))

df_ger$Index[is.na(df_ger$Index)] <- df_ger$`Index_Germany West`[is.na(df_ger$Index)]

df_ger$Index_SE <- round(sqrt(((df_ger$`Index_SE_Germany East`)^2*subd$frac[subd$reg=="Germany East"] +
                                 (df_ger$`Index_SE_Germany West`)^2*subd$frac[subd$reg=="Germany West"]) /
                                (subd$frac[subd$reg=="Germany East"] + subd$frac[subd$reg=="Germany West"])))

df_ger$Index_SE[is.na(df_ger$Index_SE)] <- df_ger$`Index_SE_Germany West`[is.na(df_ger$Index_SE)]


# Merge all data

df$Code <- as.character(df$Code)
df$Code <- str_pad(df$Code, width = 5, pad = "0")
df$Code <- paste0("sp_", df$Code)

species_all <- data.frame(df %>% group_by(Code, Species) %>% summarise(count=n()))
species_all$count <- NULL

df_bel <- merge(df_bel,species_all, by="Species")
df_ger <- merge(df_ger,species_all, by="Species")

df_all_country <- rbind.fill(droplevels(df[!(df$CountryGroup %in% c("Belgium-Brussels", "Belgium-Wallonia",
                                                                    "Germany East", "Germany West")),]),
                             data.frame(df_bel[,c("Code","Species","Year","Index","Index_SE")],CountryGroup="Belgium"),
                             data.frame(df_ger[,c("Code","Species","Year","Index","Index_SE")],CountryGroup="Germany"))

names(df_all_country)[1] <- names(species_all)[1] <- "code_sp"
names(species_all)[2] <- "name_long"

# Select the period of study

df_select_timespan <- as.data.frame(df_all_country %>% group_by(CountryGroup, Year) %>% summarise(count=n()))
df_select_timespan2 <- as.data.frame.matrix(table(df_select_timespan$CountryGroup,df_select_timespan$Year))
df_select_timespan <- dcast(df_select_timespan, CountryGroup ~ Year, value.var = "count", fun.aggregate =  mean)
df_select_timespan <- as.data.frame(df_select_timespan)
df_select_timespan[29,1] <- "Total"
df_select_timespan[29,-1] <- apply(df_select_timespan[,-1], 2, function(x){return(sum(sign(x),na.rm=T))})
plot(unlist(df_select_timespan[29,-1])~as.numeric(names(df_select_timespan[,-1])))

country_to_keep <- df_select_timespan$CountryGroup[which(!is.na(df_select_timespan$`2000`))]
country_to_keep <- country_to_keep[1:(length(country_to_keep)-1)]

df_all_country_2000 <- droplevels(df_all_country[df_all_country$CountryGroup %in% country_to_keep & df_all_country$Year >= 2000,])

# List farmland and woodland species by country

# Austria
# https://info.bml.gv.at/themen/landwirtschaft/eu-agrarpolitik-foerderungen/laendl_entwicklung/le-07-13/evaluierung/le_studien/FBI.html
# https://info.bml.gv.at/themen/landwirtschaft/eu-agrarpolitik-foerderungen/laendl_entwicklung/le-07-13/evaluierung/le_studien/Waldvogelindikator.html

species_aus_farm <- data.frame(name_long=c("Falco tinnunculus","Perdix perdix","Vanellus vanellus","Streptopelia turtur",
                      "Jynx torquilla","Lullula arborea","Alauda arvensis","Anthus trivialis",
                      "Anthus spinoletta","Saxicola rubetra","Saxicola rubicola","Oenanthe oenanthe",
                      "Turdus pilaris","Acrocephalus palustris","Curruca communis","Lanius collurio",
                      "Sturnus vulgaris","Passer montanus","Serinus serinus","Carduelis citrinella",
                      "Carduelis carduelis","Linaria cannabina","Emberiza citrinella","Emberiza calandra"))
#code_sp = paste0(substr(sub("*. ", "", toupper(species_aus_farm)),1,3),substr(sub(".* ", "", toupper(species_aus_farm)),1,3))
species_aus_farm <- merge(species_aus_farm, species_all, by="name_long", all.x=T)
species_sub <- species_aus_farm <- na.omit(species_aus_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Austria",]
species_sub <- species_aus_farm <- species_aus_farm[species_aus_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_aus_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_aus_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_aus_forest <- data.frame(name_long=c("Jynx torquilla","Lullula arborea","Anthus trivialis","Turdus pilaris",
                                       "Serinus serinus"))

species_aus_forest <- merge(species_aus_forest, species_all, by="name_long", all.x=T)
species_sub <- species_aus_forest <- na.omit(species_aus_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Austria",]
species_sub <- species_aus_forest <- species_aus_forest[species_aus_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_forest[y_forest == 0] <- NA

dfa_aus_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                         species_sub = species_aus_forest,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)



# Belgium
# from https://www.aves.be/fileadmin/Aves/Bulletins/Articles/47_1/47_1_1.pdf (Tab1 minus species with only presence absence (x in Frequency))
# and classification from PECBMS (https://pecbms.info/methods/pecbms-methods/3-multispecies-indicators/species-selection-and-classification/  Atlantic)

species_bel_farm <- data.frame(name_long=c("Alauda arvensis","Anthus pratensis","Emberiza citrinella","Hippolais polyglotta",
                                     "Lanius collurio","Motacilla flava","Passer montanus","Perdix perdix",
                                     "Saxicola torquatus","Streptopelia turtur","Curruca communis","Curruca curruca",
                                     "Vanellus vanellus"))

species_bel_farm <- merge(species_bel_farm, species_all, by="name_long", all.x=T)
species_sub <- species_bel_farm <- na.omit(species_bel_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Belgium",]
species_sub <- species_bel_farm <- species_bel_farm[species_bel_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_bel_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_bel_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)


species_bel_forest <- data.frame(name_long=c("Certhia brachydactyla","Certhia familiaris","Coccothraustes coccothraustes","Dendrocopos major",
                                       "Dryobates minor","Dryocopus martius","Erithacus rubecula","Fringilla coelebs",
                                       "Garrulus glandarius","Leiopicus medius","Lophophanes cristatus","Oriolus oriolus",
                                       "Periparus ater","Phoenicurus phoenicurus","Phylloscopus collybita","Phylloscopus sibilatrix",
                                       "Poecile montanus","Poecile palustris","Pyrrhula pyrrhula","Regulus ignicapilla",
                                       "Regulus regulus","Sitta europaea","Sylvia borin","Turdus viscivorus"))

species_bel_forest <- merge(species_bel_forest, species_all, by="name_long", all.x=T)
species_sub <- species_bel_forest <- na.omit(species_bel_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Belgium",]
species_sub <- species_bel_forest <- species_bel_forest[species_bel_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_bel_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                         species_sub = species_bel_forest,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)



# Bulgaria
# https://bspb.org/en/about-birds/common-bird-monitoring/

species_bul_farm  <- data.frame(name_long=c("Coturnix coturnix","Galerida cristata","Alauda arvensis","Motacilla flava",
                                     "Curruca communis","Lanius collurio","Sturnus vulgaris","Emberiza hortulana",
                                     "Emberiza melanocephala","Emberiza calandra"))

species_bul_farm <- merge(species_bul_farm, species_all, by="name_long", all.x=T)
species_sub <- species_bul_farm <- na.omit(species_bul_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Bulgaria",]
species_sub <- species_bul_farm <- species_bul_farm[species_bul_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_bul_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_bul_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)



species_bul_forest  <- data.frame(name_long=c("Columba palumbus","Dendrocopos major","Erithacus rubecula","Fringilla coelebs",
                                       "Garrulus glandarius","Oriolus oriolus","Parus major","Sylvia atricapilla",
                                       "Turdus philomelos"))

species_bul_forest <- merge(species_bul_forest, species_all, by="name_long", all.x=T)
species_sub <- species_bul_forest <- na.omit(species_bul_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Bulgaria",]
species_sub <- species_bul_forest <- species_bul_forest[species_bul_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_bul_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                         species_sub = species_bul_forest,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)


# Cyprus
# https://www.ebcc.info/wp-content/uploads/2020/06/bcn-30-1.pdf#page=33

species_cyp_farm  <- data.frame(name_long=c("Falco tinnunculus","Alectoris chukar","Francolinus francolinus","Coturnix coturnix",
                                     "Columba palumbus","Streptopelia turtur","Clamator glandarius","Athene noctua",
                                     "Coracias garrulus","Galerida cristata","Hirundo rustica","Oenanthe cypriaca",
                                     "Cisticola juncidis","Iduna pallida","Sylvia conspicillata","Curruca melanocephala",
                                     "Parus major","Pica pica","Corvus corone","Passer hispaniolensis",
                                     "Chloris chloris","Carduelis carduelis","Linaria cannabina","Emberiza melanocephala",
                                     "Emberiza calandra"))

species_cyp_farm <- merge(species_cyp_farm, species_all, by="name_long", all.x=T)
species_sub <- species_cyp_farm <- na.omit(species_cyp_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Cyprus",]
species_sub <- species_cyp_farm <- species_cyp_farm[species_cyp_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_cyp_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_cyp_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_cyp_forest  <- data.frame(name_long=c("Columba palumbus","Streptopelia turtur","Troglodytes troglodytes","Oenanthe cypriaca",
                                       "Cettia cetti","Iduna pallida","Curruca melanothorax","Periparus ater",
                                       "Parus major","Certhia brachydactyla","Lanius nubicus","Garrulus glandarius",
                                       "Fringilla coelebs","Serinus serinus","Chloris chloris","Carduelis carduelis",
                                       "Emberiza caesia"))

species_cyp_forest <- merge(species_cyp_forest, species_all, by="name_long", all.x=T)
species_sub <- species_cyp_forest <- na.omit(species_cyp_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Cyprus",]
species_sub <- species_cyp_forest <- species_cyp_forest[species_cyp_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_cyp_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                         species_sub = species_cyp_forest,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

# Czech Republic
# https://www.tandfonline.com/doi/abs/10.1080/00063657.2015.1048423
# https://www.tandfonline.com/doi/abs/10.1080/00063650709461481

species_cze_farm  <- data.frame(name_long=c("Hirundo rustica","Corvus corone","Curruca communis","Falco tinnunculus",
                                     "Linaria cannabina","Sturnus vulgaris","Pica pica","Carduelis carduelis",
                                     "Chloris chloris","Coloeus monedula","Vanellus vanellus","Lanius collurio",
                                     "Emberiza schoeniclus","Alauda arvensis","Passer montanus","Streptopelia turtur",
                                     "Saxicola rubetra","Columba palumbus","Emberiza citrinella"))

species_cze_farm <- merge(species_cze_farm, species_all, by="name_long", all.x=T)
species_sub <- species_cze_farm <- na.omit(species_cze_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Czech Republic",]
species_sub <- species_cze_farm <- species_cze_farm[species_cze_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_cze_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_cze_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_cze_forest  <- data.frame(name_long=c("Buteo buteo","Accipiter nisus","Columba palumbus","Streptopelia turtur",
                                     "Cuculus canorus","Picus viridis","Picus canus","Dryocopus martius",
                                     "Dendrocopos major","Dendrocoptes medius","Dryobates minor","Jynx torquilla",
                                     "Anthus trivialis","Troglodytes troglodytes","Prunella modularis","Erithacus rubecula",
                                     "Luscinia megarhynchos","Phoenicurus phoenicurus","Turdus merula","Turdus philomelos",
                                     "Turdus viscivorus","Hippolais icterina","Sylvia borin","Sylvia atricapilla",
                                     "Phylloscopus sibilatrix","Phylloscopus collybita","Phylloscopus trochilus","Regulus regulus",
                                     "Muscicapa striata","Ficedula hypoleuca","Ficedula albicollis","Aegithalos caudatus",
                                     "Poecile palustris","Poecile montanus","Periparus ater","Cyanistes caeruleus",
                                     "Parus major","Sitta europea","Certhia familiaris","Certhia brachydactyla",
                                     "Fringilla coelebs","Spinus spinus","Acanthis flammea","Pyrrhula pyrrhula",
                                     "Coccothraustes coccothraustes","Oriolus oriolus","Garrulus glandarius"))

species_cze_forest <- merge(species_cze_forest, species_all, by="name_long", all.x=T)
species_sub <- species_cze_forest <- na.omit(species_cze_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Czech Republic",]
species_sub <- species_cze_forest <- species_cze_forest[species_cze_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_cze_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                         species_sub = species_cze_forest,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)


# Denmark
# https://www.dof.dk/fakta-om-fugle/punkttaellingsprogrammet

species_den_farm  <- data.frame(name_long=c("Falco tinnunculus","Perdix perdix","Vanellus vanellus","Gallinago gallinago",
                                     "Alauda arvensis","Hirundo rustica","Anthus pratensis","Motacilla flava",
                                     "Motacilla alba","Saxicola rubetra","Oenanthe oenanthe","Turdus pilaris",
                                     "Curruca curruca","Curruca communis","Lanius collurio","Corvus frugilegus",
                                     "Corvus corone","Passer montanus","Carduelis carduelis","Linaria cannabina",
                                     "Emberiza citrinella","Emberiza calandra"))

species_den_farm <- merge(species_den_farm, species_all, by="name_long", all.x=T)
species_sub <- species_den_farm <- na.omit(species_den_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Denmark",]
species_sub <- species_den_farm <- species_den_farm[species_den_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_den_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_den_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_den_forest  <- data.frame(name_long=c("Accipiter nisus","Columba oenas","Dryocopus martius","Dendrocopos major",
                                       "Erithacus rubecula","Phoenicurus phoenicurus","Turdus viscivorus","Sylvia borin",
                                       "Phylloscopus sibilatrix","Phylloscopus collybita","Regulus regulus","Ficedula hypoleuca",
                                       "Poecile palustris","Lophophanes cristatus","Periparus ater","Sitta europaea",
                                       "Certhia familiaris","Garrulus glandarius","Corvus corax","Fringilla coelebs",
                                       "Spinus spinus","Pyrrhula pyrrhula","Coccothraustes coccothraustes"))

species_den_forest <- merge(species_den_forest, species_all, by="name_long", all.x=T)
species_sub <- species_den_forest <- na.omit(species_den_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Denmark",]
species_sub <- species_den_forest <- species_den_forest[species_den_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[,1:9][y_forest[,1:9] == 0] <- NA
y_forest[y_forest == 0] <- NA
dfa_den_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                         species_sub = species_den_forest,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)


# Estonia
# https://www.sciencedirect.com/science/article/pii/S0167880901002122#TBL2
# plutot https://kirj.ee/public/Ecology/2011/issue_2/ecol-2011-2-88-110.pdf

species_est_farm  <- data.frame(name_long=c("Ciconia ciconia","Crex crex","Tringa totanus","Delichon urbicum",
                                     "Anthus pratensis","Locustella fluviatilis","Curruca nisoria","Muscicapa striata",
                                     "Carpodacus erythrinus","Acrocephalus dumetorum","Corvus corone", "Turdus pilaris",
                                     "Carduelis carduelis","Locustella naevia","Chloris chloris","Falco subbuteo",
                                     "Vanellus vanellus","Linaria cannabina","Pica pica","Circus aeruginosus",
                                     "Acrocephalus palustris","Emberiza hortulana","Perdix perdix","Phasianus colchicus",
                                     "Coturnix coturnix","Lanius collurio","Emberiza schoeniclus","Acrocephalus scirpaceus",
                                     "Acrocephalus schoenobaenus","Alauda arvensis","Sturnus vulgaris","Columba oenas",
                                     "Passer montanus","Oenanthe oenanthe","Saxicola rubetra","Motacilla alba",
                                     "Curruca communis","Columba palumbus","Motacilla flava","Emberiza citrinella"))

species_est_farm <- merge(species_est_farm, species_all, by="name_long", all.x=T)
species_sub <- species_est_farm <- na.omit(species_est_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Estonia",]
species_sub <- species_est_farm <- species_est_farm[species_est_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_est_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_est_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_est_forest  <- data.frame(name_long=c("Luscinia luscinia","Dryocopus martius","Dryobates minor","Phylloscopus sibilatrix",
                                       "Ficedula parva","Poecile montanus","Lophophanes cristatus","Certhia familiaris",
                                       "Oriolus oriolus","Nucifraga caryocatactes","Pyrrhula pyrrhula","Coccothraustes coccothraustes",
                                       "Turdus merula","Cyanistes caeruleus","Fringilla coelebs","Cuculus canorus",
                                       "Prunella modularis","Sylvia borin","Regulus regulus","Dendrocopus major",
                                       "Parus major","Hippolais icterina","Curruca curruca","Poecile palustris",
                                       "Sitta europaea","Ficedula hypoleuca","Corvus corax","Turdus iliacus",
                                       "Erithacus rubecula",	"Spinus spinus",	"Turdus philomelos","Anthus trivialis",
                                       "Phylloscopus trochilus","Lullula arborea"))

species_est_forest <- merge(species_est_forest, species_all, by="name_long", all.x=T)
species_sub <- species_est_forest <- na.omit(species_est_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Estonia",]
species_sub <- species_est_forest <- species_est_forest[species_est_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_forest[y_forest == 0] <- NA

dfa_est_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_est_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Finland
# https://www.biodiversity.fi/ext/en/data-pages/fa8_backgroundinfo.html
# https://www.biodiversity.fi/en/habitats/forests/fo10-forest-birds

species_fin_farm  <- data.frame(name_long=c("Crex crex","Vanellus vanellus","Numenius arquata","Alauda arvensis",
                                     "Hirundo rustica","Delichon urbicum","Anthus pratensis","Saxicola rubertra",
                                     "Turdus pilaris","Curruca communis","Coloeus monedula","Sturnus vulgaris",
                                     "Passer montanus","Emberiza hortulana"))

species_fin_farm <- merge(species_fin_farm, species_all, by="name_long", all.x=T)
species_sub <- species_fin_farm <- na.omit(species_fin_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Finland",]
species_sub <- species_fin_farm <- species_fin_farm[species_fin_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_fin_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_fin_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_fin_forest  <- data.frame(name_long=c("Cuculus canorus","Jynx torquilla","Dendrocopos major","Anthus trivialis",
                                       "Troglodytes troglodytes","Prunella modularis","Luscinia luscinia","Turdus merula",
                                       "Turdus philomelos","Turdus iliacus","Curruca curruca","Sylvia borin",
                                       "Sylvia atricapilla","Phylloscopus trochilus","Muscicapa striata","Ficedula hypoleuca",
                                       "Corvus corax","Fringilla coelebs","Fringilla montifringilla","Acanthis flammea"))

species_fin_forest <- merge(species_fin_forest, species_all, by="name_long", all.x=T)
species_sub <- species_fin_forest <- na.omit(species_fin_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Finland",]
species_sub <- species_fin_forest <- species_fin_forest[species_fin_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_fin_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_fin_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# France
# https://www.vigienature.fr/fr/produire-indicateurs-partir-indices-especes-habitat-2819
# et naturefrance.fr/sites/default/files/2020-08/160513_note_methodologique_indice_stoc.pdf

species_fra_farm  <- data.frame(name_long=c("Vanellus vanellus","Buteo buteo","Falco tinnunculus","Alectoris rufa",
                                     "Perdix perdix","Coturnix coturnix","Upupa epops","Alauda arvensis",
                                     "Lullula arborea","Galerida cristata","Anthus pratensis","Anthus campestris",
                                     "Motacilla flava","Curruca communis","Saxicola torquatus","Saxicola rubetra",
                                     "Oenanthe oenanthe","Lanius collurio","Corvus frugilegus","Linaria cannabina",
                                     "Emberiza citrinella","Emberiza cirlus","Emberiza calandra","Emberiza hortulana"))

species_fra_farm <- merge(species_fra_farm, species_all, by="name_long", all.x=T)
species_sub <- species_fra_farm <- na.omit(species_fra_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "France",]
species_sub <- species_fra_farm <- species_fra_farm[species_fra_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_fra_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_fra_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_fra_forest  <- data.frame(name_long=c("Dendrocopos major","Dendrocoptes medius","Picus canus","Dryocopus martius",
                                       "Curruca melanocephala","Phylloscopus bonelli","Phylloscopus sibilatrix","Phylloscopus collybita",
                                       "Phylloscopus trochilus","Regulus regulus","Regulus ignicapilla","Sitta europaea",
                                       "Certhia brachydactyla","Certhia familiaris","Troglodytes troglodytes","Turdus philomelos",
                                       "Turdus viscivorus","Erithacus rubecula","Lophophanes cristatus","Periparus ater",
                                       "Poecile palustris","Poecile montanus","Coccothraustes coccothraustes","Pyrrhula pyrrhula"))


species_fra_forest <- merge(species_fra_forest, species_all, by="name_long", all.x=T)
species_sub <- species_fra_forest <- na.omit(species_fra_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "France",]
species_sub <- species_fra_forest <- species_fra_forest[species_fra_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_forest[y_forest == 0] <- NA

dfa_fra_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_fra_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Germany
# https://link.springer.com/article/10.1007/s10336-020-01830-4#additional-information

species_ger_farm  <- data.frame(name_long=c("Buteo buteo","Falco tinnunculus","Phasianus colchicus","Coturnix coturnix",
                                     "Sturnus vulgaris","Curruca communis","Emberiza calandra","Alauda arvensis",
                                     "Passer montanus","Streptopelia turtur","Turdus pilaris","Perdix perdix",
                                     "Anthus pratensis","Vanellus vanellus","Milvus milvus","Lanius collurio",
                                     "Saxicola rubetra","Motacilla flava","Emberiza citrinella"))

species_ger_farm <- merge(species_ger_farm, species_all, by="name_long", all.x=T)
species_sub <- species_ger_farm <- na.omit(species_ger_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Germany",]
species_sub <- species_ger_farm <- species_ger_farm[species_ger_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_ger_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_ger_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_ger_forest  <- data.frame(name_long=c("Dryocopus martius","Sylvia atricapilla","Cyanistes caeruleus","Periparus ater",
                                       "Phylloscopus collybita","Corvus corax","Lophophanes cristatus","Pyrrhula pyrrhula",
                                       "Fringilla coelebs","Oriolus oriolus","Garrulus glandarius","Certhia familiaris",
                                       "Ficedula hypoleuca","Erithacus rubecula","Regulus ignicapilla","Sylvia borin",
                                       "Regulus regulus","Dendrocopos major","Picus canus","Coccothraustes coccothraustes",
                                       "Dryobates minor","Aegithalos caudatus","Poecile palustris","Dendrocoptes medius",
                                       "Turdus viscivorus","Certhia brachydactyla","Turdus philomelos","Columba oenas",
                                       "Anthus trivialis","Poecile montanus","Phylloscopus trochilus","Troglodytes troglodytes",
                                       "Sitta europaea","Phylloscopus sibilatrix"))


species_ger_forest <- merge(species_ger_forest, species_all, by="name_long", all.x=T)
species_sub <- species_ger_forest <- na.omit(species_ger_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Germany",]
species_sub <- species_ger_forest <- species_ger_forest[species_ger_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_ger_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_ger_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Greece
# same as EU, https://www.ornithologiki.gr/en/our-work/conservation-scientific-research/bird-monitoring/1298-programma-parakoloythisis-ton-koinon-eidon-poulion-tis-elladas
# asked

# Hungary
# https://www.mme.hu/tovabbra-csokken-hazai-mezogazdasagi-teruletek-madarvilaga

species_hun_farm  <- data.frame(name_long=c("Falco tinnunculus","Perdix perdix","Coturnix coturnix","Vanellus vanellus",
                                     "Merops apiaster","Galerida cristata","Alauda arvensis","Anthus campestris",
                                     "Motacilla flava","Locustella naevia","Curruca nisoria","Curruca communis",
                                     "Lanius collurio","Lanius minor","Sturnus vulgaris","Emberiza calandra"))

species_hun_farm <- merge(species_hun_farm, species_all, by="name_long", all.x=T)
species_sub <- species_hun_farm <- na.omit(species_hun_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Hungary",]
species_sub <- species_hun_farm <- species_hun_farm[species_hun_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_hun_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_hun_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_hun_forest  <- data.frame(name_long=c("Columba oenas","Dryocopus martius","Dendrocopos major","Dendrocoptes medius",
                                       "Dryobates minor","Lullula arborea","Troglodytes troglodytes","Prunella modularis",
                                       "Erithacus rubecula","Turdus philomelos","Turdus viscivorus","Phylloscopus sibilatrix",
                                       "Phylloscopus collybita","Ficedula albicollis","Poecile palustris","Periparus ater",
                                       "Cyanistes caeruleus","Sitta europaea","Certhia brachydactyla","Garrulus glandarius",
                                       "Fringilla coelebs","Coccothraustes coccothraustes"))


species_hun_forest <- merge(species_hun_forest, species_all, by="name_long", all.x=T)
species_sub <- species_hun_forest <- na.omit(species_hun_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Hungary",]
species_sub <- species_hun_forest <- species_hun_forest[species_hun_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_forest[y_forest == 0] <- NA

dfa_hun_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_hun_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Italy
# https://www.reterurale.it/flex/cm/pages/ServeBLOB.php/L/IT/IDPagina/15032

species_ita_farm  <- data.frame(name_long=c("Alauda arvensis","Lanius collurio","Motacilla alba","Melanocorypha calandra",
                                     "Calandrella brachydactyla","Anthus campestris","Galerida cristata","Carduelis carduelis",
                                     "Corvus cornix","Motacilla flava","Pica pica","Falco tinnunculus",
                                     "Emberiza hortulana","Passer italiae","Passer montanus","Passer hispaniolensis",
                                     "Oriolus oriolus","Hirundo rustica","Saxicola torquatus","Sturnus vulgaris",
                                     "Sturnus unicolor","Emberiza calandra","Jynx torquilla","Streptopelia turtur",
                                     "Upupa epops","Luscinia megarhynchos","Chloris chloris","Serinus serinus"))

species_ita_farm <- merge(species_ita_farm, species_all, by="name_long", all.x=T)
species_sub <- species_ita_farm <- na.omit(species_ita_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Italy",]
species_sub <- species_ita_farm <- species_ita_farm[species_ita_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_ita_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_ita_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_ita_forest  <- data.frame(name_long=c("Poecile montanus","Poecile palustris","Lophophanes cristatus","Periparus ater",
                                       "Cyanistes caeruleus","Pyrrhula pyrrhula","Aegithalos caudatus","Regulus ignicapilla",
                                       "Fringilla coelebs","Garrulus glandarius","Phylloscopus bonelli","Phylloscopus collybita",
                                       "Nucifraga caryocatactes","Erithacus rubecula","Sitta europaea","Dryocopus martius",
                                       "Dendrocopos major","Certhia familiaris","Certhia brachydactyla","Regulus regulus",
                                       "Troglodytes troglodytes","Turdus viscivorus","Turdus philomelos"))


species_ita_forest <- merge(species_ita_forest, species_all, by="name_long", all.x=T)
species_sub <- species_ita_forest <- na.omit(species_ita_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Italy",]
species_sub <- species_ita_forest <- species_ita_forest[species_ita_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_ita_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_ita_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Latvia
# https://www.researchgate.net/profile/Oskars-Keiss/publication/268363185_Experiences_with_a_Baseline_Indicator_Farmland_Bird_Index_in_Latvia/links/55097a2f0cf26ff55f859259/Experiences-with-a-Baseline-Indicator-Farmland-Bird-Index-in-Latvia.pdf
# https://www.lu.lv/fileadmin/user_upload/LU.LV/Apaksvietnes/Konferences/Comm_ForestBirds_Trends_LV_AuninsA_20181205.pdf

species_lat_farm  <- data.frame(name_long=c("Ciconia ciconia","Crex crex","Vanellus vanellus","Alauda arvensis",
                                     "Anthus pratensis","Locustella naevia","Acrocephalus palustris","Saxicola rubetra",
                                     "Carduelis carduelis","Linaria cannabina","Carpodacus erythrinus","Emberiza citrinella"))

species_lat_farm <- merge(species_lat_farm, species_all, by="name_long", all.x=T)
species_sub <- species_lat_farm <- na.omit(species_lat_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Latvia",]
species_sub <- species_lat_farm <- species_lat_farm[species_lat_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_lat_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_lat_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_lat_forest  <- data.frame(name_long=c("Accipiter gentilis","Accipiter nisus","Tetrastes bonasia","Picus canus",
                                       "Dryocopus martius","Dryobates minor","Dendrocopos leucotos","Picoides tridactylus",
                                       "Turdus viscivorus","Phylloscopus sibilatrix","Regulus regulus","Ficedula parva",
                                       "Ficedula hypoleucos","Aegithalos caudatus","Poecile palustris","Poecile montanus",
                                       "Lophophanes cristatus","Periparus ater","Certhia familiaris","Nucifraga caryocatactes",
                                       "Loxia curvirostra","Pyrrhula pyrrhula","Coccothraustes coccothraustes"))


species_lat_forest <- merge(species_lat_forest, species_all, by="name_long", all.x=T)
species_sub <- species_lat_forest <- na.omit(species_lat_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Latvia",]
species_sub <- species_lat_forest <- species_lat_forest[species_lat_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA
# no because from 2005 for forest
dfa_lat_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_lat_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Lithuania
# no official forest bird indicator

species_lit_farm  <- data.frame(name_long=c("Ciconia ciconia","Crex crex","Vanellus vanellus","Alauda arvensis",
                                            "Hirundo rustica","Anthus pratensis","Motacilla flava","Saxicola rubetra",
                                            "Curruca communis","Lanius collurio","Sturnus vulgaris","Passer montanus",
                                            "Carduelis carduelis","Emberiza citrinella"))

species_lit_farm <- merge(species_lit_farm, species_all, by="name_long", all.x=T)
species_sub <- species_lit_farm <- na.omit(species_lit_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                             df_all_country_2000$CountryGroup == "Lithuania",]
species_sub <- species_lit_farm <- species_lit_farm[species_lit_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_lit_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_lit_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)



# Luxembourg
# asked

species_lux_farm  <- data.frame(name_long=c("Alauda arvensis","Saxicola torquatus","Curruca communis","Lanius collurio",
                                     "Passer montanus","Linaria cannabina","Emberiza citrinella"))

species_lux_farm <- merge(species_lux_farm, species_all, by="name_long", all.x=T)
species_sub <- species_lux_farm <- na.omit(species_lux_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Luxembourg",]
species_sub <- species_lux_farm <- species_lux_farm[species_lux_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_lux_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_lux_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_lux_forest  <- data.frame(name_long=c("Columba oenas","Dryocopus martius","Leiopicus medius","Anthus trivialis",
                                       "Turdus viscivorus","Phylloscopus sibilatrix","Regulus regulus","Regulus ignicapillus",
                                       "Sitta europaea","Certhia familiaris","Coccothraustes coccothraustes"))


species_lux_forest <- merge(species_lux_forest, species_all, by="name_long", all.x=T)
species_sub <- species_lux_forest <- na.omit(species_lux_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Luxembourg",]
species_sub <- species_lux_forest <- species_lux_forest[species_lux_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_forest[y_forest == 0] <- NA

dfa_lux_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_lux_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)


# Netherlands
# https://stats.sovon.nl/pub/publicatie/18092
# https://www.clo.nl/indicatoren/nl1618-broedvogels-van-het-bos

species_net_farm  <- data.frame(name_long=c("Perdix perdix","Streptopelia turtur","Vanellus vanellus","Limosa limosa",
                                     "Sturnus vulgaris","Haematopus ostralegus","Numenius arquata ","Tringa totanus",
                                     "Aythya fuligula","Spatula querquedula","Falco tinnunculus ","Coturnix coturnix",
                                     "Passer montanus","Spatula clypeata","Emberiza citrinella","Anas crecca",
                                     "Gallinago gallinago","Turdus viscivorus ","Emberiza calandra","Alauda arvensis ",
                                     "Motacilla flava ","Hippolais icterina ","Hirundo rustica","Saxicola rubicola",
                                     "Mareca strepera ","Curruca communis","Carduelis carduelis"))

species_net_farm <- merge(species_net_farm, species_all, by="name_long", all.x=T)
species_sub <- species_net_farm <- na.omit(species_net_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Netherlands",]
species_sub <- species_net_farm <- species_net_farm[species_net_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_net_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_net_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_net_forest  <- data.frame(name_long=c("Coccothraustes coccothraustes","Ficedula hypoleuca","Sitta europaea","Certhia brachydactyla",
                                       "Strix aluco","Buteo buteo","Phylloscopus sibilatrix","Poecile palustris",
                                       "Regulus regulus","Pyrrhula pyrrhula","Muscicapa striata","Picus viridis",
                                       "Dendrocopos major","Turdus viscivorus","Accipiter gentilis","Dryobates minor",
                                       "Loxia curvirostra","Lophophanes cristatus","Poecile montanus","Corvus corax",
                                       "Spinus spinus","Accipiter nisus","Fringilla coelebs","Regulus ignicapilla",
                                       "Oriolus oriolus","Periparus ater","Dryocopus martius"))



species_net_forest <- merge(species_net_forest, species_all, by="name_long", all.x=T)
species_sub <- species_net_forest <- na.omit(species_net_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Netherlands",]
species_sub <- species_net_forest <- species_net_forest[species_net_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_net_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_net_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Norway
# https://onlinelibrary.wiley.com/doi/full/10.1111/ibi.12896
# https://www.naturindeks.no/Themes/28

species_nor_farm  <- data.frame(name_long=c("Vanellus vanellus","Numenius arquata","Larus canus","Columba palumbus",
                                     "Turdus pilaris","Pica pica","Corvus cornix","Alauda arvensis",
                                     "Hirundo rustica","Delichon urbicum","Curruca communis","Curruca curruca",
                                     "Sturnus vulgaris","Saxicola rubetra","Passer domesticus","Motacilla alba",
                                     "Emberiza citrinella"))

species_nor_farm <- merge(species_nor_farm, species_all, by="name_long", all.x=T)
species_sub <- species_nor_farm <- na.omit(species_nor_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Norway",]
species_sub <- species_nor_farm <- species_nor_farm[species_nor_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_nor_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_nor_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_nor_forest  <- data.frame(name_long=c("Turdus viscivorus","Poecile montanus","Accipiter gentilis","Sylvia atricapilla",
                                       "Ficedula hypoleuca","Lophophanes cristatus","Hippolais icterina","Lyrurus tetrix",
                                       "Dendrocopos major","Dryocopus martius","Picus viridis","Turdus iliacus",
                                       "Garrulus glandarius","Certhia familiaris","Phoenicurus phoenicurus","Regulus regulus",
                                       "Sylvia borin","Muscicapa striata"))


species_nor_forest <- merge(species_nor_forest, species_all, by="name_long", all.x=T)
species_sub <- species_nor_forest <- na.omit(species_nor_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Norway",]
species_sub <- species_nor_forest <- species_nor_forest[species_nor_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_forest[y_forest == 0] <- NA

dfa_nor_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_nor_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Poland
# https://sdg.gov.pl/en/statistics_nat/15-1-a/
# https://agro.icm.edu.pl/agro/element/bwmeta1.element.agro-7c2d9bb3-978e-4620-bd07-62df7e877dc2

species_pol_farm  <- data.frame(name_long=c("Ciconia ciconia", "Falco tinnunculus", "Vanellus vanellus", "Limosa limosa",
                                     "Upupa epops", "Streptopelia turtur", "Alauda arvensis", "Galerida cristata",
                                     "Anthus pratensis", "Motacilla flava", "Hirundo rustica", "Saxicola rubetra",
                                     "Saxicola rubicola", "Curruca communis", "Lanius collurio", "Passer montanus",
                                     "Sturnus vulgaris", "Linaria cannabina", "Serinus serinus", "Emberiza calandra",
                                     "Emberiza citrinella", "Emberiza hortulana"))

species_pol_farm <- merge(species_pol_farm, species_all, by="name_long", all.x=T)
species_sub <- species_pol_farm <- na.omit(species_pol_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Poland",]
species_sub <- species_pol_farm <- species_pol_farm[species_pol_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_pol_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_pol_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_pol_forest  <- data.frame(name_long=c("Periparus ater", "Lophophanus cristatus", "Regulus regulus", "Erithacus rubecula",
                                       "Pyrrhula pyrrhula", "Phyloscopus sibilatrix", "Dendrocopos major", "Troglodytes troglodytes",
                                       "Certhia familiaris", "Phylloscopus colybita", "Regulus ignicapillus", "Anthus trivialis", 
                                       "Turdus philomelos", "Dryocopus martius", "Garrulus glandarius", "Turdus viscivorus",
                                       "Fringilla coelebs", "Ficedula parva", "Ficedula hypoleuca", "Sitta europea", 
                                       "Prunella modularis", "Spinus spinus", "Phylloscopus trochulus", "Poecile montanus",
                                       "Parus major", "Poecile palustris", "Turdus merula", "Lullula arborea",
                                       "Sylvia atricapilla", "Certhia brachydactyla", "Aegithalos caudatus", "Coccothraustes coccothraustes",
                                       "Columba oenas", "Phoenicurus phoenicurus"))


species_pol_forest <- merge(species_pol_forest, species_all, by="name_long", all.x=T)
species_sub <- species_pol_forest <- na.omit(species_pol_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Poland",]
species_sub <- species_pol_forest <- species_pol_forest[species_pol_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_pol_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_pol_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Portugal
# https://www.spea.pt/wp-content/uploads/2021/06/relatorio_cac_2021_vf3.pdf

species_por_farm  <- data.frame(name_long=c("Athene noctua","Bubulcus ibis","Carduelis carduelis","Chloris chloris",
                                     "Ciconia ciconia","Cisticola juncidis","Coturnix coturnix","Delichon urbicum",
                                     "Emberiza cirlus","Falco tinnunculus","Galerida cristata","Hirundo rustica",
                                     "Lanius meridionalis","Linaria cannabina","Merops apiaster","Emberiza calandra",
                                     "Milvus migrans","Passer domesticus","Pica pica","Saxicola rubicola",
                                     "Serinus serinus","Sturnus unicolor","Upupa epops"))

species_por_farm <- merge(species_por_farm, species_all, by="name_long", all.x=T)
species_sub <- species_por_farm <- na.omit(species_por_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Portugal",]
species_sub <- species_por_farm <- species_por_farm[species_por_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_por_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_por_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_por_forest  <- data.frame(name_long=c("Aegithalos caudatus","Certhia brachydactyla","Columba palumbus","Cuculus canorus",
                                       "Cyanistes caeruleus","Dendrocopos major","Erithacus rubecula","Fringilla coelebs",
                                       "Garrulus glandarius","Lanius senator","Lophophanes cristatus","Lullula arborea",
                                       "Oriolus oriolus","Parus major","Periparus ater","Picus viridis",
                                       "Sitta europaea","Streptopelia turtur","Sylvia atricapilla","Troglodytes troglodytes"))


species_por_forest <- merge(species_por_forest, species_all, by="name_long", all.x=T)
species_sub <- species_por_forest <- na.omit(species_por_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Portugal",]
species_sub <- species_por_forest <- species_por_forest[species_por_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_forest[y_forest == 0] <- NA

dfa_por_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_por_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Republic of Ireland
# https://birdwatchireland.ie/our-work/surveys-research/research-surveys/countryside-bird-survey/countryside-bird-population-indicators/
# no forest indicator but asked

species_ire_farm  <- data.frame(name_long=c("Falco tinnunculus","Phasianus colchicus","Columba oenas","Columba palumbus",
                                     "Hirundo rustica","Motacilla alba","Saxicola rubicola","Pica pica",
                                     "Coloeus monedula","Corvus frugilegus","Corvus cornix","Sturnus vulgaris",
                                     "Passer domesticus","Fringilla coelebs","Chloris chloris","Carduelis carduelis",
                                     "Linaria cannabina","Emberiza citrinella"))

species_ire_farm <- merge(species_ire_farm, species_all, by="name_long", all.x=T)
species_sub <- species_ire_farm <- na.omit(species_ire_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Republic of Ireland",]
species_sub <- species_ire_farm <- species_ire_farm[species_ire_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_ire_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_ire_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

# Romania
# asked


# Slovakia
# https://www.enviroportal.sk/indicator/detail?id=4041&print=yes
# 18 species but asked

species_slk_farm  <- data.frame(name_long=c("Alauda arvensis", "Linaria cannabina", "Carduelis carduelis", "Emberiza calandra",
                                     "Emberiza citrinella", "Falco tinnunculus", "Hirundo rustica", "Chloris chloris",
                                     "Lanius collurio", "Locustella naevia", "Motacilla flava", "Passer montanus",
                                     "Saxicola rubetra", "Saxicola torquatus", "Serinus serinus", "Streptopelia turtur",
                                     "Sturnus vulgaris", "Curruca communis", " Curruca nisoria", "Vanellus vanellus"))

species_slk_farm <- merge(species_slk_farm, species_all, by="name_long", all.x=T)
species_sub <- species_slk_farm <- na.omit(species_slk_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Slovakia",]
species_sub <- species_slk_farm <- species_slk_farm[species_slk_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_slk_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_slk_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)




# Slovenia
# https://www.sciencedirect.com/science/article/pii/S0167880920303868
# no forest bird index (confirmed)

species_sln_farm  <- data.frame(name_long=c("Acrocephalus palustris","Alauda arvensis","Anthus trivialis","Linaria cannabina",
                                     "Carduelis carduelis","Columba oenas","Columba palumbus","Emberiza calandra",
                                     "Emberiza cirlus","Emberiza citrinella","Falco tinnunculus","Galerida cristata",
                                     "Hirundo rustica","Jynx torquilla","Lanius collurio","Lullula arborea",
                                     "Luscinia megarhynchos","Motacilla flava","Passer montanus","Phoenicurus phoenicurus",
                                     "Picus viridis","Saxicola rubetra","Saxicola torquatus","Serinus serinus",
                                     "Streptopelia turtur","Sturnus vulgaris","Curruca communis","Upupa epops",
                                     "Vanellus vanellus"))

species_sln_farm <- merge(species_sln_farm, species_all, by="name_long", all.x=T)
species_sub <- species_sln_farm <- na.omit(species_sln_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Slovenia",]
species_sub <- species_sln_farm <- species_sln_farm[species_sln_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_sln_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_sln_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

# Spain
# https://www.nature.com/articles/s41598-019-45854-0
# asked

species_spa_farm  <- data.frame(name_long=c("Merops apiaster","Upupa epops","Alauda arvensis","Melanocorypha calandra",
                                     "Falco tinnunculus","Cisticola juncidis","Coturnix coturnix","Galerida cristata",
                                     "Oenanthe hispanica","Emberiza calandra","Sturnus unicolor","Sturnus vulgaris",
                                     "Pterocles orientalis","Hirundo rustica","Passer domesticus","Passer montanus",
                                     "Passer hispaniolensis","Coloeus monedula","Carduelis carduelis","Athene noctua",
                                     "Linaria cannabina","Alectoris rufa","Tetrax tetrax","Calandrella brachydactyla",
                                     "Streptopelia turtur","Pica pica"))

species_spa_farm <- merge(species_spa_farm, species_all, by="name_long", all.x=T)
species_sub <- species_spa_farm <- na.omit(species_spa_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Spain",]
species_sub <- species_spa_farm <- species_spa_farm[species_spa_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_spa_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_spa_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

# Sweden

species_swe_farm  <- data.frame(name_long= c("Falco tinnunculus","Vanellus vanellus","Alauda arvensis","Hirundo rustica",
                                      "Corvus frugilegus","Saxicola rubetra","Curruca communis","Anthus pratensis",
                                      "Motacilla flava","Lanius collurio","Sturnus vulgaris","Linaria cannabina",
                                      "Emberiza citrinella","Emberiza hortulana","Passer montanus"))

species_swe_farm <- merge(species_swe_farm, species_all, by="name_long", all.x=T)
species_sub <- species_swe_farm <- na.omit(species_swe_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Sweden",]
species_sub <- species_swe_farm <- species_swe_farm[species_swe_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_swe_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_swe_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_swe_forest  <- data.frame(name_long=c("Accipiter nisus","Tetrastes bonasia","Tringa ochropus","Columba oenas",
                                       "Dendrocopos major","Dryocopus martius","Picus viridis","Jynx torquilla",
                                       "Dryobates minor","Picoides tridactylus","Nucifraga caryocatactes","Garrulus glandarius",
                                       "Periparus ater","Lophophanes cristatus","Poecile palustris","Poecile montanus",
                                       "Sitta europaea","Certhia familiaris","Turdus viscivorus","Phoenicurus phoenicurus",
                                       "Phylloscopus collybita","Phylloscopus sibilatrix","Regulus regulus","Ficedula hypoleuca",
                                       "Anthus trivialis","Coccothraustes coccothraustes","Spinus spinus","Pyrrhula pyrrhula",
                                       "Emberiza rustica"))


species_swe_forest <- merge(species_swe_forest, species_all, by="name_long", all.x=T)
species_sub <- species_swe_forest <- na.omit(species_swe_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Sweden",]
species_sub <- species_swe_forest <- species_swe_forest[species_swe_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_swe_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_swe_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)

# Switzerland
# https://www.researchgate.net/publication/228432513_Fur_welche_Vogelarten_tragt_die_Schweiz_eine_besondere_Verantwortung/link/09e415136eddea186c000000/download
# https://www.researchgate.net/profile/Verena-Keller/publication/289958618_Swiss_bird_index_SBIR_-_Combined_population_trend_indices_for_different_groups_of_regular_breeding_birds_in_Switzerland_1990-2004/links/589b161c92851c8bb685e99c/Swiss-bird-index-SBIR-Combined-population-trend-indices-for-different-groups-of-regular-breeding-birds-in-Switzerland-1990-2004.pdf?_sg%5B0%5D=started_experiment_milestone&origin=journalDetail

species_swi_farm  <- data.frame(name_long=c("Falco tinnunculus","Ciconia ciconia","Alectoris rufa","Perdix perdix",
                                     "Crex crex","Vanellus vanellus","Tyto alba","Otus scops",
                                     "Athene noctua","Asio otus","Upupa epops","Jynx torquilla",
                                     "Lullula arborea","Alauda arvensis","Anthus pratensis","Motacilla flava",
                                     "Phoenicurus phoenicurus","Saxicola rubetra","Saxicola torquatus","Curruca communis",
                                     "Lanius minor","Lanius excubitor","Lanius senator","Coloeus monedula",
                                     "Corvus frugilegus","Emberiza cirlus","Emberiza calandra","Milvus milvus",
                                     "Buteo buteo","Turdus pilaris","Corvus corone","Coturnix coturnix",
                                     "Streptopelia turtur","Hirundo rustica","Anthus trivialis","Lanius collurio",
                                     "Sturnus vulgaris","Passer montanus","Linaria cannabina","Emberiza citrinella",
                                     "Circus pygargus"," Curruca nisoria"))


species_swi_farm <- merge(species_swi_farm, species_all, by="name_long", all.x=T)
species_sub <- species_swi_farm <- na.omit(species_swi_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Switzerland",]
species_sub <- species_swi_farm <- species_swi_farm[species_swi_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

y_farm[y_farm == 0] <- NA

dfa_swi_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_swi_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_swi_forest  <- data.frame(name_long=c("Tetrastes bonasia","Picus canus","Pernis apivorus","Lyrurus tetrix",
                                       "Tetrao urogallus","Scolopax rusticola","Glaucidium passerinum","Caprimulgus europaeus",
                                       "Dendrocoptes medius","Luscinia megarhynchos","Phylloscopus sibilatrix","Phylloscopus trochilus",
                                       "Accipiter gentilis","Accipiter nisus","Aegolius funereus","Dryocopus martius",
                                       "Picoides tridactylus","Prunella modularis","Erithacus rubecula","Turdus torquatus",
                                       "Turdus merula","Turdus philomelos","Turdus viscivorus","Phylloscopus collybita",
                                       "Regulus regulus","Regulus ignicapillus","Poecile palustris","Lophophanes cristatus",
                                       "Periparus ater","Parus major","Certhia familiaris","Nucifraga caryocatactes",
                                       "Fringilla coelebs","Carduelis citrinella","Loxia curvirostra","Pyrrhula pyrrhula",
                                       "Columba oenas","Columba palumbus","Strix aluco","Dendrocopos major",
                                       "Dryobates minor","Troglodytes troglodytes","Curruca curruca","Sylvia borin",
                                       "Sylvia atricapilla","Phylloscopus bonelli","Ficedula hypoleuca","Aegithalos caudatus",
                                       "Poecile montanus","Cyanistes caeruleus","Sitta europaea","Certhia brachydactyla",
                                       "Oriolus oriolus","Garrulus glandarius","Spinus spinus","Acanthis flammea",
                                       "Coccothraustes coccothraustes","Ficedula albicollis"))

species_swi_forest <- merge(species_swi_forest, species_all, by="name_long", all.x=T)
species_sub <- species_swi_forest <- na.omit(species_swi_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "Switzerland",]
species_sub <- species_swi_forest <- species_swi_forest[species_swi_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_swi_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_swi_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)


# United Kingdom     
# https://www.rspb.org.uk/our-work/conservation/conservation-and-sustainability/farming/near-you/farmland-bird-indicator/
# and https://www.gov.uk/government/statistics/wild-bird-populations-in-the-uk

species_uk_farm  <- data.frame(name_long=c("Passer montanus", "Emberiza calandra", "Streptopelia turtur", "Perdix perdix",
                                    "Motacilla flava", "Sturnus vulgaris", "Linaria cannabina", "Vanellus vanellus", 
                                    "Emberiza citrinella", "Alauda arvensis", "Falco tinnunculus", "Emberiza schoeniclus",
                                    "Curruca communis", "Chloris chloris", "Corvus frugilegus", "Columba oenas",  
                                    "Carduelis carduelis", "Columba palumbus", "Coloeus monedula"))

species_uk_farm <- merge(species_uk_farm, species_all, by="name_long", all.x=T)
species_sub <- species_uk_farm <- na.omit(species_uk_farm)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "United Kingdom",]
species_sub <- species_uk_farm <- species_uk_farm[species_uk_farm$code_sp %in% unique(Obs$code_sp),]
y_farm <- dcast(Obs[,c("code_sp","Index","Year")],
                code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_farm <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                     code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_farm[y_farm == 0] <- NA

dfa_uk_farm <- make_dfa(data_ts = y_farm,data_ts_se = obs_se_farm,
                         species_sub = species_uk_farm,nfac = 0,
                         mintrend = 1,maxtrend = 5,AIC = TRUE,
                         nboot = 500,silent = TRUE,control = list(),
                         se_log = FALSE,is_mean_centred = FALSE)

species_uk_forest  <- data.frame(name_long=c("Turdus merula","Cyanistes caeruleus","Pyrrhula pyrrhula","Fringilla coelebs",
                                    "Prunella modularis","Parus major","Curruca curruca","Aegithalos caudatus",
                                    "Erithacus rubecula","Turdus philomelos","Strix aluco","Troglodytes troglodytes",
                                    "Sylvia atricapilla","Phylloscopus collybita","Periparus ater","Sylvia borin",
                                    "Regulus regulus","Dendrocopos major","Picus viridis","Garrulus glandarius",
                                    "Dryobates minor","Poecile palustris","Luscinia megarhynchos","Sitta europaea",
                                    "Carduelis cabaret","Phoenicurus phoenicurus","Accipiter nisus","Muscicapa striata",
                                    "Anthus trivialis","Certhia familiaris","Poecile montana","Phylloscopus trochilus",
                                    "Ficedula hypoleuca","Phylloscopus sibilatrix","Loxia curvirostra","Spinus spinus",
                                    "Tetrao urogallus"))

species_uk_forest <- merge(species_uk_forest, species_all, by="name_long", all.x=T)
species_sub <- species_uk_forest <- na.omit(species_uk_forest)

Obs <- df_all_country_2000[df_all_country_2000$code_sp %in% species_sub$code_sp &
                        df_all_country_2000$CountryGroup == "United Kingdom",]
species_sub <- species_uk_forest <- species_uk_forest[species_uk_forest$code_sp %in% unique(Obs$code_sp),]
y_forest <- dcast(Obs[,c("code_sp","Index","Year")],
                  code_sp~Year, fun.aggregate = sum, value.var = "Index")
obs_se_forest <- dcast(Obs[,c("code_sp","Index_SE","Year")],
                       code_sp~Year, fun.aggregate = sum, value.var = "Index_SE")

#y_forest[y_forest == 0] <- NA

dfa_uk_forest <- make_dfa(data_ts = y_forest,data_ts_se = obs_se_forest,
                           species_sub = species_uk_forest,nfac = 0,
                           mintrend = 1,maxtrend = 5,AIC = TRUE,
                           nboot = 500,silent = TRUE,control = list(),
                           se_log = FALSE,is_mean_centred = FALSE)


# Species indices

## Species specialisation index (SSI)
# https://onlinelibrary.wiley.com/doi/10.1111/oik.02276

SSI <- read.csv("raw_data/SSI.csv", header=T)
SSI$Species <- paste0(SSI$Genus, sep=" ", SSI$Species)
SSI$Genus <- NULL

SSI_cornix <- SSI[SSI$Species=="Corvus corone",]
SSI_cornix$Species <- "Corvus cornix"
SSI <- rbind(SSI,SSI_cornix)
SSI$Species[SSI$Species=="Bonasa bonasia"] <- "Tetrastes bonasia"
SSI$Species[SSI$Species=="Carduelis cannabina"] <- "Linaria cannabina"
SSI$Species[SSI$Species=="Carduelis chloris"] <- "Chloris chloris"
SSI$Species[SSI$Species=="Carduelis flammea"] <- "Acanthis flammea"
SSI$Species[SSI$Species=="Carduelis spinus"] <- "Spinus spinus"
SSI$Species[SSI$Species=="Corvus monedula"] <- "Coloeus monedula"
SSI$Species[SSI$Species=="Delichon urbica"] <- "Delichon urbicum"
SSI$Species[SSI$Species=="Dendrocopos medius"] <- "Dendrocoptes medius"
SSI$Species[SSI$Species=="Dendrocopos minor"] <- "Dryobates minor"
SSI$Species[SSI$Species=="Hippolais pallida"] <- "Iduna pallida"
SSI$Species[SSI$Species=="Hirundo daurica"] <- "Cecropis daurica"
SSI$Species[SSI$Species=="Hirundo rupestris"] <- "Ptyonoprogne rupestris"
SSI$Species[SSI$Species=="Larus ridibundus"] <- "Chroicocephalus ridibundus"
SSI$Species[SSI$Species=="Miliaria calandra"] <- "Emberiza calandra"
SSI$Species[SSI$Species=="Parus ater"] <- "Periparus ater"
SSI$Species[SSI$Species=="Parus caeruleus"] <- "Cyanistes caeruleus"
SSI$Species[SSI$Species=="Parus cristatus"] <- "Lophophanes cristatus"
SSI$Species[SSI$Species=="Parus montanus"] <- "Poecile montanus"
SSI$Species[SSI$Species=="Parus palustris"] <- "Poecile palustris"
SSI$Species[SSI$Species=="Saxicola torquata"] <- "Saxicola torquatus"
SSI$Species[SSI$Species=="Serinus citrinella"] <- "Carduelis citrinella"
SSI$Species[SSI$Species=="Sylvia cantillans"] <- "Curruca cantillans"
SSI$Species[SSI$Species=="Sylvia communis"] <- "Curruca communis"
SSI$Species[SSI$Species=="Sylvia curruca"] <- "Curruca curruca"
SSI$Species[SSI$Species=="Sylvia hortensis"] <- "Curruca hortensis"
SSI$Species[SSI$Species=="Sylvia melanocephala"] <- "Curruca melanocephala"
SSI$Species[SSI$Species=="Sylvia melanothorax"] <- "Curruca melanothorax"
SSI$Species[SSI$Species=="Sylvia nisoria"] <- "Curruca nisoria"
SSI$Species[SSI$Species=="Sylvia undata"] <- "Curruca undata"
SSI$Species[SSI$Species=="Tetrao tetrix"] <- "Lyrurus tetrix"

SSI$Species[SSI$Species==""] <- "Acrocephalus arundinaceus"
SSI$Species[SSI$Species==""] <- "Acrocephalus schoenobaenus"
SSI$Species[SSI$Species==""] <- "Acrocephalus scirpaceus"
SSI$Species[SSI$Species==""] <- "Actitis hypoleucos"
SSI$Species[SSI$Species==""] <- "Alcedo atthis"
SSI$Species[SSI$Species==""] <- "Anas platyrhynchos"
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""
SSI$Species[SSI$Species==""] <- ""

## Species temperature index (STI)
# https://royalsocietypublishing.org/doi/epdf/10.1098/rspb.2008.0878
# https://onlinelibrary.wiley.com/doi/epdf/10.1111/j.1600-0587.2012.07799.x
# get bird occurence from https://ebba2.info/data-request/
# get climat data from https://surfobs.climate.copernicus.eu/dataaccess/access_eobs.php

# map of temperature
require(raster)
temp_se <- brick("raw_data/tg_ens_mean_0.1deg_reg_v23.1e.nc")

end_year <- 0

for(i in 1:(2021-1950)){
  print(i)
  beg_year <- end_year+1
  end_year <- end_year+365
  year <- paste0("temp_",1949+i)
  if(i %in% c(seq(1952,2021,4)-1949)){
    end_year<-end_year+1 # account for bisextil years
  }
  map_year <- assign(year, mean(temp_se[[(beg_year+31):(beg_year+201)]], na.rm=T)) # from February to July cf Loxia
  path <- paste0("output/temp_",1949+i,".tif")
  writeRaster(x=map_year, filename=path, overwrite=T)
}

# species occurrence

species_occ <- read.csv("raw_data/ebba2_data_occurrence_50km.csv", sep=";", header = T)

# grid of occurrence

grid_occ <- readOGR("raw_data/ebba2_grid50x50_v1.shp")
grid_occ_poly <- SpatialPolygons(grid_occ@polygons)

# temperature by grid cell and by year

for(i in 2000:2017){ # data from the same period as the survey (2000-2017)
  print(i)
  year <- paste0("temp_",i)
  path <- paste0("output/temp_",i,".tif")
  temp_year <- paste0("temp_",i,"_tr")
  value_year <- paste0("value_",i)
  mean_value_year <- paste0("mean_value_",i)
  
  map_year <- assign(year, raster(path))
  map_year_tr <- assign(temp_year, projectRaster(map_year, crs = proj4string(grid_occ)))
  
  value_year_list <- assign(value_year, extract(map_year_tr,grid_occ_poly))
  mean_value_year_list <- assign(mean_value_year, lapply(value_year_list, FUN = function(x){return(mean(x,na.rm=T))}))
  
  grid_occ@data[,ncol(grid_occ@data)+1] <- unlist(mean_value_year_list)
  names(grid_occ@data)[ncol(grid_occ@data)] <- year
  
}

grid_occ_fort <- fortify(grid_occ, region='cell50x50')
grid_occ_fort <- merge(grid_occ_fort, grid_occ@data,
                       by.x = "id", by.y="cell50x50")

ggplot(grid_occ_fort, aes(x = long, y = lat, group = id)) + 
  geom_polygon(aes(fill=temp_2000))+
  scale_fill_gradient2()

grid_occ_fort$grid_occ_mean <- apply(grid_occ_fort,1,function(x){return(mean(as.numeric(x[8:25]), na.rm=T))})

grid_mean <- data.frame(grid_occ_fort %>% group_by(id) %>% summarise(mean_temp = mean(grid_occ_mean, na.rm=T)))

# Merge species occurrence and temperature

species_occ_temp <- merge(species_occ,grid_mean,
                          by.x = "cell50x50", by.y = "id", all.x=T)

species_sti <- data.frame(species_occ_temp %>% group_by(birdlife_scientific_name) %>% summarise(STI = mean(mean_temp, na.rm=T),
                                                                                                sd_STI = sd(mean_temp, na.rm=T)))
STI <- species_sti

names(STI)[1] <- "Species"

STI$Species <- as.character(STI$Species)
STI <- STI[which(STI$Species!=""),]
STI$Species[STI$Species=="Bonasa bonasia"] <- "Tetrastes bonasia"
STI_cornix <- STI[STI$Species=="Corvus corone",]
STI_cornix$Species <- "Corvus cornix"
STI <- rbind(STI,STI_cornix)
STI$Species[STI$Species=="Corvus monedula"] <- "Coloeus monedula"
STI$Species[STI$Species=="Cyanecula svecica"] <- "Luscinia svecica"
STI$Species[STI$Species=="Larus ridibundus"] <- "Chroicocephalus ridibundus"
STI$Species[STI$Species=="Leiopicus medius"] <- "Dendrocoptes medius"
STI$Species[STI$Species=="Sylvia cantillans"] <- "Curruca cantillans"
STI$Species[STI$Species=="Sylvia communis"] <- "Curruca communis"
STI$Species[STI$Species=="Sylvia curruca"] <- "Curruca curruca"
STI$Species[STI$Species=="Sylvia hortensis"] <- "Curruca hortensis"
STI$Species[STI$Species=="Sylvia melanocephala"] <- "Curruca melanocephala"
STI$Species[STI$Species=="Sylvia melanothorax"] <- "Curruca melanothorax"
STI$Species[STI$Species=="Sylvia nisoria"] <- "Curruca nisoria"
STI$Species[STI$Species=="Sylvia undata"] <- "Curruca undata"


## Species functional index (SFI)
# https://onlinelibrary.wiley.com/doi/pdf/10.1111/j.1461-0248.2010.01493.x
# https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12266
# https://onlinelibrary.wiley.com/doi/full/10.1111/geb.12709
# https://www.jstor.org/stable/4539206#metadata_info_tab_contents

library(ade4)
library(RVAideMemoire)
library(cluster)
library(geometry)

trait <- read.table("raw_data/Life-history characteristics of European birds.txt", header = T, sep = "\t")

# Select traits

trait_selected <- trait[,c("Species",
                           "WeightU_MEAN","WingU_MEAN","BillU_MEAN","TarsusU_MEAN","TailU_MEAN","Sexual.dimorphism","Clutch_MEAN","Life.span","Age.of.independence", "Age.of.first.breeding", "Fledging.period","Egg_MASS", # Resource quantity
                           "Fish_B","Other.vertebrates_B","Carrion_B", # Vertebrate diet to aggregate
                           "Arthropods_B","Other.invertebrates_B", # Invertebrate diet to aggregate
                           "Folivore_B","Frugivore_B","Granivore_B", # Plant diet to aggregate
                           "Nest.type","Nest.building", # Nest location
                           "Association.during.nesting", "Mating.system", "Hatching", "Association.outside.the.breeding.season", "Territoriality", "Human.settlements", "Incubation.sex","Young", # Behaviour
                           "Sedentary","Short.distance.migrant","Long.distance.migrant" # Migratory status to aggregate
                           )]

trait_selected$Vertebrate_diet <- 0
trait_selected$Vertebrate_diet[trait_selected$Fish_B == 1 |
                               trait_selected$Other.vertebrates_B == 1 |
                               trait_selected$Carrion_B == 1] <- 1
trait_selected$Fish_B <- trait_selected$Other.vertebrates_B <- trait_selected$Carrion_B <- NULL

trait_selected$Invertebrate_diet <- 0
trait_selected$Invertebrate_diet[trait_selected$Arthropods_B == 1 |
                                 trait_selected$Other.invertebrates_B == 1] <- 1
trait_selected$Arthropods_B <- trait_selected$Other.invertebrates_B <- NULL

trait_selected$Plant_diet <- 0
trait_selected$Plant_diet[trait_selected$Folivore_B == 1 |
                                   trait_selected$Frugivore_B == 1 |
                                   trait_selected$Granivore_B == 1] <- 1
trait_selected$Folivore_B <- trait_selected$Frugivore_B <- trait_selected$Granivore_B  <- NULL

trait_selected$Migratory_status <- "S"
trait_selected$Migratory_status[trait_selected$Sedentary == 1 &
                                trait_selected$Short.distance.migrant == 1] <- "S,SD"
trait_selected$Migratory_status[trait_selected$Short.distance.migrant == 1 &
                                  trait_selected$Long.distance.migrant == 1] <- "SD,LD"
trait_selected$Migratory_status[trait_selected$Short.distance.migrant == 0 &
                                  trait_selected$Long.distance.migrant == 1] <- "LD"
trait_selected$Sedentary <- trait_selected$Short.distance.migrant <- trait_selected$Long.distance.migrant  <- NULL

# Correct species names

trait_selected$Species <- as.character(trait_selected$Species)
trait_selected <- trait_selected[which(trait_selected$Species!=""),]
trait_selected$Species[trait_selected$Species=="Bonasa bonasia"] <- "Tetrastes bonasia"
trait_selected_cornix <- trait_selected[trait_selected$Species=="Corvus corone",]
trait_selected_cornix$Species <- "Corvus cornix"
trait_selected <- rbind(trait_selected,trait_selected_cornix)
trait_selected$Species[trait_selected$Species=="Corvus monedula"] <- "Coloeus monedula"
trait_selected$Species[trait_selected$Species=="Cyanecula svecica"] <- "Luscinia svecica"
trait_selected$Species[trait_selected$Species=="Larus ridibundus"] <- "Chroicocephalus ridibundus"
trait_selected$Species[trait_selected$Species=="Leiopicus medius"] <- "Dendrocoptes medius"
trait_selected$Species[trait_selected$Species=="Sylvia cantillans"] <- "Curruca cantillans"
trait_selected$Species[trait_selected$Species=="Sylvia communis"] <- "Curruca communis"
trait_selected$Species[trait_selected$Species=="Sylvia curruca"] <- "Curruca curruca"
trait_selected$Species[trait_selected$Species=="Sylvia hortensis"] <- "Curruca hortensis"
trait_selected$Species[trait_selected$Species=="Sylvia melanocephala"] <- "Curruca melanocephala"
trait_selected$Species[trait_selected$Species=="Sylvia melanothorax"] <- "Curruca melanothorax"
trait_selected$Species[trait_selected$Species=="Sylvia nisoria"] <- "Curruca nisoria"
trait_selected$Species[trait_selected$Species=="Sylvia undata"] <- "Curruca undata"

trait_selected_sp <- trait_selected[which(trait_selected$Species %in% species_all$name_long),]
trait_selected_sp[,c("Sexual.dimorphism","Nest.type", "Nest.building", "Association.during.nesting",
          "Mating.system", "Hatching", "Association.outside.the.breeding.season",
          "Territoriality", "Human.settlements", "Incubation.sex",
          "Young", "Vertebrate_diet", "Invertebrate_diet",
          "Plant_diet", "Migratory_status")] <- lapply(trait_selected_sp[,c("Sexual.dimorphism","Nest.type", "Nest.building", "Association.during.nesting",
                                                                            "Mating.system", "Hatching", "Association.outside.the.breeding.season",
                                                                            "Territoriality", "Human.settlements", "Incubation.sex",
                                                                            "Young", "Vertebrate_diet", "Invertebrate_diet",
                                                                            "Plant_diet", "Migratory_status")], factor)

table.na <- apply(trait_selected_sp,2,function(x){sum(is.na(x))})

# 1st possibility: Hill&Smith

# Replacing NA by mean and scale

trait_selected_sp[,c(2:6,8:13)] <- apply(trait_selected_sp[,c(2:6,8:13)],2,function(x){x[is.na(x)] <- mean(x,na.rm=T)})
for(i in c(7,14:27)){
  trait_selected_sp[is.na(trait_selected_sp[,i]),i] <- names(which.max(table(trait_selected_sp[,i])))
}


trait_selected_sp[,c(2:6,8:13)] <- apply(trait_selected_sp[,c(2:6,8:13)],2,function(x){return(scale(x))})

# Explanatory mixed analysis 

AMix <- dudi.hillsmith(trait_selected_sp[,-1], scannf=F, nf=25) # 25 axes for 75%
MVA.synt(AMix, rows=25)
MVA.plot(AMix)
MVA.plot(AMix,"corr")
scat.cr(AMix,axis=1)

trait_hs <- AMix$li
row.names(trait_hs) <- trait_selected_sp[,1]

mat.dista <- daisy(trait_hs, metric="euclidean")
names(mat.dista) <- trait_selected_sp$Species

hc1 <- hclust(mat.dista)
hc2 <- hclust(cophenetic(hc1), method = "ave")
plot(hc2)

vect.dista <- data.frame(Species=trait_selected_sp$Species,
                         fd.obs=rowSums(as.matrix(mat.dista)),
                         fd.obs.coph=rowSums(as.matrix(cophenetic(hc1))))

vect.dista$SFI <- scale(vect.dista$fd.obs.coph, center=F)


# 2nd possibility: Gower distance matrix

mat.dist <- daisy(trait_selected_sp[,-1], metric="gower", type=list(asymm=c("Human.settlements","Vertebrate_diet","Invertebrate_diet","Plant_diet")))
mat.dist2 <- sqrt(mat.dist) # euclidean matrix
names(mat.dist2) <- trait_selected_sp$Species
vect.dist <- data.frame(Species=trait_selected_sp$Species, fd.obs=rowSums(as.matrix(mat.dist2)))

hc1 <- hclust(mat.dist2)
hc2 <- hclust(cophenetic(hc1), method = "ave")
plot(hc2)

vect.distb <- data.frame(Species=trait_selected_sp$Species,
                         fd.obs=rowSums(as.matrix(mat.dist2)),
                         fd.obs.coph=rowSums(as.matrix(cophenetic(hc1))))

vect.distb$SFI <- scale(vect.distb$fd.obs.coph, center=F)

SFI <- merge(vect.dista[,c("Species","SFI")],vect.distb[,c("Species","SFI")],
             by="Species")

SXI <- merge(SFI,STI,by="Species", all.x=T)
SXI <- merge(SXI,SSI,by="Species", all.x=T)

write.csv(SXI,"output/SXI.csv", row.names = F)

### Test model to explain clusters

data_mod <- merge(dfa_aus_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_aus_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_aus_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_aus_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))


data_mod <- merge(dfa_bel_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_bel_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_bel_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_bel_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_cze_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_cze_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_cze_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_cze_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_den_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_den_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_den_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_den_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_est_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_est_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_est_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_est_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_fin_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_fin_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_fin_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_fin_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_fra_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_fra_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_fra_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_fra_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_ger_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_ger_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_ger_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_ger_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_hun_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_hun_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_hun_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_hun_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_ire_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_ire_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_ita_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_ita_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_ita_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_ita_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_lat_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_lat_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_lit_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_lit_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_net_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_net_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_net_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_net_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_nor_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_nor_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_nor_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_nor_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_pol_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_pol_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_pol_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_pol_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_spa_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_spa_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_swe_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_swe_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_swe_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_swe_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_swi_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_swi_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_swi_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_swi_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

data_mod <- merge(dfa_uk_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_uk_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC1,data_mod$SFI.y)
cor.test(data_mod$PC1,data_mod$STI)
cor.test(data_mod$PC1,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC1","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_uk_farm$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))
data_mod <- merge(dfa_uk_forest$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))
cor.test(data_mod$PC2,data_mod$SFI.y)
cor.test(data_mod$PC2,data_mod$STI)
cor.test(data_mod$PC2,data_mod$SSI)
pcor(na.omit(data_mod[,c("PC2","SFI.y","STI","SSI")]))

# Aggragate results

list_dfa <- list(dfa_aus_farm,dfa_aus_forest,dfa_bel_farm,dfa_bel_forest,
                 dfa_cze_farm,dfa_cze_forest,dfa_den_farm,dfa_den_forest,
                 dfa_est_farm,dfa_est_forest,dfa_fin_farm,dfa_fin_forest,
                 dfa_fra_farm,dfa_fra_forest,dfa_ger_farm,dfa_ger_forest,
                 dfa_hun_farm,dfa_hun_forest,dfa_ire_farm,
                 dfa_ita_farm,dfa_ita_forest,dfa_lat_farm,
                 dfa_lit_farm,dfa_net_farm,dfa_net_forest,
                 dfa_nor_farm,dfa_nor_forest,dfa_pol_farm,dfa_pol_forest,
                 dfa_spa_farm,dfa_swe_farm,dfa_swe_forest,
                 dfa_swi_farm,dfa_swi_forest,dfa_uk_farm,dfa_uk_forest)


result_cor <- data.frame(Country = c(rep("Austria",2),rep("Belgium",2),rep("Czechia",2),rep("Denmark",2),
                                     rep("Estonia",2),rep("Finland",2),rep("France",2),rep("Germany",2),
                                     rep("Hungary",2),rep("Ireland",1),rep("Italia",2),rep("Latvia",1),
                                     rep("Lithuania",1),rep("Netherlands",2),rep("Norway",2),rep("Poland",2),
                                     rep("Spain",1),rep("Sweden",2),rep("Switzerland",2),rep("UK",2)),
                         Index = c(rep(c("FBI","WBI"),9),"FBI",rep(c("FBI","WBI"),1),
                                   "FBI","FBI",rep(c("FBI","WBI"),3),"FBI",rep(c("FBI","WBI"),3)),
                         Nb_lat_trend = NA, Nb_cluster = NA,
                         PCA1_SFI = NA, PCA1_SFI_pval = NA, PCA1_STI = NA, PCA1_STI_pval = NA,
                         PCA1_SSI = NA, PCA1_SSI_pval = NA, R2_PCA1 = NA,
                         PCA2_SFI = NA, PCA2_SFI_pval = NA, PCA2_STI = NA, PCA2_STI_pval = NA,
                         PCA2_SSI = NA, PCA2_SSI_pval = NA, R2_PCA2 = NA)

for(i in 1:length(list_dfa)){
  
  data_dfa <- list_dfa[[i]]
  
  result_cor[i,3] <- length(unique(data_dfa$data_loadings$variable))
  result_cor[i,4] <- length(unique(data_dfa$group[[1]][[1]]$group))
  
  data_mod <- merge(data_dfa$group[[1]][[1]],SXI, by.x="name_long", by.y="Species")
  
  result_cor[i,5] <- cor.test(data_mod$PC1,data_mod$SFI.y)$estimate
  result_cor[i,6] <- cor.test(data_mod$PC1,data_mod$SFI.y)$p.value
  result_cor[i,7] <- cor.test(data_mod$PC1,data_mod$STI)$estimate
  result_cor[i,8] <- cor.test(data_mod$PC1,data_mod$STI)$p.value
  result_cor[i,9] <- cor.test(data_mod$PC1,data_mod$SSI)$estimate
  result_cor[i,10] <- cor.test(data_mod$PC1,data_mod$SSI)$p.value
  result_cor[i,11] <- summary(lm(PC1~SFI.y+STI+SSI, data=data_mod))$r.squared

  result_cor[i,12] <- cor.test(data_mod$PC2,data_mod$SFI.y)$estimate
  result_cor[i,13] <- cor.test(data_mod$PC2,data_mod$SFI.y)$p.value
  result_cor[i,14] <- cor.test(data_mod$PC2,data_mod$STI)$estimate
  result_cor[i,15] <- cor.test(data_mod$PC2,data_mod$STI)$p.value
  result_cor[i,16] <- cor.test(data_mod$PC2,data_mod$SSI)$estimate
  result_cor[i,17] <- cor.test(data_mod$PC2,data_mod$SSI)$p.value
  result_cor[i,18] <- summary(lm(PC2~SFI.y+STI+SSI, data=data_mod))$r.squared
  
}

# Plot results

result_to_plot <- melt(result_cor, id.vars = c("Country","Index"),
                       measure.vars = c("PCA1_SFI","PCA1_STI","PCA1_SSI","PCA2_SFI","PCA2_STI","PCA2_SSI"))

result_to_plot$country_index <- paste0(result_to_plot$Country,result_to_plot$Index)

ggplot(data = result_to_plot, aes(y=country_index, x=variable, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal() + geom_text(aes(label=round(value,2)))

# Do the same by species
