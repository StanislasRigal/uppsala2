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
df$Species[df$Species=="Carduelis cannabina"] <- "Linaria cannabina"
df$Species[df$Species=="Carduelis chloris"] <- "Chloris chloris"
df$Species[df$Species=="Carduelis flammea"] <- "Acanthis flammea"
df$Species[df$Species=="Carduelis spinus"] <- "Spinus spinus"
df$Species[df$Species=="Corvus corone+cornix"] <- "Corvus corone"
df$Species[df$Species=="Delichon urbica"] <- "Delichon urbicum"
df$Species[df$Species=="Dendrocopos medius"] <- "Dendrocoptes medius"
df$Species[df$Species=="Dendrocopos minor"] <- "Dryobates minor"
df$Species[df$Species=="Hippolais pallida"] <- "Iduna pallida"
df$Species[df$Species=="Hirundo daurica"] <- "Cecropis daurica"
df$Species[df$Species=="Hirundo rupestris"] <- "Ptyonoprogne rupestris"
df$Species[df$Species=="Miliaria calandra"] <- "Emberiza calandra"
df$Species[df$Species=="Parus ater"] <- "Periparus ater"
df$Species[df$Species=="Parus caeruleus"] <- "Cyanistes caeruleus"
df$Species[df$Species=="Parus cristatus"] <- "Lophophanes cristatus"
df$Species[df$Species=="Parus montanus"] <- "Poecile montanus"
df$Species[df$Species=="Parus palustris"] <- "Poecile palustris"
df$Species[df$Species=="Saxicola torquata"] <- "Saxicola torquatus"
df$Species[df$Species=="Serinus citrinella"] <- "Carduelis citrinella"
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

df_all_country <- rbind.fill(droplevels(df[!(df$CountryGroup %in% c("Belgium-Brussels", "Belgium-Wallonia",
                                                                    "Germany East", "Germany West")),]),
                             data.frame(df_bel[,c("Species","Year","Index","Index_SE")],CountryGroup="Belgium"),
                             data.frame(df_ger[,c("Species","Year","Index","Index_SE")],CountryGroup="Germany"))


# List farmland and woodland species by country

# Austria
# https://info.bml.gv.at/themen/landwirtschaft/eu-agrarpolitik-foerderungen/laendl_entwicklung/le-07-13/evaluierung/le_studien/FBI.html
# https://info.bml.gv.at/themen/landwirtschaft/eu-agrarpolitik-foerderungen/laendl_entwicklung/le-07-13/evaluierung/le_studien/Waldvogelindikator.html


species_sub <- species_aus_farm <- c("Falco tinnunculus","Perdix perdix","Vanellus vanellus","Streptopelia turtur",
                                     "Jynx torquilla","Lullula arborea","Alauda arvensis","Anthus trivialis",
                                     "Anthus spinoletta","Saxicola rubetra","Saxicola rubicola","Oenanthe oenanthe",
                                     "Turdus pilaris","Acrocephalus palustris","Sylvia communis","Lanius collurio",
                                     "Sturnus vulgaris","Passer montanus","Serinus serinus","Carduelis citrinella",
                                     "Carduelis carduelis","Linaria cannabina","Emberiza citrinella","Emberiza calandra")

species_sub <- species_aus_forest <- c("Bonasa bonasia","Tetrao tetrix","Tetrao urogallus","Ciconia nigra",
                                       "Pernis apivorus","Accipiter gentilis","Accipiter nisus","Buteo buteo",
                                       "Scolopax rusticola","Columba oenas","Columba palumbus","Cuculus canorus",
                                       "Bubo bubo","Glaucidium passerinum","Strix aluco","Aegolius funereus",
                                       "Caprimulgus europaeus","Jynx torquilla","Picus canus","Picus viridis",
                                       "Dryocopus martius","Dendrocopos major","Dendrocopos medius","Dendrocopos leucotos",
                                       "Dendrocopos minor","Picoides tridactylus","Lullula arborea","Anthus trivialis",
                                       "Troglodytes troglodytes","Prunella modularis","Erithacus rubecula","Luscinia megarhynchos",
                                       "Phoenicurus phoenicurus","Turdus torquatus","Turdus merula","Turdus pilaris",
                                       "Turdus philomelos","Turdus viscivorus","Locustella fluviatilis","Hippolais icterina",
                                       "Sylvia curruca","Sylvia borin","Sylvia atricapilla","Phylloscopus bonelli",
                                       "Phylloscopus sibilatrix","Phylloscopus collybita","Phylloscopus trochilus","Regulus regulus",
                                       "Regulus ignicapilla","Muscicapa striata","Ficedula parva","Ficedula albicollis",
                                       "Ficedula hypoleuca","Aegithalos caudatus","Parus palustris","Parus montanus",
                                       "Parus cristatus","Parus ater","Parus caeruleus","Parus major",
                                       "Sitta europaea","Certhia familiaris","Certhia brachydactyla","Oriolus oriolus",
                                       "Garrulus glandarius","Nucifraga caryocatactes","Fringilla coelebs","Serinus serinus",
                                       "Carduelis chloris","Carduelis spinus","Carduelis flammea","Loxia curvirostra",
                                       "Pyrrhula pyrrhula","Coccothraustes coccothraustes")


# Belgium
# from https://www.aves.be/fileadmin/Aves/Bulletins/Articles/47_1/47_1_1.pdf (Tab1 minus species with only presence absence (x in Frequency))
# and classification from PECBMS (https://pecbms.info/methods/pecbms-methods/3-multispecies-indicators/species-selection-and-classification/  Atlantic)

species_sub <- species_bel_farm <- c("Alauda arvensis","Anthus pratensis","Emberiza citrinella","Hippolais polyglotta",
                                     "Lanius collurio","Motacilla flava","Passer montanus","Perdix perdix",
                                     "Saxicola torquatus","Streptopelia turtur","Sylvia communis","Sylvia curruca",
                                     "Vanellus vanellus")

species_sub <- species_bel_forest <- c("Certhia brachydactyla","Certhia familiaris","Coccothraustes coccothraustes","Dendrocopos major",
                                       "Dryobates minor","Dryocopus martius","Erithacus rubecula","Fringilla coelebs",
                                       "Garrulus glandarius","Leiopicus medius","Lophophanes cristatus","Oriolus oriolus",
                                       "Periparus ater","Phoenicurus phoenicurus","Phylloscopus collybita","Phylloscopus sibilatrix",
                                       "Poecile montanus","Poecile palustris","Pyrrhula pyrrhula","Regulus ignicapilla",
                                       "Regulus regulus","Sitta europaea","Sylvia borin","Turdus viscivorus")

# Bulgaria
# https://bspb.org/en/about-birds/common-bird-monitoring/

species_sub <- species_bul_farm <- c("Coturnix coturnix","Galerida cristata","Alauda arvensis","Motacilla flava",
                                     "Sylvia communis","Lanius collurio","Sturnus vulgaris","Emberiza hortulana",
                                     "Emberiza melanocephala","Miliaria calandra")

species_sub <- species_bul_forest <- c("Columba palumbus","Dendrocopos major","Erithacus rubecula","Fringilla coelebs",
                                       "Garrulus glandarius","Oriolus oriolus","Parus major","Sylvia atricapilla",
                                       "Turdus philomelos")

# Cyprus
# asked

# Czech Republic
# https://www.tandfonline.com/doi/abs/10.1080/00063657.2015.1048423
# https://www.tandfonline.com/doi/abs/10.1080/00063650709461481

species_sub <- species_cze_farm <- c("Hirundo rustica","Corvus corone","Sylvia communis","Falco tinnunculus",
                                     "Carduelis cannabina","Sturnus vulgaris","Pica pica","Carduelis carduelis",
                                     "Carduelis chloris","Corvus monedula","Vanellus vanellus","Lanius collurio",
                                     "Emberiza schoeniclus","Alauda arvensis","Passer montanus","Streptopelia turtur",
                                     "Saxicola rubetra","Columba palumbus","Emberiza citrinella")

species_sub <- species_cze_forest <- c("Buteo buteo","Accipiter nisus","Columba palumbus","Streptopelia turtur",
                                     "Cuculus canorus","Picus viridis","Picus canus","Dryocopus martius",
                                     "Dendrocopos major","Dendrocopos medius","Dendrocopos minor","Jynx torquilla",
                                     "Anthus trivialis","Troglodytes troglodytes","Prunella modularis","Erithacus rubecula",
                                     "Luscinia megarhynchos","Phoenicurus phoenicurus","Turdus merula","Turdus philomelos",
                                     "Turdus viscivorus","Hippolais icterina","Sylvia borin","Sylvia atricapilla",
                                     "Phylloscopus sibilatrix","Phylloscopus collybita","Phylloscopus trochilus","Regulus regulus",
                                     "Muscicapa striata","Ficedula hypoleuca","Ficedula albicollis","Aegithalos caudatus",
                                     "Poecile palustris","Poecile montanus","Periparus ater","Cyanistes caeruleus",
                                     "Parus major","Sitta europea","Certhia familiaris","Certhia brachydactyla",
                                     "Fringilla coelebs","Carduelis spinus","Carduelis flammea","Pyrrhula pyrrhula",
                                     "Coccothraustes coccothraustes","Oriolus oriolus","Garrulus glandarius")

# Denmark
# https://www.dof.dk/fakta-om-fugle/punkttaellingsprogrammet

species_sub <- species_den_farm <- c("Falco tinnunculus","Perdix perdix","Vanellus vanellus","Gallinago gallinago",
                                     "Alauda arvensis","Hirundo rustica","Anthus pratensis","Motacilla flava",
                                     "Motacilla alba","Saxicola rubetra","Oenanthe oenanthe","Turdus pilaris",
                                     "Sylvia curruca","Curruca communis","Lanius collurio","Corvus frugilegus",
                                     "Corvus corone","Passer montanus","Carduelis carduelis","Linaria cannabina",
                                     "Emberiza citrinella","Emberiza calandra")

species_sub <- species_den_forest <- c("Accipiter nisus","Columba oenas","Dryocopus martius","Dendrocopos major",
                                       "Erithacus rubecula","Phoenicurus phoenicurus","Turdus viscivorus","Sylvia borin",
                                       "Phylloscopus sibilatrix","Phylloscopus collybita","Regulus regulus","Ficedula hypoleuca",
                                       "Poecile palustris","Lophophanes cristatus","Periparus ater","Sitta europaea",
                                       "Certhia familiaris","Garrulus glandarius","Corvus corax","Fringilla coelebs",
                                       "Carduelis spinus","Pyrrhula pyrrhula","Coccothraustes coccothraustes")




# Estonia

# Finland

# France

# Germany

# Greece

# Hungary

# Italy

# Latvia

# Lithuania

# Luxembourg

# Netherlands

# Norway

# Poland

# Portugal

# Republic of Ireland

# Romania

# Slovakia

# Slovenia

# Spain

# Sweden

# Switzerland

# United Kingdom     
# https://www.rspb.org.uk/our-work/conservation/conservation-and-sustainability/farming/near-you/farmland-bird-indicator/
# and https://www.gov.uk/government/statistics/wild-bird-populations-in-the-uk

species_sub <- species_uk_farm <- c("Passer montanus", "Emberiza calandra", "Streptopelia turtur", "Perdix perdix",
                                    "Motacilla flava", "Sturnus vulgaris", "Linaria cannabina", "Vanellus vanellus", 
                                    "Emberiza citrinella", "Alauda arvensis", "Falco tinnunculus", "Emberiza schoeniclus",
                                    "Curruca communis", "Chloris chloris", "Corvus frugilegus", "Columba oenas",  
                                    "Carduelis carduelis", "Columba palumbus", "Coloeus monedula")

species_sub <- species_uk_forest <- c("Turdus merula","Cyanistes caeruleus","Pyrrhula pyrrhula","Fringilla coelebs",
                                    "Prunella modularis","Parus major","Sylvia curruca","Aegithalos caudatus",
                                    "Erithacus rubecula","Turdus philomelos","Strix aluco","Troglodytes troglodytes",
                                    "Sylvia atricapilla","Phylloscopus collybita","Periparus ater","Sylvia borin",
                                    "Regulus regulus","Dendrocopos major","Picus viridis","Garrulus glandarius",
                                    "Dendrocopos minor","Poecile palustris","Luscinia megarhynchos","Sitta europaea",
                                    "Carduelis cabaret","Phoenicurus phoenicurus","Accipiter nisus","Muscicapa striata",
                                    "Anthus trivialis","Certhia familiaris","Poecile montana","Phylloscopus trochilus",
                                    "Ficedula hypoleuca","Phylloscopus sibilatrix","Loxia curvirostra","Carduelis spinus",
                                    "Tetrao urogallus")



# Run the DFA cluster analysis