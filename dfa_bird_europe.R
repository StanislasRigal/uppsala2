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
# https://www.ebcc.info/wp-content/uploads/2020/06/bcn-30-1.pdf#page=33

species_sub <- species_cyp_farm <- c("Falco tinnunculus","Alectoris chukar","Francolinus francolinus","Coturnix coturnix",
                                     "Columba palumbus","Streptopelia turtur","Clamator glandarius","Athene noctua",
                                     "Coracias garrulus","Galerida cristata","Hirundo rustica","Oenanthe cypriaca",
                                     "Cisticola juncidis","Iduna pallida","Sylvia conspicillata","Sylvia melanocephala",
                                     "Parus major","Pica pica","Corvus corone","Passer hispaniolensis",
                                     "Chloris chloris","Carduelis carduelis","Linaria cannabina","Emberiza melanocephala",
                                     "Emberiza calandra")

species_sub <- species_cyp_forest <- c("Columba palumbus","Streptopelia turtur","Troglodytes troglodytes","Oenanthe cypriaca",
                                       "Cettia cetti","Hippolais pallida","Sylvia melanothorax","Periparus ater",
                                       "Parus major","Certhia brachydactyla","Lanius nubicus","Garrulus glandarius",
                                       "Fringilla coelebs","Serinus serinus","Chloris chloris","Carduelis carduelis",
                                       "Emberiza caesia")



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
# https://www.sciencedirect.com/science/article/pii/S0167880901002122#TBL2

species_sub <- species_est_farm <- c("Ciconia ciconia","Crex crex","Tringa totanus","Delichon urbicum",
                                     "Anthus pratensis","Locustella fluviatilis","Sylvia nisoria","Muscicapa striata",
                                     "Carpodacus erythrinus","Acrocephalus dumetorum","Corvus corone", "Turdus pilaris",
                                     "Carduelis carduelis","Locustella naevia","Carduelis chloris","Falco subbuteo",
                                     "Vanellus vanellus","Carduelis cannabina","Pica pica","Circus aeruginosus",
                                     "Acrocephalus palustris","Emberiza hortulana","Perdix perdix","Phasianus colchicus",
                                     "Coturnix coturnix","Lanius collurio","Emberiza schoeniclus","Acrocephalus scirpaceus",
                                     "Acrocephalus schoenobaenus","Alauda arvensis","Sturnus vulgaris","Columba oenas",
                                     "Passer montanus","Oenanthe oenanthe","Saxicola rubetra","Motacilla alba",
                                     "Sylvia communis","Columba palumbus","Motacilla flava","Emberiza citrinella")


species_sub <- species_est_forest <- c("Luscinia luscinia","Dryocopus martius","Dendrocopos minor","Phylloscopus sibilatrix",
                                       "Ficedula parva","Parus montanus","Parus cristatus","Certhia familiaris",
                                       "Oriolus oriolus","Nucifraga caryocatactes","Pyrrhula pyrrhula","Coccothraustes coccothraustes",
                                       "Turdus merula","Parus caeruleus","Fringilla coelebs","Cuculus canorus",
                                       "Prunella modularis","Sylvia borin","Regulus regulus","Dendrocopus major",
                                       "Parus major","Hippolais icterina","Sylvia curruca","Parus palustris",
                                       "Sitta europaea","Ficedula hypoleuca","Corvus corax","Turdus iliacus",
                                       "Erithacus rubecula",	"Carduelis spinus",	"Turdus philomelos","Anthus trivialis",
                                       "Phylloscopus trochilus","Lullula arborea")



# Finland
#https://www.biodiversity.fi/ext/en/data-pages/fa8_backgroundinfo.html
#https://www.biodiversity.fi/en/habitats/forests/fo10-forest-birds

species_sub <- species_fin_farm <- c("Crex crex","Vanellus vanellus","Numenius arquata","Alauda arvensis",
                                     "Hirundo rustica","Delichon urbica","Anthus pratensis","Saxicola rubertra",
                                     "Turdus pilaris","Sylvia communis","Corvus monedula","Sturnus vulgaris",
                                     "Passer montanus","Emberiza hortulana")

species_sub <- species_fin_forest <- c("Cuculus canorus","Jynx torquilla","Dendrocopos major","Anthus trivialis",
                                       "Troglodytes troglodytes","Prunella modularis","Luscinia luscinia","Turdus merula",
                                       "Turdus philomelos","Turdus iliacus","Sylvia curruca","Sylvia borin",
                                       "Sylvia atricapilla","Phylloscopus trochilus","Muscicapa striata","Ficedula hypoleuca",
                                       "Corvus corax","Fringilla coelebs","Fringilla montifringilla","Carduelis flammea")


# France
# https://www.vigienature.fr/fr/produire-indicateurs-partir-indices-especes-habitat-2819
# et naturefrance.fr/sites/default/files/2020-08/160513_note_methodologique_indice_stoc.pdf

species_sub <- species_fra_farm <- c("Vanellus vanellus","Buteo buteo","Falco tinnunculus","Alectoris rufa",
                                     "Perdix perdix","Coturnix coturnix","Upupa epops","Alauda arvensis",
                                     "Lullula arborea","Galerida cristata","Anthus pratensis","Anthus campestris",
                                     "Motacilla flava","Curruca communis","Saxicola torquatus","Saxicola rubetra",
                                     "Oenanthe oenanthe","Lanius collurio","Corvus frugilegus","Linaria cannabina",
                                     "Emberiza citrinella","Emberiza cirlus","Emberiza calandra","Emberiza hortulana")

species_sub <- species_fra_forest <- c("Dendrocopos major","Dendrocoptes medius","Picus canus","Dryocopus martius",
                                       "Curruca melanocephala","Phylloscopus bonelli","Phylloscopus sibilatrix","Phylloscopus collybita",
                                       "Phylloscopus trochilus","Regulus regulus","Regulus ignicapilla","Sitta europaea",
                                       "Certhia brachydactyla","Certhia familiaris","Troglodytes troglodytes","Turdus philomelos",
                                       "Turdus viscivorus","Erithacus rubecula","Lophophanes cristatus","Periparus ater",
                                       "Poecile palustris","Poecile montanus","Coccothraustes coccothraustes","Pyrrhula pyrrhula")


# Germany
# https://link.springer.com/article/10.1007/s10336-020-01830-4#additional-information

species_sub <- species_ger_farm <- c("Buteo buteo","Falco tinnunculus","Phasianus colchicus","Coturnix coturnix",
                                     "Sturnus vulgaris","Curruca communis","Emberiza calandra","Alauda arvensis",
                                     "Passer montanus","Streptopelia turtur","Turdus pilaris","Perdix perdix",
                                     "Anthus pratensis","Vanellus vanellus","Milvus milvus","Lanius collurio",
                                     "Saxicola rubetra","Motacilla flava","Emberiza citrinella")

species_sub <- species_ger_forest <- c("Dryocopus martius","Sylvia atricapilla","Cyanistes caeruleus","Periparus ater",
                                       "Phylloscopus collybita","Corvus corax","Lophophanes cristatus","Pyrrhula pyrrhula",
                                       "Fringilla coelebs","Oriolus oriolus","Garrulus glandarius","Certhia familiaris",
                                       "Ficedula hypoleuca","Erithacus rubecula","Regulus ignicapilla","Sylvia borin",
                                       "Regulus regulus","Dendrocopos major","Picus canus","Coccothraustes coccothraustes",
                                       "Dryobates minor","Aegithalos caudatus","Poecile palustris","Dendrocoptes medius",
                                       "Turdus viscivorus","Certhia brachydactyla","Turdus philomelos","Columba oenas",
                                       "Anthus trivialis","Poecile montanus","Phylloscopus trochilus","Troglodytes troglodytes",
                                       "Sitta europaea","Phylloscopus sibilatrix")


# Greece
# same as EU, https://www.ornithologiki.gr/en/our-work/conservation-scientific-research/bird-monitoring/1298-programma-parakoloythisis-ton-koinon-eidon-poulion-tis-elladas


# Hungary
# https://www.mme.hu/tovabbra-csokken-hazai-mezogazdasagi-teruletek-madarvilaga

species_sub <- species_hun_farm <- c("Falco tinnunculus","Perdix perdix","Coturnix coturnix","Vanellus vanellus",
                                     "Merops apiaster","Galerida cristata","Alauda arvensis","Anthus campestris",
                                     "Motacilla flava","Locustella naevia","Curruca nisoria","Curruca communis",
                                     "Lanius collurio","Lanius minor","Sturnus vulgaris","Emberiza calandra")

species_sub <- species_hun_forest <- c("Columba oenas","Dryocopus martius","Dendrocopos major","Dendrocoptes medius",
                                       "Dryobates minor","Lullula arborea","Troglodytes troglodytes","Prunella modularis",
                                       "Erithacus rubecula","Turdus philomelos","Turdus viscivorus","Phylloscopus sibilatrix",
                                       "Phylloscopus collybita","Ficedula albicollis","Poecile palustris","Periparus ater",
                                       "Cyanistes caeruleus","Sitta europaea","Certhia brachydactyla","Garrulus glandarius",
                                       "Fringilla coelebs","Coccothraustes coccothraustes")



# Italy
# https://www.reterurale.it/flex/cm/pages/ServeBLOB.php/L/IT/IDPagina/15032

species_sub <- species_ita_farm <- c("Alauda arvensis","Lanius collurio","Motacilla alba","Melanocorypha calandra",
                                     "Calandrella brachydactyla","Anthus campestris","Galerida cristata","Carduelis carduelis",
                                     "Corvus cornix","Motacilla flava","Pica pica","Falco tinnunculus",
                                     "Emberiza hortulana","Passer italiae","Passer montanus","Passer hispaniolensis",
                                     "Oriolus oriolus","Hirundo rustica","Saxicola torquatus","Sturnus vulgaris",
                                     "Sturnus unicolor","Emberiza calandra","Jynx torquilla","Streptopelia turtur",
                                     "Upupa epops","Luscinia megarhynchos","Carduelis chloris","Serinus serinus")

species_sub <- species_ita_forest <- c("Poecile montanus","Poecile palustris","Lophophanes cristatus","Periparus ater",
                                       "Cyanistes caeruleus","Pyrrhula pyrrhula","Aegithalos caudatus","Regulus ignicapilla",
                                       "Fringilla coelebs","Garrulus glandarius","Phylloscopus bonelli","Phylloscopus collybita",
                                       "Nucifraga caryocatactes","Erithacus rubecula","Sitta europaea","Dryocopus martius",
                                       "Dendrocopos major","Certhia familiaris","Certhia brachydactyla","Regulus regulus",
                                       "Troglodytes troglodytes","Turdus viscivorus","Turdus philomelos")


# Latvia
# https://www.researchgate.net/profile/Oskars-Keiss/publication/268363185_Experiences_with_a_Baseline_Indicator_Farmland_Bird_Index_in_Latvia/links/55097a2f0cf26ff55f859259/Experiences-with-a-Baseline-Indicator-Farmland-Bird-Index-in-Latvia.pdf
# https://www.lu.lv/fileadmin/user_upload/LU.LV/Apaksvietnes/Konferences/Comm_ForestBirds_Trends_LV_AuninsA_20181205.pdf

species_sub <- species_lat_farm <- c("Ciconia ciconia","Crex crex","Vanellus vanellus","Alauda arvensis",
                                     "Anthus pratensis","Locustella naevia","Acrocephalus palustris","Saxicola rubetra",
                                     "Carduelis carduelis","Linaria cannabina","Carpodacus erythrinus","Emberiza citrinella")


species_sub <- species_lat_forest <- c("Accipiter gentilis","Accipiter nisus","Bonasa bonasia","Picus canus",
                                       "Dryocopus martius","Dryobates minor","Dendrocopos leucotos","Picoides tridactylus",
                                       "Turdus viscivorus","Phylloscopus sibilatrix","Regulus regulus","Ficedula parva",
                                       "Ficedula hypoleucos","Aegithalos caudatus","Poecile palustris","Poecile montanus",
                                       "Lophophanes cristatus","Periparus ater","Certhia familiaris","Nucifraga caryocatactes",
                                       "Loxia curvirostra","Pyrrhula pyrrhula","Coccothraustes coccothraustes")


# Lithuania

# Luxembourg

# Netherlands
# https://stats.sovon.nl/pub/publicatie/18092

Patrijs
Zomertortel
Kievit
Grutto
Spreeuw
Scholekster 
Wulp 
Tureluur
Kuifeend
Zomertaling
Torenvalk 
Kwartel
Ringmus
Slobeend
Geelgors
Wintertaling
Watersnip 
Grote
Lijster 
Grauwe
Gors
Veldleeuwerik 
Gele
Kwikstaart 
Spotvogel 
Boerenzwaluw
Roodborsttapuit
Krakeend 
Grasmus
Putter



Spreeuw 
Wilde Eend 
Fitis 
Merel
Houtduif
Zwarte
Kraai
Groenling
Heggenmus
Oeverzwaluw
Tuinfluiter 
Witte Kwikstaart
Zanglijster 
Meerkoet
Pimpelmees 
Rietgors 
Bosrietzanger 
Koolmees 
Tjiftjaf 
Vink
Kleine Karekiet
Winterkoning 
Grote Bonte Specht 
Grasmus 
Zwartkop 
Kerkuil 
Putter 
Grauwe Gans



# Norway

# Poland

# Portugal

# Republic of Ireland

# Romania

# Slovakia
# https://www.enviroportal.sk/indicator/detail?id=4041&print=yes
Alauda arvensis, Carduelis cannabina, Carduelis carduelis, Emberiza calandra, Emberiza citrinella, Falco tinnunculus, Hirundo rustica, Chloris chloris, Lanius collurio, Locustella naevia, Motacilla flava, Passer montanus, Saxicola rubetra, Saxicola torquata, Serinus serinus, Streptopelia turtur, Sturnus vulgaris, Sylvia communis, Sylvia nisoria, Vanellus vanellus.

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