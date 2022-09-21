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

# Belgium

# Bulgaria           

# Cyprus

# Czech Republic

# Denmark            

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



species_sub <- species_uk_farm <- c("Passer montanus", "Emberiza calandra", "Streptopelia turtur", "Perdix perdix",
                                    "Motacilla flava", "Sturnus vulgaris", "Linaria cannabina", "Vanellus vanellus", 
                                    "Emberiza citrinella", "Alauda arvensis", "Falco tinnunculus", "Emberiza schoeniclus",
                                    "Curruca communis", "Chloris chloris", "Corvus frugilegus", "Columba oenas",  
                                    "Carduelis carduelis", "Columba palumbus", "Coloeus monedula")




# Run the DFA cluster analysis