# clear NeoXen database
library(dplyr)
library(plyr)

setwd("/Users/christiellyborges/Library/Mobile Documents/com~apple~CloudDocs/Research/Pilosa/Data/Raw-occurrence")

data=read.csv("NeoXen-Pilosa.csv", h=T)
head(data)


data_clean = data %>% distinct(SPECIES, LONG_X, LAT_Y, .keep_all= TRUE)
data_clean = data_clean[,c(5:8,13,16:18)]
head(data_clean)
names(data_clean) = c("species", "genus","family", "order", "locality",
                      "longitude", "latitude", "precision")

data_clean["source"] = "NeoXen"
data_clean$species = gsub(" ", "_", data_clean$species)
# Join data for each species #

# Bra_var
bra_var_df = read.table(paste0("SP-occur-2020/", "Bradypus_variegatus.txt"), h=T)
bra_var_df = bra_var_df[which(bra_var_df$precision==1),]
bra_var_neo = data_clean[which(data_clean$species=="Bradypus_variegatus"),]
bra_var = rbind.fill(bra_var_df, bra_var_neo)
bra_var_clean = bra_var %>% distinct(longitude, latitude, .keep_all= TRUE)
# nenhum dado repetido
# usar apenas o neoxen?

# Bra_tor
bra_tor_df = read.table(paste0("SP-occur-2020/", "Bradypus_torquatus.txt"), h=T)
bra_tor_neo = data_clean[which(data_clean$species=="Bradypus_torquatus"),]
bra_tor = rbind.fill(bra_tor_df, bra_tor_neo)
bra_tor_clean = bra_tor %>% distinct(longitude, latitude, .keep_all= TRUE)

# Bra_tri
bra_tri_df = read.table(paste0("SP-occur-2020/", "Bradypus_tridactylus.txt"), h=T)
bra_tri_df = bra_tri_df[which(bra_tri_df$precision==1),]
bra_tri_neo = data_clean[which(data_clean$species=="Bradypus_tridactylus"),]
bra_tri = rbind.fill(bra_tri_df, bra_tri_neo)
bra_tri_clean = bra_tri %>% distinct(longitude, latitude, .keep_all= TRUE)

# Cho_did
cho_did_df = read.table(paste0("SP-occur-2020/", "Choloepus_didactylus.txt"), h=T)
cho_did_df = cho_did_df[which(cho_did_df$precision==1),]
cho_did_neo = data_clean[which(data_clean$species=="Choloepus_didactylus"),]
cho_did = rbind.fill(cho_did_df, cho_did_neo)
cho_did_clean = cho_did %>% distinct(longitude, latitude, .keep_all= TRUE)

# Cho_hof
cho_hof_df = read.table(paste0("SP-occur-2020/", "Choloepus_hoffmanni.txt"), h=T)
cho_hof_df = cho_hof_df[which(cho_hof_df$precision==1),]
cho_hof_neo = data_clean[which(data_clean$species=="Choloepus_hoffmanni"),]
cho_hof = rbind.fill(cho_hof_df, cho_hof_neo)
cho_hof_clean = cho_hof %>% distinct(longitude, latitude, .keep_all= TRUE)

# Myr_tri
myr_tri_df = read.table(paste0("SP-occur-2020/", "Myrmecophaga_tridactyla.txt"), h=T)
myr_tri_df = myr_tri_df[which(myr_tri_df$precision==1),]
myr_tri_neo = data_clean[which(data_clean$species=="Myrmecophaga_tridactyla"),]
myr_tri = rbind.fill(myr_tri_df, myr_tri_neo)
myr_tri_clean = myr_tri %>% distinct(longitude, latitude, .keep_all= TRUE)

# Tam_mex
tam_mex_df = read.table(paste0("SP-occur-2020/", "Tamandua_mexicana.txt"), h=T)
tam_mex_df = tam_mex_df[which(tam_mex_df$precision==1),]
tam_mex_neo = data_clean[which(data_clean$species=="Tamandua_mexicana"),]
tam_mex = rbind.fill(tam_mex_df, tam_mex_neo)
tam_mex_clean = tam_mex %>% distinct(longitude, latitude, .keep_all= TRUE)

# Tam_tet
tam_tet_df = read.table(paste0("SP-occur-2020/", "Tamandua_tetradactyla.txt"), h=T)
tam_tet_df = tam_tet_df[which(tam_tet_df$precision==1),]
tam_tet_neo = data_clean[which(data_clean$species=="Tamandua_tetradactyla"),]
tam_tet = rbind.fill(tam_tet_df, tam_tet_neo)
tam_tet_clean = tam_tet %>% distinct(longitude, latitude, .keep_all= TRUE)

## save new raw occurrences
drive="/Users/christiellyborges/Library/Mobile Documents/com~apple~CloudDocs/Research/Pilosa/Data/SP-occur-2021/"
write.csv(bra_tor_clean, paste0(drive,"Bradypus_torquatus.csv"), row.names=F)


# usar apenas o banco de dados NeoXen
write.csv(bra_tor_neo, paste0(drive,"Bradypus_torquatus.csv"), row.names=F)
write.csv(bra_var_clean,  paste0(drive,"Bradypus_variegatus.csv"), row.names=F)
write.csv(bra_tri_neo,  paste0(drive,"Bradypus_tridactylus.csv"), row.names=F)
write.csv(myr_tri_neo,  paste0(drive,"Myrmecophaga_tridactyla.csv"), row.names=F)
write.csv(tam_mex_clean,  paste0(drive,"Tamandua_mexicana.csv"), row.names=F)
write.csv(tam_tet_neo,  paste0(drive,"Tamandua_tetradactyla.csv"), row.names=F)
write.csv(cho_did_clean,  paste0(drive,"Choloepus_didactylus.csv"), row.names=F)
write.csv(cho_hof_clean,  paste0(drive,"Choloepus_hoffmanni.csv"), row.names=F)
