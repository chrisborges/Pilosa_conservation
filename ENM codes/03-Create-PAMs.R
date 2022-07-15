#####################################################################
## The threshold code was adapted from the function coded by       ##
## Cecina Babich Morrow and available at:                          ##                         #
## <https://babichmorrowc.github.io/post/2019-04-12-sdm-threshold/>##
####################################################################

library(raster)

enm_threshold = function(files_suit, files_obs, type="mtp", dir){ 
  for(sp in 1:length(files_obs)){
    # read suitability and occurrences files
    suit_df = read.csv(files_suit[sp], h=T)
    obs = read.csv(files_obs[sp], h=T)
    obs = obs[,c(1:2)] # keep only the coordinates
    
    # transform df to raster
    suit=suit_df
    gridded(suit) = ~x+y
    suit = stack(suit)
    
    #species name
    species=gsub(".*[/]([^.]+)[_].*","\\1",files_obs)
    print(species[sp])
    
    suit_pre=suit$presente
    occPredVals <- raster::extract(suit_pre, obs)
    
    if(type=="mtp"){
      thresh = min(na.omit(occPredVals))
    } else if(type=="p10"){
      if(length(occPredVals) < 10){
        p10 <- floor(length(occPredVals) * 0.9)
      } else {
        p10 <- ceiling(length(occPredVals) * 0.9)
      }
      thresh <- rev(sort(occPredVals))[p10]
    }
    
    # data.frame
    suit_df_pam = ifelse(suit_df[,c(3:7)] >= thresh, yes=1, no=0)
    suit_df_pam = cbind(suit_df[,c(1:2)], suit_df_pam) # binding coords
    
    # save all as dataframe
    write.csv(suit_df_pam, paste0(dir,"ENM/",species[sp],"_PAM.csv"), row.names=F)
  }
  return(suit_df_pam)
  print("All PAMs were created")
}


maps_dir=("/Pilosa/Data/")
# read all suitability files
files_suit = list.files(paste0(maps_dir, "ENM"), pattern="*Suitability.csv", full.names=T)
# read all occurrence files
files_obs = list.files(paste0(maps_dir, "SP-occur-env"), pattern="*_var.csv", full.names=T)

# run binary map function
enm_threshold(files_suit, files_obs, type="p10", dir=maps_dir)

########
## create richness maps ##
files_pam = list.files(paste0(maps_dir, "ENM"), pattern="*PAM.csv", full.names=T)

bra_tor_p = read.csv(files_pam[1], h=T)
bra_tri_p = read.csv(files_pam[2], h=T)
bra_var_p = read.csv(files_pam[3], h=T)
cho_did_p = read.csv(files_pam[4], h=T)
cho_hof_p = read.csv(files_pam[5], h=T)
cyc_did_p = read.csv(files_pam[6], h=T)
cyc_dor_p = read.csv(files_pam[7], h=T)
cyc_ida_p = read.csv(files_pam[8], h=T)
cyc_tho_p = read.csv(files_pam[9], h=T)
cyc_xin_p = read.csv(files_pam[10], h=T)
myr_tri_p = read.csv(files_pam[11], h=T)
tam.mex.p = read.csv(files_pam[12], h=T)
tam_tet_p = read.csv(files_pam[13], h=T)

# coords
coords = bra_tor_p[,c(1:2)]

########
# 0k
present=cbind(bra_tor_p[,3], bra_tri_p[,3], bra_var_p[,3],
              cho_did_p[,3], cho_hof_p[,3],
              cyc_did_p[,3], cyc_dor_p[,3], cyc_ida_p[,3],
              cyc_tho_p[,3], cyc_xin_p[,3],
              myr_tri_p[,3],
              tam.mex.p[,3],tam_tet_p[,3])

head(present)
present = as.data.frame(present)
names(present)=c("Bra_tor","Bra_tri","Bra.var",
                  "Cho_did","Cho_hof",
                  "Cyc_did", "Cyc_dor", "Cyc_ida",
                  "Cyc_tho", "Cyc_xin",
                  "Myr_tri",
                  "Tam.mex", "Tam.tet")

present["Richness"] = rowSums(present)
present = cbind(coords, present)
head(present)

########
# LGM 21k
lgm=cbind(bra_tor_p[,4], bra_tri_p[,4], bra_var_p[,4],
          cho_did_p[,4], cho_hof_p[,4],
          cyc_did_p[,4], cyc_dor_p[,4], cyc_ida_p[,4],
          cyc_tho_p[,4], cyc_xin_p[,4],
          myr_tri_p[,4],
          tam.mex.p[,4],tam_tet_p[,4])

head(lgm)
lgm = as.data.frame(lgm)
names(lgm)=c("Bra_tor","Bra_tri","Bra.var",
             "Cho_did","Cho_hof",
             "Cyc_did", "Cyc_dor", "Cyc_ida",
             "Cyc_tho", "Cyc_xin",
             "Myr_tri",
             "Tam.mex", "Tam.tet")

lgm["Richness"] = rowSums(lgm)
lgm = cbind(coords, lgm)
head(lgm)

########
# rcp26
rcp26=cbind(bra_tor_p[,5], bra_tri_p[,5], bra_var_p[,5],
            cho_did_p[,5], cho_hof_p[,5],
            cyc_did_p[,5], cyc_dor_p[,5], cyc_ida_p[,5],
            cyc_tho_p[,5], cyc_xin_p[,5],
            myr_tri_p[,5],
            tam.mex.p[,5],tam_tet_p[,5])

head(rcp26)
rcp26 = as.data.frame(rcp26)
names(rcp26)=c("Bra_tor","Bra_tri","Bra.var",
               "Cho_did","Cho_hof",
               "Cyc_did", "Cyc_dor", "Cyc_ida",
               "Cyc_tho", "Cyc_xin",
               "Myr_tri",
               "Tam.mex", "Tam.tet")

rcp26["Richness"] = rowSums(rcp26)
rcp26 = cbind(coords, rcp26)
head(rcp26)

########
# rcp45
rcp45=cbind(bra_tor_p[,6], bra_tri_p[,6], bra_var_p[,6],
            cho_did_p[,6], cho_hof_p[,6],
            cyc_did_p[,6], cyc_dor_p[,6], cyc_ida_p[,6],
            cyc_tho_p[,6], cyc_xin_p[,6],
            myr_tri_p[,6],
            tam.mex.p[,6],tam_tet_p[,6])

head(rcp45)
rcp45 = as.data.frame(rcp45)
names(rcp45)=c("Bra_tor","Bra_tri","Bra.var",
               "Cho_did","Cho_hof",
               "Cyc_did", "Cyc_dor", "Cyc_ida",
               "Cyc_tho", "Cyc_xin",
               "Myr_tri",
               "Tam.mex", "Tam.tet")

rcp45["Richness"] = rowSums(rcp45)
rcp45 = cbind(coords, rcp45)
head(rcp45)

########
# rcp85
rcp85=cbind(bra_tor_p[,7], bra_tri_p[,7], bra_var_p[,7],
            cho_did_p[,7], cho_hof_p[,7],
            cyc_did_p[,7], cyc_dor_p[,7], cyc_ida_p[,7],
            cyc_tho_p[,7], cyc_xin_p[,7],
            myr_tri_p[,7],
            tam.mex.p[,7],tam_tet_p[,7])

head(rcp85)
rcp85 = as.data.frame(rcp85)
names(rcp85)=c("Bra_tor","Bra_tri","Bra.var",
               "Cho_did","Cho_hof",
               "Cyc_did", "Cyc_dor", "Cyc_ida",
               "Cyc_tho", "Cyc_xin",
               "Myr_tri",
               "Tam.mex", "Tam.tet")

rcp85["Richness"] = rowSums(rcp85)
rcp85 = cbind(coords, rcp85)
head(rcp85)

########
# all richness 2gether
allrich=cbind(present[,16], lgm[,16], rcp26[,16], rcp45[,16], rcp85[,16])
allrich=as.data.frame(allrich)
names(allrich)=c("0k","LGM","RCP_2.6","RCP_4.5","RCP_8.5")
allrich=cbind(coords, allrich)
head(allrich)

####### 
# saving it
write.csv(present, paste0(maps_dir,"ENM/","Richness_0k.csv"), row.names=F)
write.csv(lgm, paste0(maps_dir,"ENM/","Richness_21k.csv"), row.names=F)
write.csv(rcp26, paste0(maps_dir,"ENM/","Richness_rcp26.csv"), row.names=F)
write.csv(rcp45, paste0(maps_dir,"ENM/","Richness_rcp45.csv"), row.names=F)
write.csv(rcp85, paste0(maps_dir,"ENM/","Richness_rcp85.csv"), row.names=F)

write.csv(allrich, paste0(maps_dir,"ENM/","Richness_all_time.csv"), row.names=F)

