library(maptools)
library(sp)
library(raster)
library(rgdal)
library(rgeos)

main_dir=("/Pilosa/Data/")

#lendo os arquivos raster salvos a partir do hd
bio10 <- raster(paste0(main_dir,"Environmental/","bio10_AS_0k.grd"))
bio16 <- raster(paste0(main_dir,"Environmental/","bio16_AS_0k.grd"))
bio17 <- raster(paste0(main_dir,"Environmental/","bio17_AS_0k.grd"))
bio2 <- raster(paste0(main_dir,"Environmental/","bio2_AS_0k.grd"))
bio4 <- raster(paste0(main_dir,"Environmental/","bio4_AS_0k.grd"))

# clip to mainland area
shape = readOGR(dsn=paste0(main_dir,"Environmental/"), layer="Neotropic-shape")
# crop and mask
bio10.crop = crop(bio10, extent(shape))
bio10.m = mask(bio10.crop, shape)

bio16.crop = crop(bio16, extent(shape))
bio16.m = mask(bio16.crop, shape)

bio17.crop = crop(bio17, extent(shape))
bio17.m = mask(bio17.crop, shape)

bio2.crop = crop(bio2, extent(shape))
bio2.m = mask(bio2.crop, shape)

bio4.crop = crop(bio4, extent(shape))
bio4.m = mask(bio4.crop, shape)

#
clima.AS <- stack(c(bio10.m,bio16.m, bio17.m, bio2.m, bio4.m))
names(clima.AS) <- c("bio10","bio16", "bio17", "bio2", "bio4")
plot(clima.AS)

### Occurrence point data

files = list.files(paste0(main_dir, "SP-occur-2021"), pattern="*.csv", full.names=T)
sp = lapply(files, read.csv)

for(i in 1:length(files)){
  sp_df = sp[[i]]
  save_as = sp_df[1,"species"]
  print(save_as)
  xy = sp_df[,c("longitude","latitude")]
  sp_shp = SpatialPointsDataFrame(coords=xy, data=xy)
  
  plot(clima.AS$bio10)
  points(sp_df[,"longitude"], sp_df[,"latitude"])
  
  # extraindo variáveis a partir dos pontos de ocorrencia
  sp.var <- extract(clima.AS, sp_shp, cellnumbers=T)
  sp.var <- cbind(xy, sp.var)
  sp.var    #note que existem NAs na matriz e células repetidas (duas ou mais ocorrencias em uma mesma celula)
  
  duplicated(sp.var[,"cells"])
  dup <- which(duplicated(sp.var[,"cells"]) == TRUE)
  
  sp.var <- sp.var[-dup,]  #note que ainda restaram NAs
  sp.var <- na.omit(sp.var)
  nrow(sp.var)
  points(sp.var[,"longitude"], sp.var[,"latitude"], col="blue")
  
  # create buffer around occur locations and define it as study region
  # this creates a 4-decimal-degree buffer around the occurrence data
  crs(sp_shp) = "+proj=longlat +datum=WGS84 +no_defs"
  crs(sp_shp) = "+proj=utm +datum=WGS84 +no_defs" # so buffer can be in degrees
  occ_buff = raster::buffer(sp_shp, 1, doEdge=T, dissolve=T)
  crs(sp_shp) = "+proj=longlat +datum=WGS84 +no_defs" # set back
  
  plot(clima.AS[[1]])
  plot(occ_buff, add=T, col="blue")
  plot(sp_shp, add=T, col="red", pch=20)
  
  # the 'study area' created by extracting the buffer area from the raster stack
  studyarea = mask(clima.AS, occ_buff, inverse=T)
  plot(studyarea$bio10)

  #### AMOSTRANDO BACKGROUND (aleatorio)
  bg = sampleRandom(x=studyarea,
                    size=1000,
                    na.rm=T,
                    sp=T)
  
  
  plot(clima.AS[[1]])
  plot(bg,add=T) # add the background points to the plotted raster
  plot(sp_shp,add=T,col="red") # add the occurrence data to the plotted raster
  
  #
  sp.bg = as.data.frame(bg)
  head(sp.bg)
  names(sp.bg)[names(sp.bg) == "x"] = "longitude"
  names(sp.bg)[names(sp.bg) == "y"] = "latitude"
  sp.bg = sp.bg[,c(6,7,1:5)]
  head(sp.bg)
  head(sp.var)

  # save dfs
  write.table(sp.bg, paste0(main_dir,"SP-occur-env/", "Background_random-", save_as,".csv"), row.names=F, sep=",")
  write.table(sp.var, paste0(main_dir,"SP-occur-env/", save_as,"_var.csv"), row.names=F, sep=",")
  
}


