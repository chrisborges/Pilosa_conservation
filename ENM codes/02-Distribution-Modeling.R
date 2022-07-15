rm(list=ls())

library(sp)
library(raster)
library(dismo)
library(kernlab) # SVM
library(randomForest)
library(rJava)
library(vegan)
library(rgdal)
suppressMessages(library(dismo))

ENM_func = function(sp_back, sp_ocor, cross_validation){
  
  sp_back = lapply(files_back, read.csv)
  sp_ocor = lapply(files_ocor, read.csv)
  
  for(sp in 1:length(sp_back)){
    species=gsub(".*[-]([^.]+)[.].*","\\1",files_back)
    print(species[sp])
    ocor = sp_ocor[[sp]]
    back = sp_back[[sp]]
    #back = back[sample(nrow(back), nrow(ocor)),]
    back = back[,-8]
    names(back)[names(back) == "x"] = "longitude"
    names(back)[names(back) == "y"] = "latitude"
    back = back[,c(6,7,1:5)]
    ocor = ocor[,-3]
    
    #shapefile to clip
    shape = readOGR(dsn=paste0(main_dir,"Environmental/"), layer="Neotropic-shape")
    
    #####
    # loop for different AOGCMs
    Output0k = NULL
    Output21k = NULL
    Outputrcp26 = NULL
    Outputrcp45 = NULL
    Outputrcp85 = NULL
    
    #To save AUC, TSS, thresholds
    Bioclim.auc = NULL
    Maxent.auc = NULL
    SVM.auc = NULL
    GLM.auc = NULL
    Gower.auc = NULL
    RndFor.auc = NULL
    
    AOGCMs = c("CCSM", "GISS", "MIROC", "MRI")
    
    for(j in AOGCMs){
      print(j)
      # reading climatic data 
      clim_dir=paste0(main_dir,"Environmental/")
      AOGCM0k = paste(clim_dir,"be_biovar_" ,j, "_Modern(1950-1999)_selecao.csv", sep="")
      AOGCM21k = paste(clim_dir,"be_biovar_" ,j, "_LGM_Modern(1950-1999)_selecao.csv", sep="")
      AOGCMrcp26 = paste(clim_dir,"be_biovar_" ,j, "_rcp26_Modern(1950-1999)_selecao.csv", sep="")
      AOGCMrcp45 = paste(clim_dir,"be_biovar_" ,j, "_rcp45_Modern(1950-1999)_selecao.csv", sep="")
      AOGCMrcp85 = paste(clim_dir,"be_biovar_" ,j, "_rcp85_Modern(1950-1999)_selecao.csv", sep="")
      
      clima0k = read.table(AOGCM0k, h=T, sep=",")
      clima21k = read.table(AOGCM21k, h=T, sep=",")
      climarcp26 = read.table(AOGCMrcp26, h=T, sep=",")
      climarcp45 = read.table(AOGCMrcp45, h=T, sep=",")
      climarcp85 = read.table(AOGCMrcp85, h=T, sep=",")
      
      gridded(clima0k) = ~long+lat
      gridded(clima21k) = ~long+lat
      gridded(climarcp26) = ~long+lat
      gridded(climarcp45) = ~long+lat
      gridded(climarcp85) = ~long+lat
      
      clima0k.r = stack(clima0k)
      clima21k.r = stack(clima21k)
      climarcp26.r = stack(climarcp26)
      climarcp45.r = stack(climarcp45)
      climarcp85.r = stack(climarcp85)
      
      # clip to mainland area
      clima0k.r = crop(clima0k.r, extent(shape))
      clima0k.r = mask(clima0k.r, shape)
      clima21k.r = crop(clima21k.r, extent(shape))
      clima21k.r = mask(clima21k.r, shape)
      climarcp26.r = crop(climarcp26.r, extent(shape))
      climarcp26.r = mask(climarcp26.r, shape)
      climarcp45.r = crop(climarcp45.r, extent(shape))
      climarcp45.r = mask(climarcp45.r, shape)
      climarcp85.r = crop(climarcp85.r, extent(shape))
      climarcp85.r = mask(climarcp85.r, shape)
      
      # save partial results from each cross validation loop
      Bioclim.Pout0k = NULL
      Maxent.Pout0k = NULL
      SVM.Pout0k = NULL
      GLM.Pout0k = NULL
      Gower.Pout0k = NULL
      RndFor.Pout0k = NULL
      
      Bioclim.Pout21k = NULL
      Maxent.Pout21k = NULL
      SVM.Pout21k = NULL
      GLM.Pout21k = NULL
      Gower.Pout21k = NULL
      RndFor.Pout21k = NULL
      
      Bioclim.Poutrcp26 = NULL
      Maxent.Poutrcp26 = NULL
      SVM.Poutrcp26 = NULL
      GLM.Poutrcp26 = NULL
      Gower.Poutrcp26 = NULL
      RndFor.Poutrcp26 = NULL
      
      Bioclim.Poutrcp45 = NULL
      Maxent.Poutrcp45 = NULL
      SVM.Poutrcp45 = NULL
      GLM.Poutrcp45 = NULL
      Gower.Poutrcp45 = NULL
      RndFor.Poutrcp45 = NULL
      
      Bioclim.Poutrcp85 = NULL
      Maxent.Poutrcp85 = NULL
      SVM.Poutrcp85 = NULL
      GLM.Poutrcp85 = NULL
      Gower.Poutrcp85 = NULL
      RndFor.Poutrcp85 = NULL
      
      Bioclim.values = NULL
      Maxent.values = NULL
      SVM.values = NULL
      GLM.values = NULL
      Gower.values = NULL
      RndFor.values = NULL
        
      ###
      # loop para cross-validation
      for(i in 1:cross_validation){
        
        # dados de treino e teste
        sample.ocor = sample(1:nrow(ocor), round(0.75*nrow(ocor)))
        sample.back = sample(1:nrow(back), round(0.75*nrow(back)))
        
        treino = prepareData(x= clima0k.r, p= ocor[sample.ocor,1:2], b= back[sample.back,1:2])
        teste = prepareData(x=clima0k.r, p=ocor[-sample.ocor,1:2], b=back[-sample.back,1:2])
 
        ######
        ## Bioclim
        
        #ajustando o modelo
        Bioclim.model = bioclim(treino[treino[,"pb"]==1, -1])
        
        #fazendo predicoes
        Bioclim0k = predict(clima0k.r, Bioclim.model, type="response")
        Bioclim21k = predict(clima21k.r, Bioclim.model, type="response")
        Bioclimrcp26 = predict(climarcp26.r, Bioclim.model, type="response")
        Bioclimrcp45 = predict(climarcp45.r, Bioclim.model, type="response")
        Bioclimrcp85 = predict(climarcp85.r, Bioclim.model, type="response")

        #avaliando o modelo
        Bioclim.eval = evaluate(p= teste[teste[,"pb"]==1, -1], a= teste[teste[,"pb"]==0, -1], model= Bioclim.model)
        Bioclim.auc = Bioclim.eval@auc
        Bioclim.thr = threshold(Bioclim.eval)
        Bioclim_values = as.data.frame(cbind(Bioclim.auc, Bioclim.thr))
        
        ######
        ## Maxent
       
        #ajustando o modelo
        Sys.setenv(NOAWT=T)
        Maxent.model = maxent(treino[,-c(1)] , p=treino[,"pb"])
        
        #fazendo predicoes
        Maxent0k = predict(clima0k.r, Maxent.model, type="response")
        Maxent21k = predict(clima21k.r, Maxent.model, type="response")
        Maxentrcp26 = predict(climarcp26.r, Maxent.model, type="response")
        Maxentrcp45 = predict(climarcp45.r, Maxent.model, type="response")
        Maxentrcp85 = predict(climarcp85.r, Maxent.model, type="response")
        
        #avaliando o modelo
        Maxent.eval = evaluate(p=teste[teste[,"pb"]==1, -1], a = teste[teste[,"pb"]==0, -1], model=Maxent.model)
        Maxent.auc = Maxent.eval@auc
        Maxent.thr = threshold(Maxent.eval)
        Maxent_values = as.data.frame(cbind(Maxent.auc, Maxent.thr))
        
        
        ######
        ## SVM
        
        #ajustando o modelo
        SVM.model = ksvm(pb ~ bio.10+bio.16+bio.17+bio.2+bio.4, data= treino)

        #fazendo predicoes
        SVM0k = predict(clima0k.r, SVM.model, type="response")
        SVM21k = predict(clima21k.r, SVM.model, type="response")
        SVMrcp26 = predict(climarcp26.r, SVM.model, type="response")
        SVMrcp45 = predict(climarcp45.r, SVM.model, type="response")
        SVMrcp85 = predict(climarcp85.r, SVM.model, type="response")
        
        #avaliando o modelo
        SVM.eval = evaluate(p=teste[teste[,"pb"]==1, -1], a=teste[teste[,"pb"]==0, -1], model= SVM.model)
        SVM.auc = SVM.eval@auc
        SVM.thr = threshold(SVM.eval)
        SVM_values = as.data.frame(cbind(SVM.auc, SVM.thr))

        
        ######
        ## GLM
        
        #ajustando o modelo
        GLM.model = glm(pb ~ bio.10+bio.16+bio.17+bio.2+bio.4, family= binomial(link="logit"), data= treino)
        
        #fazendo predicoes
        GLM0k = predict(clima0k.r, GLM.model, type="response")
        GLM21k = predict(clima21k.r, GLM.model, type="response")
        GLMrcp26 = predict(climarcp26.r, GLM.model, type="response")
        GLMrcp45 = predict(climarcp45.r, GLM.model, type="response")
        GLMrcp85 = predict(climarcp85.r, GLM.model, type="response")
        
        #plot(GLM0k)

        #avaliando o modelo
        GLM.eval = evaluate(p= teste[teste[,"pb"]==1, -1], a= teste[teste[,"pb"]==0, -1], model= GLM.model)
        GLM.auc = GLM.eval@auc
        GLM.thr = threshold(GLM.eval)
        GLM_values = as.data.frame(cbind(GLM.auc, GLM.thr))

        
        ######
        ## Gower
        
        #ajustando o modelo
        Gower.model = domain(treino[treino[,"pb"]==1, -1])
        
        #fazendo predicoes
        Gower0k = predict(clima0k.r, Gower.model, type="response")
        Gower21k = predict(clima21k.r, Gower.model, type="response")
        Gowerrcp26 = predict(climarcp26.r, Gower.model, type="response")
        Gowerrcp45 = predict(climarcp45.r, Gower.model, type="response")
        Gowerrcp85 = predict(climarcp85.r, Gower.model, type="response")
        
        #avaliando o modelo
        Gower.eval = evaluate(p= teste[teste[,"pb"]==1, -1], a= teste[teste[,"pb"]==0, -1], model= Gower.model)
        Gower.auc = Gower.eval@auc
        Gower.thr = threshold(Gower.eval)
        Gower_values = as.data.frame(cbind(Gower.auc, Gower.thr))
        

        #########
        # Random Forest
        RndFor.model = randomForest(pb ~ bio.10+bio.16+bio.17+bio.2+bio.4, data= treino)
        
        RndFor0k = predict(clima0k.r, RndFor.model, type="response")
        RndFor21k = predict(clima21k.r, RndFor.model, type="response")
        RndForrcp26 = predict(climarcp26.r, RndFor.model, type="response")
        RndForrcp45 = predict(climarcp45.r, RndFor.model, type="response")
        RndForrcp85 = predict(climarcp85.r, RndFor.model, type="response")
        
        #avaliando o modelo
        RndFor.eval = evaluate(p=teste[teste[,"pb"]==1, -1], a=teste[teste[,"pb"]==0, -1], RndFor.model)
        RndFor.auc = RndFor.eval@auc
        RndFor.thr = threshold(RndFor.eval)
        RndFor_values = as.data.frame(cbind(RndFor.auc, RndFor.thr))
        
        #####
        # saving partial outputs for each "i" cross-validation loop
        
        Bioclim.Pout0k = cbind(Bioclim.Pout0k, values(Bioclim0k))
        Maxent.Pout0k = cbind(Maxent.Pout0k, values(Maxent0k))
        SVM.Pout0k = cbind(SVM.Pout0k, values(SVM0k))
        GLM.Pout0k = cbind(GLM.Pout0k, values(GLM0k))
        Gower.Pout0k = cbind(Gower.Pout0k, values(Gower0k))
        RndFor.Pout0k = cbind(RndFor.Pout0k, values(RndFor0k))
        
        Bioclim.Pout21k = cbind(Bioclim.Pout21k, values(Bioclim21k))
        Maxent.Pout21k = cbind(Maxent.Pout21k, values(Maxent21k))
        SVM.Pout21k = cbind(SVM.Pout21k, values(SVM21k))
        GLM.Pout21k = cbind(GLM.Pout21k, values(GLM21k))
        Gower.Pout21k = cbind(Gower.Pout21k, values(Gower21k))
        RndFor.Pout21k = cbind(RndFor.Pout21k, values(RndFor21k))
        
        Bioclim.Poutrcp26 = cbind(Bioclim.Poutrcp26, values(Bioclimrcp26))
        Maxent.Poutrcp26 = cbind(Maxent.Poutrcp26, values(Maxentrcp26))
        SVM.Poutrcp26 = cbind(SVM.Poutrcp26, values(SVMrcp26))
        GLM.Poutrcp26 = cbind(GLM.Poutrcp26, values(GLMrcp26))
        Gower.Poutrcp26 = cbind(Gower.Poutrcp26, values(Gowerrcp26))
        RndFor.Poutrcp26 = cbind(RndFor.Poutrcp26, values(RndForrcp26))
        
        Bioclim.Poutrcp45 = cbind(Bioclim.Poutrcp45, values(Bioclimrcp45))
        Maxent.Poutrcp45 = cbind(Maxent.Poutrcp45, values(Maxentrcp45))
        SVM.Poutrcp45 = cbind(SVM.Poutrcp45, values(SVMrcp45))
        GLM.Poutrcp45 = cbind(GLM.Poutrcp45, values(GLMrcp45))
        Gower.Poutrcp45 = cbind(Gower.Poutrcp45, values(Gowerrcp45))
        RndFor.Poutrcp45 = cbind(RndFor.Poutrcp45, values(RndForrcp45))
        
        Bioclim.Poutrcp85 = cbind(Bioclim.Poutrcp85, values(Bioclimrcp85))
        Maxent.Poutrcp85 = cbind(Maxent.Poutrcp85, values(Maxentrcp85))
        SVM.Poutrcp85 = cbind(SVM.Poutrcp85, values(SVMrcp85))
        GLM.Poutrcp85 = cbind(GLM.Poutrcp85, values(GLMrcp85))
        Gower.Poutrcp85 = cbind(Gower.Poutrcp85, values(Gowerrcp85))
        RndFor.Poutrcp85 = cbind(RndFor.Poutrcp85, values(RndForrcp85))

        #salvando auc, tss, thresholds
        Bioclim.values = rbind(Bioclim.values, (Bioclim_values))
        Maxent.values = rbind(Maxent.values, Maxent_values)
        SVM.values = rbind(SVM.values, (SVM_values))
        GLM.values = rbind(GLM.values, (GLM_values))
        Gower.values = rbind(Gower.values, (Gower_values))
        RndFor.values = rbind(RndFor.values, (RndFor_values))
   
      } #fecha for 'i' - cross-validation
    
      save_dir=paste0(main_dir, "ENM-output/")
      write.csv(Bioclim.Pout0k, paste0(save_dir,species[sp],"_",j,"_","Bioclim_0k.csv"), row.names=F)
      write.csv(Bioclim.Pout21k, paste0(save_dir,species[sp],"_",j,"_","Bioclim_21k.csv"), row.names=F)
      write.csv(Bioclim.Poutrcp26, paste0(save_dir,species[sp],"_",j,"_","Bioclim_rcp26.csv"), row.names=F)
      write.csv(Bioclim.Poutrcp45, paste0(save_dir,species[sp],"_",j,"_","Bioclim_rcp45.csv"), row.names=F)
      write.csv(Bioclim.Poutrcp85, paste0(save_dir,species[sp],"_",j,"_","Bioclim_rcp85.csv"), row.names=F)
      
      write.csv(Maxent.Pout0k, paste0(save_dir,species[sp],"_",j,"_","Maxent_0k.csv"), row.names=F)
      write.csv(Maxent.Pout21k, paste0(save_dir,species[sp],"_",j,"_","Maxent_21k.csv"), row.names=F)
      write.csv(Maxent.Poutrcp26, paste0(save_dir,species[sp],"_",j,"_","Maxent_rcp26.csv"), row.names=F)
      write.csv(Maxent.Poutrcp45, paste0(save_dir,species[sp],"_",j,"_","Maxent_rcp45.csv"), row.names=F)
      write.csv(Maxent.Poutrcp85, paste0(save_dir,species[sp],"_",j,"_","Maxent_rcp85.csv"), row.names=F)
      
      write.csv(SVM.Pout0k, paste0(save_dir,species[sp],"_",j,"_","SVM_0k.csv"), row.names=F)
      write.csv(SVM.Pout21k, paste0(save_dir,species[sp],"_",j,"_","SVM_21k.csv"), row.names=F)
      write.csv(SVM.Poutrcp26, paste0(save_dir,species[sp],"_",j,"_","SVM_rcp26.csv"), row.names=F)
      write.csv(SVM.Poutrcp45, paste0(save_dir,species[sp],"_",j,"_","SVM_rcp45.csv"), row.names=F)
      write.csv(SVM.Poutrcp85, paste0(save_dir,species[sp],"_",j,"_","SVM_rcp85.csv"), row.names=F)
      
      write.csv(GLM.Pout0k, paste0(save_dir,species[sp],"_",j,"_","GLM_0k.csv"), row.names=F)
      write.csv(GLM.Pout21k, paste0(save_dir,species[sp],"_",j,"_","GLM_21k.csv"), row.names=F)
      write.csv(GLM.Poutrcp26, paste0(save_dir,species[sp],"_",j,"_","GLM_rcp26.csv"), row.names=F)
      write.csv(GLM.Poutrcp45, paste0(save_dir,species[sp],"_",j,"_","GLM_rcp45.csv"), row.names=F)
      write.csv(GLM.Poutrcp85, paste0(save_dir,species[sp],"_",j,"_","GLM_rcp85.csv"), row.names=F)
      
      write.csv(Gower.Pout0k, paste0(save_dir,species[sp],"_",j,"_","Gower_0k.csv"), row.names=F)
      write.csv(Gower.Pout21k, paste0(save_dir,species[sp],"_",j,"_","Gower_21k.csv"), row.names=F)
      write.csv(Gower.Poutrcp26, paste0(save_dir,species[sp],"_",j,"_","Gower_rcp26.csv"), row.names=F)
      write.csv(Gower.Poutrcp45, paste0(save_dir,species[sp],"_",j,"_","Gower_rcp45.csv"), row.names=F)
      write.csv(Gower.Poutrcp85, paste0(save_dir,species[sp],"_",j,"_","Gower_rcp85.csv"), row.names=F)
      
      write.csv(RndFor.Pout0k, paste0(save_dir,species[sp],"_",j,"_","RndFor_0k.csv"), row.names=F)
      write.csv(RndFor.Pout21k, paste0(save_dir,species[sp],"_",j,"_","RndFor_21k.csv"), row.names=F)
      write.csv(RndFor.Poutrcp26, paste0(save_dir,species[sp],"_",j,"_","RndFor_rcp26.csv"), row.names=F)
      write.csv(RndFor.Poutrcp45, paste0(save_dir,species[sp],"_",j,"_","RndFor_rcp45.csv"), row.names=F)
      write.csv(RndFor.Poutrcp85, paste0(save_dir,species[sp],"_",j,"_","RndFor_rcp85.csv"), row.names=F)
      
      # which are the best models?
      # remove poor performing model
      # AUC < 0.7 are bad models
      if(any(Bioclim.values$Bioclim.auc < 0.7)){
        bad_mods = which(Bioclim.values$Bioclim.auc < 0.7)
        Bioclim.values = Bioclim.values[-bad_mods,]
        Bioclim.Pout0k = Bioclim.Pout0k[,-bad_mods]
        Bioclim.Pout21k = Bioclim.Pout21k[,-bad_mods]
        Bioclim.Poutrcp26 = Bioclim.Poutrcp26[,-bad_mods]
        Bioclim.Poutrcp45 = Bioclim.Poutrcp45[,-bad_mods]
        Bioclim.Poutrcp85 = Bioclim.Poutrcp85[,-bad_mods]
      }
      
      if(any(Maxent.values$Maxent.auc < 0.7)){
        bad_mods = which(Maxent.values$Maxent.auc < 0.7)
        Maxent.values = Maxent.values[-bad_mods,]
        Maxent.Pout0k = Maxent.Pout0k[,-bad_mods]
        Maxent.Pout21k = Maxent.Pout21k[,-bad_mods]
        Maxent.Poutrcp26 = Maxent.Poutrcp26[,-bad_mods]
        Maxent.Poutrcp45 = Maxent.Poutrcp45[,-bad_mods]
        Maxent.Poutrcp85 = Maxent.Poutrcp85[,-bad_mods]
      }
      
      if(any(SVM.values$SVM.auc < 0.7)){
        bad_mods = which(SVM.values$SVM.auc < 0.7)
        SVM.values = SVM.values[-bad_mods,]
        SVM.Pout0k = SVM.Pout0k[,-bad_mods]
        SVM.Pout21k = SVM.Pout21k[,-bad_mods]
        SVM.Poutrcp26 = SVM.Poutrcp26[,-bad_mods]
        SVM.Poutrcp45 = SVM.Poutrcp45[,-bad_mods]
        SVM.Poutrcp85 = SVM.Poutrcp85[,-bad_mods]
      }
      
      if(any(GLM.values$GLM.auc < 0.7)){
        bad_mods = which(GLM.values$GLM.auc < 0.7)
        GLM.values = GLM.values[-bad_mods,]
        GLM.Pout0k = GLM.Pout0k[,-bad_mods]
        GLM.Pout21k = GLM.Pout21k[,-bad_mods]
        GLM.Poutrcp26 = GLM.Poutrcp26[,-bad_mods]
        GLM.Poutrcp45 = GLM.Poutrcp45[,-bad_mods]
        GLM.Poutrcp85 = GLM.Poutrcp85[,-bad_mods]
      }
      
      
      if(any(Gower.values$Gower.auc < 0.7)){
        bad_mods = which(Gower.values$Gower.auc < 0.7)
        Gower.values = Gower.values[-bad_mods,]
        Gower.Pout0k = Gower.Pout0k[,-bad_mods]
        Gower.Pout21k = Gower.Pout21k[,-bad_mods]
        Gower.Poutrcp26 = Gower.Poutrcp26[,-bad_mods]
        Gower.Poutrcp45 = Gower.Poutrcp45[,-bad_mods]
        Gower.Poutrcp85 = Gower.Poutrcp85[,-bad_mods]
      }
      
      
      if(any(RndFor.values$RndFor.auc < 0.7)){
        bad_mods = which(RndFor.values$RndFor.auc < 0.7)
        RndFor.values = RndFor.values[-bad_mods,]
        RndFor.Pout0k = RndFor.Pout0k[,-bad_mods]
        RndFor.Pout21k = RndFor.Pout21k[,-bad_mods]
        RndFor.Poutrcp26 = RndFor.Poutrcp26[,-bad_mods]
        RndFor.Poutrcp45 = RndFor.Poutrcp45[,-bad_mods]
        RndFor.Poutrcp85 = RndFor.Poutrcp85[,-bad_mods]
      }
      
      #calculando as medias dos modelos parciais cross-validation
      
      Bioclim.Pout0k.mean = apply(Bioclim.Pout0k, 1, mean)
      Maxent.Pout0k.mean = apply(Maxent.Pout0k, 1, mean)
      SVM.Pout0k.mean = apply(SVM.Pout0k, 1, mean)
      GLM.Pout0k.mean = apply(GLM.Pout0k, 1, mean)
      Gower.Pout0k.mean = apply(Gower.Pout0k, 1, mean)
      RndFor.Pout0k.mean = apply(RndFor.Pout0k, 1, mean)
      
      Bioclim.Pout21k.mean = apply(Bioclim.Pout21k, 1, mean)
      Maxent.Pout21k.mean = apply(Maxent.Pout21k, 1, mean)
      SVM.Pout21k.mean = apply(SVM.Pout21k, 1, mean)
      GLM.Pout21k.mean = apply(GLM.Pout21k, 1, mean)
      Gower.Pout21k.mean = apply(Gower.Pout21k, 1, mean)
      RndFor.Pout21k.mean = apply(RndFor.Pout21k, 1, mean)
      
      Bioclim.Poutrcp26.mean = apply(Bioclim.Poutrcp26, 1, mean)
      Maxent.Poutrcp26.mean = apply(Maxent.Poutrcp26, 1, mean)
      SVM.Poutrcp26.mean = apply(SVM.Poutrcp26, 1, mean)
      GLM.Poutrcp26.mean = apply(GLM.Poutrcp26, 1, mean)
      Gower.Poutrcp26.mean = apply(Gower.Poutrcp26, 1, mean)
      RndFor.Poutrcp26.mean = apply(RndFor.Poutrcp26, 1, mean)
      
      Bioclim.Poutrcp45.mean = apply(Bioclim.Poutrcp45, 1, mean)
      Maxent.Poutrcp45.mean = apply(Maxent.Poutrcp45, 1, mean)
      SVM.Poutrcp45.mean = apply(SVM.Poutrcp45, 1, mean)
      GLM.Poutrcp45.mean = apply(GLM.Poutrcp45, 1, mean)
      Gower.Poutrcp45.mean = apply(Gower.Poutrcp45, 1, mean)
      RndFor.Poutrcp45.mean = apply(RndFor.Poutrcp45, 1, mean)
      
      Bioclim.Poutrcp85.mean = apply(Bioclim.Poutrcp85, 1, mean)
      Maxent.Poutrcp85.mean = apply(Maxent.Poutrcp85, 1, mean)
      SVM.Poutrcp85.mean = apply(SVM.Poutrcp85, 1, mean)
      GLM.Poutrcp85.mean = apply(GLM.Poutrcp85, 1, mean)
      Gower.Poutrcp85.mean = apply(Gower.Poutrcp85, 1, mean)
      RndFor.Poutrcp85.mean = apply(RndFor.Poutrcp85, 1, mean)
      
      Output0k = cbind(Output0k, 
                       Bioclim= Bioclim.Pout0k.mean, 
                       Maxent= Maxent.Pout0k.mean,
                       SVM= SVM.Pout0k.mean, 
                       GLM= GLM.Pout0k.mean, 
                       Gower= Gower.Pout0k.mean, 
                       RndFor=RndFor.Pout0k.mean)  
      
      Output21k = cbind(Output21k, 
                        Bioclim= Bioclim.Pout21k.mean, 
                        Maxent= Maxent.Pout21k.mean,
                        SVM= SVM.Pout21k.mean, 
                        GLM= GLM.Pout21k.mean, 
                        Gower= Gower.Pout21k.mean, 
                        RndFor=RndFor.Pout21k.mean)
      
      Outputrcp26 = cbind(Outputrcp26, 
                          Bioclim= Bioclim.Poutrcp26.mean, 
                          Maxent= Maxent.Poutrcp26.mean,
                          SVM= SVM.Poutrcp26.mean, 
                          GLM= GLM.Poutrcp26.mean, 
                          Gower= Gower.Poutrcp26.mean, 
                          RndFor=RndFor.Poutrcp26.mean)
      
      Outputrcp45 = cbind(Outputrcp45, 
                          Bioclim= Bioclim.Poutrcp45.mean, 
                          Maxent= Maxent.Poutrcp45.mean,
                          SVM= SVM.Poutrcp45.mean, 
                          GLM= GLM.Poutrcp45.mean, 
                          Gower= Gower.Poutrcp45.mean, 
                          RndFor=RndFor.Poutrcp45.mean)
      
      Outputrcp85 = cbind(Outputrcp85, 
                          Bioclim= Bioclim.Poutrcp85.mean, 
                          Maxent= Maxent.Poutrcp85.mean, 
                          SVM= SVM.Poutrcp85.mean, 
                          GLM= GLM.Poutrcp85.mean, 
                          Gower= Gower.Poutrcp85.mean, 
                          RndFor=RndFor.Poutrcp85.mean)
      
    } # close 'j' - AOGCMs loop
    
    #inserindo coordenadas geograficas aos outputs
    coords = xyFromCell(clima0k.r, 1:ncell(clima0k.r))
    
    Output0k = cbind(coords, Output0k)
    Output21k = cbind(coords, Output21k)
    Outputrcp26 = cbind(coords, Outputrcp26)
    Outputrcp45 = cbind(coords, Outputrcp45)
    Outputrcp85 = cbind(coords, Outputrcp85)
    
    # exclude NAs from outputs
    Output0k = na.omit(Output0k)
    Output21k = na.omit(Output21k)
    Outputrcp26 = na.omit(Outputrcp26)
    Outputrcp45 = na.omit(Outputrcp45)
    Outputrcp85 = na.omit(Outputrcp85)
    
    # save each result
    save_dir=paste0(main_dir, "ENM-output/")
    write.csv(Output0k, paste0(save_dir,"AOGCM_av_",species[sp],"_0k.csv"), row.names=F)
    write.csv(Output21k, paste0(save_dir,"AOGCM_av_",species[sp],"_21k.csv"), row.names=F)
    write.csv(Outputrcp26, paste0(save_dir,"AOGCM_av_",species[sp],"_rcp26.csv"), row.names=F)
    write.csv(Outputrcp45, paste0(save_dir,"AOGCM_av_",species[sp],"_rcp46.csv"), row.names=F)
    write.csv(Outputrcp85, paste0(save_dir,"AOGCM_av_",species[sp],"_rcp85.csv"), row.names=F)
    
    write.csv(Bioclim.values, paste0(save_dir, species[sp], "_Bioclim.thr.csv"), row.names=F)
    write.csv(Maxent.values, paste0(save_dir, species[sp], "_Maxent.thr.csv"), row.names=F)
    write.csv(SVM.values, paste0(save_dir, species[sp], "_SVM.thr.csv"), row.names=F)
    write.csv(GLM.values, paste0(save_dir, species[sp], "_GLM.thr.csv"), row.names=F)
    write.csv(Gower.values, paste0(save_dir, species[sp], "_Gower.thr.csv"), row.names=F)
    write.csv(RndFor.values, paste0(save_dir, species[sp], "_RndFor.thr.csv"), row.names=F)

    # Ensemble
    
    mean_suit = data.frame(presente=rowMeans(Output0k[,-c(1:2)],na.rm=T), 
                               LGM=rowMeans(Output21k[,-c(1:2)],na.rm=T), 
                               Futuro26=rowMeans(Outputrcp26[,-c(1:2)],na.rm=T), 
                               Futuro45=rowMeans(Outputrcp45[,-c(1:2)],na.rm=T), 
                               Futuro85=rowMeans(Outputrcp85[,-c(1:2)],na.rm=T))
  
    mean_thr = mean(mean(Bioclim.values$prevalence), 
                    mean(Maxent.values$prevalence),
                    mean(SVM.values$prevalence), 
                    mean(GLM.values$prevalence),
                    mean(Gower.values$prevalence), 
                    mean(RndFor.values$prevalence))
    
    # create PAMs
    mean_suit=decostand(mean_suit, method= "range", margin=2)
    sp.pam=ifelse(mean_suit >= mean_thr, yes=1,no=0)
    coords = Output0k[,(1:2)]
    sp.pam=cbind(coords, sp.pam)
    sp.suit=cbind(coords, mean_suit)
    
    write.csv(sp.suit, paste0(main_dir,"ENM/",species[sp],"_Suitability.csv"), row.names=F)
    write.csv(sp.pam, paste0(main_dir,"ENM/",species[sp],"_PAM.csv"), row.names=F)
    
    } # close sp loop
  print("ENM is done.")
} # end function


# run ENM
main_dir="/Users/christiellyborges/Library/Mobile Documents/com~apple~CloudDocs/Research/Pilosa/Data/"
files_back = list.files(paste0(main_dir, "SP-occur-env"), pattern="Background_random", full.names=T)
files_ocor = list.files(paste0(main_dir, "SP-occur-env"), pattern="*_var.csv", full.names=T)


files_back=files_back[13]
files_ocor=files_ocor[13]
cross_validation=50

ENM_func(sp_back=files_back, sp_ocor=files_ocor, cross_val=50)
 

