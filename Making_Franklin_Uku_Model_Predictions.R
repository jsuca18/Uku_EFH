setwd("/home/jsuca/Uku_EFH_Recreation/Bathymetry")
library(raster)
library(terra)
library(gbm)
library(dismo)
library(matrixStats)
library(ggplot2)
library(viridis)

#generating the layers for the prediction maps
#for some reason predictions were not behaving when using raster stack, so 
#they were switched to dataframe. 
Uku_Layers<-stack("Uku_Prediction_Layers.gri")
names(Uku_Layers)<-c("Depth", "Aspect","Slope","Rugosity","Max_WvHgt","Mean_WvHgt","Extra")
Tight<-extent(-160.3, -154.7, 18.8,22.3)

Uku_Prediction_Layers_Trimmed<-crop(Uku_Layers, Tight)
Uku_Prediction_Layers_Trimmed<-Uku_Prediction_Layers_Trimmed[[1:6]]
Uku_Prediction_Layers_Agg<-aggregate(Uku_Prediction_Layers_Trimmed, fact=10)
Uku_Prediction_DF<-as.data.frame(rasterToPoints(Uku_Prediction_Layers_Agg))
Uku_Prediction_DF<-Uku_Prediction_DF[!is.na(Uku_Prediction_DF$Depth),]
Uku_Prediction_DF<-Uku_Prediction_DF[Uku_Prediction_DF$Depth>-300 & Uku_Prediction_DF$Depth<1,,]


saveRDS(Uku_Prediction_DF, "Uku_DF_for_Predictions.rds")


setwd("/home/jsuca/Uku_EFH_Recreation")
Uku_Model_Shallow<-readRDS("Shallow_Erik_Uku_Model_Ensemble_Sandless.rds")

Model_Prediction_Estimates<-matrix(, nrow=nrow(Uku_Prediction_DF), ncol=100)

#shallow model, no backscatter/sand
for (k in 1:100){
  Model_Prediction_Estimates[,k]<-predict.gbm(Uku_Model_Shallow[[1]][[k]],Uku_Prediction_DF, 
                                              n.trees=Uku_Model_Shallow[[1]][[k]]$gbm.call$best.trees, type="response")
  print(paste("completed", k, "of 100"))}


saveRDS(Model_Prediction_Estimates, "Erik_Uku_Model_Shallow_Sandless_Full.rds")

#Deep model with backscatter
Uku_Model_Deep<-readRDS("Deep_Erik_Uku_Model_Ensemble_Sand.rds")
setwd("/home/jsuca/Uku_EFH_Recreation/MHI_backscatterSynthesis")
BackScatter<-raster("BackScatter.tif")


Uku_Prediction_DF$BackScat<-raster::extract(BackScatter, cbind(Uku_Prediction_DF$x, Uku_Prediction_DF$y))

Model_Prediction_Estimates_Deep<-matrix(, nrow=nrow(Uku_Prediction_DF), ncol=100)

for (k in 1:100){
  Model_Prediction_Estimates_Deep[,k]<-predict.gbm(Uku_Model_Deep[[1]][[k]],Uku_Prediction_DF, 
                                              n.trees=Uku_Model_Deep[[1]][[k]]$gbm.call$best.trees, type="response")
  print(paste("completed", k, "of 100"))}


saveRDS(Model_Prediction_Estimates_Deep, "Erik_Uku_Model_Deep_Sand_Full.rds")



###############################processing estimates###########################


Uku_Prediction_DF<-readRDS( "Uku_DF_for_Predictions.rds")

Model_Prediction_Estimates<-readRDS("Erik_Uku_Model_Shallow_Sandless_Full.rds")

#getting the mean, CV, and range for both shallow and deep models
Mean_Uku<-rowMeans(Model_Prediction_Estimates, na.rm=T)
Min_Uku<-rowMins(Model_Prediction_Estimates, na.rm=T)
Max_Uku<-rowMaxs(Model_Prediction_Estimates, na.rm=T)
Range_Uku<-Max_Uku-Min_Uku
CV_Uku<-rowSds(Model_Prediction_Estimates)/rowMeans(Model_Prediction_Estimates, na.rm=T)

Mean_Uku_Deep<-rowMeans(Model_Prediction_Estimates_Deep, na.rm=T)
Min_Uku_Deep<-rowMins(Model_Prediction_Estimates_Deep, na.rm=T)
Max_Uku_Deep<-rowMaxs(Model_Prediction_Estimates_Deep, na.rm=T)
Range_Uku_Deep<-Max_Uku_Deep-Min_Uku_Deep
CV_Uku_Deep<-rowSds(Model_Prediction_Estimates_Deep)/rowMeans(Model_Prediction_Estimates_Deep, na.rm=T)



df_predictions<-cbind(Uku_Prediction_DF, Mean_Uku, Range_Uku, CV_Uku, Mean_Uku_Deep, Range_Uku_Deep, CV_Uku_Deep)
df_predictions$Lat<-df_predictions$y
df_predictions$Lon<-df_predictions$x

#parsing which model to use for which depths 
for (i in 1:nrow(df_predictions)){
df_predictions$Pr1_Uku[i]<-ifelse(df_predictions$Depth[i]>= -30,df_predictions$Mean_Uku[i], df_predictions$Mean_Uku_Deep[i])
df_predictions$Pr1_Uku_Range[i]<-ifelse(df_predictions$Depth[i]>= -30,df_predictions$Range_Uku[i], df_predictions$Range_Uku_Deep[i])
df_predictions$Pr1_Uku_CV[i]<-ifelse(df_predictions$Depth[i]>= -30,df_predictions$CV_Uku[i], df_predictions$CV_Uku_Deep[i])


}


png("MHI_Example_Full_Mean.png", height=6, width=8, res=300, units="in")
ggplot()+coord_fixed(ratio = 1)+geom_raster(data= df_predictions,aes(x=Lon, y=Lat, fill=Pr1_Uku))+scale_fill_viridis_c( guide = guide_colourbar(title="Pr (1)"))+
  theme_bw()+geom_sf()+coord_sf(xlim=c(-160.3,-154.7), ylim=c(18.8, 22.3))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=14),legend.direction = "vertical", legend.box = "vertical")+ggtitle("Uku Pr(1)")
dev.off()

png("MHI_Example_Full_Range.png", height=6, width=8, res=300, units="in")
ggplot()+coord_fixed(ratio = 1)+geom_raster(data= df_predictions,aes(x=Lon, y=Lat, fill=Pr1_Uku_Range))+scale_fill_viridis_c( guide = guide_colourbar(title="Pr (1) Range"))+
  theme_bw()+geom_sf()+coord_sf(xlim=c(-160.3,-154.7), ylim=c(18.8, 22.3))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=14),legend.direction = "vertical", legend.box = "vertical")+ggtitle("Uku Pr (1) Range")
dev.off()

png("MHI_Example_Full_CV.png", height=6, width=8, res=300, units="in")
ggplot()+coord_fixed(ratio = 1)+geom_raster(data= df_predictions,aes(x=Lon, y=Lat, fill=Pr1_Uku_CV))+scale_fill_viridis_c( guide = guide_colourbar(title="Pr (1) CV"))+
  theme_bw()+geom_sf()+coord_sf(xlim=c(-160.3,-154.7), ylim=c(18.8, 22.3))+
  theme(legend.title=element_text(size=16),legend.text=element_text(size=14),legend.direction = "vertical", legend.box = "vertical")+ggtitle("Uku Pr(1) CV")
dev.off()


#save just the individual model estimates that are used
Uku_Model_Predictions<-df_predictions[,16:20]

saveRDS(Uku_Model_Predictions, "JJS_version_of_Franklin_Uku_Model.rds")

