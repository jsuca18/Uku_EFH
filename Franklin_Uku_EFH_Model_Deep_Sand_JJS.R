setwd("/home/jsuca/Uku_EFH_Recreation")
library(tidyr)
library(dplyr)
library(raster)
library(ncdf4)
library(terra)
library(matrixStats)
library(fmsb)

Full_Fish_Data<-read.csv("Uku_Bmouss_DropCam.csv")
Full_Fish_Data<-Full_Fish_Data[Full_Fish_Data$Species=="Aprion virescens" | Full_Fish_Data$Species=="No Fish in Drop",]



setwd("/home/jsuca/Uku_EFH_Recreation/Bathymetry")
MHI_Bathy<-raster("mhi_mbsyn_bathytopo_50m_v21.nc")
Aspect<-terra::terrain(MHI_Bathy, opt="aspect", unit="degrees", neighbors=8)
Slope1<-terra::terrain(MHI_Bathy, opt="slope", unit="degrees", neighbors=8)
Slope<-terra::terrain(Slope1, opt="slope", unit="degrees", neighbors=8)
Rugosity<-terra::terrain(MHI_Bathy, opt="TRI", unit="degrees", neighbors=8)

setwd("/home/jsuca/Uku_EFH_Recreation/MHI_backscatterSynthesis")

BackScatter<-raster("BackScatter.tif")


Uku_PA_Data<-Full_Fish_Data[, c(3,7,8,9,14) ]

colnames(Uku_PA_Data)<-c("Date","Latitude","Longitude","Depth","PA")
Uku_PA_Data$Depth<-as.numeric(Uku_PA_Data$Depth)*-1
Uku_PA_Data$Latitude<-as.numeric(as.character(Uku_PA_Data$Latitude))
Uku_PA_Data$Longitude<-as.numeric(as.character(Uku_PA_Data$Longitude))

Uku_PA_Data$PA[Uku_PA_Data$PA>0]<-1
Locations<-c(Uku_PA_Data$Longitude, Uku_PA_Data$Latitude)
Uku_PA_Data$Rugosity<-raster::extract(Rugosity, cbind(Uku_PA_Data$Longitude, Uku_PA_Data$Latitude))
Uku_PA_Data$Aspect<-raster::extract(Aspect, cbind(Uku_PA_Data$Longitude, Uku_PA_Data$Latitude))
Uku_PA_Data$Slope<-raster::extract(Slope, cbind(Uku_PA_Data$Longitude, Uku_PA_Data$Latitude))
Uku_PA_Data$BackScat<-raster::extract(BackScatter, cbind(Uku_PA_Data$Longitude, Uku_PA_Data$Latitude))



setwd("/home/jsuca/Longline_Projects")
source("BRT_Eval_Function_JJS.R")
setwd("/home/jsuca/Uku_EFH_Recreation")
Uku_PA_Data$Depth[Uku_PA_Data$Depth=="NaN"]<-NA

Uku_Models_Est<-list()

Uku_Models_Est<-fit.brt.n_eval_Balanced_Fixed(Uku_PA_Data, gbm.x = c(4,6:9), gbm.y = 5,lr=0.001, tc=4,family="bernoulli",nt=3300, bag.fraction = 0.5,100)

saveRDS(Uku_Models_Est,"Deep_Erik_Uku_Model_Ensemble_Sand.rds")

Model_Skill_PA<-matrix(,1,4)
Model_Evals_PA<-unlist(Uku_Models_Est[[2]])
Model_PA_Eval<-matrix(,100,2)

    for (i in 1:100){
      Model_PA_Eval[i,1]<-Model_Evals_PA[[i]]@auc
      Model_PA_Eval[i,2]<-max(Model_Evals_PA[[i]]@TPR+Model_Evals_PA[[i]]@TNR-1)
    }
    Model_Skill_PA[1,1]<-mean( Model_PA_Eval[,1])
    Model_Skill_PA[1,2]<-min( Model_PA_Eval[,1])
    Model_Skill_PA[1,3]<-mean( Model_PA_Eval[,2])
    Model_Skill_PA[1,4]<-min( Model_PA_Eval[,2])
rownames(Model_Skill_PA)<-"Uku_Deep"
colnames(Model_Skill_PA)<-c("Mean_AUC","Min_AUC","Mean_TSS","Min_TSS")
write.csv(Model_Skill_PA,"PA_Metrics_Uku_Deep_Erik_Sand.csv")

# #####Full model plot#########
`%nin%` = Negate(`%in%`)
    Fish_Models<-Uku_Models_Est[[1]]
    Species_Name<-"Uku_Deep"
    var_tested<-Fish_Models[[1]]$var.names
    #var_tested<-colnames(PA_Data_Reef_Fishes_MHI[Predictors])
    Fish_Models_Good<-Fish_Models#[Q]
    percent_contrib<-NULL#list()
    iters=length(Fish_Models_Good)
    part_plot<-list()
    part_plot<-list()
    percent_contrib<-NULL#list()
    Continuous_Preds<-which(var_tested %nin% c("Year"))
    for(q in 1:iters){                                #this was 50
      mod<-Fish_Models_Good[q][[1]]
      ###
      part_plot1<-data.frame(row.names=1:100)
      for(x in c(Continuous_Preds)){ ###
        
        pp<-plot(mod ,var_tested[x],return.grid=T) ###
        part_plot1<-cbind(part_plot1, pp) ###
      }
      
      #   ###
      part_plot[[q]]<-part_plot1 ###
      
      sum1<-summary(Fish_Models_Good[q][[1]]  , plot=F )
      sum2<-sum1[order(sum1[,1], levels = var_tested),]
      percent_contrib<-cbind(percent_contrib, sum2[,2])
      rownames(percent_contrib)<-sum1[order(sum1[,1], levels = var_tested),1]
    }
    All_percent_contribution<-cbind(rownames(percent_contrib), paste(round(rowMeans(percent_contrib),2), round(rowSds(percent_contrib),2), sep=" ± "))
    Combined_All_percent_contribution<-All_percent_contribution
    #
    #
    Mean_PA_Contributions<-as.data.frame(t(rowMeans(percent_contrib)))
    write.csv( Mean_PA_Contributions,paste0("Var_Contributions_",Species_Name,"_PA_Sand_Deep.csv"))
    
    PA_Predictors_Plot<- rbind(rep(max(Mean_PA_Contributions),length(var_tested)) , rep(0,length(var_tested)) , Mean_PA_Contributions)
    PA_Predictors_Plot[]<-sapply(PA_Predictors_Plot, as.numeric)
    par(mfrow=c(1,1))
    
    png(paste0("Radar_Chart_",Species_Name,"_PA_Sand_Deep.png"), height=6, width=6, units="in",res=300)
    radarchart(PA_Predictors_Plot,  pfcol=rgb(0.0,0.3,0.5,0.5), pcol=rgb(0.0,0.3,0.5,0.5), title=paste0(Species_Name,"_PA"))
    dev.off()
    #
    All_percent_contribution<-cbind(rownames(percent_contrib), paste(round(rowMeans(percent_contrib),2), round(rowSds(percent_contrib),2), sep=" ± "))
    #
    png(paste0("Partial_plots_",Species_Name,"_PA_Sand_Deep.png"), height=18,width=14, res=300, units="in")
    par(mfrow=c(4,2))
    mn_part_plot<-list()
    for(y in c(Continuous_Preds)){
      id<-which(colnames(part_plot[[1]])==var_tested[y])
      all1<-NULL
      all2<-NULL
      for(z in 1:iters){											 #this was 50
        all1<-rbind(all1, cbind(c(part_plot[[z]][,id])))
        all2<-rbind(all2, cbind(c(part_plot[[z]][,id+1])))
      }
      all3<-cbind(all1, all2)
      all1<-all3[order(all3[,1]),]
      #
      plot(all1, xlab=var_tested[y], col="white", ylab=paste("f(",var_tested[y], ")", sep=""),cex.axis=1.2, cex.lab=1.2)
      plx<-predict(loess(all1[,2] ~ all1[,1], span = 0.3), se=T)
      mn_part_plot[[y]]<- cbind(all1[,1], plx$fit)
      lines(all1[,1],plx$fit)
      lines(all1[,1],plx$fit - qt(0.975,plx$df)*plx$se, lty=2)#0.975
      lines(all1[,1],plx$fit + qt(0.975,plx$df)*plx$se, lty=2)
      rug(na.omit(unlist(Uku_PA_Data[var_tested[y]])))
      legend("bottomright", paste(All_percent_contribution[which(All_percent_contribution[,1]==var_tested[y]),2],"%", sep=" "), bty="n", cex=1.4)
    }
    dev.off()
    
    
   