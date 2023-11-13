# load the libries
library(sp)
library(spatstat)
library(ncf)
library(spdep)
library(data.table)
library(gam)
library(mgcv)
# set the working directoy
# setwd("~/Desktop/BIOMASS")
setwd("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation")

# delete the na values



#  Step 1 Moran's I dynamics along distance plotting
######################################################################################
#  Step 1 Moran's I dynamics along distance plotting
######################################################################################
# load the train table
rawDataTable = fread(paste("Data/SpatialCrossValidation/SDs/20221209_RemoteSensing_Biomass_Merged_sampled_dataset_1000000.csv",sep=""))[,-1] %>% na.omit()
trainVaribles = c('Aridity_Index','CHELSA_Annual_Mean_Temperature','CHELSA_Annual_Precipitation','CHELSA_Isothermality','CHELSA_Max_Temperature_of_Warmest_Month','CHELSA_Mean_Diurnal_Range','CHELSA_Mean_Temperature_of_Coldest_Quarter','CHELSA_Mean_Temperature_of_Driest_Quarter','CHELSA_Mean_Temperature_of_Warmest_Quarter','CHELSA_Mean_Temperature_of_Wettest_Quarter','CHELSA_Min_Temperature_of_Coldest_Month','CHELSA_Precipitation_Seasonality','CHELSA_Precipitation_of_Coldest_Quarter','CHELSA_Precipitation_of_Driest_Month','CHELSA_Precipitation_of_Driest_Quarter','CHELSA_Precipitation_of_Warmest_Quarter','CHELSA_Precipitation_of_Wettest_Month','CHELSA_Precipitation_of_Wettest_Quarter','CHELSA_Temperature_Annual_Range','CHELSA_Temperature_Seasonality','Depth_to_Water_Table','EarthEnvTopoMed_Eastness','EarthEnvTopoMed_Elevation','EarthEnvTopoMed_Northness','EarthEnvTopoMed_ProfileCurvature','EarthEnvTopoMed_Roughness','EarthEnvTopoMed_Slope','SG_Absolute_depth_to_bedrock','WorldClim2_SolarRadiation_AnnualMean','WorldClim2_WindSpeed_AnnualMean','EarthEnvCloudCover_MODCF_interannualSD','EarthEnvCloudCover_MODCF_intraannualSD','EarthEnvCloudCover_MODCF_meanannual','EarthEnvTopoMed_AspectCosine','EarthEnvTopoMed_AspectSine','LandCoverClass_Cultivated_and_Managed_Vegetation','Human_Disturbance','LandCoverClass_Urban_Builtup','SG_Clay_Content_0_100cm','SG_Coarse_fragments_0_100cm','SG_Sand_Content_0_100cm','SG_Silt_Content_0_100cm','SG_Soil_pH_H2O_0_100cm','cropland','grazing',"pasture","rangeland","PresentTreeCover") #

set.seed(10000)
trainTable = rawDataTable[sample(nrow(rawDataTable), 40000), ]

formulaString = as.formula(paste('densityHM','~', paste("ti(",trainVaribles,")", collapse="+",sep=""), sep=""))
# formulaString = as.formula(paste('RemoteBiomass','~', paste(trainVaribles, collapse="+",sep=""), sep=""))
gamMulti = mgcv::gam(data = trainTable, formula = formulaString) 
trainTable$gamResiduals = gamMulti$residuals

duplicateTable = trainTable
# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+init=epsg:4326")

# # calculate the Moran I changes along the spatial distance
spatialMoranDynamicsHM1 = spline.correlog(x=coordinates(duplicateTable)[,1], y=coordinates(duplicateTable)[,2],z=duplicateTable$gamResiduals, resamp=20, quiet=TRUE,latlon=T,xmax=1000,df=100,na.rm=T)
# plot(spatialMoranDynamics02,main="Remote Biomass Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1))

pdf(paste("Plots/Figure_SX_MoranI_along_distance_of_random_sampled_Remote_Sensing_Residual_with_1000_km_xmax_df_100_HM1.pdf",sep=""),width = 5, height=4)
plot(spatialMoranDynamicsHM1,main="HM1 Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1),xlab="Distance (km)",ylab="Moran's I")
abline(v=800,col="red", lty=2)
dev.off()



set.seed(10000)
trainTable = rawDataTable[sample(nrow(rawDataTable), 40000), ]

formulaString = as.formula(paste('densitySD','~', paste("ti(",trainVaribles,")", collapse="+",sep=""), sep=""))
# formulaString = as.formula(paste('RemoteBiomass','~', paste(trainVaribles, collapse="+",sep=""), sep=""))
gamMulti = mgcv::gam(data = trainTable, formula = formulaString) 
trainTable$gamResiduals = gamMulti$residuals

duplicateTable = trainTable
# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+init=epsg:4326")

# # calculate the Moran I changes along the spatial distance
spatialMoranDynamicsSD1 = spline.correlog(x=coordinates(duplicateTable)[,1], y=coordinates(duplicateTable)[,2],z=duplicateTable$gamResiduals, resamp=20, quiet=TRUE,latlon=T,xmax=1000,df=100,na.rm=T)
# plot(spatialMoranDynamics02,main="Remote Biomass Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1))

pdf(paste("Plots/Figure_SX_MoranI_along_distance_of_random_sampled_Remote_Sensing_Residual_with_1000_km_xmax_df_100_SD1.pdf",sep=""),width = 5, height=4)
plot(spatialMoranDynamicsSD1,main="SD1 Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1),xlab="Distance (km)",ylab="Moran's I")
abline(v=900,col="red", lty=2)
dev.off()


set.seed(10000)
trainTable = rawDataTable[sample(nrow(rawDataTable), 40000), ]

formulaString = as.formula(paste('densityWK','~', paste("ti(",trainVaribles,")", collapse="+",sep=""), sep=""))
# formulaString = as.formula(paste('RemoteBiomass','~', paste(trainVaribles, collapse="+",sep=""), sep=""))
gamMulti = mgcv::gam(data = trainTable, formula = formulaString) 
trainTable$gamResiduals = gamMulti$residuals

duplicateTable = trainTable
# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+init=epsg:4326")

# # calculate the Moran I changes along the spatial distance
spatialMoranDynamicsWK1 = spline.correlog(x=coordinates(duplicateTable)[,1], y=coordinates(duplicateTable)[,2],z=duplicateTable$gamResiduals, resamp=20, quiet=TRUE,latlon=T,xmax=1000,df=100,na.rm=T)
# plot(spatialMoranDynamics02,main="Remote Biomass Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1))

pdf(paste("Plots/Figure_SX_MoranI_along_distance_of_random_sampled_Remote_Sensing_Residual_with_1000_km_xmax_df_100_WK1.pdf",sep=""),width = 5, height=4)
plot(spatialMoranDynamicsWK1,main="WK1 Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1),xlab="Distance (km)",ylab="Moran's I")
abline(v=900,col="red", lty=2)
dev.off()


###############################################################################################################################
# For SD 2 models


rawDataTable = fread(paste("Data/SpatialCrossValidation/SDs/20221209_RemoteSensing_Biomass_Merged_sampled_dataset_1000000.csv",sep=""))[,-1] %>% na.omit() %>% dplyr::filter(WDPA==1|IntactForest==1) 
trainVaribles = c('Aridity_Index','CHELSA_Annual_Mean_Temperature','CHELSA_Annual_Precipitation','CHELSA_Isothermality','CHELSA_Max_Temperature_of_Warmest_Month','CHELSA_Mean_Diurnal_Range','CHELSA_Mean_Temperature_of_Coldest_Quarter','CHELSA_Mean_Temperature_of_Driest_Quarter','CHELSA_Mean_Temperature_of_Warmest_Quarter','CHELSA_Mean_Temperature_of_Wettest_Quarter','CHELSA_Min_Temperature_of_Coldest_Month','CHELSA_Precipitation_Seasonality','CHELSA_Precipitation_of_Coldest_Quarter','CHELSA_Precipitation_of_Driest_Month','CHELSA_Precipitation_of_Driest_Quarter','CHELSA_Precipitation_of_Warmest_Quarter','CHELSA_Precipitation_of_Wettest_Month','CHELSA_Precipitation_of_Wettest_Quarter','CHELSA_Temperature_Annual_Range','CHELSA_Temperature_Seasonality','Depth_to_Water_Table','EarthEnvTopoMed_Eastness','EarthEnvTopoMed_Elevation','EarthEnvTopoMed_Northness','EarthEnvTopoMed_ProfileCurvature','EarthEnvTopoMed_Roughness','EarthEnvTopoMed_Slope','SG_Absolute_depth_to_bedrock','WorldClim2_SolarRadiation_AnnualMean','WorldClim2_WindSpeed_AnnualMean','EarthEnvCloudCover_MODCF_interannualSD','EarthEnvCloudCover_MODCF_intraannualSD','EarthEnvCloudCover_MODCF_meanannual','EarthEnvTopoMed_AspectCosine','EarthEnvTopoMed_AspectSine','SG_Clay_Content_0_100cm','SG_Coarse_fragments_0_100cm','SG_Sand_Content_0_100cm','SG_Silt_Content_0_100cm','SG_Soil_pH_H2O_0_100cm',"PresentTreeCover")


set.seed(10000)
trainTable = rawDataTable[sample(nrow(rawDataTable), 40000), ]

formulaString = as.formula(paste('densityHM','~', paste("ti(",trainVaribles,")", collapse="+",sep=""), sep=""))
# formulaString = as.formula(paste('RemoteBiomass','~', paste(trainVaribles, collapse="+",sep=""), sep=""))
gamMulti = mgcv::gam(data = trainTable, formula = formulaString) 
trainTable$gamResiduals = gamMulti$residuals

duplicateTable = trainTable
# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+init=epsg:4326")

# # calculate the Moran I changes along the spatial distance
spatialMoranDynamicsHM2 = spline.correlog(x=coordinates(duplicateTable)[,1], y=coordinates(duplicateTable)[,2],z=duplicateTable$gamResiduals, resamp=20, quiet=TRUE,latlon=T,xmax=1000,df=100,na.rm=T)

pdf(paste("Plots/Figure_SX_MoranI_along_distance_of_random_sampled_Remote_Sensing_Residual_with_1000_km_xmax_df_100_HM2.pdf",sep=""),width = 5, height=4)
plot(spatialMoranDynamicsHM2,main="HM2 Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1),xlab="Distance (km)",ylab="Moran's I")
abline(v=550,col="red", lty=2)
dev.off()


set.seed(10000)
trainTable = rawDataTable[sample(nrow(rawDataTable), 40000), ]

formulaString = as.formula(paste('densitySD','~', paste("ti(",trainVaribles,")", collapse="+",sep=""), sep=""))
# formulaString = as.formula(paste('RemoteBiomass','~', paste(trainVaribles, collapse="+",sep=""), sep=""))
gamMulti = mgcv::gam(data = trainTable, formula = formulaString) 
trainTable$gamResiduals = gamMulti$residuals

duplicateTable = trainTable
# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+init=epsg:4326")

# # calculate the Moran I changes along the spatial distance
spatialMoranDynamicsSD2 = spline.correlog(x=coordinates(duplicateTable)[,1], y=coordinates(duplicateTable)[,2],z=duplicateTable$gamResiduals, resamp=20, quiet=TRUE,latlon=T,xmax=1000,df=100,na.rm=T)

pdf(paste("Plots/Figure_SX_MoranI_along_distance_of_random_sampled_Remote_Sensing_Residual_with_1000_km_xmax_df_100_SD2.pdf",sep=""),width = 5, height=4)
plot(spatialMoranDynamicsSD2,main="SD2 Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1),xlab="Distance (km)",ylab="Moran's I")
abline(v=725,col="red", lty=2)
dev.off()


set.seed(10000)
trainTable = rawDataTable[sample(nrow(rawDataTable), 40000), ]

formulaString = as.formula(paste('densityWK','~', paste("ti(",trainVaribles,")", collapse="+",sep=""), sep=""))
# formulaString = as.formula(paste('RemoteBiomass','~', paste(trainVaribles, collapse="+",sep=""), sep=""))
gamMulti = mgcv::gam(data = trainTable, formula = formulaString) 
trainTable$gamResiduals = gamMulti$residuals

duplicateTable = trainTable
# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+init=epsg:4326")

# # calculate the Moran I changes along the spatial distance
spatialMoranDynamicsWK2 = spline.correlog(x=coordinates(duplicateTable)[,1], y=coordinates(duplicateTable)[,2],z=duplicateTable$gamResiduals, resamp=20, quiet=TRUE,latlon=T,xmax=1000,df=100,na.rm=T)

pdf(paste("Plots/Figure_SX_MoranI_along_distance_of_random_sampled_Remote_Sensing_Residual_with_1000_km_xmax_df_100_WK2.pdf",sep=""),width = 5, height=4)
plot(spatialMoranDynamicsWK2,main="WK2 Residuals Moran's I",xlim=c(0,1000),ylim=c(-1,1),xlab="Distance (km)",ylab="Moran's I")
abline(v=700,col="red", lty=2)
dev.off()








