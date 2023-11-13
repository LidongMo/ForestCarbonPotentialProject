# Figure 1

library(data.table)
library(parallel)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(GSIF)
library(rgdal)
library(RColorBrewer)
library(rgeos)
library(gridExtra)
library(effects)
library(ggpubr)
library(tidyverse)
library(stringr)
library(sf)
library(sp)

####################################
########## Panel A for Figure 1
####################################
# load the table for training data
inputTable =fread("Data/GroundSourcedModel/20230126_Merged_Covariates_sampled_dataset_outliers_cleaned_for_Figure_0427.csv")[,-1] %>% dplyr::select(x,y,CrbnDns)
inputTable$mod = inputTable$CrbnDns
# change all the values larger than 500 to 500
inputTable = within(inputTable, mod[CrbnDns>=250] <- 250)

# load the biomass map for panel b
remoteBiomassImage = raster("Data/SupplementaryFigures/Figure_01/SD_Arnan_carbon_stock_density_mean_Map_Merged.tif")

 # projInfo()
# load the world border polygon
worldBorder = shapefile("Data/WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorder,c(-180, 180, -60, 84)) 
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)


# define the colour pallete
pal= colorRampPalette((c("gray80","#A1D99B" ,"#74C476", "#41AB5D" ,"#238B45", "#006D2C", "#00441B","#004521")))

# define the equal earth projection
equalEarthProj =  "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# project raster
extentInfo = projectExtent(remoteBiomassImage,equalEarthProj)
rasEqualEarth = projectRaster(remoteBiomassImage, to=extentInfo,over=T)

# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj))

# define the colour for the points
inputTable$Col = pal(100)[as.numeric(cut(inputTable$mod,breaks = 100))]
duplicateTable = inputTable
# transform the data into spatial points
coordinates(duplicateTable) = ~x+y
proj4string(duplicateTable) = CRS("+proj=longlat")
transformedPoints = spTransform(duplicateTable, to=extentInfo,CRS(equalEarthProj))


pdf("Plots/Figure_1_Panel_a_20230427.pdf",width =16, height=20)
par(mfrow=c(3,1))
# panel a
plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Plots Map",lwd=0.2)
points(transformedPoints,col=inputTable$Col,pch=20,cex=1.25)
dev.off()

# panel b 

pdf("Plots/Fig_1_Panel_b_20230427.pdf",width =16, height=20)
par(mfrow=c(3,1))
#This adds a column of color values
plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Presen biomass density map (ESA_CCI)",lwd=0.2)
plot(rasEqualEarth,col=(pal(100)),add=T,breaks=c(seq(0,250,2.5)),legend=T,maxpixels = 50000000)
dev.off()



# this is code we used to generate the legend for panel A has been made in other script by raster maps in fig 4
####################################
########## Panel B for Figure 1
####################################
tableWithScaler = fread("Data/GroundSourcedModel/20230126_Merged_Covariates_sampled_dataset_outliers_cleaned_for_Figure_0427.csv")[,-1] %>% na.omit()

# to avoid the uneven distribution of the data points, we applied spatial subsampling 
tableWithScaler$MeanScaledDensity = tableWithScaler$CrbnDns*tableWithScaler$MeanCoverScaler
tableWithScaler$MaxScaledDensity = tableWithScaler$CrbnDns*tableWithScaler$MaxCoverScaler
# get the mean of the two densities
tableWithScaler$MeanDensity = (tableWithScaler$MaxScaledDensity +tableWithScaler$MeanScaledDensity)/2
# limited number of biome 3 and 14, we did not do the boxplot for that two biomes

pal= colorRampPalette(brewer.pal(6, "Set2"))
inputTableGS = tableWithScaler[tableWithScaler$WWF_Bim<=14,] 
inputTableGS = inputTableGS%>% mutate(MeanDensity = replace(MeanDensity, WWF_Bim ==3|WWF_Bim ==14, 0)) 
inputTableGS$WWF_Bim = as.factor(inputTableGS$WWF_Bim)

# pdf("BiomeLevelEstimation/plotsForPaper/Present_bioma_biomass_boxplot.pdf")
panelC = ggplot(inputTableGS, aes(x=WWF_Bim, y=MeanDensity,group=WWF_Bim,fill=WWF_Bim)) + 
        geom_boxplot(outlier.shape = NA) +
        theme_classic()+
        scale_y_continuous(expand = c(0,0),limits = c(0, 250), breaks = 2*c(0,50,100))+
        guides(fill=FALSE) +
        scale_fill_manual(values =pal(14))+
        labs(x="Biome", y = "Forest tree carbon density (tC/ha)")+
        theme(panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.title.x = element_blank(),
        axis.text.x=element_blank(),
        axis.title.y = element_text(size = 15),
        axis.ticks.length = unit(0.2, "cm"))

####################################
########## Panel C for Figure 1
####################################
# merge all the subtables for the bootstrapping
tableWithScaler = fread("Data/GroundSourcedModel/20230126_Merged_Covariates_sampled_dataset_outliers_cleaned_for_Figure_0427.csv")[,-1]  %>% na.omit()
# to avoid the uneven distribution of the data points, we applied spatial subsampling 
tableWithScaler$MeanScaledDensity = tableWithScaler$CrbnDns*tableWithScaler$MeanCoverScaler
tableWithScaler$MaxScaledDensity = tableWithScaler$CrbnDns*tableWithScaler$MaxCoverScaler

processTable = tableWithScaler
# # tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(processTable) = ~ x + y
# allocate the projection
proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
processTable@data = tableWithScaler

set.seed(10000)
gridSubsampledPoints = sample.grid(processTable, cell.size = c(0.1,0.1), n =1)$subset


pcaDataTable  = gridSubsampledPoints@data
# define the variable 
retainVaribles = c('Aridity_Index','CHELSA_Annual_Mean_Temperature','CHELSA_Annual_Precipitation','CHELSA_Isothermality','CHELSA_Max_Temperature_of_Warmest_Month','CHELSA_Mean_Diurnal_Range','CHELSA_Mean_Temperature_of_Coldest_Quarter','CHELSA_Mean_Temperature_of_Driest_Quarter','CHELSA_Mean_Temperature_of_Warmest_Quarter','CHELSA_Mean_Temperature_of_Wettest_Quarter','CHELSA_Min_Temperature_of_Coldest_Month','CHELSA_Precipitation_Seasonality','CHELSA_Precipitation_of_Coldest_Quarter','CHELSA_Precipitation_of_Driest_Month','CHELSA_Precipitation_of_Driest_Quarter','CHELSA_Precipitation_of_Warmest_Quarter','CHELSA_Precipitation_of_Wettest_Month','CHELSA_Precipitation_of_Wettest_Quarter','CHELSA_Temperature_Annual_Range','CHELSA_Temperature_Seasonality','Depth_to_Water_Table','EarthEnvTopoMed_Eastness','EarthEnvTopoMed_Elevation','EarthEnvTopoMed_Northness','EarthEnvTopoMed_ProfileCurvature','EarthEnvTopoMed_Roughness','EarthEnvTopoMed_Slope','SG_Absolute_depth_to_bedrock','WorldClim2_SolarRadiation_AnnualMean','WorldClim2_WindSpeed_AnnualMean','EarthEnvCloudCover_MODCF_interannualSD','EarthEnvCloudCover_MODCF_intraannualSD','EarthEnvCloudCover_MODCF_meanannual','EarthEnvTopoMed_AspectCosine','EarthEnvTopoMed_AspectSine','SG_Clay_Content_0_100cm','SG_Coarse_fragments_0_100cm','SG_Sand_Content_0_100cm','SG_Silt_Content_0_100cm','SG_Soil_pH_H2O_0_100cm')

# pcaDataTable$BmssDns = log(pcaDataTable$BmssDns +1)
# kick out the na values in the data frame
# human activitiy variable table
scaledPCATable = pcaDataTable[,c('LandCoverClass_Cultivated_and_Managed_Vegetation','Human_Disturbance','LandCoverClass_Urban_Builtup','cropland','grazing','pasture','rangeland','WDPA')]
pcaResult = prcomp(scaledPCATable,center = TRUE,scale. = TRUE)
summaryPCA1 = summary(pcaResult)
# plot(pcaResult, type = "l")
library(factoextra)

panelD = fviz_pca_var(pcaResult,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             ggtheme= theme_classic()) +  
             theme(legend.position = "none",
             panel.border = element_rect(colour = "black", fill=NA, size=1.5),
             axis.ticks.length = unit(0.2, "cm"),
             axis.title = element_text(size = 15),
             axis.text = element_text(size = 12))+
             labs(title = NULL,x = str_c("PC1 (",round(summaryPCA1$importance[2,1]*100,1),"%)"),y = str_c("PC2 (",round(summaryPCA1$importance[2,2]*100,1),"%)"))+scale_y_reverse()

# combine the table with the pca results
combinedPCATableGS = data.frame(subset(pcaDataTable,select = retainVaribles),
                                RelativeDensityMeanScaler = ((pcaDataTable$MeanScaledDensity)/mean(pcaDataTable$MeanScaledDensity))-1,
                                RelativeDensityMaxScaler = ((pcaDataTable$MaxScaledDensity)/mean(pcaDataTable$MaxScaledDensity))-1,
                                PC1 = scales::rescale(pcaResult$x[,1])) 
# define the independent variables for formula making
independentVars = c('Aridity_Index','CHELSA_Annual_Mean_Temperature','CHELSA_Annual_Precipitation','CHELSA_Isothermality','CHELSA_Max_Temperature_of_Warmest_Month','CHELSA_Mean_Diurnal_Range','CHELSA_Mean_Temperature_of_Coldest_Quarter','CHELSA_Mean_Temperature_of_Driest_Quarter','CHELSA_Mean_Temperature_of_Warmest_Quarter','CHELSA_Mean_Temperature_of_Wettest_Quarter','CHELSA_Min_Temperature_of_Coldest_Month','CHELSA_Precipitation_Seasonality','CHELSA_Precipitation_of_Coldest_Quarter','CHELSA_Precipitation_of_Driest_Month','CHELSA_Precipitation_of_Driest_Quarter','CHELSA_Precipitation_of_Warmest_Quarter','CHELSA_Precipitation_of_Wettest_Month','CHELSA_Precipitation_of_Wettest_Quarter','CHELSA_Temperature_Annual_Range','CHELSA_Temperature_Seasonality','Depth_to_Water_Table','EarthEnvTopoMed_Eastness','EarthEnvTopoMed_Elevation','EarthEnvTopoMed_Northness','EarthEnvTopoMed_ProfileCurvature','EarthEnvTopoMed_Roughness','EarthEnvTopoMed_Slope','SG_Absolute_depth_to_bedrock','WorldClim2_SolarRadiation_AnnualMean','WorldClim2_WindSpeed_AnnualMean','EarthEnvCloudCover_MODCF_interannualSD','EarthEnvCloudCover_MODCF_intraannualSD','EarthEnvCloudCover_MODCF_meanannual','EarthEnvTopoMed_AspectCosine','EarthEnvTopoMed_AspectSine','SG_Clay_Content_0_100cm','SG_Coarse_fragments_0_100cm','SG_Sand_Content_0_100cm','SG_Silt_Content_0_100cm','SG_Soil_pH_H2O_0_100cm','PC1')

formulaStringMean = as.formula(paste('RelativeDensityMeanScaler', paste(independentVars, collapse="+"), sep="~"))
lmMultiMean = lm(data = combinedPCATableGS, formula = formulaStringMean) 

PC1.effMean = effect(term = "PC1", mod = lmMultiMean,focal.level=10)

effectsTableGS_Mean = data.frame(
  PC1.effMean$x,
  eff.lower = PC1.effMean$lower,
  eff.upper = PC1.effMean$upper,
  eff.fit = PC1.effMean$fit,
  Type = "MeanScaler"
)


formulaStringMax = as.formula(paste('RelativeDensityMaxScaler', paste(independentVars, collapse="+"), sep="~"))
lmMultiMax = lm(data = combinedPCATableGS, formula = formulaStringMax) 

PC1.effMax = effect(term = "PC1", mod = lmMultiMax,focal.level=10)

effectsTableGS_Max = data.frame(
  PC1.effMax$x,
  eff.lower = PC1.effMax $lower,
  eff.upper = PC1.effMax$upper,
  eff.fit = PC1.effMax$fit,
  Type = "MaxScaler"
)


effectsTableGS = rbind(effectsTableGS_Max,effectsTableGS_Mean)
# library(ppcor)
# pcorResutl = pcor(scale(combinedPCATableGS))
# pcorResutl$estimate #

# lm.with<-lm(data = combinedPCATableGS, formula = formulaString1) 
# lm.without<-update(lm.with, ~. - PC1)


panelE = ggplot(data=effectsTableGS, aes(x = PC1,group=Type)) +
        geom_ribbon(aes(ymin = eff.lower, ymax = eff.upper,fill=factor(Type)), alpha = 0.3) +
        geom_line(aes(y = eff.fit,colour=factor(Type)))+
        scale_color_manual(values=c('red', 'blue'),labels = c("Upper canopy cover", "Lower canopy cover"))+
        scale_fill_manual(values=c('red', 'blue'),labels = c("Upper canopy cover", "Lower canopy cover"))+
        theme_classic()+
        guides(colour = guide_legend(title = "Type"),fill = "none")+
        theme(legend.position = c(0.7, 0.85),
        legend.title=element_text(size=15),
        legend.text=element_text(size=12),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.title.x=element_blank(),
        axis.ticks.length = unit(0.2, "cm"))+
        ylim(c(-1.1,0.75))+
        labs(x="Human disturbance", y = "Relative carbon density (GS)")


####################################
########## Panel E for Figure 1
####################################
library(data.table)
library(parallel)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(GSIF)
library(rgdal)
library(RColorBrewer)
library(rgeos)
library(gridExtra)
library(effects)
library(ggpubr)
library(factoextra)
library(ppcor)
# module load java/1.8.0_91
setwd("~/Desktop/BIOMASS/")

# load the data table for PCA analysis
trainDataTable = fread("Data/SatelliteDerivedModel/CovariatesExtracted/20221209_RemoteSensing_Biomass_Merged_sampled_dataset_1000000.csv")[,-1] %>% mutate(RemoteBiomass= (densitySD+densityHM+densityWK)/3)
trainDataTable = na.omit(trainDataTable)
trainDataTableSub = trainDataTable[trainDataTable$WWF_Biome<=14,]
trainDataTableSub$WWF_Biome = as.factor(trainDataTableSub$WWF_Biome)
pal= colorRampPalette(brewer.pal(6, "Set2"))
# pdf("BiomeLevelEstimation/plotsForPaper/Present_bioma_biomass_boxplot.pdf")
panelF = ggplot(trainDataTableSub, aes(x=WWF_Biome, y=RemoteBiomass,group=WWF_Biome,fill=WWF_Biome)) + 
        geom_boxplot(outlier.shape = NA) +
        theme_classic()+
        scale_y_continuous(expand = c(0,0),limits = c(0, 250), breaks = 2*c(0,50,100))+
        guides(fill=FALSE) +
        scale_fill_manual(values =pal(14))+
        # scale_x_discrete(labels=c('Tropical moist','Tropical dry','Tropical coniferous','Temperate broadleaf','Temperate coniferous','Boreal','Tropical savanna','Temperate savanna','Flooded savanna','Montane grassland','Tundra','Mediterranean forest','Desert','Mangroves'))+
        theme(legend.position = "none",panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12, angle=45,vjust=0.95, hjust=1),
        axis.ticks.length = unit(0.2, "cm"))+
        labs(x="Biome", y = "Forest tree carbon density (tC/ha)")
        # +geom_jitter(shape=16, position=position_jitter(0.2))

# load the data table for PCA analysis
trainDataTable = fread("Data/SatelliteDerivedModel/CovariatesExtracted/20221209_RemoteSensing_Biomass_Merged_sampled_dataset_1000000.csv")[,-1] %>% mutate(RemoteBiomass= (densitySD+densityHM+densityWK)/3) %>% na.omit()

processTable = trainDataTable
# # tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(processTable) = ~ x + y
# allocate the projection
proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
processTable@data = trainDataTable

set.seed(10000)
gridSubsampledPoints = sample.grid(processTable, cell.size = c(0.5,0.5), n =1)$subset


pcaDataTable  = gridSubsampledPoints@data

# human activitiy variable table
scaledPCATable = pcaDataTable[,c('LandCoverClass_Cultivated_and_Managed_Vegetation','Human_Disturbance','LandCoverClass_Urban_Builtup','cropland','grazing','pasture','rangeland','WDPA')]
pcaResult = prcomp(scaledPCATable,center = TRUE,scale. = TRUE) 

# define the variable 
retainVaribles = c('LandCoverClass_Cultivated_and_Managed_Vegetation','Human_Disturbance','LandCoverClass_Urban_Builtup','cropland','grazing','pasture','rangeland','WDPA') #,"PresentTreeCover"
# human activitiy variable table
PCA_DataTable = subset(pcaDataTable,select=retainVaribles) *(-1)
pcaResultRM = prcomp(PCA_DataTable,center = TRUE,scale. = TRUE) 
# plot(pcaResult, type = "l")
summaryPCA = summary(pcaResultRM)
summaryPCA$importance[2,1]
panelG = fviz_pca_var(pcaResultRM,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE,
             ggtheme= theme_classic()) +  
             theme(legend.position = "none",
             panel.border = element_rect(colour = "black", fill=NA, size=1.5),
             axis.ticks.length = unit(0.2, "cm"),
             axis.title = element_text(size = 15),
             axis.text = element_text(size = 12))+
             labs(title = NULL,x = str_c("PC1 (",round(summaryPCA$importance[2,1]*100,1),"%)"),y = str_c("PC2 (",round(summaryPCA$importance[2,2]*100,1),"%)"))

####################################
########## Panel F for Figure 1
####################################
# define the environmental variables
envVars = c('Aridity_Index','CHELSA_Annual_Mean_Temperature','CHELSA_Annual_Precipitation','CHELSA_Isothermality','CHELSA_Max_Temperature_of_Warmest_Month','CHELSA_Mean_Diurnal_Range','CHELSA_Mean_Temperature_of_Coldest_Quarter','CHELSA_Mean_Temperature_of_Driest_Quarter','CHELSA_Mean_Temperature_of_Warmest_Quarter','CHELSA_Mean_Temperature_of_Wettest_Quarter','CHELSA_Min_Temperature_of_Coldest_Month','CHELSA_Precipitation_Seasonality','CHELSA_Precipitation_of_Coldest_Quarter','CHELSA_Precipitation_of_Driest_Month','CHELSA_Precipitation_of_Driest_Quarter','CHELSA_Precipitation_of_Warmest_Quarter','CHELSA_Precipitation_of_Wettest_Month','CHELSA_Precipitation_of_Wettest_Quarter','CHELSA_Temperature_Annual_Range','CHELSA_Temperature_Seasonality','Depth_to_Water_Table','EarthEnvTopoMed_Eastness','EarthEnvTopoMed_Elevation','EarthEnvTopoMed_Northness','EarthEnvTopoMed_ProfileCurvature','EarthEnvTopoMed_Roughness','EarthEnvTopoMed_Slope','SG_Absolute_depth_to_bedrock','WorldClim2_SolarRadiation_AnnualMean','WorldClim2_WindSpeed_AnnualMean','EarthEnvCloudCover_MODCF_interannualSD','EarthEnvCloudCover_MODCF_intraannualSD','EarthEnvCloudCover_MODCF_meanannual','EarthEnvTopoMed_AspectCosine','EarthEnvTopoMed_AspectSine','SG_Clay_Content_0_100cm','SG_Coarse_fragments_0_100cm','SG_Sand_Content_0_100cm','SG_Silt_Content_0_100cm','SG_Soil_pH_H2O_0_100cm')
# combine the table with the pca results
combinedPCATableRM = data.frame(subset(pcaDataTable,select = envVars),
                    RelativeDensitySD = pcaDataTable$densitySD/mean(pcaDataTable$densitySD)-1,
                    RelativeDensityWK = pcaDataTable$densityWK/mean(pcaDataTable$densityWK)-1,
                    RelativeDensityHM = pcaDataTable$densityHM/mean(pcaDataTable$densityHM)-1,
                    PC1 = scales::rescale(pcaResultRM$x[,1]*(-1))) # since the pca 1 shows the lower value indicates higher human distrubance,then multiply by -1 to get the reversed pca1
# define the independent variables for formula making
predictorVars = c('Aridity_Index','CHELSA_Annual_Mean_Temperature','CHELSA_Annual_Precipitation','CHELSA_Isothermality','CHELSA_Max_Temperature_of_Warmest_Month','CHELSA_Mean_Diurnal_Range','CHELSA_Mean_Temperature_of_Coldest_Quarter','CHELSA_Mean_Temperature_of_Driest_Quarter','CHELSA_Mean_Temperature_of_Warmest_Quarter','CHELSA_Mean_Temperature_of_Wettest_Quarter','CHELSA_Min_Temperature_of_Coldest_Month','CHELSA_Precipitation_Seasonality','CHELSA_Precipitation_of_Coldest_Quarter','CHELSA_Precipitation_of_Driest_Month','CHELSA_Precipitation_of_Driest_Quarter','CHELSA_Precipitation_of_Warmest_Quarter','CHELSA_Precipitation_of_Wettest_Month','CHELSA_Precipitation_of_Wettest_Quarter','CHELSA_Temperature_Annual_Range','CHELSA_Temperature_Seasonality','Depth_to_Water_Table','EarthEnvTopoMed_Eastness','EarthEnvTopoMed_Elevation','EarthEnvTopoMed_Northness','EarthEnvTopoMed_ProfileCurvature','EarthEnvTopoMed_Roughness','EarthEnvTopoMed_Slope','SG_Absolute_depth_to_bedrock','WorldClim2_SolarRadiation_AnnualMean','WorldClim2_WindSpeed_AnnualMean','EarthEnvCloudCover_MODCF_interannualSD','EarthEnvCloudCover_MODCF_intraannualSD','EarthEnvCloudCover_MODCF_meanannual','EarthEnvTopoMed_AspectCosine','EarthEnvTopoMed_AspectSine','SG_Clay_Content_0_100cm','SG_Coarse_fragments_0_100cm','SG_Sand_Content_0_100cm','SG_Silt_Content_0_100cm','SG_Soil_pH_H2O_0_100cm','PC1')

formulaStringSD = as.formula(paste('RelativeDensitySD', paste(predictorVars, collapse="+"), sep="~"))
lmMultiSD = lm(data = combinedPCATableRM, formula = formulaStringSD) 

PC1.effSD = effect(term = "PC1", mod = lmMultiSD,focal.level=10)

effectsTableSD = data.frame(PC1.effSD $x,
                          eff.lower = PC1.effSD $lower,
                          eff.upper = PC1.effSD$upper,
                          eff.fit = PC1.effSD$fit,
                          Type = 'SD')

formulaStringHM = as.formula(paste('RelativeDensityHM', paste(predictorVars, collapse="+"), sep="~"))
lmMultiHM = lm(data = combinedPCATableRM, formula = formulaStringHM) 

PC1.effHM = effect(term = "PC1", mod = lmMultiHM,focal.level=10)

effectsTableHM = data.frame(PC1.effHM $x,
                          eff.lower = PC1.effHM $lower,
                          eff.upper = PC1.effHM$upper,
                          eff.fit = PC1.effHM$fit,
                          Type = 'HM')

formulaStringWK = as.formula(paste('RelativeDensityWK', paste(predictorVars, collapse="+"), sep="~"))
lmMultiWK = lm(data = combinedPCATableRM, formula = formulaStringWK) 

PC1.effWK = effect(term = "PC1", mod = lmMultiWK,focal.level=10)

effectsTableWK = data.frame(PC1.effWK $x,
                          eff.lower = PC1.effWK$lower,
                          eff.upper = PC1.effWK$upper,
                          eff.fit = PC1.effWK$fit,
                          Type = 'WK')

effectsTableRM = rbind(effectsTableSD,effectsTableHM,effectsTableWK)


panelH = ggplot(data=effectsTableRM, aes(x = PC1,group=Type)) +
        geom_ribbon(aes(ymin = eff.lower, ymax = eff.upper,fill=factor(Type)), alpha = 0.3) +
        geom_line(aes(y = eff.fit,colour=factor(Type)))+
        scale_color_manual(values=c('red', 'blue','darkcyan'),labels = c("ESA_CCI", "Harmonized","Walker et al. 2022"))+
        scale_fill_manual(values=c('red', 'blue','darkcyan'),labels = c("ESA_CCI", "Harmonized","Walker et al. 2022"))+
        theme_classic()+
        guides(colour = guide_legend(title = "Type"),fill = "none")+
        theme(legend.position = c(0.7, 0.85),legend.title=element_text(size=15),legend.text=element_text(size=12),axis.title = element_text(size = 15),axis.text = element_text(size = 12),panel.border = element_rect(colour = "black", fill=NA, size=1.5),axis.ticks.length = unit(0.2, "cm"))+
        ylim(c(-1.4,0.75))+
        labs(x="Human disturbance", y = "Relative carbon density (SD)")


pdf("Plots/Fig_1_Panel_c_d_e_f_g_h.pdf",width =15, height=10)
library("gridExtra")
ggarrange(panelC,panelD,panelE,panelF,panelG,panelH,align="hv", widths = c(1,1,1),heights = c(1,1),ncol=3,nrow=2)

dev.off()