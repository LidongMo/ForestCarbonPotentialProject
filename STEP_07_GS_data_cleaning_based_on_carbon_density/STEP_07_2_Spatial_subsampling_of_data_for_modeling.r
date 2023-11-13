# load the library
library(data.table)
library(dplyr)
library(raster)
library(rgdal)
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


# load the data table
tableWithScaler = fread("Data/GroundSourcedModel/20230126_Merged_Covariates_sampled_dataset_outliers_cleaned.csv")[,-1]
# calculate a column with the scaled biomass density
tableWithScaler$MeanScaledDensity = tableWithScaler$CarbonDensity*tableWithScaler$MeanCoverScaler 
tableWithScaler$LogMeanScaledDensity = log(tableWithScaler$MeanScaledDensity+1)

tableWithScaler$MaxScaledDensity = tableWithScaler$CarbonDensity*tableWithScaler$MaxCoverScaler
tableWithScaler$LogMaxScaledDensity = log(tableWithScaler$MaxScaledDensity+1)

# trasfer this into shape file for GEE uploading
processTable = tableWithScaler %>% dplyr::select(x,y,LogMeanScaledDensity,LogMaxScaledDensity,MeanScaledDensity,MaxScaledDensity,PlotAre,WWF_Biome)
# write csv
write.csv(processTable,"Data/GroundSourcedModel/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution_BD_Transformed.csv")


#for Model type 1 model
# conduct the spatial subsampling
tableWithScaler= fread("Data/GroundSourcedModel/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution_BD_Transformed.csv")[,-1] %>% dplyr::select(x,y,LogMeanScaledDensity,LogMaxScaledDensity) %>% na.omit()
processTable = tableWithScaler
# # tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(processTable) = ~ x + y
# allocate the projection
proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
processTable@data = tableWithScaler

for (seeds in 0:99)
{
    set.seed(seeds)
    gridSubsampledPoints = sample.grid(processTable, cell.size = c(0.65,0.65), n =1)$subset

    writeOGR(gridSubsampledPoints,dsn="Data/GroundSourcedModel/GridSubsampledShapfiles",layer= paste("GS1_allScalers_GFBI_Gridsampled_0427_",seeds,sep=""), driver = "ESRI Shapefile",overwrite=T)
}

# upload manually to cloud storage and upload to assets by bash code below

for ((i=0; i<100; i++)); 
do
    earthengine upload table --asset_id=users/leonidmoore/ForestBiomass/GroundSourcedModel/GridSubShapefiles/GS1_allScalers_GFBI_Gridsampled_0427_$i gs://crowtherlab_gcsb_lidong/GridSubsampledGS1/GS1_allScalers_GFBI_Gridsampled_0427_$i.shp
done





#for Model type 2 model

# load the data table
tableWithScaler = fread("Data/GroundSourcedModel/20230126_Merged_Covariates_sampled_dataset_outliers_cleaned.csv")[,-1]
# calculate a column with the scaled biomass density
tableWithScaler$MeanScaledDensity = tableWithScaler$CarbonDensity*tableWithScaler$MeanCoverScaler 
tableWithScaler$LogMeanScaledDensity = log(tableWithScaler$MeanScaledDensity+1)

tableWithScaler$MaxScaledDensity = tableWithScaler$CarbonDensity*tableWithScaler$MaxCoverScaler
tableWithScaler$LogMaxScaledDensity = log(tableWithScaler$MaxScaledDensity+1)

# trasfer this into shape file for GEE uploading
processTable = tableWithScaler %>% filter(IntactForest==1|WDPA==1|Human_Disturbance<=0.01) %>% dplyr::select(x,y,LogMeanScaledDensity,LogMaxScaledDensity,MeanScaledDensity,MaxScaledDensity,PlotAre,WWF_Biome) 

# write csv
write.csv(processTable,"Data/GroundSourcedModel/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution_BD_Transformed_GS2.csv")


tableWithScaler= fread("Data/GroundSourcedModel/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution_BD_Transformed_GS2.csv")[,-1] %>% dplyr::select(x,y,LogMeanScaledDensity,LogMaxScaledDensity) %>% na.omit()
processTable = tableWithScaler
# # tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(processTable) = ~ x + y
# allocate the projection
proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
processTable@data = tableWithScaler
for (seeds in 0:99)
{
    set.seed(seeds)
    gridSubsampledPoints = sample.grid(processTable, cell.size = c(0.25,0.25), n =1)$subset

    writeOGR(gridSubsampledPoints,dsn="Data/GroundSourcedModel/GridSubsampledShapfiles",layer= paste("GS2_allScalers_GFBI_Natural_Gridsampled_0427_",seeds,sep=""), driver = "ESRI Shapefile",overwrite=T)
}

# upload manually to cloud storage and upload to assets by bash code below

for ((i=0; i<100; i++)); 
do
    earthengine upload table --asset_id=users/leonidmoore/ForestBiomass/GroundSourcedModel/GridSubShapefiles/GS2_allScalers_GFBI_Natural_Gridsampled_0427_$i gs://crowtherlab_gcsb_lidong/GridSubsampledGS2/GS2_allScalers_GFBI_Natural_Gridsampled_0427_${i}.shp
done


# prepare the data for the figure 2
rawData1 = fread(paste("Data/ForestCoverScaler/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution.csv",sep=""))[,-1] 
# adaptation of the far east russian points
farEastRussiaPlots = rawData1 %>% dplyr::filter(x>136&y>58) %>% dplyr::mutate(x= x*(-1))#229

otherPlots = rawData1 %>% dplyr::filter(x<=136|y<=58) #563852
# rbind the modified tables together
tableWithScaler = rbind(farEastRussiaPlots,otherPlots)
# add the carbon density row
tableWithScaler$CarbonDensity = tableWithScaler$BmssDns/1000*tableWithScaler$CarbonConcentration 

tableWithScaler = tableWithScaler%>% dplyr::select(x,y,CarbonDensity) %>% na.omit()
processTable = tableWithScaler
# # tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(processTable) = ~ x + y
# allocate the projection
proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
processTable@data = tableWithScaler
# write to local folder
writeOGR(processTable,dsn="Data/ForestCoverScale/ShapefileFull",layer= "Full_GFBI_AllYears_NewScaled_with_outliers_0427", driver = "ESRI Shapefile",overwrite=T)

