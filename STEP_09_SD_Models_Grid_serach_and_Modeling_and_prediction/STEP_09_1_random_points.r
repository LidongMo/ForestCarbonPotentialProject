library(raster)
library(rgdal)
library(sp)
library(dplyr)
library(data.table)

# Part 1, data preperation for Type 1 model
##########################################################################################################################

# load the forest cover map
forestCoverMap = raster("Data/SatelliteDerivedModel/HansenTreeCoverMap/Hansen_Tree_Cover_Map_merged.tif")
# do a random samppling to get 100,0000 points
randomPoints = as.data.frame(sampleRandom(forestCoverMap, 2000000, na.rm=T, xy=T))
# get 100,0000 points 
set.seed(1000)
randomPoints = randomPoints[sample(nrow(randomPoints), 1000000), ]
replicatedTable = randomPoints
coordinates(replicatedTable) = ~x+y
# allocate projection system
proj4string(replicatedTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
replicatedTable@data = randomPoints[,c("x","y")]
# write the shapefile into local folder
writeOGR(replicatedTable,dsn="Data/SatelliteDerivedModel/RandomPointsShp",layer="Random_Points_Shp_100000", driver = "ESRI Shapefile",overwrite=T)

# Part 2, data preperation for Type 1 model
##########################################################################################################################
# load the full data and filter by the WDPA and intact forest 
fullTrainTable = fread("Data/SatelliteDerivedModel/CovariatesExtracted/20221209_RemoteSensing_Biomass_Merged_sampled_dataset_1000000.csv")  %>% filter(WDPA==1|IntactForest==1) %>% dplyr::select(x,y,CHELSA_Annual_Precipitation) %>% na.omit()  #
# transfer this into shapefile
replicatedTable = fullTrainTable
coordinates(replicatedTable) = ~x+y
# allocate projection system
proj4string(replicatedTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
replicatedTable@data = fullTrainTable[,c("x","y")]


for (seeds in 0:99)
{
    set.seed(seeds)
    gridSubsampledPoints = sample.grid(replicatedTable, cell.size = c(1,1), n =1)$subset

    writeOGR(gridSubsampledPoints,dsn="Data/SatelliteDerivedModel/GridSubsampledShapfiles",layer= paste("SD2_Natural_WDPA_Intact_Gridsampled_",seeds,sep=""), driver = "ESRI Shapefile",overwrite=T)
}

# transfer the data in google cloud to google eaerh engine 
for ((i=0; i<100; i++)); 
do
    earthengine upload table --asset_id=users/leonidmoore/ForestBiomass/RemoteSensingModel/GridSampleShapefiles/SD2_Gridsubsampled_Natural_Seed_${i} gs://crowtherlab_gcsb_lidong/GridSubsampledSD2/SD2_Natural_WDPA_Intact_Gridsampled_${i}.shp
done


