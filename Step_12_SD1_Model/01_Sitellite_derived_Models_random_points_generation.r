library(raster)
library(rgdal)
library(sp)
library(dplyr)
library(data.table)
library(gam)
library(h2o)
# set working directory
setwd("~/Github_Folder")

# STEP 1 read the hansen forest cover
forestCoverMap = raster("SD_RandomPointsGeneration/HansenTreeCoverMap/Hansen_Tree_Cover_Map_merged.tif") # This tif layer was not in the folder, but you can get the map by yourself or ask help from us.
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
writeOGR(replicatedTable,dsn="SD_RandomPointsGeneration/RandomPointsShp",layer="Random_Points_Shp_100000", driver = "ESRI Shapefile",overwrite=T) # this shapefile was uploaded in the GEE assests and you can apply the authentication from us. The path on GEE is :  ee.FeatureCollection("users/leonidmoore/ForestBiomass/RemoteSensingModel/RandomPoints/Random_Points_Shp_100000")



