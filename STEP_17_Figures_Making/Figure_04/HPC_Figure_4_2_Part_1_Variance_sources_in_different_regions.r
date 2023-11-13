# load the library
library(raster)
library(RColorBrewer)
library(rgeos)
library(viridis)
library(sp)
library(rgdal)


referenceRaster = raster("Data/BiomassMergedMaps/Full_TGB_potential_Map_of_ensembled_mean_merged.tif")
# load the SD and GS maps
GS_VarianceID_Stack = stack("Data/BiomassMergedMaps/GS_Max_variance_IDs_Merged.tif")[[2]] #the second layer is the highest variance source
RM_VarianceID_Stack = stack("Data/BiomassMergedMaps/RM_Max_variance_IDs_Merged.tif")[[2]]

 # define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = raster::projectExtent(referenceRaster, equalEarthProj)
# reproject the rasters
# reproject the rasters
GS_varianceMapProjected = projectRaster(GS_VarianceID_Stack, to=extentInfo,over=T)
RM_varianceMapProjected = projectRaster(RM_VarianceID_Stack, to=extentInfo,over=T)

worldBorderInit = shapefile("Data/WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj))
# export the figure out in PDF format
pdf("Plots/Figure_4_Part_1_Variance_sourced_map.pdf",width =12, height=16)
par(mfrow=c(3,1))

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Ground-sourced Map",lwd=0.2)
plot(GS_varianceMapProjected,col=c('#440154','#414487','#2a788e','#22a884','#7ad151','#fde725'),add=T,legend=F,maxpixels = 50000000)

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Satellite-derived Map",lwd=0.2)
plot(RM_varianceMapProjected,col=c('#440154','#414487','#2a788e','#22a884','#7ad151','#fde725'),add=T,legend=F,maxpixels = 50000000)

dev.off()
