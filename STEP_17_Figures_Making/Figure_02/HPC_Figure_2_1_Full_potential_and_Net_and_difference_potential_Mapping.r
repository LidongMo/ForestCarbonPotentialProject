# load the library
library(raster)
library(RColorBrewer)
library(rgeos)
library(viridis)
library(sp)
library(rgdal)

# setwd("/Volumes/Scratch2/Lidong/BiomassEstimation/")
# setwd("~/Desktop/BIOMASS/")
# load the map tiles downloaded from GEE ang cloud
# fullPotentialDensity1 = raster("BiomassMaps/Full_TGB_potential_Map_of_ensembled_mean0000000000-0000000000.tif")
# fullPotentialDensity2 = raster("BiomassMaps/Full_TGB_potential_Map_of_ensembled_mean0000000000-0000023296.tif")
# netPotentialDensity1 = raster("BiomassMaps/Net_TGB_potential_Map_of_ensembled_mean0000000000-0000000000.tif")
# netPotentialDensity2 = raster("BiomassMaps/Net_TGB_potential_Map_of_ensembled_mean0000000000-0000023296.tif")

# fullPotentialDensity = merge(fullPotentialDensity1,fullPotentialDensity2)
# netPotentialDensity = merge(netPotentialDensity1,netPotentialDensity2)

# writeRaster(fullPotentialDensity,"BiomassMergedMaps/Full_TGB_potential_Map_of_ensembled_mean_merged.tif",overwrite=T)
# writeRaster(netPotentialDensity,"BiomassMergedMaps/Net_TGB_potential_Map_of_ensembled_mean_merged.tif",overwrite=T)

openWaterMask = raster("Data/BiomassMergedMaps/Open_Water_mask_Map.tif")
fullPotentialDensity = raster("Data/BiomassMergedMaps/Full_TGB_potential_Map_of_ensembled_mean_merged.tif")
netPotentialDensity = raster("Data/BiomassMergedMaps/Net_TGB_potential_Map_of_ensembled_mean_merged.tif")
# load the merged maps
netDifferenceMap = raster("Data/BiomassMergedMaps/Difference_in_biomass_density_of_GS_and_SD_Net.tif")

potentialCover = raster("Data/BiomassMergedMaps/Bastin_et_al_2019_Potential_Forest_Cover_Adjusted_Merged.tif")
potentialCover[potentialCover<0.05] = NA

netPotentialDensity = mask(netPotentialDensity,potentialCover)
fullPotentialDensity = mask(fullPotentialDensity,potentialCover)
netDiffMasked = mask(netDifferenceMap,potentialCover)

# mask the open water 
openWaterMask[openWaterMask<=0] = NA
netPotentialDensity = mask(netPotentialDensity,openWaterMask)

# adjust the values in the maps
netDiffMasked[netDiffMasked >= 75] = 75
netDiffMasked[netDiffMasked <= -75] = -75

fullPotentialDensity[fullPotentialDensity >= 250] = 250
netPotentialDensity[netPotentialDensity >= 150] =150 

fullPotentialDensity[fullPotentialDensity <=0] = NA
netPotentialDensity[netPotentialDensity <=0] = NA


 # define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = raster::projectExtent(fullPotentialDensity, equalEarthProj)
# reproject the rasters
# reproject the rasters
netDifferenceMapProjected = projectRaster(netDiffMasked, to=extentInfo,over=T)
fullPotentialDensityProjected = projectRaster(fullPotentialDensity, to=extentInfo,over=T)
netPotentialDensityProjected = projectRaster(netPotentialDensity, to=extentInfo,over=T)

worldBorderInit = shapefile("Data/WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj))

wwfBiome = shapefile("Data/WWF_Pol/wwf_terr_ecos.shp")
wwfBiome = crop(wwfBiome,c(-180, 180, -60, 84))
# subset shapefile
subWwfBiome = wwfBiome[wwfBiome$BIOME %in% c(7,8,9),]
aggreBiome = aggregate(subWwfBiome)
aggreBiomeS = gSimplify(aggreBiome, tol = 0.01, topologyPreserve = TRUE)
biomeEqualEarth = sp::spTransform(aggreBiomeS,to=extentInfo, CRS(equalEarthProj))

# define the colour pallete
# pal= colorRampPalette((c("gray35","gray60", "#A1D99B" ,"#74C476", "#41AB5D" ,"#238B45", "#006D2C", "#00441B","#004521")))
pal= colorRampPalette((c("gray80","#A1D99B" ,"#74C476", "#41AB5D" ,"#238B45", "#006D2C", "#00441B","#004521")))
# pal= colorRampPalette((c("#753217","#895025","#9C6F34","#AF9246" ,"#C2B75B", "#B9BF5F",  "#90A349" ,"#6D8A3A", "#4B732A", "#2E5C1C" ,"#10450F")))
# pal= colorRampPalette(brewer.pal(5, "YlGn"))

pdf("Plots/Fig_2_Restoration_Abs_Reforestation_conservation_20230427.pdf",width =12, height=16)
par(mfrow=c(3,1))

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Absolute Poltential Map",lwd=0.2)
plot(fullPotentialDensityProjected,col=pal(100),add=T,breaks=c(seq(0,250,2.5)),legend=T,maxpixels = 50000000)
plot(biomeEqualEarth,density=25,angle=45,lwd=0.5,border="gray55",add=T)

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Difference Poltential Map",lwd=0.2)
plot(netPotentialDensityProjected,col=pal(100),add=T,breaks=c(seq(0,150,1.5)),legend=T,maxpixels = 50000000)
plot(biomeEqualEarth,density=25,angle=45,lwd=0.5,border="gray55",add=T)

pal= colorRampPalette(c('#A50026', '#D73027', '#F46D43', '#FDAE61', '#FEE090', '#FFFFFF','#E0F3F8', '#ABD9E9', '#74ADD1', '#4575B4', '#313695'))

plot(worldEqualEarth,border="gray30",col="gray30",axes=T,main="Difference comparison (GSs-SDs) Net potential",lwd=0.2)
plot(netDifferenceMapProjected,col=pal(100),add=T,breaks=c(seq(-75,75,1.5)),legend=T,maxpixels = 50000000)
plot(biomeEqualEarth,density=25,angle=45,lwd=0.5,border="gray55",add=T)

dev.off()
