# load the library
library(raster)
library(RColorBrewer)
library(rgeos)
library(sp)
library(rgdal)
# setwd("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation")
# load the uncertainty maps for present and potential predictions
GS1_Max_CV = raster("Data/BiomassMergedMaps/GS1_MaxScaler_Potential_Model_Variation_Coefficient_Map_Merged.tif")
GS2_Max_CV = raster("Data/BiomassMergedMaps/GS2_MaxScaler_Potential_Model_Variation_Coefficient_Map_Merged.tif")

GS1_Mean_CV = raster("Data/BiomassMergedMaps/GS1_MeanScaler_Potential_Model_Variation_Coefficient_Map_Merged.tif")
GS2_Mean_CV = raster("Data/BiomassMergedMaps/GS2_MeanScaler_Potential_Model_Variation_Coefficient_Map_Merged.tif")

SD1_CV = raster("Data/BiomassMergedMaps/SD1_Potential_Model_Variation_Coefficient_Map_Merged.tif")
SD2_CV = raster("Data/BiomassMergedMaps/SD2_Potential_Model_Variation_Coefficient_Map_Merged.tif")

HM1_CV = raster("Data/BiomassMergedMaps/HM1_Potential_Model_Variation_Coefficient_Map_Merged.tif")
HM2_CV = raster("Data/BiomassMergedMaps/HM2_Potential_Model_Variation_Coefficient_Map_Merged.tif")

WK1_CV = raster("Data/BiomassMergedMaps/WK1_Potential_Model_Variation_Coefficient_Map_Merged.tif")
WK2_CV = raster("Data/BiomassMergedMaps/WK2_Potential_Model_Variation_Coefficient_Map_Merged.tif")

# print(quantile(GS1_Max_CV,probs=c(0.975)))
# print(quantile(GS2_Max_CV,probs=c(0.975)))

# print(quantile(GS1_Mean_CV,probs=c(0.975)))
# print(quantile(GS2_Mean_CV,probs=c(0.975)))

# print(quantile(HM1_CV,probs=c(0.975)))
# print(quantile(HM2_CV,probs=c(0.975)))

# print(quantile(SD1_CV,probs=c(0.975)))
# print(quantile(SD2_CV,probs=c(0.975)))

# print(quantile(WK1_CV,probs=c(0.975)))
# print(quantile(WK2_CV,probs=c(0.975)))

# 95% quantile at the 97.5%
# GS1_Max_CV:0.2261669 
# GS2_Max_CV:0.2310852 
# GS1_Mean_CV:0.2145732 
# GS2_Mean_CV:0.2278503 
# HM1_CV:0.2539754 
# HM2_CV:0.2515216 
# SD1_CV:0.1358976 
# SD2_CV:0.1297116 
# WK1_CV:0.1491379 
# WK2_CV:0.1527334 

#  mask the pixesl with cover less than 1%
forestCoverPotential = raster("Data/BiomassMergedMaps/Bastin_et_al_2019_Potential_Forest_Cover_Adjusted_Merged.tif")
forestCoverPotential[forestCoverPotential<=0.01] =NA
 
# since there are not too many pixels have the variation higher than 25%,we modifiy the maximum to 0.25
GS1_Max_CV[GS1_Max_CV >=0.3] =0.3
GS2_Max_CV[GS2_Max_CV >=0.3] =0.3

GS1_Mean_CV[GS1_Mean_CV >=0.3] =0.3
GS2_Mean_CV[GS2_Mean_CV >=0.3] =0.3

SD1_CV[SD1_CV >=0.3] =0.3
SD2_CV[SD2_CV >=0.3] =0.3

HM1_CV[HM1_CV >=0.3] =0.3
HM2_CV[HM2_CV >=0.3] =0.3

WK1_CV[WK1_CV >=0.3] =0.3
WK2_CV[WK2_CV >=0.3] =0.3

# mask by the potential cover
GS1_Max_CV_Map = mask(GS1_Max_CV,forestCoverPotential)
GS2_Max_CV_Map = mask(GS2_Max_CV,forestCoverPotential)

GS1_Mean_CV_Map = mask(GS1_Mean_CV,forestCoverPotential)
GS2_Mean_CV_Map = mask(GS2_Mean_CV,forestCoverPotential)

HM1_CV_Map = mask(HM1_CV,forestCoverPotential)
HM2_CV_Map = mask(HM2_CV,forestCoverPotential)

SD1_CV_Map = mask(SD1_CV,forestCoverPotential)
SD2_CV_Map = mask(SD2_CV,forestCoverPotential)

WK1_CV_Map = mask(WK1_CV,forestCoverPotential)
WK2_CV_Map = mask(WK2_CV,forestCoverPotential)


 # define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = projectExtent(WK2_CV_Map, equalEarthProj)
# reproject the rasters
GS1_Max_CV_MapProjected = projectRaster(GS1_Max_CV_Map , to=extentInfo,over=T)
GS2_Max_CV_MapProjected = projectRaster(GS2_Max_CV_Map, to=extentInfo,over=T)

GS1_Mean_CV_MapProjected = projectRaster(GS1_Mean_CV_Map , to=extentInfo,over=T)
GS2_Mean_CV_MapProjected = projectRaster(GS2_Mean_CV_Map, to=extentInfo,over=T)

HM1_CV_MapProjected = projectRaster(HM1_CV_Map , to=extentInfo,over=T)
HM2_CV_MapProjected = projectRaster(HM2_CV_Map, to=extentInfo,over=T)

SD1_CV_MapProjected = projectRaster(SD1_CV_Map , to=extentInfo,over=T)
SD2_CV_MapProjected = projectRaster(SD2_CV_Map, to=extentInfo,over=T)

HM1_CV_MapProjected = projectRaster(HM1_CV_Map , to=extentInfo,over=T)
HM2_CV_MapProjected = projectRaster(HM2_CV_Map, to=extentInfo,over=T)

WK1_CV_MapProjected = projectRaster(WK1_CV_Map , to=extentInfo,over=T)
WK2_CV_MapProjected = projectRaster(WK2_CV_Map, to=extentInfo,over=T)

# pal= colorRampPalette(brewer.pal(6, "YlGnBu"))

pal = colorRampPalette(c("#F7F7F7","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))

worldBorderInit = shapefile("WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj))

pdf("Plots/Fig_S1_Bootstrap_Uncertainties_Variation_Coefficients_all_Models.pdf",width =24, height=40)
par(mfrow=c(5,2))
plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="GS1 Max uncertainty",lwd=0.2)
plot(GS1_Max_CV_MapProjected ,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="GS2 Maxuncertainty",lwd=0.2)
plot(GS2_Max_CV_MapProjected ,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="GS1 Mean uncertainty",lwd=0.2)
plot(GS1_Mean_CV_MapProjected ,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="GS2 Mean uncertainty",lwd=0.2)
plot(GS2_Mean_CV_MapProjected ,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="HM1 uncertainty",lwd=0.2)
plot(HM1_CV_MapProjected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="HM2 uncertainty",lwd=0.2)
plot(HM2_CV_MapProjected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="SD1 uncertainty",lwd=0.2)
plot(SD1_CV_MapProjected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="SD2 uncertainty",lwd=0.2)
plot(SD2_CV_MapProjected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="WK1 uncertainty",lwd=0.2)
plot(WK1_CV_MapProjected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="WK2 uncertainty",lwd=0.2)
plot(WK2_CV_MapProjected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)
dev.off()

