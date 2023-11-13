# load the library
library(raster)
library(RColorBrewer)
library(rgeos)
library(sp)
library(rgdal)

# load the uncertainty maps for present and potential predictions
GS_Max_CV = raster("Data/BiomassMergedMaps/GS_Max_Present_AGB_Modeling_VC_Map_Merged.tif")

GS_Mean_CV = raster("DataBiomassMergedMaps/GS_Mean_Present_AGB_Modeling_VC_Map_Merged.tif")

# print(quantile(GS_Max_CV,probs=c(0.975)))

# print(quantile(GS_Mean_CV,probs=c(0.975)))

# 0.2947121 

# 0.2804561 

#  mask the pixesl with cover less than 1%
forestCoverPotential = raster("Data/BiomassMergedMaps/Bastin_et_al_2019_Potential_Forest_Cover_Adjusted_Merged.tif")
forestCoverPotential[forestCoverPotential<=0.01] =NA

# since there are not too many pixels have the variation higher than 25%,we modifiy the maximum to 0.25
GS_Max_CV[GS_Max_CV >=0.3] =0.3
GS_Mean_CV[GS_Mean_CV >=0.3] =0.3

# mask by the potential cover
GS_Max_CV_Map = mask(GS_Max_CV,forestCoverPotential)

GS_Mean_CV_Map = mask(GS_Mean_CV,forestCoverPotential)

# define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = projectExtent(GS_Mean_CV_Map, equalEarthProj)

worldBorderInit = shapefile("WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj))

GS_Max_CV_Projected= projectRaster(GS_Max_CV_Map, to=extentInfo,over=T)
GS_Mean_CV_Projected= projectRaster(GS_Mean_CV_Map, to=extentInfo,over=T)


pal = colorRampPalette(c("#F7F7F7","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))

pdf("Plots/Fig_S2_Present_AGB_Uncertainty_of_GS_Models.pdf",,width =12, height=16)
par(mfrow=c(3,1))

# plot(worldBorder,border="gray40",col="gray40",xlim=c(-180,180),ylim=c(-60,84),axes=T,main="GS Present uncertainty",lwd=0.2)
# plot(presentUncertainMap,col=rev(pal(100)),add=T,breaks=c(seq(0,0.25,0.0025)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="GS Max uncertainty",lwd=0.2)
plot(GS_Max_CV_Projected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray25",col="gray25",axes=T,main="GS Mean uncertainty",lwd=0.2)
plot(GS_Mean_CV_Projected,col=(pal(100)),add=T,breaks=c(seq(0,0.3,0.003)),legend=T,maxpixels = 50000000)

dev.off()

