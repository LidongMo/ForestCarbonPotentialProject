# load the library
library(raster)
library(RColorBrewer)
library(rgeos)
library(sp)
library(rgdal)
# load the TGb carbon density maps for present predictions

GS_Max_TGB_Density = raster("Data/BiomassMergedMaps/GS_Max_Present_TGB_Density_Map_Merged.tif")
GS_Mean_TGB_Density = raster("Data/BiomassMergedMaps/GS_Mean_Present_TGB_Density_Map_Merged.tif")

HM_TGB_Density = raster("Data/BiomassMergedMaps/HM_Present_TGB_Density_Map_Merged.tif")
SD_TGB_Density = raster("Data/BiomassMergedMaps/SD_Present_TGB_Density_Map_Merged.tif")
WK_TGB_Density = raster("Data/BiomassMergedMaps/WK_Present_TGB_Density_Map_Merged.tif")


# load the forest cover map
forestCover = raster("Data/BiomassMergedMaps/Hansen_2010_Present_Forest_Cover_Map_Merged.tif")
forestCover[forestCover>0] =1
forestCover[forestCover<=0] =NA
# mask the region inside the defined forest cover 
GS_Max_TGB_Density_Masked = mask(GS_Max_TGB_Density,forestCover)
GS_Mean_TGB_Density_Masked = mask(GS_Mean_TGB_Density,forestCover)
HM_TGB_Density_Masked = mask(HM_TGB_Density,forestCover)
SD_TGB_Density_Masked = mask(SD_TGB_Density,forestCover)
WK_TGB_Density_Masked = mask(WK_TGB_Density,forestCover)
# adjust the pixels that have density higher than 250t C/ha
GS_Max_TGB_Density_Masked[GS_Max_TGB_Density_Masked >=250] = 250
GS_Mean_TGB_Density_Masked[GS_Mean_TGB_Density_Masked >=250] = 250
HM_TGB_Density_Masked[HM_TGB_Density_Masked >=250] = 250
SD_TGB_Density_Masked[SD_TGB_Density_Masked >=250] = 250
WK_TGB_Density_Masked[WK_TGB_Density_Masked >=250] = 250
 # define the equal earth projection
equalEarthProj = "+proj=eqearth +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
# project raster to equal earth
extentInfo = projectExtent(GS_Mean_TGB_Density_Masked, equalEarthProj)

worldBorderInit = shapefile("WORLD_BORDERS/TM_WORLD_BORDERS_SIMPL-0.3.shp")
worldBorderRaw = crop(worldBorderInit,c(-180, 180, -60, 84))
# set the resolution of how the polygon plot process
worldBorder = gSimplify(worldBorderRaw, tol = 0.01, topologyPreserve = TRUE)
# reproject the world border
worldEqualEarth = spTransform(worldBorder,to=extentInfo, CRS(equalEarthProj))

GS_Max_TGB_DensityProjected= projectRaster(GS_Max_TGB_Density_Masked, to=extentInfo,over=T)
GS_Mean_TGB_DensityProjected = projectRaster(GS_Mean_TGB_Density_Masked, to=extentInfo,over=T)

HM_TGB_DensityProjected= projectRaster(HM_TGB_Density_Masked, to=extentInfo,over=T)
SD_TGB_DensityProjected = projectRaster(SD_TGB_Density_Masked, to=extentInfo,over=T)
WK_TGB_DensityProjected= projectRaster(WK_TGB_Density_Masked, to=extentInfo,over=T)


# define the colour pallete
pal= colorRampPalette((c("gray80","#A1D99B" ,"#74C476", "#41AB5D" ,"#238B45", "#006D2C", "#00441B","#004521")))

pdf("Plots/Fig_S5_Present_TGB_Biomass_Density_of_GS_and_SD_models.pdf",width =16, height=30)
par(mfrow=c(5,1))

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="GS Max TGB Density",lwd=0.2)
plot(GS_Max_TGB_DensityProjected,col=(pal(50)),add=T,breaks=c(seq(0,250,5)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="GS mean TGB Density",lwd=0.2)
plot(GS_Mean_TGB_DensityProjected,col=(pal(50)),add=T,breaks=c(seq(0,250,5)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Harmionized TGB Density",lwd=0.2)
plot(HM_TGB_DensityProjected,col=(pal(50)),add=T,breaks=c(seq(0,250,5)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="ESA-CCI TGB Density",lwd=0.2)
plot(SD_TGB_DensityProjected,col=(pal(50)),add=T,breaks=c(seq(0,250,5)),legend=T,maxpixels = 50000000)

plot(worldEqualEarth,border="gray80",col="gray80",axes=T,main="Walker TGB Density",lwd=0.2)
plot(WK_TGB_DensityProjected,col=(pal(50)),add=T,breaks=c(seq(0,250,5)),legend=T,maxpixels = 50000000)


dev.off()

