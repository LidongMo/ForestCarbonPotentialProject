library(raster)
library(stringr)
# set working directory 

# don't need to set working directory as we did this on HPC

# define the vector of map names

mapNames = c("GS_Max","GS_Mean","SD","HM","WK","GB","IPCC")


for (mp in mapNames)
{
	# generate the names of sub tiers 
	tierList = list.files(path="Data/BiomassMaps", pattern = paste(mp,"_TGB_Density*.",sep=""),full.name=T, ignore.case = TRUE)
	print(tierList)
	# # read them into a list of raster
	mergedMap = raster::merge(raster(tierList[[1]]),raster(tierList[[2]]))
	# print("done")
	writeRaster(mergedMap,paste("Data/BiomassMergedMaps/",mp,"_Present_TGB_Density_Map_Merged.tif",sep=""),overwrite=T)
	print(mp)
}

