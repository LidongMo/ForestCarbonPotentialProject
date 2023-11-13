library(raster)
library(stringr)

# don't need to run this code as we are providing the merged maps
# define the vector of map names

mapNames = c("GS_Max","GS_Mean")


for (mp in mapNames)
{
	# generate the names of sub tiers 
	tierList = list.files(path="Data/BiomassMaps", pattern = paste(mp,"_AGB_Modeling_VC_Map*.",sep=""),full.name=T, ignore.case = TRUE)
	print(tierList)
	# # read them into a list of raster
	mergedMap = raster::merge(raster(tierList[[1]]),raster(tierList[[2]]))
	# print("done")
	writeRaster(mergedMap,paste("Data/BiomassMergedMaps/",mp,"_Present_AGB_Modeling_VC_Map_Merged.tif",sep=""),overwrite=T)
	print(mp)
}

