library(raster)
library(stringr)
# don't need to run this code, the mereged maps are alreday providec

# define the vector of map names

mapNames = c("GS1_MaxScaler","GS2_MaxScaler","GS1_MeanScaler","GS2_MeanScaler","SD1","SD2","HM1","HM2","WK1","WK2")

for (mp in mapNames)
{
	# generate the names of sub tiers 
	tierList = list.files(path="Data/BiomassMaps", pattern = paste(mp,"_Potential*.",sep=""),full.name=T, ignore.case = TRUE) 
	print(tierList)
	# # read them into a list of raster
	mergedMap = raster::merge(raster(tierList[[1]]),raster(tierList[[2]]))
	# print("done")
	writeRaster(mergedMap,paste("Data/BiomassMergedMaps/",mp,"_Potential_Model_Variation_Coefficient_Map_Merged.tif",sep=""),overwrite=T)
	print(mp)
}

