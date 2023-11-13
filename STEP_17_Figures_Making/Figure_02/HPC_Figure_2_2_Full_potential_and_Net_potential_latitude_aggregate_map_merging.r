library(raster)
library(stringr)

# don't need to set working directory as we did this on HPC

# define the vector of map names
mapNames = c("GS_Max1","GS_Max2","GS_Mean1","GS_Mean2","SD1","SD2","HM1","HM2","WK1","WK2")

for (tp in c("Full","Net"))
{
	for (mp in mapNames)
	{
		# generate the names of sub tiers 
		tierList = list.files(path="Data/BiomassMaps", pattern = paste(mp,"*.",tp,"*.",sep=""),full.name=T)
		print(tierList)
		# read them into a list of raster
		mergedMap = raster::merge(raster(tierList[[1]]),raster(tierList[[2]]))
		print("done")
		writeRaster(mergedMap,paste("Data/BiomassMergedMaps/",mp,"_",tp,"_TGB_carbon_density_Map_Merged.tif",sep=""),overwrite=T)
		print(mp)
	}
}

