library(raster)
library(RColorBrewer)
library(rgeos)
library(viridis)
library(rgdal)
library(sp)
# set working directory
# load the two difference maps

netDifferenceMap01 = raster("Data/BiomassMaps/Net_TGB_potential_diff_Map_of_GS_SD0000000000-0000023296.tif")
netDifferenceMap02 = raster("Data/BiomassMaps/Net_TGB_potential_diff_Map_of_GS_SD0000000000-0000000000.tif")


# merge the tiles of all three maps 

netDifferenceMap = merge(netDifferenceMap01,netDifferenceMap02)

# write them into local folder
writeRaster(netDifferenceMap,"Data/BiomassMergedMaps/Difference_in_biomass_density_of_GS_and_SD_Net.tif",overwrite=T)

