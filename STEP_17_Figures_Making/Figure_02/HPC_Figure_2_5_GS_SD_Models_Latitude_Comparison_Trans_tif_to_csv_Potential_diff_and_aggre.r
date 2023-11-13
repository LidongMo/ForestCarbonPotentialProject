library(raster)
library(parallel)
library(rgdal)
library(stringr)
library(raster)
library(RColorBrewer)
library(rgeos)
library(viridis)
library(data.table)
library(dplyr)
library(stats)
# load the full and net potential 

# read the image for each model
netDifferenceMap = raster("Data/BiomassMergedMaps/Difference_in_biomass_density_of_GS_and_SD_Net.tif")
# transfer raster into data frame

# transfer to data frame
transferedTable = rasterToPoints(netDifferenceMap,na.rm=T)

# write table to local folder
write.csv(transferedTable,"Data/GS_SD_Compare/Net_diff_potential_carbon_density_tif_transfer_to_table.csv")

# read the image for each model
perModelTable =  fread("Data/GS_SD_Compare/Net_diff_potential_carbon_density_tif_transfer_to_table.csv")[,-1] %>% na.omit()
names(perModelTable) = c("x","y","Density")
# aggregate the table by each 0.1 degree of latitude
aggregatedTable = perModelTable %>% mutate(Type = "Net") %>%
               mutate(LatRound = round(y,1)) %>%
               group_by(LatRound) %>%
               summarise(mean = mean(Density),
                         lowCI = mean(Density) - sd(Density),
                         highCI = mean(Density) + sd(Density)) %>%as.data.frame()

# write table to local folder
write.csv(aggregatedTable,paste("Data/GS_SD_Compare/Net_diff_potential_carbon_density_latitude_agregated_table_with_SD.csv",sep=""))
print(paste("_transfered",sep=""))
