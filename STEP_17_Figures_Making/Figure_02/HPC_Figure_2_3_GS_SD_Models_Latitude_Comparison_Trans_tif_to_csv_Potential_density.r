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


mapNames = c("GS_Max1","GS_Max2","GS_Mean1","GS_Mean2","SD1","SD2","HM1","HM2","WK1","WK2")

for (tp in c("Full","Net"))
{
    for (mp in mapNames)
    {
        # # read the image for each model
    perModelMap =  raster(paste("Data/BiomassMergedMaps/",mp,"_",tp,"_TGB_carbon_density_Map_Merged.tif",sep=""))
    # transfer raster into data frame
    # replace the values equals or smaller than 0 to NA
    perModelMap[perModelMap<=0] = NA
    # transfer to data frame
    transferedTable = rasterToPoints(perModelMap)
    # write table to local folder
    write.csv(transferedTable,paste("Data/GS_SD_Compare/",mp,"_",tp,"_TGB_carbon_density_Map_transfer_to_table.csv",sep=""))
    print(paste(mp,"_",tp,"_transfered",sep=""))
    }
}
