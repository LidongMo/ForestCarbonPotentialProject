library(raster)
library(parallel)
library(rgdal)
library(stringr)
library(data.table)
library(gmodels)
library(dplyr)

# aggregate the two pixel wise data frame

mapNames = c("GS_Max1","GS_Max2","GS_Mean1","GS_Mean2","SD1","SD2","HM1","HM2","WK1","WK2")

for (tp in c("Full","Net"))
{
    for (mp in mapNames)
    {
        # read the image for each model
        perModelTable =  fread(paste("Data/GS_SD_Compare/",mp,"_",tp,"_TGB_carbon_density_Map_transfer_to_table.csv",sep=""))[,-1]
        # transfer raster into data frame
        names(perModelTable) = c("x","y","Density")
        if(mp %in%c("GS_Max1","GS_Max2","GS_Mean1","GS_Mean2"))
        {
            # add the model type logo
            perModelTable$Type = "GS"
        }else
        {
            perModelTable$Type = "SD"
        }

        # aggregate the table by each 0.1 degree of latitude
        aggregatedTable = perModelTable %>%
                   mutate(LatRound = round(y,1)) %>%
                   group_by(LatRound,Type) %>%
                   summarise(mean = mean(Density),
                             lowCI = mean(Density) - sd(Density),
                             highCI = mean(Density) + sd(Density)) %>%
                  as.data.frame()

        # write table to local folder
        write.csv(aggregatedTable,paste("Data/GS_SD_Compare/",mp,"_",tp,"_TGB_carbon_density_latitude_agregated_table_with_SD.csv",sep=""))
        print(paste(mp,"_",tp,"_aggregated",sep=""))
    }
}

