# get the 
library(data.table)
library(parallel)
library(dplyr)
library(raster)
# set the working directory by cmd 
setwd("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation/")

#################################################################
## STEP 1 datra checking
#################################################################

coverTable = fread("PlantedForestPotential/Planted_and_Present_Cover_table.csv")[,-1] %>% na.omit()
names(coverTable) = c("x","y","PlantedCover","PresentCover")
coverTable$Ratio = coverTable$PlantedCover/coverTable$PresentCover
coverTableFiltered = coverTable %>% filter(Ratio<2.21441006)

# > quantile(coverTableFiltered$Ratio,c(0.025,0.5,0.975))
#        2.5%         50%       97.5% 
# 0.009945641 0.574408084 1.994094524 
# make a plot
pdf("Data/PlantedForestPotential/Planted_Present_Ratio.pdf",height=10,width=12)
hist(coverTableFiltered$Ratio)
dev.off()


#################################################################
## STEP 2 calculate the plantation potential 
#################################################################
library(data.table)
library(parallel)
library(dplyr)
library(raster)


modelNames = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
sumGetFunc = function(mn)
{
    # read the tables of carbon stock partition
    plantedPotentialModel= read.csv(paste("PlantedForestPotential/OutputTable/",mn,"_Biome_Level_Statistics_Planted_PGB.csv",sep=""))[,-1]   %>% data.table::last() %>% dplyr::select(PlantedPotential)

    return(plantedPotentialModel)
}

outputList = lapply(modelNames,sumGetFunc)
outputTable = rbindlist(outputList)

# > outputTable
#     PlantedPotential
#  1:         8.484635
#  2:         8.898071
#  3:         9.187048
#  4:         8.945669
#  5:         7.255393
#  6:         9.451000
#  7:         4.756956
#  8:         5.263217
#  9:         7.150824
# 10:         8.636897

plantedPotentialSOC = read.csv(paste("PlantedForestPotential/OutputTable/SoilCarbon_Biome_Level_Statistics_Planted_PGB.csv",sep=""))[,-1]   %>% data.table::last() %>% dplyr::select(PlantedPotential)

# full of the plantated potential 

mean(outputTable$PlantedPotential)+plantedPotentialSOC$PlantedPotential
 # the value is 10.5
