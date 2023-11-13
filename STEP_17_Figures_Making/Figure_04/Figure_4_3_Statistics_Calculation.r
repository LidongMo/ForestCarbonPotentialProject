# Biome level, land use type level uncertainty for different parts of carbon.
# set working directory

library(data.table)
library(dplyr)
library(abind)
library(stringr)


# calculate the lower and upper boundary of the uncertainty from boot strapping modeling
modelNames = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
uncerFunc = function(mn)
{
    # load the uncertainty table for each model
    perModelUncertainty = fread(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Uncertainty_for_diff_parts_at_Biome_Level.csv",sep=""))[,-1]
    # calculate the  TGB and PGB upper and lower boundary
    # read the tables of carbon stock partition
    pgbFullPot = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Biome_Level_Statistics_with_litter.csv",sep=""))[,-1] %>% tail(1) %>% dplyr::select(AbsolutePotential)%>% as.numeric()
    soilCarbonPot = fread("Data/BiomeLevelStatistics/StatisticsForModels/SoilCarbon_Landuse_type_Uncertainty.csv")[,-1] %>% 
                      dplyr::select(AbsolutePotential) %>% tail(1) %>% as.numeric()

    fullPotentialVal = pgbFullPot +soilCarbonPot
    # get the sum of each model
    lastRow = perModelUncertainty %>% tail(1) %>% mutate(ModelName = mn,potTotal_Mean = fullPotentialVal) %>% dplyr::select(ModelName, everything())
    return(lastRow)
}

allModelUncertainTable = lapply(modelNames,uncerFunc) %>% rbindlist()
soilUncertaity = fread("Data/BiomeLevelStatistics/StatisticsForModels/SoilCarbon_Uncertainty_for_diff_parts_at_Biome_Level.csv") %>% dplyr::select(SoilCarbonLower,SoilCarbonUpper) %>% tail(1) 
# add the soil uncertainty to the all model uncertaity table
allModelUncertainTable = allModelUncertainTable %>% mutate(model=c("GS","GS","GS","GS","SD","SD","SD","SD","SD","SD"),
                                                           type = c("M1","M1","M2","M2","M1","M2","M1","M2","M1","M2"),
                                                           subtype = c("mean","max","mean","max","ESA-CCI","ESA-CCI","Walker et al.","Walker et al.","Harmonized","Harmonized"),
                                                           potSoil_Lower = soilUncertaity$SoilCarbonLower,
                                                           potSoil_Upper = soilUncertaity$SoilCarbonUpper) %>%relocate(model, type,subtype, .before = preAGB_Lower) 
# fineResultTable = rbind(allModelUncertainTable,lowerRow,upperRow)
write.csv(allModelUncertainTable,"Data/SupplementaryFigures/Figure_04/Bootstapped_uncertainty_all_models.csv")










