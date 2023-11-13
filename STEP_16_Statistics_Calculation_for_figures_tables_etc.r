# Biome level, land use type level uncertainty for different parts of carbon.
library(data.table)
library(dplyr)

##STEP 2 Total living biomass

modelNames1 = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
sumGetFunc = function(mn)
{
    # read the tables of carbon stock partition
    landuseTypePotentialTB = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Biome_Level_Statistics.csv",sep=""))[,-1] %>% 
    mutate(ConservationPotential = PresentPotential- Present) %>%data.table::last() 

    return(landuseTypePotentialTB)
}

outputList1 = lapply(modelNames1,sumGetFunc)
outputTable1 = rbindlist(outputList1)

TGB_meanVec = outputTable1 %>% apply(MARGIN =2, FUN=mean)  %>%t() %>% as.data.frame()
TGB_minVec = outputTable1 %>%  apply(MARGIN =2, FUN=min)%>%t() %>% as.data.frame()
TGB_maxVec = outputTable1 %>% apply(MARGIN =2, FUN=max)%>%t() %>% as.data.frame()


modelNames2 = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler")

outputList2 = lapply(modelNames2,sumGetFunc)
outputTable2 = rbindlist(outputList2)

TGB_meanVec_GS = outputTable2 %>% apply(MARGIN =2, FUN=mean)%>%t() %>% as.data.frame()
TGB_minVec_GS = outputTable2 %>%  apply(MARGIN =2, FUN=min)%>%t() %>% as.data.frame()
TGB_maxVec_GS = outputTable2 %>% apply(MARGIN =2, FUN=max) %>%t() %>% as.data.frame()

modelNames3 = c("SD1","SD2","WK1","WK2","HM1","HM2")

outputList3 = lapply(modelNames3,sumGetFunc)
outputTable3 = rbindlist(outputList3)

TGB_meanVec_SD = outputTable3 %>% apply(MARGIN =2, FUN=mean)  %>%t() %>% as.data.frame()
TGB_minVec_SD = outputTable3 %>%  apply(MARGIN =2, FUN=min) %>%t() %>% as.data.frame()
TGB_maxVec_SD = outputTable3 %>% apply(MARGIN =2, FUN=max)%>%t() %>% as.data.frame()


TGB_Vec = paste(TGB_meanVec,"(",TGB_minVec,"-",TGB_maxVec,")",sep="")

finalTable = rbind(outputTable1,TGB_meanVec,TGB_minVec,TGB_maxVec,TGB_meanVec_GS,TGB_minVec_GS,TGB_maxVec_GS,TGB_meanVec_SD,TGB_minVec_SD,TGB_maxVec_SD)  %>% mutate(Names = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2","Mean","Min","Max","Mean_GS","Min_GS","Max_GS","Mean_SD","Min_SD","Max_SD"))

write.csv(finalTable,"Data/BiomeLevelStatistics/StatisticsForModels/Summery_of_TGB_potential_for_different_landuse.csv")
            

##STEP3 # plant biomass
modelNames1 = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
sumGetFunc = function(mn)
{
    # read the tables of carbon stock partition
    landuseTypePotentialPGB = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Biome_Level_Statistics_with_litter.csv",sep=""))[,-1] %>% 
    mutate(ConservationPotential = PresentPotential- Present) %>%data.table::last() 

    return(landuseTypePotentialPGB)
}

outputList1 = lapply(modelNames1,sumGetFunc)
outputTable1 = rbindlist(outputList1)

TGB_meanVec = outputTable1 %>% apply(MARGIN =2, FUN=mean) %>%t() %>% as.data.frame()
TGB_minVec = outputTable1 %>%  apply(MARGIN =2, FUN=min)%>%t() %>% as.data.frame()
TGB_maxVec = outputTable1 %>% apply(MARGIN =2, FUN=max) %>%t() %>% as.data.frame()


modelNames2 = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler")

outputList2 = lapply(modelNames2,sumGetFunc)
outputTable2 = rbindlist(outputList2)

TGB_meanVec_GS = outputTable2 %>% apply(MARGIN =2, FUN=mean) %>%t() %>% as.data.frame()
TGB_minVec_GS = outputTable2 %>%  apply(MARGIN =2, FUN=min) %>%t() %>% as.data.frame()
TGB_maxVec_GS = outputTable2 %>% apply(MARGIN =2, FUN=max) %>%t() %>% as.data.frame()

modelNames3 = c("SD1","SD2","WK1","WK2","HM1","HM2")

outputList3 = lapply(modelNames3,sumGetFunc)
outputTable3 = rbindlist(outputList3)

TGB_meanVec_SD = outputTable3 %>% apply(MARGIN =2, FUN=mean) %>%t() %>% as.data.frame()
TGB_minVec_SD = outputTable3 %>%  apply(MARGIN =2, FUN=min)%>%t() %>% as.data.frame()
TGB_maxVec_SD = outputTable3 %>% apply(MARGIN =2, FUN=max) %>%t() %>% as.data.frame()


TGB_Vec = paste(TGB_meanVec,"(",TGB_minVec,"-",TGB_maxVec,")",sep="")

finalTable = rbind(outputTable1,TGB_meanVec,TGB_minVec,TGB_maxVec,TGB_meanVec_GS,TGB_minVec_GS,TGB_maxVec_GS,TGB_meanVec_SD,TGB_minVec_SD,TGB_maxVec_SD)  %>% mutate(Names = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2","Mean","Min","Max","Mean_GS","Min_GS","Max_GS","Mean_SD","Min_SD","Max_SD"))

write.csv(finalTable,"Data/BiomeLevelStatistics/StatisticsForModels/Summery_of_PGB_potential_for_different_landuse.csv")


############################################################################################
# STEP 4 TGB upper and lower boundary this is for the Figure 4
# Biome level, land use type level uncertainty for different parts of carbon.
library(data.table)
library(dplyr)

# read csv of the upper and lower boundary for each model

modelNames = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
rangeGetFunc = function(mn)
{
    # read the tables of carbon stock partition
    landuseTypePotentialRange = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Uncertainty_for_diff_parts_at_Biome_Level.csv",sep=""))[,-1] %>% tail(1) %>% mutate(potTGB_Lower = potAGB_Lower+potRoot_Lower,potTGB_Upper = potAGB_Upper+potRoot_Upper) %>% dplyr::select(potTGB_Lower,potTGB_Upper) 

    landuseTypePotentialMean = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Biome_Level_Statistics.csv",sep=""))[,-1] %>% tail(1) %>% dplyr::select(Present,AbsolutePotential) 
    if (mn %in% c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler"))
    {
        ModelType = "Ground-sourced model"
    }else
    {
         ModelType = "Satellite-derived model"
    }
    # rbind the mean and the range with min and max
    bindedTable = data.frame(landuseTypePotentialMean,landuseTypePotentialRange,ModelName =mn,Model =ModelType)
    # retrun the table 
    return(bindedTable)
}

outputList4 = lapply(modelNames,rangeGetFunc)
outputTable4 = rbindlist(outputList4)

write.csv(outputTable4,"Data/BiomeLevelStatistics/TGB_Statistics_summary_with_mean_and_range.csv")




# STEP 4 Table S2 preperation
library(data.table)
library(dplyr)


modelNames1 = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
sumGetFunc = function(mn)
{
    # read the tables of carbon stock partition
    landuseTypePotentialTB = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Biome_Level_Statistics.csv",sep=""))[,-1] %>% 
    mutate(ConservationPotential = PresentPotential- Present) # subset the table 

    return(landuseTypePotentialTB)
}

outputList = lapply(modelNames1,sumGetFunc)

allMatrix = abind::abind(outputList, along=3)
# oget the mean for each element
outputMean = apply(allMatrix, c(1,2), mean)
outputMin = apply(allMatrix, c(1,2), min)
outputMax = apply(allMatrix, c(1,2), max)

outputS2Table = rbind(outputMean,outputMin,outputMax)

write.csv(outputS2Table,"Data/BiomeLevelStatistics/TGB_Statistics_summary_Biome_and_land_use_type_for_Table_S2.csv")



# STEP 5 Table S3 preperation
library(data.table)
library(dplyr)


modelNames = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
forestTyowSumGetFunc = function(mn)
{
    # read the tables of carbon stock partition
    landuseTypePotentialTB = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Biome_Level_Statistics.csv",sep=""))[,-1] %>% 
    mutate(ConservationPotential = PresentPotential- Present) %>% mutate(ForestType = c("Tropical","Tropical","Tropical","Temperate","Temperate","Boreal","Tropical","Temperate","Tropical","Temperate","Boreal","Dryland","Dryland","Tropical","SUM"))
    aggregatedTable = aggregate(landuseTypePotentialTB %>% dplyr::select(-ForestType),by = list(landuseTypePotentialTB$ForestType),FUN = sum)  %>% arrange(factor(Group.1, levels = c("Tropical","Temperate","Boreal","Dryland","SUM"))) %>% dplyr::select(-Group.1)
    return(aggregatedTable)
}

outputList = lapply(modelNames,forestTyowSumGetFunc)

allMatrix = abind::abind(outputList, along=3)
# oget the mean for each element
outputMean = apply(allMatrix, c(1,2), mean) %>% 'row.names<-'(c("Tropical","Temperate","Boreal","Dryland","SUM"))

outputMin = apply(allMatrix, c(1,2), min) %>% 'row.names<-'(c("Tropical","Temperate","Boreal","Dryland","SUM"))

outputMax = apply(allMatrix, c(1,2), max) %>% 'row.names<-'(c("Tropical","Temperate","Boreal","Dryland","SUM"))

outputS3Table = rbind(outputMean,outputMin,outputMax)

write.csv(outputS3Table,"Data/BiomeLevelStatistics/TGB_Statistics_summary_ForestType_and_land_use_type_for_Table_S3.csv")


