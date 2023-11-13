# Different data filtering comparison by grid seach 
library(data.table)
library(parallel)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(GSIF)
library(rgdal)


# load the data with upper and lower bound information of cover
rawData1 = fread(paste("Data/ForestCoverScaler/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution.csv",sep=""))[,-1] 

# correction of the far east russian points
farEastRussiaPlots = rawData1 %>% dplyr::filter(x>136&y>58) %>% dplyr::mutate(x= x*(-1))#229

otherPlots = rawData1 %>% dplyr::filter(x<=136|y<=58) #563852
# rbind the modified tables together
rawData = rbind(farEastRussiaPlots,otherPlots)

# add the carbon density row
rawData$CarbonDensity = rawData$BmssDns/1000*rawData$CarbonConcentration #divide 1000 for the kg to t transformation
# subset the train data by per biome
biomeVector = c(1,2,4,5,6,7,8,9,10,11,12,13)

for (bm in 1:14)
{
    if (bm %in% biomeVector)
    {
        # get the traind data for each biome
        # typeLevelTrainData = na.omit(rawData[rawData$ForestType == fType,-c('.geo')])
        typeLevelTrainData = na.omit(rawData[rawData$WWF_Biome == bm,])

        # # typeLevelTrainData = rawData
        # # delete the biomass values larger than 1790
        typeLevelTrainData = typeLevelTrainData[typeLevelTrainData$CarbonDensity <1867,]
        typeLevelTrainData$LogCarbonDensity = log(typeLevelTrainData$CarbonDensity+1)
        # # # # quantileInfo = quantile(biomeLevelTrainData$BimssAr, c(0.025,0.975),na.rm=T) 
        MADValue = mad(typeLevelTrainData$LogCarbonDensity)
        # # print(MADValue)
        largeRange = 2.5*MADValue + median(typeLevelTrainData$LogCarbonDensity)
        smallRange =  median(typeLevelTrainData$LogCarbonDensity)-2.5*MADValue
        # # # subset the data frame by the quantile range
        subsetTailsDF_Part1 = typeLevelTrainData[typeLevelTrainData$LogCarbonDensity<=largeRange&typeLevelTrainData$LogCarbonDensity>=smallRange,]

        # # # subset the data frame by the quantile range
        subsetTailsDF_Part2 = typeLevelTrainData[typeLevelTrainData$LogCarbonDensity<smallRange,] %>% dplyr::filter((Human_Disturbance>0.1&PresentTreeCover<0.1))
        
        print(bm)
        subsetTailsDF = rbind(subsetTailsDF_Part1,subsetTailsDF_Part2)

        print(1-nrow(subsetTailsDF)/nrow(typeLevelTrainData))
    }else
    {
        typeLevelTrainData = na.omit(rawData[rawData$WWF_Biome == bm,])
        typeLevelTrainData = typeLevelTrainData[typeLevelTrainData$CarbonDensity <1867,]
        typeLevelTrainData$LogCarbonDensity = log(typeLevelTrainData$CarbonDensity+1)
        print(bm)
        subsetTailsDF = typeLevelTrainData
        print((1-nrow(subsetTailsDF)/nrow(typeLevelTrainData)))
    }
    

    write.csv(subsetTailsDF,paste("Data/GroundSourcedModel/CleanedSubTables/Biome_",bm,"_Sub_Sample_Table.csv",sep=""))
}
# biome number and the cleaned ratio of the data
# [1] 1
# [1] 0.02777778
# [1] 2
# [1] 0.05177112
# [1] 3
# [1] 0
# [1] 4
# [1] 0.05012255
# [1] 5
# [1] 0.02964339
# [1] 6
# [1] 0.03354969
# [1] 7
# [1] 0.0384216
# [1] 8
# [1] 0.01958079
# [1] 9
# [1] 0.04347826
# [1] 10
# [1] 0.05839416
# [1] 11
# [1] 0.0008654262
# [1] 12
# [1] 0.009319333
# [1] 13
# [1] 0.03097673
# [1] 14
# [1] 0



# this part is for figure S13
rawData = fread(paste("Data/ForestCoverScaler/20230126_Merged_Table_with_Forest_Cover_Scaler_with_1km_resolution.csv",sep=""))[,-1] 

# adaptation of the far east russian points
farEastRussiaPlots = rawData1 %>% dplyr::filter(x>136&y>58) %>% dplyr::mutate(x= x*(-1))#229

otherPlots = rawData1 %>% dplyr::filter(x<=136|y<=58) #563852
# rbind the modified tables together
rawData = rbind(farEastRussiaPlots,otherPlots)

# add the carbon density row
rawData$CarbonDensity = rawData$BmssDns/1000*rawData$CarbonConcentration #divide 1000 for the kg to t transformation
# define the forest type vector
biomeVector = 1:14
biomeNames = c("Tropical moist","Tropical dry","Temperate broadleaf","Temperate coniferous","Boreal","Tropical savanna","Temperate savanna","Flooded savanna","Montane grassland","Tundra","Mediterranean forest","Desert")
pdf(paste("Plots/Figure_S13_HistGram_of_carbon_density_for_Biome_Level_LOG_Transed_MAD_25_Filtered_20230428.pdf",sep=""))
par(mfrow=c(4,3),oma=c(0,0,3,0))

for (bm in biomeVector)
{
    if (bm %in% c(1,2,4,5,6,7,8,9,10,11,12,13))
    {
        # get the traind data for each biome
        # typeLevelTrainData = na.omit(rawData[rawData$ForestType == fType,-c('.geo')])
        typeLevelTrainData = na.omit(rawData[rawData$WWF_Biome == bm,])

        # # typeLevelTrainData = rawData
        # # delete the biomass values larger than 1867
        typeLevelTrainData = typeLevelTrainData[typeLevelTrainData$CarbonDensity <1867,]
        typeLevelTrainData$LogCarbonDensity = log(typeLevelTrainData$CarbonDensity+1)
        # # # # quantileInfo = quantile(biomeLevelTrainData$BimssAr, c(0.025,0.975),na.rm=T) 
        MADValue = mad(typeLevelTrainData$LogCarbonDensity)
        # # print(MADValue)
        largeRange = 2.5*MADValue + median(typeLevelTrainData$LogCarbonDensity)
        smallRange =  median(typeLevelTrainData$LogCarbonDensity)-2.5*MADValue
        # # # subset the data frame by the quantile range
        subsetTailsDF_Part1 = typeLevelTrainData[typeLevelTrainData$LogCarbonDensity<=largeRange&typeLevelTrainData$LogCarbonDensity>=smallRange,]

        # # # subset the data frame by the quantile range
        subsetTailsDF_Part2 = typeLevelTrainData[typeLevelTrainData$LogCarbonDensity<smallRange,] %>% dplyr::filter((Human_Disturbance>0.1&PresentTreeCover<0.1))
        
        print(bm)
        subsetTailsDF = rbind(subsetTailsDF_Part1,subsetTailsDF_Part2)
        # # make a histgram
        hist(subsetTailsDF$LogCarbonDensity,main=paste(biomeNames[which(c(1,2,4,5,6,7,8,9,10,11,12,13)==bm)],sep=""),xlab="Carbon density (log transformed)",ylab="Frequency")
        abline(v = largeRange, col="red", lwd=1, lty=2)
        abline(v = smallRange, col="red", lwd=1, lty=2)
        print(paste("Biome: ",biomeNames[which(c(1,2,4,5,6,7,8,10,11,12,13)==bm)],sep=""))


    }else
    {
        typeLevelTrainData = na.omit(rawData[rawData$WWF_Biome == bm,])
        typeLevelTrainData = typeLevelTrainData[typeLevelTrainData$CarbonDensity <1867,]
        typeLevelTrainData$LogCarbonDensity = log(typeLevelTrainData$CarbonDensity+1)
        print(bm)
        subsetTailsDF = typeLevelTrainData
        print(1-nrow(subsetTailsDF)/nrow(typeLevelTrainData))

    }
    
}

mtext("Log tranformed (filtered 2.5 MAD of Carbon density(log))", line=0, side=3, outer=TRUE, cex=1.5)
dev.off()

################################
# # merge the cleaned biome level tables into one and make a shapefile
################################

tableNameList = list.files(path = "Data/GroundSourcedModel/CleanedSubTables",pattern=".csv",full.name=T)
# load each table into a a list which has biome level table as an element
tableList = lapply(tableNameList,fread,na.strings="None")
# kick out those rows with NA's
rbindedTable = rbindlist(tableList)

# write to local folder 
write.csv(rbindedTable,"Data/GroundSourcedModel/20230126_Merged_Covariates_sampled_dataset_outliers_cleaned.csv")

# ################################
# # transfer this into shape files
# ################################

cleanedTableWithCovariates = fread("Data/GroundSourcedModel/20230126_Merged_Covariates_sampled_dataset_outliers_cleaned.csv") %>% dplyr::select(x,y,WWF_Biome,CarbonDensity,PlotAre)
# make a replicate
processTable = cleanedTableWithCovariates
# # tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(processTable) = ~ x + y
# allocate the projection
proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
processTable@data = cleanedTableWithCovariates 
# write to localfolder
writeOGR(processTable,dsn="Data/GroundSourcedModel/ShapeFilesOutliersCleaned",layer= "Full_GFBI_AllYears_with_outliers_Cleaned_0427", driver = "ESRI Shapefile",overwrite=T)








