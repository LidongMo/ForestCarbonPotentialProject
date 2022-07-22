# Different data filtering comparison by grid seach 
library(data.table)
library(parallel)
library(ggplot2)
library(raster)
library(purrr)
library(dplyr)
library(h2o)
library(GSIF)
library(rgdal)

setwd("~/Github_Folder")
biomeVector = c(1,2,4,5,6,7,8,10,11,12,13)


rawData = fread(paste("CovariatesTable/20201007_Merged_Covariates_sampled_dataset.csv",sep=""))[,-1]
# delete the na values
rawData = na.omit(rawData)
# subset the train data by per biome
for (bm in biomeVector)
# for (fType in forestTypeVector)
{
    # get the traind data for each biome
    # typeLevelTrainData = na.omit(rawData[rawData$ForestType == fType,-c('.geo')])
    typeLevelTrainData = na.omit(rawData[rawData$WWF_Biome == bm,])
    # detelete the plots from the regions with small tree cover
    # typeLevelTrainData = typeLevelTrainData[typeLevelTrainData$treecover2000 >=10,]
    # typeLevelTrainData = na.omit(rawData[rawData$WWF_Biome == bm,-c('.geo')])
    # typeLevelTrainData = rawData[rawData$WDPA == 1,]
    # typeLevelTrainData = fread(paste("TreeDensityFIltering/SubTablesCovariates/20200217_Present_Forest_sampledType_",fType,"_seed_500_Covariates_Extraction_Nonfil2.csv",sep=""))
    # typeLevelTrainData = fread(paste("TreeDensityFIltering/SubTablesCovariates/20200523_Present_Forest_sampledType_",fType,"_seed_500_Covariates_Extraction.csv",sep=""))
    
    # subsetTailsDF = typeLevelTrainData
    # subsetTailsDF = subsetTailsDF[subsetTailsDF$treecover2000 >=10,]

    # # typeLevelTrainData = rawData
    # # delete the biomass values larger than 1790
    typeLevelTrainData = typeLevelTrainData[typeLevelTrainData$BmssDns <1867,]
    typeLevelTrainData$BmssDns = log(typeLevelTrainData$BmssDns+1)
    # # # # quantileInfo = quantile(biomeLevelTrainData$BimssAr, c(0.025,0.975),na.rm=T) 
    MADValue = mad(typeLevelTrainData$BmssDns)
    # # print(MADValue)
    largeRange = 2.5*MADValue + median(typeLevelTrainData$BmssDns)
    smallRange = median(typeLevelTrainData$BmssDns) -2.5*MADValue 
    # # # subset the data frame by the quantile range
    subsetTailsDF = typeLevelTrainData[typeLevelTrainData$BmssDns<=largeRange&typeLevelTrainData$BmssDns>=smallRange,]
    #  subsetTailsDF = typeLevelTrainData[typeLevelTrainData$BmssDns<=largeRange,]
    # # subsetTailsDF = typeLevelTrainData
    subsetTailsDF$BmssDns = exp(subsetTailsDF$BmssDns)-1
    print(bm)

    print(dim(subsetTailsDF))

    write.csv(subsetTailsDF,paste("CovariatesTable/CleanedSubTables/Biome_",bm,"_Sub_Sample_Table.csv",sep=""))
}


################################
# # merge the cleaned biome level tables into one and make a shapefile
################################

tableNameList = list.files(path = "CovariatesTable/CleanedSubTables",pattern=".csv",full.name=T)
# load each table into a a list which has biome level table as an element
tableList = lapply(tableNameList,fread,na.strings="None")
# kick out those rows with NA's
rbindedTable = rbindlist(tableList)

# write to local folder 
write.csv(rbindedTable,"CovariatesTable/20211221_Merged_Covariates_sampled_dataset_outliers_cleaned.csv")
################################
# transfer this into shape files
################################

cleanedTableWithCovariates = fread("CovariatesTable/20211221_Merged_Covariates_sampled_dataset_outliers_cleaned.csv") %>% dplyr::select(x,y,BmssDns,treecover2010,WDPA,Human_Disturbance) %>% filter(Human_Disturbance<=0.05|WDPA==1)
cleanedTableWithCovariates$lgBD = log(cleanedTableWithCovariates$BmssDns+1)
# make a replicate
processTable = cleanedTableWithCovariates
# # tranform the data frame format lat lon into spatial lat long as spatial points
coordinates(processTable) = ~ x + y
# allocate the projection
proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
processTable@data = cleanedTableWithCovariates
for (seeds in 0:99)
{
    set.seed(seeds)
    gridSubsampledPoints = sample.grid(processTable, cell.size = c(0.5,0.5), n =1)$subset

    writeOGR(gridSubsampledPoints,dsn="CovariatesTable/GridSampledShapeFiles",layer= paste("GFBI_AllYears_Natural_Gridsampled_",seeds,sep=""), driver = "ESRI Shapefile",overwrite=T)
}