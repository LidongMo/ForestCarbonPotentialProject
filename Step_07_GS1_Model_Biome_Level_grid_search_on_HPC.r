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

##################################################################################
##########STEP 1 Grid Search 
##################################################################################

# 1: 'Tropical & Subtropical Moist Broadleaf Forests',
# 2: 'Tropical & Subtropical Dry Broadleaf Forests',
# 3: 'Tropical & Subtropical Coniferous Forests',
# 4: 'Temperate Broadleaf & Mixed Forests',
# 5: 'Temperate Conifer Forests',
# 6: 'Boreal Forests/Taiga',
# 7: 'Tropical & Subtropical Grasslands, Savannas & Shrublands',
# 8: 'Temperate Grasslands, Savannas & Shrublands',
# 9: 'Flooded Grasslands & Savannas',
# 10: 'Montane Grasslands & Shrublands',
# 11: 'Tundra',
# 12: 'Mediterranean Forests, Woodlands & Scrub',
# 13: 'Deserts & Xeric Shrublands',
# 14: 'Mangroves'
# setwd("/Volumes/CrowtherLabRAID/Lidong_Mo/BiomassEstimation/")

# # generate the biomeVector for the biomes which has enough data points

portNumber = 54321 
seedVector = c(500:514)
biomeVector = c(1,2,4,5,6,7,8,10,11,12,13)
for (ss in seedVector)
{
    rawData = fread(paste("CovariatesTable/20201007_Merged_Covariates_sampled_dataset.csv",sep="")) # this file was not in the github folder, but you can generate by your self from last step
    # delete the na values
    rawData = na.omit(rawData)
    # subset the train data by per biome
    for (bm in biomeVector)
    # for (fType in forestTypeVector)
    {
        # get the traind data for each biome
        # typeLevelTrainData = na.omit(rawData[rawData$ForestType == fType,-c('.geo')])
        typeLevelTrainData = na.omit(rawData[rawData$WWF_Biome == bm,])
        
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
        # show the information
        print(bm)
        print(dim(subsetTailsDF))

        # intial the h2o 
        duplicateTable = subsetTailsDF
        # tranform the data frame format lat lon into spatial lat long as spatial points
        coordinates(subsetTailsDF) = ~ x + y
        # # allocate the projection
        proj4string(subsetTailsDF) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
        subsetTailsDF@data = duplicateTable
        # # set the random seed
        set.seed(ss)
        if (bm %in% c(4,5,6,12))
        {
            gridSubsampledPoints = sample.grid(subsetTailsDF, cell.size = c(0.15,0.15), n = 1)
            gridSubsampledTable = gridSubsampledPoints$subset@data
            print(dim(gridSubsampledPoints$subset@data))
        }else
        {
            # just directly use the raw data 
            gridSubsampledTable = duplicateTable
        }
        # print the information of the grid sub sampled table
        print(dim(gridSubsampledTable))
        # duplicate a new table for the shaepfile data allocation
        processTable = gridSubsampledTable
        # # tranform the data frame format lat lon into spatial lat long as spatial points
        coordinates(processTable) = ~ x + y
        # allocate the projection
        proj4string(processTable) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
        processTable@data = gridSubsampledTable
        
        writeOGR(processTable,dsn="ShapeFilesBiome",layer=paste("Biome_",bm,"_Seed_",ss,"_FinalModel",sep=""), driver = "ESRI Shapefile",overwrite=T)
        # gridSubsampledPoints = sample.grid(subsetTailsDF, cell.size = c(0.25,0.25), n = 1)
        # gridSubsampledTable = gridSubsampledPoints$subset@data
        # print(dim(gridSubsampledPoints$subset@data))
        
        localH2O = h2o.init(ip = 'localhost', port = portNumber, nthreads= 24,max_mem_size = '240g') #portNumber[portOrder]
        trainData.hex = as.h2o(gridSubsampledTable, destination_frame = "trainData.hex")
    #  #     finalTable = subsetTails
        trainVaribles = c('Aridity_Index','CHELSA_Annual_Mean_Temperature','CHELSA_Annual_Precipitation','CHELSA_Isothermality','CHELSA_Max_Temperature_of_Warmest_Month','CHELSA_Mean_Diurnal_Range','CHELSA_Mean_Temperature_of_Coldest_Quarter','CHELSA_Mean_Temperature_of_Driest_Quarter','CHELSA_Mean_Temperature_of_Warmest_Quarter','CHELSA_Mean_Temperature_of_Wettest_Quarter','CHELSA_Min_Temperature_of_Coldest_Month','CHELSA_Precipitation_Seasonality','CHELSA_Precipitation_of_Coldest_Quarter','CHELSA_Precipitation_of_Driest_Month','CHELSA_Precipitation_of_Driest_Quarter','CHELSA_Precipitation_of_Warmest_Quarter','CHELSA_Precipitation_of_Wettest_Month','CHELSA_Precipitation_of_Wettest_Quarter','CHELSA_Temperature_Annual_Range','CHELSA_Temperature_Seasonality','Depth_to_Water_Table','EarthEnvTopoMed_Eastness','EarthEnvTopoMed_Elevation','EarthEnvTopoMed_Northness','EarthEnvTopoMed_ProfileCurvature','EarthEnvTopoMed_Roughness','EarthEnvTopoMed_Slope','SG_Absolute_depth_to_bedrock','WorldClim2_SolarRadiation_AnnualMean','WorldClim2_WindSpeed_AnnualMean','EarthEnvCloudCover_MODCF_interannualSD','EarthEnvCloudCover_MODCF_intraannualSD','EarthEnvCloudCover_MODCF_meanannual','EarthEnvTopoMed_AspectCosine','EarthEnvTopoMed_AspectSine','LandCoverClass_Cultivated_and_Managed_Vegetation','Human_Disturbance','LandCoverClass_Urban_Builtup','SG_Clay_Content_0_100cm','SG_Coarse_fragments_0_100cm','SG_Sand_Content_0_100cm','SG_Silt_Content_0_100cm','SG_Soil_pH_H2O_0_100cm','WDPA','cropland','grazing',"pasture","rangeland")
        #  define the grid for coming grid search
        RF_Parameters = list(ntrees = c(100),
                            mtries = c(3,6,9,13),
                            min_rows = c(1,3,5,10,15,20),
                            max_depth = c(10,20,30,40,50,60,70,80,100),
                            sample_rate = c(0.632))
        
        RF_Search_Criteria = list(strategy = "RandomDiscrete", 
                                max_models = 48,
                                seed=ss)

    #  # run the parameter grid
        RF_Grid = h2o.grid("randomForest", x = trainVaribles, y = "BmssDns",
                        grid_id = "RF_Grid",
                        training_frame = trainData.hex,
                        seed = ss,
                        nfolds = 10,
                        hyper_params = RF_Parameters,
                        search_criteria = RF_Search_Criteria)
        RF_GridResult = h2o.getGrid(grid_id = "RF_Grid",
                                    sort_by = "R2",
                                    decreasing = T)
        outputTable = RF_GridResult@summary_table
        # add one column information for 
        outputTable$Biome = bm
        # print(RF_GridResult@summary_table)
        write.csv(outputTable,paste("GridSearchResults_GS1/Parameter_From_Grid_for_Biomass_GridSample_Transformed_Filter_AllYears_Forest_Type_",bm,"_new_Model_Trial3_without_Vegetation_Seed_",ss,"_LOG_Ver_02.csv",sep=""))

        h2o.shutdown(prompt = F)
        portNumber = portNumber+2
    }
}

##################################################################################
##########STEP 2 Grid Search result merge
##################################################################################

# filterNames = c("Full_GFBI_AllYears")
biomeVector = c(1,2,4,5,6,7,8,10,11,12,13)
for (modelNumber in 1:15)
{
    processTable = data.frame()
    for(bm in biomeVector)
    {
        if (bm %in% c(4,5,6,12))
        {
            ss = 499+modelNumber #                                                      
            rawGridSearchTable = fread(paste("GridSearchResults_GS1/Parameter_From_Grid_for_Biomass_GridSample_Transformed_Filter_AllYears_Forest_Type_",bm,"_new_Model_Trial3_without_Vegetation_Seed_",ss,"_LOG_Ver_02.csv",sep=""))[,-1]
            returnRow = rawGridSearchTable[1,]
            
        }else
        {
            ss= 500
            rawGridSearchTable = fread(paste("GridSearchResults_GS1/Parameter_From_Grid_for_Biomass_GridSample_Transformed_Filter_AllYears_Forest_Type_",bm,"_new_Model_Trial3_without_Vegetation_Seed_",ss,"_LOG_Ver_02.csv",sep=""))[,-1]
            returnRow = rawGridSearchTable[modelNumber,]
        }
        returnRow$Biome = bm
        processTable = rbind(processTable,returnRow)

    }
    write.csv(processTable,paste("GridSearchResults_GS1/GridSearchResultsMerged/Parameter_From_Grid_for_Biomass_GridSample_Transformed_Filter_Full_GFBI_AllYears_Forest_Type_Final_Model_Parameters_LOG_Tranformed_Ver_02_Model_Number_",modelNumber,".csv",sep=""))
}

