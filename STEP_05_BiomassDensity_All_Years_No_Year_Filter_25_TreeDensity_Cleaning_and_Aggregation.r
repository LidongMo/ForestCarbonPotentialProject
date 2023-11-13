# get the 
library(data.table)
library(parallel)
library(dplyr)
library(raster)
library(rgdal)
# set the working directory
#################################################################
## STEP 1 plot siez based cleanning
#################################################################

gfbiPlotsTable = fread("Data/BiomassCalculationData/GFBI_Plots_with_All_Years_Time_Series_Records_20190412_Lidong.csv")[,-1] #this is only the subset of the data, if you need the full data, please contact me directly

# check the plot numbers which has size smaller than 0.01 ha or larger than 1 ha
kickedPlots = gfbiPlotsTable[gfbiPlotsTable$plotArea >=0.01 &gfbiPlotsTable$plotArea <=1,] #|gfbiPlotsTable$plotArea >1 gfbiPlotsTable$plotArea <0.01 &
length(unique(kickedPlots$PLT)) #[1] 114936 plots

length(unique(gfbiPlotsTable$PLT)) #[1] 1188771 

# because some plot names with some plot size lager than 0.01 then the left plot numbers are larger than the 1188771  -114936 = 1073835

# delete the plots with small plot size for invesitagetion 
# the standard is biger than 0.01ha = 100m2
currentYearPlotsTable = gfbiPlotsTable[gfbiPlotsTable$plotArea >=0.01& gfbiPlotsTable$plotArea <=1,] # 1073835 plots left with 1115129 recordings(some plots has timeseries observation)
#  in this plot size cleaning process we have 2355 plots have been dropped
# plot number drecreased from 1188771  to 1186416  
#  don't do outliers cleaning based on absolute tree number in each plot
# currentYearPlotsTable = currentYearPlotsTable[currentYearPlotsTable$treeNumber >=5&currentYearPlotsTable$treeNumber <=2000,]

pdf(paste("Plots/HistGram_of_Plot_Size_Full_Data_All_Years.pdf",sep=""))
hist(currentYearPlotsTable$plotArea,main="Hist of plot size (ha)",xlab="Plot Size",ylab="Frequency")
dev.off()

# write the cleaned the data frame into local folder
write.csv(currentYearPlotsTable,"Data/BiomassCalculationData/GFBI_Plots_Biomass_Data_with_Full_Data_Year_Recording_20230126_PlotSize_Cleaned.csv")



#################################################################
## STEP 2 forest types definition and investigation year filter
#################################################################

# load the table we generated from STEP 1
cleanedPlotsTable = fread("Data/BiomassCalculationData/GFBI_Plots_Biomass_Data_with_Full_Data_Year_Recording_20230126_PlotSize_Cleaned.csv")[,-1] # no providing of the data

WWF_Biome = raster("Data/WWFRaster/WWF_Raster_Uniformed_Extent_Uniformed.tif")
# get the cordinates table to extract biome infromation
coordinatesTable = cleanedPlotsTable[,c("LON","LAT")]
extractedBiome = extract(WWF_Biome,coordinatesTable)
# add the biome information to the table
cleanedPlotsTable$WWF_Biome = extractedBiome
# duplicate the wwf_biome column as forest type column
cleanedPlotsTable$ForestType = cleanedPlotsTable$WWF_Biome
# define the forest types
# here dont to year filter
dim(cleanedPlotsTable) # 1115129      13
# write to local folder
write.csv(cleanedPlotsTable,"Data/BiomassCalculationData/GFBI_Plots_Biomass_Data_with_Full_Data_Year_Recording_20230126_Biome_Defined.csv")  # no providing of the data

#################################################################
## STEP 3 tree density based data cleaning
#################################################################

# load the table 
currentYearPlotsTable = fread("Data/BiomassCalculationData/GFBI_Plots_Biomass_Data_with_Full_Data_Year_Recording_20200520_Biome_Defined.csv")[,-1] %>% filter(!is.na(WWF_Biome))  # no providing of the data

# define the forest type vector
biomeVector = 1:14
biomeNames = c("Tropical moist","Tropical dry","Temperate broadleaf","Temperate coniferous","Boreal","Tropical savanna","Temperate savanna","Flooded savanna","Montane grassland","Tundra","Mediterranean forest","Desert")
pdf(paste("Plots/Figure_HistGram_of_Tree_density_for_Biome_Level_LOG_Transed_MAD_25_PlotSize_Filtered_20200525_FullData_AllYear.pdf",sep=""))
par(mfrow=c(4,3),oma=c(0,0,3,0))

for (bm in biomeVector)
{
    if (bm %in% c(1,2,4,5,6,7,8,9,10,11,12,13))
    {
        # subset the data for that biome 
        biomeSubsetTable = currentYearPlotsTable[currentYearPlotsTable$WWF_Biome == bm,]
        # !!!!!!!!!!!!!!!!!!!!!!!!!!transfer the unit from lg/ha biomass to t/ha Carbon. this was not applicable with new biome level carbon concentrations
        # biomeSubsetTable$biomassDensity = biomeSubsetTable$biomassDensity/2000
        # use the highest biomass density as the standard to delete the
        biomeSubsetTable$TreePerAreea = biomeSubsetTable$treeNumber/biomeSubsetTable$plotArea
        biomeSubsetTable$TreePerAreea = log(biomeSubsetTable$TreePerAreea+1)
        # print(nrow(biomeSubsetTable))
        madValue = mad(biomeSubsetTable$TreePerAreea)
        medianValue = median(biomeSubsetTable$TreePerAreea)
        largeRange = 2.5*madValue + medianValue
        # print(largeRange)
        smallRange = medianValue-2.5*madValue
        # print(smallRange)
        # subset data frame in mad range
        madFilteredTable = biomeSubsetTable[biomeSubsetTable$TreePerAreea <= largeRange&biomeSubsetTable$TreePerAreea >= smallRange,]
        # # make a histgram
        hist(biomeSubsetTable$TreePerAreea,main=paste(biomeNames[which(c(1,2,4,5,6,7,8,9,10,11,12,13)==bm)],sep=""),xlab="Tree density(log(N)/ha)",ylab="Frequency")
        abline(v = smallRange, col="red", lwd=1, lty=2)
        abline(v = largeRange, col="red", lwd=1, lty=2)
        # hist(log(madFilteredTable$biomassDensity+1),main=paste("Forest Type ",bm,sep=""),xlab="Biomass Density(LOG)",ylab="Frequency")
        # hist(exp(madFilteredTableNew$biomassDensity)-1,main=paste("Biome",biome,sep=""),xlab="Biomass Density(t/ha)",ylab="Frequency")
        print(paste("Biome: ",biomeNames[which(c(1,2,4,5,6,7,8,9,10,11,12,13)==bm)],sep=""))
        print((nrow(madFilteredTable)/nrow(biomeSubsetTable)))
        write.csv(madFilteredTable,paste("Data/BiomassCalculationData/TreeDensityFilter_SubTables/Full_Data_all_Year_25_PlotSize_Filtered_Biome_Level_Table_for_",bm,"_Biome_NON_Log.csv",sep=""))

    }else
    {
        # subset the data for that biome 
        biomeSubsetTable = currentYearPlotsTable[currentYearPlotsTable$WWF_Biome == bm,]
        biomeSubsetTable$TreePerAreea = biomeSubsetTable$treeNumber/biomeSubsetTable$plotArea
        biomeSubsetTable$TreePerAreea = log(biomeSubsetTable$TreePerAreea+1)
        print(paste("Biome: ",bm,sep=""))
        print((nrow(biomeSubsetTable)/nrow(biomeSubsetTable)))
        write.csv(biomeSubsetTable,paste("Data/BiomassCalculationData/TreeDensityFilter_SubTables/Full_Data_all_Year_25_PlotSize_Filtered_Biome_Level_Table_for_",bm,"_Biome_NON_Log.csv",sep=""))
    }
    
}

mtext("Log tranformed (filtered 2.5 MAD of treePerArea)", line=0, side=3, outer=TRUE, cex=1.5)
dev.off()



# [1] "Biome: Tropical moist"
# [1] 0.8485257
# [1] "Biome: Tropical dry"
# [1] 0.8901186
# [1] "Biome: Temperate broadleaf"
# [1] 0.9341215
# [1] "Biome: Temperate coniferous"
# [1] 0.9523375
# [1] "Biome: Boreal"
# [1] 0.9418345
# [1] "Biome: Tropical savanna"
# [1] 0.941896
# [1] "Biome: Temperate savanna"
# [1] 0.9438949
# [1] "Biome: Montane grassland"
# [1] 0.9969136
# [1] "Biome: Tundra"
# [1] 0.9597396
# [1] "Biome: Mediterranean forest"
# [1] 0.9810319
# [1] "Biome: Desert"
# [1] 0.9628302

# merge the forest type table into one
# rbind list those tables fitered 
fileNameVector = list.files(path="Data/BiomassCalculationData/TreeDensityFilter_SubTables",pattern = "_Biome_NON_Log.csv",full.name=T)
# load the table content into a list
fileList = lapply(fileNameVector,fread)
# rbind to a data frame 
rbindedTable = rbindlist(fileList)
# write to local folder
write.csv(rbindedTable,"Data/BiomassCalculationData/GFBI_Plots_Biomass_Density_Full_Data_TreeDensity_filtered_Recording_20230126_MAD_25_AllYear.csv")  

#################################################################
## STEP 4 Aggregate to pixel level by mean
#################################################################

library(data.table)
library(parallel)
library(dplyr)
library(raster)
library(rgdal)
# read the data
filteredTable = fread("Data/BiomassCalculationData/GFBI_Plots_Biomass_Density_Full_Data_TreeDensity_filtered_Recording_20200520_MAD_25_AllYear.csv")[,-1]
# filteredTable$biomassDensity = exp(filteredTable$biomassDensity)-1
initialRaster = raster(nrows=17400, ncols=43200, xmn=-180.0001, xmx=179.9999, ymn=-90.00014,ymx=83.99986,crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0",resolution=c(0.008333333, 0.008333333))
# backgroud=0, is an option to set the na values in the data absence cells
recordsRaster = rasterize(filteredTable[,c("LON","LAT")], initialRaster, subset(filteredTable,select="biomassDensity"), fun=mean)
# recordsRaster = rasterize(filteredTable[,c("LON","LAT")], initialRaster, subset(filteredTable,select="biomassDensity"), fun=median)
# write raster into the local folder
writeRaster(recordsRaster,paste("Data/BiomassCalculationData/Rasters/Full_Data_Retained_Plots_with_TreeDensity_25_MAD_Filter_Mean_after_AllYear.tif",sep=""),overwrite=T)


#################################################################
## STEP 5 Extract biome information
#################################################################

# transfer those rasters into data frames

recordsRaster = raster("Data/BiomassCalculationData/Rasters/Full_Data_Retained_Plots_with_TreeDensity_25_MAD_Filter_Mean_after_AllYear.tif")
# get the data frame for the cell whichi have values
rasterizedDataFram = na.omit(as.data.frame(recordsRaster,xy=TRUE))
# rename the data frame
names(rasterizedDataFram) = c("x","y","BiomassDensity")
print(dim(rasterizedDataFram))
# get the coordinates which has the 

WWF_Biome = raster("Data/WWFRaster/WWF_Raster_Uniformed_Extent_Uniformed.tif")
# get the cordinates table to extract biome infromation
coordinatesTable = rasterizedDataFram[,c("x","y")]
extractedBiome = extract(WWF_Biome,coordinatesTable)
# add the biome information to the table
rasterizedDataFram$WWF_Biome = extractedBiome
# duplicate the wwf_biome column as forest type column
# rasterizedDataFram$ForestType = rasterizedDataFram$WWF_Biome
# write to local drive
write.csv(rasterizedDataFram,"Data/BiomassCalculationData/Full_Data_Aggregated_Retained_Plots_with_TreeDensity_MAD_25_Filtered_Biomass_Density_Mean_20230126_AllYear.csv")


#################################################################
## STEP 6 Transform to shape file
#################################################################


rasterizedDataFram = fread ("Data/BiomassCalculationData/Full_Data_Aggregated_Retained_Plots_with_TreeDensity_MAD_25_Filtered_Biomass_Density_Mean_20230126_AllYear.csv")[,-1]
# duplicate the table 
duplicateTable = rasterizedDataFram

# allote the coordinates
coordinates(rasterizedDataFram) = ~x+y
# allocate projection system
proj4string(rasterizedDataFram) = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
rasterizedDataFram@data = duplicateTable  
# write the shapefile into local folder
writeOGR(rasterizedDataFram,dsn="Data/Full_Shapefile",layer="GFBI_Full_MAD25_AllYear_Mean_Aggre_Latest", driver = "ESRI Shapefile",overwrite=T)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# After this step, we uploaded the shapefile into google earth engine, then all the data will be open for the following steps.

