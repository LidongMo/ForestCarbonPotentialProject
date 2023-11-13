# get the 
library(data.table)
library(parallel)
library(dplyr)
# set the working directory by cmd 
# set the R version

BiomassBiomeTPH = fread("Data/BiomassCalculationData/GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181024.csv")[,-1]
# get the biomass per area for each individual
# BiomassBiomeTPH$BiomassArea     <- BiomassBiomeTPH$Biomass*BiomassBiomeTPH$TPH
# load the wood density data table
woodDensityData = read.csv("Data/BiomassCalculationData/original_wd_of_all_spsecies_GFBI.csv") 
# get the plot names
plotNames = as.vector(unique(BiomassBiomeTPH$NewPLT))

plotLevelOperationFunc = function(pn=PlotName,inputPlot = InputPlot)
{
    
    # get the subset data for each plot by plot name
    perPlotDataFrame = inputPlot
    # get the total tree number in per plot
    treeNumber = nrow(perPlotDataFrame)
    # plot size by hectare
    plotArea = 1/(mean(perPlotDataFrame$TPH))
    # total biomass
    totalBiomass = sum(perPlotDataFrame$Biomass)
    # biomass density kg/ha
    biomassDensity = totalBiomass/plotArea
    sumBasalArea = sum((perPlotDataFrame$DBH/2)^2)*pi
    arverageBasalArea = sumBasalArea/plotArea
    treeNumberArea = treeNumber/plotArea
    # mean biomass per tree
    # meanTreeBiomass = totalBiomass/ treeNumber
    # # sd of biomass per tree
    # sdBiomass = sd(perPlotDataFrame$Biomass)
    # # sd of DBH per tree
    # sdDBH = sd(perPlotDataFrame$DBH)
    # # mean DBH in each plot
    # meanDBH = mean(perPlotDataFrame$DBH)
    # woodDenstiyAllocateFUnc = function(x)
    # {
    #     # for the individuals come from the tropical region, we will apply the computeAGB function from BIOMASS
    #     speciesBinomial = as.character(x[6])
    #     # get the wood density row of the specific species
    #     woodDensityRow = woodDensityData[woodDensityData$GFBIRawName==speciesBinomial,]
    #     woodDensity = mean(woodDensityRow$meanWD)
    #     return(woodDensity)
    # }
    # # get the wood density data for each tree in the plot
    # woodDensityVector = apply (perPlotDataFrame,1,woodDenstiyAllocateFUnc)
    # # get the mean of wood density for each plot
    # meanWoodDensity = mean(woodDensityVector)
    # # get the sd of the wood density in each plot
    # sdWoodDensity = sd(woodDensityVector)
    # get the lat and lon for that tph in that year
    Lati = unique(perPlotDataFrame$LAT)
    Long = unique(perPlotDataFrame$LON)
    yr = unique(perPlotDataFrame$Year)
    # return the information row out
    outputRow = data.frame(PLT=pn,
                          LAT=Lati,
                          YEAR=yr,
                          LON=Long,
                          treeNumber = treeNumber,
                          plotArea = plotArea,
                          totalBiomass = totalBiomass,
                          biomassDensity = biomassDensity,
                          sumBasalArea = sumBasalArea,
                          averageBasalArea = arverageBasalArea,
                          treeNumberArea = treeNumberArea)
    print(paste("--- the metadata for plot ",pn," has been calculated ---"))
    return(outputRow)   
}

# in order to get the plot with time series, we just to check the plot with multiple years records
MultipleYearPlotFind            <- function(pn)
{
    # get the per plot data 
    perPlotDataFrame = fread(paste("Data/BiomassCalculationData/PerPLT/",pn,".csv",sep=""))[,-1]
    # check is the plot has more than one years observations
    startDatFrame = data.frame()
    # get the latest year
    maxYear = max(unique(perPlotDataFrame$Year))
    # subset the max year data frame
    yearlyDataFrame = perPlotDataFrame[perPlotDataFrame$Year == maxYear,]

    if (length(unique(perPlotDataFrame$Year)) > 1)
    {
        # get the years in that plot
        yearsVector = unique(perPlotDataFrame$Year)
        # lopp by year
        for (yr in yearsVector)
        {
            # subset the yearly data frame 
            yearlyDataFrame = perPlotDataFrame[perPlotDataFrame$Year == yr,]
            startDatFrame = rbind(startDatFrame,plotLevelOperationFunc(pn,yearlyDataFrame))
        }
    }else
    {
        startDatFrame = rbind(startDatFrame,plotLevelOperationFunc(pn,yearlyDataFrame))
    } 
    return(startDatFrame) 
}

system.time(OutputList  <- mclapply(plotNames,MultipleYearPlotFind,mc.cores = 4, mc.preschedule=FALSE))
OutputDF = rbindlist(OutputList)
# write the data frame into local folder
write.csv(OutputDF,"Data/BiomassCalculationData/GFBI_Plots_with_All_Years_Time_Series_Records_20190412_Lidong.csv")

