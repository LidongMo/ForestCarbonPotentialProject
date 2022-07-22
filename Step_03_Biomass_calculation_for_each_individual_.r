# Biomass caluclation for each individual in GFBi
#biomass calculation#
#for the tropical regions we use the biomass estimation eqaution from the "BIOMASS" package which is associated with wood density, whice has already been prepared#
#for the other regions, we apply the equations trained by ourself#


library(stringr)
library(parallel)
library(dplyr)
library(BIOMASS)

# this the work dirctory setting in big computer
setwd("~/Github_Folder")
# this part is the operation on the raw data. 
# Allocte plot ID to tree level data of GFBI
# allote the WWF biome information as well

# Data_PlotID                      <- fread("/Volumes/CrowtherLabRAID/GFBI_Data/GFBI_Original_Data/Global_treelist_06222017.csv")[,-1]
# Data_Biome                       <- fread("/Volumes/CrowtherLabRAID/GFBI_Data/Tree_Level_Data/20180424_GFBI_Treelist_Data_wBiomes.csv")[,-1]

# Data_PlotID$New                  <- paste(Data_PlotID$LAT,Data_PlotID$LON,Data_PlotID$SPCD,Data_PlotID$DBH,Data_PlotID$FID,sep="_")
# Data_Biome$New                   <- paste(Data_Biome$LAT,Data_Biome$LON,Data_Biome$SPCD,Data_Biome$DBH,Data_Biome$FID,sep="_")
# MergeData                        <- merge(Data_PlotID[,c("PLT","TPH","New")],Data_Biome,by="New")
# DataBiomePLT                     <- MergeData[,-"New"]
# write.csv(DataBiomePLT,"/Volumes/CrowtherLabRAID/GFBI_Data/Tree_Level_Data/20180827_GFBI_Treelist_Data_wBiomes_PLT_TPH_Lidong.csv")
# write.csv(DataBiomePLT,"Global_treelist_Biome_PLT_TPH.csv")


WoodDensityData                  <- read.csv("BiomassCalculationData/original_wd_of_all_spsecies_GFBI.csv") #the file path should be changed if you want to apply the code on other machine#
# this the normal approach for reltively small file size csv loading#
#for the loading of the extremely large csv data frame, the most effictive approach is the function
#from the package "data.table" and the function is "fread"#
#l
#this is the data frame of the raw binomials and the TNRS corrected genus and species name#
GFBIRawAndCorrectedName          <- read.csv("BiomassCalculationData/GFBI_raw_and_tnrs_corrected_binomial_df.csv")[,-1]
# read the ecoregion equation data frame
EcoregionEquation                <- read.csv("BiomassCalculationData/pesudo_biomass_based_ecoregion_equations.csv")

# read the original GFBI data and subset as example data frame
# we used the top 50,000 lines as the data for the example 
# Original_GFBI_full             <- fread("Global_treelist_Biome_PLT_TPH.csv") [,-1]
# write.csv(Original_GFBI_full[1:50000,],"Global_treelist_Biome_PLT_TPH_Subset.csv")

Original_GFBI                    <- fread("BiomassCalculationData/Global_treelist_Biome_PLT_TPH_Subset.csv") 
# remove the rows with NAs, specifical for the plots which located out side of the WWF_Biome mask
NAomit_GFBI                      <- na.omit(Original_GFBI)
# remove the rows which has the WWF_Biome values 99 and 98
CleanedGFBI                      <- subset(NAomit_GFBI,NAomit_GFBI$WWF_Biome!=98&NAomit_GFBI$WWF_Biome!=99)
# kick out those individuals whose diameter is less than 5
FiveDeleteGFBI                   <- subset(CleanedGFBI,CleanedGFBI$DBH>=5)
# split the big data frame into samller pieces
SplitList                        <- split(FiveDeleteGFBI, (as.numeric(rownames(FiveDeleteGFBI))-1) %/% 10000)
# for (i in names(SplitList))
# {
# 	write.csv(SplitList[[as.character(i)]], file = paste("GFBI_Data_",i,"_piece.csv", sep = ""))
# }
# in order to inprove the calcualtion speed, we split the big data frame into smalle pieces and name them
lapply(names(SplitList), function(x){write.csv(SplitList[[as.character(x)]], file = paste("BiomassCalculationData/Split_After_Cleaning/Split_GFBI_Data_",x,"_piece.csv", sep = ""))})

# get the file neme list in the split data frame folder
FileList                             <- list.files(path="BiomassCalculationData/Split_After_Cleaning",pattern="*.csv")
for (fl in FileList)
# for (fl in FileList)
{
	IndividualDataDBH                <- fread(paste("BiomassCalculationData/Split_After_Cleaning/",fl,sep=""))[,-1]
	# ListElement                      <- which(FileList==fl)
	BiomassCalculation               <- function(individual)
    {
    	#got the ecoregion code for each tree#
	    EcoregionCode                <- IndividualDataDBH[individual,]$WWF_Biome
	    #there is no allometic equations trained for 14 Mangroves and 3 Tropical and Subtropical Coniferous Forests#
	    if (EcoregionCode!=3&EcoregionCode!=14)
	    {
		    #get the equation information for each specific ecoregion code#
		    EquationInformation      <- EcoregionEquation[EcoregionEquation$ID_WWF==EcoregionCode,]
		    #get the individual information for each individual in the  IndividualDataDBH data frame#
		    IndividualInformation    <- IndividualDataDBH[individual,]
		    # get the dbh value for this individual
		    dbh                      <- IndividualInformation$DBH
		    # get the eqaution for the ecoregion where the individual is#
		    EquationFormula          <- tolower(as.character(EquationInformation$Eqaution))
            #Find the "=" symbol in the equation string,the charater after "=" will be the start of the equation#
            StartPosition            <- regexpr("=",EquationFormula)[1]+1
            #Find the end point in the equation string#
            EndPostion               <- nchar(EquationFormula)
            #Get the final equation for following calculations#
            EqautionForCalculation   <- str_sub(EquationFormula, start = StartPosition, end = EndPostion)
            #calculate out the individual biomass#
            IndividualBiomass        <- exp(eval(parse(text=EqautionForCalculation)))
            OutputDataFrame          <- data.frame(select(IndividualInformation,-Year),BiomassEquationBased=IndividualBiomass,BiomassPackageBased=NA)
            #for the tropical trees, we also apply the calculation approach from BIOMASS package#
            if(EcoregionCode!=1&EcoregionCode!=2&EcoregionCode!=7)
            {
        	    # for the individuals which come from the non-tropical region, we do not apply the computeAGB function to get the BIOMASS package based biomass, then set it to NA@
        	    OutputDataFrame$BiomassPackageBased <-NA
            }else
            {
        	    # for the individuals come from the tropical region, we will apply the computeAGB function from BIOMASS
        	    SpeciesBinomial      <- as.character(IndividualInformation$SPCD)
        	    # get the wood density row of the specific species
        	    WoodDensityRow       <- WoodDensityData[WoodDensityData$GFBIRawName==SpeciesBinomial,]
        	    # compute the biomass with the computeAGB equation from BIOMASS
        	    OutputDataFrame$BiomassPackageBased <-computeAGB(dbh,WoodDensityRow$meanWD,coord=cbind(IndividualInformation$LON,IndividualInformation$LAT))*1000
            }
        }else
        # (EcoregionCode==3|14)
        {
    	    # for the individuals come from the two ecoregions which do not hava any equations to get the biomass calculation, we set the biomass value to NA
    	    IndividualInformation    <- IndividualDataDBH[individual,]
    	    OutputDataFrame          <- data.frame(select(IndividualInformation,-Year),BiomassEquationBased=NA,BiomassPackageBased=NA)
        }

        print(paste("---the ",individual,"th line of ",nrow(IndividualDataDBH),",has been calculated---",sep=""))

        return(OutputDataFrame)
    }
    # apply parallel calculation of all the individuals in the GFBI database
    # be carfull of the input for the arguement "mc.cores", Please check how many cores in your own computer
    system.time(OutputLists          <- mclapply(1:nrow(IndividualDataDBH), BiomassCalculation, mc.cores = 16, mc.preschedule=FALSE)) 
    OutputDF                         <- do.call(rbind, OutputLists)
    write.csv(OutputDF,paste("BiomassCalculationData/Split_After_Cleaning_Biomass/individual_biomass_for_",fl,sep=""))
    print(paste("-----the ",fl,"has been calculated-----",sep=""))
}



#STEP 10#
# transfer the data from the big computer to my local drive 
# scp  standardlogin@uwis-cx-dock-11-040.ethz.ch:/Users/Shared/Lidong_Mo/Split_Data_GFBI/individual_biomass_of_GFBI_data_for_ 6 _split.csv /Users/LeonidMoore/Desktop/BIOMASS/individual_biomass_of_GFBI_data_for_ 6 _split.csv

FileList           <- list.files(path="BiomassCalculationData/Split_After_Cleaning_Biomass",pattern="*.csv")
# use lapply to combine all the data into a singe object in R
ReunionDataList                 <- lapply(paste("BiomassCalculationData/Split_After_Cleaning_Biomass/",FileList,sep=""),fread)
ReunionData                     <- do.call(rbind,ReunionDataList)
# write the data cleaned biomass calculation result to the local folder
write.csv(ReunionData[,-1:-2],"BiomassCalculationData/GFBI_Biomass_Calcualtion_Result_Data_Frame_0522.csv")

# ###################################################
# STEP 11
# merge the data frame with TPH with the one just has the biomass value
# Load the table with TPH
DataBiomePLTTPH               <- fread("BiomassCalculationData/Global_treelist_Biome_PLT_TPH_Subset.csv")[,-1]
# load the biomass data
BiomassBiome                  <- fread("BiomassCalculationData/GFBI_Biomass_Calcualtion_Result_Data_Frame_0522.csv")[,-1]

# add new column to the data frame
DataBiomePLTTPH$New           <- paste(DataBiomePLTTPH$LAT,DataBiomePLTTPH$LON,DataBiomePLTTPH$SPCD,DataBiomePLTTPH$DBH,DataBiomePLTTPH$FID,sep="_")
BiomassBiome$New              <- paste(BiomassBiome$LAT,BiomassBiome$LON,BiomassBiome$SPCD,BiomassBiome$DBH,BiomassBiome$FID,sep="_")
# merge the two data frames by the "new" column
MergeData                     <- merge(DataBiomePLTTPH[,c("Year","New")],BiomassBiome,by="New")
# delete the column "New"
BiomassBiomeTPH               <- MergeData[,-"New"]

write.csv(BiomassBiomeTPH,"BiomassCalculationData/GFBI_Biomass_Calcualtion_Result_with_TPH_Data_Frame_20181008.csv")

# load the bimass data frame for plots only have 1 year record
SingleYearTreeLevelData        <- fread("BiomassCalculationData/GFBI_Biomass_Calcualtion_Result_with_TPH_Data_Frame_20181008.csv")
# do not do the NA kick out process
ReunionData                    <- SingleYearTreeLevelData
# kick out the data in biome 3 and 14
SubDataFrame                   <- ReunionData[ReunionData$WWF_Biome!=3&ReunionData$WWF_Biome!=14,]
# subset the individual biomass for the tropcial region
TropicalDataFrame              <- SubDataFrame[SubDataFrame$WWF_Biome==1|SubDataFrame$WWF_Biome==2|SubDataFrame$WWF_Biome==7,]
TropicalDataFrame              <- TropicalDataFrame[,c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","BiomassPackageBased")]
TropicalDataFrame              <- na.omit(TropicalDataFrame)
names(TropicalDataFrame)       <- c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","Biomass")
write.csv(TropicalDataFrame,"BiomassCalculationData/Tropical_Individual_Biomass_Data_Frame_20181012.csv") 

# subset the data frame for the non tropical region and allocate the name Biomass to the Biomasspackagebased column
NonTropicalDataFrame           <- SubDataFrame[SubDataFrame$WWF_Biome!=1&SubDataFrame$WWF_Biome!=2&SubDataFrame$WWF_Biome!=7,]
# select the column of biomass for the nontropical region individuals
NonTropicalDataFrame           <- NonTropicalDataFrame[,c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","BiomassEquationBased")]
NonTropicalDataFrame           <- na.omit(NonTropicalDataFrame)
# change the name for the biomass column
names(NonTropicalDataFrame)    <- c("FID","PLT","LAT","LON","Year","SPCD","DBH","TPH","WWF_Biome","Biomass")
write.csv(NonTropicalDataFrame,"Non_Tropical_Individual_Biomass_Data_Frame_20181012.csv") 
# rbind the tropical and notropical region individual biomass data 
FinalIndividualBiomass         <- rbind(TropicalDataFrame,NonTropicalDataFrame)
# write the tropical and nontropical individual biomass data to the local directory
write.csv(FinalIndividualBiomass,"BiomassCalculationData/GFBI_Biomass_Calcualtion_Result_Reunion_Data_Frame_20181012.csv")


# ######################
# ######################

library(data.table)
library(stringr)
library(parallel)
library(dplyr)
library(BIOMASS)
library(reshape2)
# because there are some duplicates in the PLT, which is not good for our calculation, We have add our own PLT to the data
# check the unique plots by the lat and lon
# genewrete a new column by the combination of lat and lon
FinalIndividualBiomass         <- fread("BiomassCalculationData/GFBI_Biomass_Calcualtion_Result_Reunion_Data_Frame_20181012.csv")[,-1]
FinalIndividualBiomass$new      <- paste(FinalIndividualBiomass$LAT,"_",FinalIndividualBiomass$LON,"_",FinalIndividualBiomass$PLT,sep="")
# find the out the unique lat and lon combination
UniqueCoordinates               <- unique(FinalIndividualBiomass$new)
# ADD NEW PLT INTO EACH COORDINATES PLT
NewPLTAllocating                <- function(uc)
# for (uc in UniqueCoordinates)
{
    # get the order ot the uc in the uniquecorrdinates
    UCOrder                     <- which(UniqueCoordinates==uc)
    # generate new plot id by the order of the uc and allocate this to the data frame
    NewPlotID                   <- paste("PLT_",str_pad(UCOrder, 8, pad = "0"),sep="")
    # subset the data in that coordinates plot
    PerCoordinate               <- FinalIndividualBiomass[FinalIndividualBiomass$new==uc,]
    PerCoordinateNew            <- data.frame(PerCoordinate,NewPLT=NewPlotID)
    PerCoordinateOutput         <- PerCoordinateNew[,c("FID","PLT","LAT","LON","Year","NewPLT","SPCD","DBH","TPH","WWF_Biome","Biomass")] 
    # print(UCOrder)
    # return(PerCoordinateOutput)
    write.csv(PerCoordinateOutput,paste("BiomassCalculationData/PerPLT/PLT_",NewPlotID,".csv",sep=""))
    print(UCOrder)
}

system.time(OutputList <- mclapply(UniqueCoordinates,NewPLTAllocating,mc.cores=12,mc.preschedule=F))
# OutputDataFrame                   <- rbindlist(OutputList)
# write.csv(OutputDataFrame,"GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181022.csv")
# the noted codes are for the classic calculation

# rbind all the data frame in the PerPLT folder
# list all the data frame in the perplt folder
DataFrameList                    <- list.files(path="BiomassCalculationData/PerPLT",pattern=".csv")
# load all the data frame into the memory
ReunionDataList                 <- lapply(paste("BiomassCalculationData/PerPLT/",DataFrameList ,sep=""),read.csv)
ReunionData                     <- rbindlist(ReunionDataList)
# write the data cleaned biomass calculation result to the local folder
fwrite(ReunionData[,-1],"BiomassCalculationData/GFBI_Biomass_Tree_Level_Biomass_New_PLT_Data_Frame_20181024.csv")