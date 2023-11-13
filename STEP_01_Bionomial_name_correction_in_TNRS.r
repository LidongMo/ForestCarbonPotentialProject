library(BIOMASS)
library(stringr)
library(ggplot2)
library(fields)
library(plyr)
library(gtools)
library(parallel)
library(dplyr)
library(RANN)
library(stringi)

# change to the path of the right dirctory

####################################################################################################
#STEP 1#
#SPECIES NAME CLEANING#

##################
###ORGINAL NAME###
##################

#Load the dbh data matrix#
dbh_path                     <- paste("Data/input_data/20180201_Tropical_DBH_Matrix_1.csv")
dbh_matrix                   <- read.csv(dbh_path)
#Load the dbh data based on individual level#
dbh_sp_path                  <- paste("Data/input_data/20180208_GFBI_DBH_Tree_List.csv")
dbh_sp_matrix                <- read.csv(dbh_sp_path)

#Load the wood density data matrix#
wd_path                      <- paste("Data/input_data/GlobalWoodDensityDatabase_data.csv")
wd_matrix                    <- read.csv(wd_path)

#Delete the dots in the binomial name from the sp_name_dbh dataframe#
sp_name_dbh                  <- gsub("\\.", " ",names(dbh_matrix)[-1:-3])
sp_name_wd                   <- as.character(wd_matrix[,3])
sp_individual_dbh            <- as.character(dbh_sp_matrix$SPCD)

sp_name_dbh_list             <- strsplit(as.character(sp_name_dbh)," ")
sp_name_wd_list              <- strsplit(as.character(sp_name_wd)," ")
sp_name_dbh_indi_list        <- strsplit(as.character(sp_individual_dbh)," ")

#There are 2 species binomial name sources, one is original,another is tnrs corrected#
#split the binomial name into columns of genus name and specie name#

orig_name_dbh_df             <- as.data.frame(t(stri_list2matrix(sp_name_dbh_list)))
orig_name_wd_df              <- as.data.frame(t(stri_list2matrix(sp_name_wd_list)))
orig_name_individual_df      <- as.data.frame(t(stri_list2matrix(sp_name_dbh_indi_list)))

#########################
###TNRS CORRECTED NAME###
#########################

#save the name list as data frame in CSV format#
dbh_df                       <- data.frame(binomial=sp_name_dbh)
write.csv(dbh_df,"Data/input_data/species_name_dbh.csv")

wd_df                        <- data.frame(binomial=sp_name_wd)
write.csv(wd_df,"Data/input_data/species_name_wd.csv")

individual_df                <- data.frame(binomial=sp_individual_dbh)
write.csv(individual_df,"Data/input_data/species_name_dbh_individual.csv")

#TNRS binomial correction on  http://tnrs.iplantcollaborative.org/index.html#

#load the tnrs matched binomial name and transfer it into 1 column data frame#
tnrs_dbh_path                <- paste("Data/input_data/tnrs_results_dbh.csv")
tnrs_sp_name_dbh             <- read.csv(tnrs_dbh_path)[-1,]

tnrs_wd_path                 <- paste("Data/input_data/tnrs_results_wd.csv")
tnrs_sp_name_wd              <- read.csv(tnrs_wd_path)[-1,]

tnrs_dbh_indi_path           <- paste("Data/input_data/tnrs_results_dbh_individual.csv")
tnrs_sp_name_dbh_indi        <- read.csv(tnrs_dbh_indi_path)[-1,]

tnrs_gfbi_path               <- paste("Data/input_data/tnrs_results_GFBI_list.csv")
tnrs_gfbi                    <- read.csv(tnrs_gfbi_path)

cor_name_dbh                 <- as.character(tnrs_sp_name_dbh[,2])
cor_name_wd                  <- as.character(tnrs_sp_name_wd[,2])
cor_individual_dbh           <- as.character(tnrs_sp_name_dbh_indi[,2])
cor_name_gfbi                <- as.character(tnrs_gfbi[,2])

#There are 2 species binomial name sources, one is original,another is tnrs corrected#
#split the binomial name into columns of genus name and specie name#
cor_name_dbh_list            <- strsplit(as.character(cor_name_dbh)," ")
cor_name_wd_list             <- strsplit(as.character(cor_name_wd)," ")
cor_name_dbh_indi_list       <- strsplit(as.character(cor_individual_dbh)," ")
cor_name_gfbi_list           <- strsplit(as.character(cor_name_gfbi)," ")


cor_name_dbh_df              <- as.data.frame(t(stri_list2matrix(cor_name_dbh_list)))
cor_name_wd_df               <- as.data.frame(t(stri_list2matrix(cor_name_wd_list)))
cor_name_individual_df       <- as.data.frame(t(stri_list2matrix(cor_name_dbh_indi_list)))
cor_name_gfbi_df             <- as.data.frame(t(stri_list2matrix(cor_name_gfbi_list)))


#combine the binomial dataframe with the TNRS coreccted names' data freame#
GFBIRawAndCorrectedName      <- cbind(tnrs_gfbi[,1],cor_name_gfbi_df[,1:2])
#name the columns#
names(GFBIRawAndCorrectedName) <- c("RawName","CorrectedGenus","CorrectedSpecies")
write.csv(GFBIRawAndCorrectedName,"Data/input_data/GFBI_raw_and_tnrs_corrected_binomial_df.csv")



####################################################################################################

#SETP 2#
#GET THE WOOD DENSITY OF EACH SPECIES#

###########################
###ORIGINAL WOOD DENSITY###
###########################

#the binomial name of the wood density database not be corrected#
origin_wd                    <- getWoodDensity(genus=orig_name_dbh_df[,1],species=orig_name_dbh_df[,2])
#this is the sample data with corrected name#
origin_wd_indi               <- getWoodDensity(genus=orig_name_individual_df[,1],species=orig_name_individual_df[,2])

############################
###CORRECTED WOOD DENSITY###
############################

#the name of the wood density database have been corrected by TNRS#
#replace the binomial in the wood density data base with tnrs corrected binomial#
correc_wd_matrix             <- wd_matrix
correc_wd_matrix             <- cbind(correc_wd_matrix[,1:2],cor_name_wd_df[,1:2],correc_wd_matrix[,4:6])
correc_wd_matrix             <- correc_wd_matrix[,-1]
names(correc_wd_matrix)      <- c("family", "genus", "species","wd","Region","Reference.Number")
write.csv(correc_wd_matrix,"Code/Resources/name_corrected_wooddensity_database.csv",row.names = FALSE)

#modiciaction the wood density data sources of the original getWoodDensity function in BIOMASS packages#
source("Code/Resources/getWoodDensity_MO.r")
# corrected_wd               <- getWoodDensity_MO(genus=cor_name_dbh_df[,1],species=cor_name_dbh_df[,2])
# corrected_wd_indi          <- getWoodDensity_MO(genus=cor_name_individual_df[,1],species=cor_name_individual_df[,2])
corrected_wd_gfbi            <- getWoodDensity(genus=cor_name_gfbi_df[,1],species=cor_name_gfbi_df[,2])
table(corrected_wd_gfbi$levelWD)

# cbind the GFBIRawAndCorrectedName data frame with this wood density data frame#
corrected_wd_gfbi            <- cbind(GFBIRawName=GFBIRawAndCorrectedName[,1],corrected_wd_gfbi)
################################
#use the original function to get the wood density value, which is the right one we are using for calculation#
write.csv(corrected_wd_gfbi,"Data/output_data/original_wd_of_all_spsecies_GFBI.csv")
table(corrected_wd_gfbi$levelWD)
################################
# Result from the old BIOMASS version before 2018
#dataset   genus species
#6849   16738    4038

# This is the result from the latest version 
# dataset   genus species 
#    6853   16696    4076 

# In the following calculation, we still use the old result that was applied in the biomass calculation