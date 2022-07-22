

library(BIOMASS)
library(stringr)
library(ggplot2)
library(fields)
library(plyr)
library(gtools)
library(parallel)
library(dplyr)
library(RANN)

setwd("~/Github_Folder/")
#allometry function database cleaning#
allometry_func_path          <- paste("input_data/allometric_equations.csv")
allom_func                   <- read.csv(allometry_func_path)
#At first check how many diffenrent Xs for the allomety functions#
X_df                         <- data.frame(table(allom_func$X))
#After check the variety of X, 55 kinds of X metrics have been found#
#We just need the DBH associated functions (there are some functions like Log10(DBH) etc) will be condidered for next manipulation#
#Keep the allometry functions based only on DBH#


sub_allom_func               <- subset(allom_func,allom_func$X=="DBH"|allom_func$X=="DBH^2"|allom_func$X=="DBH1"|allom_func$X=="Log10(DBH)")
#table the Component type of the eqautions and transfer the results to a data frame#
component_type               <-as.data.frame(table(sub_allom_func$Veg_Component))
#write the result data frame to a csv foramt file#
write.csv(component_type,"output_data/component_type.csv")

# after we acquired the biomass component type, we only retain the types for the AGB, this was manually done
selectedTypes = fread("output_data/selected_component_type.csv")[,-1]
# and then subset the allom_func by this data frame
sub_allom_func_true = sub_allom_func %>% filter(Veg_Component %in% selectedTypes$Var1)

#Furthermore, some allometry functions are build based on DBH and wood density or height. #
#However, we just need the functions only based on DBH#
#For the data operating reasons, I did the replacement of the NA values to any self defined non-NA factor, and then do subset based on this replacement#
sub_allom_func_true$Z        <- as.character(sub_allom_func_true$Z)
sub_allom_func_true$Z[is.na(sub_allom_func_true$Z)] <- "keep"
sub_allom_func_dbh           <- subset(sub_allom_func_true,sub_allom_func_true$Z=="keep")

#Delete the columns which contains useless information for our data manipulation#
#We keep the colmns named "ID_AE","Latitude","Longitude","Family","Genus","Species","X","Unit_X","Min_X","Max_X","Output","Output_TR","Unit_Y","Equation","Substitute_equation","Sample_size","Population","Location","Country","Ecoregion_WWF","Reference"#
SelectedColumn               <- c("ID_AE","Latitude","Longitude","Family","Genus","Species","X","Unit_X","Min_X","Max_X","Output","Output_TR","Unit_Y","Equation","RMSE","SEE","Corrected_for_bias","Substitute_equation","Sample_size","Veg_Component","Population","Location","Country","Ecoregion_WWF","Reference")
selected_dbh_func            <- subset(sub_allom_func_dbh,select=SelectedColumn)
sub_allom_func_biomass       <- subset(selected_dbh_func,selected_dbh_func$Output=="Biomass")
write.csv(sub_allom_func_biomass,"input_data/full_information_of_selected_allometry_functions.csv")

#Load the subset of the DBH base function data frame#
sub_allom_func_biomass           <- read.csv("input_data/full_information_of_selected_allometry_functions.csv")[,-1]

#The equations from the allometry equation database are not uniformed#
#They might have different X,differen output transformation or different input and output Unit#
#Check those dimensions first#
#The output Units are: g,kg,kg/ha,kg/tree,lb,Mg,Mg/ha#
#We have to delete the per tree or per hectre based output Unit_Y functions: kg/ha,kg/tree,Mg/ha#
sudo_allom_func                  <- subset(sub_allom_func_biomass,sub_allom_func_biomass$Unit_Y=="g"|sub_allom_func_biomass$Unit_Y=="kg"|sub_allom_func_biomass$Unit_Y=="lb"|sub_allom_func_biomass$Unit_Y=="Mg")
#delete the attributes with empty function#
sudo_allom_func                  <- subset(sudo_allom_func,sudo_allom_func$Substitute_equation!="#NAME?")
write.csv(sudo_allom_func,"input_data/sudo_allom_functions_raw.csv")
#function for bioamss calculation based each allometry equation#

#delete the equations with "none" inside, if an euqation has a "none", that means it is not a DBH only based equations#
sudo_allom_func_none              <- read.csv("input_data/sudo_allom_functions_raw.csv")
sudo_allom_func_final             <- subset(sudo_allom_func_none,!(grepl("None",sudo_allom_func_none$Substitute_equation)))
#delete the equations without aplliable min and max#
#Define the function which can delete the rows with NAs in speciefic columns#
#function is from the https://stackoverflow.com/questions/11254524/omit-rows-containing-specific-column-of-na #

SpecificNARemove                  <- function(data, desiredCols)
{
    completeVec                   <- complete.cases(data[, desiredCols])
    return(data[completeVec, ])
}

sudo_allom_func_final             <- SpecificNARemove(sudo_allom_func_final,c("Min_X","Max_X"))


#Write the function data frame with the simple information: "ID_AE","Latitude","Longitude","Family","Genus","Species","Ecoregion_WWF"#
FinalSelected                <- c("ID_AE","Latitude","Longitude","Family","Genus","Species","Substitute_equation","Unit_X","Min_X","Max_X","Output_TR","Unit_Y","Ecoregion_WWF","Veg_Component")
simple_dbh_func              <- subset(sudo_allom_func_final,select=FinalSelected)

write.csv(simple_dbh_func,"input_data/final_simple_information_of_selected_allometry_functions.csv")


#Check the minimum and maximum value of DBH in the allometry function database#
#subset the data frame, and keep the columns of X, Unit_X, Min_X and Max_X#
sudo_allom_func_final            <- read.csv("input_data/final_simple_information_of_selected_allometry_functions.csv")

sub_min_max                      <- sudo_allom_func_final[,9:12]
sub_min_max                      <- na.omit(sub_min_max)

unit_uniformed                   <- sub_min_max

#the function which can tranfor the minimum and maximum value of Min_X and Max_X, and return it#

for (i in 1:nrow(sub_min_max))
{
    if (sub_min_max[i,2]=="cm")
    {
        unit_uniformed[i,3]       <- sub_min_max[i,3]
        unit_uniformed[i,4]       <- sub_min_max[i,4]
    }else
    if (sub_min_max[i,2]=="cm x 10")
    {
        unit_uniformed[i,3]       <- sub_min_max[i,3]/10
        unit_uniformed[i,4]       <- sub_min_max[i,4]/10
    }else
    if (sub_min_max[i,2]=="in")
    {
        unit_uniformed[i,3]       <- sub_min_max[i,3]*2.54
        unit_uniformed[i,4]       <- sub_min_max[i,4]*2.54
    }else
    if (sub_min_max[i,2]=="in2")
    {
        unit_uniformed[i,3]       <- sqrt(sub_min_max[i,3])*2.54
        unit_uniformed[i,4]       <- sqrt(sub_min_max[i,4])*2.54
    }else
    if (sub_min_max[i,2]=="m")
    {
        unit_uniformed[i,3]       <- sub_min_max[i,3]
        unit_uniformed[i,4]       <- sub_min_max[i,4]
    }else
    if (sub_min_max[i,2]=="mm")
    {
        unit_uniformed[i,3]       <- sub_min_max[i,3]/10
        unit_uniformed[i,4]       <- sub_min_max[i,4]/10
    }else
    if (sub_min_max[i,2]=="mm2")
    {
        unit_uniformed[i,3]       <- sqrt(sub_min_max[i,3])/10
        unit_uniformed[i,4]       <- sqrt(sub_min_max[i,4])/10
    }
}
#write the transferomed minimum and maximum value to a csv format file#
write.csv(unit_uniformed,"input_data/uniformed_unite_min_max_x_cm.csv")

############################
#PESUDO BIOMASS CALCULATION#
############################

#Data cleaning for some not real DBH only based equations#
#Those functions have a component inside the equation string is "none"#
#Kick out those equations contain "none"#
source("Resources/sudo_biomass.r")

# get the unique 
sudo_allom_func_final = read.csv("input_data/final_simple_information_of_selected_allometry_functions.csv")[,-1] 

# sudo_allom_func_final$Ecoregion_WWF = unlist(lapply(sudo_allom_func_final$Ecoregion_WWF, as.character))

# check how many equations in the non-tropical biomes 
NonTropical_Sudo_allometirc_Eq = sudo_allom_func_final %>% dplyr::filter(Ecoregion_WWF %in% c("Temperate broadleaf and mixed forests","Temperate coniferous forests","Boreal forest/taigas","Tropical and subtropical grasslands, savannas, and shrublands","Temperate grasslands, savannas, and shrublands","Montane grasslands","Tundra","Mediterranean scrub","Deserts and xeric shrublands")) #3090 equations with duplicates in non-tropical regions

#Check the distribution of the eqautions in diverse ecoregions#
DistributionEquations = as.data.frame(table(sudo_allom_func_final$Ecoregion_WWF))
write.csv(DistributionEquations,"output_data/distribution_of_allometric_eqautions_ecoregion.csv")

#in this part we should take care of some mistakes in "()" in the fucntions, I did this manualy and modified the equations by hand#
#for (ron in 1:nrow(sudo_allom_func_final))
multi_sudo_biomass = function(ron)
{
    print(ron)
    #According to the Jenkins 2003‘s method#
    sudo_data_vector =c(rep(5:24,1),seq(25,100,5),seq(110,300,10))
    name_list = paste("DBH",sudo_data_vector,sep="")
    per_func_row = NonTropical_Sudo_allometirc_Eq[ron,]
    sudo_biomass(InputDf=per_func_row,SudoData=sudo_data_vector,NameList=name_list)
}

mycores = detectCores()
cores =mycores-2

#parallel calculation of the embed function multi_prediction#
system.time(output_list <- mclapply(1:nrow(NonTropical_Sudo_allometirc_Eq), multi_sudo_biomass, mc.cores = 4, mc.preschedule=FALSE))
#rbind the pesudo data calculation results into a data frame#

output_df= rbindlist(output_list)

#it looks we have many repeated equations, which is meaningless.Therefore, we have to kick them out! #
RepeatKickout = as.data.frame(output_df %>% distinct(Equations, .keep_all = T))
#write the calculation results in to local file as csv format#
write.csv(RepeatKickout,"output_data/primary_sudo_biomass_data_repeat_kickout.csv")

###############
#MODELING PART#
###############
#this the original pesudo biomass data withe those strange and extremly large values#
PesudoData_df                    <- read.csv("output_data/primary_sudo_biomass_data_repeat_kickout.csv")
#this is the tail kicked out data,without those strange and extremly large values#
PesudoData_df =read.csv("output_data/final_pesudo_data_repeat_kickout.csv") %>% filter(Ecoregion_WWF!="Tropical and subtropical moist broadleaf forests"&Ecoregion_WWF!="Flooded grasslands"&Ecoregion_WWF!="Tropical and subtropical dry broadleaf forests"&Ecoregion_WWF!="Tropical and subtropical coniferous forests"&Ecoregion_WWF!= 'NA'&Ecoregion_WWF!="Mangroves") %>% mutate(Ecoregion_WWF=as.character(Ecoregion_WWF))
#get the ecoregion name lists in side the pesudo biomass data frame#
ecoregions = as.vector(data.frame( table(PesudoData_df$Ecoregion_WWF))[,1])
# change the order of the ecoregions
ecoregionsOrdered = c("Temperate broadleaf and mixed forests","Temperate coniferous forests","Boreal forest/taigas","Tropical and subtropical grasslands, savannas, and shrublands","Temperate grasslands, savannas, and shrublands","Montane grasslands","Tundra","Mediterranean scrub","Deserts and xeric shrublands")


#create the x lab value for each row#
sudo_data_vector = c(rep(5:24,1),seq(25,100,5),seq(110,300,10))
#this is for the flooded grassland, just applied the tail deletion to that single region and then get the 2 equations left and get the model#
#PesudoData_df = read.csv("output_data/final_pesudo_data_repeat_kickout_flooded_grasslands.csv")
#for (eco in "Flooded grasslands")
# define a list to save the plot in the loop
plotList = list()
for (eco in ecoregionsOrdered)
{
    #par(mfrow=c(1:2))
    #scatter plot of the sudo biomass data frame#
    TempoaryDataframe = subset(PesudoData_df,PesudoData_df$Ecoregion_WWF==eco)
    print(nrow(TempoaryDataframe))
    #TempoaryDataframe                              <- TempoaryDataframe[TempoaryDataframe$Ecoregion_WWF==eco,]
    #plot(x=rep(sudo_data_vector,nrow(TempoaryDataframe)),y=unlist(TempoaryDataframe[,-1:-4]))
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    #train_data                   <- data.frame(y_value=(unlist(tt[,-1:-7])),x_value=(rep(sudo_data_vector,nrow(tt))))
    #please do not us unlist to transfer a dataframe to a vector use the "unlist" function!!!!!!!!!!!!!!!!!!!
    #train_data                   <- data.frame(y_value=as.vector(t(TempoaryDataframe[,-1:-6])),x_value=(rep(sudo_data_vector,nrow(TempoaryDataframe))))
    #plot(x=train_data$x_value,y=train_data$y_value)
    #test_df                      <- train_data_f[1:50,]
    #train_data_f                 <- test_df
    #try to train a small model for the data from the selected the ecoregion#
    #work out a two columns data frame with transformed x and y values#
    train_data = data.frame(y_value=log(as.vector(t(TempoaryDataframe[,-1:-7]))),x_value=log(rep(sudo_data_vector,nrow(TempoaryDataframe))))
    #remove the 0,NA and infinite vlaues in the train_data#
    train_data = na.omit(train_data)
    train_data = train_data[is.finite(rowSums(train_data)),]
    #remove the 0 and NA vlaues in the train_data#
    #train_data_f                     <- subset(train_data,train_data$y_value!=0)

    #delete the values which is smaller than 1, becaues it will produce some nagative and na values#
    #plot the scater plot and add the linear model line to the plot#
    #kick out the "/" in eocregion name#
    
    plotList[[which(ecoregionsOrdered==eco)]] = ggplot(train_data, aes(x = x_value,y = y_value)) + 
                                           geom_point(colour="cadetblue4",size=0.5) +
                                           geom_smooth(se = F, method = "lm",size =0.6,colour="coral3") +
                                           theme_classic() +
                                           theme(legend.position = 'none',panel.border = element_rect(color = "black",fill = NA,size = 1))+
                                           labs(x="Log(DBH)(cm)",y="Log(biomass)(kg)")+
                                           ylim(min(train_data$y_value),max(train_data$y_value))+
                                           xlim(min(train_data$x_value),max(train_data$x_value))+
                                           ggtitle(gsub("\\/", " ",eco))#delete the dash line in the boreal biome

    #save the summary of the lm models to text file#
    LM_Summary = summary(lm(train_data$y_value~train_data$x_value))
    sink(paste("plots_less/",gsub("\\/", " ",eco),"_primary_plot.txt",sep=""))
    print(LM_Summary)
    sink()
}
#STEP 9#
pdf("Plots/Figure_S1_Biome_Level_Allometric_equations_for_non_tropical_regions.pdf",height=10,width=10)
library(ggpubr)
ggarrange(plotlist=plotList,align="hv",nrow=3,ncol=3)
dev.off()

# 20180424_GFBI_Treelist_Data_wBiomes.csv column metadata:
# 'Unnamed_Column': 'Unnamed column; an index that can be removed',
# 'FID': 'Unique GFBI ID',
# 'LAT': 'Latitude in decimal degrees',
# 'LON': 'Longitude in decimal degrees',
# 'SPCD': 'Unique species code',
# 'DBH': 'Diameter at breast height of the tree (cm)',
# 'Year': 'The year the tree record was taken',
# 'WWF_Biome' : 'See the Biome Information below'

# WWF Biome Numbers:
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