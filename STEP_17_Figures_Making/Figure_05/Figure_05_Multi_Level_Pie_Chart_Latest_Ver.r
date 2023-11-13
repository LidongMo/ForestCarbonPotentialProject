library(gridExtra)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(dplyr)
library(abind)
library(ggpattern)
# define the data frame 
# load the data frame
# get the average of the TGB estimates 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# panel a was created by anther code



modelNames = c("GS1_MeanScaler","GS1_MaxScaler","GS2_MeanScaler","GS2_MaxScaler","SD1","SD2","WK1","WK2","HM1","HM2")
tableGetFunc = function(mn)
{
    # read the tables of carbon stock partition
    landuseTypePotentialTB = read.csv(paste("Data/BiomeLevelStatistics/StatisticsForModels/",mn,"_Biome_Level_Statistics_with_Litter.csv",sep=""))[,-1] %>% 
    mutate(ConservationPotential = PresentPotential- Present) %>%data.table::first(14) 

    return(landuseTypePotentialTB)
}

outputList = lapply(modelNames,tableGetFunc)
# combine the list of data frame into a 3D array
allMatrix = abind::abind(outputList, along=3)
# oget the mean for each element
outputTable = apply(allMatrix, c(1,2), mean)

# write this into local folder
write.csv(outputTable,"Code/STEP_17_Figures_Making/Figure_05/Mean_of_the_biome_level_statistics_with_landsue_type_with_Litter.csv")

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

# write this into local folder
write.csv(outputMean,"FCode/STEP_17_Figures_Making/Figure_05/Mean_of_the_biome_level_statistics_with_landsue_type.csv")


# define the data frame 1
# existing carbon stock, conservation potential and reforestation potential in the global forest
df1 = data.frame(
  group = c("A_Low-human pressure land", "B_Rangeland","C_Pasture","D_Cropland","E_Urban","F_Plantation","G_Conservation potential"),
  value = c(53.3,33.7,43.3,56.8,1.7,10.5,128.4)) #10.5 for planted 138.9 42.4% for conservation.

panelA = ggplot(df1, aes(x="", y=value, fill=group))+
                # geom_bar(width = 1, stat = "identity")+
                geom_col_pattern(aes(pattern = group),
                   pattern_color = "gray30",
                   pattern_fill = "gray30",
                   pattern_density = 0.05, 
                   pattern_angle = 30,
                   pattern_spacing =0.02,
                   pattern_size=0.5,
                   pattern_linetype=1)+
                coord_polar("y", start=0)+
                scale_fill_manual(values=c("#FEE8C8","#FDD49E","#FDBB84","#FC8D59", "#EF6548","#2B8CBE","#2B8CBE")) + #"#D95F02"
                geom_text(size=2.5,aes(label = paste(c("A_Low-human pressure land", "B_Rangeland","C_Pasture","D_Cropland","E_Urban","F_Plantation","G_Conservation potential")," ",value," Gt C(",round(value/sum(value) * 100, 1), "%)",sep=""), x = 1.3),position = position_stack(vjust = 0.5))+
                scale_pattern_manual(values = c("none","none","none","none", "none","wave","none"))+
                # scale_pattern_spacing_discrete(range = c(0., 0.1))+
                # ggpattern::scale_pattern_discrete(choices = c("none","none","none","none", "none","stripe","none"))+
                theme(legend.position = "none")


# conservation potential of different forest types
df2 = data.frame(
  group = c("A_Tropical", "B_Temperate", "C_Boreal","D_Dryland"), #
  value = c(87.1,31.7,16.7,3.5 ))

panelB = ggplot(df2, aes(x="", y=value, fill=group))+
                geom_bar(width = 1, stat = "identity")+
                coord_polar("y", start=0)+
                scale_fill_manual(values=c( "#2171B5","#4292C6","#6BAED6","#9ECAE1"))+
                geom_text(size=2.5,aes(label = paste(c("A_Tropical", "B_Temperate", "C_Boreal","D_Dryland")," ",value," Gt C(",round(value/sum(value) * 100, 1), "%)",sep=""), x = 1.3),position = position_stack(vjust = 0.5))+
                theme(legend.position = "none")

df5 = data.frame(
  group = c("A_Brazil","B_United_States","C_Russia","D_China","E_Australia","F_Canada","G_India","H_Indonesia","I_DR_Congo","J_Others"),
  value = c(39.24,28.13,22.72,21.60,15.54,12.82,12.04,8.04,7.27,160.40))

panelE = ggplot(df5, aes(x="", y=value, fill=group))+
                geom_bar(width = 1, stat = "identity")+
                coord_polar("y", start=0)+
                scale_fill_manual(values=c("#00441B","#006D2C", "#238B45" ,"#41AB5D" ,"#74C476" ,"#A1D99B", "#C7E9C0", "#E5F5E0" ,"#F7FCF5","#b5b5b5"))+
                geom_text(size=2.5,aes(label = paste(c("A_Brazil","B_United_States","C_Russia","D_China","E_Australia","F_Canada","G_India","H_Indonesia","I_DR_Congo","I_Others")," ",value," Gt C(",round(value/sum(value) * 100, 1), "%)",sep=""), x = 1.3),position = position_stack(vjust = 0.5))+
                theme(legend.position = "none")


# 
# full difference potential in different parts
df7 = data.frame(
  group = c("A_Aboveground biomass", "B_Root_Biomass","C_Dead wood and litter", "D_Soil"),
  value = c(170.7,46.0,61.8,49.3))

panelG = ggplot(df7, aes(x="", y=value, fill=group))+
                geom_bar(width = 1, stat = "identity")+
                coord_polar("y", start=0)+
                scale_fill_manual(values=c("#2B7700","#559233","#B9AE00","#8B5900"))+
                geom_text(size=2.5,aes(label = paste(c("A_Aboveground biomass", "B_Root_Biomass","C_Dead wood and litter", "D_Soil")," ",value," Gt C(",round(value/sum(value) * 100, 1), "%)",sep=""), x = 1.3),position = position_stack(vjust = 0.5))+
                theme(legend.position = "none")





pdf("Plots/Fig_05_Panel_A_B_C_D_E_202304271.pdf",width =12, height=12)

 # ggarrange(panelA,panelB,panelC,panelD,panelE,ncol = 3, nrow = 2)

grid.arrange(panelA,panelB,panelE,panelG,ncol = 2, nrow = 2,
             layout_matrix = rbind(c(1,2),
                                   c(4,3)))
dev.off()

