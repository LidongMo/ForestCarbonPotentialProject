library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(data.table)
# set working directory 

##################################################################
# STEP 1 SD 
##################################################################

mapNames = c("GS_Max1","GS_Max2","GS_Mean1","GS_Mean2","SD1","SD2","HM1","HM2","WK1","WK2")

# define a function to add the model name/Type to each data frame 
addFunc = function(mp)
{
    # read the table for each model
    perModelTable =fread(paste("Code/STEP_17_Figures_Making/Figure_02/",mp,"_Full_TGB_carbon_density_latitude_agregated_table_with_SD.csv",sep=""))[,-1] %>% mutate(Type = mp)
    # return the table
    return(perModelTable)
}
# use lapply to merge all the tbales by row
tableList = lapply(mapNames,addFunc)
aggregatedTablePotential = rbindlist(tableList) 
aggregatedTablePotential = aggregatedTablePotential %>%mutate_at(vars(c('lowCI')),~ifelse(lowCI <=0, 0, .))


# display.brewer.pal(n, name)
p1 = aggregatedTablePotential %>%ggplot(aes(x = LatRound, y= mean, color=Type, group=Type)) + 
    geom_ribbon(aes(ymin=lowCI, ymax=highCI, fill=Type),alpha=0.1,color=NA)+
    scale_fill_manual(values = c("GS_Max1" = "#009e73","GS_Max2" = "#009e73","GS_Mean1" = "#cc79a7","GS_Mean2" = "#cc79a7","HM1"="#56b4e9","HM2"= "#56b4e9","SD1"="#e69f00","SD2"="#e69f00","WK1"="#d55e00","WK2"= "#d55e00"))+
    # geom_hline(yintercept=0)+
    geom_line(aes(linetype=Type),lwd = 0.35,)+
    scale_color_manual(values = c("#009e73","#009e73","#cc79a7","#cc79a7","#56b4e9","#56b4e9","#e69f00","#e69f00","#d55e00","#d55e00"))+
    scale_linetype_manual(values=c(rep(c("solid","dashed"), 5)))+
    scale_x_continuous(expand = c(0,0),limits = c(-60, 90), breaks = c(-60,-40,-20,0,20,40,60,80,90))+
    # scale_alpha_discrete(range = c(0.2, 0.7))+
    ylab("") +
    theme_classic()+
    coord_flip() + #ylim = c(-0.4, 0.4),xlim=c(27,75)
    ylim(0,240)+
    xlab("Latitude")+
    guides( col = guide_legend(ncol=2),linetype = guide_legend())+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position=c(0.8,0.8))


mapNames = c("GS_Max1","GS_Max2","GS_Mean1","GS_Mean2","SD1","SD2","HM1","HM2","WK1","WK2")

# define a function to add the model name/Type to each data frame 
addFunc = function(mp)
{
    # read the table for each model
    perModelTable =fread(paste("Code/STEP_17_Figures_Making/Figure_02/",mp,"_Net_TGB_carbon_density_latitude_agregated_table_with_SD.csv",sep=""))[,-1] %>% mutate(Type = mp)
    # return the table
    return(perModelTable)
}
# use lapply to merge all the tbales by row
tableList = lapply(mapNames,addFunc)
aggregatedTableDiff = rbindlist(tableList) 
aggregatedTableDiff = aggregatedTableDiff %>%mutate_at(vars(c('lowCI')),~ifelse(lowCI <=0, 0, .))



p2 = aggregatedTableDiff %>%ggplot(aes(x = LatRound, y= mean, color=Type, group=Type)) + 
        geom_ribbon(aes(ymin=lowCI, ymax=highCI, fill=Type),alpha=0.1,color=NA)+
    scale_color_manual(values = c("#009e73","#009e73","#cc79a7","#cc79a7","#56b4e9","#56b4e9","#e69f00","#e69f00","#d55e00","#d55e00"))+
        geom_line(aes(linetype=Type),lwd = 0.35)+
            scale_linetype_manual(values=c(rep(c("solid","dashed"), 5)))+
    scale_fill_manual(values = c("GS_Max1" = "#009e73","GS_Max2" = "#009e73","GS_Mean1" = "#cc79a7","GS_Mean2" = "#cc79a7","HM1"="#56b4e9","HM2"= "#56b4e9","SD1"="#e69f00","SD2"="#e69f00","WK1"="#d55e00","WK2"= "#d55e00"))+
    # scale_alpha_discrete(range = c(0.2, 0.7))+
    scale_x_continuous(expand = c(0,0),limits = c(-60, 90), breaks = c(-60,-40,-20,0,20,40,60,80,90))+
    ylab("") +
    theme_classic()+
    coord_flip() + 
    ylim(0,150)+
    xlab("Latitude")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.position="none")


# read the table for each model
netModelTable =fread(paste("Code/STEP_17_Figures_Making/Figure_02/Net_diff_potential_carbon_density_latitude_agregated_table_with_SD.csv",sep=""))[,-1] %>% mutate(Type ='Net')

p3 = netModelTable %>%ggplot(aes(x = LatRound, y= mean, color=Type)) + 
        geom_ribbon(aes(ymin=lowCI, ymax=highCI, fill=Type),alpha=0.2,color=NA)+
    scale_color_manual(values = c("gray25"))+
        geom_line(aes(linetype=Type),lwd = 0.35)+
            scale_linetype_manual(values=c(rep(c("solid","dashed"), 5)))+
    scale_fill_manual(values = c("gray25"))+
    # scale_alpha_discrete(range = c(0.2, 0.7))+
    scale_x_continuous(expand = c(0,0),limits = c(-60, 90), breaks = c(-60,-40,-20,0,20,40,60,80,90))+
    ylab("") +
    geom_hline(yintercept = 0, linetype="dotted", color = "gray20", size=0.5)+
    theme_classic()+
    coord_flip() + 
    ylim(-50,50)+
    xlab("Latitude")+
    guides( col = guide_legend(ncol=1),linetype = guide_legend())+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),legend.position=c(0.8,0.8))




pdf("Plots/Fig_2_Panel_b_d_f_GS_SD_Models_Comparison_of_Potential_Diff_biomass_estimation_SD_20230427.pdf",width =7, height=33)
library("gridExtra")
ggarrange(p1,p2,p3,ncol = 1, nrow = 3,align="hv",widths = c(1),heights = c(3,3) )

dev.off() 

