library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(data.table)
setwd("~/Desktop/BIOMASS/")

##################################################################
# STEP 1 SD 
##################################################################

mapNames = c("GS_Max","GS_Mean","SD","HM","WK")

# define a function to add the model name/Type to each data frame 
addFunc = function(mp)
{
    # read the table for each model
    perModelTable =fread(paste("Code/STEP_17_Figures_Making/Figure_S06/",mp,"_TGB_Present_carbon_density_latitude_agregated_table_with_SD.csv",sep=""))[,-1] %>% mutate(Type = mp)
    # return the table
    return(perModelTable)
}
# use lapply to merge all the tbales by row
tableList = lapply(mapNames,addFunc)
aggregatedTablePresentTGB = rbindlist(tableList) 
aggregatedTablePresentTGB = aggregatedTablePresentTGB %>%mutate_at(vars(c('lowCI')),~ifelse(lowCI <=0, 0, .))


# display.brewer.pal(n, name)
p1 = aggregatedTablePresentTGB %>%ggplot(aes(x = LatRound, y= mean, color=Type, group=Type)) + 
    geom_ribbon(aes(ymin=lowCI, ymax=highCI, fill=Type),alpha=0.1,color=NA)+
    scale_fill_manual(values = c("GS_Max" = "#009e73","GS_Mean" = "#cc79a7","HM"="#56b4e9","SD"="#e69f00","WK"="#d55e00"))+
    # geom_hline(yintercept=0)+
    geom_line(aes(linetype=Type),lwd = 0.5)+
    scale_color_manual(values = c("#009e73","#cc79a7","#56b4e9","#e69f00","#d55e00"))+
    scale_linetype_manual(values=c(rep(c("solid"), 5)))+
     scale_y_continuous(expand = c(0,0),limits = c(0, 250), breaks = 5*c(0,20,40))+
    coord_cartesian(ylim=c(0,250))+
    theme_classic()+
    xlim(-56,80)+
    # ylim(0,250)+
    labs(y="Living tree biomass density (tC/ha)",x= "Latitude")+
    guides( col = guide_legend(ncol=2),linetype = guide_legend(ncol=2))+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1),
          legend.position=c(0.8,0.8),
          legend.title=element_blank(),
          legend.text=element_text(size=14),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14))


pdf("Plots/Fig_S6_TGB_Comparison_of_Present_biomass_estimation_along_latitude_SD_20230929.pdf",width =13, height=5)
library("gridExtra")
ggarrange(p1,ncol = 1, nrow = 1,align="hv",widths = c(1),heights = c(3,3) )

dev.off() 

