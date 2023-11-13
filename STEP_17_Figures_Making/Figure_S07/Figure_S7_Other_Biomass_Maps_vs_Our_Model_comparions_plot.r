# Load the libraries 
library(raster)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(readr)
library(data.table)
library(dplyr)

###################################################################################################################
#### STEP 1 make the figure
###################################################################################################################
#define an empty list to save the plot panels
plotList = list()

randomSampledTableFull = fread("Code/STEP_17_Figures_Making/Figure_S07/Random_sample_Biomass_Maps_with_different_spatial_scales_Raw_data.csv")[,-1]
rSquaredTableFull = fread("Code/STEP_17_Figures_Making/Figure_S07/Random_sample_Biomass_Maps_with_different_spatial_scales_R_Squared.csv")[,-1]

rSquaredTableFull$NumbericScale = readr::parse_number(rSquaredTableFull$scale)
rSquaredTableFull$Models = as.factor(rSquaredTableFull$Models)
rSquaredTableFiltered = rSquaredTableFull %>% 
                        filter(Models %in% c("GB","HM","SD","WK","IPCC"))%>%
                        mutate(R_Squared = round(R_Squared, 2))

multipleModelCompare = ggplot(rSquaredTableFiltered, aes(x = NumbericScale,y = R_Squared,color = factor(Models, levels=c("GB","HM","SD","WK","IPCC")))) + 
                                geom_point() +
                                geom_smooth(se = F, method = "loess",formula = y ~ log(x),size =0.6) +
                                coord_cartesian(ylim=c(0,1))+
                                 scale_color_manual(values=c("#d55e00","#009e73","#56b4e9","#e69f00","#f0e442"))+
                                theme_bw() +
                                theme(legend.position = "none",panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),axis.title = element_text(size = 16),axis.text = element_text(size = 14))+
                                labs(y="R squared",x="Spatial scale (km)")+
                                ylim(0,1)

plotList[[1]] = multipleModelCompare


# load the two tables
randomSampledTable = fread("Code/STEP_17_Figures_Making/Figure_S07Figure_S6/Random_sample_Biomass_Maps_with_different_spatial_scales_Raw_data.csv")[,-1]
rSquaredTable = fread("Code/STEP_17_Figures_Making/Figure_S07/Random_sample_Biomass_Maps_with_different_spatial_scales_R_Squared.csv")[,-1]

# add the numberic column of the scale levels
rSquaredTable$NumbericScale = readr::parse_number(rSquaredTable$scale)
rSquaredTable$Models = as.factor(rSquaredTable$Models)
# subset the r square table at 1km res
rSquaredTableSub = rSquaredTable %>% 
                                 filter(scale == "1km") %>% 
                                 filter(Models %in% c("GB","HM","SD","WK","IPCC"))%>%
                                 mutate(R_Squared = round(R_Squared, 2))

# define the r square changes along scale gradients
barP = ggplot(rSquaredTableSub, aes(fill=factor(Models, levels=c("GB","HM","SD","WK","IPCC")), y=R_Squared, x=Models)) + 
            geom_bar(position="dodge", stat="identity",show.legend = T)+
            theme_classic()+
            coord_cartesian(ylim=c(0,1))+
            scale_fill_manual(values=c("#d55e00","#009e73","#56b4e9","#e69f00","#f0e442"),labels = c("GlobBiomass", "Harmonized","ESA-CCI","Walker et al.","IPCC"))+
            guides(fill = guide_legend(title = "Models",nrow=2))+
            theme(legend.position = c(0.5,0.9),
                legend.title=element_text(size=16),
                legend.text=element_text(size=14),
                axis.title = element_text(size = 16),
                axis.text = element_text(size = 14),
                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
            ylim(0,1)+
            scale_x_discrete(limits = c("WK","GB","HM","SD","IPCC"),labels=c("Waker et al.","GlobBiomass","Harmonized","ESA-CCI","IPCC"))+
            labs(y="R squared")+
            geom_text(aes(label= round(R_Squared,2)), position=position_dodge(width=0.9), vjust=-0.25,size = 3)
plotList[[2]] = barP


pdf(paste("Plots/Figure_S7_Other_Models_vs_Our_Model_Comparison_20230929.pdf",sep=""),height=6,width=12.32)
ggpubr::ggarrange(plotlist=plotList,align="hv",widths = c(1,1),heights = c(1),labels = c("a", "b"),font.label = list(size = 32, color = "black", face = "bold", family = NULL))
dev.off()


