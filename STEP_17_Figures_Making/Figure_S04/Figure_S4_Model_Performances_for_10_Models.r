
#*******************************************#
# Different data filtering comparison by grid seach 
library(data.table)
library(parallel)
library(ggplot2)
library(purrr)
library(dplyr)
library(GSIF)
library(rgdal)
library(ggpubr)
library(Metrics)
library(RColorBrewer)

modelNames = c("GS1_MaxScaler","GS1_MeanScaler","HM1","SD1","WK1")
rSquaredTableHD = data.frame()
# do this plot on my own computer
for (mn in modelNames)
{
    # get the position of the model name in the list 
    mnPosition = which(modelNames ==mn)
    # load all the 10-fold cross validation model train and predict tables
    tableNameList = do.call(paste0, expand.grid("Data/BiomassModelPerformance/",mn,"_10_Fold_CV_PredObs_Seed_",0:99,".csv"))
    # load each table into a a list which has biome level table as an element
    tableList = lapply(tableNameList,fread)

    cvRbindedTableFull = rbindlist(tableList) %>% dplyr::select(-c("system:index",".geo","CV_Fold")) %>%setNames(c("Human_Disturbance","Train","Predict"))
    for (hd in c(0,0.05,0.1,0.2,0.4,0.6))
    {
        cvRbindedTable = cvRbindedTableFull %>% dplyr::filter(Human_Disturbance>=hd)
        percentLeft = nrow(cvRbindedTable)/nrow(cvRbindedTableFull)*100
        # subset the table by human disturbance intensity
        humanFilteredR2 = 1 - (sum((cvRbindedTable$Train-cvRbindedTable$Predict)^2)/sum((cvRbindedTable$Train-mean(cvRbindedTable$Predict))^2))
        rSquaredLevel = data.frame(ModelName = mn,R2 = humanFilteredR2,Level= hd,Percent = percentLeft )
        print(rSquaredLevel)

        rSquaredTableHD = rbind(rSquaredTableHD,rSquaredLevel)
    }

}

write.csv(rSquaredTableHD,"Data/BiomassModelPerformance/rSquared_Table_10_fold_of_all_models_at_Different_Human_disturbance_levels.csv")

rSquaredTableHD$ModelName = as.factor(rSquaredTableHD$ModelName)
rSquaredTableFiltered = rSquaredTableHD %>% mutate(R_Squared = round(R2, 2))

modelPerformanceLevel = ggplot(rSquaredTableFiltered, aes(x = Level,y = R2,color = factor(ModelName, levels=c("GS1_MaxScaler","GS1_MeanScaler","HM1","SD1","WK1")))) + 
                                geom_point() +
                                geom_line(linewidth = 1) +
                                # geom_smooth(se = F, method = "loess",formula = y ~ x,size =0.6) +
                                coord_cartesian(ylim=c(0,1),xlim=c(0,0.6))+
                                 scale_color_manual(values=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3" ,"#A6D854"),labels = c("Upper canopy cover (GS)","Lower canopy cover (GS)","Harmonized (SD)","ESA-CCI (SD)","Walker et al. (SD)"))+
                                theme_bw() +
                                theme(legend.position = c(0.7, 0.3),
                                    legend.title=element_text(size=14),
                                    legend.text=element_text(size=12),
                                    panel.border = element_rect(fill=NA,color="black", 
                                    linewidth=1, linetype="solid"),
                                    axis.title = element_text(size = 16),
                                    axis.text = element_text(size = 14))+
                                labs(y="R squared",x="Human modification",colour = "Models")

sampleLevel = ggplot(rSquaredTableFiltered, aes(x = Level,y = Percent,color = factor(ModelName, levels=c("GS1_MaxScaler","GS1_MeanScaler","HM1","SD1","WK1")))) + 
                                geom_point() +
                                geom_line(linewidth = 1) +
                                coord_cartesian(ylim=c(0,100),xlim=c(0,0.6))+
                                 scale_color_manual(values=c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3" ,"#A6D854"))+
                                theme_bw() +
                                theme(legend.position = "none",panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),axis.title = element_text(size = 16),axis.text = element_text(size = 14))+
                                labs(y="Sample size (%)",x="Human modification",colour = "Models")


pdf(paste("WritingFolder/Plots/Figure_S4_Model_performance_at_different_human_modification_level.pdf",sep=""),height=10,width=5.2)
ggpubr::ggarrange(modelPerformanceLevel,sampleLevel,align="hv",widths = c(1),heights = c(1,1),labels = c("a", "b"),ncol = 1,font.label = list(size = 32, color = "black", face = "bold", family = NULL))
dev.off()

