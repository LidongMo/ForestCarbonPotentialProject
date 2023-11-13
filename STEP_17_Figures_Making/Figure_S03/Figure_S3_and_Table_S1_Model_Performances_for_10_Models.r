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


modelNames = c("GS1_MaxScaler","GS2_MaxScaler","GS1_MeanScaler","GS2_MeanScaler","HM1","HM2","SD1","SD2","WK1","WK2")
processList = list()
rSquaredTable = data.frame()
# Panel D
#*******************************************#
# do this plot on my own computer
for (mn in modelNames)
{
    # get the position of the model name in the list 
    mnPosition = which(modelNames ==mn)
    # load all the 10-fold cross validation model train and predict tables
    tableNameList = do.call(paste0, expand.grid("Data/BiomassModelPerformance/",mn,"_10_Fold_CV_PredObs_Seed_",0:99,".csv"))
    # load each table into a a list which has biome level table as an element
    tableList = lapply(tableNameList,fread)
    if (mn %in% c("GS1_MaxScaler","GS1_MeanScaler","HM1","SD1","WK1"))
    {
        cvRbindedTable = rbindlist(tableList) %>% dplyr::select(-c("system:index",".geo","CV_Fold","Human_Disturbance")) %>%setNames(c("Train","Predict"))
    }else
    {
        # since there is no humand disturbance been used in the model,therefore we dong have the humand distubance in the file to exclude
        cvRbindedTable = rbindlist(tableList) %>% dplyr::select(-c("system:index",".geo","CV_Fold")) %>%setNames(c("Train","Predict"))
    }


    # calculate the r square for all the models
    rSquared =  1 - (sum((cvRbindedTable$Train-cvRbindedTable$Predict)^2)/sum((cvRbindedTable$Train-mean(cvRbindedTable$Predict))^2))
    print(rSquared) 
    if (mn %in% c("GS1_MaxScaler","GS1_MeanScaler","GS2_MaxScaler","GS2_MeanScaler"))
    {
        limMax = 8
        anoteX = 2
        anoteY = 6

    }else
    {
        limMax = 225
        anoteX = 50
        anoteY = 180
    }
    if (mnPosition ==1)
    {
        p1 = ggplot(data=cvRbindedTable,aes(x=Predict, y=Train))+
                stat_bin2d(data=cvRbindedTable, aes(x=Predict, y=Train),bins=200)+
                scale_fill_gradientn(colours=c("grey80","grey30","grey10","grey10"))+
                #  geom_point( ) +
                geom_smooth(method="lm",size=1,se= F,color = "darkblue")+
                theme_classic() +
                annotate(geom="text", x=anoteX, y=anoteY, label=paste("R2 = ",round(rSquared,2),sep=""),color="black",size = 6)+
                labs(x = "Prediction", y = "Observation") +
                coord_cartesian(xlim=c(0,limMax), ylim=c(0,limMax))+
                geom_abline(slope=1, intercept=0,na.rm = FALSE, show.legend = NA,  linetype="dashed",size = 0.3)+
                theme(legend.position = c(0.8, 0.3),
                    panel.border = element_rect(color = "black",fill = NA,size = 1),
                    legend.key.size = unit(0.4, 'cm'),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14))
    }else
    {
        p1 = ggplot(data=cvRbindedTable,aes(x=Predict, y=Train))+
                stat_bin2d(data=cvRbindedTable, aes(x=Predict, y=Train),bins=200)+
                scale_fill_gradientn(colours=c("grey80","grey30","grey10","grey10"))+
                #  geom_point( ) +
                geom_smooth(method="lm",size=1,se= F,color = "darkblue")+
                theme_classic() +
                annotate(geom="text", x=anoteX, y=anoteY, label=paste("R2 = ",round(rSquared,2),sep=""),color="black",size = 6)+
                labs(x = "Prediction", y = "Observation") +
                coord_cartesian(xlim=c(0,limMax), ylim=c(0,limMax))+
                geom_abline(slope=1, intercept=0,na.rm = FALSE, show.legend = NA,  linetype="dashed",size = 0.3)+
                theme(legend.position = "none",
                    panel.border = element_rect(color = "black",fill = NA,size = 1),
                    legend.key.size = unit(0.4, 'cm'),
                    axis.title = element_text(size = 16),
                    axis.text = element_text(size = 14))
    }
    rSquaredModel = data.frame(ModelName = mn,R2 = rSquared)
    rSquaredTable = rbind(rSquaredTable,rSquaredModel)
    
    processList[[mnPosition]] = p1
}


write.csv(rSquaredTable,"Data/BiomassModelPerformance/rSquared_Table_10_fold_CV_based_of_all_models.csv")



pdf("WritingFolder/Plots/Fig_S3_Panel_A_B_C_D_E_F_G_H_I_J_K_L_Model_Performance_20230427.pdf",width =6.5, height=15.5)

# library("gridExtra")
ggpubr::ggarrange(plotlist = processList,ncol = 2, nrow = 5,align="hv",widths = c(1,1,1,1,1,1),heights = c(1,1,1,1,1),labels = c("a", "b", "c","d","e","f","g","h","i","j"),font.label = list(size = 25))

dev.off()



rSquaredTable = fread("Data/BiomassModelPerformance/rSquared_Table_10_fold_CV_based_of_all_models.csv")[,-1]
modelNames = c("GS_Max1","GS_Max2","GS_Mean1","GS_Mean2","HM1","HM2","SD1","SD2","WK1","WK2")
LOOCV_R2_Table = data.frame()
for (mnn in modelNames)
{
    tableNameList = paste("Data/BiomassModelPerformance/LOOCV_Results/Model_",mnn,"_Leav_One_Out_Cross_Validation_rSquared_",0:99,".csv",sep="")
    # read the tables into a list
    rSquaredVector = lapply(tableNameList,fread) %>% rbindlist() %>% dplyr::select(R2_val) %>% unlist()
    perVal = data.frame(LOOCV_R2 = mean(rSquaredVector))
    LOOCV_R2_Table = rbind(LOOCV_R2_Table,perVal)
}
full_rSquared_Table = cbind(rSquaredTable,LOOCV_R2_Table)

write.csv(full_rSquared_Table,"Code/STEP_17_Figures_Making/Table_S01/Table_S1_rSquared_Table_10_fold_CV_and_LOOCV_based_of_all_models_20230427.csv")





