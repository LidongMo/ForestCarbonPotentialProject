
library(data.table)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(patchwork)

# plot data frame 


dataTable = read.csv("Data/BiomeLevelStatistics/TGB_Statistics_summary_with_mean_and_range.csv")[,-1] %>% mutate(Type = c("Ground1","Ground1","Ground2","Ground2","Satellite1","Satellite2","Satellite1","Satellite2","Satellite1","Satellite2")) 

names(dataTable) = c("Present","Fit","Lower","Upper","ModelName","Model","Type")
pd = position_dodge(.1) 
dataTable$Type = factor(dataTable$Type, levels = c("Ground1","Satellite1","Ground2","Satellite2"))

dataTableRep1 = dataTable %>% mutate(Fit = Fit-Present,Lower =NA,Upper = NA,Model ='Model1') %>% mutate(facet = c("GS1","GS1","GS2","GS2","RM1","RM2","RM1","RM2","RM1","RM2")) %>% mutate(NewType = c("T2","T1","T2","T1","T2","T2","T3","T3","T1","T1")) #%>% mutate(ModelName = paste(ModelName,"_1",sep=""))
dataTableRep2 = dataTable %>% mutate(Fit =Present,Model ='Model2')  %>% mutate(facet = c("GS1","GS1","GS2","GS2","RM1","RM2","RM1","RM2","RM1","RM2")) %>%mutate(NewType = c("T2","T1","T2","T1","T2","T2","T3","T3","T1","T1")) #%>% mutate(ModelName = paste(ModelName,"_2",sep=""))
bindedDataTable = rbind(dataTableRep2,dataTableRep1)


p1 = ggplot(data=bindedDataTable, aes(x=NewType, y=Fit, fill=factor(ModelName),alpha=factor(Model))) +
  # annotate("rect", fill = "gray40", alpha = 0.3, xmin = -Inf, xmax = Inf,ymin = 343.5, ymax = 468.1)+ 
  geom_bar(stat="identity",position ="stack",colour="white")+
  geom_errorbar(aes(ymin = Lower,ymax  = Upper),width = 0.5, linewidth  = 0.6,position =position_dodge(width=0.1,preserve = "single")) +
    coord_cartesian(ylim=c(0,1100))+
  facet_grid(~ facet,scales="fixed")+
  # scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("#009e73","#cc79a7","#009e73","#cc79a7","#56b4e9","#56b4e9","#e69f00","#e69f00","#d55e00","#d55e00"))+
  scale_alpha_discrete(range=c(0.7, 1))+
  scale_y_continuous(expand = c(0,0),limits = c(0, 1100), breaks = 5*c(0,40,80,120,160,200))+
  theme_classic()+ 
  labs(y="Carbon stock (Gt C) in living tree biomass",x= "Estimations")+
  guides(fill = guide_legend(title = "Models",nrow=2),alpha = "none")+
  theme(legend.position = c(0.5, 0.95),
    legend.title=element_text(size=16),
    legend.text=element_text(size=14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    strip.background = element_blank(),
    strip.text.x = element_blank())

# p2

dataTable2 = data.frame(Number = 1:20,
  biomass = c(1089.8,980.0,977.6,971.0,956.0,907.6,900.0,877.5,857.0,924.0,850.0,772.0,771.0,737.0,695.0,610.0,641.0,923.0,916.0,795.5),
  Type = "Biomass",
  ModelType = c("Ensemble","Inventory","Ensemble","Inventory","Inventory","Inventory","Inventory","Inventory","Mechanistic model","Inventory","Mechanistic model","Mechanistic model","Ensemble","Inventory","Mechanistic model","Ensemble","Mechanistic model","Mechanistic model","Ensemble","Data-driven model"))
set.seed(1)
p2 = ggplot(dataTable2, aes(x=Type, y=biomass)) + 
  geom_violin(trim=FALSE,fill= "gray75",color = "gray70")+ #mean is 857.5 Gt C
    stat_summary(fun.data=mean_sdl,fun.args = list(mult = 1),geom="pointrange", color="gray10")+
  geom_jitter(aes(colour = ModelType,shape=ModelType), position=position_jitter(0.2),size=3)+
  scale_colour_manual(values = c("Ensemble"="#FDC086", "Inventory"="#386CB0", "Mechanistic model" = "#F0027F" ,"Data-driven model" = "#BF5B17"))+
  scale_shape_manual(values=c(15,16,17,18))+
  theme_classic()+
    scale_y_continuous(expand = c(0,0),limits = c(0, 1100),breaks = 5*c(0,40,80,120,160,200))+
  theme(legend.position = c(0.5, 0.2),panel.border = element_rect(colour = "black", fill=NA, size=1.5),axis.title = element_text(size = 16),axis.text = element_text(size = 14))+
  # ylim(0,1100)+
  coord_cartesian(ylim = c(0, 1100))+
  geom_hline(yintercept = 442.5,size=1,colour="black",linetype = c("dashed"))+annotate(geom="text", label=443, x=1.5, y=443, vjust=-1,size=5,hjust=-5)
  # geom_hline(yintercept = 438.1,linetype=c("dashed"))+annotate(geom="text", label=438.1, x=1.5, y=438.1, vjust=1.5,size=5,hjust=-5)+



dataTable3 = read.csv("Data/BiomeLevelStatistics/TGB_Statistics_summary_with_mean_and_range.csv")[,-1] %>% mutate(Type = c("Ground1","Ground1","Ground2","Ground2","Satellite1","Satellite2","Satellite1","Satellite2","Satellite1","Satellite2")) %>% mutate(AbsolutePotential = AbsolutePotential - Present,potTGB_Lower = potTGB_Lower-Present,potTGB_Upper = potTGB_Upper-Present) 

names(dataTable3) = c("Present","Fit","Lower","Upper","ModelName","Model","Type")
dataTable3 = dataTable3%>% mutate(facet = c("GS1","GS1","GS2","GS2","RM1","RM2","RM1","RM2","RM1","RM2")) %>%mutate(NewType = c("T2","T1","T2","T1","T2","T2","T3","T3","T1","T1")) 

pd = position_dodge(.3) 
dataTable3$Type = factor(dataTable3$Type, levels = c("Ground1","Satellite1","Ground2","Satellite2"))



p3 = ggplot(data=dataTable3, aes(x=NewType, y=Fit, fill=factor(ModelName),alpha=factor(Model))) +
  # annotate("rect", fill = "gray40", alpha = 0.3, xmin = -Inf, xmax = Inf,ymin = 343.5, ymax = 468.1)+ 
  geom_bar(stat="identity",position ="dodge")+
  geom_errorbar(aes(ymin = Lower,ymax  = Upper),width = 0.27, linewidth = 0.6,position =position_dodge(width=0.1,preserve = "single")) +
    coord_cartesian(ylim=c(0,550))+
  facet_grid(~ facet,scales="fixed")+
  # scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values = c("#009e73","#cc79a7","#009e73","#cc79a7","#56b4e9","#56b4e9","#e69f00","#e69f00","#d55e00","#d55e00"))+
  scale_alpha_discrete(range=c(0.7, 0.7))+
  scale_y_continuous(expand = c(0,0),limits = c(0, 550), breaks = 5*c(0,50,100))+
  theme_classic()+ 
  labs(y="Carbon stock (Gt C) in living tree biomass",x= "Estimations")+
  guides(fill = guide_legend(title = "Models",nrow=2),alpha = "none")+
  theme(legend.position = "none",axis.title = element_text(size = 16),axis.text = element_text(size = 14),strip.background = element_blank(),
  strip.text.x = element_blank())


dataTable4 = data.frame(Number = 1:9,
  biomass = c(319.0,340.0,150.0,183.8,378.6,242.0,254.0,466.0,354.4),
  Type = 'Biomass',
  ModelType = c("Inventory","Inventory","Mechanistic model","Mechanistic model","Inventory","Inventory","Mechanistic model","Ensemble","Data-driven model"))


p4 = ggplot(dataTable4, aes(x=Type, y=biomass)) + 
  geom_violin(trim=FALSE,fill= "gray75",color = "gray70")+ #mean is 298.6 Gt C
    stat_summary(fun.data=mean_sdl,fun.args = list(mult = 1),
       geom="pointrange", color="gray10")+
  geom_jitter(aes(colour = ModelType,shape=ModelType), position=position_jitter(0.2),size=3)+
  scale_colour_manual(values = c("Ensemble"="#FDC086", "Inventory"="#386CB0", "Mechanistic model" = "#F0027F" ,"Data-driven model" = "#BF5B17"))+
  scale_x_discrete(labels=c('Literature'))+
  scale_shape_manual(values=c(15,16,17,18))+
  theme_classic()+
    scale_y_continuous(expand = c(0,0),limits = c(0, 550), breaks = 5*c(0,50,100))+
  theme(legend.position = "none",panel.border = element_rect(colour = "black", fill=NA, linewidth=1.5),axis.text = element_text(size = 14),axis.title = element_text(size = 16))+
  coord_cartesian(ylim = c(0, 550))
  
# ,legend.title=element_text(size=5),legend.text=element_text(size=3)


pdf("Plots/Figure_3_Model_Conservness_comparison_20230929.pdf",width =12, height=12)
p1 + p2 + p3 + p4 + 
  plot_layout(widths = c(4, 1), heights = c(2,1))+plot_annotation(tag_levels = c("a"))&theme(plot.tag = element_text(face = 'bold',size=20))
# ggarrange(p1,p2,p3,p4,ncol = 2, nrow = 2,align="v",widths= c(0.8,0.2),heights=c(0.7,0.3))#
dev.off()   


