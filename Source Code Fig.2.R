######## Source code Figure2

#load package
`%notin%` <- Negate(`%in%`)
library(dplyr)
library(ggpubr)
require("ggpubr")
require(ggrepel)
library(readr)
library(tidyverse)
library(readr)

##setwd to "fig2table1_plottingfiles"

###P53
#2f p53 dependency

p53_dep_vol <- read_csv("p53_dep_vol.csv")

fig2f_p53depgg <- ggplot(p53_dep_vol, aes(x = Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(color = Correlation, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% c(top_n(p53_dep_vol,6)$Feature, c("MDM2"))),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = "Anti-Dependency <--- Correlation ---> Co-Dependency", 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Times",size = 14, face = "bold"),
        axis.title= element_text(color='black',family="Times",size = 14, face = "bold"),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=10,face="bold"))

ggsave(fig2f_p53depgg, filename = "Fig2f_dep.svg", width = 6, height = 5)


#2g p53 drug resistance
p53_drug_vol <- read_csv("p53_drug_mut.csv")
fig2g_p53druggg <- ggplot(p53_drug_vol, aes(x = Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(color = Correlation, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% top_n(p53_drug_vol,1)$Feature),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = "Resistance <--- Correlation --->Sensitivity ", 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Times",size = 14, face = "bold"),
        axis.title= element_text(color='black',family="Times",size = 14, face = "bold"),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=10,face="bold"))

ggsave(fig2g_p53druggg, filename = "Fig2g_p53drug.svg", width = 6, height = 5)


#2a p53 mutation box plot

fig2a_p53mut <- ggplot(p53_p53mut, aes(y=p53_p53mut$`custom data  11_14.csv`, x=as.factor(p53_p53mut$`TP53 Mutation (one hot encoding) Mutation Internal 20Q4`),fill =as.factor(p53_p53mut$`TP53 Mutation (one hot encoding) Mutation Internal 20Q4`))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("TP53 hotspot mutation")+
  ylab("p53_overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8 ,size = 6)
ggsave(fig2a_p53mut, filename = "Fig2a_p53mut.svg", width = 4, height = 5)

#2j PTEN dependency
PTEN_dep <- read_csv("PTEN_Dep_voc.csv")
PTEN_drug <- read_csv("PTEN_Drug_vol.csv")
fig2j_PTENdep <- ggplot(PTEN_dep, aes(x = Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(size = abs(Correlation), color = Correlation, alpha=0.75)) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% c("RICTOR", "MARK2", "TRPM5", "PIK3CA")),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(0.8, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = "Anti-Dependency <--- Correlation ---> Co-Dependency", 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Times",size = 14, face = "bold"),
        axis.title= element_text(color='black',family="Times",size = 14, face = "bold"),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=10,face="bold"))

ggsave(fig2j_PTENdep, filename = "Fig2j_PTENdep.svg", width = 6, height = 5)


#2k PTEN Drug
fig2k_PTENdrug <- ggplot(PTEN_drug , aes(x = Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(size = abs(Correlation), color = Correlation, alpha=0.75)) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% top_n(PTEN_drug,4)$Feature),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(0.7, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = "Anti-Dependency <--- Correlation ---> Co-Dependency", 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Times",size = 14, face = "bold"),
        axis.title= element_text(color='black',family="Times",size = 14, face = "bold"),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=10,face="bold"))

ggsave(fig2k_PTENdrug, filename = "Fig2k_PTENdrug.svg", width = 6, height = 5)


#Fig2d mutations
PTEN_PTENmut <- read_csv("PTEN_PTENmut.csv")
Fig2d_PTENmut <- ggplot(PTEN_PTENmut, aes(y=PTEN_PTENmut$`custom data  11_9.csv`, x=as.factor(PTEN_PTENmut$`PTEN Mutation (one hot encoding) Mutation Internal 20Q4` ),fill=as.factor(PTEN_PTENmut$`PTEN Mutation (one hot encoding) Mutation Internal 20Q4` ))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("PTEN hotspot mutation")+
  ylab("PTEN_overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 6) 

ggsave(Fig2b_PTENmut, filename = "Fig2b_PTENmut.svg", width = 4, height = 5)


PTEN_PI3Kmut <- read_csv("PTEN_PI3Kmut.csv")

Fig2e_PI3Kmut <- ggplot(PTEN_PI3Kmut, aes(y=PTEN_PTENmut$`custom data  11_9.csv`, x=as.factor(PTEN_PI3Kmut$`PIK3CA Mutation (one hot encoding) Mutation Internal 20Q4` ),fill=as.factor(PTEN_PI3Kmut$`PIK3CA Mutation (one hot encoding) Mutation Internal 20Q4` ))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("PIK3CA hotspot mutation")+
  ylab("PTEN_overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 6) 

ggsave(Fig2e_PI3Kmut, filename = "Fig2e_PI3Kmut.svg", width = 4, height = 5)


# Fig2c p16

# dep
p16_CRISPR_vol <- read_csv("p16_CRISPR_vol.csv")

fig2h_p16dep <- ggplot(p16_CRISPR_vol, aes(x = Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(size = abs(Correlation), color = Correlation, alpha=0.75)) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% c(top_n(p16_CRISPR_vol,4)$Feature, c("CDK2"))),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = "Anti-Dependency <--- Correlation ---> Co-Dependency", 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Times",size = 14, face = "bold"),
        axis.title= element_text(color='black',family="Times",size = 14, face = "bold"),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=10,face="bold"))
ggsave(fig2h_p16dep, filename = "Fig2h_p16dep.svg", width = 6, height = 5)


#p16 RBmut


Fig2b_Rbmut <- ggplot(p16_RBmut, aes(y=p16_RBmut$`custom data  11_10.csv`, x=as.factor(p16_RBmut$`RB1 Mutation (one hot encoding) Mutation Internal 20Q4`),fill =as.factor(p16_RBmut$`RB1 Mutation (one hot encoding) Mutation Internal 20Q4`))) +  
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("RB1 hotspot mutation")+
  ylab("P16_overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 6) 

ggsave(Fig2b_Rbmut, filename = "Fig2b_RBmut.svg", width = 4, height = 5)


#p16 expression
p16_expression_vol<- read_csv("p16_expression_vol.csv",col_types = cols(Correlation = col_double(),`-log10(PValue)` = col_double()))

fig2i_p16expression <- ggplot(p16_expression_vol, aes(x = -Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(size = abs(Correlation), color = Correlation, alpha=0.75)) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% c("RB1","CCNE1","CCNE2","E2F3","CDKN2A")),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = "Sensitivity <--- Correlation ---> Resistance", 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',family="Times",size = 14, face = "bold"),
        axis.title= element_text(color='black',family="Times",size = 14, face = "bold"),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=10,face="bold"))
ggsave(fig2i_p16expression, filename = "Fig2i_p16expression.svg", width = 6, height = 5)


#### Cohen's D calculation
#Table 1

library(lsr)

tp53_master <- read.csv("~/Downloads/TP53_Master.csv")
# Based on p53_Mut
cohen_tp53 <- cohensD((tp53_master%>% filter(p53_Mut == "p53_Mut"))$TP53.overexpression_LFC, 
                      (tp53_master%>% filter(p53_Mut == "p53_Intact"))$TP53.overexpression_LFC)

# Based on p53_LOF
cohen_tp53 <- cohensD((tp53_master%>% filter(p53_LOF == "p53_LOF"))$TP53.overexpression_LFC, 
                      (tp53_master%>% filter(p53_LOF == "p53_Intact"))$TP53.overexpression_LFC)

CDKN2A_master <- read.csv("~/Downloads/p16_CDKN2A_Master.csv")

# Based on CDKN2A_LOF
cohen_CDKN2A <- cohensD((CDKN2A_master%>% filter(CDKN2A_LOF == "p16_LOF"))$p16_LFC, 
                        (CDKN2A_master%>% filter(CDKN2A_LOF == "p16_Intact"))$p16_LFC)

RB_master <- read.csv("~/Downloads/p16_RB_Master.csv")

# Based on RB_LOF
cohen_RB <- cohensD((RB_master%>% filter(RB1_Mut == "RB1_Mut"))$p16_LFC, 
                    (RB_master%>% filter(RB1_Mut == "RB1_Intact"))$p16_LFC)

pten_master <- read.csv("~/Downloads/PTEN_Master_Master.csv")

# Based on pten_mut
cohen_pten <- cohensD((pten_master%>% filter(PTEN_Mut == "PTEN_Mut"))$PTEN_LFC, 
                      (pten_master%>% filter(PTEN_Mut == "PTEN_Intact"))$PTEN_LFC)

# Based on PIK3CA_Hotspot
cohen_pten <- cohensD((pten_master%>% filter(PIK3CA_Hotspot == 1))$PTEN_LFC, 
                      (pten_master%>% filter(PIK3CA_Hotspot == 0))$PTEN_LFC)

ERK2_master <- read.csv("/Users/dkesar/Sellers Lab Dropbox/Devishi Kesar/Devishi Liang Analysis/Finalized figures and codes/Figure_3/ERK2_CohensD.csv")

# Based on p53_LOF
cohen_ERK2 <- cohensD((ERK2_master%>% filter(BRAF_hotspot == 1))$ERK2, 
                      (ERK2_master%>% filter(BRAF_hotspot == 0))$ERK2)


