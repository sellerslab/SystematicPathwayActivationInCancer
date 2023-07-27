###Extended Fig 2

library(readr)
library(tidyverse)
library(svglite)
library(ggrepel)

#TP53
P53GOF_P53Dep <- read_csv("P53GOF_P53Dep.csv")
TP53_Master <- read_csv("TP53_Master.csv")
dataR1_5_1 <- inner_join (P53GOF_P53Dep, TP53_Master, by = c("DepMap ID"))

FigureR1_5_1 <- ggplot(dataR1_5_1, aes(y=dataR1_5_1$`TP53 Gene Effect (Chronos) CRISPR (DepMap Public 22Q4+Score, Chronos)`,x = dataR1_5_1$`p53 LogFC TSGs`)) +
  geom_point(size = 4, alpha = 0.8,aes(color = factor(dataR1_5_1$p53_Mut))) +  
  labs(x = "p53_ORF_lethality", y = "p53_CRISPR_KD_Lethality") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=0), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=0),
        legend.key.size = unit(1, 'cm'))+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black")+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black")


ggsave(FigureR1_5_1, filename = "FigureR1_5_1.pdf", width = 7, height = 7)


#p16
p16GOF_CDK6dep <- read_csv("p16GOF_CDK6dep.csv")
p16_RB_Master <- read_csv("p16_RB_Master.csv")
dataR1_5_2 <- inner_join (p16GOF_CDK6dep, p16_RB_Master, by = c("DepMap ID"))

FigureR1_5_2 <- ggplot(dataR1_5_2, aes(y=dataR1_5_2$`CDK6 Gene Effect (Chronos) CRISPR (DepMap Public 22Q4+Score, Chronos)`,x = dataR1_5_2$`p16 LogFC TSGs`)) +
  geom_point(size = 4, alpha = 0.8,aes(color = factor(dataR1_5_2$RB1_Mut))) +  
  labs(x = "p53_ORF_lethality", y = "p53_CRISPR_KD_Lethality") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=0), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=0),
        legend.key.size = unit(1, 'cm'))+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black")+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black")


ggsave(FigureR1_5_2, filename = "FigureR1_5_2.pdf", width = 7, height = 7)

#PTEN
PTENGOF_PIK3CAdep <- read_csv("PTENGOF_PIK3CAdep.csv")
PTEN_Master_Master <- read_csv("PTEN_Master_Master.csv")
dataR1_5_3 <- inner_join (PTENGOF_PIK3CAdep, PTEN_Master_Master, by = c("DepMap ID"))

FigureR1_5_3 <- ggplot(dataR1_5_3, aes(y=dataR1_5_3$`PIK3CA Gene Effect (Chronos) CRISPR (DepMap Public 22Q4+Score, Chronos)`,x = dataR1_5_3$`PTEN LogFC TSGs`)) +
  geom_point(size = 4, alpha = 0.8,aes(color = factor(dataR1_5_3$PIK3CA))) +  
  labs(x = "p53_ORF_lethality", y = "p53_CRISPR_KD_Lethality") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=0), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=0),
        legend.key.size = unit(1, 'cm'))+
  geom_vline(aes(xintercept = 0), linetype = "dashed", color = "black")+
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "black")


ggsave(FigureR1_5_3, filename = "FigureR1_5_3.pdf", width = 8, height = 7)