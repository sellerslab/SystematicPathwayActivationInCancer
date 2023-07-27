#Data Loading


library(readr)
library(tidyverse)
library(ggpubr)
library(ggrepel)
#setwd to fig4_plottingfiles

#sample processing
D11_PI3K <- read_csv("D11_PI3K.csv")
PI3K_mut <- read_csv("PI3K_mut.csv")
PTEN_mut <- read_csv("PTEN_mut.csv")

New_PI3K_mut <- inner_join(PI3K_mut,PTEN_mut,by = "DepmapID") 
combined_D11_PI3K <- inner_join(New_PI3K_mut,D11_PI3K,by = "DepmapID") %>% mutate(PTEN = (PTEN_Damaging + PTEN_hotspot) )

#new
new_PI3K_all <- combined_D11_PI3K %>% mutate(Hotspot_Mutation = case_when(
  PTEN > 0  & PIK3CA == 0 ~ "PTEN_mutation",
  PIK3CA == 1 & PTEN == 0 ~ "PIK3CA_Hotspot_mutation",
  PIK3CA == 1 & PTEN > 0 ~ "PTEN+PIK3CA_mutation",
  (PTEN +PIK3CA) == 0 ~ "Double_WT")) %>% mutate (PI3K_mutation = case_when(
    (PTEN + PIK3CA + AKT1) >0  ~ "YES",
    (PTEN + PIK3CA + AKT1) == 0 ~ "NO"))


new_PI3K_all <- combined_D11_PI3K %>% mutate(Hotspot_Mutation = case_when(
  PTEN_Damaging == 1 & PIK3CA == 0 ~ "PTEN_Damaging_mutation",
  PIK3CA == 1 & PTEN_Damaging == 0~ "PIK3CA_Hotspot_mutation",
  (PTEN_Damaging +PIK3CA) >1 ~ "PTEN+PIK3CA_mutation",
  (PTEN_Damaging +PIK3CA) == 0 ~ "Double_WT")) %>% mutate (PI3K_mutation = case_when(
    (PTEN + PIK3CA + PIK3CA + PIK3CD + PIK3CG + PIK3R1 + PIK3R3 + PIK3R4 + PIK3R6 + AKT1 + AKT3) >0  ~ "YES",
    (PTEN + PIK3CA + PIK3CA + PIK3CD + PIK3CG + PIK3R1 + PIK3R3 + PIK3R4 + PIK3R6 + AKT1 + AKT3) == 0 ~ "NO")) 
  


#Fig 4a PIK3CA pathway mutations

Fig4a_PI3Kmut <- ggplot(new_PI3K_all, aes(y=new_PI3K_all$AKTe17k,x = factor(PI3K_mutation),fill = factor(PI3K_mutation))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("PI3K pathway hotspot mutation")+
  ylab("AKTE17K overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 5) 
ggsave(Fig4a_PI3Kmut, filename = "Fig4a_PI3Kmut.svg", width = 4, height = 5)

#Fig 4b double pathway mutations
PI3K_comparisons <- list(c("PTEN_mutation","PIK3CA_mutation"))



PI3K_comparisons <- list( c("Double_WT","PIK3CA_Hotspot_mutation"), 
                          c("Double_WT","PTEN_mutation"), 
                          c("Double_WT","PTEN+PIK3CA_mutation"))
new_PI3K_all

fig4b_double_mutation <- ggplot(new_PI3K_all, aes(y=new_PI3K_all$AKTe17k,x = factor(Hotspot_Mutation),fill = factor(Hotspot_Mutation))) + 
  geom_boxplot()+       
  geom_jitter(position=position_dodge(0.8),alpha = 0.5)+
  xlab("PI3K pathway hotspot mutation")+
  ylab("AKTE17K overexpression_LogFC")+
  labs(fill = "Hotspot Mutation Status")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 15),axis.title.x=element_text(size = 20),axis.title.y=element_text(size = 15),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.key.size = unit(0.7, 'cm'))+
  stat_compare_means(label = "t.test",comparisons = PI3K_comparisons)

ggsave(fig4b_double_mutation, filename = "Fig4b_double.svg", width = 8, height = 5)


# Fig 4d
protein_array <- read.csv("Protein_Array.csv")
PI3K_protein <- inner_join(new_PI3K_all,protein_array,by = "DepmapID") %>% mutate (PTEN_PIK3CA_DualMut = case_when(Hotspot_Mutation == "PTEN+PIK3CA_mutation"  ~ "YES")) %>% mutate (DMUTR = case_when(
  (lineage_1 == "Uterus" & PTEN_PIK3CA_DualMut == "YES")  ~ "YES",TRUE ~"NO"))

DM_UTR_cells <- PI3K_protein %>% filter(lineage_1 == "Uterus" & PTEN_PIK3CA_DualMut == "YES") %>% select(cell_line_display_name)

fig4d <- ggplot(PI3K_protein, aes(y=Akt_pS473,x = Akt_pT308,label = cell_line_display_name)) +
  geom_point(size = 4, alpha = 0.8,aes(color = DMUTR)) +  
  labs(x = "AKT_pS473 RPPA", y = "AKT_pT308 RPPA") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=20), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.key.size = unit(1, 'cm'))

ggsave(fig4d, filename = "Fig4d.svg", width = 9, height = 6)



