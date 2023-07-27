######## Source code Figure3

#setwd to "fig3_extended3_plottingfiles"

#Data Loading

library(readr)
library(tidyverse)
library(ggpubr)
library(ggrepel)
#sample processing
D11_MAPK <- read_csv("D11_MAPK.csv")
MAPK_mutations <- read_csv("MAPK_mutations.csv")
protein_array <- read.csv("Protein_Array.csv")
combined_D11_MAPK <- inner_join(MAPK_mutations,D11_MAPK) 

write.csv(combined_D11_MAPK,"combined_D11_MAPK.csv")
combined_D11_MAPK_anno <- combined_D11_MAPK %>% add_column(RAS_mutant = ((combined_D11_MAPK$HRAS_hotspot + combined_D11_MAPK$KRAS_hotspot + combined_D11_MAPK$NRAS_hotspot)>0.99), RTK_mutant = ((combined_D11_MAPK$EGFR_hotspot + combined_D11_MAPK$ERBB2_hotspot + combined_D11_MAPK$FGFR2_hotspot + combined_D11_MAPK$PDGFRA_hotspot)>0.99)) %>% add_column(MAPK_mutant = ((combined_D11_MAPK$HRAS_hotspot + combined_D11_MAPK$KRAS_hotspot + combined_D11_MAPK$NRAS_hotspot+combined_D11_MAPK$EGFR_hotspot + combined_D11_MAPK$ERBB2_hotspot + combined_D11_MAPK$FGFR2_hotspot + combined_D11_MAPK$PDGFRA_hotspot+combined_D11_MAPK$BRAF_hotspot+combined_D11_MAPK$NF1_hotspot)>0.99)) %>% mutate(BRAFNRAS = case_when(
  BRAF_hotspot == 1  ~ "BRAF_mutation",
  NRAS_hotspot == 1 ~ "NRAS_mutation",
  BRAF_hotspot == 0 & BRAF_hotspot == 0 ~ "Double_WT"))
  

new_MAPK_all <- combined_D11_MAPK_anno%>% mutate(Hotspot_Mutation = case_when(
    RAS_mutant == TRUE ~ "RAS_mutant",
    RTK_mutant == TRUE ~ "RTK_mutant",
    BRAF_hotspot == 1 ~ "RAF_mutant",
    (RAS_mutant + RTK_mutant + BRAF_hotspot) >1 ~"Double_mutation",
    (RAS_mutant + RTK_mutant + BRAF_hotspot) == 0 ~ "No_mutation")) %>% mutate(BRAFmut_plotting = case_when(
      BRAF_hotspot == TRUE ~ "YES",BRAF_hotspot == FALSE ~ "NO"))

MAPK_Protein <- inner_join(new_MAPK_all,protein_array,by = "DepmapID") 

my_comparisons <- list( c("No_mutation","RAF_mutant"), c("No_mutation","RAS_mutant"), c("No_mutation","RTK_mutant"))


#Figure 3a ERK2 Dependency
ERK2_dep_vol <- read_csv("ERK2_dep_vol.csv")

fig3a_ERK2dep <- ggplot(ERK2_dep_vol, aes(x = Correlation, y = `#NAME?`, label = Feature)) +
  geom_point(aes(size = abs(Correlation), color = Correlation, alpha=0.75)) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% c("DUSP4","BRAF","MAPK1","DUSP5","PTPN11","GRB2")),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL, 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',size = 24),
        axis.title= element_text(color='black',size = 24),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=24))
ggsave(fig3a_ERK2dep, filename = "Fig3a_ERK2dep.png", width = 6, height = 5)

# Figure 3b BRAF mutation
Fig3b_ERK2mut <- ggplot(new_MAPK_all, aes(y=ERK2,x = factor(BRAFmut_plotting),fill = factor(BRAFmut_plotting))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("BRAF hotspot mutation")+
  ylab("ERK2 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=28), 
        axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=28),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 8) 



ggsave(Fig3b_ERK2mut, filename = "Fig3b_ERK2mut.pdf", width = 7, height = 7)

#3c scatter

BRAF_Dep <- read_csv("BRAF_Dep.csv")
data3c <- inner_join (BRAF_Dep, new_MAPK_all)

Fig3cNew <- ggplot(data3c, aes(y=data3c$BRAF_CRISPR_KD_CERES,x = data3c$ERK2)) +
  geom_point(size = 4, alpha = 0.8,aes(color = factor(BRAFmut_plotting))) +  
  labs(x = "ERK2_ORF_lethality", y = "BRAF_CRISPR_KD_Lethality") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=0), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=0),
        legend.key.size = unit(1, 'cm'))

ggsave(Fig3cNew, filename = "Fig3cNew.svg", width = 7, height = 7)


ggplot(data3c, aes(y=data3c$BRAF_CRISPR_KD_CERES,x = data3c$ERK2)) +
  geom_point(size = 4, alpha = 0.8,aes(color = factor(BRAFNRAS))) +  
  labs(x = "ERK2_ORF_lethality", y = "BRAF_CRISPR_KD_Lethality") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=0), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=0),
        legend.key.size = unit(1, 'cm'))


#Fig3d different mutations
library(ggpubr)
F3D1 <- ggplot(new_MAPK_all, aes(y=ERK2,x = factor(Hotspot_Mutation),fill = factor(Hotspot_Mutation))) + 
  geom_boxplot()+       
  geom_jitter(position=position_dodge(0.8),alpha = 0.5)+
  xlab("ERK2")+
  ylab("LogFC")+
  labs(fill = "Hotspot Mutation Status")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12),axis.title.x=element_text(size = 24),axis.title.y = element_text(size = 16),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.key.size = unit(1, 'cm'))+
  stat_compare_means(label = "t.test",comparisons = my_comparisons)

F3D2 <- ggplot(new_MAPK_all, aes(y=MEKDD,x = factor(Hotspot_Mutation),fill = factor(Hotspot_Mutation))) + 
  geom_boxplot()+       
  geom_jitter(position=position_dodge(0.8),alpha = 0.5)+
  xlab("MEK1-DD")+
  ylab(" ")+
  labs(fill = "Hotspot Mutation Status")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12),axis.title.x=element_text(size = 24),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  stat_compare_means(label = "t.test",comparisons = my_comparisons)


F3D3 <- ggplot(new_MAPK_all, aes(y=BRAFV600E,x = factor(Hotspot_Mutation),fill = factor(Hotspot_Mutation))) + 
  geom_boxplot()+       
  geom_jitter(position=position_dodge(0.8),alpha = 0.5)+
  xlab("BRAFV600E")+
  ylab(" ")+
  labs(fill = "Hotspot Mutation Status")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12),axis.title.x=element_text(size = 24),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  stat_compare_means(label = "t.test",comparisons = my_comparisons)


F3D4 <- ggplot(new_MAPK_all, aes(y=EGFRdel19,x = factor(Hotspot_Mutation),fill = factor(Hotspot_Mutation))) + 
  geom_boxplot()+       
  geom_jitter(position=position_dodge(0.8),alpha = 0.5)+
  xlab("EGFRexon19del")+
  ylab(" ")+
  labs(fill = "Hotspot Mutation Status")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 12),axis.title.x=element_text(size = 24),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  stat_compare_means(label = "t.test",comparisons = my_comparisons)


fig3d <- ggarrange(F3D1, F3D2, F3D3, F3D4,
                    ncol = 4, nrow = 1,common.legend = TRUE)

ggsave(fig3d, filename = "Fig3d_panel.svg", width = 15, height = 5)

#fig3e protein array
fig3e <- ggplot(MAPK_Protein, aes(y=MAPK_Protein$MEK1_pS217_S221,x = factor(Hotspot_Mutation),fill = factor(Hotspot_Mutation))) + 
  geom_boxplot()+       
  geom_jitter(position=position_dodge(0.8),alpha = 0.5)+
  xlab("MEK1_pS217_S221")+
  ylab("RPPA protein level")+
  labs(fill = "Hotspot Mutation Status")+
  theme_classic()+
  theme(axis.text.y = element_text(size = 15),axis.title.x=element_text(size = 20),axis.title.y=element_text(size = 15),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15),
        legend.key.size = unit(0.7, 'cm'))+
  stat_compare_means(label = "t.test",comparisons = my_comparisons)

ggsave(fig3e, filename = "Fig3e_panel.svg", width = 8, height = 5)


#fig 3f ERK2 RPPA
ERK2RPPA <- read_csv("MEK1RPPA_ERK2.csv")

Fig3f <- ggplot(ERK2RPPA,aes(x= ERK2RPPA$`ERK2 LFC LC`,y= ERK2RPPA$`MEK1_pS217_S221 RPPA signal (log2) Protein Array`)) +
  geom_point(size = 1, alpha = 0.6) +  
  labs(x = "ERK2_overexpression_LogFC", y = "MEK1_pS217_S221 RPPA Signal") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16))+
  geom_smooth(method=lm,se = FALSE)

cor(ERK2RPPA$`ERK2 LFC LC`,ERK2RPPA$`MEK1_pS217_S221 RPPA signal (log2) Protein Array`)
ggsave(Fig3f, filename = "Fig3f.svg", width = 8, height = 5)






# figure 3 individual p values calculated
#p value not significant

Fig3b_ERK2mut <- ggplot(new_MAPK_all%>% filter(), aes(y=ERK2,x = factor(BRAFmut_plotting),fill = factor(BRAFmut_plotting))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("BRAF hotspot mutation")+
  ylab("ERK2 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16), 
        axis.text.y = element_text(size=16),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 6) 

ggsave(Fig3b_ERK2mut, filename = "Fig3b_ERK2mut.svg", width = 4, height = 5)


BRAF_1_lineage <- new_MAPK_all %>% filter(Lineage %in% c("Colorectal","Lung","Skin","Uterus","Ovary","Thyroid"))
Fig3c_ERK2mutlineage <- ggplot(BRAF_1_lineage , aes(y=ERK2,x = Lineage,fill = factor(BRAFmut_plotting))) + 
  geom_boxplot()+
  geom_jitter(position=position_dodge(0.8),alpha = 0.5)+ 
  xlab("BRAF hotspot mutation")+
  ylab("ERK2 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.x = element_text(size=16,angle=45,hjust = 1), 
        axis.text.y = element_text(size=16),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 6) 

ggsave(Fig3c_ERK2mutlineage , filename = "Fig3c_ERK2mut.svg", width = 6.5, height = 5)

#Extended 3
# NRAS plotting

D11_MAPK <- read_csv("D11_MAPK.csv")
ggplot(new_MAPK_all, aes(y=ERK2,x = factor(NRAS_hotspot),fill = factor(NRAS_hotspot))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("NRAS hotspot mutation")+
  ylab("ERK2 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=28), 
        axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=28),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 8)

ggplot(new_MAPK_all, aes(y=ERK2,x = factor(BRAFNRAS),fill = factor(BRAFNRAS))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("BRAF hotspot mutation")+
  ylab("ERK2 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=28), 
        axis.text.x = element_text(size=28), 
        axis.text.y = element_text(size=28),legend.position='none')+
  stat_compare_means(method = "t.test",label.y = 8,size = 8)


BRAFNRAS_comparisons <- list( c("Double_WT","BRAF_mutation"), c("Double_WT","NRAS_mutation"))
ggplot(new_MAPK_all, aes(y=ERK2,x = factor(BRAFNRAS),fill = factor(BRAFNRAS))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+
  ylab("ERK2 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=24), 
        axis.text.y = element_text(size=24),legend.text = element_text(size=24),        
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=15),
        legend.key.size = unit(0.7, 'cm'))+
  stat_compare_means(method = "t.test",size = 5,comparisons = BRAFNRAS_comparisons) +
  scale_x_discrete(limits = c("Double_WT","BRAF_mutation","NRAS_mutation"))

new_MAPK_all %>% filter(Lineage == "Skin")
ggplot(new_MAPK_all %>% filter(Lineage == "Skin"), aes(y=ERK2,x = factor(BRAFNRAS),fill = factor(BRAFNRAS))) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+
  ylab("ERK2 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=24), 
        axis.text.y = element_text(size=24),legend.text = element_text(size=24),        
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=15),
        legend.key.size = unit(0.7, 'cm'))+
  stat_compare_means(method = "t.test",size = 5,comparisons = BRAFNRAS_comparisons) +
  scale_x_discrete(limits = c("Double_WT","BRAF_mutation","NRAS_mutation"))