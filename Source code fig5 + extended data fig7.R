
#setwd "fig5_extendedfig7_plottingfiles"
library(readr)
library(tidyverse)
library(ggpubr)
library(ggrepel)

#Data processing
CTNNB_mut_clean <- read_csv("CTNNB_mut_clean.csv") %>% mutate(Hotspot_Mutation = case_when(
  APC_hotspotMut == 1 & CTNNB1_hotspotMut == 0 ~ "APC_mutation",
  CTNNB1_hotspotMut == 1 & APC_hotspotMut == 0 ~ "CTNNB1_mutation",
  APC_hotspotMut == 0 & CTNNB1_hotspotMut == 0 ~ "Double_WT")) %>% filter(`DepMap ID` != "ACH-000991") 

Damaging_Mutations_APC <- read_csv("Damaging_Mutations_APC.csv")

CTNNB_mut_clean1 <- inner_join(read_csv("CTNNB_mut_clean.csv"),read_csv("Damaging_Mutations_APC.csv")) %>% mutate(Hotspot_Mutation = case_when(
  APC_hotspotMut == 1 & CTNNB1_hotspotMut == 0 ~ "APC_mutation",
  CTNNB1_hotspotMut == 1 & APC_hotspotMut == 0 ~ "CTNNB1_mutation",
  APC_hotspotMut == 0 & CTNNB1_hotspotMut == 0 ~ "Double_WT")) %>% filter(`DepMap ID` != "ACH-000991") %>%
  mutate(Key_WNT_Mutation = case_when(
    APC_hotspotMut == 1 | APC_Damaging_mut == 1 ~ "APC_mutation",
    CTNNB1_hotspotMut == 1 & APC_hotspotMut == 0 ~ "CTNNB1_mutation",
    APC_hotspotMut == 0 & CTNNB1_hotspotMut == 0 ~ "Double_WT")) %>%
  mutate(APC_mutations_CRC = case_when(
    (APC_hotspotMut == 1 | APC_Damaging_mut == 1) & Lineage =="Colorectal" ~ "YES"))%>%
  mutate(APC_mutation = case_when(
    APC_hotspotMut == 1 | APC_Damaging_mut == 1 ~ "APC_mutation",
    APC_hotspotMut == 0 & APC_Damaging_mut == 0 ~ "APC_WT")) %>%
  mutate(APC_hotspot_CRC = case_when(
    (APC_hotspotMut == 1 ) & Lineage =="Colorectal" ~ "YES"))

# Fig 5a APC CTNNB1 mut
Fig5a_mut <- ggplot(CTNNB_mut_clean, aes(y=CTNNB_mut_clean$CTNNB_LFC,x = Hotspot_Mutation,fill = Hotspot_Mutation)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=2,binwidth=0.05)+ 
  xlab("WNT pathway hotspot mutation")+
  ylab("CTNNB1 overexpression_LogFC")+
  theme_classic()+
  theme(axis.title=element_text(size=16), 
        axis.text.y = element_text(size=16),legend.text = element_text(size=15),        
        axis.text.x=element_text(size=16),
        axis.ticks.x=element_blank(),
        legend.title = element_text(size=15))+
  stat_compare_means(method = "t.test",size = 5,comparisons = WNT_comparisons) +
  scale_x_discrete(limits = c("Double_WT","CTNNB1_mutation","APC_mutation"))
ggsave(Fig5a_mut, filename = "Fig5a_APCmut.pdf", width = 8, height = 5)





#Fig 5b 
CTNNB_vol_dep <- read_csv("~/Sellers Lab Dropbox/Liang Chang/Bioinfo analysis/4Q2020_PRISM016_TheScreen/DAC_Heatmaps/CTNNB_vol_dep.csv")

fig5b_CTNNB1dep <- ggplot(CTNNB_vol_dep, aes(x = Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(size = abs(Correlation), color = Correlation, alpha=0.75)) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% c("CSNK1A1","APC")),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = NULL, 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(size = 20),
        axis.title= element_text(size = 20),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=20))
ggsave(fig5b_CTNNB1dep, filename = "fig5b_CTNNB1dep.svg", width = 6, height = 5)

#Fig 5c APC depency assosiation
Depmap_APC_Co_dependency <- read_csv("Depmap APC Co-dependency.csv", 
                                     col_types = cols(`-log10(PValue)` = col_double()))

Fig5c <- ggplot(Depmap_APC_Co_dependency[-1,], aes(x = Correlation, y = `-log10(PValue)`, label = Feature)) +
  geom_point(aes(size = abs(Correlation), color = Correlation, alpha=0.75)) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    data = .%>%filter(Feature %in% c("AXIN1","CSNK1A1","CTNNBIP1","AXIN2")),# change gene name 
    aes(label = Feature),
    size = 6,
    box.padding = unit(1, "lines"),
    point.padding = unit(1, "lines"))+
  labs(title = NULL,
       x = NULL, 
       y = "-log10(p_value)")+
  theme(aspect.ratio = 1, 
        plot.title = element_text(color='black', hjust = 0.5),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.grid = element_blank(),
        axis.text = element_text(color='black',size = 20),
        axis.title= element_text(color='black',size = 20),
        legend.key = element_rect(fill = NA),
        legend.position ="none",
        text = element_text(size=20))
ggsave(Fig5c, filename = "fig5c_APCdep.svg", width = 6, height = 5)



#fig 5d
CK1_APC_CTNNB <- inner_join(APC_CK1_CRISPR,CTNNB_mut_clean,by = "DepMap ID")
Fig5d <- ggplot(CK1_APC_CTNNB, aes(y=CK1_APC_CTNNB$APC_CRISPR_KD_CERES,x = CK1_APC_CTNNB$CSNK1A1_CRISPR_KD_CERES)) +
  geom_point(alpha = 0.8,aes(color = factor(`APC any Mutation`),size = -CK1_APC_CTNNB$CTNNB_LFC)) +  
  labs(x = "CSNK1A1 CRISPR KD Effect (CERES)", y = "APC CRISPR KD Effect (CERES)") +
  theme_classic()+
  theme(aspect.ratio=1,axis.title=element_text(size=0), 
        axis.text.x = element_text(size=20), 
        axis.text.y = element_text(size=20),legend.text = element_text(size=20),
        legend.title = element_text(size=0),
        legend.key.size = unit(1, 'cm'))
ggsave(Fig5d, filename = "Fig5d.svg", width = 6, height = 6)


###### RNAseq fig 5h, Extended Fig. 7

rm(list=ls())
##################################
# Loading libraries and functions
##################################
library(msigdbr)
library(dplyr)
library(openxlsx)

volcano_plot <- function(file, name = "output", highlights = "", highlight2 = "",
                         dotSize = 1, title = "", textSize = 6, labelSize = 7){
  pdf(file= paste0(name,".pdf"),width = 9, height = 7)
  plot_file <- read.csv(file)
  plot_file$color <- ifelse(plot_file$logFC <=0 & plot_file$P.Value < 0.01,"gray",
                            ifelse(plot_file$logFC >=0 & plot_file$P.Value < 0.01,"red","black"))
  plot_file$color <- ifelse(plot_file$id %in% highlights,"mediumblue",
                            ifelse(plot_file$id %in% highlight2,"darkorange3",
                                   plot_file$color))
  highlights_wnt <- (plot_file %>% filter(P.Value<=0.01)%>% filter(logFC >= 0| id == "APC") %>% filter(id %in% highlights))$id
  highlights_e2f <- (plot_file %>% filter(P.Value<=0.01)%>% filter(abs(logFC) >= 0.3)%>% filter(id %in% highlight2))$id
  gg <- ggplot(data = plot_file, aes(x = logFC, y = -log10(P.Value), label = id)) +
    geom_point(color = plot_file$color, size = dotSize) + theme_classic() + 
    ggrepel::geom_label_repel(data = subset(plot_file, id %in% highlights_wnt), 
                              mapping = aes(x = logFC, y = -log10(P.Value), 
                                            label = id),fill = "lightblue1",
                              xlim = c(1,NA),  # <--- here
                              size = labelSize, segment.color = 'black',
                              box.padding = 0.3) + ggtitle(title)+
    ggrepel::geom_label_repel(data = subset(plot_file, id %in% highlights_e2f), 
                              mapping = aes(x = logFC, y = -log10(P.Value), 
                                            label = id), fill = "lightgoldenrodyellow",
                              xlim = c(NA,-1),  # <--- here
                              size = labelSize, segment.color = 'black',
                              box.padding = 0.3) + ggtitle(title)+
    scale_color_manual(values=plot_file$color)+ xlab("log2(FC)") + ylab("-log10(P value)") + 
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      legend.position = "none") +
    geom_hline(yintercept=2, linetype='dotted')
  print(gg)
  dev.off()
}

# APC shRNA vs WT
HALLMARK_E2F_Targets <- (msigdbr(species = "human", category = "H") %>% 
                           filter(gs_name == "HALLMARK_E2F_TARGETS"))$gene_symbol
HALLMARK_G2M_CHECKPOINT <- (msigdbr(species = "human", category = "H") %>% 
                              filter(gs_name == "HALLMARK_G2M_CHECKPOINT"))$gene_symbol
REACTOME_CELL_CYCLE_CHECKPOINTS <- (msigdbr(species = "human", category = "C2") %>% 
                                      filter(gs_name == "REACTOME_CELL_CYCLE_CHECKPOINTS"))$gene_symbol
REACTOME_MITOTIC_METAPHASE_AND_ANAPHASE <- (msigdbr(species = "human", category = "C2") %>% 
                                              filter(gs_name == "REACTOME_MITOTIC_METAPHASE_AND_ANAPHASE"))$gene_symbol
REACTOME_DNA_REPLICATION_PRE_INITIATION <- (msigdbr(species = "human", category = "C2") %>% 
                                              filter(gs_name == "REACTOME_DNA_REPLICATION_PRE_INITIATION"))$gene_symbol
set.seed(40)


imp_genes_APC <- read.xlsx("curated gene list_RNAseq heatmap.xlsx", sheet  = 1,colNames = F)
imp_genes_BCatenin <- read.xlsx("/curated gene list_RNAseq heatmap.xlsx", sheet  = 2,colNames = F)

file <- "shRNA-NT1vshRNA-APC.limma_star.csv"
set.seed(70)
volcano_plot(file, highlights = wnt_gene_list,highlight2 = imp_genes_APC$X1,name = "APC_shRNA_vs_WT_volcano")

file <- "CTNNB1vGFP.limma_star.csv"
set.seed(70)
volcano_plot(file, highlights = wnt_gene_list,highlight2 = imp_genes_BCatenin$X1,name = "BetaCatenin_shRNA_vs_GFP_volcano")


library(GSA)

create_own_hmark_gset <- GSA.read.gmt("h.all.v7.5.1.symbols.gmt")
names(create_own_hmark_gset$genesets) <- create_own_hmark_gset$geneset.names

library(clusterProfiler)
library(enrichplot)
library(readxl)
library(ggplot2)
library(stringr)
library(dplyr)
library(scales)
library(readr)


plot_enrichment <- function(file, name = "output", thres = 0.005, x_axis_label=""){
  dataset <- read.csv(file) %>% dplyr::select(id,logFC) 
  gene_list <- dataset$logFC
  names(gene_list) <- dataset$id
  gene_list <- sort(gene_list, decreasing = TRUE)
  ids <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
  dedup_ids <- ids[!duplicated(ids[c("SYMBOL")]),]
  gene_list <- gene_list[names(gene_list) %in% dedup_ids$SYMBOL]
  
  library(msigdbr)
  library(fgsea)
  genesets_list <- create_own_hmark_gset$genesets
  fgseaRes <- fgsea(genesets_list, gene_list)
  fgseaRes_save <- fgseaRes %>% arrange(padj) %>% apply(2,as.character)
  write.csv(fgseaRes_save,paste0(name,".csv"))
  fgseaRes$logpval <- ifelse(fgseaRes$NES>0, -log10(fgseaRes$pval),log10(fgseaRes$pval))
  a <- fgseaRes %>% filter(pval < thres) %>%
    arrange(logpval) %>%
    ggplot(aes(x=reorder(pathway, -logpval), 
               y=logpval, fill = ifelse(NES>0, "Blue","red"))) + geom_bar(stat="identity")+ coord_flip()+
    theme_classic() +
    ylab(label= x_axis_label)+ xlab(label ="Pathway") + 
    theme(legend.position="none",
          axis.title.x = element_text(size = 20,face = "bold"),
          axis.text.x = element_text(size = 20,face = "bold"),
          axis.text.y = element_text(size = 20,face = "bold"),
          axis.title.y = element_text(size = 20,face = "bold")) +
    scale_y_continuous(breaks = c(-10, -5, 0,5,10), limits = c(-11, 10)) 
  
  ggsave(a, filename = paste0(name, '.pdf'), device = 'pdf',width = 14, height = 4)
  
}

# APC shRNA vs WT
file <- "shRNA-NT1vshRNA-APC.limma_star.csv"
plot_enrichment(file, name = "APC", x_axis_label = "<-   Depleted     -    log10(pvalue)  -   Enriched ->")

# beta catenin vs GFP
file <- "CTNNB1vGFP.limma_star.csv"
plot_enrichment(file, name = "BetaCatenin", x_axis_label = "<-Depleted       -        log10(pvalue)      -      Enriched ->")

heatmap_fun <- function(name = "output", fontSize = 11, Treatment_l = c("shRNA-APC","shRNA-NT1","GFP","CTNNB1"),
                        use_z_score = TRUE,wnt_only = TRUE, height = 7,highlight_genes = NULL,
                        gs_name = "WNT genes"){
  if(length(highlight_genes) != 0){
    highlight_genes <- highlight_genes
  }
  else{
    if(wnt_only){
      highlight_genes <- wnt_gene_list
    }
    else{
      highlight_genes <- c(create_own_hmark_gset$genesets$KEGG_CELL_CYCLE)
    }
  }
  library(ComplexHeatmap)
  library(tidyverse)
  library(tidyr)
  library(matrixStats)
  pdf(file= paste0(name,".pdf"),height = height)
  star_gene <- read.csv("STAR_Gene_TPM.csv") %>%
    gather("SampleName","TPM",-X) %>% inner_join(read.csv("summary_reports.metasheet.csv")[,1:2]) %>%
    group_by(X,Treatment) %>%
    summarise(mean_TPM = log2(mean(TPM+1)))
  heatmap_plot <- star_gene %>% filter(X %in% highlight_genes) %>% filter(Treatment %in% Treatment_l) %>%
    mutate(type = gs_name) 
  heatmap_matrix <- heatmap_plot %>% dplyr::select(-type) %>% spread(key = Treatment, value = mean_TPM) %>%
    column_to_rownames(var = "X") %>% as.matrix() 
  if(use_z_score){
    heatmap_matrix  <- (heatmap_matrix  - rowMeans(heatmap_matrix ))/rowSds(heatmap_matrix)
    heatmap_matrix <- heatmap_matrix[!rowSums(is.na(heatmap_matrix)),]
    bs <- "Z-score"
  }else{bs <- "log2(TPM+1)"}
  
  print(Heatmap(heatmap_matrix,cluster_columns = FALSE,name = bs))
  dev.off()
}

heatmap_fun(use_z_score = TRUE,wnt_only = TRUE,name = "Z-score_WNTonly")
heatmap_fun(use_z_score = FALSE,wnt_only = TRUE,name = "TPM_WNTonly")
heatmap_fun(use_z_score = TRUE,wnt_only = FALSE,name = "Z-score_KEGGonly", height = 20)
heatmap_fun(use_z_score = FALSE,wnt_only = FALSE,name = "TPM_KEGGonly", height = 20)




