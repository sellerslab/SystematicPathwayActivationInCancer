###Extended Fig 4
library(readr)
library(tidyverse)
library(svglite)
library(ggrepel)

#pGDSC1_ERK2
GDSC1_ERK2 <- read_csv("GDSC1_ERK2.csv")

GDSC1_ERK2_fig <- ggplot(GDSC1_ERK2, aes(x = Cor, y = -log10(PValue), label = GDSC1_ERK2$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
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

ggsave(GDSC1_ERK2_fig, filename = "GDSC1_ERK2_fig.pdf", width = 6, height = 5)

#pGDSC1_ERK2_label
GDSC1_ERK2_figlabel <- ggplot(GDSC1_ERK2, aes(x = Cor, y = -log10(PValue), label = GDSC1_ERK2$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  labs(title = NULL,
       x = "Anti-Dependency <--- Correlation ---> Co-Dependency", 
       y = "-log10(p_value)")+
  geom_label_repel(
    aes(label = special_label),
    size = 4,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.5, "lines"),force =3,nudge_x = 1,
    nudge_y = 1,
    force = 1,
    max.iter = 1000,min.segment.length = 0.5)+
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

ggsave(GDSC1_ERK2_figlabel, filename = "GDSC1_ERK2_figlabel.pdf", width = 8, height = 6)



# GDSC1_PTEN
GDSC1_PTEN <- read_csv("GDSC1_PTEN.csv")

GDSC1_PTEN_fig <- ggplot(GDSC1_PTEN, aes(x = Cor, y = -log10(PValue), label = GDSC1_PTEN$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    aes(label = special_label),
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

ggsave(GDSC1_PTEN_fig, filename = "GDSC1_PTEN_fig.pdf", width = 6, height = 5)

# GDSC1_p53
GDSC1_p53 <- read_csv("GDSC1_p53.csv")

GDSC1_p53_fig <- ggplot(GDSC1_p53 , aes(x = Cor, y = -log10(PValue), label = GDSC1_p53$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    aes(label = special_label),
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

ggsave(GDSC1_p53_fig, filename = "GDSC1_p53_fig.pdf", width = 6, height = 5)


# GDSC1_p16
GDSC1_p16 <- read_csv("GDSC1_p16.csv")

GDSC1_p16_fig <- ggplot(GDSC1_p16 , aes(x = Cor, y = -log10(PValue), label = GDSC1_p16$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    aes(label = special_label),
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

ggsave(GDSC1_p16_fig, filename = "GDSC1_p16_fig.pdf", width = 6, height = 5)


# Score_PTEN
SangerCeres_PTEN <- read_csv("SangerCeres_PTEN.csv")

SangerCeres_PTEN_fig <- ggplot(SangerCeres_PTEN , aes(x = Cor, y = -log10(PValue), label = SangerCeres_PTEN$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    aes(label = special_label),
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

ggsave(SangerCeres_PTEN_fig, filename = "SangerCeres_PTEN_fig.pdf", width = 6, height = 5)

# Score_p16
SangerCeres_p16 <- read_csv("SangerCeres_p16.csv")

SangerCeres_p16_fig <- ggplot(SangerCeres_p16, aes(x = Cor, y = -log10(PValue), label = SangerCeres_p16$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    aes(label = special_label),
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

ggsave(SangerCeres_p16_fig, filename = "SangerCeres_p16_fig.pdf", width = 6, height = 5)

# Score_p53
SangerCeres_p53 <- read_csv("SangerCeres_p53.csv")

SangerCeres_p53_fig <- ggplot(SangerCeres_p53, aes(x = Cor, y = -log10(PValue), label = SangerCeres_p53$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    aes(label = special_label),
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

ggsave(SangerCeres_p53_fig, filename = "SangerCeres_p53_fig.pdf", width = 6, height = 5)

# Score_ERK2
SangerCeres_ERK2 <- read_csv("SangerCeres_ERK2.csv")

SangerCeres_ERK2_fig <- ggplot(SangerCeres_ERK2, aes(x = Cor, y = -log10(PValue), label = SangerCeres_ERK2$label)) +
  geom_point(aes(color = Cor, alpha=0.75),size = 3) +
  scale_radius(limits = c(0,1), 
               range = c(2.5,12))+
  scale_y_sqrt() +
  scale_color_gradient2(low="blue", mid = "gray50", high="red", 
                        limits = c(-1,1))+
  geom_label_repel(
    aes(label = special_label),
    size = 4,
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

ggsave(SangerCeres_ERK2_fig, filename = "SangerCeres_ERK2_fig.pdf", width = 6, height = 5)
