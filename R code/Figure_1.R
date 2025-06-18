#Preparation----
##Install required packages----

#List of packages
pkg_list <- c("ggplot2", "dplyr", "readxl", "smatr", "cowplot", 
              "ggpubr", "ggpointdensity", "ggpmisc", "emmeans", 
              "multcomp", "ggridges", "patchwork", "segmented", 
              "viridis", "Hmisc")

#Function to check and install missing packages
check_and_install <- function(pkg){
  if (!require(pkg, character.only = TRUE)){
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

#Apply the function to the list of packages
sapply(pkg_list, check_and_install)

##Load the required packages----

library(ggplot2)
library(dplyr)
library(readxl)
library(smatr)
library(cowplot)
library(ggpubr)
library(ggpointdensity)
library(ggpmisc)
library(emmeans)
library(multcomp)
library(ggridges)
library(patchwork)
library(segmented)
library(viridis)
library(Hmisc)

##Set options----

options(scipen = 999) #prevent scientific notation on plots

##Import data----

df_full <- read_excel(list.files(path = "/", pattern = "xylem_scaling_data\\.xlsx$", recursive = TRUE, full.names = TRUE)[1])

#Analyses----

#Reorganize organs in correct order

Organ_order<- rev(c("Leaf", "Twig", "Branch", "Trunk", "Coarse root", "Fine root", "Very fine root"))

df_full$Organ<- factor(df_full$Organ, levels = Organ_order)

df_full$Organ_order_num<- as.numeric(factor(df_full$Organ, levels=c("Leaf", "Twig", "Branch", "Trunk", "Coarse root", "Fine root", "Very fine root")))

#Fit ANOVA model for d~Organ

ANOVA_FIG1_A<- lm(log10(DAVG)~Organ, data = df_full)

summary.aov(ANOVA_FIG1_A)

#Post-hoc test

ANOVA_FIG1_A_emmeans<- emmeans(ANOVA_FIG1_A, specs = "Organ")
ANOVA_FIG1_A_pairs<- pairs(ANOVA_FIG1_A_emmeans, adjust="sidak")
ANOVA_FIG1_A_CLD<- cld(object = ANOVA_FIG1_A_emmeans,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

ANOVA_FIG1_A_pairs

ANOVA_FIG1_A_CLD

10^(ANOVA_FIG1_A_CLD$emmean)

#Segmented regression d~Organ

ANOVA_FIG1_A<- lm(log10(DAVG)~Organ_order_num, data = df_full)

ANOVA_seg_FIG1_A <- segmented(ANOVA_FIG1_A, seg.Z = ~Organ_order_num)

summary(ANOVA_seg_FIG1_A)

#Fit ANOVA model for t/b^2~organ

ANOVA_FIG1_B<- lm(log10(t_b_2)~Organ, data = df_full)

summary.aov(ANOVA_FIG1_B)

ANOVA_FIG1_B_emmeans<- emmeans(ANOVA_FIG1_B, specs = "Organ")
ANOVA_FIG1_B_pairs<- pairs(ANOVA_FIG1_B_emmeans)
ANOVA_FIG1_B_CLD<- cld(object = ANOVA_FIG1_B_emmeans,
                       adjust = "sidak",
                       Letters = letters,
                       alpha = 0.05)

ANOVA_FIG1_B_pairs

ANOVA_FIG1_B_CLD

10^(ANOVA_FIG1_B_CLD$emmean)

#Segmented regression (t/b)^2~Organ

ANOVA_FIG1_B<- aov(log10(t_b_2)~Organ_order_num, data = df_full)

ANOVA_seg_FIG1_B <- segmented(ANOVA_FIG1_B, seg.Z = ~Organ_order_num)

summary(ANOVA_seg_FIG1_B)

#Figure----

CLD_df_d <- data.frame(Organ = ANOVA_FIG1_A_CLD$Organ, 
                       Letters = trimws(ANOVA_FIG1_A_CLD$.group),
                       emmean = 10^(ANOVA_FIG1_A_CLD$emmean))

Fig1_d <- ggplot(df_full, aes(x = DAVG, y = Organ)) + 
  geom_density_ridges2(scale=2, alpha=0.2, rel_min_height = 0.00002, quantile_lines = T, quantile_fun = mean, color="#0d0887", fill="#0d0887")+
  geom_text(data = CLD_df_d, aes(y = Organ, x = emmean, label = Letters), hjust = 1.5, vjust = -1, size = 5, color = "#cc4778") +
  theme_classic2()+
  xlab("Xylem conduit diameter (µm)")+
  theme(axis.title.y = element_blank(),
        text = element_text(size = 15),
        aspect.ratio = 1)+
  scale_x_log10(limits= c(1, 80), expand = c(0, 0.02))

Fig1_d

CLD_df_t_b_2 <- data.frame(Organ = ANOVA_FIG1_B_CLD$Organ, 
                           Letters = trimws(ANOVA_FIG1_B_CLD$.group),
                           emmean = 10^(ANOVA_FIG1_B_CLD$emmean))

Fig1_tb2<- ggplot(df_full, aes(x = t_b_2, y = Organ)) + 
  geom_density_ridges2(scale=2, alpha=0.2, rel_min_height = 0.00002, quantile_lines = T, quantile_fun = mean, color="#0d0887", fill="#0d0887")+
  geom_text(data = CLD_df_t_b_2, aes(y = Organ, x = emmean, label = Letters), hjust = 2.5, vjust = -1.5, size = 5, color = "#cc4778") +
  theme_classic2()+
  xlab("Thickness-to-span ratio (dimensionless)")+
  theme(axis.title.y = element_blank(),
        text = element_text(size = 15),
        aspect.ratio = 1)+
  scale_x_log10(limits= c(0.001, 100), expand = c(0, 0.02))

Fig1_tb2

Fig1_full <- Fig1_d + 
  labs(tag = 'A') +
  Fig1_tb2 + 
  labs(tag = 'B') +
  plot_layout(guides = 'collect') & 
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.1, 0.98))

ggsave("Figure_1.pdf", plot=Fig1_full, width = 25, height = 10, units = "cm", device = cairo_pdf)


