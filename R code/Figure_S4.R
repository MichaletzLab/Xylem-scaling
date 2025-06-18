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

#Subset stem data 

df_full_stems <- df_full[which(df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_stems$Organ<- factor(df_full_stems$Organ, levels=(c("Twig", "Branch", "Trunk"))) 

#Figure----

data_plot_full_D_d<- ggplot(df_full_stems, aes(x=Diam_org, y=DAVG))+
  geom_pointdensity(alpha=1, adjust =0.5) + 
  scale_color_viridis_c(option = "viridis", name= "Point density", alpha = 0.7)+
  annotate("rect", xmin = 0.001, xmax = Inf, ymin = 1, ymax = Inf,
           fill = "white", alpha = 0.3) +
  stat_ma_line(linewidth=1, method = "SMA", color= "black") +
  scale_y_log10()+
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1), labels= c(0.001, 0.01, 0.1, 1))+
  xlab("Stem diameter (m)")+
  ylab("Xylem diameter (µm)")+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    aspect.ratio = 1
  )

ggsave("Figure_S4.JPG", plot = data_plot_full_D_d, width = 12, height = 10, units = "cm", dpi = 600)

