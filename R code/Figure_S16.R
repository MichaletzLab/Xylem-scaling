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

#Subset data

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'| df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk"))) 

#Figure----

FigS16<- ggplot(aes(x=DAVG, y=CWTALL), data=df_full_aboveground)+
  geom_pointdensity(alpha=0.05, adjust =0.5) + 
  scale_color_viridis_c(option = "viridis", name= "Point density")+
  scale_y_continuous(breaks=c(1,3,5,7,9,11,13), limits=c(1,14), expand = c(0,0))+
  scale_x_continuous(breaks=c(0,10,20,30,40,50,60), limits=c(1,60), expand = c(0,0))+
  geom_line(aes(x=DAVG, y=(((sqrt(0.009))*DAVG)/2)), linewidth=1, linetype="dotdash")+
  xlab("Xylem conduit diameter (µm)")+
  ylab("Cell wall thickness (µm)")+
  theme_classic()+
  theme(text = element_text(size = 12),
        title = element_text(size = 14),
        aspect.ratio = 1)

ggsave("Figure_S16.JPG", plot=FigS16, width = 12, height = 10, dpi = 600, units = "cm")