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

#Calculate means

df_cats<- df_full[1:4]

df_cats_summarized <- df_cats %>%
  group_by(ID) %>%
  slice(1)

df_full_vars_plusID<- df_full[c(1,5:19)]

df_DAVG_ID<- df_full[c(1,18)]

df_DhAVG<- df_full %>%
  group_by(ID) %>%
  summarise(DhAVG = sum(DAVG^5)/sum(DAVG^4)) #calculate hydraulic diameter

df_means<-  df_full_vars_plusID %>%
  group_by(ID) %>%
  summarise_all(.funs = funs(mean(., na.rm = TRUE)))

df_means<- merge(df_means, df_DhAVG, by= "ID")
df_means<- merge(df_cats_summarized, df_means, by= "ID", all=FALSE)

#Subset data

df_means_stems <- df_means[which(df_means$Organ=='Twig'| df_means$Organ=='Branch'| df_means$Organ=='Trunk'), ]
df_means_stems$Organ<- factor(df_means_stems$Organ, levels=(c("Twig", "Branch", "Trunk")))

#Analyses----

summary(sma(data = df_means_stems, L_stem_tip~Diam_org, slope.test = 2/3, method = "SMA", log = "XY"))

#Figure----

Fig_S15 <- ggplot(df_means_stems, aes(x=Diam_org, y=L_stem_tip)) +
  geom_point(size=2, alpha=0.5, aes(color=Organ)) +  
  stat_ma_line(method = "SMA", aes(color=Organ)) + 
  scale_color_manual(values = c("#f0f921", "#cc4778", "#0d0887"))+
  scale_y_log10(limits = c(0.00001,100), breaks=c(0.0001, 0.01, 1, 100), labels = c(0.0001, 0.01, 1, 100)) +
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1), labels = c(0.001, 0.01, 0.1, 1)) +
  xlab("Stem diameter (m)") +
  ylab("Distance from stem tip (m)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    legend.position = "none",
    aspect.ratio = 1
  )

ggsave("Figure_S15.JPG", plot=Fig_S15, width = 12, height = 10, dpi = 600, units = "cm")