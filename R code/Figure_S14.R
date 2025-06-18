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

#Subset root data

df_means_roots <- df_means[which(df_means$Organ=='Coarse root'| df_means$Organ=='Fine root'| df_means$Organ=='Very fine root'), ]
df_means_roots$Organ<- factor(df_means_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

df_full_roots <- df_full[which(df_full$Organ=='Coarse root'| df_full$Organ=='Fine root'| df_full$Organ=='Very fine root'), ]
df_full_roots$Organ<- factor(df_full_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

#Analyses----

#Define the number of bins you want
num_bins_root <- 11 #This results in n>100 in each bin

#Calculate the logarithmic range of the stem
min_root <- log10(min(df_full_roots$Diam_org))
max_root <- round(log10(max(df_full_roots$Diam_org)), 2)

#Calculate the bin size
bin_size <- (max_root - min_root) / num_bins_root

#Generate the bins
log_bins_root <- seq(min_root, max_root, by = bin_size)

#Transform to original scale
bins_root <- 10^(log_bins_root)

#Determine minimum sample size across bins

summary_root <- df_full_roots %>%
  mutate(root_class = cut(Diam_org, breaks = bins_root, right = FALSE, include.lowest = TRUE, labels = paste(bins_root[-length(bins_root)], "-", bins_root[-1]))) %>%
  group_by(root_class) %>%
  summarise(original_count = n())

min_root_count<- min(summary_root$original_count)

####Add path classes to full dataset----

df_full_roots <- df_full_roots %>%
  mutate(root_class = cut(Diam_org, breaks = bins_root, right = FALSE, include.lowest = TRUE, labels = paste(bins_root[-length(bins_root)], "-", bins_root[-1])))

####Generate a sample dataset for plotting----
set.seed(999) #Set seed for reproducibility

sampled_data_sample_d_D_roots <- df_full_roots %>%
  group_by(root_class) %>%
  slice_sample(n = min_root_count, replace = FALSE) %>%
  ungroup()

#Figure----

#hydraulic diameter plot

data_plot_means_d_D_roots<- ggplot(df_means_roots, aes(x=Diam_org, y=DhAVG))+
  geom_point(size=1, color="#0d0887", alpha = 0.5)+
  stat_ma_line(method = "SMA", color="black")+
  scale_y_log10(breaks=c(3,10,30), limits=c(1,50))+
  scale_x_log10(breaks=c(0.0003, 0.003, 0.03), labels=c(0.0003, 0.003, 0.03), limits=c(0.0002, 0.03))+
  xlab("Root diameter (m)")+
  ylab("Xylem diameter (µm)")+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size = 14),
    axis.title.x=element_blank(),
    aspect.ratio = 1
  )

#Subsampled plot

sampled_data_plot_d_D_roots<- ggplot(sampled_data_sample_d_D_roots, aes(x=Diam_org, y=DAVG))+
  geom_point(size=1, color="#0d0887", alpha = 0.5)+
  stat_ma_line(method = "SMA", color="black")+
  scale_y_log10(breaks=c(3,10,30), limits=c(1,50))+
  scale_x_log10(breaks=c(0.0003, 0.003, 0.03), labels=c(0.0003, 0.003, 0.03), limits=c(0.0002, 0.03))+
  xlab("Root diameter (m)")+
  ylab("Xylem diameter (µm)")+
  theme_classic()+
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    aspect.ratio = 1)

#Combined plot

Fig_S14<- ggarrange(data_plot_means_d_D_roots, sampled_data_plot_d_D_roots, hjust=-2, vjust = 1.3, labels = c("A", "B"), nrow=1)

#Create an empty plot with the common x-axis label
empty_plot <- ggplot()  + 
  labs(x = "Root diameter (m)")+
  theme(axis.title.x = element_text(size= 15))

Fig_S14_fin <- ggarrange(Fig_S14, empty_plot, ncol = 1, nrow = 2, heights = c(1, 0.06))

ggsave("Figure_S14.JPG", plot=Fig_S14_fin, width = 22, height = 11, dpi = 600, units = "cm")