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

df_full_roots <- df_full[which(df_full$Organ=='Coarse root'| df_full$Organ=='Fine root'| df_full$Organ=='Very fine root'), ]
df_full_roots$Organ<- factor(df_full_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

#Analyses----

###First, calculate the bootstrap sample size for each individual and L----

#Create new dataframe
df_full_aboveground_d_L_ind<- df_full_aboveground %>%
  dplyr::select(Species, Individual, ID, Organ, L, DAVG)

#Since MoHe1 is missing core data, substitute coarse root data

MoHe1_cr_data <- df_full_roots[df_full_roots$ID == "MoHe1_cr", ]

MoHe1_cr_data <- MoHe1_cr_data %>%
  dplyr::select(Species, Individual, ID, Organ, L, DAVG)

MoHe1_cr_data$L<- 5.574

df_full_aboveground_d_L_ind<-  rbind(df_full_aboveground_d_L_ind, MoHe1_cr_data)

#Exclude MoHe1 due to missing stem core data

df_full_aboveground_d_L_ind <- df_full_aboveground_d_L_ind %>%
  filter(Individual != "MoHe1")

#Group by Individual and L, then count the number of DAVG observations
df_counts_d_L <- df_full_aboveground_d_L_ind %>%
  group_by(Individual, L) %>%
  summarise(count_DAVG = n()) %>%
  ungroup()

#For each Individual, find the L with the greatest number of DAVG observations
df_max_counts_d_L <- df_counts_d_L %>%
  group_by(Individual) %>%
  filter(count_DAVG == max(count_DAVG)) %>%
  ungroup()

#Join the max counts back to the original counts dataframe
df_diff_d_L <- df_counts_d_L %>%
  left_join(df_max_counts_d_L, by = "Individual", suffix = c("", "_max")) %>%
  mutate(bootstrap_count_d_L = count_DAVG_max - count_DAVG) %>%
  dplyr::select(Individual, L, bootstrap_count_d_L)

#Merge counts with full dataset
df_full_aboveground_d_L_ind<-  merge(df_full_aboveground_d_L_ind, df_diff_d_L, by = c("Individual", "L"), all.x = TRUE)

###Create new dataframe containing bootstrapped data for each individual at each L----

#Create a function to perform bootstrapping for each group:
bootstrap_samples_d_L <- function(data) {
  n <- data$bootstrap_count_d_L[1]
  data %>%
    slice_sample(n = n, replace = TRUE)
}

#Bootstrap the data
set.seed(999) #set seed for reproducibility

bootstrapped_data_d_L_ind <- df_full_aboveground_d_L_ind %>%
  group_by(Individual, L) %>%
  group_modify(~ bootstrap_samples_d_L(.x)) %>%
  ungroup()

#Combine the original and bootstrapped data and clean it up
df_full_aboveground_d_L_ind_final <- rbind(df_full_aboveground_d_L_ind, bootstrapped_data_d_L_ind)

df_full_aboveground_d_L_ind_final<- df_full_aboveground_d_L_ind_final %>%
  dplyr::select(Species,Individual, ID, Organ, L, DAVG)

###Fit SMA model to each individual and extract coefficients----

#Fit the model
d_L_ind_model<- sma(data= df_full_aboveground_d_L_ind_final, DAVG~L*Individual, log = "XY")

d_L_ind_model_summary<- as.data.frame(d_L_ind_model$groupsummary)

#Calculate mean slope

mean(d_L_ind_model_summary$Slope)

#Calculate 95% CI

d_L_ind_model_summary$Slope

smean.cl.boot(d_L_ind_model_summary$Slope)

#Fitting SMA models to each species

d_L_species_model<- sma(data= df_full_aboveground_d_L_ind_final, DAVG~L*Species, log = "XY")

summary(d_L_species_model)

d_L_species_model_summary<- as.data.frame(d_L_species_model$groupsummary)

d_L_species_model_summary

#Figure----

df_full_aboveground_d_L_ind_final$Species<- as.factor(df_full_aboveground_d_L_ind_final$Species)

bootstrapped_data_plot_d_L_ind <- ggplot(df_full_aboveground_d_L_ind_final, aes(x=L, y=DAVG)) +
  geom_point(aes(group=Species, color=Species), alpha=0.025) + 
  scale_color_manual(values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F"))+
  stat_ma_line(aes(group=Individual, color= Species), se=FALSE, linewidth=0.5, alpha=1, method = "SMA") +
  scale_y_log10(breaks=c(3,10,30), limits=c(1,80)) +
  scale_x_log10(breaks=c(0.01,0.1, 1, 10), labels =c(0.01,0.1, 1, 10)) +
  xlab("Distance from leaf tip (m)") +
  ylab("Xylem conduit diameter (µm)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size=14),
    axis.title.x = element_text(size=14),
    legend.position = "right",
    aspect.ratio = 1
  )

ggsave("Figure_S12.JPG", plot = bootstrapped_data_plot_d_L_ind, width = 12, height = 10, units = "cm", dpi = 600)