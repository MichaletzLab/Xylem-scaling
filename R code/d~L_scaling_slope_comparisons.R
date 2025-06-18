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

#Subset aboveground data

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'| df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk"))) 

df_means_aboveground <- df_means[which(df_means$Organ=='Leaf'| df_means$Organ=='Twig'| df_means$Organ=='Branch'| df_means$Organ=='Trunk'), ]
df_means_aboveground$Organ<- factor(df_means_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk")))

#Analyses----

#Define the number of bins you want
num_bins_path <- 8  #This results in each bin having n>100 

#Calculate the logarithmic range of the path
min_path <- log10(min(df_full_aboveground$L))
max_path <- round(log10(max(df_full_aboveground$L)), 2)

#Calculate the bin size
bin_size_path <- (max_path - min_path) / num_bins_path

#Generate the bins
log_bins_path <- seq(min_path, max_path, by = bin_size_path)

#Transform to original scale
bins_path <- 10^(log_bins_path)

####Determine minimum sample size across bins----

summary_path <- df_full_aboveground %>%
  mutate(path_class = cut(L, breaks = bins_path, right = FALSE, include.lowest = TRUE, labels = paste(bins_path[-length(bins_path)], "-", bins_path[-1]))) %>%
  group_by(path_class) %>%
  summarise(original_count = n())

min_path_count<- min(summary_path$original_count)

####Add path classes to full dataset----

df_full_aboveground <- df_full_aboveground %>%
  mutate(path_class = cut(L, breaks = bins_path, right = FALSE, include.lowest = TRUE, labels = paste(bins_path[-length(bins_path)], "-", bins_path[-1])))

####Generate a sample dataset for plotting----

set.seed(999) #Set seed to 999 for reproducibility

sampled_data_sample_d_L <- df_full_aboveground %>%
  group_by(path_class) %>%
  slice_sample(n = min_path_count, replace = FALSE) %>%
  ungroup()

####Fit SMA model to sample dataset----

subsampled_scaling_d_L<- sma(data=sampled_data_sample_d_L, DAVG~L, log = "XY", method = "SMA")

####Determine maximum sample size across bins----

summary_path <- df_full_aboveground %>%
  mutate(path_class = cut(L, breaks = bins_path, right = FALSE, include.lowest = TRUE, labels = paste(bins_path[-length(bins_path)], "-", bins_path[-1]))) %>%
  group_by(path_class) %>%
  summarise(original_count = n())

max_path_count<- max(summary_path$original_count)

####Classify paths into logarithmic bins and determine bootstrap sample size----

summary_path <- df_full_aboveground %>%
  mutate(path_class = cut(L, breaks = bins_path, right = FALSE, include.lowest = TRUE, labels = paste(bins_path[-length(bins_path)], "-", bins_path[-1]))) %>%
  group_by(path_class) %>%
  summarise(original_count = n(),
            bootstrap_count = if_else(n() < max_path_count, max_path_count - n(), 0)) %>%
  ungroup() 

####Add path classes to full dataset----

df_full_aboveground <- df_full_aboveground %>%
  mutate(path_class = cut(L, breaks = bins_path, right = FALSE, include.lowest = TRUE, labels = paste(bins_path[-length(bins_path)], "-", bins_path[-1])))

####Generate sample dataset for plotting----

set.seed(999) #Set seed to "999" for reproducibility

bootstrapped_data_sample_d_L <- df_full_aboveground %>%
  group_by(path_class) %>%
  do({
    #Get the bootstrap count for the current group
    current_path_class <- unique(.$path_class)
    current_bootstrap_count <- summary_path$bootstrap_count[summary_path$path_class == current_path_class]
    
    #Create a flag for the original data as FALSE
    original_data <- mutate(., Resampled = FALSE)
    
    #Bind the original data with the bootstrapped samples
    #and flag the bootstrapped data as TRUE
    bind_rows(
      original_data,
      mutate(slice_sample(., n = current_bootstrap_count, replace = TRUE), Resampled = TRUE)
    )
  }) %>%
  ungroup()

####Fit model to bootstrapped dataset----

d_L_mod_boot<- sma(DAVG~L, data=bootstrapped_data_sample_d_L, method = "SMA", log = "XY")

####Test to see if slopes differ between means, sampled, and bootstrapped datasets----

#First, clean up datasets and prepare to merge into one

#Means dataset

df_means_aboveground_clean<- df_means_aboveground %>% dplyr::select(L, DhAVG)
df_means_aboveground_clean <- rename(df_means_aboveground_clean, DAVG = DhAVG)
df_means_aboveground_clean <- df_means_aboveground_clean %>%
  mutate(Dataset = "Means")

#Subsampled dataset

sampled_data_sample_d_L_clean<- sampled_data_sample_d_L %>% dplyr::select(L, DAVG)
sampled_data_sample_d_L_clean <- sampled_data_sample_d_L_clean %>%
  mutate(Dataset = "Subsampled")

#Bootstrapped dataset

bootstrapped_data_sample_d_L_clean<- bootstrapped_data_sample_d_L %>% dplyr::select(L, DAVG)
bootstrapped_data_sample_d_L_clean <- bootstrapped_data_sample_d_L_clean %>%
  mutate(Dataset = "Bootstrapped")

#Merge into one dataset

d_L_means_sampled_bootstrapped<- rbind(df_means_aboveground_clean, sampled_data_sample_d_L_clean, bootstrapped_data_sample_d_L_clean)
d_L_means_sampled_bootstrapped$Dataset<- as.factor(d_L_means_sampled_bootstrapped$Dataset)

d_L_means_sampled_bootstrapped<- as.data.frame(d_L_means_sampled_bootstrapped)

#Test for common slope

d_L_fits_common_slope<- sma(DAVG~L*Dataset, data = d_L_means_sampled_bootstrapped, log = "XY", method = "SMA", multcomp = T)

summary(d_L_fits_common_slope)