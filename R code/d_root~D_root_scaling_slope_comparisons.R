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

df_full_roots <- df_full[which(df_full$Organ=='Coarse root'| df_full$Organ=='Fine root'| df_full$Organ=='Very fine root'), ]
df_full_roots$Organ<- factor(df_full_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

df_means_roots <- df_means[which(df_means$Organ=='Coarse root'| df_means$Organ=='Fine root'| df_means$Organ=='Very fine root'), ]
df_means_roots$Organ<- factor(df_means_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

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

bins_root

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

#Determine maximum sample size across bins

summary_root <- df_full_roots %>%
  mutate(root_class = cut(Diam_org, breaks = bins_root, right = FALSE, include.lowest = TRUE, labels = paste(bins_root[-length(bins_root)], "-", bins_root[-1]))) %>%
  group_by(root_class) %>%
  summarise(original_count = n())

max_root_count<- max(summary_root$original_count)

#Classify organs into logarithmic bins and determine bootstrap sample size

summary_root <- df_full_roots %>%
  mutate(root_class = cut(Diam_org, breaks = bins_root, right = FALSE, include.lowest = TRUE, labels = paste(bins_root[-length(bins_root)], "-", bins_root[-1]))) %>%
  group_by(root_class) %>%
  summarise(original_count = n(),
            bootstrap_count = if_else(n() < max_root_count, max_root_count - n(), 0)) %>%
  ungroup() 

#Add organ classes to full dataset

df_full_roots <- df_full_roots %>%
  mutate(root_class = cut(Diam_org, breaks = bins_root, right = FALSE, include.lowest = TRUE, labels = paste(bins_root[-length(bins_root)], "-", bins_root[-1])))

#Generate sample dataset for plotting

set.seed(999) #Set seed to "999" for reproducibility

bootstrapped_data_sample_d_D_roots <- df_full_roots %>%
  group_by(root_class) %>%
  do({
    #Get the bootstrap count for the current group
    current_root_class <- unique(.$root_class)
    current_bootstrap_count <- summary_root$bootstrap_count[summary_root$root_class == current_root_class]
    
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

##Test to see if slopes differ between means, sampled, and bootstrapped datasets----

#First, clean up datasets and prepare to merge into one

#Means dataset

df_means_roots_clean<- df_means_roots %>% dplyr::select(Diam_org, DhAVG)
df_means_roots_clean <- rename(df_means_roots_clean, DAVG = DhAVG)
df_means_roots_clean <- df_means_roots_clean %>%
  mutate(Dataset = "means_roots")

#Sampled w/o replacement dataset

sampled_data_sample_d_D_roots_clean<- sampled_data_sample_d_D_roots %>% dplyr::select(Diam_org, DAVG)
sampled_data_sample_d_D_roots_clean <- sampled_data_sample_d_D_roots_clean %>%
  mutate(Dataset = "Subsampled")

#Bootstrapped dataset

bootstrapped_data_sample_d_D_roots_clean<- bootstrapped_data_sample_d_D_roots %>% dplyr::select(Diam_org, DAVG)
bootstrapped_data_sample_d_D_roots_clean <- bootstrapped_data_sample_d_D_roots_clean %>%
  mutate(Dataset = "Bootstrapped")

#Merge into one dataset

d_D_roots_means_sampled_bootstrapped<- rbind(df_means_roots_clean, sampled_data_sample_d_D_roots_clean, bootstrapped_data_sample_d_D_roots_clean)
d_D_roots_means_sampled_bootstrapped$Dataset<- as.factor(d_D_roots_means_sampled_bootstrapped$Dataset)

d_D_roots_means_sampled_bootstrapped<- as.data.frame(d_D_roots_means_sampled_bootstrapped)

#Test for common slope

d_D_roots_fits_common_slope<- sma(DAVG~Diam_org*Dataset, data = d_D_roots_means_sampled_bootstrapped, log = "XY", method = "SMA", multcomp = T)

summary(d_D_roots_fits_common_slope)