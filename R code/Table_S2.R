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
library(ggrastr)

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

df_means_aboveground <- df_means[which(df_means$Organ=='Leaf'| df_means$Organ=='Twig'| df_means$Organ=='Branch'| df_means$Organ=='Trunk'), ]
df_means_aboveground$Organ<- factor(df_means_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk")))

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'| df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk"))) 

df_full_stems <- df_full[which(df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_stems$Organ<- factor(df_full_stems$Organ, levels=(c("Twig", "Branch", "Trunk"))) 

df_means_stems <- df_means[which(df_means$Organ=='Twig'| df_means$Organ=='Branch'| df_means$Organ=='Trunk'), ]
df_means_stems$Organ<- factor(df_means_stems$Organ, levels=(c("Twig", "Branch", "Trunk")))

df_means_roots <- df_means[which(df_means$Organ=='Coarse root'| df_means$Organ=='Fine root'| df_means$Organ=='Very fine root'), ]
df_means_roots$Organ<- factor(df_means_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

df_full_roots <- df_full[which(df_full$Organ=='Coarse root'| df_full$Organ=='Fine root'| df_full$Organ=='Very fine root'), ]
df_full_roots$Organ<- factor(df_full_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

#Analyses----

##d~L scaling----

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

##d_stem~D_stem scaling----

#Define the number of bins you want
num_bins_stem <- 15 #This results in n>100 in each bin

#Calculate the logarithmic range of the stem
min_stem <- log10(min(df_full_stems$Diam_org))
max_stem <- round(log10(max(df_full_stems$Diam_org)), 2)

#Calculate the bin size
bin_size_stem <- (max_stem - min_stem) / num_bins_stem

#Generate the bins
log_bins_stem <- seq(min_stem, max_stem, by = bin_size_stem)

#Transform to original scale
bins_stem <- 10^(log_bins_stem)

#Determine minimum sample size across bins

summary_stem <- df_full_stems %>%
  mutate(org_class = cut(Diam_org, breaks = bins_stem, right = FALSE, include.lowest = TRUE, labels = paste(bins_stem[-length(bins_stem)], "-", bins_stem[-1]))) %>%
  group_by(org_class) %>%
  summarise(original_count = n())

min_stem_count<- min(summary_stem$original_count)

####Add path classes to full dataset----

df_full_stems <- df_full_stems %>%
  mutate(org_class = cut(Diam_org, breaks = bins_stem, right = FALSE, include.lowest = TRUE, labels = paste(bins_stem[-length(bins_stem)], "-", bins_stem[-1])))

####Generate a sample dataset for plotting----
set.seed(999) #Set seed for reproducibility

sampled_data_sample_d_D_stem <- df_full_stems %>%
  group_by(org_class) %>%
  slice_sample(n = min_stem_count, replace = FALSE) %>%
  ungroup()

#Determine maximum sample size across bins

summary_org <- df_full_stems %>%
  mutate(org_class = cut(Diam_org, breaks = bins_stem, right = FALSE, include.lowest = TRUE, labels = paste(bins_stem[-length(bins_stem)], "-", bins_stem[-1]))) %>%
  group_by(org_class) %>%
  summarise(original_count = n())

max_stem_count<- max(summary_org$original_count)

#Classify organs into logarithmic bins and determine bootstrap sample size

summary_org <- df_full_stems %>%
  mutate(org_class = cut(Diam_org, breaks = bins_stem, right = FALSE, include.lowest = TRUE, labels = paste(bins_stem[-length(bins_stem)], "-", bins_stem[-1]))) %>%
  group_by(org_class) %>%
  summarise(original_count = n(),
            bootstrap_count = if_else(n() < max_stem_count, max_stem_count - n(), 0)) %>%
  ungroup() 

#Add organ classes to full dataset

df_full_stems <- df_full_stems %>%
  mutate(org_class = cut(Diam_org, breaks = bins_stem, right = FALSE, include.lowest = TRUE, labels = paste(bins_stem[-length(bins_stem)], "-", bins_stem[-1])))

#Generate sample dataset for plotting

set.seed(999) #Set seed to "999" for reproducibility

bootstrapped_data_sample_d_D_stem <- df_full_stems %>%
  group_by(org_class) %>%
  do({
    #Get the bootstrap count for the current group
    current_org_class <- unique(.$org_class)
    current_bootstrap_count <- summary_org$bootstrap_count[summary_org$org_class == current_org_class]
    
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

##d_root~D_root scaling----

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

##SMA & OLS fits for all key scaling relationships----

###d ~ L (OLS)----

#d ~ L (bootstrapped)

sma(data=bootstrapped_data_sample_d_L, DAVG~L, log = "XY", method = "SMA")
sma(data=bootstrapped_data_sample_d_L, DAVG~L, log = "XY", method = "OLS")

#d~L (subsampled)

sma(data=sampled_data_sample_d_L, DAVG~L, log = "XY", method = "SMA")
sma(data=sampled_data_sample_d_L, DAVG~L, log = "XY", method = "OLS")

#dh~l

sma(data=df_means_aboveground, DhAVG~L, log = "XY", method = "SMA")
sma(data=df_means_aboveground, DhAVG~L, log = "XY", method = "OLS")

###d ~ Dstem (OLS)----

#d ~ D (bootstrapped)

sma(data=bootstrapped_data_sample_d_D_stem, DAVG~Diam_org, log = "XY", method = "SMA")
sma(data=bootstrapped_data_sample_d_D_stem, DAVG~Diam_org, log = "XY", method = "OLS")

#d ~ D (subsampled)

sma(data=sampled_data_sample_d_D_stem, DAVG~Diam_org, log = "XY", method = "SMA")
sma(data=sampled_data_sample_d_D_stem, DAVG~Diam_org, log = "XY", method = "OLS")

#dh ~ D

sma(data=df_means_stems, DhAVG~Diam_org, log = "XY", method = "SMA")
sma(data=df_means_stems, DhAVG~Diam_org, log = "XY", method = "OLS")

###d ~ Droot (OLS)----

#d ~ D (bootstrapped)

sma(data=bootstrapped_data_sample_d_D_roots, DAVG~Diam_org, log = "XY", method = "SMA")
sma(data=bootstrapped_data_sample_d_D_roots, DAVG~Diam_org, log = "XY", method = "OLS")

#d ~ D (subsampled)

sma(data=sampled_data_sample_d_D_roots, DAVG~Diam_org, log = "XY", method = "SMA")
sma(data=sampled_data_sample_d_D_roots, DAVG~Diam_org, log = "XY", method = "OLS")

#dh ~ D

sma(data=df_means_roots, DhAVG~Diam_org, log = "XY", method = "SMA")
sma(data=df_means_roots, DhAVG~Diam_org, log = "XY", method = "OLS")