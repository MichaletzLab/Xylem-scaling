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

#Subset stem data

df_full_stems <- df_full[which(df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_stems$Organ<- factor(df_full_stems$Organ, levels=(c("Twig", "Branch", "Trunk"))) 

df_means_stems <- df_means[which(df_means$Organ=='Twig'| df_means$Organ=='Branch'| df_means$Organ=='Trunk'), ]
df_means_stems$Organ<- factor(df_means_stems$Organ, levels=(c("Twig", "Branch", "Trunk")))

#Analyses----

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

##Test to see if slopes differ between means, sampled, and bootstrapped datasets----

#First, clean up datasets and prepare to merge into one

#Means dataset

df_means_stems_clean<- df_means_stems %>% dplyr::select(Diam_org, DhAVG)
df_means_stems_clean <- rename(df_means_stems_clean, DAVG = DhAVG)
df_means_stems_clean <- df_means_stems_clean %>%
  mutate(Dataset = "means_stems")

#Sampled w/o replacement dataset

sampled_data_sample_d_D_stem_clean<- sampled_data_sample_d_D_stem %>% dplyr::select(Diam_org, DAVG)
sampled_data_sample_d_D_stem_clean <- sampled_data_sample_d_D_stem_clean %>%
  mutate(Dataset = "Subsampled")

#Bootstrapped dataset

bootstrapped_data_sample_d_D_stem_clean<- bootstrapped_data_sample_d_D_stem %>% dplyr::select(Diam_org, DAVG)
bootstrapped_data_sample_d_D_stem_clean <- bootstrapped_data_sample_d_D_stem_clean %>%
  mutate(Dataset = "Bootstrapped")

#Merge into one dataset

d_D_stem_means_sampled_bootstrapped<- rbind(df_means_stems_clean, sampled_data_sample_d_D_stem_clean, bootstrapped_data_sample_d_D_stem_clean)
d_D_stem_means_sampled_bootstrapped$Dataset<- as.factor(d_D_stem_means_sampled_bootstrapped$Dataset)

d_D_stem_means_sampled_bootstrapped<- as.data.frame(d_D_stem_means_sampled_bootstrapped)

#Test for common slope

d_D_stem_fits_common_slope<- sma(DAVG~Diam_org*Dataset, data = d_D_stem_means_sampled_bootstrapped, log = "XY", method = "SMA", multcomp = T)

summary(d_D_stem_fits_common_slope)
