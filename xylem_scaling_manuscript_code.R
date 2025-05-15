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

#Load the required packages

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

##Working directory----

wd<- #INSERT PATH TO WORKING DIRECTORY (TO SAVE FIGURES)

setwd(wd)

##Data----
###Import full dataset----

path_to_file<- #INSERT LOCAL PATH TO XYLEM_SCALING_DATA.XLSX

df_full <- read_excel(path_to_file)

options(scipen = 999) #prevent scientific notation on plots

###Calculate averages as separate dataset----

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

###Subset full dataset----
#Subset aboveground data for means dataset

df_means_aboveground <- df_means[which(df_means$Organ=='Leaf'| df_means$Organ=='Twig'| df_means$Organ=='Branch'| df_means$Organ=='Trunk'), ]
df_means_aboveground$Organ<- factor(df_means_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk")))

#Subset aboveground data for full dataset

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'| df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk"))) 

#Subset stem data 

df_full_stems <- df_full[which(df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_stems$Organ<- factor(df_full_stems$Organ, levels=(c("Twig", "Branch", "Trunk"))) 

df_means_stems <- df_means[which(df_means$Organ=='Twig'| df_means$Organ=='Branch'| df_means$Organ=='Trunk'), ]
df_means_stems$Organ<- factor(df_means_stems$Organ, levels=(c("Twig", "Branch", "Trunk")))

#Subset roots data for means dataset

df_means_roots <- df_means[which(df_means$Organ=='Coarse root'| df_means$Organ=='Fine root'| df_means$Organ=='Very fine root'), ]
df_means_roots$Organ<- factor(df_means_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

#Subset roots data for full dataset

df_full_roots <- df_full[which(df_full$Organ=='Coarse root'| df_full$Organ=='Fine root'| df_full$Organ=='Very fine root'), ]
df_full_roots$Organ<- factor(df_full_roots$Organ, levels=(c('Coarse root', 'Fine root', 'Very fine root')))

#Clean up datasets

df_full_clean<- df_full %>% dplyr::select(Diam_org, L, DAVG)
df_full_aboveground_clean <- df_full_aboveground %>% dplyr::select(Diam_org, L, DAVG)
df_full_roots_clean<- df_full_roots %>% dplyr::select(Diam_org, L, DAVG)

#Analyses----
#Differences in conduit diameter and t/b^2 between organs (Fig. 1)----

#Reorganize organs in correct order

Organ_order<- rev(c("Leaf", "Twig", "Branch", "Trunk", "Coarse root", "Fine root", "Very fine root"))

df_full$Organ<- factor(df_full$Organ, levels = Organ_order)

df_full$Organ_order_num<- as.numeric(factor(df_full$Organ, levels=c("Leaf", "Twig", "Branch", "Trunk", "Coarse root", "Fine root", "Very fine root")))

df_means$Organ<- factor(df_means$Organ, levels = Organ_order)

df_means$Organ_order_num<- as.numeric(factor(df_means$Organ, levels=c("Leaf", "Twig", "Branch", "Trunk", "Coarse root", "Fine root", "Very fine root")))

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

##Segmented regression d~Organ

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

##Segmented regression (t/b)^2~Organ

ANOVA_FIG1_B<- aov(log10(t_b_2)~Organ_order_num, data = df_full)

ANOVA_seg_FIG1_B <- segmented(ANOVA_FIG1_B, seg.Z = ~Organ_order_num)

summary(ANOVA_seg_FIG1_B)

#d~L scaling----

##Fit models to empirical (i.e., not resampled) data (d~L)----

###Fit SMA model to aboveground organs averaged for each organ (or sample)----

mean_scaling_d_L<- sma(data=df_means_aboveground, DhAVG~L, log = "XY", method = "SMA")

mean_scaling_d_L

###Fit SMA model to aboveground organs in full dataset----

full_scaling_d_L<- sma(data=df_full_aboveground, DAVG~L, log = "XY", method = "SMA")

summary(full_scaling_d_L)

##Fitting models using resampled data (d~L)----

###Fitting models using the subsampling method (d~L)----

####Define log bins for path length-----

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

bins_path

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

subsampled_scaling_d_L

####Function to perform sampling, fit SMA model, and extract coefficients----
run_iteration <- function(data) {
  #Random sampling of d across each L bin
  sampled_data <- data %>%
    group_by(path_class) %>%
    slice_sample(n = min_path_count, replace = FALSE) %>%
    ungroup()
  
  #Fit a SMA model to sampled_data
  sampled_scaling <- sma(data = sampled_data, DAVG ~ L, log = "XY", method = "SMA")
  
  #Extract coefficients
  sampled_scaling_summary <- sampled_scaling$groupsummary
  sampled_slope <- sampled_scaling_summary$Slope
  sampled_y_int <- sampled_scaling_summary$Int
  sampled_slope_lowCI <- sampled_scaling_summary$Slope_lowCI
  sampled_slope_highCI <- sampled_scaling_summary$Slope_highCI
  
  #Return a list of coefficients
  return(list(
    slope = sampled_slope,
    y_intercept = sampled_y_int,
    slope_lowCI = sampled_slope_lowCI,
    slope_highCI = sampled_slope_highCI
  ))
}

####Specify the numbers of iterations----
num_iterations_list <- c(1, 3, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000)

####Create an empty dataframe to store the results----
sampled_results_df <- data.frame(
  iteration = integer(),
  slope = double(),
  y_intercept = double(),
  slope_lowCI = double(),
  slope_highCI = double(),
  num_iterations = integer()
)

####Loop through different numbers of iterations----
for (num_iterations in num_iterations_list) {
  #Loop through iterations
  for (i in 1:num_iterations) {
    #Set seed for this iteration based on iteration index
    set.seed(i)
    #Run the iteration and extract coefficients
    iteration_results <- run_iteration(df_full_aboveground)
    
    #Store the results in the dataframe
    sampled_results_df <- rbind(sampled_results_df, c(i, iteration_results$slope, iteration_results$y_intercept, iteration_results$slope_lowCI, iteration_results$slope_highCI, num_iterations))
  }
}

####Rename the columns in sampled_results_df----
colnames(sampled_results_df) <- c("iteration", "slope", "y_intercept", "slope_lowCI", "slope_highCI", "num_iterations")

####Calculate mean and standard error for each variable at each number of iterations----
sampled_summary_df <- sampled_results_df %>%
  group_by(num_iterations) %>%
  dplyr::summarize(
    mean_slope = mean(slope),
    mean_y_intercept = mean(y_intercept),
    mean_slope_lowCI = mean(slope_lowCI),
    mean_slope_highCI = mean(slope_highCI),
    se_slope = sd(slope) / sqrt(n()),
    se_y_intercept = sd(y_intercept) / sqrt(n()),
    se_slope_lowCI = sd(slope_lowCI) / sqrt(n()),
    se_slope_highCI = sd(slope_highCI) / sqrt(n()),
    CV_slope = (sd(slope)/mean(slope))*100,
    CV_highCI = (sd(slope_highCI)/mean(slope_highCI))*100,
    CV_lowCI = (sd(slope_lowCI)/mean(slope_lowCI))*100
  )

####Rename the sampled dataset & stat summary and store----

sampled_results_df_d_L<- sampled_results_df
sampled_summary_df_d_L<- sampled_summary_df

###Fitting models using the bootstrapping method (d~L)----

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

d_L_mod_boot

#Function to perform sampling, fit SMA model, and extract coefficients

run_iteration <- function(data) {
  
  #Perform the bootstrapping method
  bootstrapped_data <- data %>%
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
  
  #Fit a SMA model to bootstrapped_data
  bootstrapped_scaling <- sma(data = bootstrapped_data, DAVG ~ L, log = "XY", method = "SMA")
  
  #Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_slope <- bootstrapped_scaling_summary$Slope
  bootstrapped_y_int <- bootstrapped_scaling_summary$Int
  bootstrapped_slope_lowCI <- bootstrapped_scaling_summary$Slope_lowCI
  bootstrapped_slope_highCI <- bootstrapped_scaling_summary$Slope_highCI
  
  #Return a list of coefficients
  return(list(
    slope = bootstrapped_slope,
    y_intercept = bootstrapped_y_int,
    slope_lowCI = bootstrapped_slope_lowCI,
    slope_highCI = bootstrapped_slope_highCI
  ))
}

####Specify the numbers of iterations----
num_iterations_list <- c(1,3,5,10,50,100,500,1000)

####Create an empty dataframe to store the results----
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  slope = double(),
  y_intercept = double(),
  slope_lowCI = double(),
  slope_highCI = double(),
  num_iterations = integer()
)

####Loop through different numbers of iterations----
for (num_iterations in num_iterations_list) {
  #Loop through iterations
  for (i in 1:num_iterations) {
    #Set seed for this iteration based on iteration index
    set.seed(i)
    #Run the iteration and extract coefficients
    iteration_results <- run_iteration(df_full_aboveground)
    
    #Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(i, iteration_results$slope, iteration_results$y_intercept, iteration_results$slope_lowCI, iteration_results$slope_highCI, num_iterations))
  }
}

####Rename the columns in bootstrapped_results_df----
colnames(bootstrapped_results_df) <- c("iteration", "slope", "y_intercept", "slope_lowCI", "slope_highCI", "num_iterations")

####Calculate mean and standard error for each variable at each number of iterations----
bootstrapped_summary_df <- bootstrapped_results_df %>%
  group_by(num_iterations) %>%
  dplyr::summarize(
    mean_slope = mean(slope),
    mean_y_intercept = mean(y_intercept),
    mean_slope_lowCI = mean(slope_lowCI),
    mean_slope_highCI = mean(slope_highCI),
    se_slope = sd(slope) / sqrt(n()),
    se_y_intercept = sd(y_intercept) / sqrt(n()),
    se_slope_lowCI = sd(slope_lowCI) / sqrt(n()),
    se_slope_highCI = sd(slope_highCI) / sqrt(n()),
    CV_slope = (sd(slope)/mean(slope))*100,
    CV_highCI = (sd(slope_highCI)/mean(slope_highCI))*100,
    CV_lowCI = (sd(slope_lowCI)/mean(slope_lowCI))*100
  )

####Rename the bootstrapped dataset & stat summary and store----

bootstrapped_results_df_d_L<- bootstrapped_results_df
bootstrapped_summary_df_d_L<- bootstrapped_summary_df

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

#d~D scaling in stems----

##Fit models to empirical (i.e., not resampled) data (d~D)----

###Fit SMA model to data averaged for each sample----

mean_scaling<- sma(data=df_means_stems, DhAVG~Diam_org, log = "XY", method = "SMA")

mean_scaling

###Fit SMA model to stems organs in full dataset----

full_scaling<- sma(data=df_full_stems, DAVG~Diam_org, log = "XY", method = "SMA")

full_scaling

###Fitting models using the subsampling method (D~d)----

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

bins_stem

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

####Fit SMA model to sample dataset----

sma(data=sampled_data_sample_d_D_stem, DAVG~Diam_org, log = "XY", method = "SMA")

####Function to perform sampling, fit SMA model, and extract coefficients----
run_iteration <- function(data) {
  #Random sampling of d across each D bin
  sampled_data <- data %>%
    group_by(org_class) %>%
    slice_sample(n = min_stem_count, replace = FALSE) %>%
    ungroup()
  
  #Fit a SMA model to sampled_data
  sampled_scaling <- sma(data = sampled_data, DAVG ~ Diam_org, log = "XY", method = "SMA")
  
  #Extract coefficients
  sampled_scaling_summary <- sampled_scaling$groupsummary
  sampled_slope <- sampled_scaling_summary$Slope
  sampled_y_int <- sampled_scaling_summary$Int
  sampled_slope_lowCI <- sampled_scaling_summary$Slope_lowCI
  sampled_slope_highCI <- sampled_scaling_summary$Slope_highCI
  
  #Return a list of coefficients
  return(list(
    slope = sampled_slope,
    y_intercept = sampled_y_int,
    slope_lowCI = sampled_slope_lowCI,
    slope_highCI = sampled_slope_highCI
  ))
}

####Specify the numbers of iterations----
num_iterations_list <- c(1, 3, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000)

####Create an empty dataframe to store the results----
sampled_results_df <- data.frame(
  iteration = integer(),
  slope = double(),
  y_intercept = double(),
  slope_lowCI = double(),
  slope_highCI = double(),
  num_iterations = integer()
)

####Loop through different numbers of iterations----
for (num_iterations in num_iterations_list) {
  #Loop through iterations
  for (i in 1:num_iterations) {
    #Set seed for this iteration based on iteration index
    set.seed(i)
    #Run the iteration and extract coefficients
    iteration_results <- run_iteration(df_full_stems)
    
    #Store the results in the dataframe
    sampled_results_df <- rbind(sampled_results_df, c(i, iteration_results$slope, iteration_results$y_intercept, iteration_results$slope_lowCI, iteration_results$slope_highCI, num_iterations))
  }
}

####Rename the columns in sampled_results_df----
colnames(sampled_results_df) <- c("iteration", "slope", "y_intercept", "slope_lowCI", "slope_highCI", "num_iterations")

####Calculate mean and standard error for each variable at each number of iterations----
sampled_summary_df <- sampled_results_df %>%
  group_by(num_iterations) %>%
  dplyr::summarize(
    mean_slope = mean(slope),
    mean_y_intercept = mean(y_intercept),
    mean_slope_lowCI = mean(slope_lowCI),
    mean_slope_highCI = mean(slope_highCI),
    se_slope = sd(slope) / sqrt(n()),
    se_y_intercept = sd(y_intercept) / sqrt(n()),
    se_slope_lowCI = sd(slope_lowCI) / sqrt(n()),
    se_slope_highCI = sd(slope_highCI) / sqrt(n()),
    CV_slope = (sd(slope)/mean(slope))*100,
    CV_highCI = (sd(slope_highCI)/mean(slope_highCI))*100,
    CV_lowCI = (sd(slope_lowCI)/mean(slope_lowCI))*100
  )

####Rename the sampled dataset & stat summary and store----

sampled_results_df_d_D_stem<- sampled_results_df
sampled_summary_df_d_D_stem<- sampled_summary_df

###Fitting models using the bootstrapping method (d~D)----

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

#Fit model to bootstrapped dataset

sma(DAVG~Diam_org, data=bootstrapped_data_sample_d_D_stem, method = "SMA", log = "XY")  

#Function to perform sampling, fit SMA model, and extract coefficients

run_iteration <- function(data) {
  
  #Perform the bootstrapping method
  bootstrapped_data <- data %>%
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
  
  #Fit a SMA model to bootstrapped_data
  bootstrapped_scaling <- sma(data = bootstrapped_data, DAVG ~ Diam_org, log = "XY", method = "SMA")
  
  #Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_slope <- bootstrapped_scaling_summary$Slope
  bootstrapped_y_int <- bootstrapped_scaling_summary$Int
  bootstrapped_slope_lowCI <- bootstrapped_scaling_summary$Slope_lowCI
  bootstrapped_slope_highCI <- bootstrapped_scaling_summary$Slope_highCI
  
  #Return a list of coefficients
  return(list(
    slope = bootstrapped_slope,
    y_intercept = bootstrapped_y_int,
    slope_lowCI = bootstrapped_slope_lowCI,
    slope_highCI = bootstrapped_slope_highCI
  ))
}

####Specify the numbers of iterations----
num_iterations_list <- c(1,3,5,10,50,100,500,1000)

####Create an empty dataframe to store the results----
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  slope = double(),
  y_intercept = double(),
  slope_lowCI = double(),
  slope_highCI = double(),
  num_iterations = integer()
)

####Loop through different numbers of iterations----
for (num_iterations in num_iterations_list) {
  #Loop through iterations
  for (i in 1:num_iterations) {
    #Set seed for this iteration based on iteration index
    set.seed(i)
    #Run the iteration and extract coefficients
    iteration_results <- run_iteration(df_full_stems)
    
    #Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(i, iteration_results$slope, iteration_results$y_intercept, iteration_results$slope_lowCI, iteration_results$slope_highCI, num_iterations))
  }
}

####Rename the columns in bootstrapped_results_df----
colnames(bootstrapped_results_df) <- c("iteration", "slope", "y_intercept", "slope_lowCI", "slope_highCI", "num_iterations")

####Calculate mean and standard error for each variable at each number of iterations----
bootstrapped_summary_df <- bootstrapped_results_df %>%
  group_by(num_iterations) %>%
  dplyr::summarize(
    mean_slope = mean(slope),
    mean_y_intercept = mean(y_intercept),
    mean_slope_lowCI = mean(slope_lowCI),
    mean_slope_highCI = mean(slope_highCI),
    se_slope = sd(slope) / sqrt(n()),
    se_y_intercept = sd(y_intercept) / sqrt(n()),
    se_slope_lowCI = sd(slope_lowCI) / sqrt(n()),
    se_slope_highCI = sd(slope_highCI) / sqrt(n()),
    CV_slope = (sd(slope)/mean(slope))*100,
    CV_highCI = (sd(slope_highCI)/mean(slope_highCI))*100,
    CV_lowCI = (sd(slope_lowCI)/mean(slope_lowCI))*100
  )

####Rename the bootstrapped dataset & stat summary and store----

bootstrapped_results_df_d_D_stem<- bootstrapped_results_df
bootstrapped_summary_df_d_D_stem<- bootstrapped_summary_df

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

#d~D scaling in roots----

##Fit models to empirical (i.e., not resampled) data (d~D)----

###Fit SMA model to data averaged for each sample----

mean_scaling<- sma(data=df_means_roots, DhAVG~Diam_org, log = "XY", method = "SMA")

mean_scaling

###Fit SMA model to roots organs in full dataset----

full_scaling<- sma(data=df_full_roots, DAVG~Diam_org, log = "XY", method = "SMA")

full_scaling

###Fitting models using the subsampling method (D~d)----

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

####Fit SMA model to sample dataset----

sma(data=sampled_data_sample_d_D_roots, DAVG~Diam_org, log = "XY", method = "SMA")

####Function to perform sampling, fit SMA model, and extract coefficients----
run_iteration <- function(data) {
  #Random sampling of d across each D bin
  sampled_data <- data %>%
    group_by(root_class) %>%
    slice_sample(n = min_root_count, replace = FALSE) %>%
    ungroup()
  
  #Fit a SMA model to sampled_data
  sampled_scaling <- sma(data = sampled_data, DAVG ~ Diam_org, log = "XY", method = "SMA")
  
  #Extract coefficients
  sampled_scaling_summary <- sampled_scaling$groupsummary
  sampled_slope <- sampled_scaling_summary$Slope
  sampled_y_int <- sampled_scaling_summary$Int
  sampled_slope_lowCI <- sampled_scaling_summary$Slope_lowCI
  sampled_slope_highCI <- sampled_scaling_summary$Slope_highCI
  
  #Return a list of coefficients
  return(list(
    slope = sampled_slope,
    y_intercept = sampled_y_int,
    slope_lowCI = sampled_slope_lowCI,
    slope_highCI = sampled_slope_highCI
  ))
}

####Specify the numbers of iterations----
num_iterations_list <- c(1, 3, 5, 10, 50, 100, 500, 1000, 5000, 10000, 50000)

####Create an empty dataframe to store the results----
sampled_results_df <- data.frame(
  iteration = integer(),
  slope = double(),
  y_intercept = double(),
  slope_lowCI = double(),
  slope_highCI = double(),
  num_iterations = integer()
)

####Loop through different numbers of iterations----
for (num_iterations in num_iterations_list) {
  #Loop through iterations
  for (i in 1:num_iterations) {
    #Set seed for this iteration based on iteration index
    set.seed(i)
    #Run the iteration and extract coefficients
    iteration_results <- run_iteration(df_full_roots)
    
    #Store the results in the dataframe
    sampled_results_df <- rbind(sampled_results_df, c(i, iteration_results$slope, iteration_results$y_intercept, iteration_results$slope_lowCI, iteration_results$slope_highCI, num_iterations))
  }
}

####Rename the columns in sampled_results_df----
colnames(sampled_results_df) <- c("iteration", "slope", "y_intercept", "slope_lowCI", "slope_highCI", "num_iterations")

####Calculate mean and standard error for each variable at each number of iterations----
sampled_summary_df <- sampled_results_df %>%
  group_by(num_iterations) %>%
  dplyr::summarize(
    mean_slope = mean(slope),
    mean_y_intercept = mean(y_intercept),
    mean_slope_lowCI = mean(slope_lowCI),
    mean_slope_highCI = mean(slope_highCI),
    se_slope = sd(slope) / sqrt(n()),
    se_y_intercept = sd(y_intercept) / sqrt(n()),
    se_slope_lowCI = sd(slope_lowCI) / sqrt(n()),
    se_slope_highCI = sd(slope_highCI) / sqrt(n()),
    CV_slope = (sd(slope)/mean(slope))*100,
    CV_highCI = (sd(slope_highCI)/mean(slope_highCI))*100,
    CV_lowCI = (sd(slope_lowCI)/mean(slope_lowCI))*100
  )

####Rename the sampled dataset & stat summary and store----

sampled_results_df_d_D_roots<- sampled_results_df
sampled_summary_df_d_D_roots<- sampled_summary_df

###Fitting models using the bootstrapping method (d~D)----

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

#Fit model to bootstrapped dataset

sma(DAVG~Diam_org, data=bootstrapped_data_sample_d_D_roots, method = "SMA", log = "XY")  

#Function to perform sampling, fit SMA model, and extract coefficients

run_iteration <- function(data) {
  
  #Perform the bootstrapping method
  bootstrapped_data <- data %>%
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
  
  #Fit a SMA model to bootstrapped_data
  bootstrapped_scaling <- sma(data = bootstrapped_data, DAVG ~ Diam_org, log = "XY", method = "SMA")
  
  #Extract coefficients
  bootstrapped_scaling_summary <- bootstrapped_scaling$groupsummary
  bootstrapped_slope <- bootstrapped_scaling_summary$Slope
  bootstrapped_y_int <- bootstrapped_scaling_summary$Int
  bootstrapped_slope_lowCI <- bootstrapped_scaling_summary$Slope_lowCI
  bootstrapped_slope_highCI <- bootstrapped_scaling_summary$Slope_highCI
  
  #Return a list of coefficients
  return(list(
    slope = bootstrapped_slope,
    y_intercept = bootstrapped_y_int,
    slope_lowCI = bootstrapped_slope_lowCI,
    slope_highCI = bootstrapped_slope_highCI
  ))
}

####Specify the numbers of iterations----
num_iterations_list <- c(1,3,5,10,50,100,500,1000)

####Create an empty dataframe to store the results----
bootstrapped_results_df <- data.frame(
  iteration = integer(),
  slope = double(),
  y_intercept = double(),
  slope_lowCI = double(),
  slope_highCI = double(),
  num_iterations = integer()
)

####Loop through different numbers of iterations----
for (num_iterations in num_iterations_list) {
  #Loop through iterations
  for (i in 1:num_iterations) {
    #Set seed for this iteration based on iteration index
    set.seed(i)
    #Run the iteration and extract coefficients
    iteration_results <- run_iteration(df_full_roots)
    
    #Store the results in the dataframe
    bootstrapped_results_df <- rbind(bootstrapped_results_df, c(i, iteration_results$slope, iteration_results$y_intercept, iteration_results$slope_lowCI, iteration_results$slope_highCI, num_iterations))
  }
}

####Rename the columns in bootstrapped_results_df----
colnames(bootstrapped_results_df) <- c("iteration", "slope", "y_intercept", "slope_lowCI", "slope_highCI", "num_iterations")

####Calculate mean and standard error for each variable at each number of iterations----
bootstrapped_summary_df <- bootstrapped_results_df %>%
  group_by(num_iterations) %>%
  dplyr::summarize(
    mean_slope = mean(slope),
    mean_y_intercept = mean(y_intercept),
    mean_slope_lowCI = mean(slope_lowCI),
    mean_slope_highCI = mean(slope_highCI),
    se_slope = sd(slope) / sqrt(n()),
    se_y_intercept = sd(y_intercept) / sqrt(n()),
    se_slope_lowCI = sd(slope_lowCI) / sqrt(n()),
    se_slope_highCI = sd(slope_highCI) / sqrt(n()),
    CV_slope = (sd(slope)/mean(slope))*100,
    CV_highCI = (sd(slope_highCI)/mean(slope_highCI))*100,
    CV_lowCI = (sd(slope_lowCI)/mean(slope_lowCI))*100
  )

####Rename the bootstrapped dataset & stat summary and store----

bootstrapped_results_df_d_D_roots<- bootstrapped_results_df
bootstrapped_summary_df_d_D_roots<- bootstrapped_summary_df

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

##Comparison of slopes and normalization constants between stems and roots----

#Add "organ" IDs for different datasets

bootstrapped_data_sample_d_D_stem_clean <- bootstrapped_data_sample_d_D_stem_clean %>%
  mutate(Organ = "Stems")

bootstrapped_data_sample_d_D_roots_clean <- bootstrapped_data_sample_d_D_roots_clean %>%
  mutate(Organ = "Roots")

#Merge datasets into one 

d_D_stems_roots_bootstrapped<- rbind(bootstrapped_data_sample_d_D_stem_clean, bootstrapped_data_sample_d_D_roots_clean)

#Slope test 

d_D_stems_roots_common_slope<- sma(DAVG~Diam_org*Organ, data = d_D_stems_roots_bootstrapped, log = "XY", method = "SMA")

d_D_stems_roots_common_slope

#One-tailed p-value

d_D_stems_roots_common_slope_onetail_p<- (0.000000000000000222/2)

d_D_stems_roots_common_slope_onetail_p

#norm_constant test 

d_D_stems_roots_common_norm_constant<- sma(DAVG~Diam_org+Organ, data = d_D_stems_roots_bootstrapped, log = "XY", method = "SMA")

d_D_stems_roots_common_norm_constant

#One-tailed p-value

d_D_stems_roots_common_norm_constant_onetail_p<- (0.000000000000000222/2)

d_D_stems_roots_common_norm_constant_onetail_p

#Hydraulic model----

#First, extract slope and y_intercept for scaling relationship between d and L

d_L_mod_boot_summary<- d_L_mod_boot$groupsummary

slope<- d_L_mod_boot_summary$Slope
y_int<- d_L_mod_boot_summary$Int

##Species & Individual characteristics----

###Callitropsis nootkatensis----

c_noot_mod_rup<- 44.1 #MPa, green, Ross 2021

c_noot_ind<- c("YeCe1",	"YeCe2","YeCe3","YeCe4","YeCe5")
c_noot_tree_h<- c(12.235,	12.712,	13.177,	10.936,	13.148)
c_noot_ID<- data.frame(c_noot_ind, c_noot_tree_h)

###Picea sitchensis----

p_sit_mod_rup<- 39.3 #MPa, green, Ross 2021

p_sit_ind<- c("SiSp1","SiSp2","SiSp3","SiSp4","SiSp5")
p_sit_tree_h<- c(15.41,	12.679,	13.413,	12.805,	14.362)
p_sit_ID<- data.frame(p_sit_ind, p_sit_tree_h)

###Tsuga heterophylla----

t_het_mod_rup<- 45.5 #MPa, green, Ross 2021

t_het_ind<- c("WeHe1","WeHe2","WeHe3","WeHe4","WeHe5","WeHe6")
t_het_tree_h<- c(3.866,	17.344,	9.733,	22.471,	12.475,	2.459)
t_het_ID<- data.frame(t_het_ind, t_het_tree_h)

###Tsuga mertensiana----

t_mer_mod_rup<- 43.4 #MPa, green, Ross 2021

t_mer_ind<- c("MoHe1","MoHe2","MoHe3","MoHe4","MoHe5")
t_mer_tree_h<- c(5.356,	6.295,	9.493,	10.996,	10.226)
t_mer_ID<- data.frame(t_mer_ind, t_mer_tree_h)

###Thuja plicata----

t_pli_mod_rup<- 35.9 #MPa, green, Ross 2021

t_pli_ind<- c("WeRe1","WeRe2","WeRe3","WeRe4","WeRe5", "WeRe6")
t_pli_tree_h<- c(8.78,	4.886,	13.072,	23.215,	27.437,	2.626)
t_pli_ID<- data.frame(t_pli_ind, t_pli_tree_h)

###Tree characteristics----

stem_scaling_eqn = function(stem){  
  return(10^(slope*log10(stem)+y_int)) #general scaling equation
}

gen_mod_rup<- sum(c_noot_mod_rup, p_sit_mod_rup, t_mer_mod_rup, t_het_mod_rup, t_pli_mod_rup)/5 #average for all species 

all_tree_h<- c(c_noot_tree_h, p_sit_tree_h, t_het_tree_h, t_mer_tree_h, t_pli_tree_h)
all_ind<- c(c_noot_ind, p_sit_ind, t_het_ind, t_mer_ind, t_pli_ind)
all_ID<- data.frame(all_ind, all_tree_h)

###Predict t_b_2_crit----

#Loop start
output_list <- list()

for (tree_h in all_tree_h) {
  
  ####Predicting average conduit (d) diameter along full path length----
  
  L_seq<- seq(0.0001,tree_h, by=0.0001)
  d_seq<- stem_scaling_eqn(L_seq)
  
  
  ####Assembling conduit data into one dataframe----
  
  df<- data.frame(L=L_seq, d=d_seq)
  
  ####redicting water potential gradient----
  
  hag_pos_eqn = function(d, l) {
    eta = 1e-9 #MPa.s
    return(128*eta*l/(pi*d^4))
  }
  
  psi_leaf<- -4.7 #MPa, minimum field measurement 
  psi_soil<- -1.5 #MPa, wilting point ##change after field soil WP values come in!!##
  
  delta_psi_full = psi_soil-psi_leaf
  
  inc<- 0.0001
  
  df <- df %>%
    mutate(
      r = hag_pos_eqn(d, inc),
      R = sum(r),
      r_cum = cumsum(r),
      Flow = delta_psi_full / R,
      delta_psi = Flow * r_cum,
      psi = psi_leaf + delta_psi
    )
  
  ####Adding water potential at the tip of the leaf----
  
  general_leaf_tip<- as.data.frame(matrix(nrow = 1, ncol = 0))
  general_leaf_tip = general_leaf_tip %>% mutate(L=0.0001,
                                                 d=NA,
                                                 r=NA,
                                                 R=NA,
                                                 r_cum=NA,
                                                 Flow=NA,
                                                 delta_psi=NA,
                                                 psi=psi_leaf)
  
  df<- rbind(general_leaf_tip,df)
  ####Implosion model----
  
  t_b_2_crit_eqn = function(p_def) {
    W=  gen_mod_rup#MPa, species specific modulus of rupture from wood database
    beta = 0.25 #coefficient which depends upon the ratio of conduit width to length - should be variable ideally 
    return((p_def*beta)/W)
  }
  
  df = df %>% mutate(t_b_2_crit = t_b_2_crit_eqn(abs(psi)))
  
  output_list[[as.character(tree_h)]] <- df
  
}

#Loop end

##Clean up model output dataframe----

general_df <- bind_rows(output_list, .id = "tree_h")

all_ID <- all_ID %>%
  rename(tree_h=all_tree_h)

general_df<- merge(general_df, all_ID, by = "tree_h", all.x = TRUE)

general_df <- general_df %>%
  rename(ind=all_ind)

#Make a column where ID is combined with L with an underscore separating the values

general_df$ind_L<- paste(general_df$ind, general_df$L, sep = "_")

#Make the same column for datasheet containing empirical data

df_full_aboveground$ind_L<- paste(df_full_aboveground$Individual, df_full_aboveground$L_trunc, sep = "_")

#Now, merge the model outputs with empirical data based on the "ind_L" column

general_df_full <- merge(df_full_aboveground, general_df, by = "ind_L", all.x = TRUE)

#Convert tree_h from character to numeric

general_df$tree_h<- as.numeric(general_df$tree_h)
general_df_full$tree_h<- as.numeric(general_df_full$tree_h)

#Add relative path to data

general_df_full<- general_df_full %>% mutate(rel_path = 100*(L_trunc/tree_h))

###Add xylem types to dataset----

#First, categorize the conduits by type of xylem (primary vs secondary)

general_df_full_xylem <- general_df_full %>%
  mutate(Xylem_type = case_when(
    Organ %in% c("Leaf", "Twig") ~ "Primary xylem",
    Organ %in% c("Branch", "Trunk") ~ "Secondary xylem",
  ))

#Then, split the dataset by xylem types

general_df_full_primary_xylem <- subset(general_df_full_xylem, Xylem_type == "Primary xylem")
general_df_full_secondary_xylem <- subset(general_df_full_xylem, Xylem_type == "Secondary xylem")

#In primary xylem dataset, add column called Xylem add populate with "PX" (primary xylem)

general_df_full_primary_xylem<- general_df_full_primary_xylem %>%
  mutate(Xylem = "PX")

#In secondary xylem, create a column called Xylem and populate depending on whether its earlywood or latewood via Mork's Index

general_df_full_secondary_xylem<- general_df_full_secondary_xylem %>%
  mutate(Xylem = ifelse(general_df_full_secondary_xylem$RTSR > 1, "LW", "EW"))

#Finally, combine the two datasets

general_df_full_all_xylem <- rbind(general_df_full_primary_xylem, general_df_full_secondary_xylem)

##Estimate % of imploded conduits in primary xylem, earlywood, latewood, etc.

implosion_summary <- general_df_full_all_xylem %>%
  group_by(Xylem) %>%
  dplyr::summarize(Count_imploded = sum(t_b_2_crit >= t_b_2),
            Count_safe = sum(t_b_2_crit < t_b_2),
            Percent_imploded = (Count_imploded/(Count_imploded+Count_safe))*100)

##Estimates of collapsed conduits using species-specific modulus of rupture----

#Create df with modulus of rupture values

mod_rup_df <- data.frame(
  Species = c("C.nootkatensis", "P.sitchensis", "T.plicata", "T.heterophylla", "T.mertensiana"),
  modulus_of_rupture = c(c_noot_mod_rup, p_sit_mod_rup, t_pli_mod_rup, t_het_mod_rup, t_mer_mod_rup)
)

#Merge above dataframe with dataframe containing thickness to span

general_df_full_all_xylem<- merge(general_df_full_all_xylem, mod_rup_df, by = "Species", all.x = TRUE)

#Calculate t_b_2_crit using species specific modulus of rupture

general_df_full_all_xylem<- general_df_full_all_xylem %>%
  mutate(t_b_2_crit_species = (abs(psi)*0.25)/modulus_of_rupture)

#Summarize imploded conduits

implosion_summary_species <- general_df_full_all_xylem %>%
  group_by(Xylem) %>%
  dplyr::summarize(Count_imploded = sum(t_b_2_crit_species >= t_b_2),
            Count_safe = sum(t_b_2_crit_species < t_b_2),
            Percent_imploded = (Count_imploded/(Count_imploded+Count_safe))*100)

##Model comparison for relationship between t/b^2 and relative path length----

#Linear model

t_b_2_rel_path_lin_model<- lm(data=general_df_full, t_b_2~rel_path)

summary.lm(t_b_2_rel_path_lin_model)

#Logarithmic model

t_b_2_rel_path_log_model<- lm(data=general_df_full, t_b_2~log10(rel_path))

summary.lm(t_b_2_rel_path_log_model)

#Power law

t_b_2_rel_path_power_law<- lm(data=general_df_full, log10(t_b_2)~log10(rel_path))

summary.lm(t_b_2_rel_path_power_law)

#Exponential

t_b_2_rel_path_exp_model<- lm(data=general_df_full, log(t_b_2)~rel_path)

summary.lm(t_b_2_rel_path_exp_model)

#Compare models via AIC

AIC(t_b_2_rel_path_lin_model)
AIC(t_b_2_rel_path_log_model)
AIC(t_b_2_rel_path_power_law)
AIC(t_b_2_rel_path_exp_model)

##Safety margins and factors (Table S4)----

#Create new dataframe with t/b^2 and t/b^2 crit

safety<- general_df_full %>% 
  dplyr::select(ID, Species, Organ, L_trunc, H_wleaves,RTSR, psi, t_b_2, t_b_2_crit)

#Calculate safety margins and factors for individual tracheids

safety<- safety %>%
  mutate(P_crit = (gen_mod_rup/0.25)*t_b_2,
         safety_factor = (P_crit)/abs(psi),
         rel_path = 100*(L_trunc/H_wleaves))

##Add xylem type

safety <- safety %>%
  mutate(Xylem = case_when(
    Organ %in% c("Leaf", "Twig") ~ "PX",
    RTSR > 1 ~ "LW",
    TRUE ~ "EW"
  ))

##Summarize safety margins and factors

safety_summary_all<- safety %>%
  dplyr::summarize(
    safety_factor_median = median(safety_factor),
    safety_factor_5th_perc = quantile(safety_factor, probs = 0.05),
    safety_factor_95th_perc = quantile(safety_factor, probs = 0.95),
    safety_factor_IQR = IQR(safety_factor),
    safety_factor_25th_perc = quantile(safety_factor, probs = 0.25),
    safety_factor_75th_perc = quantile(safety_factor, probs = 0.75),
  )

safety_summary_by_xylem <- safety %>%
  group_by(Xylem) %>%
  dplyr::summarize(
    t_b_2_5th_perc = quantile(t_b_2, probs = 0.05),
    t_b_2_95th_perc = quantile(t_b_2, probs = 0.95),
    safety_factor_median = median(safety_factor),
    safety_factor_5th_perc = quantile(safety_factor, probs = 0.05),
    safety_factor_95th_perc = quantile(safety_factor, probs = 0.95),
    safety_factor_IQR = IQR(safety_factor),
    safety_factor_25th_perc = quantile(safety_factor, probs = 0.25),
    safety_factor_75th_perc = quantile(safety_factor, probs = 0.75)
  )

safety_summary_by_xylem_species <- safety %>%
  group_by(Xylem, Species) %>%
  dplyr::summarize(
    safety_factor_median = median(safety_factor),
    safety_factor_5th_perc = quantile(safety_factor, probs = 0.05),
    safety_factor_95th_perc = quantile(safety_factor, probs = 0.95),
    safety_factor_IQR = IQR(safety_factor),
    safety_factor_25th_perc = quantile(safety_factor, probs = 0.25),
    safety_factor_75th_perc = quantile(safety_factor, probs = 0.75)
  )

###Absolute mean, minima and maxima----

#Mean

mean(safety$safety_margin)
mean(safety$safety_factor)

#Min

min(safety$safety_margin)
min(safety$safety_factor)

#Max

max(safety$safety_margin)
max(safety$safety_factor)


#Supplementary analyses (for PCE publication)----
##OLS fits for all scaling relationships (Table S2)----

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

##Model comparison for key scaling relationships (Table S3)----

###d ~ L----

bootstrapped_data_sample_d_L$log_L<- log10(bootstrapped_data_sample_d_L$L) 
bootstrapped_data_sample_d_L$log_DAVG<- log10(bootstrapped_data_sample_d_L$DAVG)

#Linear fit to log-transformed variables (power law)

d_L_power_law<- lm(data=bootstrapped_data_sample_d_L, log_DAVG~log_L)

summary.lm(d_L_power_law)

#Logarithmic fit to log-transformed variables

d_L_log_model<- lm(data=bootstrapped_data_sample_d_L, DAVG~log_L)

summary.lm(d_L_log_model)

#Exponential fit to log-transformed variables

d_L_exp_model<- lm(data=bootstrapped_data_sample_d_L, log_DAVG~L)

summary.lm(d_L_exp_model)

#Piecewise fit

d_L_seg_model <- segmented(d_L_power_law)

summary(d_L_seg_model)

#Quadratic model to detect curvature
d_L_quad_model <- lm(log_DAVG ~ log_L + I(log_L^2), data = bootstrapped_data_sample_d_L)
summary(d_L_quad_model)

#AIC comparisons

AIC(d_L_power_law)
AIC(d_L_log_model)
AIC(d_L_exp_model)
AIC(d_L_seg_model)
AIC(d_L_quadratic)

###d ~ D_stem----

bootstrapped_data_sample_d_D_stem$log_DAVG<- log10(bootstrapped_data_sample_d_D_stem$DAVG)
bootstrapped_data_sample_d_D_stem$log_D<- log10(bootstrapped_data_sample_d_D_stem$Diam_org)

#Linear fit to log-transformed variables (power law)

d_D_stem_power_law<- lm(data=bootstrapped_data_sample_d_D_stem, log_DAVG~log_D)

summary.lm(d_D_stem_power_law)

#Logarithmic fit to log-transformed variables

d_D_stem_log_model<- lm(data=bootstrapped_data_sample_d_D_stem, DAVG~log_D)

summary.lm(d_D_stem_log_model)

#Exponential fit to log-transformed variables

d_D_stem_exp_model<- lm(data=bootstrapped_data_sample_d_D_stem, log_DAVG~Diam_org)

summary.lm(d_D_stem_exp_model)

#Piecewise fit

d_D_stem_seg_model <- segmented(d_D_stem_power_law)

summary(d_D_stem_seg_model)

#Quadratic model to detect curvature
d_D_stem_quad_model <- lm(log_DAVG ~ log_D + I(log_D^2), data = bootstrapped_data_sample_d_D_stem)
summary(d_D_stem_quad_model)

#AIC comparisons

AIC(d_D_stem_power_law)
AIC(d_D_stem_log_model)
AIC(d_D_stem_exp_model)
AIC(d_D_stem_seg_model)
AIC(d_D_stem_quadratic)

###d ~ D_roots----

bootstrapped_data_sample_d_D_roots$log_DAVG<- log10(bootstrapped_data_sample_d_D_roots$DAVG)
bootstrapped_data_sample_d_D_roots$log_D<- log10(bootstrapped_data_sample_d_D_roots$Diam_org)

#Linear fit to log-transformed variables (power law)

d_D_roots_power_law<- lm(data=bootstrapped_data_sample_d_D_roots, log_DAVG~log_D)

summary.lm(d_D_roots_power_law)

#Logarithmic fit to log-transformed variables

d_D_roots_log_model<- lm(data=bootstrapped_data_sample_d_D_roots, DAVG~log_D)

summary.lm(d_D_roots_log_model)

#Exponential fit to log-transformed variables

d_D_roots_exp_model<- lm(data=bootstrapped_data_sample_d_D_roots, log_DAVG~Diam_org)

summary.lm(d_D_roots_exp_model)

#Piecewise fit

d_D_roots_seg_model <- segmented(d_D_roots_power_law)

summary(d_D_roots_seg_model)

#Quadratic model to detect curvature
d_D_roots_quad_model <- lm(log_DAVG ~ log_D + I(log_D^2), data = bootstrapped_data_sample_d_D_roots)
summary(d_D_roots_quad_model)

#AIC comparisons

AIC(d_D_roots_power_law)
AIC(d_D_roots_log_model)
AIC(d_D_roots_exp_model)
AIC(d_D_roots_seg_model)
AIC(d_D_roots_quadratic)

##Fitting SMA models to individual trees and calculating mean slope (Fig. S12)----
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

#Manuscript figures----

##Figure 1----

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

ggsave("Figure_1.SVG", plot=Fig1_full, width = 25, height = 10, dpi = 600, units = "cm")

##Figure 2----

#Set density range

density_range<- c(0,180000)

#plot d~L using bootstrapped data

bootstrapped_sample_plot_d_L <- ggplot(bootstrapped_data_sample_d_L, aes(x=L, y=DAVG)) +
  geom_pointdensity(alpha=1, adjust =0.5) + 
  scale_color_viridis_c(option = "viridis", name= "Point density", alpha = 0.7, limits = density_range)+
  annotate("rect", xmin = 0.004, xmax = Inf, ymin = 1, ymax = Inf,
           fill = "white", alpha = 0.3) +
  stat_ma_line(method = "SMA", color="black") +
  scale_y_log10(breaks=c(3,10,30), limits=c(1,80)) +
  scale_x_log10(breaks=c(0.01,0.1, 1, 10), labels =c(0.01,0.1, 1, 10), limits = c(0.004, 30)) +
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

#Plot d~D using bootstrapped data

bootstrapped_sample_plot_d_D_stem <- ggplot(data = bootstrapped_data_sample_d_D_stem, aes(x=Diam_org, y=DAVG)) +
  geom_pointdensity(alpha=1, adjust =0.5) + 
  scale_color_viridis_c(option = "viridis", name= "Point density", alpha = 0.7, limits = density_range)+
  annotate("rect", xmin = 0.001, xmax = Inf, ymin = 1, ymax = Inf,
           fill = "white", alpha = 0.3) +
  stat_ma_line(data = bootstrapped_data_subsample_d_D_stem, aes(x=Diam_org, y=DAVG), method = "SMA", color= "black") +
  scale_y_log10(breaks=c(3,10,30), limits=c(1,80)) +
  scale_x_log10(breaks=c(0.001,0.01,0.1, 1), labels=c(0.001,0.01,0.1, 1), limits=c(0.001, 1)) +
  xlab("Stem diameter (m)") +
  ylab("Xylem conduit diameter (µm)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=14),
    legend.position = "right",
    aspect.ratio = 1
  )

#Plot d~Droot using bootstrapped data

bootstrapped_sample_plot_d_D_roots <- ggplot(data=bootstrapped_data_sample_d_D_roots, aes(x=Diam_org, y=DAVG)) +
  geom_pointdensity(alpha=1, adjust =0.5) + 
  scale_color_viridis_c(option = "viridis", name= "Point density", alpha = 0.7, limits = density_range)+
  annotate("rect", xmin = 0.00015, xmax = Inf, ymin = 1, ymax = Inf,
           fill = "white", alpha = 0.3) +
  stat_ma_line(linewidth=1, method = "SMA", color= "black") +
  scale_y_log10(breaks=c(3,10,30), limits=c(1,80)) +
  scale_x_log10(breaks=c(0.0003, 0.003,0.03), labels=c(0.0003, 0.003,0.03), limits=c(0.00015,0.03)) +
  xlab("Root diameter (m)") +
  ylab("Xylem conduit diameter (µm)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size=14),
    legend.position = "right",
    aspect.ratio = 1
  )

Fig_2 <- bootstrapped_sample_plot_d_L + 
  labs(tag = 'A') +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.1, 1)) +  
  bootstrapped_sample_plot_d_D_stem + 
  labs(tag = 'B') +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.03, 1)) + 
  bootstrapped_sample_plot_d_D_roots + 
  labs(tag = 'C') +
  theme(plot.tag = element_text(size = 14, face = "bold"),
        plot.tag.position = c(0.02, 1)) + 
  plot_layout(guides = 'collect') & theme(legend.position = 'right',
                                          legend.text = element_text(size = 12),
                                          legend.title = element_text(size=14))


ggsave("Figure_2.JPG", plot=Fig_2, width = 28, height = 10, dpi = 600, units = "cm")

##Figure 3----

implosion_plot<- ggplot(data=general_df_full, aes(x=rel_path, y=t_b_2)) +
  geom_pointdensity(alpha=1, adjust = 0.5) + 
  scale_color_viridis_c(option = "viridis", name= "Point density", alpha = 0.7)+
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.002, ymax = Inf,
           fill = "white", alpha = 0.3) +
  geom_smooth(aes(x=rel_path, y=t_b_2), method="lm", color="black", linewidth=1, linetype="solid", se=TRUE) +
  geom_smooth(aes(x=rel_path, y=t_b_2_crit), method="loess", color="black", linewidth=1, linetype="dotdash", se=FALSE) +
  scale_y_log10()+
  xlab("Relative position from leaf tip to base (%)") +
  ylab("Thickness-to-span ratio (dimensionless)")+
  theme_classic() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size=14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size=14),
    aspect.ratio = 1,
    legend.margin = margin(t = 0, unit = "lines"),
    legend.position = "right"
  ) +
  annotate("text", x = 50, y = 0.005, label = "Critical collapse limit")

ggsave("Figure_3.JPG", plot=implosion_plot, width = 14, height = 14, dpi = 600, units = "cm")

##Figure S3----

data_plot_full_d_L <- ggplot(df_full_aboveground, aes(x=L, y=DAVG)) +
  geom_pointdensity(alpha=1, adjust =0.5) + 
  scale_color_viridis_c(option = "viridis", name= "Point density", alpha = 0.7)+
  annotate("rect", xmin = 0.002, xmax = Inf, ymin = 1, ymax = Inf,
           fill = "white", alpha = 0.3) +
  stat_ma_line(linewidth=1, method = "SMA", color= "black") +
  scale_y_log10(breaks=c(3,10,30)) +
  scale_x_log10(limits=c(0.001, 30), breaks=c(0.01, 0.1, 1, 10), labels= c(0.01, 0.1, 1, 10)) +
  xlab("Distance from leaf tip (m)") +
  ylab("Xylem diameter (µm)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    aspect.ratio = 1
  )

ggsave("Figure_S3.jpg", plot = data_plot_full_d_L, width = 12, height = 10, units = "cm", dpi = 600)

##Figure S4----

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

##Figure S5----

#Change num_iterations from numeric to factor
sampled_summary_df_d_L$num_iterations<- as.factor(sampled_summary_df_d_L$num_iterations)

sampled_slope_plot<- ggplot(sampled_summary_df_d_L, aes(x = num_iterations, y = mean_slope)) +
  geom_line(color = "blue", linewidth = 1, group=1) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), color = "blue", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Slope Estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )


sampled_lowCI_plot<- ggplot(sampled_summary_df_d_L, aes(x = num_iterations, y = mean_slope_lowCI)) +
  geom_line(color = "red", linewidth = 1, group=1) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_lowCI - se_slope_lowCI, ymax = mean_slope_lowCI + se_slope_lowCI), color = "red", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Lower CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

sampled_highCI_plot<- ggplot(sampled_summary_df_d_L, aes(x = num_iterations, y = mean_slope_highCI)) +
  geom_line(color = "green4", linewidth = 1, group=1) +
  geom_point(color = "green4", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_highCI - se_slope_highCI, ymax = mean_slope_highCI + se_slope_highCI), color = "green4", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Upper CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

sampled_plot_slope_CI_full_d_L<- ggarrange(sampled_slope_plot, sampled_lowCI_plot, sampled_highCI_plot, nrow=1)

sampled_plot_slope_CI_full_d_L

ggsave("Figure_S5.SVG", plot = sampled_plot_slope_CI_full_d_L, width = 27, height = 10, units = "cm", dpi = 600)

##Figure S6----

bootstrapped_summary_df_d_L$num_iterations<- as.factor(bootstrapped_summary_df_d_L$num_iterations)

bootstrapped_slope_plot<- ggplot(bootstrapped_summary_df_d_L, aes(x = num_iterations, y = mean_slope)) +
  geom_line(color = "blue", linewidth = 1, group=1) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), color = "blue", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Slope Estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_lowCI_plot<- ggplot(bootstrapped_summary_df_d_L, aes(x = num_iterations, y = mean_slope_lowCI)) +
  geom_line(color = "red", size = 1, group=1) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_lowCI - se_slope_lowCI, ymax = mean_slope_lowCI + se_slope_lowCI), color = "red", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Lower CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_highCI_plot<- ggplot(bootstrapped_summary_df_d_L, aes(x = num_iterations, y = mean_slope_highCI)) +
  geom_line(color = "green4", size = 1, group=1) +
  geom_point(color = "green4", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_highCI - se_slope_highCI, ymax = mean_slope_highCI + se_slope_highCI), color = "green4", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Upper CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_plot_slope_CI_full_d_L<- ggarrange(bootstrapped_slope_plot, bootstrapped_lowCI_plot, bootstrapped_highCI_plot, nrow=1)

bootstrapped_plot_slope_CI_full_d_L

ggsave("Figure_S6.SVG", plot = bootstrapped_plot_slope_CI_full_d_L, width = 27, height = 10, units = "cm", dpi = 600)

##Figure S7----

#Change num_iterations from numeric to factor
sampled_summary_df_d_D_stem$num_iterations<- as.factor(sampled_summary_df_d_D_stem$num_iterations)

sampled_slope_plot<- ggplot(sampled_summary_df_d_D_stem, aes(x = num_iterations, y = mean_slope)) +
  geom_line(color = "blue", size = 1, group=1) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), color = "blue", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Slope Estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )


sampled_lowCI_plot<- ggplot(sampled_summary_df_d_D_stem, aes(x = num_iterations, y = mean_slope_lowCI)) +
  geom_line(color = "red", size = 1, group=1) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_lowCI - se_slope_lowCI, ymax = mean_slope_lowCI + se_slope_lowCI), color = "red", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Lower CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

sampled_highCI_plot<- ggplot(sampled_summary_df_d_D_stem, aes(x = num_iterations, y = mean_slope_highCI)) +
  geom_line(color = "green4", size = 1, group=1) +
  geom_point(color = "green4", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_highCI - se_slope_highCI, ymax = mean_slope_highCI + se_slope_highCI), color = "green4", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Upper CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

sampled_plot_slope_CI_full_d_D<- ggarrange(sampled_slope_plot, sampled_lowCI_plot, sampled_highCI_plot, nrow=1)

sampled_plot_slope_CI_full_d_D

ggsave("Figure_S7.SVG", plot = sampled_plot_slope_CI_full_d_D, width = 27, height = 10, units = "cm", dpi = 600)

##Figure S8----

bootstrapped_summary_df_d_D_stem$num_iterations<- as.factor(bootstrapped_summary_df_d_D_stem$num_iterations)

bootstrapped_slope_plot<- ggplot(bootstrapped_summary_df_d_D_stem, aes(x = num_iterations, y = mean_slope)) +
  geom_line(color = "blue", linewidth = 1, group=1) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), color = "blue", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Slope Estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_lowCI_plot<- ggplot(bootstrapped_summary_df_d_D_stem, aes(x = num_iterations, y = mean_slope_lowCI)) +
  geom_line(color = "red", size = 1, group=1) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_lowCI - se_slope_lowCI, ymax = mean_slope_lowCI + se_slope_lowCI), color = "red", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Lower CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_highCI_plot<- ggplot(bootstrapped_summary_df_d_D_stem, aes(x = num_iterations, y = mean_slope_highCI)) +
  geom_line(color = "green4", size = 1, group=1) +
  geom_point(color = "green4", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_highCI - se_slope_highCI, ymax = mean_slope_highCI + se_slope_highCI), color = "green4", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Upper CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_plot_slope_CI_full_d_D<- ggarrange(bootstrapped_slope_plot, bootstrapped_lowCI_plot, bootstrapped_highCI_plot, nrow=1)

bootstrapped_plot_slope_CI_full_d_D

ggsave("Figure_S8.SVG", plot = bootstrapped_plot_slope_CI_full_d_D, width = 27, height = 10, units = "cm", dpi = 600)

##Figure S9----

#Change num_iterations from numeric to factor
sampled_summary_df_d_D_roots$num_iterations<- as.factor(sampled_summary_df_d_D_roots$num_iterations)

sampled_slope_plot<- ggplot(sampled_summary_df_d_D_roots, aes(x = num_iterations, y = mean_slope)) +
  geom_line(color = "blue", size = 1, group=1) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), color = "blue", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Slope Estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )


sampled_lowCI_plot<- ggplot(sampled_summary_df_d_D_roots, aes(x = num_iterations, y = mean_slope_lowCI)) +
  geom_line(color = "red", size = 1, group=1) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_lowCI - se_slope_lowCI, ymax = mean_slope_lowCI + se_slope_lowCI), color = "red", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Lower CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

sampled_highCI_plot<- ggplot(sampled_summary_df_d_D_roots, aes(x = num_iterations, y = mean_slope_highCI)) +
  geom_line(color = "green4", size = 1, group=1) +
  geom_point(color = "green4", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_highCI - se_slope_highCI, ymax = mean_slope_highCI + se_slope_highCI), color = "green4", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Upper CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

sampled_plot_slope_CI_full_d_D_roots<- ggarrange(sampled_slope_plot, sampled_lowCI_plot, sampled_highCI_plot, nrow=1)

ggsave("Figure_S9.SVG", plot = sampled_plot_slope_CI_full_d_D_roots, width = 27, height = 10, units = "cm", dpi = 600)

##Figure S10----

bootstrapped_summary_df_d_D_roots$num_iterations<- as.factor(bootstrapped_summary_df_d_D_roots$num_iterations)

bootstrapped_slope_plot<- ggplot(bootstrapped_summary_df_d_D_roots, aes(x = num_iterations, y = mean_slope)) +
  geom_line(color = "blue", linewidth = 1, group=1) +
  geom_point(color = "blue", size = 3) +
  geom_errorbar(aes(ymin = mean_slope - se_slope, ymax = mean_slope + se_slope), color = "blue", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Slope Estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_lowCI_plot<- ggplot(bootstrapped_summary_df_d_D_roots, aes(x = num_iterations, y = mean_slope_lowCI)) +
  geom_line(color = "red", size = 1, group=1) +
  geom_point(color = "red", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_lowCI - se_slope_lowCI, ymax = mean_slope_lowCI + se_slope_lowCI), color = "red", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Lower CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_highCI_plot<- ggplot(bootstrapped_summary_df_d_D_roots, aes(x = num_iterations, y = mean_slope_highCI)) +
  geom_line(color = "green4", size = 1, group=1) +
  geom_point(color = "green4", size = 3) +
  geom_errorbar(aes(ymin = mean_slope_highCI - se_slope_highCI, ymax = mean_slope_highCI + se_slope_highCI), color = "green4", width = 0.2) +
  labs(x = "Number of Iterations",
       y = "Upper CI estimate") +
  theme_classic()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

bootstrapped_plot_slope_CI_full_d_D_roots<- ggarrange(bootstrapped_slope_plot, bootstrapped_lowCI_plot, bootstrapped_highCI_plot, nrow=1)

ggsave("Figure_S10.SVG", plot = bootstrapped_plot_slope_CI_full_d_D_roots, width = 27, height = 10, units = "cm", dpi = 600)

##Figure S11 ----

#hydraulic diameter plot

data_plot_means_d_L<- ggplot(df_means_aboveground, aes(x=L, y=DhAVG))+
  geom_point(size=1, color= "#0d0887", alpha = 0.5)+
  stat_ma_line(method = "SMA", color="black")+
  scale_y_log10(breaks=c(3,10,30), limits=c(1,50))+
  scale_x_log10(breaks=c(0.01,0.1, 1, 10), labels =c(0.01,0.1, 1, 10))+
  xlab("Distance from leaf tip (m)")+
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

sampled_data_plot_d_L<- ggplot(sampled_data_sample_d_L, aes(x=L, y=DAVG))+
  geom_point(size=1, color="#0d0887", alpha = 0.5)+
  stat_ma_line(method = "SMA", color="black")+
  scale_y_log10(breaks=c(3,10,30), limits=c(1,50))+
  scale_x_log10(breaks=c(0.01,0.1, 1, 10), labels =c(0.01,0.1, 1, 10))+
  xlab("Distance from leaf tip (m)")+
  ylab("Xylem diameter (µm)")+
  theme_classic()+
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    aspect.ratio = 1)

#Combined plot

Fig_S11<- ggarrange(data_plot_means_d_L, sampled_data_plot_d_L, hjust=-2, vjust = 1.3, labels = c("A", "B"), nrow=1)

Fig_S11

#Create an empty plot with the common x-axis label
empty_plot <- ggplot()  + 
  labs(x = "Distance from leaf tip (m)")+
  theme(axis.title.x = element_text(size= 15))

empty_plot

Fig_S11_fin <- ggarrange(Fig_S11, empty_plot, ncol = 1, nrow = 2, heights = c(1, 0.06))

Fig_S11_fin

ggsave("FigS11_fin.JPG", plot=Fig_S11_fin, width = 20, height = 10, dpi = 600, units = "cm")

##Figure S12----

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

##Figure S13----

#hydraulic diameter plot

data_plot_means_d_D_stem<- ggplot(df_means_stems, aes(x=Diam_org, y=DhAVG))+
  geom_point(size=1, color="#0d0887", alpha = 0.5)+
  stat_ma_line(method = "SMA", color="black")+
  scale_y_log10(breaks=c(3,10,30), limits=c(1,50))+
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1), labels=c(0.001, 0.01, 0.1, 1))+
  xlab("Stem diameter (m)")+
  ylab("Xylem diameter (µm)")+
  theme_classic()+
  theme(
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    axis.title.y = element_text(size = 14),
    axis.title.x=element_blank(),
    aspect.ratio = 1
  )

data_plot_means_d_D_stem

#Subsampled plot

sampled_data_plot_d_D_stem<- ggplot(sampled_data_sample_d_D_stem, aes(x=Diam_org, y=DAVG))+
  geom_point(size=1, color="#0d0887", alpha = 0.5)+
  stat_ma_line(method = "SMA", color="black")+
  scale_y_log10(breaks=c(3,10,30), limits=c(1,50))+
  scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1), labels=c(0.001, 0.01, 0.1, 1))+
  xlab("Stem diameter (m)")+
  ylab("Xylem diameter (µm)")+
  theme_classic()+
  theme(
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    axis.text.x = element_text(size=12),
    axis.text.y = element_text(size=12),
    aspect.ratio = 1)


sampled_data_plot_d_D_stem

#Combined plot

Fig_S13<- ggarrange(data_plot_means_d_D_stem, sampled_data_plot_d_D_stem, hjust=-2, vjust = 1.3, labels = c("A", "B"), nrow=1)

Fig_S13

#Create an empty plot with the common x-axis label
empty_plot <- ggplot()  + 
  labs(x = "Stem diameter (m)")+
  theme(axis.title.x = element_text(size= 15))

empty_plot

Fig_S13_fin <- ggarrange(Fig_S13, empty_plot, ncol = 1, nrow = 2, heights = c(1, 0.06))

Fig_S13_fin

ggsave("FigS13_fin.JPG", plot=Fig_S13_fin, width = 22, height = 11, dpi = 600, units = "cm")

##Figure S14 ----

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

data_plot_means_d_D_roots

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


sampled_data_plot_d_D_roots

#Combined plot

Fig_S14<- ggarrange(data_plot_means_d_D_roots, sampled_data_plot_d_D_roots, hjust=-2, vjust = 1.3, labels = c("A", "B"), nrow=1)

Fig_S14

#Create an empty plot with the common x-axis label
empty_plot <- ggplot()  + 
  labs(x = "Root diameter (m)")+
  theme(axis.title.x = element_text(size= 15))

empty_plot

Fig_S14_fin <- ggarrange(Fig_S14, empty_plot, ncol = 1, nrow = 2, heights = c(1, 0.06))

Fig_S14_fin

ggsave("FigS14_fin.JPG", plot=Fig_S14_fin, width = 22, height = 11, dpi = 600, units = "cm")

##Figure S15----

##Line fits

summary(sma(data = df_means_aboveground, L~Diam_org, slope.test = 2/3, method = "SMA", log = "XY"))

##Plot

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

Fig_S15

ggsave("Fig_S15.JPG", plot=Fig_S15, width = 12, height = 10, dpi = 600, units = "cm")

##Figure S16----

FigS16<- ggplot(aes(x=DAVG, y=CWTALL), data=general_df_full)+
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

ggsave("FigS16.JPG", plot=FigS16, width = 12, height = 10, dpi = 600, units = "cm")