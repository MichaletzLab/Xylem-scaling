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

#Subset root data

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

#Figure----

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

ggsave("Figure_S10.SVG", plot = bootstrapped_plot_slope_CI_full_d_D_roots, width = 27, height = 10, units = "cm")