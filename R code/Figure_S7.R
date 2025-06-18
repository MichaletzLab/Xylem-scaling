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

#Subset stem data

df_full_stems <- df_full[which(df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_stems$Organ<- factor(df_full_stems$Organ, levels=(c("Twig", "Branch", "Trunk"))) 

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

#Figure----

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

ggsave("Figure_S7.SVG", plot = sampled_plot_slope_CI_full_d_D, width = 27, height = 10, units = "cm")