#Preparation----
##Install required packages----

#List of packages
pkg_list <- c("ggplot2", "dplyr", "readxl", "smatr", "cowplot", 
              "ggpubr", "ggpointdensity", "ggpmisc", "emmeans", 
              "multcomp", "ggridges", "patchwork", "segmented", 
              "viridis", "Hmisc", "ggrastr")

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

#Subset data

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'| df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk"))) 

df_full_stems <- df_full[which(df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_stems$Organ<- factor(df_full_stems$Organ, levels=(c("Twig", "Branch", "Trunk"))) 

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

#Fit model to bootstrapped dataset

sma(DAVG~L, data=bootstrapped_data_sample_d_L, method = "SMA", log = "XY")

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

#Figure----

#Set density range

density_range<- c(0,180000)

#plot d~L using bootstrapped data

bootstrapped_sample_plot_d_L <- ggplot(bootstrapped_data_sample_d_L, aes(x=L, y=DAVG)) +
  rasterise(geom_pointdensity(alpha=1, adjust =0.5), dpi = 300) + 
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
  rasterise(geom_pointdensity(alpha=1, adjust =0.5), dpi = 300) + 
  scale_color_viridis_c(option = "viridis", name= "Point density", alpha = 0.7, limits = density_range)+
  annotate("rect", xmin = 0.001, xmax = Inf, ymin = 1, ymax = Inf,
           fill = "white", alpha = 0.3) +
  stat_ma_line(data = bootstrapped_data_sample_d_D_stem, aes(x=Diam_org, y=DAVG), method = "SMA", color= "black") +
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
  rasterise(geom_pointdensity(alpha=1, adjust =0.5), dpi = 300) + 
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

ggsave("Figure_2.pdf", plot=Fig_2, width = 28, height = 10, units = "cm", device = cairo_pdf)