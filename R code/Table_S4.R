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

#Subset data

df_full_aboveground <- df_full[which(df_full$Organ=='Leaf'| df_full$Organ=='Twig'| df_full$Organ=='Branch'| df_full$Organ=='Trunk'), ]
df_full_aboveground$Organ<- factor(df_full_aboveground$Organ, levels=(c("Leaf", "Twig", "Branch", "Trunk"))) 

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

##Hydraulic model----

#First, extract slope and y_intercept for scaling relationship between d and L

d_L_mod_boot_summary<- d_L_mod_boot$groupsummary

slope<- d_L_mod_boot_summary$Slope
y_int<- d_L_mod_boot_summary$Int

###Species & Individual characteristics----

#Callitropsis nootkatensis

c_noot_mod_rup<- 44.1 #MPa, green, Ross 2021

c_noot_ind<- c("YeCe1",	"YeCe2","YeCe3","YeCe4","YeCe5")
c_noot_tree_h<- c(12.235,	12.712,	13.177,	10.936,	13.148)
c_noot_ID<- data.frame(c_noot_ind, c_noot_tree_h)

#Picea sitchensis

p_sit_mod_rup<- 39.3 #MPa, green, Ross 2021

p_sit_ind<- c("SiSp1","SiSp2","SiSp3","SiSp4","SiSp5")
p_sit_tree_h<- c(15.41,	12.679,	13.413,	12.805,	14.362)
p_sit_ID<- data.frame(p_sit_ind, p_sit_tree_h)

#Tsuga heterophylla

t_het_mod_rup<- 45.5 #MPa, green, Ross 2021

t_het_ind<- c("WeHe1","WeHe2","WeHe3","WeHe4","WeHe5","WeHe6")
t_het_tree_h<- c(3.866,	17.344,	9.733,	22.471,	12.475,	2.459)
t_het_ID<- data.frame(t_het_ind, t_het_tree_h)

#Tsuga mertensiana

t_mer_mod_rup<- 43.4 #MPa, green, Ross 2021

t_mer_ind<- c("MoHe1","MoHe2","MoHe3","MoHe4","MoHe5")
t_mer_tree_h<- c(5.356,	6.295,	9.493,	10.996,	10.226)
t_mer_ID<- data.frame(t_mer_ind, t_mer_tree_h)

#Thuja plicata

t_pli_mod_rup<- 35.9 #MPa, green, Ross 2021

t_pli_ind<- c("WeRe1","WeRe2","WeRe3","WeRe4","WeRe5", "WeRe6")
t_pli_tree_h<- c(8.78,	4.886,	13.072,	23.215,	27.437,	2.626)
t_pli_ID<- data.frame(t_pli_ind, t_pli_tree_h)

#Tree characteristics

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

##Safety margins and factors----

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

mean(safety$safety_factor)

#Min

min(safety$safety_factor)

#Max

max(safety$safety_factor)