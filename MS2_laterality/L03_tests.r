#This code runs permutation tests to test whether honey buzzards have a preference for left or right bank angles during different life cycle stage
#Elham Nourani PhD.
#July 22. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)
library(viridis)
library(lme4)

#-------------------------------------------------------------------------------------
# STEP1: annotate data with life cycle stage
#-------------------------------------------------------------------------------------
#open data with laterality index calculated

#this file contains laterality index calculated per day. in the wintering ground and some points pre-migration, there was only one IMU recorded per day (n_bursts_in_day). filter as necessary
laterality_8sec <- readRDS("laterality_index_per_day.rds")

#use Ellen's life cycle stage estimations
life_cycle <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/EHB_metadata_mig_dates_Fin22.csv") %>% 
  mutate_at(c("migration_start", "migration_end"), as.Date) %>% 
  filter(nest_country == "Finland")

laterality_8sec_LS <- laterality_8sec %>%
  mutate(unique_date = as.Date(unique_date)) %>% 
  full_join(life_cycle %>% select(ring_ID, migration_start, migration_end), by = c("individual_local_identifier" = "ring_ID")) %>% 
  rowwise() %>% 
  mutate(life_stage = case_when(
    between(unique_date, migration_start, migration_end) ~ "migration",
    unique_date < migration_start ~ "post-fledging",
    unique_date > migration_end ~ "wintering",
    TRUE ~ NA_character_
  )) %>% 
  ungroup() %>% 
  mutate(life_stage = case_when(
    is.na(life_stage) & unique_date < as.Date("2023-09-15") ~ "post-fledging",
    is.na(life_stage) & unique_date > as.Date("2023-10-20") ~ "migration",
    is.na(life_stage) & between(unique_date, as.Date("2023-09-15"), as.Date("2023-10-20")) ~ "wintering",
    TRUE ~ life_stage
  )) %>%
  filter(n_bursts_in_day >= 3) %>% #filter for number of IMU bursts per day.  
  as.data.frame()

#-------------------------------------------------------------------------------------
# STEP2: some plotting
#-------------------------------------------------------------------------------------

#### ridgelines for laterality in banking angle
ggplot(laterality_8sec_LS, aes(x = laterality_bank, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()

#### ridgelines for laterality in heading
ggplot(laterality_8sec_LS, aes(x = laterality_heading, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()

########old stuff
X11(width = 12, height = 12)
ggplot(data = pre_winter, aes(x = unique_date, y = laterality_bank)) +
  geom_point() +
  geom_smooth( method = "gam") +
  facet_wrap(vars(individual_local_identifier))

ggplot(pre_winter, aes(x = unique_date, y = laterality_bank, group = individual_local_identifier)) +
  geom_point(alpha = 0.5) +  # Add points with some transparency
  geom_smooth(method = "loess", se = FALSE) +  # Add smooth curves for each individual
  labs(title = "Daily Laterality Index",
       x = "Date",
       y = "Laterality Index") +
  theme_minimal() +
  theme(legend.position = "right") +
  scale_color_viridis_d() 

#look at a linear model
library(lme4)

# Assuming `individual` is a factor indicating individual subjects
model <- lmer(laterality_bank ~ (1|individual_local_identifier), data = pre_winter)
model <- lmer(laterality_bank ~ unique_date + (1|individual_local_identifier), data = pre_winter)
summary(model)

library(rptR)

rpt_model <- rpt(laterality_bank ~ (1|individual_local_identifier), grname = "individual_local_identifier", data = pre_winter, datatype = "Gaussian")
summary(rpt_model)


####gam

pre_winter <- pre_winter %>% 
  mutate(unique_date = as.Date(unique_date),
         individual_local_identifier = as.factor(individual_local_identifier))

gam_model <- gam(laterality_bank ~ s(unique_date, by = individual_local_identifier, bs = "fs", k = 10) + individual_local_identifier, 
                 data = pre_winter, method = "REML")



 
