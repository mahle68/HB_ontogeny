#This code runs tests for whether honey buzzards have a preference for left or right bank angles during different life cycle stage
#at the scale of each 8-second burst
#Elham Nourani PhD.
#July 24. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)
library(viridis)
library(lme4)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")


#-------------------------------------------------------------------------------------
# STEP1: annotate data with life cycle stage
#-------------------------------------------------------------------------------------
#open data with laterality index calculated

#this file contains laterality index calculated per burst. NOT filtered for circling flight
laterality_1sec <- readRDS("laterality_index_per_8sec_burst.rds")

#open meta-data to calculate day since tagging
meta_data <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/EHB_metadata - Sheet1.csv") %>% 
  mutate(deployment_dt_utc = as.POSIXct(deployment_dt_utc, tz = "UTC"))

laterality_1sec_days <- laterality_1sec %>%
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  full_join(meta_data %>% select(ring_ID, deployment_dt_utc), by = c("individual_local_identifier" = "ring_ID")) %>% 
  rowwise() %>% 
  mutate(days_since_tagging = floor(difftime(start_timestamp, deployment_dt_utc, unit = "days"))) %>% 
  ungroup() %>% 
  as.data.frame()

saveRDS(laterality_1sec_days, "laterality_index_per_8sec_burst_days_since.rds")

# #use Ellen's life cycle stage estimations
# life_cycle <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/EHB_metadata_mig_dates_Fin22.csv") %>% 
#   mutate_at(c("migration_start", "migration_end"), as.Date) %>% 
#   filter(nest_country == "Finland")
# 
# #create new columns for life cycle stage AND day and week since tagging (proxy for age) 
# (start_time <- Sys.time())
# laterality_1sec_LS <- laterality_1sec %>%
#   full_join(life_cycle %>% select(ring_ID, deployment_dt_utc, migration_start, migration_end), by = c("individual_local_identifier" = "ring_ID")) %>% 
#   rowwise() %>% 
#   mutate(life_stage = case_when(
#     between(start_timestamp, migration_start, migration_end) ~ "migration",
#     start_timestamp < migration_start ~ "post-fledging",
#     start_timestamp > migration_end ~ "wintering",
#     TRUE ~ NA_character_
#   ),
#   day_since_tagging = floor(difftime(start_timestamp, as.POSIXct(deployment_dt_utc), unit = "days"))) %>% 
#   ungroup() %>% 
#   mutate(life_stage = case_when(
#     is.na(life_stage) & start_timestamp < as.Date("2023-09-15") ~ "post-fledging",
#     is.na(life_stage) & start_timestamp > as.Date("2023-10-20") ~ "migration",
#     is.na(life_stage) & between(start_timestamp, as.Date("2023-09-15"), as.Date("2023-10-20")) ~ "wintering",
#     TRUE ~ life_stage
#   )) %>%
#   as.data.frame()
# Sys.time() - start_time #9 minutes
# 
# saveRDS(laterality_1sec_LS, "laterality_index_per_8sec_burst_LS.rds")

# #filter for circling flight only
# laterality_circling <- laterality_1sec_LS %>% 
#   #mutate(burst_dur2 = difftime(end_timestamp, start_timestamp)) %>% ## OR use the n of records for this. that would be the number of rows
#   #filter(burst_dur2 > 6.5 &  #remove bursts that are shorter than 6.5 seconds
#   filter(n_records > 6.5 &  #remove bursts that are shorter than 6.5 seconds
#            cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45) %>%  #only keep thermalling flight. use a threshold of 45 degrees in 8 seconds.
#   as.data.frame()
  

#-------------------------------------------------------------------------------------
# STEP2: some plotting
#-------------------------------------------------------------------------------------

#### ridgelines for laterality in banking angle
ggplot(laterality_circling, aes(x = laterality_bank, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()

ggplot(laterality_circling, aes(x = laterality_bank, y = individual_local_identifier, fill = life_stage)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("post-fledging" = "yellow", "migration" = "#c4c4fc", "wintering" = "red")) +
  theme_minimal() +
  labs(fill = "Life Stage")


#### ridgelines for laterality in heading
ggplot(laterality_circling, aes(x = laterality_heading, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()


