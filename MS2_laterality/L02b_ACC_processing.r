#This code calculates hourly and daily sumamries for vedba for European honey buzzards and calculates daily summaries for exploring migration performance in Safi et al 2025
#Elham Nourani PhD.
#Dec 6. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)
library(viridis)


#open acc data prepared in 01b_imu_processing.r
#this is mean odba and vedba for each row (1.2 seconds) of data
acc_g <- readRDS("all_acceleration_g_dba_apr_24.rds") %>% 
  mutate(unique_date = as.Date(timestamp),
         ind_day = paste0(individual_local_identifier, "_", as.character(unique_date)))


#subset for days of migration (as identified in L04c_migr_metrics.r) and summarize for each day
gps_1hr <- readRDS("gps_data_LS_migr_metrics.rds") 

#extract days of migration
migr_days <-gps_1hr %>% 
  distinct(ind_day)


acc_migr <- acc_g %>% 
  filter(ind_day %in% migr_days$ind_day) %>%  #only keep acc data for migration days 
  mutate(dt_1hr = round_date(timestamp, "1 hour")) %>%  #the hourly subset will be used just to pick days that have 12-14 hours of data. 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  mutate(hrly_mean_vedba = mean(VeDBA, na.rm = T),
         hrly_max_vedba = max(VeDBA, na.rm = T),
         hrly_min_vedba = min(VeDBA, na.rm = T),
         hrly_IQR_vedba = quantile(VeDBA, prob = 0.75) - quantile(VeDBA, prob = 0.25)) %>%  #calculate interquartile range... this includes non-flight data tooo.... so maybe not so useful
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  mutate(daily_mean_vedba = mean(VeDBA, na.rm = T),
         daily_max_vedba = max(VeDBA, na.rm = T),
         daily_min_vedba = min(VeDBA, na.rm = T),
         daily_IQR_vedba = quantile(VeDBA, prob = 0.75) - quantile(VeDBA, prob = 0.25)) %>%  #calculate interquartile range... this includes non-flight data tooo.... so maybe not so useful
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  slice(1) #just keep one row per hour


#append to the gps data summarized for migration period
migr_hourly <- gps_1hr %>% #this is an sf object. keep the lat and long
  mutate(location_long = st_coordinates(.)[1],
         location_lat = st_coordinates(.)[2]) %>% 
  select(individual_id, deployment_id, tag_id, study_id, individual_local_identifier, tag_local_identifier, individual_taxon_canonical_name, location_long, location_lat,
         ind_day, migration_start, migration_end, first_exploration, life_stage, dt_1hr, geometry, height_msl, hrly_step_length, daily_distance, daily_avg_speed, daily_avg_altitude) %>% 
  full_join(acc_migr %>% dplyr::select(individual_local_identifier, ind_day, dt_1hr, hrly_mean_vedba, hrly_max_vedba, hrly_min_vedba, hrly_IQR_vedba,
                                       daily_mean_vedba, daily_max_vedba, daily_min_vedba, daily_IQR_vedba)) %>% 
  as.data.frame()

saveRDS(migr_hourly, file = "hourly_migr_metrics_gps_vedba.rds") #still need to add summarized values for yaw, pitch and potentially vertical speed to this.

