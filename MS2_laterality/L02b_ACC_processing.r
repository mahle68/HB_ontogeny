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

saveRDS(acc_migr, file = "migr_vedba.rds")


