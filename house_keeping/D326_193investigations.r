#The curious case of bird D326_193: is it dead?
#Elham Nourani PhD.
#Nov 20, 2023. Konstanz, DE

library(tidyverse)
library(lubridate)
library(move2)
library(sf)
library(mapview)
library(terra)
library(lwgeom)
library(rgl)
library(parallel)
library(viridis)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

#STEP 1: download gps data for all individuals (whole study) -----------------------------

creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
HB_id <- movebank_get_study_id("European Honey Buzzard_Finland")

all_sensors <- movebank_retrieve(study_id = 2201086728, sensor_type_id = c("gps", "acceleration", "orientation"), individual_local_identifier = "D329_012", 
                         timestamp_start = as.POSIXct("2023-10-15 00:00:00"),
                         timestamp_end = as.POSIXct("2023-11-01 00:00:00"),
                         entity_type = "event",  attributes = "all") 


gps <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "gps", individual_local_identifier = "D329_012", 
                                 timestamp_start = as.POSIXct("2023-10-15 00:00:00"),
                                 timestamp_end = as.POSIXct("2023-11-01 00:00:00"),
                                 entity_type = "event",  attributes = "all") 


acc <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "acceleration", individual_local_identifier = "D329_012", 
                                 timestamp_start = as.POSIXct("2023-10-15 00:00:00.000"),
                                 timestamp_end = as.POSIXct("2023-11-01 00:00:00.000"),
                                 entity_type = "event",  attributes = "all") 



or <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "orientation", individual_local_identifier = "D329_012", 
                                 timestamp_start = as.POSIXct("2023-10-15 00:00:00"),
                                 timestamp_end = as.POSIXct("2023-11-01 00:00:00"),
                                 entity_type = "event",  attributes = "all") 



#look at n of imu data for a high-res burst: 
one_burst <- gps %>%  filter(between(timestamp, as.POSIXct("2023-10-16 10:36:00", tz = "UTC"), as.POSIXct("2023-10-16 10:42:00", tz = "UTC"))) %>% 
  as.data.frame()

or_burst <- or %>% filter(between(timestamp, as.POSIXct("2023-10-16 10:36:00", tz = "UTC"), as.POSIXct("2023-10-16 10:42:00", tz = "UTC"))) %>% 
  as.data.frame()

acc_burst <- acc %>% filter(between(timestamp, as.POSIXct("2023-10-16 10:36:00", tz = "UTC"), as.POSIXct("2023-10-16 10:42:00", tz = "UTC"))) %>% 
  as.data.frame()

#calculate the n of unique entries for each second of IMU


acc <- acc %>%
  mutate( eobs_acceleration_g = purrr::map_chr(
      strsplit(as.character(eobs_accelerations_raw), " "),
      ~ as.character(unlist(.x) %>% as.numeric() %>% g_transform()) %>% str_c(collapse = " ")
    ))



mapview(gps_sf[2000:2700,])


############ look at the merged gps_IMU data

or_w_gps <- readRDS("GPS_matched_orientation_Nov23_2023birds.rds")
sample_ind_or <- or_w_gps[[3]]

acc_w_gps <- readRDS("GPS_matched_ACC_Nov23.rds")
sample_ind_acc <- acc_w_gps[[28]]

#high-res burst with concurrent gps and IMU

hr_burst_or <- sample_ind_or  %>% filter(between(timestamp, as.POSIXct("2023-10-16 10:36:00", tz = "UTC"), as.POSIXct("2023-10-16 10:42:00", tz = "UTC"))) %>%
  as.data.frame()

hr_burst_acc <- sample_ind_acc  %>% filter(between(timestamp, as.POSIXct("2023-10-16 10:36:00", tz = "UTC"), as.POSIXct("2023-10-16 10:42:00", tz = "UTC"))) %>%
  as.data.frame()

#IMU with no concurrent GPS

lr_burst_or <- sample_ind_or %>% filter(between(timestamp, as.POSIXct("2023-10-16 13:20:00", tz = "UTC"), as.POSIXct("2023-10-16 13:40:00", tz = "UTC"))) %>%
  as.data.frame()

lr_burst_acc <- sample_ind_acc %>% filter(between(timestamp, as.POSIXct("2023-10-16 13:20:00", tz = "UTC"), as.POSIXct("2023-10-16 13:40:00", tz = "UTC"))) %>%
  as.data.frame()

#conclusion: 10 sec of data per row for mag and quat (both for lr and hr). ACC is fine though (20 sec per row)

#BUT, how was it in the old data?
#open orientation for the old data
old_or <- readRDS("GPS_matched_orientation_Nov23.rds")[[1]] #here too. there are 10 seconds of data per row


