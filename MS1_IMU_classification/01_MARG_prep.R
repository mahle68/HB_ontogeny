#This code downloads honey buzzard IMU data, transforms the units, and calculates ACC related features. Then it finds the neares GPS point to the IMU recordings
#the focus is on two individuals at the exploration phase: 
#Elham Nourani PhD.
#Feb 15. 2023. Konstanz, DE.

library(move2)
#library(moveACC)
library(tidyverse)
library(lubridate)
library(mapview)
library(parallel)
library(sf)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#download movebank data

inds <- c("D320_475", "D324_512")

creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
HB_id <- movebank_get_study_id("European Honey Buzzard_Finland")

movebank_retrieve(study_id = 2201086728, entity_type= "tag_type")

mag <- movebank_retrieve(study_id = 2201086728, sensor_type_id = c("magnetometer", "orientation"), 
                         individual_local_identifier = inds,  entity_type = "event",  attributes = "all", 
                         timestamp_start = as.POSIXct("2022-09-01 00:00:00"), timestamp_end = as.POSIXct("2022-10-30 00:00:00"))

# STEP 1: download data and convert units to g -------------------------------------------------
#download raw acc values
acc <- getMovebankNonLocationData(EHB_FN_id, sensorID = 2365683, login = creds, removeDuplicatedTimestamps = T)

#write a ftn for g conversion
g_transform <- function(x) {
  (x-2048)*(1/1024) #slope according to the eobs manual
}

#extract values for each axis as a numeric vector to convert to g

#make sure the values fall between -2 and 2
#i have the same slope and intercept for all axes, so don't separate into 3 axis and just convert all raw values to g.
acc_g <- lapply(split(acc,seq(nrow(acc))), function(burst){
  
  g_data <- unlist(strsplit(burst %>% dplyr::select("eobs_accelerations_raw") %>%  pull(), " ")) %>% 
    as.numeric() %>% 
    g_transform()
  
  burst <- burst %>% 
    mutate(eobs_acceleration_g = str_c(as.character(g_data), collapse = " "))

}) %>% 
  reduce(rbind)

saveRDS(acc_g, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/acc_2inds.rds")


acc_g <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/acc_2inds.rds")

# STEP 2: estimate metrics -------------------------------------------------
acc_g <- readRDS("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/acc_2inds.rds") %>% 
  rename(accelerationTransformed = eobs_acceleration_g)

lapply(split(acc_g, acc_g$individual_local_identifier), function(x){
  
  #wing beat frequency and odba
  wave_acc <- ACCwave(x, transformedData = T, showProgress = T)
  
  #flight vs non-flight assignment based on amplitude and odba.
  wave_behav_acc <- WingBeatsSelection(wave_acc,forclustering = c("amplitude","odbaAvg"), minbeat = 0, maxbeat = max(wave_acc$beatsSec)) %>% 
    rename(flight_status = behavior, # the flight/non-flight column name is behavior. rename to flight_status
           individual_local_identifier = individualID,
           event_id = event.id,
           tag_id = tagID)
  
  #append to original acc data
  acc_metrics <- wave_behav_acc %>% 
    full_join(x, by = c("timestamp", "individual_local_identifier", "event_id")) %>% 
    rename(eobs_acceleration_g = accelerationTransformed)
  
  saveRDS(acc_metrics, file = paste0("Pritish_collab_IMU/animal_",unique(acc_metrics$tag_local_identifier),"_classifiedAcc.rds"))
})

#STEP 2: find nearest GPS fix to each ACC burst -------------------------------------------------

# List all ACC files
acc_ls <- lapply(list.files("Pritish_collab_IMU", pattern = "animal_", full.names = TRUE), readRDS)

# List of individual IDs
#IDs <- sapply(acc_files, function(x) str_sub(x, -22, -19))

# Open GPS data (downloaded in 02_GPS_segmentation.R)
gps_ls <- readRDS("gps_raw_2inds.rds") %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  as.list()

names(gps_ls) <- c("D320_475", "D324_512")

# Define a function to find the closest GPS information and associate it with ACC data
find_closest_gps <- function(acc_data, gps_data, time_tolerance = 10 * 60) {
  map_df(1:nrow(acc_data), function(h) {
    acc_row_time <- acc_data[h, "timestamp"]
    gps_sub <- gps_data %>%
      filter(between(timestamp, acc_row_time - time_tolerance, acc_row_time + time_tolerance))
    
    if (nrow(gps_sub) >= 1) {
      time_diff <- abs(difftime(gps_sub$timestamp, acc_row_time, units = "secs"))
      min_diff <- which.min(time_diff)
      acc_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps")] <- gps_sub[min_diff, c("timestamp", "location_long", "location_lat", "height_above_ellipsoid")]
    } else {
      acc_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps")] <- NA
    }
    return(acc_data[h, ])
  })
}

# Create a list of data frames with ACC data and associated GPS information
(b <- Sys.time())
acc_w_gps <- map2(acc_ls, gps_ls, ~ find_closest_gps(acc_data = .x, gps_data = .y))
Sys.time() - b # 16 minutes

saveRDS(acc_w_gps, "GPS_matched_ACC.rds")
# Now you have acc_w_gps, which is a list of data frames, where each data frame contains ACC data with associated GPS information.

#add a column comparing acc and gps timestamps. then save one file per individual. also limit to migratory season!! 1.09 - 30.10
lapply(acc_w_gps, function(x){
  x2 <- x %>% 
    filter(between(timestamp, as.POSIXct("2022-09-01 00:00:00", tz = "UTC"), as.POSIXct("2022-10-31 00:00:00", tz = "UTC"))) %>% 
    mutate(acc_gps_timediff_sec = if_else(is.na(timestamp_closest_gps), NA, difftime(timestamp, timestamp_closest_gps, units = "secs") %>%  as.numeric())) %>% 
    select(-tag_id.x) %>% 
    rename(tag_id = tag_id.y)
  
  write.csv(x2, 
            file = paste0("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/matched_gps_acc/",
                          x2$individual_local_identifier[1], "_acc_w_gps.csv"))
})


