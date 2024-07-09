#This code downloads honey buzzard IMU data, transforms the units of ACC and processes the quaternion data. 
#Then it finds the nearest GPS point to the IMU recordings (modified from 01b_MARG_prep23.r)
#Elham Nourani PhD.
#Apr 15. 2024. Konstanz, DE. 

library(move2)
library(tidyverse)
library(lubridate)
library(mapview)
library(parallel)
library(sf)
library(circular)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#source functions for wind direction
source("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/lap_paper/AnEnvIPaper/data_prep/EnvironmentalData/airspeed_windsupport_crosswind.R")

#STEP 1: download all IMU data -------------------------------------------------

creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
HB_id <- movebank_get_study_id("European Honey Buzzard_Finland")

mag <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "magnetometer", 
                         entity_type = "event",  attributes = "all")

quat <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "orientation", 
                          entity_type = "event",  attributes = "all")

acc <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "acceleration", 
                         entity_type = "event",  attributes = "all")

#put the quat and mag together
or <- mag %>% 
  full_join(quat, by = c("individual_local_identifier", "tag_local_identifier", "study_id", "tag_id", "individual_taxon_canonical_name", "timestamp",
                         "eobs_start_timestamp", "eobs_key_bin_checksum", "data_decoding_software", "individual_id", "import_marked_outlier", "visible")) %>% 
  as.data.frame()

rm(mag,quat)

saveRDS(or, "all_orientation_apr15_24.rds")
saveRDS(acc, "all_acceleration_apr15_24.rds")


#STEP 2: data processing: ACC (convert units to g) -------------------------------------------------

#write a ftn for g conversion
g_transform <- function(x) {
  (x-2048)*(1/1024) #slope according to the eobs manual
}

#extract values for each axis as a numeric vector to convert to g

#make sure the values fall between -2 and 2
#I have the same slope and intercept for all axes, so don't separate into 3 axis and just convert all raw values to g.
acc <- readRDS("all_acceleration_apr15_24.rds")

(st <- Sys.time())
acc_g <- acc %>%
  mutate(
    eobs_acceleration_g = purrr::map_chr(
      strsplit(as.character(eobs_accelerations_raw), " "),
      ~ as.character(unlist(.x) %>% as.numeric() %>% g_transform()) %>% str_c(collapse = " ")
    )
  ) %>% 
  as.data.frame()
Sys.time()-st #2.2 minutes

#saveRDS(acc_g, "2023_birds_acc_g_Nov23.rds")

saveRDS(acc_g, "all_acceleration_g_apr_24.rds") #not yet written

#STEP 3: data processing: Quaternions (calculate Euler angles) -------------------------------------------------

#source the imu conversion functions
source("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/HB_ontogeny/MS1_IMU_classification/00_imu_diy.r")

#open orientation dataset, subset for data with high sampling frequency (this should be post-fledging, migration, and wintering data for laterality tests)
or <- readRDS("all_orientation_apr15_24.rds") %>% 
  group_by(individual_local_identifier) %>% 
  arrange(timestamp, .by_group = T) %>% 
  mutate(timelag = c(0, diff(timestamp)),  #calculate time lag 
         imu_burst_id = cumsum(timelag > 1)) %>% #assign a unique burst ID every time timelag is larger than 1 second. The IDs will be unique only within one ind's data
  ungroup() %>% 
  as.data.frame() 

#subset or by the length of the bursts. Only keep those that are longer than 3 second
burst_size <- or %>% 
  group_by(individual_local_identifier, imu_burst_id) %>% 
  summarize(burst_length = n(), #This is not very informative... because the length is not in seconds... it's better to calculate the burst duration instead
            burst_duration = sum(timelag[-1])) %>% #the first value of timelag in each burst is a large value indicating timediff between this and the previous burst. remove it
  ungroup()
  
or_hfreq <- or %>% 
  right_join(burst_size %>% filter(burst_duration > 3)) #there are bursts with length of 0-3

rm(or);gc()

#calculate pitch, roll, and yaw for the whole dataset
(st_time <- Sys.time())
or_angles <- or_hfreq %>%
  rowwise() %>% 
  mutate(
    pitch = process_quaternions(orientation_quaternions_raw, ~ get.pitch(.x, type = "eobs")),
    yaw = process_quaternions(orientation_quaternions_raw, ~ get.yaw(.x, type = "eobs")),
    roll = process_quaternions(orientation_quaternions_raw, ~ get.roll(.x, type = "eobs"))
  ) %>% 
  ungroup() %>%  #make sure this doesnt mess up the next step
  #convert to degrees (make sure values are from -180 to +180)
  mutate(across(c("pitch", "yaw", "roll"), 
                ~ purrr::map_chr(
                  strsplit(as.character(.x), " "),
                  ~ as.character(unlist(.) %>% as.numeric() %>% deg()) %>% str_c(collapse = " ") #make sure the values are SIGNED
                ),
                .names = "{col}_deg"
  )) %>% 
  as.data.frame()
Sys.time() - st_time #1.83 hours for the whole HB dataset

saveRDS(or_angles, file = "quat_angles_apr24.rds")


#modify the data to have one row per second. this step is not necessary, unless I want to merge the acc and orientation
#Useful for summarizing the quat angles for each second. This way I can later summarize the values for each burst, depending on the burst length of interest

or_angles <- readRDS("quat_angles_apr24.rds")

(st_time <- Sys.time())
or_seconds <- or_angles%>% 
  group_by(individual_local_identifier, floor_date(timestamp, "second")) %>% 
  #create a long character string containing all values for each second
  summarize(across(c(roll, pitch, yaw, roll_deg, pitch_deg, yaw_deg), ~paste(., collapse = " ")), 
            #retain the first value of other important columns
            across(c(study_id, individual_taxon_canonical_name, orientation_quaternions_sampling_frequency,
                     tag_local_identifier, timestamp, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
            .groups = "keep") %>%
  ungroup() %>%
  as.data.frame()

Sys.time() - st_time #10 hrs

saveRDS(or_seconds, file = "quat_angles_secs_apr24.rds")


#For each axis (roll, pitch, yaw), calculate: mean, mean_abs, max, max_abs, min, min_abs, sd, sum, sum_abs

or_seconds <- readRDS("quat_angles_secs_apr24.rds") 


or_burst_ls <- or_seconds %>%
  group_by(individual_local_identifier) %>% 
  arrange(timestamp, .by_group = T) %>% 
  #the burst assignments are correct, but burst lengths are different depending on the year, etc. for the sake of uniformity, assign new burst ids that break up the data into 8 seconds
  mutate(timelag = c(0, diff(timestamp)),  #calculate time lag 
         burst_id_8sec = cumsum(timelag > 8),
         unique_burst_id = paste0(individual_local_identifier, "_", burst_id_8sec)) %>% 
  ungroup() %>% 
  as.data.frame()

or_burst_ls <- split(or_burst_ls, or_burst_ls$unique_burst_id)

rm(or_seconds); gc()

(start_t <- Sys.time())
or_burst_summaries <- lapply(or_burst_ls, function(x){
  
  #convert character strings to numeric values
  num_ls <- list(
    roll = strings_to_numeric(x$roll_deg),
    pitch = strings_to_numeric(x$pitch_deg),
    yaw = strings_to_numeric(x$yaw_deg)
  )
  
  #calculate summary statistics for the 8-second burst.
  angle_s <- angle_summaries(yaw = num_ls$yaw,
                             pitch = num_ls$pitch,
                             roll = num_ls$roll) %>% 
    mutate(start_timestamp = head(x$timestamp,1),
           end_timestamp = tail(x$timestamp, 1)) %>% 
    bind_cols(x %>% select(c(study_id, individual_local_identifier, individual_taxon_canonical_name, orientation_quaternions_sampling_frequency,
                            tag_local_identifier, tag_id, imu_burst_id, unique_burst_id)) %>% slice(1)) %>% 
    as.data.frame() #aggregated dataframe (with one row) for one burst
  
}) %>% 
  bind_rows()
Sys.time() - start_t #1.03 hours

saveRDS(or_burst_summaries, file = "matched_GPS_IMU/quat_summaries_8secs_Jul24.rds")


# STEP 4: find nearest GPS fix to each quat/mag burst -------------------------------------------------

#open gps data: segmented and annotated
gps_ls <- list.files("gps_seg_apr24", full.names = T) %>% #these are sf files
  map(readRDS) %>% 
  bind_rows() %>% 
  mutate(row_id = row_number(), #assign a row id to be able to cbind the annotated data with the original data later on
         location_long = st_coordinates(.)[,1],
         location_lat = st_coordinates(.)[,2],
         row_id = row_number()) %>% #add row id for ease of merging with annotation data later on
  st_drop_geometry() %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

#open orientation data, with one row per second
or_ls <- readRDS("quat_angles_secs_apr24.rds") %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

# Define a function to find the closest GPS information and associate it with orientation data
find_closest_gps <- function(or_data, gps_data, time_tolerance = 10 * 60) {
  map_df(1:nrow(or_data), function(h) {
    or_row_time <- or_data[h, "timestamp"]
    gps_sub <- gps_data %>%
      filter(between(timestamp, or_row_time - time_tolerance, or_row_time + time_tolerance))
    
    if (nrow(gps_sub) >= 1) {
      time_diff <- abs(difftime(gps_sub$timestamp, or_row_time, units = "secs"))
      min_diff <- which.min(time_diff)
      or_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps", "ground_speed_closest_gps", 
                   "heading_closest_gps", "row_id","flight_type_sm2", "flight_type_sm3", "track_flight_seg_id")] <- 
        gps_sub[min_diff, c("timestamp", "location_long", "location_lat", "height_above_ellipsoid", "ground_speed", "heading", "row_id", 
                            "flight_clust_sm2", "flight_clust_sm3", "track_flight_seg_id")]
    } else {
      or_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps", "ground_speed_closest_gps", 
                   "heading_closest_gps", "row_id", "flight_type_sm2", "flight_type_sm3", "track_flight_seg_id")] <- NA
    }
    return(or_data[h, ])
  })
}

# Create a list of data frames with or data and associated GPS information
(b <- Sys.time())
or_w_gps <- map2(or_ls, gps_ls, ~ find_closest_gps(or_data = .x, gps_data = .y))
Sys.time() - b # 3.8 hrs for orientation with one sec per row

#saveRDS(or_w_gps, "GPS_matched_orientation_Nov23.rds") # a list of data frames
saveRDS(or_w_gps, "GPS_matched_orientation_Apr24.rds")


#add a column comparing or and gps timestamps. then save one file per individual. also limit to migratory season!! 1.09 - 30.10
lapply(or_w_gps, function(x){
  x2 <- x %>% 
    mutate(imu_gps_timediff_sec = if_else(is.na(timestamp_closest_gps), NA, difftime(timestamp, timestamp_closest_gps, units = "secs") %>%  as.numeric()))
  
  saveRDS(x2, 
            file = paste0("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_gps_quat/",
                          x2$individual_local_identifier[1], "_quat_w_gps.rds"))
})


# STEP 5: find nearest GPS fix to each ACC burst -------------------------------------------------

#list all acc files
acc <- readRDS("all_acceleration_g_apr_24.rds") %>% 
  as.data.frame()

acc_ls <- split(acc, acc$individual_local_identifier) #make sure this is a list of dataframes, not tibbles. make sure the order of individuals in the acc and gps lists are the same

# Define a function to find the closest GPS information and associate it with ACC data
find_closest_gps <- function(acc_data, gps_data, time_tolerance = 10 * 60) {
  map_df(1:nrow(acc_data), function(h) {
    acc_row_time <- acc_data[h, "timestamp"]
    gps_sub <- gps_data %>%
      filter(between(timestamp, acc_row_time - time_tolerance, acc_row_time + time_tolerance))
    
    if (nrow(gps_sub) >= 1) {
      time_diff <- abs(difftime(gps_sub$timestamp, acc_row_time, units = "secs"))
      min_diff <- which.min(time_diff)
      acc_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps", "ground_speed_closest_gps", "heading_closest_gps", "row_id", 
                    "flight_type_sm2", "flight_type_sm3", "track_flight_seg_id")] <- 
        gps_sub[min_diff, c("timestamp", "location_long", "location_lat", "height_above_ellipsoid", "ground_speed", "heading", "row_id", 
                            "flight_clust_sm2", "flight_clust_sm3", "track_flight_seg_id")]
    } else {
      acc_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps", "ground_speed_closest_gps", "heading_closest_gps", "row_id", 
                    "flight_type_sm2", "flight_type_sm3", "track_flight_seg_id")] <- NA
    }
    return(acc_data[h, ])
  })
}

# Create a list of data frames with ACC data and associated GPS information
(b <- Sys.time())
acc_w_gps <- map2(acc_ls, gps_ls, ~ find_closest_gps(acc_data = .x, gps_data = .y))
Sys.time() - b # 3.5 hours

saveRDS(acc_w_gps, "GPS_matched_ACC_May24_allbirds.rds")

# Now you have acc_w_gps, which is a list of data frames, where each data frame contains ACC data with associated GPS information.


#add a column comparing acc and gps timestamps. then save one file per individual. also limit to migratory season!! 1.09 - 30.10
lapply(acc_w_gps, function(x){
  x2 <- x %>% 
    mutate(acc_gps_timediff_sec = if_else(is.na(timestamp_closest_gps), NA, difftime(timestamp, timestamp_closest_gps, units = "secs") %>%  as.numeric()))
  
  write.csv(x2, 
            file = paste0("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_gps_acc/",
                          x2$individual_local_identifier[1], "_acc_w_gps.csv"))
})

############ Check one file from 2023 to make sure the matching worked ###############

sample <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_orientation/D329_012_quat_mag_w_gps.csv") #2023 bird
sample2 <-read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_orientation/D323_154_quat_mag_w_gps.csv") #2023 bird
