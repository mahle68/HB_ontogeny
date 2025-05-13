# script for preparing the data for Safi et al 2025.
# Elham Nourani, PhD. enourani@ab.mgp.de

#This code downloads honey buzzard IMU data, transforms the units of ACC and processes the quaternion data. 
#Then it finds the nearest GPS point to the IMU recordings

library(move2)
library(tidyverse)
library(lubridate)
library(mapview)
library(parallel)
library(sf)
library(circular)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#source functions for wind direction
source("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/lap_paper/AnEnvIPaper/data_prep/EnvironmentalData/airspeed_windsupport_crosswind.R")
#source the imu conversion functions
source("00_imu_functions.r")


#---------------------------------------------------------------------------------
## Step 0: Download all IMU data                                             #####
#---------------------------------------------------------------------------------

#Make sure you have set up your credentials following Move2 tutorials.
#Data used in this study was downloaded on 15.04.2024. Use this timestamp as the cutoff point when downloading the data to reproduce the results of this study.

#creds <- movebank_store_credentials(username = "your_user_name", rstudioapi::askForPassword())

mag <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "magnetometer", 
                         entity_type = "event",  attributes = "all",
                         timestamp_end = as.POSIXct("2024-04-15 00:00:00"))

quat <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "orientation", 
                          entity_type = "event",  attributes = "all",
                          timestamp_end = as.POSIXct("2024-04-15 00:00:00"))

acc <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "acceleration", 
                         entity_type = "event",  attributes = "all",
                         timestamp_end = as.POSIXct("2024-04-15 00:00:00"))

#put the quat and mag together
or <- mag %>% 
  full_join(quat, by = c("individual_local_identifier", "tag_local_identifier", "study_id", "tag_id", "individual_taxon_canonical_name", "timestamp",
                         "eobs_start_timestamp", "eobs_key_bin_checksum", "data_decoding_software", "individual_id", "import_marked_outlier", "visible")) %>% 
  as.data.frame()

#remove the originial files to save RAM
rm(mag, quat)

#saveRDS(or, "all_orientation_apr15_24.rds")
#saveRDS(acc, "all_acceleration_apr15_24.rds")


#---------------------------------------------------------------------------------
## Step 1: Process the ACC                                                   #####
#---------------------------------------------------------------------------------

#### ----------------------- convert units to g 

#write a ftn for g conversion
g_transform <- function(x) {
  (x-2048)*(1/1024) #slope according to the eobs manual
}

#make sure the values fall between -2 and 2
#For my tag generation, I have the same slope and intercept for all axes, so there is no need to separate into 3 axis before converting all raw values to g.
acc_g <- acc %>%
  mutate(
    eobs_acceleration_g = purrr::map_chr(
      strsplit(as.character(eobs_accelerations_raw), " "),
      ~ as.character(unlist(.x) %>% as.numeric() %>% g_transform()) %>% str_c(collapse = " ")
    )
  ) %>% 
  as.data.frame()

#check the number of acc records. we have 1.2 seconds of data per row
acc_n <- acc_g %>%
  rowwise() %>%
  mutate(acc_g_length = length(strings_to_numeric(eobs_acceleration_g)),
         acc_raw_length = length(strings_to_numeric(eobs_accelerations_raw))
  )

#### ----------------------- calculate vedba for each row of data (1.2 seconds)
acc_g <- readRDS("all_acceleration_g_apr_24.rds") %>% 
  rowwise() %>% 
  mutate(DBA_ls = list(DBA_ftn(eobs_acceleration_g)),
    VeDBA = DBA_ls[1],
    ODBA = DBA_ls[2],
    acc_g_length = length(strings_to_numeric(eobs_acceleration_g)),
    acc_raw_length = length(strings_to_numeric(eobs_accelerations_raw))) %>%
  ungroup() %>%
  select(-DBA_ls) %>% 
  as.data.frame()

#---------------------------------------------------------------------------------
## Step 2: Process the Quaternions (calculate Euler angles)                  #####
#---------------------------------------------------------------------------------

#subset the orientation dataframe for data with high sampling frequency (this should be post-fledging, migration, and wintering data for laterality tests.)
or <- or %>% 
  group_by(individual_local_identifier) %>% 
  arrange(timestamp, .by_group = T) %>% 
  mutate(timelag = c(0, diff(timestamp)),  #calculate time lag 
         imu_burst_id = cumsum(timelag > 1)) %>% #assign a unique burst ID every time timelag is larger than 1 second. The IDs will be unique only within one individual's data
  ungroup() %>% 
  as.data.frame()

#subset orientation by the length of the bursts. Only keep those that are longer than 3 second.
#first calculate burst durations
burst_size <- or %>% 
  group_by(individual_local_identifier, imu_burst_id) %>% 
  summarize(burst_duration = sum(timelag[-1])) %>% #the first value of timelag in each burst is a large value indicating timediff between this and the previous burst. remove it
  ungroup()

#remove short bursts
or_hfreq <- or %>% 
  right_join(burst_size %>% filter(burst_duration > 3)) #there are bursts with length of 0-3

#remove the original orientation dataframe and clean up RAM
rm(or);gc()

#calculate pitch, roll, and yaw for the whole dataset

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


saveRDS(or_angles, file = "quat_angles_apr24.rds")

###### summarize data over each second
#modify the data to have one row per second. this step is not absolutely necessary
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

#for each second, calculate the mean, min, max, sd, cumsum for each axis
or_seconds <- readRDS("quat_angles_secs_apr24.rds") 

(start_t <- Sys.time())
seconds_summaries <- or_seconds %>%
  rowwise() %>%
  mutate(
    angle_summary = list(angle_summaries(
      yaw = strings_to_numeric(yaw_deg),
      pitch = strings_to_numeric(pitch_deg),
      roll = strings_to_numeric(roll_deg)
    ))
  ) %>%
  ungroup() %>%
  unnest(angle_summary) %>% 
  as.data.frame()

Sys.time() - start_t # 3.7 hrs

saveRDS(seconds_summaries, file = "quat_summaries_1sec_Jul24.rds")


# STEP 4: find nearest GPS fix to each quat/mag burst -------------------------------------------------

#open gps data: segmented and annotated
gps_ls <- list.files("gps_seg_apr24", full.names = T) %>% #these are sf files
  map(readRDS) %>% 
  bind_rows() %>% 
  mutate(row_id = row_number(), #assign a row id to be able to cbind the annotated data with the original data later on
         location_long = st_coordinates(.)[,1],
         location_lat = st_coordinates(.)[,2]) %>% 
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


########## STEP 4.1 : bind one-second orientation data with summarized values per second with the matched GPS-Orientation data (Jul. 11. 2024) ##########

#orientation summarized per second: 
seconds_summaries <- readRDS("quat_summaries_1sec_Jul24.rds")

#GPS-matched orientation:
or_w_gps <- readRDS("GPS_matched_orientation_Apr24.rds") %>% 
  bind_rows() %>% 
  select(-`floor_date(timestamp, "second")`)

#bind the two based on individual ID and imu_burst_id

or_summaries_w_gps <- or_w_gps %>% 
  full_join(seconds_summaries, 
            by = intersect(names(or_w_gps), names(seconds_summaries)))

saveRDS(or_summaries_w_gps, file = "matched_GPS_IMU/GPS_matched_or_w_summaries_Jul24.rds")

########## STEP 4.2 : summarize quat data over 8-second bursts (Jul. 17. 2024) ##########

#assign a unique id for each 8-second burst (most of the honey buzzard IMU was collected in 8-second bursts). break up longer bursts into 8 second subsets
or_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_Jul24.rds") %>% 
  drop_na(individual_local_identifier)

#Add a new column with the number of rows of each imu_burst. We have roughly one second of data per row. 
burst_n <- or_summaries_w_gps %>% 
  group_by(individual_local_identifier, imu_burst_id) %>% 
  mutate(n_of_rows = n()) %>% 
  ungroup() %>% 
  as.data.frame()

#Add a new column with new burst IDs for bursts that are longer than 9 seconds. For the rest, keep the original burst IDs 
or_seconds_8_sec_id <- burst_n %>% 
  filter(n_of_rows > 9) %>% 
  group_by(individual_local_identifier, imu_burst_id) %>% 
  mutate(burst_id_8sec = paste0(individual_local_identifier, "_", imu_burst_id, "_", ceiling(row_number() / 8))) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  #bind it back to the rest of the data
  full_join(burst_n %>% filter(n_of_rows <= 9)) %>% 
  mutate(burst_id_8sec = ifelse(is.na(burst_id_8sec), 
                                paste0(individual_local_identifier, "_", imu_burst_id, "_", 1), burst_id_8sec)) %>% #assign unique 8-sec burst ID to bursts that were less than 9 seconds long to begin with. use 1 because they each have one sub-burst.
  as.data.frame()

saveRDS(or_seconds_8_sec_id, "matched_GPS_IMU/GPS_matched_or_w_summaries_8secIDs_Jul24.rds") #also write this to file, so I can easily access the raw values for each 8sec-burst just in case


# calculate summaries of angles for each 8-sec burst
(start_time <- Sys.time())
or_8sec_summaries <- or_seconds_8_sec_id %>% 
  group_by(burst_id_8sec) %>% #this grouping variable has unique values for individuals and bursts. 
  arrange(timestamp, .by_group = T) %>% 
  #aggregate summarized angles over 8 seconds
  summarize(across(c(study_id, individual_taxon_canonical_name, individual_local_identifier, orientation_quaternions_sampling_frequency,
                     tag_local_identifier, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
            across(c(yaw_sd, roll_sd, pitch_sd,
                     yaw_mean, roll_mean, pitch_mean,
                     yaw_min, roll_min, pitch_min,
                     yaw_max, roll_max, pitch_max), 
                   ~mean(., na.rm = T),
                   .names = "mean_{.col}"), #rename the column, because these values are the mean of the sd, max, min values across the 8 seconds
            across(c(cumulative_yaw, cumulative_roll, cumulative_pitch),
                   ~sum(., na.rm = T),
                   .names = "{.col}_8sec"),
            across(c(timestamp, timestamp_closest_gps, location_lat_closest_gps, location_long_closest_gps, height_above_ellipsoid_closest_gps, ground_speed_closest_gps, heading_closest_gps), 
                   ~head(.,1),
                   .names = "start_{.col}"),
            across(c(flight_type_sm2, flight_type_sm3, track_flight_seg_id),
                   ~Mode(.), #make sure the Mode function is defined
                   .names = "{.col}_Mode"),
            end_timestamp = tail(timestamp,1),
            .groups = "keep") %>% 
  ungroup %>% 
  as.data.frame()
(Sys.time() - start_time) #17.3 minutes

saveRDS(or_8sec_summaries, "matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_Jul24.rds") #this doesn't have the NA individual


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
