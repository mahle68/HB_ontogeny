#This code downloads honey buzzard ACC data, transforms the units, and calculates ACC related features.
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
setwd("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")
setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")

#setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")
source("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/HB_ontogeny/ACCtoGPS_ftns.R") #matching gps and acc or vice versa

# focus on two individuals:
creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

#check for sensors available. Acc = 2365683; Mag = 77740402; Orientation: 819073350
getMovebankSensors(EHB_FN_id,login = creds)

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
#match acc to gps based on closest timestamp
# # STEP 2: stitch the gps and acc together -------------------------------------------------
# 
# source("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/HB_ontogeny/ACCtoGPS_ftns.R")
# acc_files <- list.files("ACCsegmentation_Mar23", full.names = T)
# gps_files <- list.files("GPSsegmentation_Mar23/classifiedData", full.names = T)
# 
# #list of individual IDs
# IDs <- sapply(acc_files, function(x) substr(x, 30,37))
# 
# 
# lapply(IDs, function(ind){
#   
#   acc_df <- readRDS(grep(ind, acc_files, fixed = T, value = T))
#   load(grep(ind, gps_files, fixed = T, value = T)) #HRdf_smooth
#   
#  #save these objects as rds in the previous steps
#   gps_df <- HRdf_smooth; rm(HRdf_smooth)
# 
#   # Create new ACC event column (using the tag id and the burst id assigned during the acc classification)
#   acc_df$new_accEventId <- as.numeric(paste0(acc_df$tag_id.x, acc_df$burstID))
#   
#   if(nrow(acc_df)==length(unique(acc_df$new_accEventId))){
#     #format timestamps
#     # gps_df$timestamp <- as.POSIXct(as.character(gps_df$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
#     # acc_df$timestamp <- as.POSIXct(as.character(acc_df$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
#     # Order both datasets by timestamp
#     acc_df <- acc_df[order(acc_df$timestamp),]
#     gps_df <- gps_df[order(gps_df$timestamp),]
#     # Extract column names to associate to the GPS dataset
#     accColsToAssociate <- c("beatsSec","amplitude","odbaAvg","odbaMedian","accAxes","numberSamplesPerAxis","burstDurationSecs",
#                             "samplingFreqPerAxis", "eobs_accelerations_raw", "accelerationTransformed", "flight_status")
#     
#     # Define time tolerance (e.g. 5 mins in seconds) to look for the closest ACC information and associate it to the gps data
#     timeTolerance <- 5*60 
#     accGps <- ACCtoGPS(GPSdata = gps_df, ACCdata = acc_df, accEventCol="new_accEventId",
#                        timeTolerance = timeTolerance,
#                        ColsToAssociate = accColsToAssociate) #columns of the acc data that you want to associate to the gps data
#     # Save the resulting dataset (merging gps and acc) for each individual
#     saveRDS(accGps, file = paste0("ACC_GPS_behav_class/",ind,"_gps&acc_behavClass.rds"))
#   }else{warning("ACC event id are not unique!")}
#   #return(accGps)
# })
# 
# # Bind all individuals and save final dataset
# accGps_df <- lapply(list.files("ACC_GPS_behav_class", full.names = T), readRDS) %>% 
#   reduce(rbind)
# 
# saveRDS(accGps_df, file = "ACC_GPS_behav_class/two_inds_gps&acc_behavClass.rds")
# 
# #exploration
# data <- readRDS( "ACC_GPS_behav_class/two_inds_gps&acc_behavClass.rds")
# 
# 
# #c("D324_512", "D320_475")
# 
# sample <- data %>% 
#   arrange(tag_local_identifier, timestamp) %>% 
#   filter(tag_local_identifier == "9551" & burstID %in% c(2200:2500)) %>% 
#   st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) %>% 
#   select(!c("accelerationTransformed", "eobs_accelerations_raw")) %>% 
#   mutate(row_id = row_number())
# 
# mapview(sample, zcol = "flightClust_smooth3")
# 
# mapview(sample, zcol = "flight_status")
# 
# #### prepare to send to Pritish (update Jul.13.23)
# animal_9551_classifiedAcc <- readRDS("~/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/animal_9551_classifiedAcc.rds")
# 
# data_smpl <- data %>% 
#   rename(GPS_flight_class = flightClust_smooth3) %>% 
#   dplyr::select(c(c("location_lat", "location_long", "GPS_flight_class", "accelerationTransformed"),intersect(names(animal_9551_classifiedAcc), names(data))))
# 
# write.csv(data_smpl, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/matched_gps_acc/two_inds_gps&acc_behavClass.csv")

