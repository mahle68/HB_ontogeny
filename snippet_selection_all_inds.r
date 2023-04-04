#This script allows for selection of relevant IMU snippets to use for exploration of IMU
#in collab with Pritish and Ellen
#March 29. 2023
#Elham Nourani, PhD. Konstanz, DE.

library(tidyverse)
library(lubridate)
library(sf)
library(mapview)


wgs <- "+proj=longlat +datum=WGS84 +no_defs"

setwd("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")

#everything is done on two individuals for now: c("D324_512", "D320_475")

## STEP 1: open data  ----------------------------------------------------------------

#segmented GPS data from 02_GPS_segmentation
load("GPSsegmentation_Mar23/classifiedData/animal_D320_475_classifiedBursts_df.rdata") #HRdf_smooth
ind1 <- HRdf_smooth
load("GPSsegmentation_Mar23/classifiedData/animal_D324_512_classifiedBursts_df.rdata")

gps <- ind1 %>% 
  bind_rows(HRdf_smooth)

gps_sf <- gps %>% 
  filter(location_lat <= 60) %>%  #filter out breeding season movements
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs)
  
#open acc data (from 01a_ACC_prep.R)
acc_g <- readRDS("Pritish_collab_IMU/animal_8986_classifiedAcc.rds")

#download mag and quat
creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

mag <- getMovebankNonLocationData(EHB_FN_id, sensorID = 77740402, animalName = "D320_475" , login = creds, removeDuplicatedTimestamps = T)

quat <- getMovebankNonLocationData(EHB_FN_id, sensorID = 819073350, animalName = "D320_475", login = creds, removeDuplicatedTimestamps = T)


## STEP 2: select snippets ----------------------------------------------------------------

#select snippets based on visualizations, either using mapview or Firetail

snippet_ls <- list(
  circling1 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-09-28 11:33:47", tz = "UTC"), as.POSIXct("2022-09-28 11:40:29", tz = "UTC"))) %>%
    mutate(snippet_id = "circling1"),
  
  circling2 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-09-28 11:26:11", tz = "UTC"), as.POSIXct("2022-09-28 11:28:06", tz = "UTC"))) %>%
    mutate(snippet_id = "circling2"),
  circling3 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-10-06 06:53:55", tz = "UTC"), as.POSIXct("2022-10-06 06:55:12", tz = "UTC"))) %>%
    mutate(snippet_id = "circling3"),
  circling4 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-10-06 07:02:06", tz = "UTC"), as.POSIXct("2022-10-06 07:05:58", tz = "UTC"))) %>%
    mutate(snippet_id = "circling4"),
  circling5 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-10-06 07:46:28", tz = "UTC"), as.POSIXct("2022-10-06 07:52:01", tz = "UTC"))) %>%
    mutate(snippet_id = "circling5"),
  circling6 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-10-07 12:42:17", tz = "UTC"), as.POSIXct("2022-10-07 12:45:51", tz = "UTC"))) %>%
    mutate(snippet_id = "circling6"),
  gliding1 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-09-25 11:43:16", tz = "UTC"), as.POSIXct("2022-09-25 11:45:18", tz = "UTC"))) %>%
    mutate(snippet_id = "gliding1"),
  gliding2 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-09-25 11:17:29", tz = "UTC"), as.POSIXct("2022-09-25 11:18:53", tz = "UTC"))) %>%
    mutate(snippet_id = "gliding2"),
  gliding3 = gps_sf %>% 
    filter(tag_local_identifier == 8986 & 
             between(timestamp, as.POSIXct("2022-09-25 11:08:55", tz = "UTC"), as.POSIXct("2022-09-25 11:09:55", tz = "UTC"))) %>%
    mutate(snippet_id = "gliding3")
)

#mapview(gliding1, zcol = "flightClust_smooth3")

## STEP 3: Match the acc, mag and quat ----------------------------------------------------------------


lapply(snippet_ls, function(snippet){
  
  #create a directory for the snippet
  
  dir.create(paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Mar31/tag8986_", 
                     unique(snippet$snippet_id)))
  
  m <- mag %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$timestamp), max(snippet$ timestamp))) %>% 
    dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "magnetic_fields_raw"))
  
  write.csv(m, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Mar31/tag8986_", 
                   unique(snippet$snippet_id), "/mag.csv"))
  a <- acc_g %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$timestamp), max(snippet$ timestamp))) %>% 
    dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "flight_status", "odbaAvg", "eobs_acceleration_g"))
  
  write.csv(a, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Mar31/tag8986_", 
                   unique(snippet$snippet_id), "/acc.csv"))
  
  q <- quat %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$timestamp), max(snippet$ timestamp))) %>% 
    dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "orientation_quaternions_raw"))
  
  write.csv(q, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Mar31/tag8986_", 
                   unique(snippet$snippet_id), "/quat.csv"))
  
})



