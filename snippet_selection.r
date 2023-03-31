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

#open the segmented GPS data from 02_GPS_segmentation

load("GPSsegmentation_Mar23/classifiedData/animal_D320_475_classifiedBursts_df.rdata") #HRdf_smooth
ind1 <- HRdf_smooth
load("GPSsegmentation_Mar23/classifiedData/animal_D324_512_classifiedBursts_df.rdata")

gps <- ind1 %>% 
  bind_rows(HRdf_smooth)

gps_sf <- gps %>% 
  filter(location_lat <= 60) %>%  #filter out breeding season movements
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs)
  
#select snippets based on visualizations, either using mapview or Firetail

circling1 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-09-28 11:33:47", tz = "UTC"), as.POSIXct("2022-09-28 11:40:29", tz = "UTC")))

circling2 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-09-28 11:26:11", tz = "UTC"), as.POSIXct("2022-09-28 11:28:06", tz = "UTC")))

circling3 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-10-06 06:53:55", tz = "UTC"), as.POSIXct("2022-10-06 06:55:12", tz = "UTC")))

circling4 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-10-06 07:02:06", tz = "UTC"), as.POSIXct("2022-10-06 07:05:58", tz = "UTC")))

circling5 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-10-06 07:46:28", tz = "UTC"), as.POSIXct("2022-10-06 07:52:01", tz = "UTC")))

circling6 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-10-07 12:42:17", tz = "UTC"), as.POSIXct("2022-10-07 12:45:51", tz = "UTC")))




mapview(gliding1, zcol = "flightClust_smooth3")


gliding1 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-09-25 11:43:16", tz = "UTC"), as.POSIXct("2022-09-25 11:45:18", tz = "UTC")))


gliding2 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-09-25 11:17:29", tz = "UTC"), as.POSIXct("2022-09-25 11:18:53", tz = "UTC")))

gliding3 <- gps_sf %>% 
  filter(tag_local_identifier == 8986 & 
           between(timestamp, as.POSIXct("2022-09-25 11:08:55", tz = "UTC"), as.POSIXct("2022-09-25 11:09:55", tz = "UTC")))

