

library(tidyverse)
library(lubridate)
library(move2)
library(sf)
library(mapview)
library(terra)


wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")


#STEP 1: download gps data for the past n weeks -----------------------------

creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
HB_id <- movebank_get_study_id("European Honey Buzzard_Finland")

#timerange... past n weeks
n_weeks <- 3

end <- as.character(today() - 1) %>% 
  str_remove_all(pattern = "-") %>% 
  str_c("100000000")

start <- as.character(today() - (n_weeks * 7)) %>% 
  str_remove_all(pattern = "-") %>% 
  str_c("100000000")


gps <- movebank_retrieve(study_id = HB_id, sensor_type_id = "gps", #download data for wintering 
                         entity_type = "event",  attributes = "all",
                         timestamp_start = start,
                         timestamp_end = end)


#plot

gps_sf <- gps %>% 
  #remove the points at (0,0) ... there are 54 of them!!
  filter(!(location_lat == 0 & location_long == 0)) %>%
  drop_na("location_long") %>% 
  st_as_sf(coords = c("location_long", "location_lat"), crs = "EPSG:4326") 

mapview(gps_sf, zcol = "individual_local_identifier")
