#This code is for exploring the flight metrics calculated using the IMU
#Elham Nourani PhD.
#Nov 27. 2023. Konstanz, DE.

library(move2)
library(tidyverse)
library(lubridate)
library(mapview)
library(parallel)
library(sf)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")


#STEP 1: open flight metrics -------------------------------------------------
#these were prepared by Pritish.
flight_files <- list.files("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_acc/Flight metrics/",
                           patter = ".csv", full.names = T)

one_ind <- read.csv(flight_files[[1]]) %>% 
  mutate(dmyr = dmy(str_sub(t_quat, 1,12)),  #convert the day-month-year to POSIXct
         timestamp = as.POSIXct(paste0(dmyr, " ", str_sub(t_quat,13,21)), tz= "UTC")) 

one_ind_sf <- one_ind %>% 
  drop_na("location_long_closest_gps") %>%
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)


#look at a couple of individuals at once
sample <- lapply(flight_files[1:4], read.csv) %>% 
  bind_rows() %>% 
  mutate(dmyr = dmy(str_sub(t_quat, 1,12)),  #convert the day-month-year to POSIXct
         timestamp = as.POSIXct(paste0(dmyr, " ", str_sub(t_quat,13,21)), tz= "UTC"),
         soar_tf = if_else(between(netHeadChange, -45, 45), 0, 1))
         #soar_tf = ifelse(abs(netHeadChange) >= 90, "1", "0")) #1:soaring; 0:not soaring

sample_sf <- sample %>% 
  drop_na("location_long_closest_gps") %>%
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)

mapview(sample_sf, zcol = "soar_tf", alpha = 0)
mapview(sample_sf, zcol = "propFlap", alpha = 0) 

mapview(sample_sf, zcol = "maxTilt", alpha = 0) 
