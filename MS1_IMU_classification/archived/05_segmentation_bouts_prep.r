#This code searches for continous bouts of each fligth behavior from the gps-segmented trajectories of honey buzzards tagged in 2023 to validate the algorithm segments flight types for the European honey buzzards. based on my code: A01_GPS_prep.r
#Elham Nourani PhD.
#Feb 26, 2024. Konstanz, DE

library(tidyverse)
library(sf)
library(mapview)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

# STEP 1: open segmented data (prepared in 01_c_GPS_prep.r) ####
birds_23 <- c("D329_012", "D329_013", "D329_014", "D329_015", "D326_193", "D326_192")

#calculate duration of each behavioral segment
gps_seg <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/all_gps_seg_ann_Nov2023.rds") %>% 
  filter(individual_local_identifier %in% birds_23) %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  group_by(individual_local_identifier, flight_seg_id_sm3, flight_clust_sm3) %>% 
  mutate(duration = sum(time_lag_sec)) %>% 
  filter(n() == 236.0) %>% 
  ungroup()
