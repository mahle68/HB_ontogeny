#preparing data for the request by Soeren Greule (June 2024)
#requirements: 
#1-min resolution. 50 km buffer around breeding site. both Swiss and Finnish birds

#14.6.2024 Konstanz, DE.
#Elham Nourani


library(tidyverse)
library(move2)
library(sf)
library(mapview)
library(lubridate)


wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

#open Finnish data and subset for period of interest

finn_1min <- readRDS("data/all_gps_apr15_24.rds") %>% 
  drop_na("location_long") %>% 
  filter(location_lat > 59.5) %>% 
  #subsample for every minute
  mutate(dt_1min = round_date(timestamp, "1 minute")) %>% 
  group_by(individual_local_identifier, dt_1min) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(c(individual_local_identifier, individual_taxon_canonical_name, eobs_horizontal_accuracy_estimate, eobs_speed_accuracy_estimate, eobs_temperature,
         gps_satellite_count, ground_speed, heading, height_above_ellipsoid,  location_lat, location_long, timestamp)) %>% 
  as.data.frame()

fnn_sf <- finn_1min %>%  
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs)



#download data for swiss birds
creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
HB_id <- movebank_get_study_id("European Honey Buzzard Switzerland")


sw_1min <- movebank_retrieve(study_id = 2235142320, sensor_type_id = "gps",   #download data for all individuals 
                         entity_type = "event",  attributes = "all",
                         timestamp_end = as.POSIXct("2023-10-15 00:00:00")) %>% 
  drop_na("location_long") %>% 
  filter(between(location_lat, 45.7, 48.5)) %>% 
  mutate(dt_1min = round_date(timestamp, "1 minute")) %>% 
  group_by(individual_local_identifier, dt_1min) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(c(individual_local_identifier, individual_taxon_canonical_name, eobs_horizontal_accuracy_estimate, eobs_speed_accuracy_estimate, eobs_temperature,
           gps_satellite_count, ground_speed, heading, height_above_ellipsoid,  location_lat, location_long, timestamp)) %>% 
  as.data.frame()


sw_sf <- sw_1min %>%  
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs)



#folder to write the file:
data <- finn_1min %>% 
  bind_rows(sw_1min) %>% 
  as.data.frame()

write.csv(data, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/collabs/Soeren_greule_2024/honey_buzzards_ch_fn.csv",
          row.names = F)
