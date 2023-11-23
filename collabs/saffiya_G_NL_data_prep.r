#This code downloads honey buzzard GPS data and filters it for tracks over Africa, for collaboration with Saffiya Ginel et al, Wageningen
#Elham Nourani PhD.
#Nov. 23. 2023. Konstanz, DE.

library(move2)
library(tidyverse)
library(lubridate)
library(mapview)
library(sf)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#STEP 1: download all GPS data -------------------------------------------------

creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
HB_id <- movebank_get_study_id("European Honey Buzzard_Finland")


gps <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "gps", 
                         entity_type = "event",  attributes = "all")


#STEP 2: filter for Africa -------------------------------------------------

Africa <- st_read("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
  filter(CONTINENT == "Africa") %>% 
  st_union()

gps_a <- gps %>% 
  drop_na("location_lat") %>% 
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) %>% 
  st_intersection(Africa)

gps_df <- gps_a %>% 
  mutate(location_long = st_coordinates(.)[,1],
         location_lat = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  select("location_long", "location_lat", "timestamp", "individual_local_identifier", "tag_local_identifier", "individual_taxon_canonical_name","eobs_battery_voltage",
         "eobs_horizontal_accuracy_estimate", "eobs_speed_accuracy_estimate", "eobs_temperature", "gps_satellite_count", "ground_speed", "heading",
         "height_above_ellipsoid") %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier),
         tag_local_identifier = as.character(tag_local_identifier)) %>% 
  rename(eobs_battery_voltage_mV = eobs_battery_voltage,
         eobs_horizontal_accuracy_estimate_m = eobs_horizontal_accuracy_estimate,
         eobs_speed_accuracy_estimate_ms = eobs_speed_accuracy_estimate,
         eobs_temperature_C = eobs_temperature,
         ground_speed_ms = ground_speed,
         heading_deg = heading,
         height_above_ellipsoid_m = height_above_ellipsoid) %>% 
  as.data.frame()

write.csv(gps_df, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Collaborations/HB_wageningen/to_share/honey_buzzards_gps_nourani.csv", row.names = F)


#STEP 3: plot to check -------------------------------------------------

ggplot() + geom_sf(data = Africa) +
  geom_sf(data=gps_a) +
  theme_void() 

#STEP 4: prep meta-data -------------------------------------------------

raw <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Collaborations/HB_wageningen/EHB_metadata_raw.csv") %>% 
  filter(ring_ID %in% gps_df$individual_local_identifier) %>% 
  rename(individual_local_identifier = ring_ID,
         tag_local_identifier = eobs_ID) %>%
  mutate(nest_country = "Finland") %>% 
  select(-c(11:25)) #28 individual

write.csv(raw, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Collaborations/HB_wageningen/to_share/honey_buzzards_metadata_nourani.csv", row.names = F)

