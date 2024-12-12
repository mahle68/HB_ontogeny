# script for making the map for Safi et al 2025
# Elham Nouani, PhD. 
# 11.12.2024, Konstanz, DE

library(tidyverse)
library(sf)
library(mapview)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

# STEP 1: map tracking points ----------------------------------

#open metadata to extract deployment date (the tracking data has some undeployed points in it)
deployment <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/EHB_metadata - Sheet1.csv") %>% 
  filter(nest_country == "Finland") %>% 
  mutate(deployment_dt_utc = as.POSIXct(deployment_dt_utc)) %>% 
  dplyr::select(ring_ID, deployment_dt_utc) %>% 
  rename(individual_local_identifier = ring_ID)

#life-cycle stages from L03a_tests_per_day.r AND append deployment
life_cycle <- readRDS("updated_life_cycle_nov24.rds") %>% 
  full_join(deployment, by = "individual_local_identifier")

#rewrite the life cycle file, to include deployment
saveRDS(life_cycle, file = "updated_life_cycle_nov24.rds")


data <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24.rds") 

cleaned_gps <- data %>% 
  #remove the points at (0,0) 
  filter(!(location_lat == 0 & location_long == 0)) %>% 
  mutate(unique_date = as.Date(timestamp)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, deployment_dt_utc, first_exploration, migration_start, migration_end), by = "individual_local_identifier") %>% 
  #remove data before deployment
  filter(timestamp >= deployment_dt_utc) %>% 
  #hourly subset to make the next step go faster!
  group_by(individual_local_identifier, yday(timestamp), hour(timestamp)) %>% #subset to hourly
  slice(1) %>% 
  ungroup() %>% 
  group_by(individual_local_identifier) %>% 
  rowwise() %>% 
  mutate(life_stage = case_when(
    between(unique_date, migration_start, migration_end) ~ "migration",
    between(unique_date, first_exploration, migration_start) ~ "post-fledging",
    unique_date < first_exploration ~ "pre-fledging",
    unique_date > migration_end ~ "wintering",
    TRUE ~ NA_character_
  )) %>% 
  ungroup() %>% 
  arrange(individual_local_identifier, timestamp) %>% 
  as.data.frame()

saveRDS(cleaned_gps, "cleaned_gps_for_laterality_map.rds")


wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

#open the continent boundaries layer
world <- st_read("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
  st_crop(xmin = -17, xmax = 43, ymin = -35.6, ymax = 67) %>%
  st_union()

#create a rectangle to be the oceans
Polycoords <- data.frame(long = c(-17,43),
                         lat = c(-35.6,67))

pol <- st_polygon(
  list(
    cbind(
      Polycoords$lon[c(1,2,2,1,1)], 
      Polycoords$lat[c(1,1,2,2,1)])
  )
) %>% 
  st_sfc(crs = wgs)




ggplot() +
  geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
  geom_sf(data = world, fill = "white", col = "white") +
  geom_path(data = subset(cleaned_gps, life_stage == "migration"), 
            aes(x = location_long, y = location_lat, 
                group = individual_local_identifier), 
            linewidth = .5, lineend = "round", color = "#df4035", linetype = "solid") +
  geom_path(data = subset(cleaned_gps, life_stage == "post-fledging"), 
            aes(x = location_long, y = location_lat, 
                group = individual_local_identifier), 
            linewidth = .5, lineend = "round", color = "#df4035", linetype = "dotted") +
  geom_path(data = subset(cleaned_gps, life_stage == "wintering"), 
            aes(x = location_long, y = location_lat, 
                group = individual_local_identifier), 
            linewidth = .5, lineend = "round", color = "#df4035", linetype = "dotted") +
  theme_void()
