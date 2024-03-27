#This code filters honey buzzard tracks for sea-crossing sections. then segment flight types.
#follows Martina's code for segmentation: https://github.com/kamransafi/GoldenEagles/blob/main/WP3_Soaring_Ontogeny/MS1_soaring_skills/script1_GPSsegmentation_goldenEagles_newMarch2023.R
# I have renamed turn_angle_smooth (from martina's code) to turn_angle_cum, because it is the cumulative sum and not simply a smoothed value.
#Elham Nourani PhD.
#Aug 2. 2023. Canberra, AU

library(tidyverse)
library(lubridate)
library(sf)
library(mapview)
library(terra)
library(lwgeom)
library(rgl)
library(viridis)

#options(rgl.printRglwidget = TRUE)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

# STEP 1: open all GPS data and filter for latitudes #####
gps <- readRDS("data/all_gps_nov_6_23.rds") %>%
  as.data.frame() %>%
  filter(between(location_lat,53.5,61.05) | between(location_lat, 28.9,46.49))

#open coastline layer
coastlines <- st_read("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>%
  st_crop(xmin = -17, xmax = 43, ymin = -35.6, ymax = 67) %>%
  st_union()
  #poly2nb(st_make_valid(shp))

#include a negative buffer, to retain the last few points before sea-crossing
#coastlines_buffer <- coastlines %>% 
#  st_buffer(dist = -0.1)

gps_sf <- gps %>%
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) 

#extract gps points that fall over the sea and the buffer zone
#over_sea_w_buffer <- gps_sf  %>%
#  st_difference(coastlines_buffer)

#extract points that only fall over the sea
over_sea <- gps_sf  %>%
  st_difference(coastlines) %>% 
  mutate(location_long = st_coordinates(.)[,1],
         location_lat = st_coordinates(.)[,2])

saveRDS(over_sea, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/HB_sea_tracks.rds")

# STEP 2: calculate flight altitude #####

geo <- rast("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/EGM96_us_nga_egm96_15.tif")

sea <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/HB_sea_tracks.rds") %>% 
  #st_drop_geometry() %>% 
  select(individual_local_identifier, location_lat, location_long, timestamp, eobs_battery_voltage, eobs_horizontal_accuracy_estimate, 
         ground_speed, heading, height_above_ellipsoid, tag_local_identifier, individual_taxon_canonical_name) %>% 
  #st_transform(crs = crs(geo)) %>% 
  extract(x = geo, y = ., method = "simple", bind = T) %>% 
  st_as_sf() %>% 
  mutate(height_msl = height_above_ellipsoid - geoid_undulation)

#look at one sample
one <- sea %>% 
  filter(local_identifier == "D324_512")

# STEP 2: associate with flight segmentation (GPS) #####



# STEP 3: associate with IMU metrics #####

flight_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_acc/Flight metrics/",
                           patter = ".csv", full.names = T)



# STEP 4: prepare for annotation #####



#there are NAs in the heigh columns. create a row id column. then only use the complete columns to annotate. row id will help with merging afterwards
sea_df <- sea %>% 
  as.data.frame() %>% 
  mutate(row_id = row_number())

saveRDS(sea_df, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_tracks_df.rds")
  

ann_df <- sea_df %>% 
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
  select(location_lat, location_long, row_id, timestamp)

colnames(ann_df)[c(2,1)] <- c("location-long","location-lat") 

write.csv(ann_df, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_to_annotate.csv", row.names = F)

# plotting code to be used while developing the code ------------------------------------------------------------------------
#2d interactive plotting
mapview(b, zcol = "flight_clust_sm3")

#look at flight altitude
ggplot() +
  geom_path(data = b, aes(x = st_coordinates(b)[,1], y = height_msl), color = "gray") +
  geom_point(data = b, aes(x = st_coordinates(b)[,1], y = height_msl, color = factor(flight_clust_sm3)),
             size = 2.5, shape = 16) +
  scale_color_viridis(discrete=TRUE, alpha = 0.4) +
  theme_linedraw()
