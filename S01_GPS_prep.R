#This code filters honey buzzard tracks for sea-crossing sections. then segment flight types.
#follows Martina's code for segmentation: https://github.com/kamransafi/GoldenEagles/blob/main/WP3_Soaring_Ontogeny/MS1_soaring_skills/script1_GPSsegmentation_goldenEagles_newMarch2023.R
#Elham Nourani PhD.
#Aug 2. 2023. Canberra, AU

library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(mapview)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
#open all GPS data and filter for latitudes
gps <- readRDS("/media/enourani/Ellham's HDD/Elham_EHB/all_GPS_Apr4.rds") %>%
  as.data.frame() %>%
  filter(between(location_lat,53.5,61.05) | between(location_lat, 28.9,46.49))

#open coastline layer
coastlines <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/World_Continents.shp") %>%
  st_crop(xmin = -17, xmax = 43, ymin = -35.6, ymax = 67) %>%
  st_union()
  #poly2nb(st_make_valid(shp))

gps_sf <- gps %>%
  st_as_sf(coords = c("location_long.1", "location_lat.1"), crs = wgs) 

dd <- gps_sf  %>%
  st_difference(coastlines)

saveRDS(dd, file = "HB_sea_tracks.rds")
