library(tidyverse)
library(sf)
library(mapview)
library(rnaturalearth)

##all data
Cr_gps <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24.rds") %>%  #2022-09-25 07:56:41.0000 to 2024-04-15 10:02:33.0000
  drop_na(individual_local_identifier, location_lat) %>% 
  st_as_sf(coords = c("location_long", "location_lat"), crs = "EPSG:4326") %>% 
  st_crop(xmin = 11, xmax = 21, ymin = 39, ymax = 47)



#reduce res to 15 min
Cr_lres <- Cr_gps %>% 
  #subsample for every minute
   mutate(dt_15min = round_date(timestamp, "15 minutes")) %>% 
   group_by(individual_local_identifier, dt_15min) %>% 
   slice(1) %>% 
   ungroup()

mapview(Cr_lres, zcol = "individual_local_identifier")  


Cr_df<- Cr_lres %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) %>% 
  select(c(individual_local_identifier, individual_taxon_canonical_name, eobs_horizontal_accuracy_estimate, eobs_speed_accuracy_estimate, eobs_temperature,
           gps_satellite_count, ground_speed, heading, height_above_ellipsoid, location_lat, location_long, timestamp)) %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  st_drop_geometry() %>% 
  as.data.frame()

write.csv(Cr_df, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/collabs/Croatia_2025/honey_buzzards_Adriatic.csv",
          row.names = F)

##annotated gps data
# setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")
# 
# the_one <- readRDS("R_files/gps_seg_apr24/animal_D329_015_classified_bursts.rds") %>% 
#   filter(year(timestamp) == 2023) %>% 
#   st_crop(xmin = 0, xmax = 15, ymin = 42, ymax = 50)
# 
# 
# 
# mapview(the_one, zcol = "flight_clust_sm3")
