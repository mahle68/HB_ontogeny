#plot Pritish's flapping classification based on acc
#Elham Nourani, PhD. Konstanz, Germany
#Aug. 25. 2023


library(tidyverse)
library(sf)
library(mapview)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

#open files for two sample individuals
two_inds <- list.files("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_acc",
                       pattern = "flap.csv", full.names = T) %>% 
  map(read.csv) %>% 
  bind_rows() %>% 
  drop_na(location_long_closest_gps) %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs) %>% 
  mutate(flapping = as.factor(flap_indicator))

#compare to the high-res gps segmentation
gps_seg <- str_subset(list.files("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/GPS_seg_Aug23/classified_data", full.names = T),
           pattern = paste0(unique(two_inds$individual_local_identifier), collapse = '|')) %>%  #add the OR sign in between the two names!
  map(readRDS) %>% 
  bind_rows()


mapview(gps_seg, zcol = "flight_clust_sm3", alpha = 0) + 
  mapview(two_inds %>% filter(flapping == 1), color = "red") +
  mapview(two_inds %>% filter(flapping == 0), color = "black")
