#make a map for tweeting and Carla's instagram post
#Elham Nourani, PhD. Oct. 24, 2022
#Konstanz. DE

library(tidyverse)
library(move)
library(moveVis)
library(lubridate)
library(sf)
library(terra)
library(rgdal)
library(mapview)
library(plotKML)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

#### make a gif #####
#download data from Movebank
creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())

EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

data <- getMovebankData(EHB_FN_id, timestamp_start = "20220901000000000", 
                        removeDuplicatedTimestamps = T, login = creds)

#save as a dataframe to send to Patrik
data_df <- data %>% 
  as.data.frame() %>% 
  dplyr::select(c("timestamp", "height_above_ellipsoid", "location_lat", "location_long", "sensor_type", "tag_local_identifier", "trackId", "sex", "ring_id", "taxon_canonical_name"))

write.csv(data_df, "/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/collabs/Patrik_wintering_range_23/data_to_send/HB_gps_Jan23.csv")

#hourly subset using moveVis
mv_hr <- align_move(data, res = 180, unit = "mins")

#assign colors to dead birds
mv_hr$colour <- ifelse(mv_hr$tag_local_identifier %in% c("9557", "9558"), "firebrick1", "cornflowerblue")

save(mv_hr, file = "/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/HB_hrly_for_map.RData")

mv_hr <- mv_hr[, -49]

load("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/HB_hrly_for_map.RData")

ext <- extent(-12.48047, 39.72656, -35.6, 65.29347)

#create frames
frames <- frames_spatial(mv_hr,  map_service = "mapbox", map_type = "satellite", path_legend = F, ext = ext, tail_length = 22,
                         path_colours = c(rep("cornflowerblue", 8), "firebrick1", rep("cornflowerblue", 15),  "firebrick1"),
                         map_token = "pk.eyJ1IjoibWFobGU2OCIsImEiOiJjbDltZmQ3dW8wMDAwdHJudzR2aTB4cXFwIn0.3-y9DoHl8I4j0-6QJZAQ4w") %>% 
  add_labels(x = "", y = "") %>% 
  add_timestamps(type = "label")

#plot one frame to see
frames[[298]]

animate_frames(frames, out_file = "/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/moveVis/HB_map_jan23.mp4", 
               height = 700, width = 500, res = 100, fps = 70)


#### make a still map. Apr 28.2022 #####
#data is on the external

data <- readRDS("/media/enourani/Ellham's HDD/Elham_EHB/all_GPS_Apr4.rds") %>% 
  as.data.frame() %>% 
  group_by(tag_local_identifier, yday(timestamp), hour(timestamp)) %>% #subset to hourly
  slice(1) %>% 
  ungroup() %>% 
  arrange(tag_local_identifier, timestamp)

saveRDS(data, "/media/enourani/Ellham's HDD/Elham_EHB/hourly_gps.rds")

#open the continent boundaries layer
world <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
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


#plot
png("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/figs/migration_map.png", res = 300, units = "in", height = 10, width = 8)
ggplot() +
  geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
  geom_sf(data = world, fill = "white", col = "white") +
  geom_path(data = data, aes(x = location_long, y = location_lat, group = tag_local_identifier), linewidth = .5, lineend = "round", color = "#df4035") +
  theme_void()
dev.off()


#################### repeat the above map, but with coordinates. For BLS8 ################

data <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_nov_6_23.rds") %>% 
  group_by(individual_local_identifier, yday(timestamp), hour(timestamp)) %>% #subset to hourly
  slice(1) %>% 
  ungroup() %>% 
  arrange(individual_local_identifier, timestamp)

ggplot() +
  geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
  geom_sf(data = world, fill = "white", col = "white") +
  geom_path(data = data, aes(x = location_long, y = location_lat, group = tag_local_identifier), linewidth = .5, lineend = "round", color = "#df4035") +
  theme_void()
