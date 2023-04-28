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

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

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

#world <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
#  st_crop(xmin = -15, xmax = 42, ymin = -35.6, ymax = 66) %>%
#  st_union()


#download natural earth layer from here: https://www.naturalearthdata.com/downloads/10m-raster-data/10m-cross-blend-hypso/
#code below follows this tutorial loosely: https://downwithtime.wordpress.com/2013/12/04/naturalearthdata-and-r-in-ggplot2/
#it uses the raster package. eventually rewrite for terra

background <- stack("/home/enourani/ownCloud/Work/GIS_files/natural_earth/HYP_HR_SR_OB_DR.tif") %>% 
  crop(extent(-16, 43, -35.6, 67)) # crop to the extend of interest

#convert into a dataframe
rast.table <- data.frame(xyFromCell(background, 1:ncell(background)),
                         getValues(background/255))
#add the rgb values to the table
rast.table$rgb <- with(rast.table, rgb(HYP_HR_SR_OB_DR_1,
                                       HYP_HR_SR_OB_DR_2,
                                       HYP_HR_SR_OB_DR_3,
                                       1))

png("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/figs/wind_field_stills_Dec22.png"), 
    height = 5, width = 7, units = "in", res = 300)
print(final_plot)
dev.off()

#plot
ggplot() +
  geom_tile(data = rast.table, aes(x = x, y = y, fill = rgb)) +
#geom_sf(data = world, fill = "grey20", col = "grey20") +
geom_path(data = data, aes(x = location_long, y = location_lat, group = tag_local_identifier), linewidth = .5, lineend = "round", color = "black") +
  theme_void()
  
