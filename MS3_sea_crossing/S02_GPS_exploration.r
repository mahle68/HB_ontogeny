#Explore flight repertoir over the open sea
#follows from S01_GPS_prep.r where I segment the data
#Elham Nourani PhD.
#Aug 24. 2023. Konstanz, DE


library(tidyverse)
library(sf)
library(mgcv)
library(mapview)
library(ggridges)
library(move2)
library(units)

source("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/functions.R")

# STEP 0: Open segmented data -------------------------------------------------------------------------------------------------------------------------------------------------------
sea_sf <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/
                     HB_ontogeny_eobs/git_repository/R_files/GPS_seg_Aug23/classified_data", 
                       full.names = T) %>% 
  map(readRDS) %>% 
  bind_rows() %>% 
  mutate(row_id = row_number())

# STEP 1: Env-Data annotation prep ---------------------------------------------------------------------------------------------------------------------------------------------------
sea_df <- sea_sf %>% 
  mutate('location-lat' = st_coordinates(.)[,2],
         'location-long' = st_coordinates(.)[,1],
         timestamp = paste0(as.character(timestamp), ".000")) %>% 
  select('location-long', 'location-lat', "timestamp", "row_id")

write.csv(sea_df, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_to_annotate.csv", row.names = F)

#open annotated files
p950_level <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/sea_annotations/HB_sea_to_annotate.csv-8533425805309564468/HB_sea_to_annotate.csv-8533425805309564468.csv") %>% 
  rename(geopotential = ECMWF.ERA5.PL.Geopotential,
         u_950 = ECMWF.ERA5.PL.U.Wind,
         v_950 = ECMWF.ERA5.PL.V.Wind)  %>% 
  select(- c("location.long", "location.lat", "timestamp"))

surface_level <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/sea_annotations/HB_sea_to_annotate.csv-1397362286642345513/HB_sea_to_annotate.csv-1397362286642345513.csv") %>% 
  rename(blh = ECMWF.ERA5.SL.Boundary.Layer.Height,
         sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u_10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.,
         v_10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.) %>% 
  mutate(delta_t = sst - t2m) %>% 
  select(- c("location.long", "location.lat", "timestamp"))

sea_ann <- sea_sf %>% 
  full_join(p950_level, by = "row_id") %>% 
  full_join(surface_level, by = "row_id") %>% 
  mutate(wind_speed_950 = sqrt(u_950^2 + v_950^2),
         wind_speed_10m = sqrt(u_10m^2 + v_10m^2))
  
saveRDS(sea_ann, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/sea_gps_seg_ann.rds")


# STEP 3: Data exploration ---------------------------------------------------------------------------------------------------------------------------------------------------

sea_sf <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds")
  
  ggplot(sea_sf, aes(x = local_identifier, fill = flight_clust_sm3)) +
  geom_bar()

ggplot(sea_ann, aes(x = flight_clust_sm3, y = wind_speed_950)) +
  geom_boxplot() +
  theme_linedraw()

ggplot(sea_ann, aes(x = flight_clust_sm3, y = delta_t)) +
  geom_boxplot() +
  theme_linedraw()

ggplot(sea_ann, aes(x = flight_clust_sm3, y = wind_speed_10m)) +
  geom_boxplot() +
  theme_linedraw()

#proportion of each flight type
sea_sf %>% 
  group_by(flight_clust_sm3) %>% 
  reframe(ratio = n()/nrow(sea_sf))


# STEP 4: modeling soaring OR flight type as a function of env variables  -----------------------------------------------------------------------------------------------------------------

sea_df <- sea_sf %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) %>% 
  st_drop_geometry() %>% 
  na.omit(flight_clust_sm3)

#multinomial gam
m <- gam(flight_clust_sm3 ~ scale(wind_speed_950) + scale(delta_t) +
           s(location_lat, location_long),
         family = multinom(K = 5), data = sea_df) #k = 1+unique(sea_df$flight_clust_sm3)

# STEP 5: Distribution plots of delta t and wind support (for BLS8) ---------------------------------------------------------------------------------------------------------------------------------------------------

sea_sf <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds") %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) %>% 
#  st_drop_geometry() %>% 
  mutate(lat_region = factor(ifelse(location_lat > 50, "Baltic", "Mediterranean"), levels = c("Mediterranean","Baltic")))


#calculate wind support
mv <- sea_sf %>% 
  mt_as_move2(time_column = "timestamp",
              track_id_column = "local_identifier") %>% 
  mutate(heading = set_units(mt_azimuth(.), "degrees") %>% as.numeric()) %>% #change units to deg
  drop_na(heading) %>% 
  mutate(ws = wind_support(u = u_10m, v = v_10m, heading = heading),
         cw = cross_wind(u = u_10m, v = v_10m, heading = heading),
         wind_speed = sqrt(u_10m^2 + v_10m^2)) %>% 
  as("sf")

# sea_sf <- sea_sf %>% 
#   as("sf") %>% 
#   mutate(ws = wind_support(u = u, v = v, heading = heading),
#          cw = cross_wind(u = u, v = v, heading = heading),
#          wind_speed = sqrt(u^2 + v^2)) %>% 
#   mutate(airspeed = sqrt((gr_speed - ws)^2 + (cw)^2),)


#viz for sanity check:
sea_sf_hr <-  sea_sf %>% 
  group_by(local_identifier, hour(timestamp)) %>% 
  slice(1)

custom_colors <- c("darkviolet", "dodgerblue2")

#### delta t plot

X11(width = 8, height = 7)
delta_t <- ggplot(sea_sf, aes(x = delta_t, y = lat_region, color = lat_region, fill = lat_region)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5, bandwidth = 0.2) +
  labs(y = "", x = expression(Delta*" T (°C)")) +
  scale_color_manual(values = custom_colors, guide = FALSE) +
  scale_fill_manual(values = custom_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,3),
                  xlim = c(-6,10.5)) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 15)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20))) 

ggsave(plot = delta_t, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/delta_t.png", 
       height = 7, width = 8, dpi = 300)


#### wind support plot

X11(width = 8, height = 7)
wspt <- ggplot(mv, aes(x = ws, y = lat_region, color = lat_region, fill = lat_region)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5, bandwidth = 0.2) +
  labs(y = "", x = "Wind support (m/s)") +
  scale_color_manual(values = custom_colors, guide = FALSE) +
  scale_fill_manual(values = custom_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,3),
                  xlim = c(-10, 15)) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 15)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20))) 


ggsave(plot = wspt, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/wspt.png", 
       height = 7, width = 8, dpi = 300)


