#Explore flight repertoire over the open sea
#follows from S01_GPS_prep.r where I segment the data and S02_GPS_exploration.r where env annotation was done
#this code also includes my plotting attempts for BLS8 Tokyo

#Elham Nourani PhD.
#Feb 16. 2024. Konstanz, DE


library(tidyverse)
library(lubridate)
library(sf)
library(mgcv)
library(mapview)
library(ggridges)
library(viridis)
library(hrbrthemes)
library(rnaturalearth)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

# STEP 1: Open annotated gps data ---------------------------------------------------------------------------------------------
sea_df <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds") %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) %>% 
  st_drop_geometry() %>% 
  as.data.frame()

sea_sf <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds")

sf_ls <- split(sea_sf, sea_sf$local_identifier)

mapview(sf_ls[[1]], zcol = "flight_clust_sm3")


base <- world <- ne_coastline(scale = 'medium', returnclass = 'sf')

#color palettes: mako is blue, magma is red

# ws <- ggplot(data = base) +
#   geom_sf(col = "gray", fill = "gray") +
#   coord_sf(xlim = c(-10, 38), ylim = c(3, 64), expand = FALSE) +
#   geom_point(data = sea_df, aes(x = location_long, y = location_lat, col = flight_clust_sm3)) +
#   scale_colour_viridis(option = "magma", na.value = "white", name = "m/s", alpha = 0.7) +
#   theme_linedraw() +
#   scale_x_continuous(breaks = c(0,30)) +
#   scale_y_continuous(breaks = c(10,30,50)) +
#   theme(axis.text = element_text(size = 10, colour = 1),
#         legend.text = element_text(size = 10, colour = 1),
#         legend.title = element_text(size = 10, colour = 1),
#         legend.position = "right",
#         legend.background = element_rect(colour = NULL, fill = "white"))+
#   labs(x = NULL, y = NULL, title = "Wind support") +
#   facet_wrap(.~year, nrow = 1)


# STEP 2: Open IMU metrics ---------------------------------------------------------------------------------------------

#for each individual, filter the IMU metrics for the time periods of sea-crossing (Baltic and Med seas)
flight_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_acc/Flight metrics",
                           patter = ".csv", full.names = T)

flight_ls <- lapply(sf_ls, function(x){ #for each individual
  
  #assign track ID to points that are max 15 min apart. use this to match with the IMU locations, but also allow for flexibility to include the IMU between GPS bursts
  
  x <- x %>% 
    arrange(timestamp) %>% 
    mutate(time_lag_min = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "mins") %>% as.numeric()),
           sea_track_id = base::cumsum(time_lag_min > 15)) #increase id by one, every time time_lag is more than 15 mins
  
  #open file that contains data for this individual
  flight_for_ind <- flight_files[[str_which(flight_files,unique(x$local_identifier))]] %>% 
    read.csv() %>% 
    mutate(dmyr = dmy(str_sub(t_quat, 1,12)),  #convert the day-month-year to POSIXct
           timestamp = as.POSIXct(paste0(dmyr, " ", str_sub(t_quat, 13, 21)), tz= "UTC")) 
  
  #for each track, extract the IMU data between the start and end time of the track
  lapply(split(x, x$sea_track_id), function(seg){
    
    
    flight_for_seg <- flight_for_ind %>% 
      filter(between(timestamp, min(seg$timestamp), max(seg$timestamp))) %>% 
      mutate(individual_local_identifier = unique(seg$local_identifier),
             sea_track_id = unique(seg$sea_track_id))
    
  }) %>% 
    bind_rows()
  
})

saveRDS(flight_ls, "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_IMU.rds")


# STEP 4: Exploration ---------------------------------------------------------------------------------------------

flight_df <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_IMU.rds") %>% 
  bind_rows()

#proportion of flapping to non-flapping
flight_df %>% 
  group_by(round(propFlap)) %>% 
  reframe(ratio = n()/nrow(flight_df)) #0.57


ggplot(flight_df, aes(x = yday(timestamp), y = stdBank_nonFlap, color = individual_local_identifier)) +
  geom_point() + 
  facet_wrap(~individual_local_identifier)

#density plots of flight metrics for each region
flight_df <- flight_df %>% 
  mutate(lat_region = factor(ifelse(location_lat_closest_gps > 50, "Baltic", "Mediterranean"), levels = c("Mediterranean","Baltic")))

#proportion of circular soaring out of all non-flapping flights
non_flap <- flight_df %>% 
  filter(propFlap < 0.6) %>% 
  mutate(soar_tf = if_else(between(netHeadChange, -45, 45), 0, 1)) #%>% 
#group_by(soar_tf) %>% 
#reframe(ratio = n()/nrow(non_flap)) 

custom_colors <- c("darkviolet", "dodgerblue2")

X11(width = 8, height = 8)
prop_flap <- ggplot(flight_df, aes(x = propFlap, y = lat_region, color = lat_region, fill = lat_region)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5) +
  labs(y = "", x = "Proportion of flapping") +
  scale_color_manual(values = custom_colors, guide = FALSE) +
  scale_fill_manual(values = custom_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,3)) +
  #xlim(0,1) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 15)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20))) 


ggsave(plot = prop_flap, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/prop_flap_sea.jpeg", 
       height = 8, width = 8, dpi = 300)

#--------------
X11(width = 8, height = 8)
numcircl <- ggplot(non_flap, aes(x = netHeadChange, y = lat_region, color = lat_region, fill = lat_region)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5) +
  labs(y = "", x = "N of circles (in 8 seconds)") +
  scale_color_manual(values = custom_colors, guide = FALSE) +
  scale_fill_manual(values = custom_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,3)) +
  #xlim(c(-30, 12)) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 15)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20)))

ggsave(plot = numcircl, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/numcircl.jpeg", 
       height = 8, width = 8, dpi = 300)

## explore the proportion of circling vs other passive soaring
non_flap_soaring <- non_flap %>% 
  filter(meanTilt >= 0) #only upward movement

numcircl <- ggplot(non_flap_soaring, aes(x = numCircles, y = lat_region, color = lat_region, fill = lat_region)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5) +
  labs(y = "", x = "N of circles (in 8 seconds)") +
  scale_color_manual(values = custom_colors, guide = FALSE) +
  scale_fill_manual(values = custom_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,3)) +
  xlim(c(-0.08, 1)) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 15)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20)))

ggsave(plot = numcircl, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/numcircl_when_passive.png", 
       height = 8, width = 9, dpi = 300)



# --------------
#std bank angle
X11(width = 8, height = 7)
stdBank <- ggplot(non_flap, aes(x = stdBank, y = lat_region, color = lat_region, fill = lat_region)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5) +
  labs(y = "", x = "Std bank angle (deg)") +
  scale_color_manual(values = custom_colors, guide = FALSE) +
  scale_fill_manual(values = custom_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,3)) +
  #xlim(c(-360, 360)) +
  #xlim(-0.1, 1.5) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 15)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20)))

ggsave(plot = stdBank, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/stdbank.jpeg", 
       height = 8, width = 9, dpi = 300)


# --------------
#std body tilt
X11(width = 8, height = 7)
stdTilt <- ggplot(non_flap, aes(x = stdTilt, y = lat_region, color = lat_region, fill = lat_region)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5) +
  labs(y = "", x = "Std body tilt (deg)") +
  scale_color_manual(values = custom_colors, guide = FALSE) +
  scale_fill_manual(values = custom_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,3)) +
  #xlim(c(-360, 360)) +
  #xlim(-0.1, 1.5) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 15)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20)))

ggsave(plot = stdTilt, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/stdtilt.jpeg", 
       height = 8, width = 9, dpi = 300)


### MAKE A MAP of proportion of flapping for the whole migration ####################

#open IMU-matched with GPS data... as of Feb 24: Pritish has deleted many of the files, so I only have 7 to work with
all_flight_df <- list.files("/home/mahle68/Desktop/matched_imu_gps_dec17/matched_gps_acc/Flight metrics/",
                            patter = ".csv", full.names = T) %>% 
  lapply(read.csv) %>% 
  bind_rows() %>% 
  mutate(dmyr = dmy(str_sub(t_quat, 1,12)),  #convert the day-month-year to POSIXct
         timestamp = as.POSIXct(paste0(dmyr, " ", str_sub(t_quat, 13, 21)), tz= "UTC")) 


#ss <- all_flight_df %>% 
#  st_as_sf(coords = c("lon", "lat"), crs = "+proj=longlat +datum=WGS84 +no_defs")

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

#the original code is in migration_map.R
data <-readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds")
  

X11(width = 13, height = 15)
(flappiong_map <- ggplot()+
    geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
    geom_sf(data = world, fill = "white", col = "white") +
    coord_sf(xlim = c(-17, 38), ylim = c(2, 65), expand = FALSE) +
    geom_path(data = all_flight_df, aes(x = location_long_closest_gps, y = location_lat_closest_gps, col = propFlap), 
              lwd = 3, lineend = "round") +
    scale_colour_viridis(option = "cividis", na.value = "white", direction = -1,
                         name = "Proportion\nof flapping", alpha = 0.7) +
    theme_void() +
    scale_x_continuous(breaks = c(0,20)) +
    scale_y_continuous(breaks = c(20,40,60)) +
    theme(axis.text = element_text(size = 16, colour = 1),
          legend.text = element_text(size = 16, colour = 1), 
          legend.title = element_text(size = 20, colour = 1, margin = margin(b = 10)),
          legend.position = c(.92,.10),
          legend.key.width = unit(1, "cm"),
          legend.key.height = unit(0.7, "cm"),
          legend.background = element_rect(colour = "white", fill = "white")) + 
          labs(x = NULL, y = NULL, title = "")
)

ggsave(plot = flappiong_map, 
       filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/flapping_map.png", 
       height = 15, width = 13, dpi = 300)

################## zoom in on the circling flight bouts ##############

#med sea: ind D324_510
data <-readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds") %>% 
  filter(local_identifier == "D324_510") %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) 

circling <- ggplot() +
  geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
  #geom_sf(data = world, fill = "white", col = "white") +
  geom_path(data = data, aes(x = location_long, y = location_lat, group = tag_local_identifier), linewidth = .5, 
            lineend = "round", color = "#df4035") +
  xlim(c(20.755, 20.77)) +
  ylim(c(37.13, 37.14)) +
  labs(y = "Latitude", x = "Longitude") +
  theme_classic()

ggsave(plot = circling, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/circling_map_med.png", 
       height = 6.3, width = 7.5, dpi = 300)

#baltic sea D321_349
data_b <-readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds") %>% 
  filter(local_identifier == "D321_349") %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1])

circling_b <- ggplot() +
  geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
  #geom_sf(data = world, fill = "white", col = "white") +
  geom_path(data = data_b, aes(x = location_long, y = location_lat, group = tag_local_identifier), linewidth = .5, 
            lineend = "round", color = "#df4035") +
  xlim(c(23.68, 23.78)) +
  ylim(c(59.4, 59.44)) +
  labs(y = "Latitude", x = "Longitude") +
  theme_classic()

ggsave(plot = circling_b, 
       filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/conferences/BLS8_tokyo_2023/presentation_prep/figs/circling_map_baltic.png", 
       height = 6.3, width = 7.5, dpi = 300)

