#This code explores the laterality of honey buzzards during circling flight. This script focuses on exploring these using the quaternion-derived features
#Elham Nourani PhD.
#Apr 24. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")
wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")



#open quat data-in original resolution: prepared in 01b_imu_processing.r
#quat_gps_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_quat",
#                             pattern = "rds", full.names = T)

#open quat data- aggregated for every second, matched with gps-informed flight classification: prepared in 01b_imu_processing.r
or_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_Jul24.rds")


#-------------------------------------------------------------------------------------
# STEP1: what range of angles represents circling flight?
#-------------------------------------------------------------------------------------
#use overlapping GPS and Quat data to calculate this. requires matching gps and quat data (do it in 01b_imu_processing.r)

#look at a sample individual from 2023
smpl <- or_summaries_w_gps %>% 
  filter(individual_local_identifier == "D329_013" &
           imu_burst_id %in% c(2672, 2674, 2676)
  )


smpl_burst_sf <- smpl %>% 
  drop_na("location_lat_closest_gps") %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)

mapview(smpl_burst_sf , zcol = "flight_type_sm3")
mapview(smpl_burst_sf , zcol = "cumulative_yaw", alpha = 0)


one_burst <- smpl_burst_sf %>% 
  filter(imu_burst_id == 2672)

#2023 data for further investigation:
data_2023 <- or_summaries_w_gps %>% filter(year(timestamp) == 2023)

#summary plots for cumulative yaw: only simultaneous GPS+IMU
X11()
ggplot(data_2023 %>% drop_na(flight_type_sm3), aes(x = abs(cumulative_yaw), y = flight_type_sm3, color = flight_type_sm3, fill = flight_type_sm3)) + 
  stat_density_ridges(jittered_points = T, rel_min_height = .01,
                      point_shape = "|", point_size = 2, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5, show.legend = F) +
  labs(y = "", x = "absolute cumulative yaw per second") +
  coord_cartesian(xlim = c(-180,180)) +
  theme_minimal() 

#try histograms for more clarity on 0s
X11()
ggplot(or_summaries_w_gps %>% drop_na(flight_type_sm3) %>% filter(year(timestamp) == 2023),
       aes(x = cumulative_yaw, 
           y = flight_type_sm3, height = stat(density))) + 
  geom_density_ridges(
    stat = "binline", bins = 100, scale = 0.98, alpha = 0.3,
    draw_baseline = FALSE
  ) +
  scale_x_continuous(limits = c(-180, 180)) +
  labs(y = "", x = "Cumulative yaw per second") +
  theme_minimal()

X11()
ggplot(or_summaries_w_gps %>% drop_na(flight_type_sm3), aes(x = pitch_mean, y = flight_type_sm3, color = flight_type_sm3, fill = flight_type_sm3)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, size = 1.5,
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5, show.legend = F) +
  labs(y = "", x = "Mean pitch per second") +
  coord_cartesian(xlim = c(-90, 90)) +
  theme_minimal() 

#there are a lot of 0s in the distributions. explore the gps flight classification when cumulative_yaw is 0
zeros <- or_summaries_w_gps %>% 
  filter(cumulative_yaw == 0 &
           flight_type_sm3 %in% c("circular_soaring", "shallow_circular_soaring"))


zeros_sf <- zeros %>% 
  drop_na("location_lat_closest_gps") %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)

mapview(zeros_sf)

#-------------------------------------------------------------------------------------
# STEP2: during circling, what is bank angle like?
#-------------------------------------------------------------------------------------
#use roll in original resolution (non-aggregated) for all individuals within circling bouts OR take the mode of roll



#-------------------------------------------------------------------------------------
# STEP3: Calculate degree of handedness
#-------------------------------------------------------------------------------------