#This code explores the laterality of honey buzzards during circling flight. This script focuses on exploring these using the quaternion-derived features
#Elham Nourani PhD.
#Apr 24. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)
library(viridis)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")
wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

birds_2023 <- c("D329_013", "D329_012", "D329_015", "D329_014", "D326_193", "D326_192")

#open quat data-in original resolution: prepared in 01b_imu_processing.r
#quat_gps_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_quat",
#                             pattern = "rds", full.names = T)

#open quat data- aggregated for every second, matched with gps-informed flight classification: prepared in 01b_imu_processing.r
or_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_Jul24.rds") %>% 
  drop_na(individual_local_identifier)


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

#2023 data for further investigation to only have simultaneous GPS+IMU
data_2023 <- or_summaries_w_gps %>% 
  filter(individual_local_identifier %in% birds_2023)

#summary plots for cumulative yaw: try histograms for more clarity on 0s
X11()
ggplot(data_2023 %>% drop_na(flight_type_sm3),
       aes(x = cumulative_yaw, 
           y = flight_type_sm3, height = stat(density))) + 
  geom_density_ridges(
    stat = "binline", bins = 500, scale = 0.98, alpha = 0.3,
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

#distribution of cumulative yaw at each flight category
circling_2023 <- data_2023 %>% 
  filter(flight_type_sm3 %in% c("shallow_circular_soaring", "circular_soaring"))

data_2023 %>% 
  #combine the two types of circling
  mutate(circling_vs_linear = ifelse(flight_type_sm3 %in% c("shallow_circular_soaring", "circular_soaring"), 
                                     "circling", flight_type_sm3)) %>% 
  group_by(circling_vs_linear) %>%
  #group_by(flight_type_sm3) %>%
  summarize(mean = mean(cumulative_yaw),
            min = min(cumulative_yaw),
            max = max(cumulative_yaw),
            sd = sd(cumulative_yaw),
            quant_90 = quantile(cumulative_yaw, probs = 0.90),
            quant_10 = quantile(cumulative_yaw, probs = 0.10),
            abs_quant_10 = quantile(abs(cumulative_yaw), probs = 0.10),
            abs_quant_25 = quantile(abs(cumulative_yaw), probs = 0.25),
            abs_quant_75 = quantile(abs(cumulative_yaw), probs = 0.75),
            abs_quant_90 = quantile(abs(cumulative_yaw), probs = 0.90),
            abs_mean = mean(abs(cumulative_yaw)),
            abs_median = median(abs(cumulative_yaw)),
            abs_min = min(abs(cumulative_yaw))) %>% 
  as.data.frame()
 
#take 21 degrees (abs_median of circling, when both circling types are combined)

data_circling <- data_2023 %>% 
  mutate(circling_ys_no = ifelse(cumulative_yaw > 14 | cumulative_yaw < -14, "circling", "not_circling"))

circling_ys_no_sf <- data_circling %>% 
  arrange(individual_local_identifier, timestamp) %>% 
  drop_na("location_lat_closest_gps") %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)

mapview(circling_ys_no_sf[1:5000,], zcol = "circling_ys_no")

#-------------------------------------------------------------------------------------
# STEP2: during circling, what is bank angle like?
#-------------------------------------------------------------------------------------
#use cumulative roll for all individuals within circling bouts

# #based on summaries for each second
# or_summaries_circling <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_Jul24.rds") %>% 
#   mutate(circling_ys_no = ifelse(cumulative_yaw > 16 | cumulative_yaw < -16, "circling", "not_circling")) %>% 
#   filter(circling_ys_no == "circling")
# 
# ggplot() +
#   geom_point(data = or_summaries_circling, aes(x = timestamp, y = cumulative_roll)) +
#   facet_wrap(~individual_local_identifier)


#based on summaries for each 8-second burst


#or_8sec_summaries <- readRDS("quat_summaries_8secs_Jul24.rds")
#cumulative yaw in the above file seems wrong. for now, calculate the cumulative yaw for the 8 second bursts using the summaries for each second
or_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_Jul24.rds") %>% 
  drop_na(individual_local_identifier)

## add a burst ID for 8-second bursts (same code was used in MS1_IMU_classification/01b_imu_processing.r)
or_8sec_summaries <- or_summaries_w_gps %>% 
  group_by(individual_local_identifier) %>% 
  arrange(timestamp, .by_group = T) %>% 
  #the burst assignments are correct, but burst lengths are different depending on the year, etc. for the sake of uniformity, assign new burst ids that break up the data into 8 seconds
  mutate(timelag = c(0, diff(timestamp)),  #calculate time lag 
         burst_id_8sec = cumsum(timelag > 8),
         unique_burst_id_8sec = paste0(individual_local_identifier, "_", burst_id_8sec)) %>% 
  ungroup() #%>%
  #aggregate summarized angles over 8 seconds
  group_by(unique_burst_id_8sec) %>% #this variable already has individual ID in it. So no need to group by ind_id AND burst_id
  summarize(across(c(study_id, individual_taxon_canonical_name, individual_local_identifier, orientation_quaternions_sampling_frequency,
                     tag_local_identifier, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
            across(c(yaw_sd, roll_sd, pitch_sd,
                     yaw_mean, roll_mean, pitch_mean,
                     yaw_min, roll_min, pitch_min,
                     yaw_max, roll_max, pitch_max), 
                   ~mean(., na.rm = T),
                   .names = "mean_{.col}"), #rename the column, because these values are the mean of the sd, max, min values across the 8 seconds
            across(c(cumulative_yaw, cumulative_roll, cumulative_pitch),
                   ~sum(., na.rm = T),
                   .names = "{.col}_8sec"),
            across(c(timestamp, location_lat_closest_gps, location_long_closest_gps, height_above_ellipsoid_closest_gps), 
                   ~head(.,1),
                   .names = "start_{.col}"),
            across(c(flight_type_sm2, flight_type_sm3, track_flight_seg_id),
                   ~Mode(.), #make sure the Model function is defined
                   .names = "{.col}_Mode"),
            .groups = "keep") %>% 
  ungroup %>% 
  as.data.frame()

saveRDS(or_8sec_summaries, "matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_Jul24.rds") #this doesn't have the NA individual

#assign circling vs non-circling to values of yaw
or_8sec_circling_ys_no <- or_8sec_summaries %>% 
  mutate(circling = ifelse(cumulative_yaw_8sec >= 20 | cumulative_yaw_8sec <= -20, "yes", "no"))

#look at the distribution of circling vs not, on the map

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


X11(width = 13, height = 15)
(circling_p <- ggplot()+
  geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
  geom_sf(data = world, fill = "white", col = "white") +
  coord_sf(xlim = c(-17, 38), ylim = c(2, 65), expand = FALSE) +
  geom_point(data = or_8sec_circling_ys_no %>% arrange(individual_local_identifier, start_timestamp) %>% drop_na(start_location_long_closest_gps), 
            aes(x = start_location_long_closest_gps, 
                y = start_location_lat_closest_gps, 
                col = circling), size = 1.2)
) 
  
  
#check the circling assignments with the 1Hz GPS data for a subset
#use smpl_burst_sf from above

smpl_8sec <- or_8sec_circling_ys_no %>% 
  filter(between(start_timestamp, head(smpl_burst_sf$timestamp, 1), tail(smpl_burst_sf$timestamp, 1))) %>% 
  drop_na(start_location_long_closest_gps) %>% 
 # filter(track_flight_seg_id_Mode %in% smpl_burst_sf$track_flight_seg_id) %>% 
  st_as_sf(coords = c("start_location_long_closest_gps", "start_location_lat_closest_gps"), crs = wgs)

mapview(smpl_burst_sf, alpha = 0, cex = 0.5) + mapview(smpl_8sec, zcol = "circling")

#-------------------------------------------------------------------------------------
# STEP3: Calculate degree of handedness
#-------------------------------------------------------------------------------------

#look at the relationship between yaw and roll. they should be highly correlated
or_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_Jul24.rds")

#scatterplot of cumulative yaw and cumulative roll
ggplot() +
  geom_point(data = or_summaries_w_gps, aes(x = cumulative_yaw, y = cumulative_roll))







