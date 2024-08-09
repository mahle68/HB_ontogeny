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

#source the imu conversion functions
source("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/HB_ontogeny/MS1_IMU_classification/00_imu_diy.r")

birds_2023 <- c("D329_013", "D329_012", "D329_015", "D329_014", "D326_193", "D326_192")

#open quat data-in original resolution: prepared in 01b_imu_processing.r
#quat_gps_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_quat",
#                             pattern = "rds", full.names = T)

#open quat data- aggregated for every second, matched with gps-informed flight classification: prepared in 01b_imu_processing.r
or_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8secIDs_Jul24.rds") %>% #this version also has the IDs assigned to 8-sec bursts.
  drop_na(individual_local_identifier)

#-------------------------------------------------------------------------------------
# STEP1: what range of angles represents circling flight?
#-------------------------------------------------------------------------------------
#use overlapping GPS and Quat data to calculate this. requires matching gps and quat data (do it in 01b_imu_processing.r)

#look at a sample individual from 2023
smpl <- or_summaries_w_gps %>% 
  filter(individual_local_identifier == "D329_013" &
  #         imu_burst_id %in% c(2672, 2674, 2676, 2677, 2678)
    imu_burst_id %in% c(4927:5817))


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
# STEP2: during circling, what is bank angle like? (using 8-sec summaries)
#-------------------------------------------------------------------------------------
#use cumulative roll for all individuals within circling bouts

#open summaries for 8-second bursts
or_8sec_summaries <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_Jul24.rds")


#assign circling vs non-circling to values of yaw
or_8sec_circling_ys_no <- or_8sec_summaries %>% 
  mutate(circling = ifelse(cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45, "yes", "no"))

#look at the distribution of circling vs not, on the map

world <- st_read("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
  st_crop(xmin = -17, xmax = 43, ymin = -35.6, ymax = 67) %>%
  st_union()

#create a rectangle to be the oceans
Polycoords <- data.frame(long = c(-17,43),
                         lat = c(-50,67))

pol <- st_polygon(
  list(
    cbind(
      Polycoords$lon[c(1,2,2,1,1)], 
      Polycoords$lat[c(1,1,2,2,1)])
  )
) %>% 
  st_sfc(crs = wgs)


X11(width = 10, height = 15)
(circling_p <- ggplot()+
  geom_sf(data = pol, fill = "powderblue", col = "powderblue") +
  geom_sf(data = world, fill = "white", col = "white") +
  coord_sf(xlim = c(-17, 38), ylim = c(2, 65), expand = FALSE) +
  geom_point(data = or_8sec_circling_ys_no %>% arrange(individual_local_identifier, start_timestamp) %>% drop_na(start_location_long_closest_gps), 
            aes(x = start_location_long_closest_gps, 
                y = start_location_lat_closest_gps, 
                col = circling), size = 1.2)
)
  
### check the circling assignments with the 1Hz GPS data for a subset 
#use smpl_burst_sf from above

smpl_8sec <- or_8sec_circling_ys_no %>% 
  filter(individual_local_identifier == smpl_burst_sf$individual_local_identifier &
         imu_burst_id %in% smpl_burst_sf$imu_burst_id) %>% 
  #filter(between(start_timestamp, head(smpl_burst_sf$timestamp, 1), tail(smpl_burst_sf$timestamp, 1))) %>% 
  drop_na(start_location_long_closest_gps) %>% 
  st_as_sf(coords = c("start_location_long_closest_gps", "start_location_lat_closest_gps"), crs = wgs)

mapview(smpl_burst_sf, zcol = "flight_type_sm3", alpha = 0, cex = 3) + mapview(smpl_8sec, zcol = "circling")

#try another sample
smpl_burst_sf2 <- or_summaries_w_gps %>% 
  filter(individual_local_identifier == "D326_192" & imu_burst_id %in% c(757:756)) %>% 
  drop_na("location_lat_closest_gps") %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)

smpl_8sec <- or_8sec_circling_ys_no %>% 
  filter(between(start_timestamp, head(smpl_burst_sf2$timestamp, 1), tail(smpl_burst_sf2$timestamp, 1))) %>% 
  drop_na(start_location_long_closest_gps) %>% 
  st_as_sf(coords = c("start_location_long_closest_gps", "start_location_lat_closest_gps"), crs = wgs)

mapview(smpl_burst_sf2, zcol = "flight_type_sm3", alpha = 0, cex = 3) + mapview(smpl_8sec, zcol = "circling")




### explore the distribution of yaw and roll ####
#scatterplot of cumulative yaw and cumulative roll
ggplot() +
  geom_point(data = or_summaries_w_gps, aes(x = cumulative_yaw, y = cumulative_roll))

#-------------------------------------------------------------------------------------
# STEP3: Calculate degree of handedness for each day
#-------------------------------------------------------------------------------------

#look at the relationship between yaw and roll. they should be highly correlated
or_8sec_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_Jul24.rds")

#---------------------- calculate laterality for each day
#assign circling vs non-circling to values of yaw
laterality_8sec <- or_8sec_summaries_w_gps %>% #this has one row per 8-9-sec burst 
  mutate(burst_dur2 = difftime(end_timestamp, start_timestamp)) %>% 
  filter(burst_dur2 > 6.5 &  #remove bursts that are shorter than 6.5 seconds
           cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45) %>%  #only keep thermalling flight. use a threshold of 45 degrees in 8 seconds.
  mutate(bank_direction = ifelse(cumulative_roll_8sec < 0, "left",
                                 ifelse(cumulative_roll_8sec > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_8sec < 0, "left",
                                    ifelse(cumulative_yaw_8sec > 0, "right", "straight")),
         unique_date = as.character(date(start_timestamp)),
         yday = yday(start_timestamp)) %>% #this will make it easier to filter later for life cycle stages 
  group_by(individual_local_identifier, unique_date) %>% #group by individual ID and day
  summarise(n_bursts_in_day = n(),
            yday = head(yday, 1),
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            heading_left = sum(heading_direction == "left"),
            heading_right = sum(heading_direction == "right"),
            heading_straight = sum(heading_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
            laterality_heading = (heading_right - heading_left)/(heading_right + heading_left),
            .groups = 'drop') %>% 
  as.data.frame()

saveRDS(laterality_8sec, file = "laterality_index_per_day.rds")

##explore
X11(width = 12, height = 12)
ggplot(data = laterality_8sec, aes(x = unique_date, y = laterality_bank)) +
  geom_point() +
  geom_smooth( method = "loess") +
  facet_wrap(vars(individual_local_identifier))

### one day has > 1000 bursts. explore why
long_day <- laterality_8sec %>%  filter(n_bursts_in_day > 1000)
long_day_all <- or_8sec_summaries_w_gps %>%  filter(as.character(date(start_timestamp)) == long_day$unique_date & individual_local_identifier == long_day$individual_local_identifier)

long_day_sf <- long_day_all %>% 
  drop_na( "start_location_long_closest_gps" ) %>% 
  st_as_sf(coords = c("start_location_long_closest_gps", "start_location_lat_closest_gps"), crs = wgs) # seems legit! lol

### there are lots of days with one observation

short_day <- laterality_8sec %>%  filter(n_bursts_in_day == 1) %>% 
  mutate(unique_ind_day = paste0(individual_local_identifier, "_", unique_date))


short_day_all <- or_8sec_summaries_w_gps %>%
  mutate(unique_ind_day = paste0(individual_local_identifier, "_", as.character(date(start_timestamp)))) %>% 
  filter(unique_ind_day %in% short_day$unique_ind_day)

short_day_sf <- short_day_all %>% 
  drop_na("start_location_long_closest_gps" ) %>% 
  st_as_sf(coords = c("start_location_long_closest_gps", "start_location_lat_closest_gps"), crs = wgs) # seems legit! lol

mapview(short_day_sf) #mainly during the wintering OR the breeding season.


#-------------------------------------------------------------------------------------
# STEP3: Calculate degree of handedness for each 8-second burst
#-------------------------------------------------------------------------------------

#open quat data- aggregated for every second, matched with gps-informed flight classification: prepared in 01b_imu_processing.r
or_summaries_w_gps <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8secIDs_Jul24.rds") %>% #this version also has the IDs assigned to 8-sec bursts.
  drop_na(individual_local_identifier)

#---------------------- calculate laterality for each 8-second burst
#this is not filtered for circling/not OR the length of bursts

(st_time <- Sys.time())
laterality_1sec <- or_summaries_w_gps %>% #this has one row per seconds
  mutate(bank_direction = ifelse(cumulative_roll < 0, "left",
                                 ifelse(cumulative_roll > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw < 0, "left",
                                    ifelse(cumulative_yaw > 0, "right", "straight"))) %>% #this will make it easier to filter later for life cycle stages
  group_by(burst_id_8sec) %>% #group by 8-sec burst IDs, which are unique across individuals
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
            across(c(timestamp, timestamp_closest_gps, location_lat_closest_gps, location_long_closest_gps, height_above_ellipsoid_closest_gps, ground_speed_closest_gps, heading_closest_gps), 
                   ~head(.,1),
                   .names = "start_{.col}"),
            across(c(flight_type_sm2, flight_type_sm3, track_flight_seg_id),
                   ~Mode(.), #make sure the Mode function is defined
                   .names = "{.col}_Mode"),
            end_timestamp = tail(timestamp,1),
            n_records = n(),
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            heading_left = sum(heading_direction == "left"),
            heading_right = sum(heading_direction == "right"),
            heading_straight = sum(heading_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
            laterality_heading = (heading_right - heading_left)/(heading_right + heading_left),
            .groups = "keep") %>% 
  ungroup() %>% 
  as.data.frame()
Sys.time() - st_time # 8 minutes
  #this is basically the same as GPS_matched_or_w_summaries_8sec_Jul24.rds in 01b_imu_processing.r, but with laterality index added

saveRDS(laterality_1sec, file = "laterality_index_per_8sec_burst.rds")

#-------------------------------------------------------------------------------------
# STEP4: Calculate degree of handedness for each second
#-------------------------------------------------------------------------------------
#the code needs to be modified. it is really just a copy paste from above
# #open quat data- in original resolution, with one second per row (not aggregated) but with euler angles calculated: prepared in 01b_imu_processing.r
# or_angles_raw <- readRDS("quat_angles_secs_apr24.rds") %>% 
#   drop_na(individual_local_identifier)
# 
# (st_time <- Sys.time())
# laterality_1sec <- or_summaries_w_gps %>% #this has one row per seconds
#   mutate(bank_direction = ifelse(cumulative_roll < 0, "left",
#                                  ifelse(cumulative_roll > 0, "right", "straight")),
#          heading_direction = ifelse(cumulative_yaw < 0, "left",
#                                     ifelse(cumulative_yaw > 0, "right", "straight"))) %>% #this will make it easier to filter later for life cycle stages
#   group_by(burst_id_8sec) %>% #group by 8-sec burst IDs, which are unique across individuals
#   summarize(across(c(study_id, individual_taxon_canonical_name, individual_local_identifier, orientation_quaternions_sampling_frequency,
#                      tag_local_identifier, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
#             across(c(yaw_sd, roll_sd, pitch_sd,
#                      yaw_mean, roll_mean, pitch_mean,
#                      yaw_min, roll_min, pitch_min,
#                      yaw_max, roll_max, pitch_max), 
#                    ~mean(., na.rm = T),
#                    .names = "mean_{.col}"), #rename the column, because these values are the mean of the sd, max, min values across the 8 seconds
#             across(c(cumulative_yaw, cumulative_roll, cumulative_pitch),
#                    ~sum(., na.rm = T),
#                    .names = "{.col}_8sec"),
#             across(c(timestamp, timestamp_closest_gps, location_lat_closest_gps, location_long_closest_gps, height_above_ellipsoid_closest_gps, ground_speed_closest_gps, heading_closest_gps), 
#                    ~head(.,1),
#                    .names = "start_{.col}"),
#             across(c(flight_type_sm2, flight_type_sm3, track_flight_seg_id),
#                    ~Mode(.), #make sure the Mode function is defined
#                    .names = "{.col}_Mode"),
#             end_timestamp = tail(timestamp,1),
#             n_records = n(),
#             bank_left = sum(bank_direction == "left"),
#             bank_right = sum(bank_direction == "right"),
#             bank_straight = sum(bank_direction == "straight"),
#             heading_left = sum(heading_direction == "left"),
#             heading_right = sum(heading_direction == "right"),
#             heading_straight = sum(heading_direction == "straight"),
#             laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
#             laterality_heading = (heading_right - heading_left)/(heading_right + heading_left),
#             .groups = "keep") %>% 
#   ungroup() %>% 
#   as.data.frame()
# Sys.time() - st_time # 8 minutes
# #this is basically the same as GPS_matched_or_w_summaries_8sec_Jul24.rds in 01b_imu_processing.r, but with laterality index added
# 
# saveRDS(laterality_1sec, file = "laterality_index_per_8sec_burst.rds")


