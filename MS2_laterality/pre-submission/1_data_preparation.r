# script for testing the hypotheses for Safi et al 2025.
# Elham Nouani, PhD. 
#17.09.2024, Konstanz, DE


#PNAS figure guidelines: width for 1 column = 3.42, 1.5 columns = 4.5 in, 2 columns = 7 in; max height = 9 in . Text should be at least 6-8 pts

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)
library(viridis)
library(lme4)
library(mgcv)
library(mgcViz)
library(INLA)
library(terra)
library(performance)
library(corrr)
library(gridExtra)
library(patchwork)
library(rptR)
library(xtable) #for exporting latex tables
library(ggh4x) # devtools::install_github("teunbrand/ggh4x") #allows modifying colors of facets in ggplot

setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#source the imu conversion functions
source("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/HB_ontogeny/MS1_IMU_classification/00_imu_diy.r")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#---------------------------------------------------------------------------------
## Step 0: Life cycle descriptive summaries                                  #####
#---------------------------------------------------------------------------------
#to be reported in the paper

#individuals with incomplete tracks

incomplete <- c("D329_015", "D326_193", "D324_513", "D320_474", "D299_270", "D299_269", "D225_236")

#life-cycle stages from L03a_tests_per_day.r
life_cycle <- readRDS("updated_life_cycle_nov24.rds") %>% 
  mutate(age_at_first_expl = (first_exploration - as.Date(deployment_dt_utc.x)) + 30,
         in_natal = migration_start - first_exploration,
         migr_dur = ifelse(individual_local_identifier %in% incomplete, NA, migration_end - migration_start)) %>% 
  summarize(in_natal_avg = mean(in_natal, na.rm = T),
            in_natal_min = min(in_natal, na.rm = T),
            in_natal_max = max(in_natal, na.rm = T),
            in_natal_sd = sd(in_natal, na.rm = T),
            migr_dur_avg = min(migr_dur, na.rm = T),
            migr_dur_min = mean(migr_dur, na.rm = T),
            migr_dur_max = max(migr_dur, na.rm = T),
            migr_dur_sd = sd(migr_dur, na.rm = T))

#tracking duration
#open data from step 4 below
filtered_w_LI <- readRDS("thinned_laterality_w_gps_wind_all_filters.rds")

mean(filtered_w_LI$days_since_tagging)
sd(filtered_w_LI$days_since_tagging)
IQR(filtered_w_LI$days_since_tagging)

#---------------------------------------------------------------------------------
## Step 1: Summarize per 8-sec burst & calc bank angle                       #####
#---------------------------------------------------------------------------------

#start with data that was summarized for each second in MS1_IMU_classification/01b_imu_processing.r
#this is also already matched with GPS

one_sec <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8secIDs_Jul24.rds")

# calculate summaries of angles for each 8-sec burst
(start_time <- Sys.time())
eight_sec <- one_sec %>% 
  #calculate bank angle for each 1/20th of a second, then calculate the mean (make sure roll and pitch are in radians)
  mutate(
    bank_angle_rad_mean = map2_dbl(
      roll, pitch,
      ~ mean(asin(sin(strings_to_numeric(.x)) * cos(strings_to_numeric(.y))))
    ),
    bank_angle_deg_mean = bank_angle_rad_mean * 180 / pi) %>%
  #calc laterality direction
  mutate(bank_direction = ifelse(bank_angle_deg_mean < 0, "left",
                                    ifelse(bank_angle_deg_mean > 0, "right", "straight")),
         roll_direction = ifelse(cumulative_roll < 0, "left",
                                 ifelse(cumulative_roll > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw < 0, "left",
                                    ifelse(cumulative_yaw > 0, "right", "straight"))) %>% # save up to this stage to compare second-scale laterality with the roll laterality from before: saveRDS(sec, file = "one_sec_BA_laterality.rds")
  #summarize for each 8 second burst
  group_by(burst_id_8sec) %>% #this grouping variable has unique values for individuals and bursts. 
  arrange(timestamp, .by_group = T) %>% 
  #aggregate summarized angles over 8 seconds
  summarize(across(c(study_id, individual_taxon_canonical_name, individual_local_identifier, orientation_quaternions_sampling_frequency,
                     tag_local_identifier, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
            across(c(yaw_sd, roll_sd, pitch_sd,
                     yaw_mean, roll_mean, pitch_mean,
                     yaw_min, roll_min, pitch_min,
                     yaw_max, roll_max, pitch_max,
                     bank_angle_deg_mean), 
                   ~mean(., na.rm = T),
                   .names = "mean_{.col}"), #rename the column, because these values are the mean of the sd, max, min values across the 8 seconds
            across(c(cumulative_yaw, cumulative_roll, cumulative_pitch, yaw_mean, roll_mean),
                   ~sum(., na.rm = T),
                   .names = "{.col}_sum_8sec"), #in the old version, this was "{.col}_8sec"
            across(c(timestamp, timestamp_closest_gps, location_lat_closest_gps, location_long_closest_gps, height_above_ellipsoid_closest_gps, ground_speed_closest_gps, heading_closest_gps), 
                   ~head(.,1),
                   .names = "start_{.col}"),
            across(c(flight_type_sm2, flight_type_sm3, track_flight_seg_id),
                   ~Mode(.), #make sure the Mode function is defined
                   .names = "{.col}_Mode"),
            end_timestamp = tail(timestamp,1),
            n_records = n(), #n of observations for this burst
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            roll_left = sum(roll_direction == "left"),
            roll_right = sum(roll_direction == "right"),
            roll_straight = sum(roll_direction == "straight"),
            heading_left = sum(heading_direction == "left"),
            heading_right = sum(heading_direction == "right"),
            heading_straight = sum(heading_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left + bank_straight),
            laterality_roll = (roll_right - roll_left)/(roll_right + roll_left + roll_straight),
            laterality_heading = (heading_right - heading_left)/(heading_right + heading_left + heading_straight),
            bank_angle_deg_sd = sd(bank_angle_deg_mean, na.rm = T), #sd of bank angle over the 8 seconds
            bank_angle_deg_var = var(bank_angle_deg_mean, na.rm = T), #variance of bank angle over the 8 seconds
            .groups = "keep") %>% 
  ungroup %>% 
  as.data.frame()
(Sys.time() - start_time) #12 minutes (on labtop)

saveRDS(eight_sec, file = "matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_BA.rds")


#sanity check for bank angle calculation -------------------------------------------------------------------------------------------

eight_sec <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_BA.rds")

ggplot(eight_sec, aes(x = cumulative_yaw_sum_8sec, y = mean_bank_angle_deg_mean)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", col = "blue") +  # Add linear model
  labs(title = "bank angle calculated for each 1/20 of second, then averaged for each second, then averaged for 8-seconds"
  )

ggplot(eight_sec, aes(x = mean_roll_mean, y = mean_bank_angle_deg_mean)) +
  geom_point() +  # Add points
  geom_smooth(method = "lm", col = "blue") +  # Add linear model
  labs(title = "bank angle calculated for each 1/20 of second, then averaged for each second, then averaged for 8-seconds"
  )

#---------------------------------------------------------------------------------
## Step 3: Environmental annotation                                          #####
#---------------------------------------------------------------------------------

#### ----------------------- use the 8-sec data and calculate handedness for each 8 second burst
#### Prepare the data first, add filtering steps afterwards

laterality_circling <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8sec_BA.rds") %>% 
  #FILTER: remove short bursts 
  filter(n_records >= 6) %>% 
  group_by(individual_local_identifier) %>% 
  mutate(individual_local_identifier = as.factor(individual_local_identifier),
         individual_local_identifier2 = individual_local_identifier,
         individual_local_identifier3 = individual_local_identifier) %>% #dublicate individual ID to use for inla random effects specification
  #create a column for circling status based on yaw
  mutate(circling_status = case_when(
    between(cumulative_yaw_sum_8sec, -10, 10) ~ "straight",
    cumulative_yaw_sum_8sec >= 45 | cumulative_yaw_sum_8sec <= -45 ~ "circling",
    .default = "shallow circling"
  )) %>% 
  mutate(abs_cum_yaw = abs(cumulative_yaw_sum_8sec)) %>% 
  as.data.frame()

#### ----------------------- environmental annotation

#6940 rows don't have a gps location associated with them.... so, redo GPS-matching and increase the time window to one hour, to match that of the env. data
#open raw gps data, matched with wind in L04a_env_annotation.r

gps_ls <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24_wind.rds") %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  arrange(individual_local_identifier) %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

#create a list of one element for each individual
or_ls <- laterality_circling %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  arrange(individual_local_identifier) %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

# Define a function to find the closest GPS information and associate it with orientation data
find_closest_gps <- function(or_data, gps_data, time_tolerance = 60 * 60) {
  map_df(1:nrow(or_data), function(h) {
    or_row_time <- or_data[h, "start_timestamp"]
    gps_sub <- gps_data %>%
      filter(between(timestamp, or_row_time - time_tolerance, or_row_time + time_tolerance))
    
    if (nrow(gps_sub) >= 1) {
      time_diff <- abs(difftime(gps_sub$timestamp, or_row_time, units = "secs"))
      min_diff <- which.min(time_diff)
      or_data[h, c("timestamp_closest_gps_raw", "location_long_closest_gps_raw", "location_lat_closest_gps_raw", "height_above_ellipsoid_closest_gps_raw", "ground_speed_closest_gps_raw", 
                   "heading_closest_gps",  "u_900", "v_900", "wind_speed")] <- 
        gps_sub[min_diff, c("timestamp", "location_long", "location_lat", "height_above_ellipsoid", "ground_speed", "heading", "u_900", "v_900", "wind_speed")]
    } else {
      or_data[h, c("timestamp_closest_gps_raw", "location_long_closest_gps_raw", "location_lat_closest_gps_raw", "height_above_ellipsoid_closest_gps_raw", "ground_speed_closest_gps_raw", 
                   "heading_closest_gps", "u_900", "v_900", "wind_speed")] <- NA
    }
    return(or_data[h, ])
  })
}

#make sure the order of individuals is the same in the two lists
# Create a list of data frames with or data and associated GPS information
(b <- Sys.time())
or_w_gps <- map2(or_ls, gps_ls, ~ find_closest_gps(or_data = .x, gps_data = .y))
Sys.time() - b # 32 mins

#add a column comparing or and gps timestamps. then save one file per individual.
or_w_gps_df <- lapply(or_w_gps, function(x){
  x2 <- x %>% 
    mutate(imu_gps_timediff_sec = if_else(is.na(timestamp_closest_gps_raw), NA, difftime(start_timestamp, timestamp_closest_gps_raw, units = "mins") %>%  as.numeric()))
  x2
}) %>% 
  bind_rows()

sum(is.na(or_w_gps_df$location_long_closest_gps_raw)) #406
sum(is.na(or_w_gps_df$start_location_long_closest_gps)) #6940
#post ba analysis, there are way more rows. why?

#there are still many rows with no assigned gps
no_gps <- or_w_gps_df %>% 
  filter(is.na(timestamp_closest_gps_raw))

saveRDS(or_w_gps_df, file = "annotated_gps_w_wind_BA.rds")

#### ----------------------- add life-cycle stage

### ALSO ADD DAY SINCE TAGGING HERE
#life-cycle stages from L03a_tests_per_day.r
life_cycle <- readRDS("updated_life_cycle_nov24.rds")

or_w_gps_df_LS <- or_w_gps_df %>% 
  #select(-deployment_dt_utc) %>% #remove deployment in the dataset and just use the one in life_cycle, which is from the metadata
  mutate(unique_date = as.Date(start_timestamp)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, deployment_dt_utc, first_exploration, migration_start, migration_end), by = "individual_local_identifier") %>% 
  #remove data before deployment
  filter(start_timestamp >= deployment_dt_utc) %>% 
  group_by(individual_local_identifier) %>% 
  rowwise() %>% 
  mutate(life_stage = case_when(
    between(unique_date, migration_start, migration_end) ~ "migration",
    between(unique_date, first_exploration, migration_start) ~ "post-fledging",
    unique_date < first_exploration ~ "pre-fledging",
    unique_date > migration_end ~ "wintering",
    TRUE ~ NA_character_),
    days_since_tagging = floor(difftime(start_timestamp, deployment_dt_utc, unit = "days")) %>% as.numeric()) %>% 
  ungroup() %>% 
  as.data.frame()

saveRDS(or_w_gps_df_LS, file = "laterality_w_gps_wind_LS_no_filter_BA.rds")

#---------------------------------------------------------------------------------
## Step 2: Data filtering                                                    #####
#---------------------------------------------------------------------------------

#open prepared data that has not been filtered yet.
#levels of filtering:
#1) only keep circling flight (MUST do before LI calculations)
#2) thin the data, so that the consecutive bursts are at least 2 min apart.. (MUST do before LI calculations)
#only keep data after fledging
#only keep first 300 days


or_w_gps_flt <- readRDS("laterality_w_gps_wind_LS_no_filter_BA.rds") %>% 
  #FILTER: remove non-circling bursts
  filter(circling_status == "circling") %>% 
  #calculate time lag between remaining bursts
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = TRUE) %>% 
  mutate(time_lag_sec = if_else(row_number() == 1, 0, 
                                difftime(start_timestamp, lag(start_timestamp), units = "secs") %>% as.numeric())) %>% 
  ungroup() %>% 
  #FILTER: thin the data, so that the consecutive bursts are at least 2 min apart.. for consistency
  filter(time_lag_sec >= 120) %>% 
  #FILTER: filter out pre-fledging
  filter(life_stage != "pre-fledging") %>% 
  as.data.frame()


#-------------------------------------------------------------------------
## Step 2: Calculate laterality index for days, life stages, and ind #####
#-------------------------------------------------------------------------

filtered_w_LI <- or_w_gps_flt %>% 
  mutate(laterality_dir_8sec = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  group_by(individual_local_identifier) %>% 
  arrange(days_since_tagging, .by_group = T) %>% 
  mutate(weeks_since_tagging = ceiling(days_since_tagging/7),  #not all individuals have a week 1. 
         bank_direction = ifelse(mean_bank_angle_deg_mean < 0, "left",
                                 ifelse(mean_bank_angle_deg_mean > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_sum_8sec < 0, "left",
                                    ifelse(cumulative_yaw_sum_8sec > 0, "right", "straight"))) %>% 
  ungroup() %>% 
  
  #------------------laterality index for each day
  group_by(individual_local_identifier, days_since_tagging) %>% 
  mutate(total_n = n(),
         bank_left = sum(bank_direction == "left"),
         bank_right = sum(bank_direction == "right"),
         bank_straight = sum(bank_direction == "straight"),
         laterality_bank_day = (bank_right - bank_left)/(bank_right + bank_left + bank_straight)) %>%
  mutate(laterality_dir_day = case_when(
    between(laterality_bank_day, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank_day, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank_day, -0.25, 0.25) ~ "ambidextrous")) %>% 
  ungroup() %>% 
  select(-c(total_n, bank_left, bank_right, bank_straight)) %>% 
  
  #------------------laterality index for life stage
  group_by(individual_local_identifier, life_stage) %>% 
  mutate(total_n = n(),
         bank_left = sum(bank_direction == "left"),
         bank_right = sum(bank_direction == "right"),
         bank_straight = sum(bank_direction == "straight"),
         laterality_bank_stage = (bank_right - bank_left)/(bank_right + bank_left + bank_straight)) %>% 
  mutate(laterality_dir_stage = case_when(
    between(laterality_bank_stage, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank_stage, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank_stage, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  ungroup() %>% 
  select(-c(total_n, bank_left, bank_right, bank_straight)) %>% 
  
  #------------------laterality index for each individual. pooling all data for each ind
  group_by(individual_local_identifier) %>% 
  mutate(total_n = n(),
         bank_left = sum(bank_direction == "left"),
         bank_right = sum(bank_direction == "right"),
         bank_straight = sum(bank_direction == "straight"),
         laterality_bank_ind = (bank_right - bank_left)/(bank_right + bank_left + bank_straight)) %>% 
  mutate(laterality_dir_ind = case_when(
    between(laterality_bank_ind, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank_ind, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank_ind, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  ungroup() %>% 
  select(-c(total_n, bank_left, bank_right, bank_straight)) %>% 
  as.data.frame()

#re-order laterality direction
filtered_w_LI$laterality_dir_8sec <- factor(filtered_w_LI$laterality_dir_8sec, levels = c("right_handed", "ambidextrous", "left_handed"))
filtered_w_LI$laterality_dir_day <- factor(filtered_w_LI$laterality_dir_day, levels = c("right_handed", "ambidextrous", "left_handed"))
filtered_w_LI$laterality_dir_stage <- factor(filtered_w_LI$laterality_dir_stage, levels = c("right_handed", "ambidextrous", "left_handed"))

saveRDS(filtered_w_LI, file = "thinned_laterality_w_gps_wind_all_filters_BA.rds")


