# script for preparing the data for Safi et al 2025.
# Elham Nourani, PhD. enourani@ab.mgp.de

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

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#source the imu conversion functions. do i need these here?
#source("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/HB_ontogeny/MS1_IMU_classification/00_imu_diy.r")

#create a funciton for mode
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
## Step 1: Summarize per 8-sec burst                                         #####
#---------------------------------------------------------------------------------

#Start with data that was summarized for each second in 02_imu_processing.r
#the aim of this step is to calculate summaries of the quaternion-derived angles for each 8-second burst
#and to determine whether the individual was circling during a burst (based on total yaw) and the 
#laterality index (based on total roll)


#### ----------------------- Assign a unique burst id to each 8-second burst 
seconds_summaries <- readRDS("quat_summaries_1sec_Jul24.rds") %>% 
  drop_na(individual_local_identifier)

#Add a new column with the number of rows of each imu_burst. We have roughly one second of data per row.
burst_n <- seconds_summaries %>%
  group_by(individual_local_identifier, imu_burst_id) %>%
  mutate(n_of_rows = n()) %>%
  ungroup() %>%
  as.data.frame()

#Add a new column with new burst IDs for bursts that are longer than 9 seconds. For the rest, keep the original burst IDs
or_seconds_8_sec_id <- burst_n %>%
  filter(n_of_rows > 9) %>%
  group_by(individual_local_identifier, imu_burst_id) %>%
  mutate(burst_id_8sec = paste0(individual_local_identifier, "_", imu_burst_id, "_", ceiling(row_number() / 8))) %>%
  ungroup() %>%
  as.data.frame() %>%
  #bind it back to the rest of the data
  full_join(burst_n %>% filter(n_of_rows <= 9)) %>%
  mutate(burst_id_8sec = ifelse(is.na(burst_id_8sec),
                                paste0(individual_local_identifier, "_", imu_burst_id, "_", 1), burst_id_8sec)) %>% #assign unique 8-sec burst ID to bursts that were less than 9 seconds long to begin with. use 1 because they each have one sub-burst.
  as.data.frame()

#### ----------------------- Calculate summaries of angles and the laterality index for each 8-sec burst

(start_time <- Sys.time())
eight_sec <- or_seconds_8_sec_id %>% 
  #calculate direction of banking based on the roll angle. this will assign the direction for each second (each row)
  mutate(bank_direction = ifelse(cumulative_roll < 0, "left",
                                 ifelse(cumulative_roll > 0, "right", "straight"))) %>% 
  group_by(burst_id_8sec) %>% #this grouping variable has unique values for individuals and bursts. 
  arrange(timestamp, .by_group = T) %>% 
  #aggregate summarized angles over 8 seconds
  summarize(across(c(study_id, individual_taxon_canonical_name, individual_local_identifier, orientation_quaternions_sampling_frequency,
                     tag_local_identifier, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
            across(c(yaw_sd, roll_sd, pitch_sd,
                     yaw_mean, roll_mean, pitch_mean,
                     yaw_min, roll_min, pitch_min,
                     yaw_max, roll_max, pitch_max), 
                   ~mean(., na.rm = T),
                   .names = "mean_{.col}"), #rename the column, because these values are the mean of the sd, max, min values across the 8 seconds
            across(c(cumulative_yaw, cumulative_roll, cumulative_pitch, yaw_mean, roll_mean),
                   ~sum(., na.rm = T),
                   .names = "{.col}_sum_8sec"),
            start_timestamp = head(timestamp, 1),
            end_timestamp = tail(timestamp, 1),
            n_records = n(),
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
            .groups = "keep") %>% 
  ungroup %>% #now each row corresponds to one 8-second burst
  #FILTER: remove short bursts. some bursts are shorter than 8 seconds. use 6 seconds as a cut-off
  filter(n_records >= 6) %>% 
  #group_by(individual_local_identifier) %>%
  #calculate circling based on total yaw angle in the 8-second brust
  mutate(circling_status = case_when(
    between(cumulative_yaw_sum_8sec, -10, 10) ~ "straight",
    cumulative_yaw_sum_8sec >= 45 | cumulative_yaw_sum_8sec <= -45 ~ "circling",
    .default = "shallow circling"
  )) %>% 
  #calculate laterality index for each burst
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  #create a binary variable for lateralised vs ambidextrous flight in the burst
  mutate(laterality_bi = ifelse(laterality_dir == "ambidextrous", 0, 1), 
         abs_cum_yaw = abs(cumulative_yaw_sum_8sec)) %>% 
  as.data.frame()
(Sys.time() - start_time) #45 minutes

#saveRDS(eight_sec, file = "matched_GPS_IMU/GPS_matched_or_w_summaries_8sec.rds")

#---------------------------------------------------------------------------------
## Step 2: Environmental annotation                                          #####
#---------------------------------------------------------------------------------

#to match wind conditions and IMU data, the intermediate step was to annotate the GPS locations with wind speed (02_wind_annotation.r)
#in this step, the closest GPS point to each IMU burst will be identified and the wind speed associated with that point will now be associated with the IMU burst


#open gps data, matched with wind in 02_wind_annotation.r (a list with one element per year) and create one dataframe
#gps_ls <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24_wind.rds") %>% 
gps_ls <- ann_ls %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  arrange(individual_local_identifier) %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

#create a list with one element for each individual
#or_ls <- laterality_circling %>% 
or_ls <- eight_sec %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  arrange(individual_local_identifier) %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

# Define a function to find the closest GPS information and associate it with orientation data. Please note that the column names are hard-coded. modify to match your data.
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
identical(sapply(gps_ls, function(x) x$individual_local_identifier[1]), sapply(or_ls, function(x) x$individual_local_identifier[1])) #should return TRUE

# Run the function on the list of gps and list of imu. This will create a list of data frames with or data and associated GPS information
(b <- Sys.time())
or_w_gps <- map2(or_ls, gps_ls, ~ find_closest_gps(or_data = .x, gps_data = .y))
Sys.time() - b # 23 mins

#add a column comparing or and gps timestamps.
or_w_gps_df <- lapply(or_w_gps, function(x){
  x2 <- x %>% 
    mutate(imu_gps_timediff_sec = if_else(is.na(timestamp_closest_gps_raw), NA, difftime(start_timestamp, timestamp_closest_gps_raw, units = "mins") %>%  as.numeric()))
  x2
}) %>% 
  bind_rows()

sum(is.na(or_w_gps_df$location_long_closest_gps_raw)) #8573 there are still some rows that don't get assigned a gps location

#remove the rows with no assigned gps
saveRDS(or_w_gps_df, file = "annotated_gps_w_wind_public_prep.rds")

#---------------------------------------------------------------------------------
## Step 3: Life stage annotation                                             #####
#---------------------------------------------------------------------------------

#life-cycle stages from L03a_tests_per_day.r
life_cycle <- readRDS("updated_life_cycle_nov24.rds")

#add age as days since tagging, and assign life stage based on HMM analysis
or_w_gps_df_LS <- or_w_gps_df %>% 
  mutate(unique_date = as.Date(start_timestamp)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, deployment_dt_utc, first_exploration, migration_start, migration_end), by = "individual_local_identifier") %>% 
  #remove data before deployment
  filter(start_timestamp >= deployment_dt_utc) %>% 
  group_by(individual_local_identifier) %>% 
  rowwise() %>% 
  mutate(days_since_tagging = floor(difftime(start_timestamp, deployment_dt_utc, unit = "days")),
         life_stage = case_when(
           between(unique_date, migration_start, migration_end) ~ "migration",
           between(unique_date, first_exploration, migration_start) ~ "post-fledging",
           unique_date < first_exploration ~ "pre-fledging",
           unique_date > migration_end ~ "wintering",
           TRUE ~ NA_character_
         )) %>% 
  ungroup() %>% 
  as.data.frame()

#saveRDS(or_w_gps_df_LS, file = "laterality_w_gps_wind_LS_no_filter.rds")

#---------------------------------------------------------------------------------
## Step 4: Data filtering                                                    #####
#---------------------------------------------------------------------------------

#open prepared data that has not been filtered yet.
#levels of filtering:
#1) only keep circling flight (MUST do before LI calculations at the daily, life stage, and individual levels)
#2) thin the data, so that the consecutive bursts are at least 2 min apart. (MUST do before LI calculations)
#3)only keep data after fledging


or_w_gps_flt <- or_w_gps_df_LS%>% 
  
  #FILTER: remove non-circling bursts
  filter(circling_status == "circling") %>% 
  #calculate time lag between remaining bursts
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = TRUE) %>% 
  mutate(time_lag_sec = if_else(row_number() == 1, 0, 
                                difftime(start_timestamp, lag(start_timestamp), units = "secs") %>% as.numeric())) %>% 
  ungroup() %>% 
  
  #FILTER: thin the data, so that the consecutive bursts are at least 2 min apart
  filter(time_lag_sec >= 120) %>% 
  
  #FILTER: filter out pre-fledging
  filter(life_stage != "pre-fledging") %>% 
  as.data.frame()


#---------------------------------------------------------------------------------
## Step 5: Calculate laterality index for days, life stages, and inds        #####
#---------------------------------------------------------------------------------

filtered_w_LI <- or_w_gps_flt %>% 
  group_by(individual_local_identifier) %>% 
  arrange(days_since_tagging, .by_group = T) %>% 
  mutate(weeks_since_tagging = ceiling(days_since_tagging/7),  #not all individuals have a week 1. 
         bank_direction = ifelse(mean_roll_mean < 0, "left",
                                 ifelse(mean_roll_mean > 0, "right", "straight"))) %>% 
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
  
  #------------------laterality index for each life stage
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
filtered_w_LI$laterality_dir_day <- factor(filtered_w_LI$laterality_dir_day, levels = c("right_handed", "ambidextrous", "left_handed"))
filtered_w_LI$laterality_dir_stage <- factor(filtered_w_LI$laterality_dir_stage, levels = c("right_handed", "ambidextrous", "left_handed"))

saveRDS(filtered_w_LI, file = "thinned_laterality_w_gps_wind_all_filters2_public_prep.rds")

#saveRDS(filtered_w_LI, file = "thinned_laterality_w_gps_wind_all_filters2.rds")


