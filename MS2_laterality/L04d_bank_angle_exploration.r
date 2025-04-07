## figuring out bank angle!! calculating bank angle at the original resolution (20 Hz)
#Feb 12.2025. Konstanz, DE
#Elham Nourani, enourani@ab.mpg.de



library(tidyverse)
library(sf)
library(mapview)

setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")
source("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/HB_ontogeny/MS1_IMU_classification/00_imu_diy.r")


## step 1: bank angle function--------------------------------------------------------------

bank_angle <- function(roll_rad, pitch_rad) {
  tibble(roll_rad = roll_rad, pitch_rad = pitch_rad) %>%
    mutate(
      bank_rad = asin(sin(roll_rad) * cos(pitch_rad)),
      bank_deg = bank_rad * 180 / pi,
      bank_diff = bank_deg - lag(bank_deg, default = first(bank_deg)),
      bank_diff_circ = case_when(
        bank_diff > 180 ~ bank_diff - 360,
        bank_diff < -180 ~ bank_diff + 360,
        TRUE ~ bank_diff
      ),
      cumulative_bank_circ = cumsum(bank_diff_circ),
      cumulative_bank = cumsum(bank_diff)
    ) %>%
    summarise(
      bank_sd = sd(bank_deg),
      bank_mean = mean(bank_deg),
      bank_max = max(bank_deg),
      bank_min = min(bank_deg),
      cumulative_bank = last(cumulative_bank),
      cumulative_bank_circ = last(cumulative_bank_circ)
    )
}

## step 2: calc bank angle for a subset of data in original resolution--------------------------------------------------------------
#input a sample of data (same as the one used in L04_full_workflow.r)

#start with data that was summarized for each second in MS1_IMU_classification/01b_imu_processing.r
#this is also already matched with GPS
one_sec <- readRDS("matched_GPS_IMU/GPS_matched_or_w_summaries_8secIDs_Jul24.rds")

sample_b <- one_sec %>% 
  filter(individual_local_identifier == "D329_013" &
           imu_burst_id %in% c(2672, 2674, 2676, 2677, 2678)) %>% 
  #calculate bank angle for the original resolution of 20 HZ
  rowwise() %>%
  mutate(
    bank_angle = list(bank_angle(
      pitch = strings_to_numeric(pitch),
      roll = strings_to_numeric(roll)
    ))
  ) %>%
  ungroup() %>%
  unnest(bank_angle) %>% 
  #add a variable just for the direction of bank
  mutate(bank_direction_cumulative = case_when(
    cumulative_bank > 0 ~ "positive",
    cumulative_bank < 0 ~ "negative",
    cumulative_bank == 0 ~ "straight"
  ),
  bank_direction_circ_cumulative = case_when(
    cumulative_bank_circ > 0 ~ "positive",
    cumulative_bank_circ < 0 ~ "negative",
    cumulative_bank_circ == 0 ~ "straight"
  ),
  bank_direction_mean = case_when(
    bank_mean > 0 ~ "positive",
    bank_mean < 0 ~ "negative",
    bank_mean == 0 ~ "straight"
  )) %>% 
  as.data.frame()


plot(sample_b$cumulative_bank, sample_b$cumulative_bank_circ) ## these two are identical. at the second scale, the range of cumulative bank doesnt exceed 180 degrees anyway

#then implement the angle correction at the 8 second scale?
#try it. because summing up the variables, i get a range of -400 to 400 degrees

sample_b_8sec <- sample_b %>% 
  group_by(burst_id_8sec) %>% #this grouping variable has unique values for individuals and bursts. 
  arrange(timestamp, .by_group = T) %>% 
  mutate(bank_diff_8sec = bank_mean - lag(bank_mean, default = first(bank_mean)), #calc differences between mean bank across 8 sec
  bank_diff_8sec = case_when(
    bank_diff_8sec > 180 ~ bank_diff_8sec - 360,
    bank_diff_8sec < -180 ~ bank_diff_8sec + 360,
    TRUE ~ bank_diff_8sec
  ),
  cumulative_bank_8sec = cumsum(bank_diff_8sec)
  ) %>% 
  summarize(across(c(study_id, individual_taxon_canonical_name, individual_local_identifier, orientation_quaternions_sampling_frequency,
                     tag_local_identifier, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
            across(c(timestamp, timestamp_closest_gps, location_lat_closest_gps, location_long_closest_gps, height_above_ellipsoid_closest_gps, ground_speed_closest_gps, heading_closest_gps), 
                   ~head(.,1),
                   .names = "start_{.col}"),
            across(c(flight_type_sm2, flight_type_sm3, track_flight_seg_id),
                   ~Mode(.), #make sure the Mode function is defined
                   .names = "{.col}_Mode"),
            end_timestamp = tail(timestamp,1),
            cum_bank_8sec = last(cumulative_bank_8sec),
            .groups = "keep")
  
## step 3: plotting--------------------------------------------------------------

#create an sf object
sample_sf <- sample_b %>% 
  drop_na("location_lat_closest_gps") %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)

#cumulative bank values
mapview(sample_sf , zcol = "cumulative_bank") 

#direction when using cumulative sum
mapview(sample_sf , zcol = "bank_direction_cumulative") 

#direction when using the mean
mapview(sample_sf , zcol = "bank_direction_mean") 
#mean bank values
mapview(sample_sf , zcol = "bank_mean") 


#compare to roll
#direction when using the mean
mapview(sample_sf , zcol = "cumulative_roll")

plot(sample_sf$roll_mean, sample_sf$bank_mean)

plot(sample_sf$cumulative_bank, sample_sf$bank_mean)
plot(sample_sf$cumulative_roll, sample_sf$roll_mean)

#mean bank works better than cumulative bank.



### yaw and roll

mapview(sample_sf , zcol = "cumulative_yaw") 

