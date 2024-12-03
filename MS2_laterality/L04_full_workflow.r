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

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#---------------------------------------------------------------------------------
## Step 1: Data prep                                                         #####
#---------------------------------------------------------------------------------

#### ----------------------- use the 8-sec data and calculate handedness for each 8 second burst
#### Prepare the data first, add filtering steps afterwards

laterality_circling <- readRDS("laterality_index_per_8sec_burst_days_since.rds") %>% 
  #FILTER: remove short bursts 
  filter(n_records >= 6) %>% 
  group_by(individual_local_identifier) %>% 
  #arrange(start_timestamp, .by_group = TRUE) %>% 
  #calculate the time difference between consecutive bursts  ### DO THIS AFTER FITLERING FOR CIRCLING
  #mutate(time_lag_sec = if_else(row_number() == 1, 0, 
  #                              difftime(start_timestamp, lag(start_timestamp), units = "secs") %>% as.numeric())) %>% 
  #ungroup() %>% 
  mutate(days_since_tagging = as.numeric(days_since_tagging),
         individual_local_identifier = as.factor(individual_local_identifier),
         individual_local_identifier2 = individual_local_identifier,
         individual_local_identifier3 = individual_local_identifier) %>% #dublicate individual ID to use for inla random effects specification
  mutate(circling_status = case_when(
    between(cumulative_yaw_8sec, -10, 10) ~ "straight",
    cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45 ~ "circling",
    .default = "shallow circling"
  )) %>% 
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  mutate(laterality_bi = ifelse(laterality_dir == "ambidextrous", 0, 1), #create a binary variable for handedness vs not
         abs_cum_yaw = abs(cumulative_yaw_8sec)) %>% 
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
or_ls <- laterality_circling_thin %>% 
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
Sys.time() - b # 2.5 mins

#add a column comparing or and gps timestamps. then save one file per individual.
or_w_gps_df <- lapply(or_w_gps, function(x){
  x2 <- x %>% 
    mutate(imu_gps_timediff_sec = if_else(is.na(timestamp_closest_gps_raw), NA, difftime(start_timestamp, timestamp_closest_gps_raw, units = "mins") %>%  as.numeric()))
  x2
}) %>% 
  bind_rows()

sum(is.na(or_w_gps_df$location_long_closest_gps_raw)) #406
sum(is.na(or_w_gps_df$start_location_long_closest_gps)) #6940


#there are still many rows with no assigned gps
no_gps <- or_w_gps_df %>% 
  filter(is.na(timestamp_closest_gps_raw))


#### ----------------------- add life-cycle stage
#life-cycle stages from L03a_tests_per_day.r
life_cycle <- readRDS("updated_life_cycle_nov24.rds")

or_w_gps_df_LS <- or_w_gps_df %>% 
  mutate(unique_date = as.Date(start_timestamp)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, migration_start, migration_end, first_exploration), by = "individual_local_identifier") %>% 
  group_by(individual_local_identifier) %>% 
  rowwise() %>% 
  mutate(life_stage = case_when(
    between(unique_date, migration_start, migration_end) ~ "migration",
    between(unique_date, first_exploration, migration_start) ~ "post-fledging",
    unique_date < first_exploration ~ "pre-fledging",
    unique_date > migration_end ~ "wintering",
    TRUE ~ NA_character_
  )) %>% 
  ungroup() %>% 
  as.data.frame()

saveRDS(or_w_gps_df_LS, file = "laterality_w_gps_wind_LS_no_filter.rds")

#FILTER: remove pre-fledging
#or_w_gps_df_LS_pf <- or_w_gps_df_LS %>% 
#  filter(life_stage != "pre-fledging")

#saveRDS(or_w_gps_df_LS_pf, file = "thinned_laterality_w_gps_wind_4min_LS_PF.rds") #only post-fledging data

#---------------------------------------------------------------------------------
## Step 2: Data filtering                                                    #####
#---------------------------------------------------------------------------------

#open prepared data that has not been filtered yet.
#levels of filtering:
#1) only keep circling flight (MUST do before LI calculations)
#2) thin the data, so that the consecutive bursts are at least 2 min apart.. (MUST do before LI calculations)
#only keep data after fledging
#only keep first 300 days


or_w_gps_flt <- readRDS("laterality_w_gps_wind_LS_no_filter.rds") %>% 
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


#-----------------------------------------------------------------------
## Step 2: Calculate laterality index for days, (weeks,) and life stages #####
#-----------------------------------------------------------------------

filtered_w_LI <- or_w_gps_flt %>% 
  group_by(individual_local_identifier) %>% 
  arrange(days_since_tagging, .by_group = T) %>% 
  mutate(weeks_since_tagging = ceiling(days_since_tagging/7),  #not all individuals have a week 1. 
         bank_direction = ifelse(mean_roll_mean < 0, "left",
                                 ifelse(mean_roll_mean > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_8sec < 0, "left",
                                    ifelse(cumulative_yaw_8sec > 0, "right", "straight"))) %>% 
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
  as.data.frame()

#re-order laterality direction
filtered_w_LI$laterality_dir_day <- factor(filtered_w_LI$laterality_dir_day, levels = c("left_handed", "ambidextrous", "right_handed"))
filtered_w_LI$laterality_dir_stage <- factor(filtered_w_LI$laterality_dir_stage, levels = c("left_handed", "ambidextrous", "right_handed"))

saveRDS(filtered_w_LI, file = "thinned_laterality_w_gps_wind_all_filters.rds")

#-----------------------------------------------------------------------
## Step 3: How does laterality vary with latitude                  #####
#-----------------------------------------------------------------------

#### ----------------------- average start and end of migration
population_migr_period <- filtered_w_LI %>% 
  filter(life_stage == "migration") %>% 
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = T) %>% 
  summarize(migration_start = min(days_since_tagging),
            migration_end = max(days_since_tagging)) %>% 
  ungroup() %>% 
  summarize(avg_migration_start = mean(migration_start) %>% round(), 
            avg_migration_end = mean(migration_end) %>% round())

#avg_migration_start 38 
#avg_migration_end 69

#### ----------------------- plots for latitude
#bin latitude into zones
filtered_w_LI_LZ <- filtered_w_LI %>%
  mutate(lat_zone = cut(location_lat_closest_gps_raw, breaks = seq(-30, 70, by = 10), include.lowest = TRUE))

#plot distribution of bank angles at each latitudinal zone
ggplot(filtered_w_LI_LZ, aes(x = mean_roll_mean, y = lat_zone), height = stat(density)) +
  geom_density_ridges(stat = "binline", bins = 200, scale = 0.98, alpha = 0.8, draw_baseline = F)

#re-order laterality direction
filtered_w_LI_LZ$laterality_dir <- factor(filtered_w_LI_LZ$laterality_dir, levels = c("left_handed", "ambidextrous", "right_handed"))



#2D plot for latitudinal zones
X11(width = 3.5, height = 7)
(d_d <- ggplot(data = filtered_w_LI_LZ, aes(x = laterality_dir_day, y = location_lat_closest_gps_raw)) +
    stat_bin2d(aes(fill = ..density..), bins = 50) + #one bar per week
    scale_fill_viridis(option = "plasma") +
    labs(x = "", y = "Latitude") +
    scale_x_discrete(labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    #start and end of migration on average
    #geom_vline(xintercept = c(population_migr_period$avg_migration_start,population_migr_period$avg_migration_end),
    #           linetype = "dashed",  color = "white", linewidth = 0.5, show.legend = TRUE) +
    ggtitle("Density of laterality categories \nbased on laterality index calculated \nfor each day for each individual.") +
    theme_minimal() +
    theme(text = element_text(size = 9)))

ggsave(plot = d_d, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/daily_handedness_lat_2d.jpg", 
       device = "jpg", width = 3.5, height = 7)


#### ----------------------- plots for age

data_300 <- filtered_w_LI_LZ %>% 
  filter(days_since_tagging <= 300)

#2D plot for days since tagging
X11(width = 7, height = 2.5)
(d_2 <- ggplot(data = data_300, aes(x = days_since_tagging, y = laterality_bank_day)) +
    stat_bin2d(aes(fill = stat(density)), geom = "raster", bins = 300/7) + #one bar per week
    #stat_density_2d(aes(fill = after_stat(density)), geom = "polygon") +
    scale_fill_viridis(option = "plasma") +
    labs(x = "Days since tagging", y = "") +
    scale_y_discrete(labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    #start and end of migration on average
    geom_vline(xintercept = c(population_migr_period$avg_migration_start,population_migr_period$avg_migration_end),
               linetype = "dashed",  color = "white", linewidth = 0.5, show.legend = TRUE) +
    ggtitle("Density of handedness categories based on laterality index calculated for \neach day for each individual. \nWhite dashed lines show the migration period.") +
    theme_minimal() +
    theme(text = element_text(size = 9)))

ggsave(plot = d_2, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/daily_handedness_2d.jpg", 
       device = "jpg", width = 7, height = 2.5)

#---------------------------------------------------------------------
## Step 4: Is laterality more likely when the task is difficult? #####
#---------------------------------------------------------------------
#logistic regression: handedness ~ tightness of circles * age . maybe only for individuals with laterality

or_w_gps_df_LS_pf <- readRDS("thinned_laterality_w_gps_wind_4min_LS_PF.rds") %>% 
  
  
  #### ----------------------- filter for day since tagging and z-transform
  
  quantile(or_w_gps_df_LS_pf$days_since_tagging, probs = 0.9) #284 ... 

circling_data <- or_w_gps_df_LS_pf %>% 
  #filter for days since tagging. only keep the 
  filter(days_since_tagging < 300) %>% 
  #make sure to do all the filtering before scaling the variables!!!!
  mutate_at(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
              "abs_cum_yaw", "wind_speed", "days_since_tagging", "location_lat_closest_gps_raw"),
            list(z = ~scale(.))) %>% 
  as.data.frame()

saveRDS(circling_data, file = "thinned_laterality_for_modeling3.rds")

#### ----------------------- look at multi-collinearity

circling_data %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
                  "abs_cum_yaw", "days_since_tagging", "wind_speed", "location_lat_closest_gps_raw")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #sd of roll and pitch = 0.64; sd of yaw and cumulative yaw = 0.7. age and latitude = -0.773

#### ----------------------- exploratory plot
ggplot(circling_data, aes(x = factor(laterality_bi), y = abs(cumulative_yaw_8sec))) +
  geom_boxplot() +
  labs(x = "Laterality", y = "Absolute Cumulative Yaw (8 sec)") +
  theme_minimal()

ggplot(circling_data, aes(x = factor(laterality_dir), y = days_since_tagging)) +
  geom_boxplot() +
  labs(x = "Laterality", y = "") +
  theme_minimal()


#compare n and s hemisphere
ggplot(circling_data, aes(x = factor(laterality_dir), y = location_lat_closest_gps_raw)) +
  geom_boxplot() +
  labs(x = "Laterality", y = "latitude") +
  theme_minimal()

ggplot(circling_data, aes(x = days_since_tagging, y = abs(cumulative_yaw_8sec))) +
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(circling_data, aes(x = days_since_tagging, y = location_lat_closest_gps_raw)) +
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(circling_data, aes(x = days_since_tagging, y = laterality_bank)) +
  geom_point() + 
  geom_smooth(method = "lm")


#### ----------------------- model: binomial logistic regression with inla

circling_data <- readRDS("thinned_laterality_for_modeling3.rds")

#create new data: make predictions for values that fall on a regular grid for optimal visualization 

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations. The max of the variables might be outliers, so use the 90% quantile instead
grd_pitch_yaw <- expand.grid(x = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = max(circling_data$mean_pitch_mean, na.rm = T),  length.out = 50),
                             y = seq(from = min(circling_data$abs_cum_yaw, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(mean_pitch_mean = x,
         abs_cum_yaw = y)  %>% 
  mutate(wind_speed = attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "pitch_yaw")

grd_wind_yaw <- expand.grid(x = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                            y = seq(from = min(circling_data$abs_cum_yaw, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(wind_speed = x,
         abs_cum_yaw = y)  %>% 
  mutate(mean_pitch_mean = attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "wind_yaw")

grd_wind_pitch <- expand.grid(x = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                              y = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = max(circling_data$mean_pitch_mean, na.rm = T),  length.out = 50)) %>% 
  rename(wind_speed = x,
         mean_pitch_mean = y)  %>% 
  mutate(abs_cum_yaw = attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "wind_pitch")


#merge all together
grd_all <- bind_rows(grd_pitch_yaw, grd_wind_yaw, grd_wind_pitch) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation for consistency
  mutate(mean_pitch_mean_z = (mean_pitch_mean - attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:scale'),
         abs_cum_yaw_z = (abs_cum_yaw - attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:scale'),
         wind_speed_z = (wind_speed - attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:scale'),
         days_since_tagging_z = (days_since_tagging - attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:scale'),
         individual_local_identifier = sample(circling_data$individual_local_identifier, nrow(.), replace = T) %>% as.factor(),
         laterality_bi = NA) 


#bin the age variable and append the new data
data <- circling_data %>% 
  dplyr::select(c(intersect(colnames(grd_all), colnames(.)), "life_stage", "location_lat_closest_gps_raw_z")) %>% 
  full_join(grd_all) %>% 
  mutate(individual_local_identifier2 = individual_local_identifier, #repeat individual ID column to be used in the model formula for random effects
         individual_local_identifier3 = individual_local_identifier,
         age_group = inla.group(days_since_tagging_z, n = 100, method = "quantile"), #age will be included as a smooth term
         lat_group = inla.group(location_lat_closest_gps_raw_z, n = 10, method = "quantile"))  


#re-order life stage, so that post-fledging is first
data$life_stage <- factor(data$life_stage, levels = c("post-fledging", "migration", "wintering"))


saveRDS(data, file = "thinned_laterality_for_modeling_w_new_data3.rds")

m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z * wind_speed_z + 
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                 f(individual_local_identifier3, wind_speed_z, model = "iid") + 
                 f(age_group, model = "rw1"),
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z * wind_speed_z + 
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                 f(individual_local_identifier3, wind_speed_z, model = "iid") + 
                 f(lat_group, model = "rw1"),
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z * wind_speed_z + days_since_tagging_z +
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                 f(individual_local_identifier3, wind_speed_z, model = "iid"),
               data = circling_data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

# m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z * wind_speed_z + life_stage +
#                  f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
#                  f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
#                  f(individual_local_identifier3, wind_speed_z, model = "iid"),
#                data = data, family = "binomial",
#                control.compute = list(cpo = TRUE),
#                control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#### model evaluation -----------------------

#look at residuals
residuals <- data$laterality_bi - m_inla$summary.fitted.values$mean
plot(residuals)

#calculate variance explained: McFadden's pseudo R²
log_likelihood_full <- m_inla$mlik[1, 1]

null_model <- inla(laterality_bi ~ 1, family = "binomial",
                   control.compute = list(cpo = TRUE),
                   control.predictor = list(link = 1, compute = TRUE),
                   data = data)

log_likelihood_null <- null_model$mlik[1, 1]

pseudo_r_squared <- 1 - (log_likelihood_full / log_likelihood_null) #0.006

#model validation metrics
eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), # 0.52
                   Mlik = as.numeric(m_inla$mlik[,1][2])) # -13888

#### coefficients plot -----------------------------------------------------------------------------

# posterior means of coefficients
graph <- as.data.frame(summary(m_inla)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)

#remove weeks since dispersal
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

#graph$Factor_n <- as.numeric(graph$Factor)

#plot the coefficients
X11(width = 7, height = 3)
(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#39568CFF", size = 2)  +
    labs(x = "Estimate", y = "") +
    #scale_y_discrete(labels = rev(c("Intercept", "Average pitch", "Absolute cumulative yaw"))) +
    #scale_y_discrete(labels = rev(c("Intercept", "Average pitch", "Absolute cumulative yaw",
    #                                "Wind speed", "Average pitch: Absolute cumulative yaw", "Average pitch: Wind speed",
    #                                "Average cumulative yaw: Wind speed", "Average pitch : Average cumulative yaw: \nWind speed "))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#39568CFF", linewidth = 0.5) +
    ggtitle("a") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

ggsave(plot = coefs, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_coeffs_4min.pdf", 
       device = "pdf", width = 7, height = 4, dpi = 600)

#### ind_specific coefficients plot -----------------------------------------------------------------------------

#extract handedness for each individual. to use for coloring
handedness <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir)

#extract random effects,  ID is for the unique individuals
#mean_pitch_mean_z
random_effects_pitch <- m_inla$summary.random$individual_local_identifier

#abs_cum_yaw_z
random_effects_yaw <- m_inla$summary.random$individual_local_identifier2

#wind_speed_z
random_effects_wind <- m_inla$summary.random$individual_local_identifier3

#extract unique individual IDs from original data
ind_IDs <- unique(data$individual_local_identifier)

pitch <- random_effects_pitch %>% 
  mutate(coef = mean + graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Upper),
         variable = "mean_pitch_mean_z")

yaw <- random_effects_yaw %>% 
  mutate(coef = mean + graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Upper),
         variable = "abs_cum_yaw_z")

wind <- random_effects_wind %>% 
  mutate(coef = mean + graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "wind_speed_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "wind_speed_z") %>% pull(Upper),
         variable = "wind_speed_z")

three_vars <- bind_rows(pitch, yaw, wind) %>%
  full_join(handedness, by = c("ID" = "individual_local_identifier")) %>% 
  mutate(v_line = case_when(
    variable == "mean_pitch_mean_z" ~ graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate),
    variable == "abs_cum_yaw_z" ~ graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate),
    variable == "wind_speed_z" ~ graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate)
  ))


#re-order the individuals based on the laterality index
three_vars$ID <- reorder(three_vars$ID, desc(three_vars$laterality_dir))

X11(width = 9, height = 5)
(coefs_inds <- ggplot(three_vars, aes(x = coef, y = ID, color = laterality_dir)) +
    geom_vline(data = filter(three_vars, variable == "mean_pitch_mean_z"), 
               aes(xintercept = graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75",size = 0.5) + 
    geom_vline(data = filter(three_vars, variable == "abs_cum_yaw_z"), 
               aes(xintercept = graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", size = 0.5) +  
    geom_vline(data = filter(three_vars, variable == "wind_speed_z"), 
               aes(xintercept = graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", size = 0.5) +  
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed"),
                       name = "Handedness") +
    scale_y_discrete(labels = ind_IDs) +
    labs(x = "Estimate", y = "Individual ID") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) + #increase distance between x-axis values and title
    facet_wrap(~ variable, scales = "free_x", labeller = as_labeller(c(
      "mean_pitch_mean_z" = "Average pitch",
      "abs_cum_yaw_z" = "Absolute cumulative yaw",
      "wind_speed_z" = "Wind speed"
    )) # Separate panels for each variable
    ))

ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_coeffs_inds_4min.pdf", 
       device = "pdf", width = 9, height = 5, dpi = 600)


#### plot the smooth term -----------------------------------------------------------------------------

# Extract the summary of the smooth term 
smooth_effects <- m_inla$summary.random$age_group %>% 
  #back transform the age values
  mutate(age = ID * attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:scale') +
           attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:center') )

# Plot the smooth term.... HIGHLIGHT MIGRATION period
X11(width = 3.42, height = 3)
(s <- ggplot(smooth_effects, aes(x = age, y = mean)) +
    geom_line(color = "#39568CFF", linewidth = 0.3) +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "#39568CFF", alpha = 0.12) +
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    labs(y = "Effect Size", x = "Days since tagging") +
    ggtitle("b") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))
)


#combine the two plots for the linear and the non-linear terms
X11(width = 9, height = 3)
model_output_p <- grid.arrange(coefs, s, nrow = 1, widths = c(0.65, 0.35))

ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_output_5min.pdf", 
       device = "pdf", width = 9, height = 3, dpi = 600)


#### plot the interaction terms -----------------------------------------------------------------------------

## yaw:pitch----------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "pitch_yaw")

preds <- data.frame(yaw = data[na_rows,"abs_cum_yaw"],
                    pitch = data[na_rows,"mean_pitch_mean"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 7, height = 2)
pred_py <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                       values = c(0, 0.5, 0.6, 0.7, 1),
                       limits = c(0, 1),
                       na.value = "white",
                       name = "Probability of Laterality") +
  guides(fill = guide_colourbar(title.vjust = .95)) + # the legend title needs to move up a bit
  labs(x = "Absolute cumulative yaw", y = "Average pitch") +
  ggtitle("a") +
  theme_classic() +
  theme(plot.margin = margin(0, 15, 0, 0, "pt"),
        text = element_text(size = 11),
        legend.direction="horizontal",
        legend.position = "bottom",
        legend.key.width=unit(.7,"cm"),
        legend.key.height=unit(.25,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 5))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))


## wind:yaw----------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "wind_yaw")

preds <- data.frame(yaw = data[na_rows,"abs_cum_yaw"],
                    wind = data[na_rows,"wind_speed"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 7, height = 2)
pred_wy <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                       values = c(0, 0.5, 0.6, 0.7, 1),
                       limits = c(0, 1),
                       na.value = "white",
                       name = "Probability of Laterality") +
  guides(fill = guide_colourbar(title.vjust = .95)) + # the legend title needs to move up a bit
  labs(x = "Absolute cumulative yaw", y = "Wind speed") +
  ggtitle("b") +
  theme_classic() +
  theme(plot.margin = margin(0, 15, 0, 0, "pt"),
        text = element_text(size = 11),
        legend.direction="horizontal",
        legend.position = "bottom",
        legend.key.width=unit(.7,"cm"),
        legend.key.height=unit(.25,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 5))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))


##wind:pitch----------------------------------------------------------------

#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "wind_pitch")

preds <- data.frame(pitch = data[na_rows,"mean_pitch_mean"],
                    wind = data[na_rows,"wind_speed"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 7, height = 2)
pred_wp <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                       values = c(0, 0.5, 0.6, 0.7, 1),
                       limits = c(0, 1),
                       na.value = "white",
                       name = "Probability of Laterality") +
  guides(fill = guide_colourbar(title.vjust = .95)) + # the legend title needs to move up a bit
  labs(x = "Average pitch", y = "Wind speed") +
  ggtitle("c") +
  theme_classic() +
  theme(plot.margin = margin(0, 15, 0, 0, "pt"),
        text = element_text(size = 11),
        legend.direction="horizontal",
        legend.position = "bottom",
        legend.key.width=unit(.7,"cm"),
        legend.key.height=unit(.25,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 5))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))


#combine the three plots for the coefficients and the interaction term--------------
X11(width = 4.5, height = 6.5)
combined <- pred_py + pred_wy + pred_wp & theme(legend.position = "bottom")
(p <- combined + plot_layout(guides = "collect", nrow = 3))

ggsave(plot = p, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_interactions.pdf", 
       device = "pdf", width = 4.5, height = 6.5, dpi = 600)

#-----------------------------------------------------------------------------------------------------------------------
## Step 5.1: Does laterality help with better performance when individuals are not experienced? Flight performance #####
#-----------------------------------------------------------------------------------------------------------------------
#flight performance ~ level of handedness * age

#open data
circling_data <- readRDS("thinned_laterality_for_modeling.rds")

#I have mean of sd of yaw, pitch, and roll. consider calculating the maxes for this analysis.... will need to go back to the original
#code where I summarized the Quat for each 8-sec burst

#### ----------------------- look at multi-collinearity
circling_data %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
                  "abs_cum_yaw", "days_since_tagging", "wind_speed")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #sd of roll and pitch = 0.65; sd of yaw and cumulative yaw = 0.7

#model each flight performance separately, but because pitch and roll are correlated, maybe only do pitch and yaw? ...
#also because laterality is calculated using roll anyway.... maybe also control for thermal strength


#### ----------------------- model: regression with inla

set.seed(777)
#create new data to predict the interaction between age and laterality
#create a new dataset to use for making predictions for interaction of mean pitch and cumulative yaw
new_data <- expand.grid(wind_speed = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                        laterality_bi = c(0, 1), #laterality values
                        days_since_tagging = seq(from = min(circling_data$days_since_tagging, na.rm = T), to = max(circling_data$days_since_tagging, na.rm = T),  length.out = 50)) %>% 
  mutate(#set the dependent variables to NA. keep both the z-transormed and original columns. The original will be used for log-transformation of the response
    mean_yaw_sd_z = NA, 
    mean_pitch_sd_z = NA, 
    mean_yaw_sd = NA, 
    mean_pitch_sd = NA, 
    #calculate z scores
    wind_speed_z = (wind_speed - attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:scale'),
    days_since_tagging_z = (days_since_tagging - attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:scale'),
    #randomly assign individual IDs
    individual_local_identifier = sample(circling_data$individual_local_identifier, nrow(.), replace = T))

#bin the age variable and append the new data
data <- circling_data %>% 
  #only keep columns that exist in the new_data
  dplyr::select(intersect(colnames(new_data), colnames(.))) %>% 
  full_join(new_data) %>% 
  mutate(individual_local_identifier2 = individual_local_identifier, #repeat individual ID column to be used in the model formula for random effects
         individual_local_identifier3 = individual_local_identifier,
         age_group = inla.group(days_since_tagging_z, n = 100, method = "quantile"),
         laterality_bi = as.factor(laterality_bi),
         individual_local_identifier = as.factor(individual_local_identifier))

saveRDS(data, file = "data_flight_performance_models.rds")

#### run the models one at a time! for pitch and yaw separately

#long-transform instead of z-transform, to make sure response values are positive
m_inla_p <- inla(log(mean_pitch_sd) ~ 1 + laterality_bi * age_group * wind_speed_z +
                   f(individual_local_identifier, age_group, model = "iid"),
                 data = data,
                 control.compute = list(cpo = TRUE),
                 control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

### the predictions are too large.... use the gamma family to restrict the predictions to be positive, but not too large
m_inla_y <- inla(log(mean_yaw_sd) ~ 1 + laterality_bi * age_group * wind_speed_z +
                   f(individual_local_identifier, age_group, model = "iid"),
                 data = data,
                 #family = "gamma",
                 control.compute = list(cpo = TRUE),
                 control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted


####  model evaluation -----------------------

#model validation metrics
data.frame(CPO = mean(m_inla_p$cpo$cpo, na.rm = T), # 0.49 with log transformation
           Mlik = as.numeric(m_inla_p$mlik[,1][2])) # -17933.39

data.frame(CPO = mean(m_inla_y$cpo$cpo, na.rm = T), # 0.32 #improves when predicting the log-transformed variable
           Mlik = as.numeric(m_inla_y$mlik[,1][2])) # -24833.32

#### pitch_model: coefficients plot -----------------------------------------------------------------------------

# posterior means of coefficients
graph_p <- as.data.frame(summary(m_inla_p)$fixed)
colnames(graph_p)[which(colnames(graph_p)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph_p)[which(colnames(graph_p)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph_p)[which(colnames(graph_p)%in%c("mean"))]<-c("Estimate")

graph_p$Factor <- rownames(graph_p)

#remove weeks since dispersal
VarOrder <- rev(unique(graph_p$Factor))
VarNames <- VarOrder

graph_p$Factor <- factor(graph_p$Factor, levels = VarOrder)
levels(graph_p$Factor) <- VarNames


#after log-transforming the response, the intercept becomes too big....  remove it for easier visualization. report it in the caption of the figure
graph_p <- graph_p[-1,]

#plot the coefficients
X11(width = 7, height = 3)
(coefs_p <- ggplot(graph_p, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#8a2be2ff", size = 2)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Laterality", "Days since tagging", 
                                    "Wind speed", "Laterality: Days since tagging", "Laterality: Wind speed",
                                    "Days since tagging: Wind speed", "Laterality: Days since tagging:\nWind speed "))) + 
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#8a2be2ff", linewidth = 0.5) +
    theme_classic() +
    ggtitle("a") +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

ggsave(plot = coefs_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_pitch_model_ws_coeffs.pdf", 
       device = "pdf", width = 7, height = 3, dpi = 600)


## conclusions: horizontal wobble (sd_yaw) goes down with age and with laterality,
## vertical wobble (sd_pitch) goes up with age and laterality, BUT there is a negative interaction term
## TO DO: 1) plot the interaction term. 2) figure out if the response should be transformed or not. 3) add ind level
## random slope on age and plot ind-specific differences (but this isn't really the focus here... so maybe just add
## inds as a random effect on the slope, just to control for indivdiual variation. )


#### pitch_model: ind_specific coefficients plot -----------------------------------------------------------------------------

#extract handedness for each individual. to use for coloring
handedness <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir)

#extract random effects,  ID is for the unique individuals
#age
random_effects_age <- m_inla_p$summary.random$individual_local_identifier

#extract unique individual IDs from original data
ind_IDs <- unique(data$individual_local_identifier)

age <- random_effects_age %>% 
  mutate(coef = mean + graph_p %>% filter(Factor == "age_group") %>% pull(Estimate),
         lower = .[,4] +  graph_p %>% filter(Factor == "age_group") %>% pull(Lower),
         upper = .[,6] +  graph_p %>% filter(Factor == "age_group") %>% pull(Upper),
         variable = "age_group") %>% 
  full_join(handedness, by = c("ID" = "individual_local_identifier")) %>% 
  mutate(v_line =  graph_p %>% filter(Factor == "age_group") %>% pull(Estimate))

#re-order the individuals based on the laterality index
age$ID <- reorder(age$ID, desc(age$laterality_dir))

X11(width = 7, height = 5)
(coefs_inds_p <- ggplot(age, aes(x = coef, y = ID, color = laterality_dir)) +
    geom_vline(data = age, 
               aes(xintercept = graph_p %>% filter(Factor == "age_group") %>% pull(Estimate)), linetype = "dashed", color = "gray75",size = 0.5) + 
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#20A387FF", "right_handed" = "#B8DE29FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed"),
                       name = "Handedness") +
    scale_y_discrete(labels = ind_IDs) +
    labs(x = "Estimate", y = "") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))  #increase distance between x-axis values and title
)

ggsave(plot = coefs_inds_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_pitch_model_ws_coeffs_inds.pdf", 
       device = "pdf", width = 7, height = 5, dpi = 600)


#### pitch_model: interaction plot -----------------------------------------------------------------------------

#extract information for rows that had NAs as response variables
na_rows <- which(is.na(data$mean_pitch_sd_z))

preds_p <- data.frame(days_since_tagging = data[na_rows,"days_since_tagging"],
                      wind = data[na_rows,"wind_speed"],
                      preds = m_inla_p$summary.fitted.values[na_rows,"mean"]) %>% 
  #back-transform the log values
  mutate(preds = exp(preds)) %>% 
  #back-transform the z-scores back to sd_pitch values
  #mutate(preds = preds * attr(circling_data[,colnames(circling_data) == "mean_pitch_sd_z"],'scaled:center') +
  #         attr(circling_data[,colnames(circling_data) == "mean_pitch_sd_z"],'scaled:scale')) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 7, height = 2)
pred_dw <- preds_p %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                       values = c(0, 0.5, 0.7, 4.7),
                       limits = c(0, 5),
                       na.value = "white",
                       name = "Pitch \nwobble") +
  guides(fill = guide_colourbar(title.vjust = .95)) + # the legend title needs to move up a bit
  labs(x = "Days since tagging", y = "Wind speed") +
  ggtitle("b") +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        text = element_text(size = 11),
        legend.key.width=unit(.25,"cm"),
        legend.key.height=unit(.4,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 9),
        legend.margin = margin(0, 0, 0, 0), # Reduce margin around the legend
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 5))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))


#combine the two plots for the coefficients and the interaction term
X11(width = 7, height = 2.3)
model_output_p <- grid.arrange(coefs_p, pred_dw, nrow = 1, widths = c(0.6, 0.4))

ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_pitch_model_ws_output.pdf", 
       device = "pdf", width = 7, height = 2.3, dpi = 600)



#### yaw_model: coefficients plot -----------------------------------------------------------------------------

# posterior means of coefficients
graph_y <- as.data.frame(summary(m_inla_y)$fixed)
colnames(graph_y)[which(colnames(graph_y)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph_y)[which(colnames(graph_y)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph_y)[which(colnames(graph_y)%in%c("mean"))]<-c("Estimate")

graph_y$Factor <- rownames(graph_y)

#remove weeks since dispersal
VarOrder <- rev(unique(graph_y$Factor))
VarNames <- VarOrder

graph_y$Factor <- factor(graph_y$Factor, levels = VarOrder)
levels(graph_y$Factor) <- VarNames


#after log-transforming the response, the intercept becomes too big....  remove it for easier visualization. report it in the caption of the figure
graph_y <- graph_y[-1,]

#plot the coefficients
X11(width = 7, height = 3)
(coefs_y <- ggplot(graph_y, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#8a2be2ff", size = 2)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Laterality", "Days since tagging", 
                                    "Wind speed", "Laterality: Days since tagging", "Laterality: Wind speed",
                                    "Days since tagging: Wind speed", "Laterality: Days since tagging:\nWind speed "))) + 
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#8a2be2ff", linewidth = 0.5) +
    theme_classic() +
    ggtitle("a") +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

ggsave(plot = coefs_y, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_ws_coeffs.pdf", 
       device = "pdf", width = 7, height = 3, dpi = 600)

## conclusions: horizontal wobble (sd_yaw) goes down with age and with laterality,
## vertical wobble (sd_pitch) goes up with age and laterality, BUT there is a negative interaction term
## TO DO: 1) plot the interaction term. 2) figure out if the response should be transformed or not. 3) add ind level
## random slope on age and plot ind-specific differences (but this isn't really the focus here... so maybe just add
## inds as a random effect on the slope, just to control for indivdiual variation. )


#### yaw_model: ind_specific coefficients plot -----------------------------------------------------------------------------

#extract handedness for each individual. to use for coloring
handedness <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir)

#extract random effects,  ID is for the unique individuals
#age
random_effects_age <- m_inla_y$summary.random$individual_local_identifier

#extract unique individual IDs from original data
ind_IDs <- unique(data$individual_local_identifier)

age <- random_effects_age %>% 
  mutate(coef = mean + graph_y %>% filter(Factor == "age_group") %>% pull(Estimate),
         lower = .[,4] +  graph_y %>% filter(Factor == "age_group") %>% pull(Lower),
         upper = .[,6] +  graph_y %>% filter(Factor == "age_group") %>% pull(Upper),
         variable = "age_group") %>% 
  full_join(handedness, by = c("ID" = "individual_local_identifier")) %>% 
  mutate(v_line =  graph_y %>% filter(Factor == "age_group") %>% pull(Estimate))

#re-order the individuals based on the laterality index
age$ID <- reorder(age$ID, desc(age$laterality_dir))

X11(width = 7, height = 5)
(coefs_inds_y <- ggplot(age, aes(x = coef, y = ID, color = laterality_dir)) +
    geom_vline(data = age, 
               aes(xintercept = graph_y %>% filter(Factor == "age_group") %>% pull(Estimate)), linetype = "dashed", color = "gray75",size = 0.5) + 
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#20A387FF", "right_handed" = "#B8DE29FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed"),
                       name = "Handedness") +
    scale_y_discrete(labels = ind_IDs) +
    labs(x = "Estimate", y = "") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))  #increase distance between x-axis values and title
)

ggsave(plot = coefs_inds_y, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_ws_coeffs_inds.pdf", 
       device = "pdf", width = 7, height = 5, dpi = 600)

#### yaw_model: interaction plot -----------------------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(is.na(data$mean_yaw_sd_z))

preds_y <- data.frame(days_since_tagging = data[na_rows,"days_since_tagging"],
                      wind = data[na_rows,"wind_speed"],
                      preds = m_inla_y$summary.fitted.values[na_rows,"mean"]) %>% 
  #back-transform the log values
  mutate(preds = exp(preds)) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 7, height = 2)
(pred_dw <- preds_y %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = preds)) +
    scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                         values = c(0, 3.5, 4, 20, 270),
                         limits = c(0, 5),
                         na.value = "white",
                         name = "Yaw \nwobble") +
    guides(fill = guide_colourbar(title.vjust = .95)) + # the legend title needs to move up a bit
    labs(x = "Days since tagging", y = "Wind speed") +
    ggtitle("b") +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          text = element_text(size = 11),
          legend.key.width=unit(.25,"cm"),
          legend.key.height=unit(.4,"cm"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 9),
          legend.margin = margin(0, 0, 0, 0), # Reduce margin around the legend
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) + # increase distance between x-axis values and title
    scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
    scale_y_continuous(expand = c(0, 0))
)

#combine the two plots for the coefficients and the interaction term
X11(width = 7, height = 2.3)
model_output_y <- grid.arrange(coefs_y, pred_dw, nrow = 1, widths = c(0.6, 0.4))

ggsave(plot = model_output_y, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_ws_output.pdf", 
       device = "pdf", width = 7, height = 2.3, dpi = 600)

#-----------------------------------------------------------------------------------------------------------------------
## Step 5.2: Does laterality help with better performance when individuals are not experienced? migration performance #####
#-----------------------------------------------------------------------------------------------------------------------

### Overall migration performance -----------------------------------------------------------------------------
#migration performance ~ level of handedness  (one row per individual)

#extract handedness for each individual. to use for coloring
handedness <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir) %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier))

#data from Ellen
migr_dates <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/HB_dates_Fin_all_M2.rds")
migr_info <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/mig_sum.rds") %>% 
  full_join(handedness) %>% 
  drop_na(individual_local_identifier) %>% 
  mutate(laterality_bi = ifelse(laterality_dir == "ambidextrous", 0, 1) %>% as.factor()) #create a binary variable for handedness vs not)


#exploratory plots
ggplot(migr_info, aes(x = laterality_dir, y = avg_speed_kmh)) +
  geom_boxplot() +
  geom_point() +
  labs(x = "Laterality", y = "avg_speed_kmh") +
  theme_classic()

ggplot(migr_info, aes(x = laterality_bi, y = avg_speed_kmh)) +
  geom_boxplot() +
  geom_point() +
  labs(x = "Laterality", y = "avg_speed_kmh") +
  theme_classic()


### Daily migration performance -----------------------------------------------------------------------------

#extract handedness for each individual
handedness <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir) %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier))

#open data from Ellen and append handedness
migr_d <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/daily_mig_metrics_wobble.rds") %>% 
  full_join(handedness) %>% 
  drop_na(individual_local_identifier) %>% 
  rename(Overall_ind_laterality_dir = laterality_dir) %>% 
  mutate(daily_mean_vert_speed = as.numeric(daily_mean_vert_speed),
         Overall_ind_laterality_bi = ifelse(Overall_ind_laterality_dir == "ambidextrous", 0, 1) %>% as.factor(), #create a binary variable for handedness vs not)
         ind_day = paste0(individual_local_identifier, "_", date)) %>% 
  as.data.frame()


#exploratory plots
migr_long <- migr_d %>%
  pivot_longer(cols = c(daily_vedba, daily_distance_km, daily_speed_kmh, 
                        daily_mean_vert_speed, daily_sd_vert_speed, 
                        daily_sd_turn_angle, daily_mean_turn_angle),
               names_to = "measure", values_to = "value")

ggplot(migr_long, aes(x = laterality_dir, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 0.2) +
  labs(x = "Laterality", y = "Value") +
  theme_classic() +
  facet_wrap(~ measure, scales = "free_y")

ggplot(migr_long, aes(x = laterality_bi, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 0.2) +
  labs(x = "Laterality", y = "Value") +
  theme_classic() +
  facet_wrap(~ measure, scales = "free_y")

### Annotate daily migration with handedness and wind -----------------------------------------------------

#open laterality data for circling
circling_data <- readRDS("thinned_laterality_for_modeling.rds")

#migration days
migration_days <- migr_d %>% 
  drop_na(date) %>% 
  distinct(ind_day) %>% 
  pull(ind_day)

#group by unique day, calculate max and mean wind speed, calculate the laterality index
migr_d_l_w <- circling_data %>% 
  mutate(day = date(start_timestamp),
         ind_day = paste0(individual_local_identifier, "_", day)) %>% 
  #subset for days in the migration data 
  filter(ind_day %in% migration_days) %>% 
  group_by(individual_local_identifier, day) %>% 
  summarize(daily_max_wind = max(wind_speed, na.rm = T),
            daily_mean_wind = mean(wind_speed, na.rm = T),
            daily_max_cum_yaw = max(abs_cum_yaw, na.rm = T),
            daily_mean_cum_yaw = mean(abs_cum_yaw, na.rm = T),
            daily_max_mean_pitch = max(mean_pitch_mean, na.rm = T),
            daily_mean_mean_pitch = mean(mean_pitch_mean, na.rm = T),
            #calculate laterality for this day ...
            daily_sum_bank_left = sum(bank_left),
            daily_sum_bank_right = sum(bank_right),
            daily_sum_bank_straight = sum(bank_straight),
            daily_laterality_index = (daily_sum_bank_right - daily_sum_bank_left)/(daily_sum_bank_right + daily_sum_bank_left),
            mode_laterality = getmode(laterality_dir),
            ind_day = head(ind_day, 1), #could have used this for grouping too
            .groups = "drop") %>% 
  mutate(daily_laterality_dir = case_when(
    between(daily_laterality_index, 0.25, 1.0) ~ "right_handed",
    between(daily_laterality_index, -1.0, -0.25) ~ "left_handed",
    between(daily_laterality_index, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  #bind to the daily migration data
  left_join(migr_d %>% drop_na(date)) %>% 
  as.data.frame() 

saveRDS(migr_d_l_w, file = "data_migration_performance_models.rds")


#exploratory plots
migr_long <- migr_d_l_w %>%
  pivot_longer(cols = c(daily_vedba, daily_distance_km, daily_speed_kmh, 
                        daily_mean_vert_speed, daily_sd_vert_speed, 
                        daily_sd_turn_angle, daily_mean_turn_angle, daily_max_wind, daily_mean_wind),
               names_to = "measure", values_to = "value")

ggplot(migr_long, aes(x = mode_laterality, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 0.2) +
  labs(x = "Laterality", y = "Value") +
  theme_classic() +
  facet_wrap(~ measure, scales = "free_y")


### Model migration performance as a function of laterality -----------------------------------------------------
# 
# ggplot(migr_d_l_w, aes(x = daily_vedba, y = daily_max_wind, color = mode_laterality)) +
#   geom_point() + 
#   geom_smooth(method = "loess")

#z-transform the variables
data <- migr_d_l_w %>% 
  #make sure to do all the filtering before scaling the variables!!!!
  mutate_at(c(3:6, 17:23),
            list(z = ~scale(.))) %>% 
  as.data.frame()

#check for autocorrelation
data %>% 
  dplyr::select(c(26:36)) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #correlated: mean and max wind, max cum yaw and mean cum yaw, daily distance an daily speed

#model
m_inla <- inla(daily_vedba ~ 1 + mode_laterality + daily_max_wind +
                 f(individual_local_identifier, model = "iid"), #put the random effect only on the slope
               data = data,
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#model distance instead of speed. they are correlated, but the distance model performs better. with no transformation of the predictor
m_inla <- inla(daily_distance_km ~ 1 + mode_laterality + daily_max_wind + 
                 f(individual_local_identifier, model = "iid"),
               data = data,
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

m_inla <- inla(daily_mean_cum_yaw ~ 1 + mode_laterality + daily_max_wind +
                 f(individual_local_identifier, model = "iid"),
               data = data,
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

m_inla <- inla(daily_mean_mean_pitch ~ 1 + mode_laterality + daily_max_wind +
                 f(individual_local_identifier,  model = "iid"),
               data = data,
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted


# vertical speed was calculated using the gps.... so think about whether it should be included or not
m_inla <- inla(log(daily_mean_vert_speed) ~ 1 + mode_laterality + daily_max_wind_z +
                 f(individual_local_identifier, daily_max_wind_z, model = "iid"),
               data = data,
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted




#### coefficients plot -----------------------------------------------------------------------------

# posterior means of coefficients
graph <- as.data.frame(summary(m_inla)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)

#remove weeks since dispersal
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

#graph$Factor_n <- as.numeric(graph$Factor)

#plot the coefficients
X11(width = 4.5, height = 1.6)
(coefs_y <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#8a2be2ff", size = 2)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Intercept", "Laterality left-handed", "Laterality rigt-handed",
                                    "Maximum wind speed"))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#8a2be2ff", linewidth = 0.5) +
    #ggtitle("Daily distance (km)") +
    #ggtitle("VeDBA (g)") +
    #ggtitle("Average pitch (degrees)") +
    ggtitle("Average cumulative yaw (degrees)") +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 10, 0, "pt"),
          text = element_text(size = 9),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

#put all plots together
X11(width = 4.5, height = 6.5)
combined <- coefs_v + coefs_d + coefs_y + coefs_p 
(p <- combined + plot_layout(axis_titles = "collect", nrow = 4))

ggsave(plot = p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/migration_model_coeffs.pdf", 
       device = "pdf", width = 4.5, height = 6.5, dpi = 600)


############################### plotting efforts & proof of laterality ################################################################


#---------------------------------------------------------------------------------
## Step 2: Repeatability over time                                           #####
#---------------------------------------------------------------------------------

#FILTER: remove pre-fledging
post_fledging <- or_w_gps_df_LS %>% 
  filter(life_stage != "pre-fledging")

lapply(split(post_fledging, post_fledging$individual_local_identifier), function(ind){
  
  #repeatability for the whole dataset
  rpt(laterality_bank ~ (1 | days_since_tagging), grname = "days_since_tagging", data = ind, 
      datatype = "Gaussian", nboot = 0, npermut = 0) #
  
})

rpt(laterality_bank ~ (1 | individual_local_identifier), grname = "individual_local_identifier", data = post_fledging, datatype = "Gaussian", 
    nboot = 0, npermut = 0) #no repeatability

rpt(laterality_bi ~ (1 | life_stage), grname = "life_stage", data = or_w_gps_df_LS, datatype = "Binary", 
    nboot = 0, npermut = 0) #no repeatability

rpt(cumulative_roll_8sec ~ (1 | individual_local_identifier) + (1 | life_stage), grname = c("individual_local_identifier","life_stage"),
    data = or_w_gps_df_LS,
    datatype = "Gaussian", 
    nboot = 0, npermut = 0) #no repeatability

rpt(laterality_bi ~ (1 | individual_local_identifier) + (1 | life_stage), grname = c("individual_local_identifier","life_stage"),
    data = or_w_gps_df_LS,
    datatype = "Binary", 
    nboot = 0, npermut = 0) #no repeatability

#only for one individual
rpt(laterality_bi ~  (1 | days_since_tagging), grname = "days_since_tagging",
    data = or_w_gps_df_LS %>% filter(individual_local_identifier == "D163_696" & life_stage == "wintering"),
    datatype = "Binary", 
    nboot = 0, npermut = 0) #no repeatability

#only pre-migration

rpt(laterality_bi ~ (1 | days_since_tagging), grname = "days_since_tagging", data = or_w_gps_df_LS %>% filter(life_stage == "post-fledging"), 
    datatype = "Binary", nboot = 0, npermut = 0)

rpt(laterality_bi ~ (1 | days_since_tagging), grname = "days_since_tagging", data = or_w_gps_df_LS %>% filter(life_stage == "migration"), 
    datatype = "Binary", nboot = 0, npermut = 0)

rpt(laterality_bi ~ (1 | days_since_tagging), grname = "days_since_tagging", data = or_w_gps_df_LS %>% filter(life_stage == "wintering"), 
    datatype = "Binary", nboot = 0, npermut = 0)


#try cohen's Kappa
library(irr)

# Function to calculate Cohen's Kappa for each individual
calculate_kappa <- function(individual_data) {
  # Create a table of laterality assignments
  laterality_table <- table(individual_data$days_since_tagging, individual_data$laterality_dir)
  
  # Calculate Cohen's Kappa
  kappa_result <- kappa2(laterality_table)
  
  return(kappa_result$value)
}

# Apply the function to each individual
kappa_values <- by(data, data$Individual, calculate_kappa)

#---------------------------------------------------------------------------------
## Step 2: Plot laterality for each latitudinal zone                         #####
#---------------------------------------------------------------------------------

#also add a map of the data. color migration section (use same color in step 3 to shade migration period)



#---------------------------------------------------------
## Step 3: Plot laterality over time                 #####
#---------------------------------------------------------

### DAILY LI
#data: daily laterality index, filtered for circling only.  annotated with day since tagging and life cycle stage

#open meta-data to calculate day since tagging
meta_data <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/EHB_metadata - Sheet1.csv") %>% 
  mutate(deployment_dt_utc = as.POSIXct(deployment_dt_utc, tz = "UTC")) %>% 
  filter(nest_country != "Switzerland")

#daily LI calculated in L02_QUAT_exploration.r (L 192~). already filtered for circling flight. add days since tagging.
LI_days <- readRDS("laterality_index_per_day.rds") %>%
  mutate(individual_local_identifier = as.character(individual_local_identifier),
         date = as.Date(unique_date)) %>% 
  full_join(meta_data %>% select(ring_ID, deployment_dt_utc), by = c("individual_local_identifier" = "ring_ID")) %>% 
  rowwise() %>% 
  mutate(days_since_tagging = floor(difftime(date, deployment_dt_utc, unit = "days")) %>% as.numeric()) %>% 
  ungroup() %>% 
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  filter(between(days_since_tagging, 1, 270)) %>% 
  arrange(individual_local_identifier,unique_date ) %>% 
  as.data.frame()


#Plot 2D densities of different directions of handedness

#re-order laterality direction
LI_days$laterality_dir <- factor(LI_days$laterality_dir, levels = c("left_handed", "ambidextrous", "right_handed"))


X11(width = 7, height = 2.5)
d_d <- ggplot(data = LI_days, aes(x = days_since_tagging, y = laterality_dir)) +
  stat_bin2d(aes(fill = ..density..), bins = 270/7) + #one bar per week
  scale_fill_viridis(option = "plasma") +
  labs(x = "Days since tagging", y = "") +
  scale_y_discrete(labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
  #start and end of migration on average
  geom_vline(xintercept = c(population_migr_period$avg_migration_start,population_migr_period$avg_migration_end),
             linetype = "dashed",  color = "white", linewidth = 0.5, show.legend = TRUE) +
  ggtitle("Density of handedness categories based on laterality index calculated for \neach day for each individual. \nWhite dashed lines show the migration period.") +
  theme_minimal() +
  theme(text = element_text(size = 9))

ggsave(plot = d_d, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/daily_handedness_2d.jpg", 
       device = "jpg", width = 7, height = 2.5)

### BANK ANGLES IN EACH DAY

#this one is weird
b_d <- ggplot(data = or_w_gps_df %>% filter(days_since_tagging <= 270), aes(x = days_since_tagging, y = laterality_dir)) +
  stat_bin2d(aes(fill = ..density..), bins = 270/7) + #one bar per week
  scale_fill_viridis(option = "plasma") +
  labs(x = "Days since tagging", y = "") +
  scale_y_discrete(labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
  #start and end of migration on average
  geom_vline(xintercept = c(population_migr_period$avg_migration_start,population_migr_period$avg_migration_end),
             linetype = "dashed",  color = "white", linewidth = 0.5, show.legend = TRUE) +
  theme_minimal() +
  theme(text = element_text(size = 11))

ggsave(plot = b_d, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/daily_handedness_bursts_2d.jpg", 
       device = "jpg", width = 7, height = 2.5)



## old stuff..............................................................................................................................
laterality_1sec_days <- laterality_1sec_days %>% 
  mutate(circling_status = as.factor(circling_status))

X11(width = 7, height = 2.8) 
(p <- ggplot(data = laterality_1sec_days, aes(x = mean_roll_mean, fill = circling_status)) +
    geom_histogram(aes(after_stat(density)), binwidth = 1) + #plot the probability
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_x_continuous(breaks =seq(-90, 90, by = 10), 
                       labels = seq(-90, 90, by = 10),
                       limits = c(-92, 92)) +
    scale_fill_manual(values = c("straight" = "#238A8DFF", "circling" = "#8a2be2ff", "shallow circling" = "#DCE319FF"),
                      labels = rev(c("Straight", "Shallow circling", "Circling"))) +
    labs(x = "Bank angle (deg) averaged over each 8 second burst",
         y = "Count",
         fill = "Circularity") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

ggsave(plot = p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distribution.pdf", 
       device = "pdf", width = 7, height = 2.8, dpi = 600)


##QQ plot style

# Extract data for each circling status
straight <- laterality_1sec_days$mean_roll_mean[laterality_1sec_days$circling_status == "straight"]
circling <- laterality_1sec_days$mean_roll_mean[laterality_1sec_days$circling_status == "circling"]
shallow <- laterality_1sec_days$mean_roll_mean[laterality_1sec_days$circling_status == "shallow circling"]

# Calculate ECDFs
ecdfS <- ecdf(straight)
ecdfC <- ecdf(circling)
ecdfSh <- ecdf(shallow)

# Create Q-Q data for circling and shallow circling
qqC <- data.frame(
  ecdf_straight = ecdfS(circling),
  ecdf_circling = ecdfC(circling)
)

qqSh <- data.frame(
  ecdf_straight = ecdfS(shallow),
  ecdf_shallow_circling = ecdfSh(shallow)
)

# Plot using ggplot2
X11(width = 7, height = 2.8) 
qq <- ggplot() +
  geom_line(data = qqC, aes(x = ecdf_straight, y = ecdf_circling, color = "Circling"), linetype = "dashed", linewidth = 0.5) +
  geom_line(data = qqSh, aes(x = ecdf_straight, y = ecdf_shallow_circling, color = "Shallow Circling"), linetype = "dashed", linewidth = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Circling" = "#8a2be2ff", "Shallow Circling" = "#238A8DFF")) +
  labs(x = "ECDF Straight", y = "ECDF Comparison", 
       #title = "Q-Q Plot of Circling and Shallow Circling vs. Straight",
       color = "Circularity") +
  theme_classic() +
  theme(text = element_text(size = 11),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 9),
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 5)))

ggsave(plot = qq, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distribution_qq.pdf", 
       device = "pdf", width = 7, height = 2.8, dpi = 600)


#-------------------------------------------------------
## Step 3.2: Is there laterality? Individual-level #####
#-------------------------------------------------------

#calculate laterality index for each individual. overall.
#based on Bennison et al: 
#right-side bias (LI = 1.0 to 0.25), a left side bias (LI = −1 to −0.25), no bias (LI = −0.25 to 0.25)

#large circles only
circling_only <- laterality_1sec_days %>% 
  filter(circling_status == "circling") %>% #only keep thermaling flight (>45 deg)
  select(-c("n_records", "bank_left", "bank_right", "bank_straight", "heading_left", #remove laterality columns for the burst-specific indices, because i'll calculate this for each individual
            "heading_right", "heading_straight", "laterality_bank", "laterality_heading"))

#calculate overall laterality index for each individual
ind_laterlaity_index <- circling_only %>%
  mutate(bank_direction = ifelse(mean_roll_mean < 0, "left",
                                 ifelse(mean_roll_mean > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_8sec < 0, "left",
                                    ifelse(cumulative_yaw_8sec > 0, "right", "straight"))) %>% 
  group_by(individual_local_identifier) %>% 
  summarise(total_n = n(),
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            heading_left = sum(heading_direction == "left"),
            heading_right = sum(heading_direction == "right"),
            heading_straight = sum(heading_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
            laterality_heading = (heading_right - heading_left)/(heading_right + heading_left),
            .groups = 'drop') %>% 
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  as.data.frame()

#add laterality as a new column to the data
circling_w_LI <- circling_only %>% 
  full_join(ind_laterlaity_index, by = "individual_local_identifier") %>% 
  mutate(laterality_dir = factor(laterality_dir)) %>% 
  mutate(laterality_dir = fct_relevel(laterality_dir, c("left_handed", "ambidextrous", "right_handed"))) 


saveRDS(circling_w_LI, file = "circling_w_LI_population.rds")

#re-order the individuals based on the laterality index
circling_w_LI$individual_local_identifier <- reorder(circling_w_LI$individual_local_identifier, desc(circling_w_LI$laterality_dir))

X11(width = 7, height = 7)
(p_inds <- ggplot(data = circling_w_LI, aes(x = mean_roll_mean, y = individual_local_identifier, 
                                            color = laterality_dir, fill = laterality_dir),
                  height = stat(density)) +
    geom_density_ridges(stat = "binline", bins = 150, scale = 0.98, alpha = 0.8, draw_baseline = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_x_continuous(breaks = seq(-90, 90, by = 10), 
                       labels = seq(-90, 90, by = 10),
                       limits = c(-92, 92)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.03))) +  # Add extra space above the highest level
    scale_fill_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    labs(x = "Bank angle (deg) averaged over each 8 second burst",
         y = "Individual ID",
         fill = "Handedness",
         color = "Handedness") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)
ggsave(plot = p_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distr_circling_inds.pdf", 
       device = "pdf", width = 7, height = 7, dpi = 600)


