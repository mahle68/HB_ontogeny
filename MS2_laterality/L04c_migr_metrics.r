# script for calculating migration metrics for Safi et al 2025.
# Elham Nouani, PhD. 
# 05.11.2024, Konstanz, DE

library(move2)
library(sf)

#open migration details to add life stage to gps data
#life-cycle stages from L03a_tests_per_day.r
life_cycle <- readRDS("updated_life_cycle_nov24.rds")

#open gps data. make a unique ID for each ind-day combo. for speed, etc. calculations.
gps_no_winter <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24.rds") %>%  #2022-09-25 07:56:41.0000 to 2024-04-15 10:02:33.0000
  drop_na(individual_local_identifier, location_lat) %>% #remove NA individuals and NA locations.
  mutate(unique_date = as.Date(timestamp),
         yr = year(timestamp),
         mn = month(timestamp),
         dy = day(timestamp),
         hr = hour(timestamp),
         unique_hr = paste(yr,mn,dy,hr, sep = "_"),
         closest_hr = round(timestamp, units = "hours") %>% as.character()) %>% 
  mutate(ind_day = paste0(individual_local_identifier, "_", unique_date)) %>% 
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
  filter(!(life_stage %in% c("wintering", "pre-fledging"))) %>% 
  as.data.frame()

saveRDS(gps_no_winter, file = "gps_data_LS_no_winter.rds")


#filter the data for 15 minutes. then add burs IDs for consecutive points 15 min apart and calculate distances/speeds for only these bursts.
#define a thinning threshold
thinning_th <- 15*60

gps_15min <- gps_no_winter %>% 
  st_as_sf(coords = c("location_long", "location_lat"), crs = "EPSG:4326") %>% 
  mutate(dt_15min = round_date(timestamp, "15 minutes")) %>% 
  group_by(individual_local_identifier, unique_date, dt_15min) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  arrange(timestamp, .by_group = T) %>% 
  mutate(time_lag_sec = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "secs") %>% as.numeric()),
         burst_id = cumsum(time_lag_sec > thinning_th)) %>%  #increase burst id by one, every time time_lag is more than thinning_th
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date, burst_id) %>% 
  filter(n() >= 2, .preserve = T) %>% # Remove groups smaller than 2
  mutate(step_length = if_else(row_number() == 1, 0, st_distance(geometry, lag(geometry), by_element = TRUE) %>% as.numeric())) %>% #meters
  ungroup() %>% 
  as.data.frame()

#calculate speed using the step lengths..... right?

################# hourly attempt

gps_1hr <- gps_no_winter %>% 
  filter(life_stage == "migration") %>%  
  st_as_sf(coords = c("location_long", "location_lat"), crs = "EPSG:4326") %>% 
  mutate(dt_1hr = round_date(timestamp, "1 hour")) %>% 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  #summarize(n_hrs_in_day = n()) #most days have 12-14 hours of data
  filter(n() >= 12, .preserve = T) #keep days with >= 12 hours of data

#how many days per ind remains?
gps_1hr %>% 
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  slice(1) %>% 
  ungroup() %>% 
  group_by(individual_local_identifier) %>% 
  summarize(n()) %>% 
  as.data.frame() #all individuals have a good number of days left (20-30)

