#environmental annotation of imu data being used for the laterality paper
#sub-step in L04_full_workflow.r 
#14.10.2024 Konstanz, Germany
#Elham Nourani. enourani@ab.mpg.de


#open data filtered and matched with GPS in L04_full_workflow_r

or_w_gps_df <- readRDS("thinned_laterality_w_gps.rds") %>% 
  drop_na(location_long_closest_gps_raw) %>% 
  mutate(row_id = row_number(),
         #create grouping factor to split the df into 5 groups for annotation request based on latitude and time (each request can't be bigger than 10GB)
         lat_group = case_when(
           location_lat_closest_gps_raw >= 45 ~ "one",
           between(location_lat_closest_gps_raw, 27, 45) ~ "two",
           between(location_lat_closest_gps_raw, 9, 27) ~ "three",
           between(location_lat_closest_gps_raw, -9, 9) ~ "four",
           location_lat_closest_gps_raw <= -9 ~ "five"
         ))
#create grouping factor to split the df into 5 groups for annotation request (each request can't be bigger than 10GB)
#annotation_group = c(rep(1:4, each = ceiling(nrow(or_w_gps_df)/5)), rep(5, nrow(or_w_gps_df) - ceiling(nrow(or_w_gps_df)/5)*4)))



#try my luck with movebank #-----------------------------------------------------------------
env_request <- or_w_gps_df %>% 
  select(row_id, individual_local_identifier, location_long_closest_gps_raw, location_lat_closest_gps_raw, timestamp_closest_gps_raw, lat_group) %>% 
  mutate(timestamp_closest_gps_raw = paste0(as.character(timestamp_closest_gps_raw), ".000")) %>% 
  rename('location-long' = location_long_closest_gps_raw,
         'location-lat' = location_lat_closest_gps_raw,
         timestamp = timestamp_closest_gps_raw)

#save each group separately as a csv
env_request %>% 
  group_by(lat_group) %>%
  group_walk(~ write_csv(.x, paste0("laterality_annotation_group_", .y$lat_group, ".csv")))

#groups 4 and 5 failed for wind annotations. Break them up further based on year (they are the only zones that each have data for 3 years)
leftovers <- or_w_gps_df %>% 
  filter(lat_group %in% c("four", "five")) %>% 
  mutate(year = year(timestamp_closest_gps_raw),
         new_groups = case_when(
           lat_group == "four" & year == 2022 ~ "four_1",
           lat_group == "four" & year == 2023 ~ "four_2",
           lat_group == "four" & year == 2024 ~ "four_3",
           lat_group == "five" & year == 2022 ~ "five_1",
           lat_group == "five" & year == 2023 ~ "five_2",
           lat_group == "five" & year == 2024 ~ "five_3",
         ))  %>% 
  select(row_id, individual_local_identifier, location_long_closest_gps_raw, location_lat_closest_gps_raw, timestamp_closest_gps_raw, new_groups, lat_group, year) %>% 
  mutate(timestamp_closest_gps_raw = paste0(as.character(timestamp_closest_gps_raw), ".000")) %>% 
  rename('location-long' = location_long_closest_gps_raw,
         'location-lat' = location_lat_closest_gps_raw,
         timestamp = timestamp_closest_gps_raw)

leftovers %>% 
  group_by(new_groups) %>%
  group_walk(~ write_csv(.x, paste0("laterality_annotation_group_", .y$new_groups, ".csv")))


#groups 4-1, 4-2, 5-1, 5-2 failed for wind annotations. Break them up further based on year

leftovers2 <-  or_w_gps_df %>% 
  mutate(year = year(timestamp_closest_gps_raw),
         month = month(timestamp_closest_gps_raw)) %>% 
  filter(lat_group %in% c("four", "five") & year %in% c(2022, 2023)) %>% 
  mutate(new_groups = case_when(
    lat_group == "four" & year == 2022 & month %in% c(1:6) ~ "four_1_1",
    lat_group == "four" & year == 2022 & month %in% c(7:12) ~ "four_1_2",
    lat_group == "four" & year == 2023 & month %in% c(1:6) ~ "four_2_1",
    lat_group == "four" & year == 2023 & month %in% c(7:12) ~ "four_2_2",
    lat_group == "five" & year == 2022 & month %in% c(1:6) ~ "five_1_1",
    lat_group == "five" & year == 2022 & month %in% c(7:12) ~ "five_1_2",
    lat_group == "five" & year == 2023 & month %in% c(1:6) ~ "five_2_1",
    lat_group == "five" & year == 2023 & month %in% c(7:12) ~ "five_2_2"
  ))  %>% 
  select(row_id, individual_local_identifier, location_long_closest_gps_raw, location_lat_closest_gps_raw, timestamp_closest_gps_raw, new_groups, lat_group, year) %>% 
  mutate(timestamp_closest_gps_raw = paste0(as.character(timestamp_closest_gps_raw), ".000")) %>% 
  rename('location-long' = location_long_closest_gps_raw,
         'location-lat' = location_lat_closest_gps_raw,
         timestamp = timestamp_closest_gps_raw)


leftovers2 %>% 
  group_by(new_groups) %>%
  group_walk(~ write_csv(.x, paste0("laterality_annotation_group_", .y$new_groups, ".csv")))



# 4_1_2, 4_2_1, and 5_1_2, 5_2_1 failed
leftovers3_rows <-  leftovers2 %>% 
  filter(new_groups %in% c("four_2_1", "four_1_2", "five_1_2", "five_2_1"))

leftovers3 <- or_w_gps_df %>% 
  filter(row_id %in% leftovers3_rows$row_id) %>% 
  mutate(year = year(timestamp_closest_gps_raw),
         month = month(timestamp_closest_gps_raw),
         hour = hour(timestamp_closest_gps_raw),
         new_lat_groups = case_when(
           location_lat_closest_gps_raw >= 5.5 ~ "new_4",
           between(location_lat_closest_gps_raw, 1.5, 5.5) ~ "new_5",
           between(location_lat_closest_gps_raw, -2, 1.5) ~ "new_6",
           between(location_lat_closest_gps_raw, -5.5, -2) ~ "new_7",
           between(location_lat_closest_gps_raw, -9, -5.5) ~ "new_8",
           between(location_lat_closest_gps_raw, -12.5, -9) ~ "new_9",
           between(location_lat_closest_gps_raw, -16, -12.5) ~ "new_10",
           between(location_lat_closest_gps_raw, -19.5, -16) ~ "new_11",
           between(location_lat_closest_gps_raw, -23, -19.5) ~ "new_12",
           location_lat_closest_gps_raw <= -23 ~ "new_13")) %>% 
  #long_group = ifelse(location_long_closest_gps_raw >= 15, 2, 1)) %>% #also use longitude to group the rows :/
  mutate(new_groups = paste0(year, "_", new_lat_groups)) %>% 
  select(row_id, individual_local_identifier, location_long_closest_gps_raw, location_lat_closest_gps_raw, timestamp_closest_gps_raw, new_groups, new_lat_groups, year, month, hour) %>% 
  mutate(timestamp_closest_gps_raw = paste0(as.character(timestamp_closest_gps_raw), ".000")) %>% 
  rename('location-long' = location_long_closest_gps_raw,
         'location-lat' = location_lat_closest_gps_raw,
         timestamp = timestamp_closest_gps_raw)

leftovers3 %>% 
  group_by(new_groups) %>%
  group_walk(~ write_csv(.x, paste0("laterality_annotation_group_", .y$new_groups, ".csv")))



# append annotated data to original dataset  

#read in all annotation files for wind............... there are many NAN values!!!!!!!!

wind_data <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/laterality_annotations/wind", pattern = ".csv", full.names = T) %>% 
  map(read.csv)


#download directly from cds #-----------------------------------------------------------------

#specify time periods for which data should be downloaded: not all years have data for all months, so specify this and submit the requests separately
hours <- unique(hour(or_w_gps_df$timestamp_closest_gps_raw)) %>% str_pad(2,"left","0")
yr_mn <- data.frame(yr = c(rep(2022, 5), rep(2023, 12), rep(2024, 4)),
                    mnth = c(8:12, 1:12, 1:4) %>% str_pad(2,"left","0"))  
days <- days <- c(1:31) %>% str_pad(2,"left","0")

# 
# 
# ## cant set up the API. downloading from the gui
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2022"],
#   "month": [
#     "08", "09", "10",
#     "11", "12"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "03:00", "04:00", "05:00",
#     "06:00", "07:00", "08:00",
#     "09:00", "10:00", "11:00",
#     "12:00", "13:00", "14:00",
#     "15:00", "16:00", "17:00",
#     "18:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
# 
# 
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2023"],
#   "month": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "04:00", "05:00", "06:00",
#     "07:00", "08:00", "09:00",
#     "10:00", "11:00", "12:00",
#     "13:00", "14:00", "15:00",
#     "16:00", "17:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
# 
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2024"],
#   "month": [
#     "01", "02", "03",
#     "04"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "06:00", "07:00", "08:00",
#     "09:00", "10:00", "11:00",
#     "12:00", "13:00", "14:00",
#     "15:00", "16:00", "17:00",
#     "18:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
