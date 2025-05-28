# script for preparing the data for Safi et al 2025.
# Elham Nourani, PhD. enourani@ab.mgp.de

#This code downloads honey buzzard IMU data, transforms the units of ACC and processes the quaternion data. 

library(move2)
library(tidyverse)
library(lubridate)

#setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#source the imu conversion functions
source("00_imu_functions.r")

#---------------------------------------------------------------------------------
## Step 0: Download all IMU data                                             #####
#---------------------------------------------------------------------------------

#Make sure you have set up your credentials following Move2 tutorials.
#Data used in this study was downloaded on 15.04.2024. Use this timestamp as the cutoff point when downloading the data to reproduce the results of this study.

#creds <- movebank_store_credentials(username = "your_user_name", rstudioapi::askForPassword())

mag <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "magnetometer", 
                         entity_type = "event",  attributes = "all",
                         timestamp_end = as.POSIXct("2024-04-15 00:00:00"))

quat <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "orientation", 
                          entity_type = "event",  attributes = "all",
                          timestamp_end = as.POSIXct("2024-04-15 00:00:00"))

acc <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "acceleration", 
                         entity_type = "event",  attributes = "all",
                         timestamp_end = as.POSIXct("2024-04-15 00:00:00"))

#put the quat and mag together
or <- mag %>% 
  full_join(quat, by = c("individual_local_identifier", "tag_local_identifier", "study_id", "tag_id", "individual_taxon_canonical_name", "timestamp",
                         "eobs_start_timestamp", "eobs_key_bin_checksum", "data_decoding_software", "individual_id", "import_marked_outlier", "visible")) %>% 
  as.data.frame()

#remove the originial files to save RAM
rm(mag, quat)

#saveRDS(or, "all_orientation_apr15_24.rds")
#saveRDS(acc, "all_acceleration_apr15_24.rds")


#---------------------------------------------------------------------------------
## Step 1: Process the ACC                                                   #####
#---------------------------------------------------------------------------------

#### ----------------------- convert units to g 

#write a ftn for g conversion
g_transform <- function(x) {
  (x-2048)*(1/1024) #slope according to the eobs manual
}

#make sure the values fall between -2 and 2
#For my tag generation, I have the same slope and intercept for all axes, so there is no need to separate into 3 axis before converting all raw values to g.
acc_g <- acc %>%
  mutate(
    eobs_acceleration_g = purrr::map_chr(
      strsplit(as.character(eobs_accelerations_raw), " "),
      ~ as.character(unlist(.x) %>% as.numeric() %>% g_transform()) %>% str_c(collapse = " ")
    )
  ) %>% 
  as.data.frame()

#check the number of acc records. we have 1.2 seconds of data per row
acc_n <- acc_g %>%
  rowwise() %>%
  mutate(acc_g_length = length(strings_to_numeric(eobs_acceleration_g)),
         acc_raw_length = length(strings_to_numeric(eobs_accelerations_raw))
  )

#### ----------------------- calculate vedba for each row of data (1.2 seconds)
acc_g <- readRDS("all_acceleration_g_apr_24.rds") %>% 
  rowwise() %>% 
  mutate(DBA_ls = list(DBA_ftn(eobs_acceleration_g)),
    VeDBA = DBA_ls[1],
    ODBA = DBA_ls[2],
    acc_g_length = length(strings_to_numeric(eobs_acceleration_g)),
    acc_raw_length = length(strings_to_numeric(eobs_accelerations_raw))) %>%
  ungroup() %>%
  select(-DBA_ls) %>% 
  as.data.frame()

#save the file. This data will be used in script 03b_data_prep_days.r to calculate daily VeDBA
saveRDS(acc_g, "your_path/your_processed_acc_data.rds")

#saveRDS(acc_g, "all_acceleration_g_dba_apr_24.rds")


#---------------------------------------------------------------------------------
## Step 2: Process the Quaternions (calculate Euler angles)                  #####
#---------------------------------------------------------------------------------

#subset the orientation dataframe for data with high sampling frequency (this should be post-fledging, migration, and wintering data for laterality tests.)
or <- or %>% 
  group_by(individual_local_identifier) %>% 
  arrange(timestamp, .by_group = T) %>% 
  mutate(timelag = c(0, diff(timestamp)),  #calculate time lag 
         imu_burst_id = cumsum(timelag > 1)) %>% #assign a unique burst ID every time timelag is larger than 1 second. The IDs will be unique only within one individual's data
  ungroup() %>% 
  as.data.frame()

#subset orientation by the length of the bursts. Only keep those that are longer than 3 second.
#first calculate burst durations
burst_size <- or %>% 
  group_by(individual_local_identifier, imu_burst_id) %>% 
  summarize(burst_duration = sum(timelag[-1])) %>% #the first value of timelag in each burst is a large value indicating timediff between this and the previous burst. remove it
  ungroup()

#remove short bursts
or_hfreq <- or %>% 
  right_join(burst_size %>% filter(burst_duration > 3)) #there are bursts with length of 0-3

#remove the original orientation dataframe and clean up RAM
rm(or);gc()

#calculate pitch, roll, and yaw for the whole dataset

or_angles <- or_hfreq %>%
  rowwise() %>% 
  mutate(
    pitch = process_quaternions(orientation_quaternions_raw, ~ get.pitch(.x, type = "eobs")),
    yaw = process_quaternions(orientation_quaternions_raw, ~ get.yaw(.x, type = "eobs")),
    roll = process_quaternions(orientation_quaternions_raw, ~ get.roll(.x, type = "eobs"))
  ) %>% 
  ungroup() %>%  
  #convert to degrees (make sure values are from -180 to +180)
  mutate(across(c("pitch", "yaw", "roll"), 
                ~ purrr::map_chr(
                  strsplit(as.character(.x), " "),
                  ~ as.character(unlist(.) %>% as.numeric() %>% deg()) %>% str_c(collapse = " ") #make sure the values are SIGNED
                ),
                .names = "{col}_deg"
  )) %>% 
  as.data.frame()


#saveRDS(or_angles, file = "quat_angles_apr24.rds")

#---------------------------------------------------------------------------------
## Step 3: Summarize the data over each second                               #####
#---------------------------------------------------------------------------------

#IMU data is collected at 20 Hz frequency and there are 1 or 1.2 seconds of data per row. summarize these so that there is one second of data per row.
#This step is not absolutely necessary. Also it is very specific to how e-obs data is stored (i.e. a long character string containing all quat or acc measurements for each burst of 1 or 1.2 seconds)
#It is useful for summarizing the quaternion-derived angles for each second. This way I can later summarize the values for each burst, depending on the burst length of interest.

#or_angles <- readRDS("quat_angles_apr24.rds")

or_seconds <- or_angles%>% 
  group_by(individual_local_identifier, floor_date(timestamp, "second")) %>% 
  #create a long character string containing all values for each second
  summarize(across(c(roll, pitch, yaw, roll_deg, pitch_deg, yaw_deg), ~paste(., collapse = " ")), 
            #retain the first value of other important columns
            across(c(study_id, individual_taxon_canonical_name, orientation_quaternions_sampling_frequency,
                     tag_local_identifier, timestamp, tag_id, imu_burst_id, burst_duration), ~head(.,1)),
            .groups = "keep") %>%
  ungroup() %>%
  as.data.frame()

#saveRDS(or_seconds, file = "quat_angles_secs_apr24.rds")

#for each second, calculate the mean, min, max, sd, cumsum for each axis
seconds_summaries <- or_seconds %>%
  rowwise() %>%
  mutate(
    angle_summary = list(angle_summaries(
      yaw = strings_to_numeric(yaw_deg),
      pitch = strings_to_numeric(pitch_deg),
      roll = strings_to_numeric(roll_deg)
    ))
  ) %>%
  ungroup() %>%
  unnest(angle_summary) %>% 
  as.data.frame()

saveRDS(seconds_summaries, file = "quat_summaries_1sec_Jul24.rds")


#sessionInfo()
# R version 4.5.0 (2025-04-11)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 22.04.5 LTS
# 
# locale:
# [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8      
# [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Europe/Berlin
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] xtable_1.8-4    patchwork_1.3.0 INLA_24.12.11   Matrix_1.7-3    terra_1.8-42    mapview_2.11.2  sf_1.0-20       move2_0.4.4     lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
# [13] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1    ggplot2_3.5.2   tidyverse_2.0.0
# 
# loaded via a namespace (and not attached):
# [1] tidyselect_1.2.1        farver_2.1.2            filelock_1.0.3          fastmap_1.2.0           leaflet_2.2.2           digest_0.6.37           timechange_0.3.0        lifecycle_1.0.4        
# [9] magrittr_2.0.3          compiler_4.5.0          rlang_1.1.6             tools_4.5.0             yaml_2.3.10             askpass_1.2.1           labeling_0.4.3          htmlwidgets_1.6.4      
# [17] bit_4.6.0               sp_2.2-0                classInt_0.4-11         curl_6.2.2              RColorBrewer_1.1-3      pkgload_1.4.0           KernSmooth_2.23-26      withr_3.0.2            
# [25] grid_4.5.0              stats4_4.5.0            e1071_1.7-16            leafem_0.2.4            scales_1.4.0            cli_3.6.5               crayon_1.5.3            ragg_1.4.0             
# [33] generics_0.1.3          rstudioapi_0.17.1       tzdb_0.5.0              DBI_1.2.3               sodium_1.4.0            proxy_0.4-27            splines_4.5.0           assertthat_0.2.1       
# [41] parallel_4.5.0          base64enc_0.1-3         vctrs_0.6.5             jsonlite_2.0.0          hms_1.1.3               bit64_4.6.0-1           corrr_0.4.4             systemfonts_1.2.3      
# [49] crosstalk_1.2.1         jquerylib_0.1.4         units_0.8-7             keyring_1.3.2           glue_1.8.0              codetools_0.2-19        leaflet.providers_2.0.0 stringi_1.8.7          
# [57] gtable_0.3.6            raster_3.6-32           pillar_1.10.2           rappdirs_0.3.3          htmltools_0.5.8.1       satellite_1.0.5         R6_2.6.1                fmesher_0.3.0          
# [65] textshaping_1.0.1       vroom_1.6.5             lattice_0.22-5          png_0.1-8               class_7.3-23            MatrixModels_0.5-4      Rcpp_1.0.14             pkgconfig_2.0.3 
