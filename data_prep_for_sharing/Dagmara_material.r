#This code prepares one csv file for Dagmara's testing purposes. follows from 01_c_GPS_prep.r
#Elham Nourani PhD.
#Apr 2, 2024. Konstanz, DE

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")


#prep a csv for the bird that flew over Switzerland. ID == "D329_015"
gps <- readRDS("R_files/all_gps_seg_ann_Nov2023.rds")

ch_bird <- gps %>% 
  filter(individual_local_identifier == "D329_015") %>% 
  select(flight_clust_sm3, turn_angle, vert_speed, heading)

#open csv file with the flight metrics
ch_bird <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_acc/Flight metrics/D329_015_flightMetrics.csv") %>% 
  select(individual_local_identifier, t_gps, lon, lat, ground_speed_closest_gps, heading_closest_gps, 22:27) #the format is similar to eobs imu. Each row has a character string, with values separated by " "

#create a flat table. did not select summary variables for each burs, so that there would be one value per row.
ch_flat <- ch_bird %>% 
  mutate_at(vars(-1), strsplit(., " "))

#create an empty data.frame
ch_flat <- data.frame(matrix(NA, nrow = 0, ncol = ncol(ch_bird)-1))
colnames(ch_flat) <- names(ch_bird)[-1]

for(i in nrow(ch_bird))
  
acc_g <- acc %>%
  mutate(
    eobs_acceleration_g = purrr::map_chr(
      strsplit(as.character(eobs_accelerations_raw), " "),
      ~ as.character(unlist(.x) %>% as.numeric() %>% g_transform()) %>% str_c(collapse = " ")
    )
  )


################ prep a golden eagle file for now!! #######################

one_eagle <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/segmented_eagle_tracks/299029866_classifiedBursts_df.rds") %>% 
  select(5:8, 10:11, 14:17,22:24, 35) %>% 
  select(-ground.speed) %>% 
  rename(ground.speed = gr.speed,
         flight.category = flightClust_smooth3)
  
write.csv(one_eagle, "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Mentoring_Teaching/Dagmara_ETH_MSC_2024/golden_eagle_flight_classes.csv", row.names = F)

