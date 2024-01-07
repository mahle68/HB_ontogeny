#look into stationary bursts to check the acc
#Elham Nourani, PhD. Konstanz, Germany
#Jan. 5. 2024


library(tidyverse)
library(sf)
library(mapview)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#function to calculate vectorial magnitude
vect_mag <- function(x, y, z) {
  sqrt(x^2 + y^2 + z^2)
}

## identify bursts with ground speed of 0, based on the most recent matching of acc and gps

#open matched gps and acc data (prepared in 01b_MARG_prep.r)
acc_w_gps <- readRDS("GPS_matched_ACC_Nov23_allbirds.rds")

gr_speed <- purrr::map(acc_w_gps, ~.x %>% drop_na(ground_speed_closest_gps)) %>% bind_rows() #~ creates a function, x is a placeholder for the current element

#calculate the vectorial magnitude of the g-transformed acc
#indices for elements that correspond to each axis
axis_1_i <- seq(1,length(sample), 3)
axis_2_i <- seq(2,length(sample), 3)
axis_3_i <- seq(3,length(sample), 3)

gr_speed <- gr_speed[1:2,] %>% 
  mutate(vect_mag_cc = purrr::map_chr(
    strsplit(eobs_acceleration_g, " ") %>% unlist() %>% as.numeric(),
    #~ unlist(.x) %>% as.numeric() %>% sqrt(sum())
    ~ as.character(sqrt(sum(.x[axis_1_i]^2, .x[axis_2_i]^2,.x[axis_3_i]^2 ))) %>% str_c(collapse = " ")
  ))

 gr_speed <- apply(gr_speed, 1, function(x){
   
   acc_g <- strsplit(x$eobs_acceleration_g, " ") %>% unlist() %>%  as.numeric()
   ax1 <- acc_g[axis_1_i]
   ax2 <- acc_g[axis_2_i]
   ax3 <- acc_g[axis_3_i]
   
   vectorial_mag <- mapply(vect_mag, ax1, ax2, ax3)
   
   x$vectorial_mag <- as.character(vectorial_mag) %>%  str_c(collapse = " ")
   
 })

#gr_speed <- lapply(gr_speed %>% 
#  mutate(vect_mag_acc = )
str_split(stationary$eobs_acceleration_g[1], " ") %>%  unlist() %>% as.numeric()

axis_1_i <- seq(1,length(sample), 3)
axis_2_i <- seq(2,length(sample), 3)
axis_3_i <- seq(3,length(sample), 3)

stationary <- gr_speed %>% 
  filter(ground_speed_closest_gps == 0)

#calculate the vectorial magnitude of the g-transformed acc

str_subset(stationary$eobs_acceleration_g[1], " ")
axis_1_i <- seq(1,length(sample), 3)
axis_2_i <- seq(2,length(sample), 3)
axis_3_i <- seq(3,length(sample), 3)



#################### copied stuff
#open files for two sample individuals
two_inds <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_acc",
                       pattern = "flap.csv", full.names = T) %>% 
  map(read.csv) %>% 
  bind_rows() %>% 
  drop_na(location_long_closest_gps) %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs) %>% 
  mutate(flapping = as.factor(flap_indicator))

#compare to the high-res gps segmentation
gps_seg <- str_subset(list.files("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/GPS_seg_Aug23/classified_data", full.names = T),
                      pattern = paste0(unique(two_inds$individual_local_identifier), collapse = '|')) %>%  #add the OR sign in between the two names!
  map(readRDS) %>% 
  bind_rows()


mapview(gps_seg, zcol = "flight_clust_sm3", alpha = 0) + 
  mapview(two_inds %>% filter(flapping == 1), color = "red") +
  mapview(two_inds %>% filter(flapping == 0), color = "black")
