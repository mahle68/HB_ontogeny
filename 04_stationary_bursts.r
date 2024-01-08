#look into stationary bursts to check the acc
#Elham Nourani, PhD. Konstanz, Germany
#Jan. 5. 2024


library(tidyverse)
library(sf)
library(mapview)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#function to calculate vector magnitude
vect_mag <- function(x, y, z) {
  sqrt(x^2 + y^2 + z^2)
}

## identify bursts with ground speed of 0, based on the most recent matching of acc and gps

#open matched gps and acc data (prepared in 01b_MARG_prep.r)
acc_w_gps <- readRDS("GPS_matched_ACC_Nov23_allbirds.rds")

gr_speed <- purrr::map(acc_w_gps, ~.x %>% drop_na(ground_speed_closest_gps)) %>% bind_rows() #~ creates a function, x is a placeholder for the current element

#calculate the vector magnitude of the g-transformed acc
#indices for elements that correspond to each axis
axis_1_i <- seq(1,length(sample), 3)
axis_2_i <- seq(2,length(sample), 3)
axis_3_i <- seq(3,length(sample), 3)

#calculate vector magnitude
gr_speed_vec_mag <- gr_speed  %>%
  mutate(vector_mag = map_chr(eobs_acceleration_g, ~{
    acc_g <- strsplit(.x, " ") %>% unlist() %>% as.numeric()
    ax1 <- acc_g[axis_1_i]
    ax2 <- acc_g[axis_2_i]
    ax3 <- acc_g[axis_3_i]
    vector_mag <- vect_mag(ax1, ax2, ax3)
    as.character(vector_mag) %>% str_c(collapse = " ")
  }))


#extract stationary bursts based on ground speed
 stationary <- gr_speed %>% 
   filter(ground_speed_closest_gps == 0) 

#plot relationship between ground speed and accelerometry
hist(str_split(stationary$vector_mag, " ") %>% unlist() %>% as.numeric(), xlab = "Vector magnitude", ylab = "Frequency", main = "")

#plot the histogram of all vector magnitudes
#2022
#stationary22 <- stationary %>% 
#  filter(year(timestamp) == 2022)
#hist(str_split(stationary22$vector_mag, " ") %>% unlist() %>% as.numeric())

#2023
#stationary23 <- stationary %>% 
#  filter(year(timestamp) == 2023)
#hist(str_split(stationary23$vector_mag, " ") %>% unlist() %>% as.numeric())


#convert to sf object and plot on map
#stationary_sf <- stationary %>% 
# st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"), crs = wgs)
