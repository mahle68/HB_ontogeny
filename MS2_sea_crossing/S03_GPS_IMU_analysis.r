#Explore flight repertoire over the open sea
#follows from S01_GPS_prep.r where I segment the data and S02_GPS_exploration.r where env annotation was done
#Elham Nourani PhD.
#Feb 16. 2024. Konstanz, DE


library(tidyverse)
library(sf)
library(mgcv)
library(mapview)
library(ggridges)
library(viridis)
library(hrbrthemes)

# STEP 1: Open annotated gps data ---------------------------------------------------------------------------------------------
sea_df <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds") %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) %>% 
  st_drop_geometry() %>% 
  as.data.frame()

sea_sf <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds")

sf_ls <- split(sea_sf, sea_sf$local_identifier)

mapview(sf_ls[[1]], zcol = "flight_clust_sm3")

# STEP 2: Open IMU metrics ---------------------------------------------------------------------------------------------

#for each individual, filter the IMU metrics for the time periods of sea-crossing (Baltic and Med seas)
flight_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/matched_gps_acc/Flight metrics",
                           patter = ".csv", full.names = T)

flight_ls <- lapply(sf_ls, function(x){ #for each individual
  
  #assign track ID to points that are max 15 min apart. use this to match with the IMU locations, but also allow for flexibility to include the IMU between GPS bursts
  
  x <- x %>% 
    arrange(timestamp) %>% 
    mutate(time_lag_min = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "mins") %>% as.numeric()),
           sea_track_id = base::cumsum(time_lag_min > 15)) #increase id by one, every time time_lag is more than 15 mins
  
  #open file that contains data for this individual
  flight_for_ind <- flight_files[[str_which(flight_files,unique(x$local_identifier))]] %>% 
    read.csv() %>% 
    mutate(dmyr = dmy(str_sub(t_quat, 1,12)),  #convert the day-month-year to POSIXct
           timestamp = as.POSIXct(paste0(dmyr, " ", str_sub(t_quat, 13, 21)), tz= "UTC")) 
  
  #for each track, extract the IMU data between the start and end time of the track
  lapply(split(x, x$sea_track_id), function(seg){
    
    
    flight_for_seg <- flight_for_ind %>% 
      filter(between(timestamp, min(seg$timestamp), max(seg$timestamp))) %>% 
      mutate(individual_local_identifier = unique(seg$local_identifier),
             sea_track_id = unique(seg$sea_track_id))
    
  }) %>% 
    bind_rows()
  
})

saveRDS(flight_ls, "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_IMU.rds")


# STEP 4: Exploration ---------------------------------------------------------------------------------------------

flight_df <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_IMU.rds") %>% 
  bind_rows()

ggplot(flight_df, aes(x = yday(timestamp), y = stdBank_nonFlap, color = individual_local_identifier)) +
  geom_point() + 
  facet_wrap(~individual_local_identifier)

#density plots of flight metrics for each region
flight_df <- flight_df %>% 
  mutate(lat_region = ifelse(location_lat_closest_gps > 50, "Baltic", "Mediterranean"))

ggplot(flight_df, aes(x = propFlap, y = lat_region, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Temp. [F]", option = "C") +
  labs(title = 'Temperatures in Lincoln NE in 2016') +
  theme_ipsum() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )

 #OR

ggplot(flight_df, aes(x = propFlap, color = lat_region)) +
  geom_density(lwd = 1.2, linetype = 1) + 
  geom_density(alpha = 0.7) + 
  scale_fill_manual(values = cols)


# ---------------------------------------------------------------------------------------------------------------------------------------------------


  
  ggplot(sea_sf, aes(x = local_identifier, fill = flight_clust_sm3)) +
  geom_bar()

ggplot(sea_ann, aes(x = flight_clust_sm3, y = wind_speed_950)) +
  geom_boxplot() +
  theme_linedraw()

ggplot(sea_ann, aes(x = flight_clust_sm3, y = delta_t)) +
  geom_boxplot() +
  theme_linedraw()

ggplot(sea_ann, aes(x = flight_clust_sm3, y = wind_speed_10m)) +
  geom_boxplot() +
  theme_linedraw()

#proportion of each flight type
sea_sf %>% 
  group_by(flight_clust_sm3) %>% 
  reframe(ratio = n()/nrow(sea_sf))


# STEP 4: modeling soaring OR flight type as a function of env variables  -----------------------------------------------------------------------------------------------------------------

sea_df <- sea_sf %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) %>% 
  st_drop_geometry() %>% 
  na.omit(flight_clust_sm3)

#multinomial gam
m <- gam(flight_clust_sm3 ~ scale(wind_speed_950) + scale(delta_t) +
           s(location_lat, location_long),
         family = multinom(K = 5), data = sea_df) #k = 1+unique(sea_df$flight_clust_sm3)
