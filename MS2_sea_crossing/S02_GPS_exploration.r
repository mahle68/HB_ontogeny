#Explore flight repertoir over the open sea
#follows from S01_GPS_prep.r where I segment the data
#Elham Nourani PhD.
#Aug 24. 2023. Konstanz, DE


library(tidyverse)
library(sf)
library(mgcv)

# STEP 0: Open segmented data -------------------------------------------------------------------------------------------------------------------------------------------------------
sea_sf <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/GPS_seg_Aug23/classified_data", 
                       full.names = T) %>% 
  map(readRDS) %>% 
  bind_rows() %>% 
  mutate(row_id = row_number())

# STEP 1: Env-Data annotation prep ---------------------------------------------------------------------------------------------------------------------------------------------------
sea_df <- sea_sf %>% 
  mutate('location-lat' = st_coordinates(.)[,2],
         'location-long' = st_coordinates(.)[,1],
         timestamp = paste0(as.character(timestamp), ".000")) %>% 
  select('location-long', 'location-lat', "timestamp", "row_id")

write.csv(sea_df, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_to_annotate.csv", row.names = F)


#open annotated files
p950_level <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/sea_annotations/HB_sea_to_annotate.csv-8533425805309564468/HB_sea_to_annotate.csv-8533425805309564468.csv") %>% 
  rename(geopotential = ECMWF.ERA5.PL.Geopotential,
         u_950 = ECMWF.ERA5.PL.U.Wind,
         v_950 = ECMWF.ERA5.PL.V.Wind)  %>% 
  select(- c("location.long", "location.lat", "timestamp"))

surface_level <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/sea_annotations/HB_sea_to_annotate.csv-1397362286642345513/HB_sea_to_annotate.csv-1397362286642345513.csv") %>% 
  rename(blh = ECMWF.ERA5.SL.Boundary.Layer.Height,
         sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u_10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.U.Component.,
         v_10m = ECMWF.ERA5.SL.Wind..10.m.above.Ground.V.Component.) %>% 
  mutate(delta_t = sst - t2m) %>% 
  select(- c("location.long", "location.lat", "timestamp"))

sea_ann <- sea_sf %>% 
  full_join(p950_level, by = "row_id") %>% 
  full_join(surface_level, by = "row_id") %>% 
  mutate(wind_speed_950 = sqrt(u_950^2 + v_950^2),
         wind_speed_10m = sqrt(u_10m^2 + v_10m^2))
  
saveRDS(sea_ann, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/sea_gps_seg_ann.rds")


# STEP 3: Data exploration ---------------------------------------------------------------------------------------------------------------------------------------------------

sea_sf <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds")
  
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
