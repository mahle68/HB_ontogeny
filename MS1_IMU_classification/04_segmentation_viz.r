#This code visualizes the segmented trajectories of honey buzzards tagged in 2023 to validate the algorithm segments flight types for the European honey buzzards. based on my code: A01_GPS_prep.r
#Elham Nourani PhD.
#Feb 13, 2024. Konstanz, DE

library(tidyverse)
library(sf)
library(mapview)
library(plotly)
library(viridis)

#library(lubridate)

#library(terra)
#library(lwgeom)
#library(rgl)
#library(parallel)
#library(viridis)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

# STEP 1: open segmented data (prepared in 01_c_GPS_prep.r) ####
birds_23 <- c("D329_012", "D329_013", "D329_014", "D329_015", "D326_193", "D326_192")
gps_seg <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/all_gps_seg_ann_Nov2023.rds") %>% 
  filter(individual_local_identifier %in% birds_23) %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier))

# STEP 2: for each individual, visualize the bursts ####

birds_ls <- split(gps_seg, gps_seg$individual_local_identifier)

#generate random seeds for reproducibility of selected bursts for visualization
set.seed(20)
rnd_seeds <- sample(1:100, length(birds_ls))

#create a color palette that matches mapview
pal <- c("circular_soaring" = viridis(4)[[1]],
         "gliding" = viridis(4)[[4]],
         #"gliding" = "red",
         "linear_soaring" = viridis(4)[[2]],
         "shallow_circular_soaring" = viridis(4)[[3]],
         "other" = "grey") 

# pal <- c("circular_soaring" = "darkviolet",
#          "gliding" = "gold1",
#          "linear_soaring" = "limegreen",
#          "shallow_circular_soaring" = "firebrick1",
#          "other" = "grey") 




plots_3d_ls <- lapply(c(1:length(birds_ls)), function(x){
  
  
  #create an sf object for 2D plotting within the next lapply call. reduce the res for faster plotting
  track_sf <- birds_ls[[x]] %>% 
    group_by(yday(timestamp), hour(timestamp)) %>% 
    slice(1) %>% 
    st_as_sf(coords = c("location_long", "location_lat"), crs = "epsg:4326")
  
  
  bursts_4min_df <- birds_ls[[x]] %>% 
    group_by(burst_id) %>% 
    filter(n() == 236.0) %>% 
    ungroup()
  
  #randomly select n bursts
  set.seed(rnd_seeds[x])
  bursts_rnd <- sample(unique(bursts_4min_df$burst_id), 20, replace = F)
  
  select_bursts <- bursts_4min_df %>% 
    filter(burst_id %in% bursts_rnd)
  
  burst_ls <- split(select_bursts, select_bursts$burst_id)
  
  
  plots_3d <- lapply(burst_ls, function(burst){
    
    
    #burst$flight_clust_sm3 <- as.factor(burst$flight_clust_sm3, 
    #          levels = c("circular_soaring", "linear_soaring", "shallow_circular_soaring", "gliding"))
    
    burst$flight_clust_sm <- factor(burst$flight_clust_sm3, 
                                        levels = c("circular_soaring", "linear_soaring", "shallow_circular_soaring", "gliding"))
    
    
    
    p <- plot_ly(burst, x = ~location_long, y = ~location_lat, z = ~height_above_ellipsoid, 
                 colors = pal, color = ~flight_clust_sm,
                 type = "scatter3d", mode = "markers",
                 marker = list(size = 8),
                 showlegend = T,
                 text = ~paste0("vert_speed: ", round(vert_speed, 2),"\n",
                              "turn_angle: ",round(turn_angle_cum, 2),"\n",
                              "ground_speed: ",round(ground_speed,2),"\n",
                              "timestamp: ", timestamp)) %>% 
      layout(title = paste0("Bird ",unique(burst$individual_local_identifier), " - burst ", unique(burst$burst_id)),
             scene = list(xaxis = list(title = "Longitude"), 
                        yaxis = list(title = "Latitude"), 
                        zaxis = list(title = "Height above ellipsoid (m)")))
   
     return(p)
    
    #2d plot
    burst_sf <- burst %>%
      st_as_sf(coords = c("location_long", "location_lat"), crs = "epsg:4326")
    
    mapview(track_sf, alpha = 0, size = 0.2) + 
      mapview(burst_sf, zcol = "flight_clust_sm3", alpha = 0, color = pal)
    
  })
  

  
  plots_3d
  
})








