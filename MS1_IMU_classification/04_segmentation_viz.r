#This code visualizes the segmented trajectories of honey buzzards tagged in 2023 to validate the algorithm segments flight types for the European honey buzzards. based on my code: A01_GPS_prep.r
#Elham Nourani PhD.
#Feb 13, 2024. Konstanz, DE

library(tidyverse)
library(sf)
library(mapview)

library(lubridate)

library(terra)
library(lwgeom)
library(rgl)
library(parallel)
library(viridis)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

# STEP 1: open segmented data (prepared in 01_c_GPS_prep.r) ####
birds_23 <- c("D329_012", "D329_013", "D329_014", "D329_015", "D326_193", "D326_192")
gps_seg <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/all_gps_seg_ann_Nov2023.rds") %>% 
  filter(individual_local_identifier %in% birds_23) %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier))

# STEP 2: for each individual, visualize the bursts ####

birds_ls <- split(gps_seg, gps_seg$individual_local_identifier)

lapply(birds_ls, function(x){
  
  x_4mins <- x %>% 
    group_by(track_flight_seg_id) %>% 
    group_size()
    filter(n() == 5)
    
    x<- x %>% 
    st_as_sf(coords = c("location_long", "location_lat"), crs = "epsg:4326")
  #extract 4 min-long bursts
  
  
  mapview(x[10000:10500,], zcol = "flight_clust_sm3")
  
  soar <- x %>%  filter(burst_id == "545")
  mapview(soar, zcol = "flight_clust_sm3")
  
  p <- plot_ly(soar, x=~location.long, y=~location.lat, z=~height.above.ellipsoid, colors=pal,
               type="scatter3d", mode="lines", name=~burstID,
               line = list(color = 'darkgrey'),
               showlegend = F) %>%
    add_trace(data = b, x=~location.long, y=~location.lat, z=~height.above.ellipsoid, 
              colors=pal, color=~flightClust_smooth3,
              type="scatter3d", mode="markers",
              text=~paste0("vert.speed: ", round(vert.speed,2),"\n",
                           "turn.angle: ",round(turnAngle_smooth,2),"\n",
                           "ground.speed: ",round(gr.speed,2)),
              marker = list(size = 5),
              showlegend = T, inherit = F) %>%
    layout(title = paste0("Animal ",unique(b$individual.local.identifier)," - burst ",unique(b$burstID)),
           scene=list(xaxis = list(title = "Longitude"), 
                      yaxis = list(title = "Latitude"), 
                      zaxis = list(title = "Height above ellipsoid (m)"),
                      aspectmode = "manual", aspectratio=list(x=aspects["x"], y=aspects["y"], z=aspects["z"]))) #preserve long/lat aspect ratio
  print(p)
  readline(prompt="Press [enter] to continue")
  
  
  ###
  # 3D plotting interactive
  ## Investigate specific bursts in 3D (these plots are not exported)
  someBursts <- c("6","32", "2983", '3026', "3071", '6975', "7100", "7164", "3023", "3029", "3163")
  # interactive 3D plots with plot_ly just for a few bursts
  lapply(burst_ls_class_smooth[someBursts], function(b){
    animalID <- unique(b$individual.local.identifier)
    b <- b[order(b$timestamp),]
    # calculate aspect ratios for plot along 3 axes
    rangeLong <- max(b$location.long)-min(b$location.long)
    rangeLat <- max(b$location.lat)-min(b$location.lat)
    rangeLat_m <- (rangeLat*111.139)*1000 #transform range in metres (approximated) to compare with elevation
    rangeElev <- max(b$height.above.ellipsoid)-min(b$height.above.ellipsoid)
    ratioLat <- 1
    ratioLong <- rangeLong/rangeLat
    ratioElev <- rangeElev/rangeLat_m
    aspects <- c(x=ratioLong, y=ratioLat, z=ratioElev)
    aspects <- aspects/aspects[which.max(aspects)]
    # set colors
    pal <- c("circular soaring"="red",
             "linear soaring"="darkgreen",
             "gliding"="blue",
             "other"="grey")
    # interactive 3d plot
    p <- plot_ly(b, x=~location.long, y=~location.lat, z=~height.above.ellipsoid, colors=pal,
                 type="scatter3d", mode="lines", name=~burstID,
                 line = list(color = 'darkgrey'),
                 showlegend = F) %>%
      add_trace(data = b, x=~location.long, y=~location.lat, z=~height.above.ellipsoid, 
                colors=pal, color=~flightClust_smooth3,
                type="scatter3d", mode="markers",
                text=~paste0("vert.speed: ", round(vert.speed,2),"\n",
                             "turn.angle: ",round(turnAngle_smooth,2),"\n",
                             "ground.speed: ",round(gr.speed,2)),
                marker = list(size = 5),
                showlegend = T, inherit = F) %>%
      layout(title = paste0("Animal ",unique(b$individual.local.identifier)," - burst ",unique(b$burstID)),
             scene=list(xaxis = list(title = "Longitude"), 
                        yaxis = list(title = "Latitude"), 
                        zaxis = list(title = "Height above ellipsoid (m)"),
                        aspectmode = "manual", aspectratio=list(x=aspects["x"], y=aspects["y"], z=aspects["z"]))) #preserve long/lat aspect ratio
    print(p)
    readline(prompt="Press [enter] to continue")
  })
  
  ###
  
  
})





