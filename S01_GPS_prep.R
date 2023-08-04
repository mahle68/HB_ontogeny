#This code filters honey buzzard tracks for sea-crossing sections. then segment flight types.
#follows Martina's code for segmentation: https://github.com/kamransafi/GoldenEagles/blob/main/WP3_Soaring_Ontogeny/MS1_soaring_skills/script1_GPSsegmentation_goldenEagles_newMarch2023.R
#Elham Nourani PhD.
#Aug 2. 2023. Canberra, AU

library(tidyverse)
library(lubridate)
library(sf)
library(mapview)
library(terra)
library(lwgeom)
#library(move2)
#library(units)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

# STEP 1: open all GPS data and filter for latitudes #####
gps <- readRDS("/media/enourani/Ellham's HDD/Elham_EHB/all_GPS_Apr4.rds") %>%
  as.data.frame() %>%
  filter(between(location_lat,53.5,61.05) | between(location_lat, 28.9,46.49))

#open coastline layer
coastlines <- st_read("/home/enourani/ownCloud/Work/GIS_files/continent_shapefile/World_Continents.shp") %>%
  st_crop(xmin = -17, xmax = 43, ymin = -35.6, ymax = 67) %>%
  st_union()
  #poly2nb(st_make_valid(shp))

gps_sf <- gps %>%
  st_as_sf(coords = c("location_long.1", "location_lat.1"), crs = wgs) 

dd <- gps_sf  %>%
  st_difference(coastlines)

saveRDS(dd, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_tracks.rds")

# STEP 2: calculate flight altitude #####

geo <- rast("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/GE_ontogeny_of_soaring/R_files/EGM96_us_nga_egm96_15.tif")

sea <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_tracks.rds") %>% 
  #st_drop_geometry() %>% 
  select(local_identifier, location_lat, location_long, eobs_start_timestamp, timestamp, eobs_battery_voltage, eobs_horizontal_accuracy_estimate, ground_speed, heading, height_above_ellipsoid,
         tag_local_identifier, sex, taxon_canonical_name, timestamp_start, timestamp_end) %>% 
  #st_transform(crs = crs(geo)) %>% 
  extract(x = geo, y = ., method = "simple", bind = T) %>% 
  st_as_sf() %>% 
  mutate(height_msl = height_above_ellipsoid - EGM96_us_nga_egm96_15)


#look at one sample
one <- sea %>% 
  filter(local_identifier == "D324_512")

# STEP 3: prepare for annotation #####
#there are NAs in the heigh columns. create a row id column. then only use the complete columns to annotate. row id will help with merging afterwards
sea_df <- sea %>% 
  as.data.frame() %>% 
  mutate(row_id = row_number())

saveRDS(sea_df, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_tracks_df.rds")
  

ann_df <- sea_df %>% 
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) %>% 
  select(location_lat, location_long, row_id, timestamp)

colnames(ann_df)[c(2,1)] <- c("location-long","location-lat") 

write.csv(ann_df, file = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_to_annotate.csv", row.names = F)

# STEP 4: flight segmentation #####

sea_ls <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/HB_sea_tracks_df.rds") %>% 
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) %>% 
  group_split(local_identifier)
  #mt_as_move2(coords = c("location_long", "location_lat"), time_column = "timestamp", track_id_column = "local_identifier", 
  #            track_attributes = c("taxon_canonical_name", "sex", "tag_local_identifier"))
# 
# rm(sea_sf)
# gc(gc())
# 
# sea_ls <- split(sea_sf, sea_sf$local_identifier)
  
#mv_ls <- split(sea_mv, mt_track_id(sea_mv))

#dir.create("GPS_seg_Apr23/classifiedData/")

#define variables for segmentation
min_res_tl <- 2 # 1 to max 2 sec timelag
min_burst_d <- 30 # we want bursts of at least 30 secs
sw_vs <- 2 #smoothing window of 5 seconds (< min burst duration, 2 before 2 after each loc) for vertical speed for later classification
sw_th <- 12 #smoothing window of 12 seconds for thermalling behavior
circl_deg <- 250 #degrees of rotation to be achieved in the time defined by sw_th*2
min_behav_d <- 5 #minimum duration in seconds of a specific behavior, when less than this and if in between two segments of a different behavior it will be incorporated in the previous and following segment 
min_therm_d <- 30 #minimum duration for a circling event to be considered as thermalling


(b <- Sys.time())

lapply(sea_ls, function(ind){
  
  # STEP 1: calculate track variables -------------------------------------------------------------------------
  ind_id <- unique(ind$local_identifier)
  
  ind <- ind %>% 
    group_by(timestamp) %>% #remove duplicated timestamps
    slice(1) %>% 
    ungroup() %>% 
    arrange(timestamp) %>% 
    mutate(time_lag_sec = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "secs") %>%  as.numeric()),
           burst_id = cumsum(time_lag_sec > min_res_tl),
           azimuth_rad = c(NA, st_geod_azimuth(.)),
           step_length = st_distance(.),
           azimuth =  (azimuth_rad * 180)/pi,
           azimuth_positive = if_else(azimuth >= 0, azimuth, azimuth + 360)) %>% #assign one burst id to rows of consecutive 1Hz gps
    group_by(burst_id) %>% 
    mutate(altitude_diff = if_else(row_number() == 1, NA, lag(height_above_ellipsoid,1) - height_above_ellipsoid),
           vert_speed = if_else(row_number() == 1, NA, altitude_diff/time_lag_sec), #m/s
           turn_angle = if_else(row_number() == 1 | row_number() == n(), NA,
                                ((lag(azimuth_positive) - azimuth_positive + 180) %% 360) - 180),
           ground_speed = if_else(row_number() == 1, NA, step_length/time_lag_sec))
  

   
    #   time_lag_sec = mt_time_lags(.) %>% set_units("sec") %>% drop_units(),
    #        altitude_diff = lag(height_above_ellipsoid,1) - height_above_ellipsoid ,
    #        vert_speed = altitude_diff/time_lag_sec, #m/s
    #        azimuth = mt_azimuth(.),
    #        turn_angle = ifelse(row_number() == 1 | row_number == nrow(mv), NA,)
    #        #turn_angle_base = calculate_turning_angle(),
    #        turn_angle_rad = mt_turnangle(.) - set_units(pi, "rad"),
    #        turn_angle =  turn_angle_rad %>% set_units("degrees") %>% drop_units(),
    #        step_length = mt_distance(.),
    #        ground_speed = step_length/ti(.) %>% drop_units()) %>% 
    # # with cumsum R assigns the same value to all consecutive locations for which the condition is false (timelag <= 1 sec)
    # mutate(burst_id = cumsum(time_lag_sec > min_res_tl)) #create burst IDs for bursts of continous 1 sec data
    # 
  
  # mv_b <- mv %>% 
  #   as.data.frame() %>% 
  #   group_by(timestamp) %>% #remove duplicated timestamps
  #   slice(1) %>% 
  #   ungroup() %>% 
  #   arrange(timestamp) %>% 
  #   mutate(time_lag_sec = mt_time_lags(.) %>% set_units("sec") %>% drop_units(),
  #          altitude_diff = lag(height_above_ellipsoid,1) - height_above_ellipsoid ,
  #          vert_speed = altitude_diff/time_lag_sec, #m/s
  #          azimuth = mt_azimuth(.),
  #          turn_angle = ifelse(row_number() == 1 | row_number == nrow(mv), NA,)
  #          #turn_angle_base = calculate_turning_angle(),
  #          turn_angle_rad = mt_turnangle(.) - set_units(pi, "rad"),
  #          turn_angle =  turn_angle_rad %>% set_units("degrees") %>% drop_units(),
  #          step_length = mt_distance(.),
  #          ground_speed = mt_speed(.) %>% drop_units()) %>% 
  #   # with cumsum R assigns the same value to all consecutive locations for which the condition is false (timelag <= 1 sec)
  #   mutate(burst_id = cumsum(time_lag_sec > min_res_tl)) #create burst IDs for bursts of continous 1 sec data
  # 
  #mapview(st_as_sf(mv_b %>%  as.data.frame(), crs = wgs), zcol = "azimuth")
  #lwgeom::st_geod_azimuth()
  
  #set_units(pi, "rad"
  
  bursts_to_keep <- ind %>% 
    st_drop_geometry() %>% 
    group_by(burst_id) %>% 
    summarize (freq = n()) %>% 
    filter(freq >= min_burst_d) %>% 
    pull(burst_id)
  
  ind <- ind %>% 
    filter(burst_id %in% bursts_to_keep)
  
  
  if(nrow(ind) > 0){
    
    df_hr_bursts <- ind %>% 
      mutate(location_long = st_coordinates(.)[,1],
             location_lat = st_coordinates(.)[,2]) %>% 
      as.data.frame()
    
    # Split each individual dataframe by burst ID
    burst_ls_corr <- split(df_hr_bursts, df_hr_bursts$burst_id)
    
    #set the first value of timelag to NA
    burst_ls_corr <- lapply(burst_ls_corr, function(x){
      x <- x %>% 
        mutate(time_lag_sec = replace(time_lag_sec, 1, NA))
      return(x)
      })
    
    #only keep bursts with nrows >= min_burst_d. 
    burst_ls_corr_sub <- burst_ls_corr %>% 
      keep(function(x) nrow(x) >= min_burst_d)
    
    # Keep only bursts with min_burst_d (30 of smoothing window will be NA) 
    #burst_ls_corr_sub <- burst_ls_corr[which(sapply(burst_ls_corr, nrow) >= min_burst_d)]
    
    # Compute smoothed turning angle separately for each burst
    df_hr <- lapply(burst_ls_corr_sub, function(b){ 
      
      b <- b %>%
        mutate(vert_speed_smooth = ifelse(row_number() <= sw_vs | row_number() > n() - sw_vs,
                                          NA,
                                          map_dbl((sw_vs + 1):(n() - sw_vs), ~ { #Within the ifelse() call, we use map_dbl() from the purrr package to loop through the row indices for which smoothing can be performed ((sw_vs + 1):(n() - sw_vs))
                                            #Inside the map_dbl() call, we define an anonymous function that calculates the mean of vert_speed values within the range of (i - sw_vs) to (i + sw_vs) for each valid row index i. The na.rm = TRUE argument ensures that NA values are removed before calculating the mean.
                                            i <- .x
                                            mean(vert_speed[(i - sw_vs):(i + sw_vs)], na.rm = TRUE)
                                          })),
               turn_angle_smooth = ifelse(row_number() <= sw_th | row_number() > n() - sw_th,
                                          NA,
                                          map_dbl((sw_th + 1):(n() - sw_th), ~ { 
                                            i <- .x
                                            max(abs(cumsum(turn_angle[(i - sw_th):(i + sw_th)])))
                                          })))
      return(b)
      
    }) %>% 
      bind_rows() %>% 
      as.data.frame()
  
    
    # Classify soaring only based on vertical speed
    df_hr <- df_hr %>% 
      drop_na(vert_speed_smooth) %>% 
      mutate(kmean_cluster = kmeans(vert_speed_smooth, 2)$cluster) #get two clusters
    
    #find the cluster that matches soaring (i.e. the cluster with the higher mean vertical speed)
    soar_id <- df_hr %>%
      group_by(kmean_cluster) %>% #take the mean of smoothed vertical speed for each kmeans cluster
      summarise(mean_vert_speed_smooth = mean(vert_speed_smooth)) %>%
      ungroup() %>%
      filter(mean_vert_speed_smooth == max(mean_vert_speed_smooth)) %>%
      pull(kmean_cluster) 
    
    df_hr <- df_hr %>% 
      mutate(soar_clust = if_else(kmean_cluster == soar_id, "soar", "glide"),
             flight_clust = if_else(soar_clust == "glide", "gliding",
                                   if_else(turn_angle_smooth >= circl_deg, "circular_soaring",
                                           if_else(turn_angle_smooth < circl_deg, "linear_soaring",
                                                   "other")))
             ) 
    
    
    
    sf_s <- st_as_sf(df_hr, coords = c("location_long", "location_lat"), crs = wgs)
    mapview(sf_s, zcol = "flight_clust")
    
    dd <- df_hr %>%  slice(532:542) %>% 
      st_as_sf( coords = c("location_long", "location_lat"), crs = wgs)
    cumsum_angles <- abs(cumsum(dd$turn_angle))
    
    #
    # soar_clust <- rep("glide", length(kmeanV$cluster))
    # 
    # soar_clust[which(kmeanV$cluster==soar_id)] <- "soar"
    # 
    # df_hr$soar_clust <- factor(soar_clust, levels=c("soar","glide"))  
    
    # Now classify thermalling only based on turning angle (cumulated to a 25 s time window in previous step)
    # df_hr$flight_clust <- "other"
    # df_hr$flight_clust[which(df_hr$soar_clust=="soar" & df_hr$turn_angle_smooth >= circl_deg)] <- "circular soaring" #complete 150 degrees in 15 sec
    # df_hr$flight_clust[which(df_hr$soar_clust=="soar" & df_hr$flight_clust != "circular soaring")] <- "linear soaring"
    # df_hr$flight_clust[which(df_hr$soar_clust=="glide")] <- "gliding"
    # df_hr$flight_clust <- factor(df_hr$flight_clust, levels=c("circular soaring","linear soaring","gliding","other"))
    # 
    # Add some steps of smoothing based on duration of behaviors:
    burst_ls_class <- split(df_hr, df_hr$burst_id)
    burst_ls_class_smooth <- lapply(burst_ls_class, function(b){
      # We assign a unique ID to each consecutive flight segment based on a rule (the ID increases if the class of each obs is different from the previous one)
      print(unique(b$burst_id))
      b <- b[order(b$timestamp),]
      b$flightNum <- c(0, cumsum(b$flight_clust[-1] != b$flight_clust[-nrow(b)]))
      
      # we calculate the duration and unique class of each behavioral segment
      behavDuration <- merge(aggregate(timelag.sec~flightNum, data=b, FUN=sum), 
                             aggregate(flight_clust~flightNum, data=b, FUN=unique), by="flightNum")
      # we create a new category ID, which is the same as the original one, 
      # unless the duration of the segment is < 5 sec and the behav before and after are the same, in which case it all becomes one segment
      behavDuration$flightNum_smooth <- behavDuration$flightNum
      behavDuration$flight_clust_smooth <- behavDuration$flight_clust
      if(nrow(behavDuration)>2){
        for(i in 2:(nrow(behavDuration)-1)){
          if(behavDuration$timelag.sec[i] <= min_behav_d & behavDuration$flight_clust_smooth[i-1] == behavDuration$flight_clust_smooth[i+1]){
            behavDuration$flightNum_smooth[c(i, i+1)] <- behavDuration$flightNum_smooth[i-1]
            behavDuration$flight_clust_smooth[c(i, i+1)] <- behavDuration$flight_clust_smooth[i-1]
          }}
      }
      if(nrow(behavDuration)==2){
        if(any(behavDuration$timelag.sec <= min_behav_d)){
          longestBehav <- which.max(behavDuration$timelag.sec)
          behavDuration$flightNum_smooth <- behavDuration$flightNum_smooth[longestBehav]
          behavDuration$flight_clust_smooth <- behavDuration$flight_clust_smooth[longestBehav]
        }
      }
      b <- merge(b, behavDuration[c("flightNum","flightNum_smooth","flight_clust_smooth")], by="flightNum", all.x=T)
      
      # recalculate segment duration based on smoothed classification and reclassify as linear soaring all circling that lasts < 30 seconds
      behavDuration_smooth <- merge(aggregate(timelag.sec~flightNum_smooth, data=b, FUN=sum), 
                                    aggregate(flight_clust_smooth~flightNum_smooth, data=b, FUN=unique), by="flightNum_smooth")
      behavDuration_smooth$flightNum_smooth2 <- behavDuration_smooth$flightNum_smooth
      behavDuration_smooth$flight_clust_smooth2 <- behavDuration_smooth$flight_clust_smooth
      if(nrow(behavDuration_smooth)>2){
        for(i in 2:(nrow(behavDuration_smooth)-1)){
          if(behavDuration_smooth$flight_clust_smooth2[i] == "circular soaring"){
            if(behavDuration_smooth$timelag.sec[i] <= min_therm_d & behavDuration_smooth$flight_clust_smooth2[i-1] == "linear soaring" & behavDuration_smooth$flight_clust_smooth2[i+1] == "linear soaring"){
              behavDuration_smooth$flightNum_smooth2[c(i, i+1)] <- behavDuration_smooth$flightNum_smooth2[i-1]
              behavDuration_smooth$flight_clust_smooth2[c(i, i+1)] <- behavDuration_smooth$flight_clust_smooth2[i-1]
            }}}
      }
      b <- merge(b, behavDuration_smooth[c("flightNum_smooth","flightNum_smooth2","flight_clust_smooth2")], by="flightNum_smooth", all.x=T)
      # finally check the classification of the time window at the start and end of the track
      # these first and last points can only be classified as either gliding or linear, as their classification was only based on vertical speed but not turning angle
      # so if at the start or end there is linear soaring, but they are preceded or followed by circluar soaring, they become circular
      behavDuration_smooth2 <- merge(aggregate(timelag.sec~flightNum_smooth2, data=b, FUN=sum), 
                                     aggregate(flight_clust_smooth2~flightNum_smooth2, data=b, FUN=unique), by="flightNum_smooth2")
      behavDuration_smooth2$flightNum_smooth3 <- behavDuration_smooth2$flightNum_smooth2
      behavDuration_smooth2$flight_clust_smooth3 <- behavDuration_smooth2$flight_clust_smooth2
      if(nrow(behavDuration_smooth2) >= 2){
        if(behavDuration_smooth2$timelag.sec[1] <= sw_th & 
           behavDuration_smooth2$flight_clust_smooth2[1] == "linear soaring" & behavDuration_smooth2$flight_clust_smooth2[2] == "circular soaring"){
          behavDuration_smooth2$flightNum_smooth3[1] <- behavDuration_smooth2$flightNum_smooth3[2]
          behavDuration_smooth2$flight_clust_smooth3[1] <- "circular soaring"
        }
        if(behavDuration_smooth2$timelag.sec[nrow(behavDuration_smooth2)] <= sw_th & 
           behavDuration_smooth2$flight_clust_smooth2[nrow(behavDuration_smooth2)] == "linear soaring" & behavDuration_smooth2$flight_clust_smooth2[nrow(behavDuration_smooth2)-1] == "circular soaring"){
          behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)] <- behavDuration_smooth2$flightNum_smooth3[nrow(behavDuration_smooth2)-1]
          behavDuration_smooth2$flight_clust_smooth3[nrow(behavDuration_smooth2)] <- "circular soaring"
        }}
      # merge with burst
      b <- merge(b, behavDuration_smooth2[c("flightNum_smooth2","flightNum_smooth3","flight_clust_smooth3")], by="flightNum_smooth2", all.x=T)
      
      # Assign unique ID to the behavioural segment based on the final smoothest classification
      b$track_flight_id <- paste0(unique(b$individual_id),"_",unique(b$burst_id),"_segm_",b$flightNum_smooth3) 
      
      return(b) #return each classified and smoothed burst to a list
    }) # to make it run in parallel use llply (instead of lapply) with .parallel=T
    
    # Rbind all bursts and save classified and smoothed dataframe per individual
    df_hr_smooth <- as.data.frame(rbindlist(burst_ls_class_smooth))
    save(df_hr_smooth, file = paste0("GPS_seg_Apr23/classifiedData/animal_",animalID,"_classifiedBursts_df.rdata"))
  }
})

Sys.time() - b #1.7 hrs

#STEP 3: plotting ------------------------------------------------------------------------

burst_ls_class_smooth <- split(df_hr_smooth, df_hr_smooth$burst_id)
randomBursts <- sample(1:length(burst_ls_class_smooth), 100)

lapply(burst_ls_class_smooth[randomBursts], function(b){
  #b=burst_ls_class_smooth[["7006"]]
  
  print(unique(b$burst_id))
  animalID <- unique(b$local_identifier)
  
  b <- b[order(b$timestamp),]
  # cbind(as.character(b$flight_clust),as.character(b$flight_clust_smooth),as.character(b$flight_clust_smooth2), as.character(b$flight_clust_smooth3))
  
  # calculate aspect ratio for plot
  rangeLong <- max(b$location_long)-min(b$location_long)
  rangeLat <- max(b$location_lat)-min(b$location_lat)
  
  # Plot results
  png(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burst_id),".png"))
  par(mfrow=c(1,2))
  plot(b$timestamp, b$height_above_ellipsoid, type="l", col="darkgrey", lwd=2)
  points(b$timestamp, b$height_above_ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flight_clust], pch=19)
  legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
  
  plot(b$location_long, b$location_lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2)
  points(b$location_long, b$location_lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flight_clust], pch=19)
  dev.off()
  
  png(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burst_id),"_smooth.png"))
  par(mfrow=c(1,2))
  plot(b$timestamp, b$height_above_ellipsoid, type="l", col="darkgrey", lwd=2)
  points(b$timestamp, b$height_above_ellipsoid, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flight_clust_smooth3], pch=19)
  legend("topright", c("circular soaring","linear soaring","gliding","other"), col=alpha(c("red","darkgreen","blue","grey"),0.7), bty="n", pch=19, cex=0.7)
  
  plot(b$location_long, b$location_lat, asp=rangeLong/rangeLat, type="l", col="darkgrey", lwd=2)
  points(b$location_long, b$location_lat, col=alpha(c("red","darkgreen","blue","grey"),0.7)[b$flight_clust_smooth3], pch=19)
  dev.off()
  
  # this 3D plots are only useful when interactive, exporting them doesn not make much sense
  # plot3d(b[,c("location_long","location_lat","height_above_ellipsoid")], type="l", col="darkgrey")
  # points3d(b[,c("location_long","location_lat","height_above_ellipsoid")], col=c("red","darkgreen","blue","grey")[b$flight_clust_smooth3], size=5)
  # aspect3d(x=rangeLong/rangeLat, y=1, z=1)
  # snapshot3d(paste0("newGPSsegmentation_March2023/segmentationPlots/",animalID,"_burst",unique(b$burst_id),"_smooth3D.png"), width = 600, height = 600)
})  

