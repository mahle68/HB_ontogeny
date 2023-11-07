#This code segments flight types for the European honey buzzards. based on my code: A01_GPS_prep.r
#Elham Nourani PhD.
#Nov 6, 2023. Konstanz, DE

library(tidyverse)
library(lubridate)
library(sf)
library(mapview)
library(terra)
library(lwgeom)
library(rgl)
library(parallel)
library(viridis)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

#STEP 1: download gps data for all individuals (whole study) -----------------------------

creds <- movebank_store_credentials(username = "mahle68", rstudioapi::askForPassword())
HB_id <- movebank_get_study_id("European Honey Buzzard_Finland")

movebank_retrieve(study_id = 2201086728, entity_type= "tag_type")

gps <- movebank_retrieve(study_id = 2201086728, sensor_type_id = "gps",   #download data for all individuals 
                         entity_type = "event",  attributes = "all")

saveRDS(gps, file = "data/all_gps_nov_6_23.rds")

#STEP 2: segmentation based on 1 Hz GPS -----------------------------

#convert to sf object
inds_ls <- gps %>% 
  drop_na("location_long") %>% 
  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) %>% 
  group_split(individual_local_identifier)


#define variables for segmentation
min_res_tl <- 2 # 1 to max 2 sec time lag
min_burst_d <- 30 # we want bursts of at least 30 secs
sw_vs <- 2 #smoothing window of 5 seconds (> min burst duration, 2 before 2 after each loc) for vertical speed for later classification
sw_th <- 12 # HALF of the smoothing window of interest. 25 seconds for the entire thermaling behavior. 1 is row i, sw_th rows before and sw_th rows after
circl_deg <- 250 #degrees of rotation to be achieved in the time defined by sw_th*2+1 for full circular soaring
circl_deg_sh <- 80 #degrees of rotation to be achieved in the time defined by sw_th*2+1 for shallow thermal soaring
min_behav_d <- 7 #minimum duration in seconds of a specific behavior, when less than this and if in between two segments of a different behavior it will be incorporated in the previous and following segment 
min_therm_d <- 20 #minimum duration for a circling event to be considered as thermalling

#calculate rotation per second required for each circling mode
ta_per_sec_circl <- circl_deg/((2*sw_th)+1)
ta_per_sec_shallow <- circl_deg_sh/((2*sw_th)+1)

#calculate cumulative turning angles for all individuals, to then find good thresholds for separating shallow soaring and thermal soaring!


mycl <- makeCluster(detectCores()-10, export = list("inds_ls")) #only two CPUs because of RAM usage

#clusterExport(mycl, c("mv", "hr", "tolerance", "n_alt","wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(tidyverse)
})


(st_t <- Sys.time())

lapply(inds_ls, function(ind){
  
  #calculate track variables
  
  #everything is calculated from one point compared to its previous
  ind <- ind %>% 
    group_by(timestamp) %>% #remove duplicated timestamps
    slice(1) %>% 
    ungroup() %>% 
    arrange(timestamp) %>% 
    mutate(time_lag_sec = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "secs") %>% as.numeric()),
           burst_id = cumsum(time_lag_sec > min_res_tl), #increase burst id by one, every time time_lag is more than min_res_tl
           azimuth_rad = c(NA, st_geod_azimuth(.)),
           step_length = if_else(row_number() == 1, 0, st_distance(geometry, lag(geometry), by_element = TRUE) %>% as.numeric()), #meters
           azimuth =  (azimuth_rad * 180)/pi,
           azimuth_positive = if_else(azimuth >= 0, azimuth, azimuth + 360)) %>% #assign one burst id to rows of consecutive 1Hz gps
    group_by(burst_id) %>% 
    mutate(altitude_diff = if_else(row_number() == 1, NA, height_above_ellipsoid - lag(height_above_ellipsoid,1)),
           vert_speed = if_else(row_number() == 1, NA, altitude_diff/time_lag_sec), #m/s
           turn_angle = if_else(row_number() == 1 | row_number() == n(), NA,
                                ((lag(azimuth_positive) - azimuth_positive + 180) %% 360) - 180),
           ground_speed = if_else(row_number() == 1, NA, step_length/time_lag_sec)) %>% 
    ungroup()
  
  #garbage collection to free up RAM
  gc()
  
  #calculate nrow of each burst
  bursts_to_keep <- ind %>% 
    st_drop_geometry() %>% 
    group_by(burst_id) %>% 
    summarize (freq = n()) %>% 
    filter(freq >= min_burst_d) %>% 
    pull(burst_id)
  
  #keep bursts with more rows than min_burst_d 
  ind <- ind %>% 
    filter(burst_id %in% bursts_to_keep)
  
  gc()
  
  if(nrow(ind) > 0){
    
    #split the individual by burst_id
    burst_ls_corr <- ind %>% 
      group_split(burst_id)  
    
    # Compute smoothed turning angle separately for each burst
    # apply a moving window of sw_th*2+1 s to calculate the absolute cumulative sum of the turning angles (hereafter cumulative turning angle) and 
    # a moving window of sw_vs*2+1 s to calculate the average vertical speed
    df_hr <- lapply(burst_ls_corr, function(b){ 
      b <- b %>%
        mutate(vert_speed_smooth = ifelse(row_number() <= sw_vs | row_number() > n() - sw_vs,
                                          NA,
                                          map_dbl((sw_vs + 1):(n() - sw_vs), ~ { #Within the ifelse() call, we use map_dbl() from the purrr package to loop through the row indices for which smoothing can be performed ((sw_vs + 1):(n() - sw_vs))
                                            #Inside the map_dbl() call, we define an anonymous function that calculates the mean of vert_speed values within the range of (i - sw_vs) to (i + sw_vs) for each valid row index i. The na.rm = TRUE argument ensures that NA values are removed before calculating the mean.
                                            i <- .x
                                            mean(b %>% slice((i - sw_vs):(i + sw_vs)) %>% pull(vert_speed), na.rm = T)
                                          })),
               turn_angle_cum = ifelse(row_number() <= sw_th | row_number() > n() - sw_th,
                                       NA,
                                       map_dbl((sw_th + 1):(n() - sw_th), ~ { 
                                         i <- .x
                                         max(abs(cumsum(b %>% slice((i - sw_th):(i + sw_th)) %>% drop_na(turn_angle) %>% pull(turn_angle))))
                                       })))
      return(b)
    }) %>% 
      bind_rows()
    
    # classify soaring only based on vertical speed using k-means clustering
    df_hr <- df_hr %>% 
      drop_na(vert_speed_smooth) %>% 
      mutate(kmean_cluster = kmeans(vert_speed_smooth, 2)$cluster) #get two clusters
    
    #find the cluster id that matches soaring (i.e. the cluster with the higher mean vertical speed)
    soar_id <- df_hr %>%
      group_by(kmean_cluster) %>% #take the mean of smoothed vertical speed for each kmeans cluster
      summarise(mean_vert_speed_smooth = mean(vert_speed_smooth)) %>%
      ungroup() %>%
      filter(mean_vert_speed_smooth == max(mean_vert_speed_smooth)) %>%
      pull(kmean_cluster) 
    
    # classify flight based on cumulative rotation at each point AND the clustering based on vertical speed
    df_hr <- df_hr %>% 
      mutate(soar_clust = if_else(kmean_cluster == soar_id, "soar", "glide"),
             flight_clust = if_else(is.na(turn_angle_cum), NA,
                                    if_else(soar_clust == "glide", "gliding",
                                            if_else(turn_angle_cum >= circl_deg, "circular_soaring",
                                                    if_else(turn_angle_cum >= circl_deg_sh, "shallow_circular_soaring",
                                                            if_else(turn_angle_cum < circl_deg_sh, "linear_soaring",
                                                                    "other")))))
      ) 
    
    #visual check
    #mapview(df_hr, zcol = "flight_clust")
    
    # Add some steps of smoothing based on duration of behaviors:
    burst_ls_class <- df_hr %>% 
      group_split(burst_id)
    
    bursts_smooth <- lapply(burst_ls_class, function(b){ #sample:  b <- burst_ls_class[[5]]
      
      ##--------------------- SMOOTHING level 0: assign classification of the time window at the start and end of the segment (Martina does this as the third smoothing stage)
      # these first and last points can only be classified as either gliding or linear, as their classification was only based on vertical speed but not turning angle
      #calculate sum of turn_angles. If they are classified as circular soaring, but are straighter than circl_deg/((2*sw_th)+1) or circl_deg_sh/((2*sw_th)+1), assign linear soaring
      
      b <- b %>% 
        mutate(na_index = if_else(!is.na(flight_clust), NA, 
                                  if_else(row_number() %in% 1:(sw_th-2), "head", "tail"))) %>% #give a separate ID to the rows of NAs at the beginning and end of the segment
        group_by(na_index) %>% 
        mutate(duration = sum(time_lag_sec),
               turn_angle_sum = abs(sum(turn_angle, na.rm = T))) %>% 
        mutate(flight_clust = if_else(!is.na(flight_clust), flight_clust, #if the row already had a flight cluster assigned to it, keep the assignment
                                      if_else(soar_clust == "glide", "gliding",
                                              if_else(soar_clust == "soar" & (turn_angle_sum/duration) >= ta_per_sec_circl, "circular_soaring",
                                                      if_else(soar_clust == "soar" & between((turn_angle_sum/duration), ta_per_sec_shallow, ta_per_sec_circl), "shallow_circular_soaring", 
                                                              if_else(soar_clust == "soar" & (turn_angle_sum/duration) < ta_per_sec_shallow, "linear_soaring", "linear_soaring")))))) %>% 
        ungroup() %>% 
        dplyr::select(-c(duration, turn_angle_sum, na_index))
      
      ##--------------------- SMOOTHING level 1: based on duration
      #Point-based smoothing and assigning flight segment ID
      b <- b %>% 
        arrange(timestamp) %>% 
        #Before getting into flight segments, do a point-wise smoothing.
        #specifically, deal with the situations where for example there are two points within a soaring bout with 2 different assignments
        #if the previous point and the point after next are circular soaring, assign circular soaring
        mutate(flight_clust_sm = if_else(row_number() %in% c(1, 2, nrow(.) - 1, nrow(.)), flight_clust, #keep the first and last points' assignments the same. otherwise they become NA
                                         if_else(lag(flight_clust, 1) == "circular_soaring" & lead(flight_clust, 2) == "circular_soaring", "circular_soaring",
                                                 #similarly, if the previous point and the point after next are shallow circular soaring, assign shallow circular soaring
                                                 if_else(lag(flight_clust, 1) == "shallow_circular_soaring" & lead(flight_clust, 2) == "shallow_circular_soaring", "shallow_circular_soaring",
                                                         #if the point before previous point and the point after are circular soaring, assign circular soaring
                                                         if_else(lag(flight_clust, 2) == "circular_soaring" & lead(flight_clust, 1) == "circular_soaring", "circular_soaring",
                                                                 #similarly, if the point before previous point and the point after are shallow circular soaring, assign shallow circular soaring
                                                                 if_else(lag(flight_clust, 2) == "shallow_circular_soaring" & lead(flight_clust, 1) == "shallow_circular_soaring", "shallow_circular_soaring",
                                                                         flight_clust)))))) %>% 
        #assign a unique flight segment id to consecutive rows with the same flight category
        mutate(flight_clust_diff = if_else(c(0, diff(as.numeric(as.factor(flight_clust_sm)))) == 0, 0, 1), 
               flight_seg_id = cumsum(coalesce(flight_clust_diff,0))) #coalesce ignores the NAs
      
      #Segment-based smoothing. First calculate the duration and unique class of each behavioral segment   
      behav_duration <- b %>% 
        group_by(flight_seg_id) %>% 
        reframe(timestamp_start = head(timestamp,1),
                duration = sum(time_lag_sec),
                flight_clust_sm = unique(flight_clust_sm)) %>% 
        ungroup() %>% 
        arrange(timestamp_start) %>% 
        #First, if the assignment is linear soaring, duration is <5 sec, and the behav before and after are shallow soaring, assign shallow soaring
        mutate(flight_clust_sm1 = if_else(flight_clust_sm == "linear_soaring" & duration <= 17 & #the very shallow circles are 19 seconds long. so if only the start and end are classified as shallow soaring, the whole bout should be classified as such
                                            lag(flight_clust_sm, 1) == "shallow_circular_soaring" &
                                            lead(flight_clust_sm, 1) == "shallow_circular_soaring", "shallow_circular_soaring", flight_clust_sm)) %>% 
        #Next, more generally, unless the duration of the segment is < 5 sec and the behav before and after are the same, in which case it all becomes one segment
        #keep the assignments for the first and last two points the same. otherwise they will get NAs based on the next rule
        mutate(flight_clust_sm1 = if_else(row_number() %in% c(1,nrow(.)), flight_clust_sm, 
                                          #if the previous and next flight clust are the same and the current segment is less than min_behav_d sec, assign the prev/next flight clust
                                          if_else(lag(flight_clust_sm1, 1) == lead(flight_clust_sm1,1) & duration <= min_behav_d, lag(flight_clust_sm1), flight_clust_sm1)))
      
      #merge the new clusters with the original data
      b <- b %>% 
        left_join(behav_duration %>% select(flight_seg_id, flight_clust_sm, flight_clust_sm1), by = c("flight_seg_id","flight_clust_sm")) %>% 
        #calculate flight_seg_id for the new clustering
        mutate(flight_clust_diff_sm1 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm1)))) == 0, 0, 1), #assign a unique flight segment ID to consecutive rows with the same flight category
               flight_seg_id_sm1 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm1), NA,
                                           cumsum(coalesce(flight_clust_diff_sm1, 0))))
      #visual sanity check
      #mapview(b, zcol = "flight_clust_sm1")
      
      ##--------------------- SMOOTHING level 2: REPEAT smoothing based on segment duration
      #Segment-based smoothing. First calculate the duration and unique class of each behavioral segment   
      behav_duration2 <- b %>% 
        group_by(flight_seg_id_sm1) %>% 
        reframe(timestamp_start = head(timestamp,1),
                duration = sum(time_lag_sec),
                flight_clust_sm1 = unique(flight_clust_sm1)) %>% 
        ungroup() %>% 
        arrange(timestamp_start) %>% 
        #unless the duration of the segment is < 5 sec and the behav before and after are the same, in which case it all becomes one segment
        #keep the assignments for the first and last two points the same. otherwise they will get NAs based on the next rule
        mutate(flight_clust_sm2 = if_else(row_number() %in% c(1,nrow(.)), flight_clust_sm1, 
                                          #if the previous and next flight clust are the same and the current segment is less than min_behav_d sec, assign the prev/next flight clust
                                          if_else(lag(flight_clust_sm1, 1) == lead(flight_clust_sm1,1) & duration <= min_behav_d, lag(flight_clust_sm1), flight_clust_sm1)))
      
      #merge the new clusters with the original data
      b <- b %>% 
        left_join(behav_duration2 %>% select(flight_seg_id_sm1, flight_clust_sm1, flight_clust_sm2), by = c("flight_seg_id_sm1","flight_clust_sm1")) %>% 
        #calculate flight_seg_id for the new clustering
        mutate(flight_clust_diff_sm2 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm2)))) == 0, 0, 1), #assign a unique flight segment ID to consecutive rows with the same flight category
               flight_seg_id_sm2 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm2), NA,
                                           cumsum(coalesce(flight_clust_diff_sm2, 0))))
      
      
      ##--------------------- SMOOTHING level 3: based on turning angle of circular soaring bouts (Martina's code did this based on circling duration)
      #calculate sum of turn_angles. If they are classified as circular soaring, but are straighter than circl_deg/((2*sw_th)+1) or circl_deg_sh/((2*sw_th)+1), assign linear soaring
      #alternatively, if the assignment is linear soaring, but the rate of rotation is more than circl_deg/((2*sw_th)+1) or circl_deg_sh/((2*sw_th)+1), assign the appropriate soaring class
      
      behav_ta <- b %>% 
        group_by(flight_seg_id_sm2, flight_clust_sm2) %>% 
        reframe(duration = sum(time_lag_sec),
                turn_angle_sum = abs(sum(turn_angle, na.rm = T))) %>% 
        mutate(flight_clust_sm3 = if_else(flight_clust_sm2 == "circular_soaring" & (turn_angle_sum/duration) >= ta_per_sec_circl, "circular_soaring",
                                          if_else(flight_clust_sm2 == "circular_soaring" & (turn_angle_sum/duration) < ta_per_sec_shallow, "linear_soaring", 
                                                  if_else(flight_clust_sm2 == "linear_soaring" & (turn_angle_sum/duration) > ta_per_sec_circl, "circular_soaring",
                                                          if_else(flight_clust_sm2 == "linear_soaring" & (turn_angle_sum/duration) > ta_per_sec_shallow, "shallow_circular_soaring",
                                                                  flight_clust_sm2)))))
      
      #merge the new clusters with the original data
      b <- b %>% 
        left_join(behav_ta %>% select(flight_seg_id_sm2, flight_clust_sm2, flight_clust_sm3), by = c("flight_seg_id_sm2","flight_clust_sm2")) %>% 
        #calculate flight_seg_id for the new clustering
        mutate(flight_clust_diff_sm3 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm3)))) == 0, 0, 1), #assign a unique flight segment ID to consecutive rows with the same flight category
               flight_seg_id_sm3 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm3), NA,
                                           cumsum(coalesce(flight_clust_diff_sm3,0))))
      
      #visual sanity check
      #mapview(b, zcol = "flight_clust_sm3")
      
      # Assign unique ID to the behavioral segment based on the final smoothest classification
      b <- b %>% 
        mutate(track_flight_seg_id = paste(unique(b$individual_local_identifier), "burst", unique(b$burst_id), "seg", b$flight_seg_id_sm3, sep = "_")) 
      
      #return each classified and smoothed burst to a list
      return(b) 
    }) %>%
      #rbind all bursts and save classified and smoothed dataframe per individual
      bind_rows()
    
    saveRDS(bursts_smooth, file = paste0("R_files/GPS_seg_Nov23/animal_", 
                                         unique(ind$individual_local_identifier), "_classified_bursts.rds")) #this is an sf file
    gc()
  } 
})

Sys.time() - st_t # 8 hours for 31 individuals

stopCluster(mycl) 

