#This script allows for selection of relevant IMU snippets to use for exploration of IMU
#in collab with Pritish and Ellen
#March 29. 2023
#Elham Nourani, PhD. Konstanz, DE.

library(tidyverse)
library(lubridate)
library(sf)
library(mapview)
library(parallel)

wgs <- "+proj=longlat +datum=WGS84 +no_defs"

setwd("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")
setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")

#write a ftn for g conversion
g_transform <- function(x) {
  (x-2048)*(1/1024) #slope according to the eobs manual
}


## STEP 1: open data  ----------------------------------------------------------------

#segmented GPS data from 02_GPS_segmentation
seg_files <- list.files("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/GPS_seg_Apr23/classifiedData", full.names = T)

#for each individual, pick 4 soaring and 4 gliding snippets
snippet_ls <- lapply(seg_files, function(ind){
  
  load(ind) #later on change the original code to save this as a rds file
  gps <- HRdf_smooth; rm(HRdf_smooth); gc(gc())
  
  #filter out non-migratory flight
  gps_migr <- gps %>% 
    filter(between(location_lat,5,60))  #& #filter out breeding and wintering season movements
             #flightClust_smooth3 != "linear soaring") %>%  #filter out linear soaring for now. some thermal soaring got in there. i need to modify martina's code for this sp
  
  #find bursts with one behavior #group consecutive bursts with the same behavior together. find a group with the largest n of burstIDs. this helps find a group with the the same behavior in consecutive bursts. easier to find many IMU recordings for the period 
  burst_groups <- gps_migr %>% 
    group_by(burstID) %>% 
    summarize(n_behavs = n_distinct(flightClust_smooth3),
              behav = head(flightClust_smooth3,1)) %>% #this only makes sense for groups with one unique behavior! 
    filter(n_behavs == 1) %>% 
    group_by(behav) %>% 
    mutate(id_diff = ifelse(row_number() == 1, 0,
                            burstID - lag(burstID,1))) %>%
    mutate(burst_group = cumsum(ifelse(id_diff == 1, 0,1))) %>% 
    ungroup()
  
  #identify the longest burst groups
  long_groups <- burst_groups %>% 
    group_by(behav, burst_group) %>% 
    summarize(n = n()) %>% 
    arrange(desc(n)) %>% 
    slice(1:4)  %>% #extract the 4 burst groups with the largest n of bursts
    mutate(snippet_id = paste0(behav, "_", burst_group)) %>% 
    ungroup()
  
  burstIDs <- burst_groups %>% 
    inner_join(long_groups, by = c("behav", "burst_group")) %>% #join with the burst_groups to get the burstID back
    dplyr::select(c("burstID", "snippet_id"))
  
  snps <- gps_migr %>% 
    inner_join(burstIDs, by = c("burstID"))
  
  #explore
  #snps_sf <- snps %>% 
  #  st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) 
    
  #mapview(snps_sf, zcol = "snippet_id")
  
  snps
  
})

saveRDS(snippet_ls, file = "IMU_snippet_ls_all_inds_per_burst.rds")

## STEP 2: download IMU  ----------------------------------------------------------------

#Dont forget to convert the acc and quat after matching with the gps snippets
creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

acc <- getMovebankNonLocationData(EHB_FN_id, sensorID = 2365683, login = creds, removeDuplicatedTimestamps = T)
mag <- getMovebankNonLocationData(EHB_FN_id, sensorID = 77740402, login = creds, removeDuplicatedTimestamps = T)
quat <- getMovebankNonLocationData(EHB_FN_id, sensorID = 819073350, login = creds, removeDuplicatedTimestamps = T)

#write as a list
IMU <- list(acc = acc,
            mag = mag,
            quat = quat)

saveRDS(IMU, file = "IMU_for_all_Apr5_23.rds") #this is all raw values. no conversions have been done


## STEP 3: Match the acc, mag and quat ----------------------------------------------------------------

IMU <- readRDS("IMU_for_all_Apr5_23.rds")
snippet_ls <- readRDS("IMU_snippet_ls_all_inds_per_burst.rds") #one list per individual. with all the snippets in it

# 
# mycl <- makeCluster(5) #the number of CPUs to use (adjust this based on your machine)
# 
# clusterExport(mycl, c("snippet_ls", "IMU", "g_transform")) #define the variable that will be used within the ParLapply call. these will be copied to all the cpus
# 
# clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
#   library(tidyverse)
#   library(lubridate)
# })

(b <- Sys.time())
#parLapply(mycl, snippet_ls, function(ind){
lapply(snippet_ls, function(ind){  
  
  for(i in unique(ind$snippet_id)){
    
    snp <- ind %>% 
      filter(snippet_id == i) %>% 
      as.data.frame() %>% 
      dplyr::select(c("tag_local_identifier", "local_identifier", "timestamp", "location_long", "location_lat",  "height_above_ellipsoid", "burstID", "vert.speed", "turn.angle", "gr.speed", "snippet_id", "flightClust_smooth3"))
    
    #visualize:
    # snp_sf <- snp %>% 
    #   st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) 
    # mapview(snp_sf, zcol = "snippet_id")
    
    #create a directory to save the files for this ind-snippet
    
    path <- paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Apr6/", unique(snp$tag_local_identifier), 
                   "_", unique(snp$snippet_id))
    
    
    #extract IMU---------------------------------------------------------------------------- 
    
    #MAG
    m <- IMU$mag %>% 
      filter(individual_local_identifier == unique(snp$local_identifier)) %>% 
      filter(dplyr::between(timestamp, min(snp$timestamp)-20, max(snp$ timestamp)+20)) %>% #adding the time ensures that we get the 8-sec IMU window before and after the gps burst  
      dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "magnetic_fields_raw"))
    
    if(nrow(m) == 0){ message(paste0("no IMU data for tag ", unique(snp$tag_local_identifier), " ",unique(snp$snippet_id), "!"))
    }else{
      
      #now that we know there is matching IMU, create a new directory and save the gps and IMU
      dir.create(path)
      write.csv(snp, paste0(path, "/gps.csv"))
      
      write.csv(m, paste0(path, "/mag.csv"))
      
      #ACC
      a <- IMU$acc %>% 
        filter(individual_local_identifier == unique(snp$local_identifier)) %>% 
        filter(dplyr::between(timestamp, min(snp$timestamp)-20, max(snp$ timestamp)+20))
      
      #if(nrow(a) == 0){ message(paste("no acc data for tag", unique(snp$tag_local_identifier), unique(snp$snippet_id), "!", sep = " "))}
      
      
      #convert the acc to units of g:
      a_g <- lapply(split(a,seq(nrow(a))), function(burst){
        
        g_data <- unlist(strsplit(burst %>% dplyr::select("eobs_accelerations_raw") %>%  pull(), " ")) %>% 
          as.numeric() %>% 
          g_transform()
        
        burst <- burst %>% 
          mutate(eobs_acceleration_g = str_c(as.character(g_data), collapse = " "))
        
      }) %>% 
        reduce(rbind) %>% 
        dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "eobs_acceleration_g"))
      
      write.csv(a_g, paste0(path, "/acc.csv"))
      
      #QUAT
      q <- IMU$quat %>% 
        filter(individual_local_identifier == unique(snp$local_identifier)) %>% 
        filter(dplyr::between(timestamp, min(snp$timestamp)-20, max(snp$ timestamp)+20)) %>% 
        dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "orientation_quaternions_raw"))
      
      #if(nrow(q) == 0){ message(paste("no quat data for tag", unique(snp$tag_local_identifier), unique(snp$snippet_id), "!", sep = " "))}
      
      #convert the quaternions:
      q_c <- lapply(split(q,seq(nrow(q))), function(burst){
        
        raw_data <- unlist(strsplit(burst %>% dplyr::select("orientation_quaternions_raw") %>%  pull(), " ")) %>% 
          as.numeric()
        
        #extract values for each axis. 
        quat_df <- data.frame(w_raw = raw_data[seq(1,length(raw_data),4)], #create one long character string from the 10 values
                              x_raw = raw_data[seq(2,length(raw_data),4)], 
                              y_raw = raw_data[seq(3,length(raw_data),4)],
                              z_raw = raw_data[seq(4,length(raw_data),4)]) %>% 
          mutate(r = sqrt(x_raw^2 + y_raw^2 + z_raw^2), #these lines are from the eobs manual p. 110. according to Pritish: skip these except for the division of w_raw by 32768. ALso, divide x, y, and z by 128
                 w = w_raw/32768) %>% 
          mutate(s = ifelse(r != 0, sqrt(1-w^2)/r, 0)) %>% 
          mutate(x = s*x_raw,
                 y = s*y_raw,
                 z = s*z_raw)
        
        
        burst <- burst %>% 
          mutate(quat_w = str_c(as.character(quat_df$w), collapse = " "),
                 quat_x = str_c(as.character(quat_df$x), collapse = " "),
                 quat_y = str_c(as.character(quat_df$y), collapse = " "),
                 quat_z = str_c(as.character(quat_df$z), collapse = " "))
        
      }) %>% 
        reduce(rbind)
      
      write.csv(q_c, paste0(path, "/quat.csv"))
      
    }
  }
})

Sys.time() - b #3.8 min

#stopCluster()
