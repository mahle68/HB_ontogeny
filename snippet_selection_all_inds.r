#This script allows for selection of relevant IMU snippets to use for exploration of IMU
#in collab with Pritish and Ellen
#March 29. 2023
#Elham Nourani, PhD. Konstanz, DE.

library(tidyverse)
library(lubridate)
library(sf)
library(mapview)


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
    select(c("burstID", "snippet_id"))
  
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
                      

lapply(snippet_ls, function(ind){
  
  snp <- ind %>% 
    mutate(timelag = as.numeric(difftime(lead(timestamp, 1), timestamp, units = "secs"))) %>% 
    mutate(burst_id = cumsum(ifelse(timelag == 1,0,1)))

  #create a directory for the snippet
  
  dir.create(paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Apr2/", 
                     unique(ind$snippet_id)))
  
  m <- IMU$mag %>% 
    filter(individual_local_identifier == unique(ind$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(ind$timestamp), max(ind$ timestamp))) %>% 
    dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "magnetic_fields_raw"))
  
  write.csv(m, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Mar31/tag8986_", 
                   unique(ind$snippet_id), "/mag.csv"))
  #add the acc conversion here:
  a <- acc_g %>% 
    filter(individual_local_identifier == unique(ind$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(ind$timestamp), max(ind$ timestamp))) %>% 
    dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "flight_status", "odbaAvg", "eobs_acceleration_g"))
  
  write.csv(a, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Mar31/tag8986_", 
                   unique(ind$snippet_id), "/acc.csv"))
  
  #add quat conversion here:
  q <- quat %>% 
    filter(individual_local_identifier == unique(ind$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(ind$timestamp), max(ind$ timestamp))) %>% 
    dplyr::select(c("tag_local_identifier", "timestamp", "eobs_start_timestamp", "sensor_type", "orientation_quaternions_raw"))
  
  write.csv(q, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/Pritish_collab_IMU/IMU_snippets_Mar31/tag8986_", 
                   unique(ind$snippet_id), "/quat.csv"))
  
})



