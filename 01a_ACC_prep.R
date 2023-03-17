#This code downloads honey buzzard ACC data, transforms the units, and calculates ACC related features.
#Elham Nourani PhD.
#Feb 15. 2023. Konstanz, DE.

library(move)
library(moveACC)
library(tidyverse)
library(lubridate)
library(mapview)
library(parallel)
library(sf)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/")

# focus on one individual: D324-512
creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

#check for sensors available. Acc = 2365683; Mag = 77740402; Orientation: 819073350
getMovebankSensors(EHB_FN_id,login = creds)

# STEP 1: download data and convert units to g -------------------------------------------------

acc_g <- getMovebankNonLocationData(EHB_FN_id, sensorID = 2365683, animalName = c("D324_512", "D320_475"), login = creds, removeDuplicatedTimestamps = T) %>% 
  TransformRawACC(units = "g")

# STEP 2: estimate metrics -------------------------------------------------

lapply(split(acc_g, acc_g$individual_local_identifier), function(x){
  
  animalID <- unique(x$individual_local_identifier)
  
  #wing beat frequency and odba
  wave_acc <- ACCwave(x, transformedData = T, showProgress = F)
  
  #flight vs non-flight assignment based on amplitude and odba.
  wave_behav_acc <- WingBeatsSelection(wave_acc,forclustering = c("amplitude","odbaAvg"), minbeat = 0, maxbeat = max(wave_acc$beatsSec)) %>% 
    rename(flight_status = behavior, # the flight/non-flight column name is behavior. rename to flight_status
           individual_local_identifier = individualID,
           event_id = event.id,
           tag_id = tagID)
  
  #append to original acc data
  acc_metrics <- wave_behav_acc %>% 
    full_join(x, by = c("timestamp", "individual_local_identifier", "event_id"))
  
  saveRDS(acc_metrics, file = paste0("ACCsegmentation_Mar23/animal_",animalID,"_classifiedAcc.rds"))
})

# STEP 2: stitch the gps and acc together -------------------------------------------------

source("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/git_repository/HB_ontogeny/ACCtoGPS_ftns.R")
acc_files <- list.files("ACCsegmentation_Mar23", full.names = T)
gps_files <- list.files("GPSsegmentation_Mar23/classifiedData", full.names = T)

#list of individual IDs
IDs <- sapply(acc_files, function(x) substr(x, 30,37))


lapply(IDs, function(ind){
  
  acc_df <- readRDS(grep(ind, acc_files, fixed = T, value = T))
  load(grep(ind, gps_files, fixed = T, value = T)) #HRdf_smooth
  
 #save these objects as rds in the previous steps
  gps_df <- HRdf_smooth; rm(HRdf_smooth)

  # Create new ACC event column (using the tag id and the burst id assigned during the acc classification)
  acc_df$new_accEventId <- as.numeric(paste0(acc_df$tag_id.x, acc_df$burstID))
  
  if(nrow(acc_df)==length(unique(acc_df$new_accEventId))){
    #format timestamps
    # gps_df$timestamp <- as.POSIXct(as.character(gps_df$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    # acc_df$timestamp <- as.POSIXct(as.character(acc_df$timestamp), format="%Y-%m-%d %H:%M:%OS", tz="UTC")
    # Order both datasets by timestamp
    acc_df <- acc_df[order(acc_df$timestamp),]
    gps_df <- gps_df[order(gps_df$timestamp),]
    # Extract column names to associate to the GPS dataset
    accColsToAssociate <- c("beatsSec","amplitude","odbaAvg","odbaMedian","accAxes","numberSamplesPerAxis","burstDurationSecs",
                            "samplingFreqPerAxis", "eobs_accelerations_raw", "accelerationTransformed", "flight_status")
    
    # Define time tolerance (e.g. 5 mins in seconds) to look for the closest ACC information and associate it to the gps data
    timeTolerance <- 5*60 
    accGps <- ACCtoGPS(GPSdata = gps_df, ACCdata = acc_df, accEventCol="new_accEventId",
                       timeTolerance = timeTolerance,
                       ColsToAssociate = accColsToAssociate) #columns of the acc data that you want to associate to the gps data
    # Save the resulting dataset (merging gps and acc) for each individual
    saveRDS(accGps, file = paste0("ACC_GPS_behav_class/",ind,"_gps&acc_behavClass.rds"))
  }else{warning("ACC event id are not unique!")}
  #return(accGps)
})

# Bind all individuals and save final dataset
accGps_df <- lapply(list.files("ACC_GPS_behav_class", full.names = T), readRDS) %>% 
  reduce(rbind)

saveRDS(accGps_df, file = "ACC_GPS_behav_class/two_inds_gps&acc_behavClass.rds")
