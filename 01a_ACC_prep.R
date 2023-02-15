#This code downloads honey buzzard ACC data, transforms the units, and calculates ACC related features.
#Elham Nourani PhD.
#Feb 15. 2023. Konstanz, DE.

library(move)
library(tidyverse)
library(lubridate)
library(mapview)
library(parallel)
library(sf)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

# focus on one individual: D324-512
creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

#check for sensors available. Acc = 2365683; Mag = 77740402; Orientation: 819073350
getMovebankSensors(EHB_FN_id,login = creds)

# STEP 1: download data and convert units to g -------------------------------------------------

acc <- getMovebankNonLocationData(EHB_FN_id, sensorID = 2365683, animalName = "D324_512", login = creds, removeDuplicatedTimestamps = T) %>% 
  TransformRawACC(units = "g")

# STEP 2: estimate metrics -------------------------------------------------

#wing beat frequency and odba
wave_acc <- ACCwave(acc_g, transformedData = T, showProgress = F)

#flight vs non-flight assignment based on amplitude and odba.
wave_behav_acc <- WingBeatsSelection(wave_acc,forclustering = c("amplitude","odbaAvg"), minbeat = 0, maxbeat = max(wave_acc$beatsSec)) %>% 
  rename(flight_status = behavior) # the flight/non-flight column name is behavior. rename to flight_status

#append to original acc data
acc_metrics <- acc %>% 
  full_join(wave_behav_acc, by = c("timestamp", 
                             "individual_local_identifier" = "individualID",
                             "event_id" = "event.id",
                             "tag_id" = "tagID"))

