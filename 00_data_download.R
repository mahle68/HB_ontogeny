#This code downloads honey buzzard GPS and IMU data and matches them together.
#Elham Nourani PhD.
#Feb 7. 2023. Konstanz, DE.

library(move)
library(tidyverse)
library(lubridate)
library(mapview)

#download gps and IMU from movebank. focus on one individual: D324-512

creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

#check for sensors available. Acc = 2365683; Mag = 77740402; Orientation: 819073350
getMovebankSensors(EHB_FN_id,login = creds)

GPS <- getMovebankData(EHB_FN_id,  animalName = "D324_512", removeDuplicatedTimestamps = T, login = creds)

IMU <- list(
acc = getMovebankNonLocationData(EHB_FN_id, sensorID = 2365683, animalName = "D324_512", login = creds, removeDuplicatedTimestamps = T),
mag = getMovebankNonLocationData(EHB_FN_id, sensorID = 77740402, animalName = "D324_512", login = creds, removeDuplicatedTimestamps = T),
quat = getMovebankNonLocationData(EHB_FN_id, sensorID = 819073350, animalName = "D324_512", login = creds, removeDuplicatedTimestamps = T)
)

#checks -------------------------------------------------------------

#make sure sampling frequencies are consistent
lapply(IMU, function(x) x %>%  dplyr::select(contains("sampling_frequency")) %>%  distinct()) # all should be 20. as in 20 Hz sampling frequency

#look at burst duration. as length of the imu raw data/number of axes
lapply(IMU, function(x){
  x %>%  
    dplyr::select(contains("_raw") | contains("axes") | contains("sampling_frequency")) %>%  
    mutate(n_raw_values = sapply(strsplit(.[[1]], " "), length), #each row of raw values is a character string, with spaces in between the values.
           n_axes =  nchar(.[[2]])) %>%  #the axis values are also a string. e.g "XYZ"
    mutate(samples_per_axis = n_raw_values/n_axes,
           burst_duration_sec = samples_per_axis/.[[3]]) %>% 
    distinct(burst_duration_sec) #acc burs duration is 1.2 seconds. mag and quat are 0.5 seconds
    #distinct(samples_per_axis) #acc = 24; mag = 15; quat = 20 
}) 

#Create one dataframe with separate columns for each axis of each sensor -------------------------------------------------------------

lapply(IMU, function(x){
  
  sensor <- substr(x$sensor_type,1,3)[1]
  
    x %>% 
    mutate(x = )
  
  
})



