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

#STEP 1: checks -------------------------------------------------------------

#make sure sampling frequencies are consistent
lapply(IMU, function(x) x %>%  dplyr::select(contains("sampling_frequency")) %>%  distinct()) # all should be 20. as in 20 Hz sampling frequency


###fix this code. orientation has more than 30 values per row for some reason. it should be 30. this is taken into account in the next step.
#look at burst duration. as length of the imu raw data/number of axes
lapply(IMU, function(x){
  x %>%  
    dplyr::select(contains("_raw") | contains("sampling_frequency")) %>%  
    mutate(n_raw_values = sapply(strsplit(.[[1]], " "), length), #each row of raw values is a character string, with spaces in between the values.
           n_axes = 3) %>%  #the axis values are also a string. e.g "XYZ" only acc has the column that specifies the axes. but I know it's 3 for all sensors...
    mutate(samples_per_axis = n_raw_values/n_axes,
           burst_duration_sec = samples_per_axis/.[[2]]) %>% 
    distinct(burst_duration_sec) #acc = 1.2; mag = 0.5; quat = 0.66 BUT should be 10
    #distinct(samples_per_axis) #acc = 24; mag = 10; quat = 13.33 BUT should be 10 
}) 
#notes:
#the quaternion raw data have 40 entries instead of 30. every 4th number is a sort of index that should be removed.

#STEP 2: Create one dataframe with separate columns for each axis of each sensor -------------------------------------------------------------

#test the code with a sample
#sample <- lapply(IMU, function(x) x %>% slice(1:100))

b <- Sys.time()
axes_separated <- lapply(sample, function(x){
  
  sensor <- substr(x$sensor_type,1,3)[1]
  
  #the quaternion raw data have 40 entries instead of 30. every 4th number is a sort of index that should be removed.
  if(sensor == "Ori"){
    
    axes <- lapply(split(x,seq(nrow(x))), function(y){ #every 4th element is some sort of an index.... dont include
      raw_data <- unlist(strsplit(y %>% dplyr::select(contains("_raw")) %>%  pull(), " "))
      
      y %>% 
        mutate(x_axis_ = raw_data[seq(2,length(raw_data),4)] %>%  str_c(collapse = " "), #create one long character string from the 10 values
               y_axis_ = raw_data[seq(3,length(raw_data),4)]%>%  str_c(collapse = " "),
               z_axis_ = raw_data[seq(4,length(raw_data),4)]%>%  str_c(collapse = " ")) %>% 
        rename_if(stringr::str_detect(names(.), "_axis_"), ~paste(sensor, ., sep = "_"))
      
    }) %>% 
      reduce(rbind)
    
  } else {
    
    axes <- lapply(split(x,seq(nrow(x))), function(y){ #every 4th element is some sort of an index.... dont include
      raw_data <- unlist(strsplit(y %>% dplyr::select(contains("_raw")) %>%  pull(), " "))
      
      y %>% 
        mutate(x_axis_ = raw_data[seq(1,length(raw_data),3)] %>%  str_c(collapse = " "), #create one long character string from the 10 values
               y_axis_ = raw_data[seq(2,length(raw_data),3)]%>%  str_c(collapse = " "),
               z_axis_ = raw_data[seq(3,length(raw_data),3)]%>%  str_c(collapse = " ")) %>% 
        rename_if(stringr::str_detect(names(.), "_axis_"), ~paste(sensor, ., sep = "_"))
      
    }) %>% 
      reduce(rbind)
  }
  
  axes
})

Sys.time() - b


