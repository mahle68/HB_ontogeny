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
