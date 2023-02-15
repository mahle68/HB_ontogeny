#This code downloads honey buzzard Mag data, transforms the units, and calculates Mag related features.
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

mag <- getMovebankNonLocationData(EHB_FN_id, sensorID = 77740402, animalName = "D324_512", login = creds, removeDuplicatedTimestamps = T)


#necessary steps based on eobs manual (p. 112):
#1) find the center of the sphere
#2) subtract the center from the raw vectors to get the direction of the magnetic field
# “the Earth's field ranges between approximately 25 and 65 µT”, depending on the location on earth. Sensitivity is 0.15 MicroTesla. 
#Therefore, the radius of the sphere of the raw values is between 25/0.15=167 and 65/0.15=433

#metrics to estimate based on Hannah's paper:
#

#Richard has a script for compass distortion correction: Gundog.Compass.R