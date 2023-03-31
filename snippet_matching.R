


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

# STEP 1: download IMU -------------------------------------------------

mag <- getMovebankNonLocationData(EHB_FN_id, sensorID = 77740402, animalName = c("D324_512", "D320_475"), login = creds, removeDuplicatedTimestamps = T)

acc_g <- getMovebankNonLocationData(EHB_FN_id, sensorID = 2365683, animalName = c("D324_512", "D320_475"), login = creds, removeDuplicatedTimestamps = T) %>% 
  TransformRawACC(units = "g")

quat <- getMovebankNonLocationData(EHB_FN_id, sensorID = 819073350, animalName = c("D324_512", "D320_475"), login = creds, removeDuplicatedTimestamps = T)

# STEP 2: open gps snippets -------------------------------------------------

circling <- readRDS("circling_snippets.rds")
gliding <- readRDS("gliding_snippets.rds")

# STEP 3: match IMU to gps -------------------------------------------------

lapply(c(1: length(circling)), function(snippet_id){
  
  snippet <- circling[[snippet_id]]
  
  m <- mag %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$ acc_closest_timestamp), max(snippet$ acc_closest_timestamp)))
  
  a <- acc_g %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$ acc_closest_timestamp), max(snippet$ acc_closest_timestamp)))
  
  q <- quat %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$ acc_closest_timestamp), max(snippet$ acc_closest_timestamp)))
  
  #put everything together
  df <- m %>% 
    full_join(q, by = "timestamp") %>% 
    full_join(a, by = "timestamp") %>% 
    mutate(snippet = paste0("circling_", snippet_id))
  
  
  write.csv(df, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/IMU_snippets/", "circling_", snippet_id , ".csv") ,row.names = F)
  
})



lapply(c(1: length(gliding)), function(snippet_id){
  
  snippet <- gliding[[snippet_id]]
  
  m <- mag %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$ acc_closest_timestamp), max(snippet$ acc_closest_timestamp)))
  
  a <- acc_g %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$ acc_closest_timestamp), max(snippet$ acc_closest_timestamp)))
  
  q <- quat %>% 
    filter(individual_local_identifier == unique(snippet$local_identifier)) %>% 
    filter(dplyr::between(timestamp, min(snippet$ acc_closest_timestamp), max(snippet$ acc_closest_timestamp)))
  
  #put everything together
  df <- m %>% 
    full_join(q, by = "timestamp") %>% 
    full_join(a, by = "timestamp") %>% 
    mutate(snippet = paste0("gliding_", snippet_id))
  
  
  write.csv(df, paste0("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/git_repository/R_files/IMU_snippets/", "gliding", snippet_id , ".csv") ,row.names = F)
  
})
