#identify birds that have arrived at wintering grounds
#October 10, 2022. Konstanz, DE
#Elham Nourani, PhD.

#assume wintering as going below 15 degrees latitude

library(tidyverse)
library(lubridate)
library(move)
library(sf)
library(rgdal)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
#load movebank login
load("/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/login.RData")

#download data for the day:
creds <- login
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)
EHB_CH_id <- getMovebankID("European Honey Buzzard Switzerland", creds)

#timerange
start <- as.character(today()-1) %>% 
  str_remove_all(pattern = "-") %>% 
  str_c("100000000")

end <- as.character(today()) %>% 
  str_remove_all(pattern = "-") %>% 
  str_c("100000000")


# finland -------------------------------------------------------
fn_today <- getMovebankData(EHB_FN_id,
                        removeDuplicatedTimestamps = T, login = creds, 
                        timestamp_start = start, timestamp_end = end)

wintering_sf <- fn_today %>% 
  as.data.frame() %>% 
  filter(location_lat < 15) %>% 
  group_by(tag_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  st_as_sf(coodrs = c("location_long.1", "location_lat.1"), crs = wgs)


wintering_IDs <- fn_today %>% 
  as.data.frame() %>% 
  filter(location_lat < 15) %>% 
  distinct(tag_local_identifier)
  
# switzerland -----------------------------------------------------
ch_today <- getMovebankData(EHB_CH_id,
                            removeDuplicatedTimestamps = T, login = creds, 
                            timestamp_start = start, timestamp_end = end)

wintering_IDs <- ch_today %>% 
  as.data.frame() %>% 
  filter(location_lat < 15) %>% 
  distinct(tag_local_identifier)
