#prepare gps tracks to send to Malta
#Oct 10, 2022
#Elham Nourani, PhD


library(tidyverse)
library(move)
library(lubridate)
library(sf)
library(rgdal)
library(mapview)
library(plotKML)

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")

creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)
data <- getMovebankData(EHB_FN_id, animalName = c("D324_513", "D225_232"),
                        removeDuplicatedTimestamps = T, login = creds)#,
                        #timestamp_start = "20221001000000000", 
                        #timestamp_end = "20221003000000000")

dead_gps <- data %>% 
  as.data.frame() %>% 
  filter(local_identifier == "D324_513") %>% 
  arrange(desc(timestamp)) %>% 
  slice(1:30) %>% 
  #st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) %>% 
  summarize(lat = median(location_lat),
            long = median(location_long))
  

#create hrly res tracks to send to BirdLife Malta


ind_9558 <- data[data$tag_local_identifier == "9558" & data$timestamp > as.Date("2022-08-24") & data$timestamp < as.Date("2022-10-09"), ] %>% 
  as.data.frame() %>% 
  mutate(dt = date(timestamp),
         hr = hour(timestamp)) %>% 
  group_by(local_identifier, dt, hr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  rename(ring_number = "local_identifier") %>% 
  dplyr::select(c("location_lat", "location_long", "timestamp", "tag_local_identifier", "ring_number"))

write.csv(ind_9558, file = "/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/kml_files/ind_9558.csv")

coordinates(ind_9558) <- ~ location_long + location_lat
proj4string(ind_9558) <- wgs
mapview(ind_9558)

writeOGR(ind_9558, dsn="/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/kml_files/ind_9558.kml", 
         layer= "ind_9558", driver="KML")



##########
#this ind has one extra point over malta that wont be included if I only do an hrly subset. so extract it here and add it later
malta_pt <- data[data$tag_local_identifier == "9543" & data$location_lat > 35.82783 & data$location_lat < 35.86234,] %>% 
  as.data.frame() %>% 
  slice(1)

ind_9543 <- data[data$tag_local_identifier == "9543" & data$timestamp > as.Date("2022-08-24"),] %>% 
  as.data.frame() %>% 
  mutate(dt = date(timestamp),
         hr = hour(timestamp)) %>% 
  group_by(local_identifier, dt, hr) %>% 
  slice(1) %>% 
  ungroup() %>% 
  full_join(malta_pt) %>% 
  arrange(timestamp) %>%
  rename(ring_number = "local_identifier") %>% 
  dplyr::select(c("location_lat", "location_long", "timestamp", "tag_local_identifier", "ring_number"))


write.csv(ind_9543, file = "/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/kml_files/ind_9543.csv")


coordinates(ind_9543) <- ~ location_long + location_lat
proj4string(ind_9543) <- wgs

writeOGR(ind_9543, dsn="/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/kml_files/ind_9543.kml", 
         layer= "ind_9543", driver="KML")


##########
#send original res data to Nicholas

creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)
data <- getMovebankData(EHB_FN_id, animalName = "D324_513",
                        removeDuplicatedTimestamps = T, login = creds,
                        timestamp_start = "20221007000000000", 
                        timestamp_end = "20221008170000000")


malta <- data[data$location_lat < 36.10460] %>% 
  as.data.frame() %>% 
  dplyr::select(c("location_lat", "location_long", "timestamp", "local_identifier")) %>% 
  rename(ring_unmber = "local_identifier")


coordinates(malta) <- ~ location_long + location_lat
proj4string(malta) <- wgs

writeOGR(malta, dsn="/home/enourani/ownCloud/Work/Projects/HB_ontogeny_eobs/R_files/kml_files/D324_513.kml", 
         layer= "malta", driver="KML")

