#find gps locations of the dead birds
#5.10.2022 Konstanz, DE.
#Elham Nourani, PhD.


library(tidyverse)
library(move)

creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_CH_id <- getMovebankID("European Honey Buzzard Switzerland", creds)
data <- getMovebankData(EHB_CH_id, animalName = c("R6589 (9711)", "R6583 (9709)"),
                        removeDuplicatedTimestamps = T, login = creds,
                        timestamp_start = "20221001000000000", 
                        timestamp_end = "20221003000000000")

GPS <- as.data.frame(data) %>% 
  group_by(tag_local_identifier) %>% 
  slice(1) %>% 
  select(c("location_lat","location_long")) %>% 
  as.data.frame()


