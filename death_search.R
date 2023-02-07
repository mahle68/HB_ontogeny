#automatic detection of deaths
#10.10.2022, Konstanz, DE
#Elham Nourani, PhD

library(tidyverse)
library(lubridate)
library(move)
library(moveACC)


#download data for the day:

creds <- movebankLogin(username = "mahle68", rstudioapi::askForPassword())
EHB_FN_id <- getMovebankID("European Honey Buzzard_Finland", creds)

#timerange
start <- as.character(today()-1) %>% 
  str_remove_all(pattern = "-") %>% 
  str_c("100000000")

end <- as.character(today()) %>% 
  str_remove_all(pattern = "-") %>% 
  str_c("100000000")

acc_data <- getMovebankNonLocationData(EHB_FN_id, sensorID = 2365683, login = creds,
                                       removeDuplicatedTimestamps = T,
                                       timestamp_start = start, timestamp_end = end)

acc_m2 <- TransformRawACC(acc_data, units = "g")

#get odba
odba <- ACCstats(df = acc_m2) %>% 
  group_by(individualID) %>% 
  summarize(min_odba = min(odbaAvg),
            avg_odba = mean(odbaAvg),
            max_odba = max(odbaAvg))

### swiss birds ------------------------------------------------------------------------------

EHB_CH_id <- getMovebankID("European Honey Buzzard Switzerland", creds)

acc_data <- getMovebankNonLocationData(EHB_CH_id,sensorID = 2365683, login = creds,
                                       removeDuplicatedTimestamps = T,
                                       timestamp_start = start, timestamp_end = end)

acc_m2 <- TransformRawACC(acc_data, units="g")

#get odba
odba <- ACCstats(df=acc_m2) %>% 
  group_by(individualID) %>% 
  summarize(min_odba = min(odbaAvg),
            avg_odba = mean(odbaAvg),
            max_odba = max(odbaAvg))






