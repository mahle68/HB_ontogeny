#This code explores the laterality of honey buzzards during circling flight. This script focuses on exploring these using the quaternion-derived features
#Elham Nourani PhD.
#Apr 24. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)


setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#open quat data-in original resolution: prepared in 01b_imu_processing.r

#open quat data- aggregated for every 8 sec-burst: prepared in 01b_imu_processing.r
burst_agg <- readRDS("quat_summaries_8secs_apr24.rds")

#-------------------------------------------------------------------------------------
# STEP1: what range of angles represents circling flight?
#-------------------------------------------------------------------------------------
#use overlapping GPS and Quat data to calculate this. requires matching gps and quat data (do it in 01b_imu_processing.r)



#-------------------------------------------------------------------------------------
# STEP2: during circling, what is bank angle like?
#-------------------------------------------------------------------------------------
#use roll in original resolution (non-aggregated) for all individuals within circling bouts OR take the mode of roll



#-------------------------------------------------------------------------------------
# STEP3: Calculate degree of handedness
#-------------------------------------------------------------------------------------