#script for calculating Euler angles from quaterions provided by eobs.
#Apr. 08.2024. Elham Nourani, PhD.


#library(tidyverse)

# STEP 1: write some functions #####

#define a function to process quaternions in eobs format (eobs specific) 
process_quaternions <- function(quaternion_string, ftn) {
  #define a function to split up a vector of multiple quaternions into a list with each element as a vector of 4 floats (eobs specific) 
  vector_to_quat_ls <- function(x) {
    n <- length(x)
    split(x, rep(1:(n/4), each = 4))
  }
  
  quaternions <- str_split(quaternion_string, " ")[[1]] %>%
    as.numeric()
  
  result <- vector_to_quat_ls(quaternions) %>%
    map(ftn) %>%
    unlist() %>%
    as.character()
  
  return(str_c(result, collapse = " "))
}

#function to calculate summary statistics for a numeric vector 
 string_to_numeric <- function(x) {str_split(x, " ")[[1]] %>% 
   as.numeric()
 }

strings_to_numeric <- function(angle_strings, ftn) { #input is multiple rows of data
  
  #convert the multiple strings to one (for summarizing)
  numeric_vec <- unlist(map(angle_strings, string_to_numeric))
  
  return(numeric_vec)
  #numeric_vec <- str_split(angle_string, " ")[[1]] %>%
  #  as.numeric()
  
  #result <- numeric_vec %>%
  #  map(ftn) %>%
  #  unlist() %>%
  #  as.character()
  
  #return(str_c(result, collapse = " "))
}

#function to calculate the statistical mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


angle_summaries <- function(yaw, pitch, roll) {
  tibble(yaw = yaw, pitch = pitch, roll = roll) %>%
    mutate(
      yaw_diff = yaw - lag(yaw, default = first(yaw)),
      yaw_diff = case_when(
        yaw_diff > 180 ~ yaw_diff - 360,
        yaw_diff < -180 ~ yaw_diff + 360,
        TRUE ~ yaw_diff
      ),
      pitch_diff = pitch - lag(pitch, default = first(pitch)),
      roll_diff = roll - lag(roll, default = first(roll)),
      roll_diff = case_when(
        roll_diff > 180 ~ roll_diff - 360,
        roll_diff < -180 ~ roll_diff + 360,
        TRUE ~ roll_diff
      ),
      cumulative_yaw = cumsum(yaw_diff),
      cumulative_pitch = cumsum(pitch_diff),
      cumulative_roll = cumsum(roll_diff)
    ) %>%
    summarise(
      yaw_sd = sd(yaw),
      yaw_mean = mean(yaw),
      yaw_max = max(yaw),
      yaw_min = min(yaw),
      pitch_sd = sd(pitch),
      pitch_mean = mean(pitch),
      pitch_max = max(pitch),
      pitch_min = min(pitch),
      roll_sd = sd(roll),
      roll_mean = mean(roll),
      roll_max = max(roll),
      roll_min = min(roll),
      cumulative_yaw = last(cumulative_yaw),
      cumulative_pitch = last(cumulative_pitch),
      cumulative_roll = last(cumulative_roll)
    )
}


#the following functions are from Kami's IMU_conversion.r

#define a function for converting raw quaternion values mathematically from integers to floats (based on eobs manual) 
.convertEobs <- function(x){
  if(x[1]==-32768){
    x[2:4] <- 0 # corresponds to scalar part of Quaternion 1.0 or -1.0
  }
  # STEP 1: Convert raw quaternion values from integers to floats based on eobs manual
  # Calculate r
  r <- sqrt(x[2]^2 + x[3]^2 + x[4]^2)
  # Calculate the scalar value (w)
  qw <- x[1] / 32768  # The denominator 32768 is used for normalization from a 16 bit signed integer memory size
  # Calculate s
  if (r != 0) {
    s <- sqrt(1 - qw^2) / r
  } else {
    s <- 0
  }
  return(as.numeric(c(qw, s * x[2], s * x[3], s * x[4])))
}

#define a function to Calculate pitch angle from quaternion
get.pitch <- function(x, type=c("eobs", "quaternion")) {
  if(length(x) != 4){
    stop("Improper quaternion passed to function")
  }
  if(any(is.na(x))){
    pitchAngle <- NA
  } else {
    if(type == "eobs"){
      quat <- .convertEobs(x)
    } else {
      quat <- x
    }
    pitchAngle <- asin(2 * (quat[1]*quat[2] + quat[3]*quat[4]))
  }
  return(pitchAngle)
}

#define a function to calculate roll angle from quaternion
get.roll <- function(x, type=c("eobs", "quaternion")) {
  if(length(x) != 4){
    stop("Improper quaternion passed to function")
  }
  if(any(is.na(x))){
    rollAngle <- NA
  } else {
    if(type == "eobs"){
      quat <- .convertEobs(x)
    } else {
      quat <- x
    }
    rollAngle <- -atan2(2 * (quat[1] * quat[3] - quat[2] * quat[4]), 1.0 - 2.0 * (quat[2]^2 + quat[3]^2))
  }
  return(rollAngle)
}


#define a function to Calculate yaw angle from quaternion
get.yaw <- function(x, type=c("eobs", "quaternion")) {
  if(length(x) != 4){
    stop("Improper quaternion passed to function")
  }
  if(any(is.na(x))){
    yawAngle <- NA
  } else {
    if(type == "eobs"){
      quat <- .convertEobs(x)
    } else {
      quat <- x
    }
    yawAngle <- -1*atan2(2.0*(quat[2]*quat[3] - quat[1]*quat[4]) , 1-2*(quat[2]^2 + quat[4]^2) )
  }
  return(yawAngle)
}

# # STEP 2: apply to eobs data #####
# #open sample data for one individual
# sample_data <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_orientation/D163_696_quat_mag_w_gps.csv")
# 
# #calculate pitch, roll, and yaw
# sample_data_quat <- sample_data %>%
#   mutate(
#     pitch = process_quaternions(orientation_quaternions_raw, ~ get.pitch(.x, type = "eobs")),
#     yaw = process_quaternions(orientation_quaternions_raw, ~ get.yaw(.x, type = "eobs")),
#     roll = process_quaternions(orientation_quaternions_raw, ~ get.roll(.x, type = "eobs"))
#   )
