#script for calculating Euler angles from quaterions provided by eobs.
#Apr. 08.2024. Elham Nourani, PhD.


library(tidyverse)

# STEP 1: write some functions #####

#order of quaternion values in eobs format: w,x,y,z

#define a function to split up a vector of multiple quaternions into a list with each element as a vector of 4 floats (eobs specific) 
vector_to_quat_ls <- function(x) split(x, cut(seq_along(x), length(x)/4, labels = FALSE)) 

#define a function to process quaternions in eobs format (eobs specific) 
process_quaternions <- function(quaternion_string, func) {
  strsplit(quaternion_string, " ") %>%
    unlist() %>%
    as.numeric() %>%
    vector_to_quat_ls() %>%
    map(func) %>%
    unlist(use.names = FALSE) %>%
    as.character() %>%
    str_c(collapse = " ")
}

process_quaternions <- function(quaternion_string, func) {
  #define a function to split up a vector of multiple quaternions into a list with each element as a vector of 4 floats (eobs specific)
  vector_to_quat_ls <- function(x) {
    n <- length(x)
    split(x, rep(1:(n/4), each = 4))
  }
  
  quaternions <- strsplit(quaternion_string, " ")[[1]]
  result <- unlist(lapply(vector_to_quat_ls(as.numeric(quaternions)), func))
  paste(result, collapse = " ")
}

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

# STEP 2: apply to eobs data #####
#open sample data for one individual
sample_data <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_orientation/D163_696_quat_mag_w_gps.csv")

#calculate pitch, roll, and yaw
sample_data_quat <- sample_data %>%
  mutate(
    pitch = process_quaternions(orientation_quaternions_raw, ~ get.pitch(.x, type = "eobs")),
    yaw = process_quaternions(orientation_quaternions_raw, ~ get.yaw(.x, type = "eobs")),
    roll = process_quaternions(orientation_quaternions_raw, ~ get.roll(.x, type = "eobs"))
  )
