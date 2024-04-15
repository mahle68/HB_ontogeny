#script for calculating Euler angles from quaterions provided by eobs.
#Apr. 08.2024. Elham Nourani, PhD.


library(tidyverse)

#order of values in eobs format: w,x,y,z
raw_quat <- c(32496, 112, 64, -27) 

# STEP 1: convert raw quaternion values mathematically from integers to floats (based on eobs manual) #####

# Calculate r
r <- sqrt(raw_quat[2]^2 + raw_quat[3]^2 + raw_quat[4]^2)

#the scalar value
qw <- raw_quat[1]/32768


# Calculate s
if (r != 0) {
  s <- sqrt(1 - qw^2) / r
} else {
  s <- 0
}


q <- c(qw = qw, 
       qx = s * raw_quat[2], 
       qy = s * raw_quat[3], 
       qz = s * raw_quat[3])


# STEP 2: calculate roll, pitch and yaw (based on eobs manual, translated from C to R) #####

# Calculate pitch (x-axis rotation)
pitch <- asin(2 * (qw * qx + qy * qz))

# Calculate roll (y-axis rotation)
roll <- -atan2(2.0 * (qw * qy - qx * qz), 1.0 - 2.0 * (qx^2 + qy^2))

# Calculate yaw (not provided by eobs. taken from chatGPT's suggestion)
yaw <- atan2(2 * (qw * qz + qx * qy), qw^2 + qx^2 - qy^2 - qz^2)



#### next steps:
# prepare Quat data for the whole time period of interest
# write a function for calculating pitch roll and yaw for all the dataset: pitch, roll, and yaw for each quat value and the whole 8-sec burst (sd, mean, cumulative). 
# acc calculations for probability of flapping


#write a function to split up a vector of multiple quaternions into a list 
vector_to_quat_ls <- function(x) split(x, cut(seq_along(x), length(x)/4, labels = FALSE)) 

#write a function to transform the raw quat values from integer to float 
calculate_quat_float <- function(raw_quat) { #a numeric vector of 4 values: w,x,y,z
    
    # Calculate r
    r <- sqrt(raw_quat[2]^2 + raw_quat[3]^2 + raw_quat[4]^2)
    
    # Calculate scalar value
    qw <- raw_quat[1] / 32768
    
    # Calculate s
    if (r != 0) {
      s <- sqrt(1 - qw^2) / r
    } else {
      s <- 0
    }
    
    # Calculate quaternion components
    q <- c(qw, s * raw_quat[2], s * raw_quat[3], s * raw_quat[4])
    
    return(q)
  
}

# Calculate pitch angle from quaternion
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

# Calculate roll angle from quaternion
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

#open sample data for one individual
sample_data <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_orientation/D163_696_quat_mag_w_gps.csv") #%>% 
#dplyr::select(orientation_quaternions_raw)

#calculate pitch, roll, and yaw
sample_data_quat <- sample_data %>%
  mutate(
    quaternions_float = purrr::map(
      #split the character string including the quaternions into a numeric vector, then into a list with one element per quaternion
      strsplit(orientation_quaternions_raw, " ") %>% unlist() %>% as.numeric() %>% 
        vector_to_quat_ls(),
      ~ calculate_quat_float(.x)) 
    %>% unlist(use.names = F) %>% as.character() %>% str_c(collapse = " "),
    pitch = purrr::map(
      #split the character string including the quaternions into a numeric vector, then into a list with one element per quaternion
      strsplit(quaternions_float, " ") %>% unlist() %>% as.numeric() %>% 
        vector_to_quat_ls(),
      ~ get.pitch(.x, type = "quaternion")) %>% 
      %>% unlist(use.names = F) %>% as.character() %>% str_c(collapse = " "),
    yaw =  purrr::map(
      #split the character string including the quaternions into a numeric vector, then into a list with one element per quaternion
      strsplit(quaternions_float, " ") %>% unlist() %>% as.numeric() %>% 
        vector_to_quat_ls(),
      ~ get.yaw(.x, type = "quaternion")) %>% 
      %>% unlist(use.names = F) %>% as.character() %>% str_c(collapse = " "),
    roll =  purrr::map(
      #split the character string including the quaternions into a numeric vector, then into a list with one element per quaternion
      strsplit(quaternions_float, " ") %>% unlist() %>% as.numeric() %>% 
        vector_to_quat_ls(),
      ~ get.roll(.x, type = "quaternion")) %>% 
      %>% unlist(use.names = F) %>% as.character() %>% str_c(collapse = " ")
  )





#calculate pitch, roll, and yaw
acc_g <- acc %>%
  mutate(
    quaternions_float = purrr::map_chr(
      strsplit(as.character(orientation_quaternions_raw), " "),
      ~ as.character(unlist(.x) %>% as.numeric() %>% g_transform()) %>% str_c(collapse = " ")
    )
  )
Sys.time()-st 


