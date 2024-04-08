#script for calculating Euler angles from quaterions provided by eobs.
#Apr. 08.2024. Elham Nourani, PhD.


library(tidyverse)


#open sample quat data
sample_quat <- read.csv("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_orientation/D163_696_quat_mag_w_gps.csv") %>% 
  dplyr::select(orientation_quaternions_raw)

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



########## playing around with various packages #################################
library(orientlib) #devtools::install_github("dmurdoch/orientlib")

#open sample quat data
sample_quat <- read.csv("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_GPS_IMU/matched_gps_orientation/D163_696_quat_mag_w_gps.csv") %>% 
  dplyr::select(orientation_quaternions_raw)

one_quat <- c(32496, 112, 64, -27)


# Sample quaternion data (replace with your actual data)
quat1 <- orientlib::quaternion(one_quat)  # Example quaternion (unit quaternion for no rotation)

# Convert quaternion to rotation matrix
rotation_matrix <- rotmatrix(quat1)


rotation_matrix <- quat2rot(one_quat)


euler_anlges <- rot2eul(rotation_matrix)

# Extract bank angle from rotation matrix
bank_angle <- atan2(rotation_matrix[2, 1], rotation_matrix[1, 1]) * (180 / pi)  # Convert radians to degrees

# Print the bank angle
print(bank_angle)

