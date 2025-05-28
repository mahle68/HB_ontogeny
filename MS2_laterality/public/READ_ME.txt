
The repository includes R scripts for reproducing the results of Safi et al 2025 "Ontogenetic Loss of Lateralisation Improves Migratory Flight Performance"

The following scripts are included
00_imu_functions.r includes the functions for coverting quaternion data to pitch, yaw, and roll angles. Plus additional functions specifically for processing data from biologging devices manufactured by e-obs GmbH
01_imu_processing.r includes the workflow for processing accelerometry and quaternion data
02_wind_annotation.r includes the workflow for downloading gps data, downloading the corresponding wind data, and annotating the gps data with wind speed
03a_data_prep_bursts.r
03b_data_prep_days.r
04_data_analysis.r
05_plots.r
000_bonus_gps_processing.r 

The following datasets are being made available through the Edmond repository xyz
"updated_life_cycle_nov24.rds" includes the meta-data for each individual, including the timing of start of different life stages
"thinned_laterality_w_gps_wind_all_filters2_public_prep.rds" includes the dataset ready for modeling
"EGM96_us_nga_egm96_15.tif"