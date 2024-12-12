#environmental annotation of imu data being used for the laterality paper
#sub-step in L04_full_workflow.r 
#14.10.2024 Konstanz, Germany
#Elham Nourani. enourani@ab.mpg.de

library(tidyverse)
library(lubridate)
library(terra)
library(sf)
library(ncdf4)
library(parallel)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#annotating data on movebank failed. even when smaller chunks of data were annotated and returned, many rows had NaN values.
#Data download using CDS API also failed because I didn't manage to set up the new api key
#So, data was downloaded from the CDS website, separately for each year

#-----------------------------------------------------------------------------
## Step 1: prepare tracking data #####
#-----------------------------------------------------------------------------

#open data filtered and matched with GPS in L04_full_workflow_r

#or_w_gps_df <- readRDS("thinned_laterality_w_gps.rds") #2022-08-20 17:16:14.0000 to 2024-04-15 10:00:46.0000 

all_gps_apr <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24.rds") %>%  #2022-09-25 07:56:41.0000 to 2024-04-15 10:02:33.0000
  drop_na(individual_local_identifier, location_lat) %>% #remove NA individuals and NA locations.
  mutate(yr = year(timestamp),
         mn = month(timestamp),
         dy = day(timestamp),
         hr = hour(timestamp),
         unique_hr = paste(yr,mn,dy,hr, sep = "_"),
         closest_hr = round(timestamp, units = "hours") %>% as.character()) %>% 
  as.data.frame()

all_gps_ls <- split(all_gps_apr, all_gps_apr$yr)

#---------------------------------------------------
## Step 2: extract wind u and v for each point #####
#---------------------------------------------------

#list nc files
nc_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/laterality_annotations/WIND_FROM_CDS",
                       pattern = ".nc", full.names = T)

output_path <- "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr24_wind_annotated/"

(st_time <- Sys.time())
lapply(all_gps_ls, function(x){
  
  #extract the year
  yr <- x$yr[1]
  
  ########### read in wind data for this year #####
  #open nc data for this year
  u <- nc_files[grep(yr, nc_files)] %>% 
    rast("u")
  
  v <- nc_files[grep(yr, nc_files)] %>% 
    rast("v")
  
  #for some reason I can't extract time using the terra::time() function. Extract the time from the layer names
  times <- names(u) %>%
    str_split("time=") %>%
    map_chr(2)
  
  #reference time according to the nc file
  ref_date <- ymd_hms("1970-01-01 00:00:00")
  
  #convert the hours into date + hour
  timestamp <- ref_date + seconds(times)
  
  #rename the layers based on the converted timestamps. these are the same for both u and v. but keep the u and v in the names for clarity when extracting values
  #from the two layers that would otherwise have the same name (ie. the timestamp)
  names(u) <- paste0("u_900_", timestamp)
  names(v) <- paste0("v_900_", timestamp)
  
  
  ########### split the tracking data into unique hours #####
  
  x_ls <- split(x, x$closest_hr)
  
  # Define the number of cores to use
  num_cores <- detectCores() - 10 #run on 2 cores
  
  wind_this_yr <- mclapply(x_ls, function(y){ #for each hour
    
    #extract the unique hourc
    unique_hr <- y$closest_hr[1]
    
    # Check whether there is any wind data for this hour
    if (any(str_detect(names(u), unique_hr))) {
      # Extract the corresponding rasters
      wind <- c(
        u[[str_which(names(u), unique_hr)]],
        v[[str_which(names(v), unique_hr)]]
      )
      #names(wind) <- c("u_900", "v_900")
      
      # Convert tracking data to SpatVector
      y_vec <- vect(y, geom = c("location_long", "location_lat"), crs = "EPSG:4326")
      
      # Extract values for each point and append directly to y
      extracted_wind <- extract(x = wind, y = y_vec, method = "bilinear", bind = FALSE, ID = FALSE)
      colnames(extracted_wind) <- c("u_900", "v_900")
      
      # Append extracted values to y
      y_df <- y %>%
        bind_cols(as.data.frame(extracted_wind))
      
    } else {
      # If there are no matching wind data for this hour, create y_df with NA values
      y_df <- y %>%
        mutate(u_900 = as.numeric(NA),
               v_900 = as.numeric(NA))
    }
    
    #saveRDS(y_df, paste0("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr24_wind_annotated/by_hr_2022/wind_ann", unique_hr, ".rds" ))
    
    rm(wind, y)
    
    y_df
    
  }, mc.cores = num_cores) %>% 
    bind_rows()
  
  
  
  saveRDS(wind_this_yr, file = paste0(output_path,"gps_annotated_", yr, ".rds")) #had some issues with 2022 so had to do it in two batches.
  
})
Sys.time() - st_time #1.8 hours for three years


#---------------------------------------------------------------
## Step 3: put all files together and calculate wind speed #####
#---------------------------------------------------------------

ann_ls <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr24_wind_annotated", full.names = T) %>% 
  map(readRDS) %>% 
  map(bind_rows) %>% 
  bind_rows() %>% 
  select(1:48) %>% 
  mutate(wind_speed =  sqrt(u_900^2 + v_900^2))

saveRDS(ann_ls, file = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24_wind.rds")

#----------------------------------------------------- 
# ## cant set up the API. downloading from the gui
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2022"],
#   "month": [
#     "07", "08", "09", "10",
#     "11", "12"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "03:00", "04:00", "05:00",
#     "06:00", "07:00", "08:00",
#     "09:00", "10:00", "11:00",
#     "12:00", "13:00", "14:00",
#     "15:00", "16:00", "17:00",
#     "18:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
# 
# 
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2023"],
#   "month": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "04:00", "05:00", "06:00",
#     "07:00", "08:00", "09:00",
#     "10:00", "11:00", "12:00",
#     "13:00", "14:00", "15:00",
#     "16:00", "17:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
# 
# import cdsapi
# 
# dataset = "reanalysis-era5-pressure-levels"
# request = {
#   "product_type": ["reanalysis"],
#   "variable": [
#     "u_component_of_wind",
#     "v_component_of_wind"
#   ],
#   "year": ["2024"],
#   "month": [
#     "01", "02", "03",
#     "04"
#   ],
#   "day": [
#     "01", "02", "03",
#     "04", "05", "06",
#     "07", "08", "09",
#     "10", "11", "12",
#     "13", "14", "15",
#     "16", "17", "18",
#     "19", "20", "21",
#     "22", "23", "24",
#     "25", "26", "27",
#     "28", "29", "30",
#     "31"
#   ],
#   "time": [
#     "06:00", "07:00", "08:00",
#     "09:00", "10:00", "11:00",
#     "12:00", "13:00", "14:00",
#     "15:00", "16:00", "17:00",
#     "18:00"
#   ],
#   "pressure_level": ["900"],
#   "data_format": "netcdf",
#   "download_format": "zip",
#   "area": [64, -11, -27, 36]
# }
# 
# client = cdsapi.Client()
# client.retrieve(dataset, request).download()
