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
## Step 1: decide on the unique hours for which data should be retrieved #####
#-----------------------------------------------------------------------------

#open data filtered and matched with GPS in L04_full_workflow_r

#or_w_gps_df <- readRDS("thinned_laterality_w_gps.rds") #2022-08-20 17:16:14.0000 to 2024-04-15 10:00:46.0000 

all_gps_apr <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24.rds") %>%  #2022-09-25 07:56:41.0000 to 2024-04-15 10:02:33.0000
  filter(!is.na(individual_local_identifier)) %>% 
  mutate(yr = year(timestamp),
         mn = month(timestamp),
         dy = day(timestamp),
         hr = hour(timestamp),
         unique_hr = paste(yr,mn,dy,hr, sep = "_"),
         closest_hr = round(timestamp, units="hours") %>% as.character())

# unique_hrs <- all_gps_apr  %>% 
#   group_by(closest_hr) %>% 
#   slice(1) #7327 unique hours; 7311 closest hour

all_gps_ls <- split(all_gps_apr, all_gps_apr$yr)

#convert everything to SpatVector for ease of extraction later on. this was not a good idea. PC crashed when trying to split one year of data into a list of hours
#all_gps_ls <- lapply(all_gps_ls, function(yr) vect(yr, geom = c("location_long", "location_lat"), crs =  "EPSG:4326"))

#---------------------------------------------------
## Step 2: extract wind u and v for each point #####
#---------------------------------------------------

#list nc files
nc_files <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/laterality_annotations/WIND_FROM_CDS",
                       pattern = ".nc", full.names = T)

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
  num_cores <- detectCores() - 4
  
  #for each hour
  (st_time <- Sys.time())
  wind_this_hr <- lapply(x_ls, function(y){
    
    #extract the unique hourc
    unique_hr <- y$closest_hr[1]
    
    #extract the corresponding rasters
    wind <- c(
      u[[grep(unique_hr, names(u))]],
      v[[grep(unique_hr, names(v))]]
    )
    #also could already annotate with wind speed. but having raw u and v could be helpful later for sure
    #wind_speed <- sqrt(u[[grep(unique_hr, names(u))]]^2 + v[[grep(unique_hr, names(v))]]^2)
    
    if(nrow(wind) > 0 ){ #if there is matching env data for this hour, extract the wind data for each location.
      
      #convert tracking data to SpatVector
      y <- vect(y, geom = c("location_long", "location_lat"), crs =  "EPSG:4326")
      
      #extract values for each point
      y_df <- y %>% 
        extract(x = wind, y = ., method = "bilinear", bind = T) %>% 
        data.frame(., geom(.)) %>% 
        dplyr::select(-c("geom", "part", "hole")) %>% 
        rename_with(~ str_sub(.x, start = -25, end = -21), matches("u_900|v_900")) %>%  #remove timestamp from the column name
        rename(location_lat = y,
               location_long = x)
      
    } else { #if there are no matching wind data for this hour, still create y_df, but with NA values for u and v
      
      y_df <- y %>% 
        mutate(u_900 = NA,
               v_900 = NA)
      
    }
    
    return(y_df)
    
  }) %>% 
    bind_rows()
  
  Sys.time() - st_time
  
  
})
#-------------------------------
### or
#extract manually using ncdf4

#to extract meta-data
nc_yr <- nc_open(nc_files[grep(yr, nc_files)])

#extract lon and lat
lat <- ncvar_get(nc_yr,'latitude')
nlat <- dim(lat) 
lon <- ncvar_get(nc_yr,'longitude')
nlon <- dim(lon) 

#extract the time
t <- ncvar_get(nc_yr, "valid_time")
nt <- dim(t)

#reference time according to the nc file
ref_date <- ymd_hms("1970-01-01 00:00:00")

#convert the hours into date + hour
timestamp <- ref_date + seconds(t)

#put everything in a large df
row_names <- expand.grid(lon, lat, timestamp)

#extract data for variables
var_df <- lapply(vname, function(i){
  df <- data.frame(as.vector(ncvar_get(nc, i)))
  colnames(df) <- paste0("data_", i)
  df
}) %>% 
  bind_cols() %>% 
  bind_cols(row_names)

colnames(var_df) <- c(vname, "lon", "lat", "date_time") #set column names





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
#     "08", "09", "10",
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
