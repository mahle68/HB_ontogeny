#Explore environmental conditions during, before, and after flight over the sea
#Elham Nourani PhD.
#Feb 20. 2024. Konstanz, DE


library(tidyverse)
library(lubridate)
library(sf)
library(mapview)
library(ncdf4)
library(ecmwfr)
library(terra)
library(oce)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/")

# STEP 1: extract dates and geographic extent ----------------------------------------------------------------------------#####
sea_ann <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds")


dates <- sea_ann %>% 
  mutate(location_lat = st_coordinates(.)[,2],
         location_long = st_coordinates(.)[,1]) %>% 
  st_drop_geometry() %>% 
  mutate(mnth = month(timestamp),
         dy = day(timestamp),
         hr = hour(timestamp),
         yr = year(timestamp),
         yday = yday(timestamp),
         region = ifelse(location_lat > 48, "Baltic", "Mediterranean")) %>% 
  group_by(yr, region, mnth, dy, hr) %>% 
  slice(1) %>% 
  select(region, timestamp, yr, yday, mnth, dy, hr) %>% 
  ungroup()

dates %>% 
  group_by(region) %>% 
  summarize(min_yday = min(timestamp), 
            max_yday = max(timestamp))

#based on the exploration of the dates, use the following for each region: Baltic: Aug 1 - Sep ; Med: Sep - Oct


tmp_sp_extents <- data.frame(
  region = c("Blatic", "Mediterranean"),
  max_lat = c(61.31,48.16), #dont change the orders here
  min_lon = c(9.66, -8.87),
  min_lat = c(53.12, 29.61),
  max_lon = c(31.72, 40.34),
  min_month = c(8,9),
  max_month = c(9,10)
)


# STEP 2: request data ----------------------------------------------------------------------------#####

#wf_set_key(user = "27732",
#           key = "xxx",
#           service = "cds")

output_path <- "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/ERA5_sea/"

vars <- c("100m_u_component_of_wind", "100m_v_component_of_wind", "10m_u_component_of_wind",
          "10m_v_component_of_wind", "2m_temperature,sea_surface_temperature",
          "surface_pressure")

hours <- c(0:23) %>% str_pad(2, "left", "0") %>% paste0(":00")
days <- c(1:31) %>% str_pad(2,"left","0")
year <- "2022"

lapply(c("Blatic", "Mediterranean"), function(x){
  
  months <- tmp_sp_extents %>%
    filter(region == x) %>% 
    select("min_month", "max_month") %>% 
    as.numeric() %>% 
    str_pad(2,"left","0")
  
  area <- tmp_sp_extents %>%
    filter(region == x) %>% 
    select(2:5) %>% 
    as.numeric
  
  lapply(months, function(mn){
    
    request <- list(
      "dataset_short_name" = "reanalysis-era5-single-levels",
      "product_type"   = "reanalysis",
      "variable"       = vars,
      "year"           = year,
      "month"          = mn,
      "day"            = days,
      "time"           = hours,
      "area"           = area, #order: N,W,S,E
      "format"         = "netcdf",
      "target"         = paste0(year, "_", mn, ".nc"))
    
    wf_request(user = "27732",
               request = request,
               transfer = TRUE,
               path = output_path,
               verbose = TRUE)
    
  })
  
})

# STEP 3: extract data ----------------------------------------------------------------------------#####

file_list <- list.files("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/ERA5_sea",
                        pattern = ".nc", full.names = TRUE)
vname <- c("u100", "v100", "u10","v10", "t2m", "sst", "sp")


sea_env_ls <- lapply(file_list, function(x){
  
  nc <- nc_open(x)
  
  #extract lon and lat
  lat <- ncvar_get(nc,'latitude')
  nlat <- dim(lat) 
  lon <- ncvar_get(nc,'longitude')
  nlon <- dim(lon) 
  
  #extract the time
  t <- ncvar_get(nc, "time")
  nt <- dim(t)
  
  #convert the hours into date + hour
  timestamp <- as_datetime(c(t*60*60),origin = "1900-01-01")
  
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
  
  #remove points over land (NAs)
  #df <- var_df %>%
  #na.omit() %>% #let's keep points over land as well. otherwise, points where significant wave height is NA, i.e. land, will be deleted
  #  mutate(yday = yday(date_time),
  #         hour = hour(date_time),
  #         year = year(date_time)) %>%
  #  data.frame()
  
  #save(df,file = paste0("/home/enourani/ownCloud/Work/Projects/seabirds_and_storms/paper_prep/wind_fields/", 
  #                      unlist(strsplit(unlist(strsplit(x, ".nc")), "./"))[2], ".RData"))
  
  var_df
}) 

# STEP 4: Plot ----------------------------------------------------------------------------#####

#open entire GPS dataset (for plotting I need the sections of the tracks before and after the sea), subset with the time period of the sea_env data
gps <- readRDS("data/all_gps_nov_6_23.rds") %>% 
  filter(between(timestamp, min(sea_env_ls[[1]]$date_time), max(sea_env_ls[[2]]$date_time))) %>% 
  mutate(unique_hour = paste(str_pad(yday(timestamp), 3, "left", "0"), str_pad(hour(timestamp), 2, "left", "0"), sep = "_"))

#open base layer
world <- st_read("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
  #st_crop(xmin = -62, xmax = 11, ymin = -51, ymax = -25) %>%
  st_union()

#make one frame for each sea. one set of plots for each variable (wind and delta t)
wind_path <- "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/sea_crossing/BLS8_plots/wind_fileds/"
delta_t_path <- "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/sea_crossing/BLS8_plots/delta_t_fields/"

lapply(sea_env_ls, function(sea){
  
  sea <- sea %>% 
    filter(hour(date_time) %in% c(6:16)) %>% #only look at the day time
    mutate(wind_speed_10 = sqrt(u10^2 + v10^2), #m/s ,
           wind_speed_100 = sqrt(u100^2 + v100^2), #m/s ,
           delta_t = sst - t2m,
           unique_hour = paste(str_pad(yday(date_time), 3, "left", "0"), str_pad(hour(date_time), 2, "left", "0"), sep = "_")) 
  
  sf_use_s2(FALSE)
  
  region <- world %>% 
    st_crop(xmin = min(sea$lon) - 0.1, xmax = max(sea$lon) + 0.1, ymin = min(sea$lat) - 0.12, ymax = max(sea$lat) + 0.12)
  
  sea_name <- ifelse(max(sea$lat) > 60, "Baltic", "Mediterranean")
  
  ############ wind fields
  
  for(i in unique(sea$unique_hour)){
    
    plot <- ggplot() +
      geom_sf(data = region, fill = "grey75", col = "black", lwd = .6) +
      geom_raster(data = sea %>% filter(unique_hour == i), aes(x = lon, y = lat, fill = wind_speed_10)) +
      geom_segment(data = sea %>% filter(unique_hour == i), 
                   aes(x = lon, xend = lon+u10/10, y = lat, 
                       yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3) +
      geom_point(data = gps %>%  filter(unique_hour == i), aes(x = location_long, y = location_lat), 
                 size = 3, colour = "gold1") +
      coord_sf(xlim = range(sea$lon), ylim =  range(sea$lat)) +
      scale_fill_gradientn(colors = alpha(oce::oceColorsPalette(120), alpha = .7), limits = c(0,18), 
                           na.value = "white", name = "Speed\n (m/s)") +
      theme_bw()+
      theme(axis.text = element_text(size = 12, colour = 1),
            legend.text = element_text(size = 10, colour = 1), 
            legend.title = element_text(size = 12, colour = 1),
            legend.position = c(.92,.17),
            legend.background = element_rect(colour = 1, fill = "white"))+
      labs(x = NULL, y = NULL, title = sea %>%  filter(unique_hour == i) %>% .$date_time %>% .[1] %>% paste0(" UTC"))
    
    ggsave(plot = plot, filename = paste0(wind_path, sea_name,"_", i, ".jpeg"), 
           height = 8, width = 12, dpi = 300)
    
  }
  
  ############ delta_t fields
  
  for(i in unique(sea$unique_hour)){
    
    plot <- ggplot() +
      geom_sf(data = region, fill = "grey75", col = "black", lwd = .6) +
      geom_raster(data = sea %>% filter(unique_hour == i), aes(x = lon, y = lat, fill = wind_speed_10)) +
      geom_segment(data = sea %>% filter(unique_hour == i), 
                   aes(x = lon, xend = lon+u10/10, y = lat, 
                       yend = lat+v10/10), arrow = arrow(length = unit(0.12, "cm")), size = 0.3) +
      geom_point(data = gps %>%  filter(unique_hour == i), aes(x = location_long, y = location_lat), 
                 size = 3, colour = "gold1") +
      coord_sf(xlim = range(sea$lon), ylim =  range(sea$lat)) +
      scale_fill_gradientn(colors = alpha(oce::oceColorsPalette(120), alpha = .7), limits = c(0,18), 
                           na.value = "white", name = "Speed\n (m/s)") +
      theme_bw()+
      theme(axis.text = element_text(size = 12, colour = 1),
            legend.text = element_text(size = 10, colour = 1), 
            legend.title = element_text(size = 12, colour = 1),
            legend.position = c(.92,.17),
            legend.background = element_rect(colour = 1, fill = "white"))+
      labs(x = NULL, y = NULL, title = sea %>%  filter(unique_hour == i) %>% .$date_time %>% .[1] %>% paste0(" UTC"))
    
    ggsave(plot = plot, filename = paste0(wind_path, sea_name,"_", i, ".jpeg"), 
           height = 8, width = 12, dpi = 300)
    
  }
  
  
  
})








