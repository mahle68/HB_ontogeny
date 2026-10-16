#This code filters honey buzzard tracks for sea-crossing sections. then segment flight types.
#follows Martina's code for segmentation: https://github.com/kamransafi/GoldenEagles/blob/main/WP3_Soaring_Ontogeny/MS1_soaring_skills/script1_GPSsegmentation_goldenEagles_newMarch2023.R
# I have renamed turn_angle_smooth (from martina's code) to turn_angle_cum, because it is the cumulative sum and not simply a smoothed value.
#Elham Nourani PhD.
#Aug 2. 2023. Canberra, AU.
#Being updated with Claude Code assistant, Oct 16. 2026 Lausanne, CH.

library(tidyverse)
library(lubridate)
library(sf)
library(move2)
library(mapview)
library(terra)
library(lwgeom)
library(rgl)
library(viridis)
library(moveWindSpeed)

#options(rgl.printRglwidget = TRUE)

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")
setwd("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/R_files/")


#------------------------------------------------------------------------#####
#STEP 1: download all GPS data, filter for latitude, and clean up -------#####
#------------------------------------------------------------------------#####

gps <- movebank_download_study("European Honey Buzzard_Finland",
                               sensor_type_id = "gps",
                               remove_movebank_outliers = TRUE)


saveRDS(gps, file = "raw_gps_Oct26.rds")

#------------------------------------------------------------------------#####
# STEP 2: open all GPS data and filter for latitudes --------------------#####
#------------------------------------------------------------------------#####

gps <- gps %>%
  filter(between(st_coordinates(.)[,2], 53.5, 61.05) | between(st_coordinates(.)[,2], 28.9, 46.49))

cleaned_gps <- gps |> 
  mt_filter_unique(criterion = "first") |> 
    mt_clean_track(mass = 1, mode = "flying", remove = TRUE) #max mass in kg


#open coastline layer
coastlines <- st_read("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/GIS_files/continent_shapefile/World_Continents.shp") %>%
  st_crop(xmin = -17, xmax = 45, ymin = -35.6, ymax = 67) %>%
  st_union()
  #poly2nb(st_make_valid(shp))

#include a negative buffer, to retain the last few points before sea-crossing
#coastlines_buffer <- coastlines %>% 
#  st_buffer(dist = -0.1)

# gps_sf <- gps %>%
#   st_as_sf(coords = c("location_long", "location_lat"), crs = wgs) 

#extract gps points that fall over the sea and the buffer zone
#over_sea_w_buffer <- gps_sf  %>%
#  st_difference(coastlines_buffer)

#extract points that only fall over the sea
over_sea <- cleaned_gps  %>%
  mt_as_event_attribute(., names(mt_track_data(.))) %>% #after this point we lose the move2 class, so keep the attributes
  st_difference(coastlines) 

saveRDS(over_sea, file = "cleaned_gps_over_sea.rds")

#------------------------------------------------------------------------#####
# STEP 3: calculate flight altitude -------------------------------------#####
#------------------------------------------------------------------------#####

geo <- rast("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/GIS_files/EGM96_us_nga_egm96_15.tif")

over_sea <- over_sea %>%
  units::drop_units() %>%
  extract(x = geo, y = ., method = "simple", bind = T) %>%
  st_as_sf() %>%
  mutate(height_msl = height_above_ellipsoid - geoid_undulation)


#look at one sample
one <- over_sea %>% 
  filter(individual_local_identifier == "D225_236")

mapview(one, zcol = "height_msl", alpha = 0)


#------------------------------------------------------------------------#####
# STEP 4: track segmentation --------------------------------------------#####
#------------------------------------------------------------------------#####

# adapted from MS2_laterality/public/000_bonus_gps_processing.r

min_res_tl <- 2
min_burst_d <- 30
sw_vs <- 2
sw_th <- 12
circl_deg <- 50
min_behav_d <- 7
min_therm_d <- 20

ta_per_sec_circl <- circl_deg/((2*sw_th)+1)

df_hr <- over_sea %>%
  filter(individual_local_identifier == "D225_236") %>% 
  group_split(individual_local_identifier) %>%
  lapply(function(ind){

    ind <- ind %>%
      group_by(timestamp) %>%
      slice(1) %>%
      ungroup() %>%
      arrange(timestamp) %>%
      mutate(time_lag_sec = if_else(row_number() == 1, 0, difftime(timestamp, lag(timestamp), units = "secs") %>% as.numeric()),
             burst_id = cumsum(time_lag_sec > min_res_tl),
             azimuth_rad = c(NA, st_geod_azimuth(.)),
             step_length = if_else(row_number() == 1, 0, st_distance(geometry, lag(geometry), by_element = TRUE) %>% as.numeric()),
             azimuth = (azimuth_rad * 180)/pi,
             azimuth_positive = if_else(azimuth >= 0, azimuth, azimuth + 360)) %>%
      group_by(burst_id) %>%
      mutate(altitude_diff = if_else(row_number() == 1, NA, height_msl - lag(height_msl,1)),
             vert_speed = if_else(row_number() == 1, NA, altitude_diff/time_lag_sec),
             turn_angle = if_else(row_number() == 1 | row_number() == n(), NA,
                                  ((lag(azimuth_positive) - azimuth_positive + 180) %% 360) - 180),
             ground_speed = if_else(row_number() == 1, NA, step_length/time_lag_sec)) %>%
      ungroup()

    bursts_to_keep <- ind %>%
      st_drop_geometry() %>%
      group_by(burst_id) %>%
      summarize(freq = n()) %>%
      filter(freq >= min_burst_d) %>%
      pull(burst_id)

    ind <- ind %>%
      filter(burst_id %in% bursts_to_keep)

    if(nrow(ind) == 0) return(NULL)

    ind %>%
      group_split(burst_id) %>%
      lapply(function(b){
        b <- b %>%
          mutate(vert_speed_smooth = ifelse(row_number() <= sw_vs | row_number() > n() - sw_vs, NA,
                                            map_dbl((sw_vs + 1):(n() - sw_vs), ~ {
                                              i <- .x
                                              mean(b %>% slice((i - sw_vs):(i + sw_vs)) %>% pull(vert_speed), na.rm = T)
                                            })),
                 turn_angle_cum = ifelse(row_number() <= sw_th | row_number() > n() - sw_th, NA,
                                         map_dbl((sw_th + 1):(n() - sw_th), ~ {
                                           i <- .x
                                           max(abs(cumsum(b %>% slice((i - sw_th):(i + sw_th)) %>% drop_na(turn_angle) %>% pull(turn_angle))))
                                         })))
        return(b)
      }) %>%
      bind_rows()
  }) %>%
  bind_rows() %>%
  drop_na(vert_speed_smooth)

df_hr <- df_hr %>%
  mutate(kmean_cluster = kmeans(vert_speed_smooth, 2)$cluster)

soar_id <- df_hr %>%
  group_by(kmean_cluster) %>%
  summarise(mean_vert_speed_smooth = mean(vert_speed_smooth)) %>%
  ungroup() %>%
  filter(mean_vert_speed_smooth == max(mean_vert_speed_smooth)) %>%
  pull(kmean_cluster)

df_hr <- df_hr %>%
  mutate(soar_clust = if_else(kmean_cluster == soar_id, "soar", "glide"),
         flight_clust = if_else(is.na(turn_angle_cum), NA,
                                if_else(soar_clust == "glide", "gliding",
                                        if_else(turn_angle_cum >= circl_deg, "circular_soaring", "linear_soaring"))))

over_sea_classified <- df_hr %>%
  group_split(individual_local_identifier, burst_id) %>%
  lapply(function(b){

      b <- b %>%
        mutate(na_index = if_else(!is.na(flight_clust), NA,
                                  if_else(row_number() %in% 1:(sw_th-2), "head", "tail"))) %>%
        group_by(na_index) %>%
        mutate(duration = sum(time_lag_sec),
               turn_angle_sum = abs(sum(turn_angle, na.rm = T))) %>%
        mutate(flight_clust = if_else(!is.na(flight_clust), flight_clust,
                                      if_else(soar_clust == "glide", "gliding",
                                              if_else(soar_clust == "soar" & (turn_angle_sum/duration) >= ta_per_sec_circl, "circular_soaring", "linear_soaring")))) %>%
        ungroup() %>%
        dplyr::select(-c(duration, turn_angle_sum, na_index))

      b <- b %>%
        arrange(timestamp) %>%
        mutate(flight_clust_sm = if_else(row_number() %in% c(1, 2, nrow(.) - 1, nrow(.)), flight_clust,
                                         if_else(lag(flight_clust, 1) == "circular_soaring" & lead(flight_clust, 2) == "circular_soaring", "circular_soaring",
                                                 if_else(lag(flight_clust, 2) == "circular_soaring" & lead(flight_clust, 1) == "circular_soaring", "circular_soaring",
                                                         flight_clust)))) %>%
        mutate(flight_clust_diff = if_else(c(0, diff(as.numeric(as.factor(flight_clust_sm)))) == 0, 0, 1),
               flight_seg_id = cumsum(coalesce(flight_clust_diff,0)))

      behav_duration <- b %>%
        group_by(flight_seg_id) %>%
        reframe(timestamp_start = head(timestamp,1),
                duration = sum(time_lag_sec),
                flight_clust_sm = unique(flight_clust_sm)) %>%
        ungroup() %>%
        arrange(timestamp_start) %>%
        mutate(flight_clust_sm1 = flight_clust_sm) %>%
        mutate(flight_clust_sm1 = if_else(row_number() %in% c(1,nrow(.)), flight_clust_sm,
                                          if_else(lag(flight_clust_sm1, 1) == lead(flight_clust_sm1,1) & duration <= min_behav_d, lag(flight_clust_sm1), flight_clust_sm1)))

      b <- b %>%
        left_join(behav_duration %>% select(flight_seg_id, flight_clust_sm, flight_clust_sm1), by = c("flight_seg_id","flight_clust_sm")) %>%
        mutate(flight_clust_diff_sm1 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm1)))) == 0, 0, 1),
               flight_seg_id_sm1 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm1), NA,
                                           cumsum(coalesce(flight_clust_diff_sm1, 0))))

      behav_duration2 <- b %>%
        group_by(flight_seg_id_sm1) %>%
        reframe(timestamp_start = head(timestamp,1),
                duration = sum(time_lag_sec),
                flight_clust_sm1 = unique(flight_clust_sm1)) %>%
        ungroup() %>%
        arrange(timestamp_start) %>%
        mutate(flight_clust_sm2 = if_else(row_number() %in% c(1,nrow(.)), flight_clust_sm1,
                                          if_else(lag(flight_clust_sm1, 1) == lead(flight_clust_sm1,1) & duration <= min_behav_d, lag(flight_clust_sm1), flight_clust_sm1)))

      b <- b %>%
        left_join(behav_duration2 %>% select(flight_seg_id_sm1, flight_clust_sm1, flight_clust_sm2), by = c("flight_seg_id_sm1","flight_clust_sm1")) %>%
        mutate(flight_clust_diff_sm2 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm2)))) == 0, 0, 1),
               flight_seg_id_sm2 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm2), NA,
                                           cumsum(coalesce(flight_clust_diff_sm2, 0))))

      behav_ta <- b %>%
        group_by(flight_seg_id_sm2, flight_clust_sm2) %>%
        reframe(duration = sum(time_lag_sec),
                turn_angle_sum = abs(sum(turn_angle, na.rm = T))) %>%
        mutate(flight_clust_sm3 = if_else(flight_clust_sm2 == "circular_soaring" & (turn_angle_sum/duration) < ta_per_sec_circl, "linear_soaring",
                                          if_else(flight_clust_sm2 == "linear_soaring" & (turn_angle_sum/duration) >= ta_per_sec_circl, "circular_soaring",
                                                  flight_clust_sm2)))

      b <- b %>%
        left_join(behav_ta %>% select(flight_seg_id_sm2, flight_clust_sm2, flight_clust_sm3), by = c("flight_seg_id_sm2","flight_clust_sm2")) %>%
        mutate(flight_clust_diff_sm3 = if_else(c(NA, diff(as.numeric(as.factor(flight_clust_sm3)))) == 0, 0, 1),
               flight_seg_id_sm3 = if_else(row_number() == 1 | is.na(flight_clust_diff_sm3), NA,
                                           cumsum(coalesce(flight_clust_diff_sm3,0))))

      b <- b %>%
        mutate(track_flight_seg_id = paste(unique(b$individual_local_identifier), "burst", unique(b$burst_id), "seg", b$flight_seg_id_sm3, sep = "_"))

      return(b)
  }) %>%
  bind_rows()

saveRDS(over_sea_classified, file = "over_sea_classified.rds")

thermal_bouts <- over_sea_classified %>%
  filter(flight_clust_sm3 == "circular_soaring") %>%
  group_by(track_flight_seg_id) %>%
  mutate(bout_duration = sum(time_lag_sec)) %>%
  ungroup() %>%
  filter(bout_duration >= min_therm_d) %>%
  rename(thermal_id = track_flight_seg_id)

saveRDS(thermal_bouts, file = "thermal_bouts_over_sea.rds")

#------------------------------------------------------------------------#####
# STEP 5: calculate wind speed from data --------------------------------#####
#------------------------------------------------------------------------#####

# per thermal bout, following Weinzierl et al. 2016 (Ecology & Evolution, doi: 10.1002/ece3.2585)
# moveWindSpeed::getWindEstimate() implements their maximum-likelihood wind estimator.
# Pennycuick morphology constants for the honey buzzard are from the "buzzard" row of the
# multi-species MORPH table in 04_wind.R (buzzard_finland -> "buzzard": m, b, S, alt).

deg2rad <- pi/180
R_earth <- 6371000
g <- 9.81
PHI <- 0
MIN_FIX <- 8
WIND_MAX <- 40
GS_MIN <- 2
EXT_MIN <- 3

morph_m <- 0.754
morph_b <- 1.35
morph_S <- 0.22
morph_alt <- 1000
k_ind <- 1.1
Cdpro <- 0.014
Cdb <- 0.10
rho <- 1.225 * exp(-morph_alt/8500)

sinkrate <- function(V, phi) {
  n <- 1/cos(pmin(phi, 60*deg2rad)); L <- n*morph_m*g; Sb <- 0.00813*morph_m^0.666
  Di <- k_ind*L^2/(0.5*rho*V^2*pi*morph_b^2); Dpar <- 0.5*rho*V^2*Sb*Cdb; Dpro <- 0.5*rho*V^2*morph_S*Cdpro
  V*(Di + Dpar + Dpro)/(morph_m*g)
}

thermal_bouts_xy <- thermal_bouts %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  group_by(thermal_id) %>%
  arrange(timestamp, .by_group = TRUE) %>%
  mutate(lon0 = first(lon), lat0 = first(lat),
         x_m = (lon - lon0) * deg2rad * R_earth * cos(lat0 * deg2rad),
         y_m = (lat - lat0) * deg2rad * R_earth,
         ts_n = as.numeric(timestamp),
         dt = ts_n - lag(ts_n),
         vx = (x_m - lag(x_m))/dt,
         vy = (y_m - lag(y_m))/dt) %>%
  mutate(vx = if_else(!is.finite(dt) | dt <= 0 | dt > 2, NA_real_, vx),
         vy = if_else(!is.finite(dt) | dt <= 0 | dt > 2, NA_real_, vy)) %>%
  ungroup()

wind_est <- thermal_bouts_xy %>%
  group_by(thermal_id) %>%
  group_modify(~{
    d <- .x
    G <- as.matrix(d[is.finite(d$vx) & is.finite(d$vy), c("vx","vy")])
    gs_mean <- mean(d$ground_speed, na.rm = TRUE)
    ext <- max(diff(range(d$x_m, na.rm = TRUE)), diff(range(d$y_m, na.rm = TRUE)))
    if(nrow(G) < MIN_FIX || !is.finite(gs_mean) || gs_mean < GS_MIN || !is.finite(ext) || ext < EXT_MIN) return(tibble())
    e <- tryCatch(moveWindSpeed::getWindEstimate(G, PHI), error = function(z) NULL)
    if(is.null(e) || !all(is.finite(e$windEst)) || sqrt(sum(e$windEst^2)) > WIND_MAX) return(tibble())
    air <- sqrt((G[,1]-e$windEst[1])^2 + (G[,2]-e$windEst[2])^2)
    airspeed <- mean(air)
    if(!is.finite(airspeed) || airspeed > WIND_MAX) return(tibble())
    tibble(u_wind = e$windEst[1], v_wind = e$windEst[2], airspeed = airspeed, air_sd = sd(air), n_fix_wind = nrow(G))
  }) %>%
  ungroup()

bout_summary <- thermal_bouts_xy %>%
  group_by(thermal_id, individual_local_identifier) %>%
  summarise(bout_duration = first(bout_duration),
            turn_angle_sum = sum(abs(turn_angle), na.rm = TRUE),
            climb = (last(height_msl) - first(height_msl))/first(bout_duration),
            .groups = "drop") %>%
  mutate(omega = (turn_angle_sum * deg2rad) / bout_duration)

bout_wind <- bout_summary %>%
  inner_join(wind_est, by = "thermal_id") %>%
  mutate(wind_speed = sqrt(u_wind^2 + v_wind^2),
         wind_dir = (atan2(v_wind, u_wind) * 180/pi) %% 360,
         bank = atan(airspeed * omega / g),
         load_factor = 1/cos(bank),
         sink_rate = sinkrate(airspeed, bank),
         updraft = climb + sink_rate) %>%
  select(-individual_local_identifier, -bout_duration)

thermal_wind <- thermal_bouts %>%
  left_join(bout_wind, by = "thermal_id")

saveRDS(thermal_wind, file = "thermal_wind_estimates_over_sea.rds")


#sanity check


