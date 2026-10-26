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
library(plotly)

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

# net-360 signed-rotation + validity-gate + gap-based thermal grouping,
# ported from 03_segment.R (laterality/thermalling project). Solar-day tagging and the
# glide-segment extraction from that script are unrelated to thermal detection and are omitted.
# One deliberate deviation: 03_segment.R projects local x/y from the bird's first-ever fix; we
# only have the over-sea subset of the track, so the origin here is the first over-sea fix instead.

min_res_tl <- 2
sw_tr <- 2
turn_lo <- 4
glide_k <- 5L
therm_gap <- 250
min_therm_d <- 20
ext_hard <- 3
ext_slow <- 10
av_slow <- 18
net_lo <- 300
gs_min <- 2

#accumulates signed turn angle per fix to detect completed circles (net 360 deg); also records why an
#incomplete rotation ended (glide-reset vs ran out of fixes), to make "other" fixes diagnosable
label_circles <- function(turn, str) {
  n <- length(turn)
  cyc <- rep(NA_integer_, n)
  rotating <- rep(FALSE, n)
  acc <- 0; straight <- 0L; k <- 0L; completed <- integer(0); reset_ks <- integer(0)
  for (i in seq_len(n)) {
    rot <- is.finite(str[i]) && str[i] >= turn_lo
    if (!rot) {
      straight <- straight + 1L
      if (straight == glide_k) { reset_ks <- c(reset_ks, k); acc <- 0; k <- k + 1L }
    } else {
      straight <- 0L
      rotating[i] <- TRUE
      acc <- acc + (if (is.finite(turn[i])) turn[i] else 0)
      cyc[i] <- k
      if (abs(acc) >= 360) { completed <- c(completed, k); k <- k + 1L; acc <- acc - sign(acc) * 360 }
    }
  }
  behav <- if_else(!is.na(cyc) & cyc %in% completed, "circling", if_else(rotating, "other", "glide"))
  other_reason <- if_else(behav != "other", NA_character_,
                          if_else(cyc %in% reset_ks, "incomplete_glide_reset", "incomplete_burst_end"))
  tibble(behav = behav, cyc = cyc, other_reason = other_reason)
}

#linearly detrends a position vector against time (removes net drift), closed-form least squares
detr <- function(v, tc, sst) if (sst <= 0) v - mean(v) else (v - mean(v)) - (sum(tc * (v - mean(v))) / sst) * tc

df_hr <- over_sea %>%
  filter(individual_local_identifier == "D225_236") %>%
  mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
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
             azimuth_positive = if_else(azimuth >= 0, azimuth, azimuth + 360),
             lon0 = first(lon), lat0 = first(lat),
             x_m = (lon - lon0) * (pi/180) * 6371000 * cos(lat0 * pi/180),
             y_m = (lat - lat0) * (pi/180) * 6371000) %>%
      group_by(burst_id) %>%
      mutate(altitude_diff = if_else(row_number() == 1, NA, height_msl - lag(height_msl,1)),
             vert_speed = if_else(row_number() == 1, NA, altitude_diff/time_lag_sec),
             turn_angle = if_else(row_number() == 1 | row_number() == n(), NA,
                                  ((lag(azimuth_positive) - azimuth_positive + 180) %% 360) - 180),
             ground_speed = if_else(row_number() == 1, NA, step_length/time_lag_sec)) %>%
      ungroup()

    ind %>%
      group_split(burst_id) %>%
      lapply(function(b){
        if(nrow(b) < 3) {
          b$behav <- "other"; b$cyc <- NA_real_; b$other_reason <- "burst_too_short"
          return(b)
        }
        b <- b %>%
          mutate(turn_rate = abs(turn_angle)/pmax(time_lag_sec, 1))
        b <- b %>%
          mutate(turn_rate_smooth = ifelse(row_number() <= sw_tr | row_number() > n() - sw_tr, NA,
                                           map_dbl((sw_tr + 1):(n() - sw_tr), ~ {
                                             i <- .x
                                             mean(b %>% slice((i - sw_tr):(i + sw_tr)) %>% pull(turn_rate), na.rm = T)
                                           })))
        lab <- label_circles(b$turn_angle, b$turn_rate_smooth)
        b$behav <- lab$behav
        b$other_reason <- lab$other_reason
        b$cyc <- if_else(is.na(lab$cyc), NA_real_, lab$cyc + unique(b$burst_id) * 100000)
        return(b)
      }) %>%
      bind_rows()
  }) %>%
  bind_rows()

#per-circle summary (every attempted rotation, completed or not) + the geometric validity gate
circles <- df_hr %>%
  st_drop_geometry() %>%
  filter(!is.na(cyc)) %>%
  group_by(individual_local_identifier, cyc) %>%
  group_modify(~{
    d <- .x
    ts_n <- as.numeric(d$timestamp)
    tc <- ts_n - mean(ts_n)
    sst <- sum(tc^2)
    px <- detr(d$x_m, tc, sst)
    py <- detr(d$y_m, tc, sst)
    tibble(burst_id = d$burst_id[1],
           behav = d$behav[1],
           t_start = min(d$timestamp), t_end = max(d$timestamp),
           dur = as.numeric(difftime(max(d$timestamp), min(d$timestamp), units = "secs")),
           net = sum(d$turn_angle, na.rm = TRUE),
           alt0 = first(d$height_msl), alt1 = last(d$height_msl),
           gs_mean = mean(d$ground_speed, na.rm = TRUE),
           ext = max(diff(range(px)), diff(range(py))))
  }) %>%
  ungroup() %>%
  mutate(ang_vel = abs(net)/dur,
         valid = behav == "circling" & !(ext < ext_hard) & !(ext < ext_slow & ang_vel < av_slow) &
           abs(net) >= net_lo & gs_mean >= gs_min)

#valid circles within therm_gap (and the same burst) of each other form one thermal; climb/bout_net_climb
#follow 03_segment.R exactly (duration-weighted mean of per-circle climb, excluding inter-circle gap time)
valid_circles <- circles %>%
  filter(valid) %>%
  arrange(individual_local_identifier, t_start) %>%
  group_by(individual_local_identifier) %>%
  mutate(gap_prev = as.numeric(difftime(t_start, lag(t_end), units = "secs")),
         new_therm = is.na(gap_prev) | gap_prev > therm_gap | burst_id != lag(burst_id),
         therm_idx = cumsum(new_therm),
         track_flight_seg_id = sprintf("%s_t%04d", individual_local_identifier, therm_idx)) %>%
  ungroup() %>%
  group_by(track_flight_seg_id) %>%
  mutate(nth_circle_in_thermal = row_number(),
         climb = (alt1 - alt0)/dur,
         bout_net_climb = sum(climb * dur)/sum(dur)) %>%
  ungroup() %>%
  mutate(circle_id = sprintf("%s_c%02d", track_flight_seg_id, nth_circle_in_thermal)) %>%
  select(individual_local_identifier, cyc, track_flight_seg_id, circle_id, climb, bout_net_climb)

over_sea_classified <- df_hr %>%
  left_join(circles %>% select(individual_local_identifier, cyc, circle_net = net, circle_dur = dur,
                                circle_ang_vel = ang_vel, circle_ext = ext, circle_gs_mean = gs_mean,
                                circle_valid = valid),
            by = c("individual_local_identifier", "cyc")) %>%
  left_join(valid_circles, by = c("individual_local_identifier", "cyc")) %>%
  mutate(flight_clust = case_when(
    behav == "circling" & circle_valid ~ "circular_soaring",
    behav == "glide" ~ "gliding",
    TRUE ~ "other"
  ),
  other_reason = case_when(
    flight_clust != "other" ~ NA_character_,
    behav == "circling" & !circle_valid ~ "completed_invalid_geometry",
    TRUE ~ other_reason
  ))


#sanity check
mapview(over_sea_classified, zcol ="flight_clust", alpha= 0)
mapview(over_sea_classified, zcol ="circle_net", alpha= 0)
mapview(over_sea_classified, zcol ="other_reason", alpha= 0)

saveRDS(over_sea_classified, file = "over_sea_classified.rds")

thermal_bouts <- over_sea_classified %>%
  filter(flight_clust == "circular_soaring") %>%
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

#Pennycuick circling sink rate for honey buzzard morphology, given airspeed and bank angle
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
            climb = first(bout_net_climb),
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

#3D trajectory colored by flight type, local planar coords (undistorted, unlike lon/lat)
traj_3d <- over_sea_classified %>%
  filter(individual_local_identifier == "D225_236" & timestamp >= as.Date("2022-10-03")) %>%
  st_drop_geometry()

plot_ly(traj_3d, x = ~x_m, y = ~y_m, z = ~height_msl, color = ~height_msl,
        colors = c("firebrick", "steelblue", "grey70"),
        type = "scatter3d", mode = "lines+markers",
        marker = list(size = 2), line = list(width = 2)) %>%
  layout(scene = list(aspectmode = "data",
                       xaxis = list(title = "x (m)"),
                       yaxis = list(title = "y (m)"),
                       zaxis = list(title = "height MSL (m)")))

#zoom in on one thermal bout to see the circling helix, colored by wind speed
one_thermal <- thermal_wind %>%
  filter(thermal_id == first(thermal_id)) %>%
  st_drop_geometry()

plot_ly(one_thermal, x = ~x_m, y = ~y_m, z = ~height_msl, color = ~wind_speed,
        type = "scatter3d", mode = "lines+markers",
        marker = list(size = 3), line = list(width = 2)) %>%
  layout(scene = list(aspectmode = "data",
                       xaxis = list(title = "x (m)"),
                       yaxis = list(title = "y (m)"),
                       zaxis = list(title = "height MSL (m)")))

#static height-vs-latitude profile for the Oct 3rd 14:30-15:30 window
chunk <- over_sea_classified %>%
  st_drop_geometry() %>%
  filter(month(timestamp) == 10, day(timestamp) == 3,
         between(hour(timestamp)*60 + minute(timestamp), 14*60+30, 15*60+30)) %>%
  arrange(timestamp) %>%
  mutate(seg_change = flight_clust != lag(flight_clust) | row_number() == 1,
         seg_id = cumsum(coalesce(seg_change, TRUE))) %>%
  group_by(seg_id) %>%
  mutate(seg_duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "secs"))) %>%
  ungroup()

ggplot(chunk, aes(x = timestamp, y = height_msl, color = turn_angle, group = seg_id)) +
  geom_path() +
  geom_point(size = 1) +
  scale_color_gradient2(name = "turn angle\n(deg)", low = "firebrick", mid = "grey85", high = "steelblue", midpoint = 0) +
  labs(x = "timestamp", y = "height above MSL (m)") +
  theme_linedraw()

