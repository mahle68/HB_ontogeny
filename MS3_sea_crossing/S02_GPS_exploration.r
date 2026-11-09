#exploration of GPS data for soaring behavior of honey buzzards over the open sea
#5.11.2026, Lausanne, CH
#elham.nourani@unil.ch

library(tidyverse)
library(sf)
library(mgcv)
library(mapview)

setwd("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/R_files/")

#------- STEP 1: Flight behavior exploration ---------------------------------------------------#####

#plot climb rate as a function of latitude, highlight latitudinal zones that correspond to the Baltic, the Mediterranean, and the Sahara desert

gps_classified <- readRDS("gps_classified_all_Oct26.rds")

#not restricted to circling bouts: the current circling definition in S01 misses some non-circling soaring, so use the per-fix vertical speed across all flight types instead of the bout-level climb rate for now
one <- gps_classified %>%
  filter(individual_local_identifier == "D225_236")


sample_fixes <- gps_classified %>%
  st_drop_geometry() %>%
  filter(individual_local_identifier == "D225_236")

#GPS bursts (1 Hz fixes) form burst_id groups with multiple rows; isolated regular tracking fixes are single-row groups
burst_sizes <- sample_fixes %>%
  count(burst_id, name = "burst_n")

sample_fixes_vs <- sample_fixes %>%
  left_join(burst_sizes, by = "burst_id") %>%
  filter(burst_n >= 3, vert_speed > 0)

#approximate latitudinal bands for the three barriers along this bird's route
lat_zones <- tibble(zone = factor(c("Sahara", "Mediterranean", "Baltic"), levels = c("Sahara", "Mediterranean", "Baltic")),
                     xmin = c(15, 30, 53.5),
                     xmax = c(30, 46.49, 61.05))

ggplot() +
  geom_rect(data = lat_zones, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = zone), alpha = 0.15) +
  geom_point(data = sample_fixes_vs, aes(x = lat, y = vert_speed), size = 0.5, alpha = 0.15) +
  geom_smooth(data = sample_fixes_vs, aes(x = lat, y = vert_speed), color = "black", se = FALSE) +
  scale_fill_manual(values = c("Baltic" = "steelblue", "Mediterranean" = "darkorange", "Sahara" = "goldenrod")) +
  scale_x_reverse() +
  labs(x = "latitude", y = "vertical speed (m/s)", fill = "", title = "D225_236: vertical speed across latitude") +
  theme_minimal()

#plot flight altitude above sea level as a function of latitude, highlight latitudinal zones that correspond to the Baltic, the Mediterranean, and the Sahara desert

ggplot() +
  geom_rect(data = lat_zones, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = zone), alpha = 0.15) +
  geom_point(data = sample_fixes, aes(x = lat, y = height_agl), size = 0.5, alpha = 0.15) +
  geom_smooth(data = sample_fixes, aes(x = lat, y = height_agl), color = "black", se = FALSE) +
  scale_fill_manual(values = c("Baltic" = "steelblue", "Mediterranean" = "darkorange", "Sahara" = "goldenrod")) +
  scale_x_reverse() +
  labs(x = "latitude", y = "height above ground (m)", fill = "", title = "D225_236: flight altitude across latitude") +
  theme_minimal()

#------- STEP 2: Flight ecology ---------------------------------------------------#####

#gam to predict climb rate as a function of delta-T and wind (potentially side-wind)

over_sea_classified <- readRDS("over_sea_classified.rds")

annotations <- readRDS("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds") %>%
  st_drop_geometry()

over_sea_ann <- over_sea_classified %>%
  st_drop_geometry() %>%
  filter(flight_clust == "circular_soaring") %>%
  left_join(annotations %>% select(c(local_identifier, timestamp, delta_t, u_10m, v_10m, wind_speed_10m)),
            by = c("individual_local_identifier" = "local_identifier", "timestamp")) %>%
  drop_na(delta_t, u_10m, v_10m) #keep only bird-years that have corresponding env annotations

#wind support (tailwind component) and cross-wind (side-wind component) relative to flight direction
heading_rad <- over_sea_ann$flight_direction * pi/180

over_sea_ann <- over_sea_ann %>%
  mutate(wind_support = u_10m * sin(heading_rad) + v_10m * cos(heading_rad),
         cross_wind = -u_10m * cos(heading_rad) + v_10m * sin(heading_rad))

#one row per thermal bout: mean env conditions over the bout, paired with its climb rate and altitude
bout_env <- over_sea_ann %>%
  group_by(individual_local_identifier, track_flight_seg_id) %>%
  summarise(climb_rate = first(bout_net_climb),
            height_agl = mean(height_agl, na.rm = TRUE),
            delta_t = mean(delta_t, na.rm = TRUE),
            wind_support = mean(wind_support, na.rm = TRUE),
            cross_wind = mean(cross_wind, na.rm = TRUE),
            .groups = "drop")

#te() rather than s(x,y): delta-T and cross-wind are on very different scales, which te()'s separate marginal bases handle better than an isotropic s()
m_climb <- gam(climb_rate ~ te(delta_t, cross_wind) + s(wind_support) +
                 s(individual_local_identifier, bs = "re"),
               data = bout_env)

m_alt <- gam(height_agl ~ te(delta_t, cross_wind) + s(wind_support) +
               s(individual_local_identifier, bs = "re"),
             data = bout_env)

summary(m_climb)
plot(m_climb, pages = 1, shade = TRUE)

summary(m_alt)
plot(m_alt, pages = 1, shade = TRUE)

#presentation heatmap: predicted climb rate over the delta-T x cross-wind surface, at the mean wind support
#and with the individual random effect excluded (population-level prediction, not bird-specific)
pred_grid <- expand.grid(
  delta_t = seq(min(bout_env$delta_t), max(bout_env$delta_t), length.out = 100),
  cross_wind = seq(min(bout_env$cross_wind), max(bout_env$cross_wind), length.out = 100)
) %>%
  mutate(wind_support = mean(bout_env$wind_support, na.rm = TRUE),
         individual_local_identifier = bout_env$individual_local_identifier[1]) #arbitrary level; excluded below

pred_grid$climb_rate <- predict(m_climb, newdata = pred_grid, type = "response",
                                 exclude = "s(individual_local_identifier)")

ggplot(pred_grid, aes(x = delta_t, y = cross_wind, fill = climb_rate)) +
  geom_raster() +
  scale_fill_viridis_c(name = "predicted\nclimb rate (m/s)") +
  labs(x = expression(Delta*T~(degree*C)), y = "cross-wind (m/s)",
       title = "Predicted climb rate: delta-T x cross-wind interaction") +
  theme_minimal()
