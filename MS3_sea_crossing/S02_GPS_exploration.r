#exploration of GPS data for soaring behavior of honey buzzards over the open sea
#5.11.2026, Lausanne, CH
#elham.nourani@unil.ch

library(tidyverse)
library(sf)
library(mgcv)
library(mapview)
library(patchwork)

setwd("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/R_files/")

fig_dir <- "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/conferences/Wynn_Nourani_symposium"
pres_text_size <- 16

#------- STEP 1: Maps ---------------------------------------------------#####
#one map showing all the trajectories to highlight the propensity of sea-crossing in different age groups. age group should be year since tagging since everyone was tagged in the nest.
#One panel per age group

gap_break_hr <- 48 #tracking gaps longer than this are not connected with a line segment, to avoid drawing a straight "teleport" line across a real data gap

all_gps <- readRDS("cleaned_gps_with_altitude_Oct26.rds") %>%
  mutate(age_group = factor(floor(as.numeric(difftime(timestamp, deploy_on_timestamp, units = "days")) / 365.25))) %>%
  arrange(individual_local_identifier, timestamp) %>%
  group_by(individual_local_identifier) %>%
  mutate(gap_hr = as.numeric(difftime(timestamp, lag(timestamp), units = "hours")),
         track_seg = cumsum(row_number() == 1 | coalesce(gap_hr > gap_break_hr, TRUE) | age_group != lag(age_group))) %>%
  ungroup()

age_tracks <- all_gps %>%
  group_by(individual_local_identifier, age_group, track_seg) %>%
  filter(n() >= 2) %>%
  summarise(do_union = FALSE, .groups = "drop") %>%
  st_cast("LINESTRING")

all_trajectories_map <- ggplot() +
  geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", scale = 50), linewidth = 0.3) +
  geom_sf(data = age_tracks, aes(color = individual_local_identifier, group = interaction(individual_local_identifier, track_seg)), linewidth = 0.6) +
  coord_sf(xlim = c(-17, 53), ylim = c(-35.6, 70)) +
  scale_x_continuous(breaks = seq(-20, 50, by = 20)) +
  facet_wrap(~ age_group, nrow = 1, labeller = labeller(age_group = ~ paste("Year", .x, "since tracking"))) +
  labs(color = "Individual", title = "All trajectories by age since tagging") +
  theme_linedraw(base_size = pres_text_size) +
  theme(legend.position = "none", axis.text = element_text(size = rel(0.55)))

ggsave(file.path(fig_dir, "all_trajectories_by_age_group.png"), all_trajectories_map, width = 22, height = 6, dpi = 300)

#map 2 the same style, but only with the focal individual "D225_236". color the points over the mediterranean and sahara based on the colors that I've chosen below.

gps_classified <- readRDS("gps_classified_all_Oct26.rds")

gis_dir <- "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/GIS_files"

landscape_zones <- st_read(file.path(gis_dir, "world_oceans/goas_v01.shp")) %>%
  filter(name %in% c("Mediterranean Region", "Baltic Sea")) %>%
  transmute(zone = recode(name, "Mediterranean Region" = "Mediterranean", "Baltic Sea" = "Baltic")) %>%
  bind_rows(
    st_read(file.path(gis_dir, "ne_geography_regions/ne_10m_geography_regions_polys.shp")) %>%
      filter(NAME == "SAHARA") %>%
      transmute(zone = "Sahara")
  ) %>%
  mutate(zone = factor(zone, levels = c("Sahara", "Mediterranean", "Baltic")))

gps_classified <- gps_classified %>%
  st_join(landscape_zones)

one <- gps_classified %>%
  filter(individual_local_identifier == "D225_236")

one_bbox <- st_bbox(one)
zoom_pad_lon <- 6
zoom_pad_lat <- 1.5

focal_map <- ggplot() +
  geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", scale = 50), linewidth = 0.3) +
  geom_sf(data = one, aes(color = zone), size = 1, alpha = 0.6) +
  coord_sf(xlim = c(one_bbox["xmin"] - zoom_pad_lon, one_bbox["xmax"] + zoom_pad_lon),
           ylim = c(one_bbox["ymin"] - zoom_pad_lat, one_bbox["ymax"] + zoom_pad_lat)) +
  scale_x_continuous(breaks = seq(10, 40, by = 10)) +
  scale_color_manual(values = c("Baltic" = "#756bb1", "Mediterranean" = "#2b8cbe", "Sahara" = "#fdae6b"), na.value = "grey60",
                     breaks = c("Sahara", "Mediterranean", "Baltic")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(color = "", title = "D225_236") +
  theme_linedraw(base_size = pres_text_size)

ggsave(file.path(fig_dir, "D225_236_trajectory_by_zone.png"), focal_map, width = 10, height = 8, dpi = 300)

#------- STEP 2: Flight behavior exploration ---------------------------------------------------#####

#plot climb rate as a function of latitude, highlight latitudinal zones that correspond to the Baltic, the Mediterranean, and the Sahara desert

sample_fixes <- one %>%
  st_drop_geometry()

vs_outlier_cap <- 15

sample_fixes_vs <- sample_fixes %>%
  filter(time_lag_sec == 1, vert_speed > 0, vert_speed <= vs_outlier_cap)

lat_axis_limits <- c(10, max(sample_fixes$lat, na.rm = TRUE))

vert_speed_plot <- ggplot() +
  geom_point(data = sample_fixes_vs, aes(x = lat, y = vert_speed, color = zone), size = 0.5, alpha = 0.3) +
  #geom_smooth(data = sample_fixes_vs, aes(x = lat, y = vert_speed), color = "black", se = FALSE) +
  scale_color_manual(values = c("Baltic" = "#756bb1", "Mediterranean" = "#2b8cbe", "Sahara" = "#fdae6b"), na.value = "grey80",
                     breaks = c("Sahara", "Mediterranean", "Baltic")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_x_reverse(limits = rev(lat_axis_limits)) +
  labs(x = "latitude", y = "vertical speed (m/s)", color = "", title = "D225_236: vertical speed across latitude") +
  theme_minimal(base_size = pres_text_size)

ggsave(file.path(fig_dir, "D225_236_vertical_speed_by_latitude.png"), vert_speed_plot, width = 10, height = 7, dpi = 300)

#plot flight altitude above sea level as a function of latitude, colored by the zone each point falls over

flight_altitude_plot <- ggplot() +
  geom_point(data = sample_fixes, aes(x = lat, y = flight_height, color = zone), size = 0.5, alpha = 0.15) +
  scale_color_manual(values = c("Baltic" = "#756bb1", "Mediterranean" = "#2b8cbe", "Sahara" = "#fdae6b"), na.value = "grey80",
                     breaks = c("Sahara", "Mediterranean", "Baltic")) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  scale_x_reverse(limits = rev(lat_axis_limits)) +
  labs(x = "latitude", y = "height above ground or sea (m)", color = "", title = "D225_236: flight altitude across latitude") +
  theme_minimal(base_size = pres_text_size)

ggsave(file.path(fig_dir, "D225_236_flight_altitude_by_latitude.png"), flight_altitude_plot, width = 10, height = 7, dpi = 300)

#------- STEP 3: Flight ecology ---------------------------------------------------#####

#gam to predict climb rate as a function of delta-T and wind (potentially side-wind)

over_sea_classified <- readRDS("over_sea_classified.rds")

annotations <- readRDS("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/R_files/sea_gps_seg_ann.rds") %>%
  st_drop_geometry()

over_sea_ann <- over_sea_classified %>%
  st_drop_geometry() %>%
  filter(flight_clust == "circular_soaring") %>%
  left_join(annotations %>% dplyr::select(c(local_identifier, timestamp, delta_t, u_10m, v_10m, wind_speed_10m)),
            by = c("individual_local_identifier" = "local_identifier", "timestamp")) %>%
  drop_na(delta_t, u_10m, v_10m) #keep only bird-years that have corresponding env annotations

#wind support (tailwind component) and cross-wind (side-wind component) relative to flight direction
heading_rad <- over_sea_ann$flight_direction * pi/180

over_sea_ann <- over_sea_ann %>%
  mutate(wind_support = u_10m * sin(heading_rad) + v_10m * cos(heading_rad),
         cross_wind = u_10m * cos(heading_rad) - v_10m * sin(heading_rad))

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
png(file.path(fig_dir, "gam_climb_rate_smooths.png"), width = 1600, height = 1200, res = 150)
plot(m_climb, pages = 1, shade = TRUE, cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.6)
dev.off()

summary(m_alt)
png(file.path(fig_dir, "gam_height_agl_smooths.png"), width = 1600, height = 1200, res = 150)
plot(m_alt, pages = 1, shade = TRUE, cex.lab = 1.6, cex.axis = 1.3, cex.main = 1.6)
dev.off()

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

heatmap_plot <- ggplot(pred_grid, aes(x = delta_t, y = cross_wind, fill = climb_rate)) +
  geom_raster() +
  scale_fill_viridis_c(name = "predicted\nclimb rate (m/s)") +
  labs(x = expression(Delta*T~(degree*C)), y = "cross-wind (m/s)",
       title = "Predicted climb rate: delta-T x cross-wind interaction") +
  theme_minimal(base_size = pres_text_size)

ggsave(file.path(fig_dir, "predicted_climb_rate_heatmap.png"), heatmap_plot, width = 9, height = 7, dpi = 300)

#------- STEP 4: Flight behavior in detail ---------------------------------------------------#####

#zoom into one part of D225_236's trajectory and look at flight height and climb rate (vertical speed)

vs_color_lim <- 5

vs_color_scale <- scale_color_gradient2(name = "vertical\nspeed (m/s)", low = "tomato", mid = "white", high = "royalblue", midpoint = 0,
                                         limits = c(-vs_color_lim, vs_color_lim), oob = scales::squish,
                                         breaks = seq(-vs_color_lim, vs_color_lim, 2),
                                         guide = guide_colorbar(barheight = unit(6, "cm"), barwidth = unit(0.6, "cm")))

chunk <- sample_fixes %>%
  filter(month(timestamp) == 10, day(timestamp) == 3,
         between(hour(timestamp)*60 + minute(timestamp), 14*60+30, 15*60+30)) %>%
  arrange(timestamp) %>%
  mutate(seg_change = flight_clust != lag(flight_clust) | row_number() == 1,
         seg_id = cumsum(coalesce(seg_change, TRUE))) %>%
  group_by(seg_id) %>%
  mutate(seg_duration = as.numeric(difftime(max(timestamp), min(timestamp), units = "secs")),
         cum_turn_deg = cumsum(replace_na(abs(turn_angle), 0))) %>%
  ungroup()

height_time_plot <- ggplot(chunk, aes(x = timestamp, y = height_msl, group = seg_id)) +
  geom_path(aes(color = vert_speed)) +
  geom_point(aes(color = vert_speed), size = 1) +
  vs_color_scale +
  labs(x = "timestamp", y = "height above mean sea level (m)", title = "D225_236: profile view, Oct 3 14:30-15:30") +
  theme_linedraw(base_size = pres_text_size)

ggsave(file.path(fig_dir, "D225_236_height_vs_time_chunk.png"), height_time_plot, width = 10, height = 7, dpi = 300)

#same chunk, top-down view: longitude vs latitude, same vert_speed coloring as the height-vs-time plot

topdown_plot <- ggplot(chunk, aes(x = lon, y = lat, group = seg_id)) +
  geom_path(aes(color = vert_speed)) +
  geom_point(aes(color = vert_speed), size = 1) +
  vs_color_scale +
  labs(x = "longitude", y = "latitude", title = "D225_236: top-down view, Oct 3 14:30-15:30") +
  theme_linedraw(base_size = pres_text_size)

ggsave(file.path(fig_dir, "D225_236_topdown_chunk.png"), topdown_plot, width = 9, height = 8, dpi = 300)

#same top-down view, zoomed in to 15:00-15:02

chunk_zoom <- chunk %>%
  filter(between(hour(timestamp)*60 + minute(timestamp), 15*60+0, 15*60+2))

topdown_zoom_plot <- ggplot(chunk_zoom, aes(x = lon, y = lat, group = seg_id)) +
  geom_path(aes(color = vert_speed)) +
  geom_point(aes(color = vert_speed), size = 1) +
  vs_color_scale +
  labs(x = "longitude", y = "latitude", title = "D225_236: top-down view, zoomed (15:00-15:02)") +
  theme_linedraw(base_size = pres_text_size)

ggsave(file.path(fig_dir, "D225_236_topdown_chunk_zoom.png"), topdown_zoom_plot, width = 9, height = 8, dpi = 300)

#same chunk, top-down views colored by wind conditions: cross-wind, wind support, wind speed, and delta-T

annotations <- readRDS("sea_gps_seg_ann.rds") %>%
  st_drop_geometry()

chunk_wind <- chunk %>%
  left_join(annotations %>% select(local_identifier, timestamp, delta_t, u_10m, v_10m, wind_speed_10m),
            by = c("individual_local_identifier" = "local_identifier", "timestamp"))

heading_rad <- chunk_wind$flight_direction * pi/180

chunk_wind <- chunk_wind %>%
  mutate(wind_support = u_10m * sin(heading_rad) + v_10m * cos(heading_rad),
         cross_wind = u_10m * cos(heading_rad) - v_10m * sin(heading_rad))

wind_topdown_plot <- function(data, var, legend_name, palette = "viridis", legend_position = "right") {
  ggplot(data, aes(x = lon, y = lat, group = seg_id)) +
    geom_path(aes(color = .data[[var]])) +
    geom_point(aes(color = .data[[var]]), size = 2) +
    scale_color_viridis_c(name = legend_name, option = palette,
                          guide = guide_colorbar(barwidth = unit(if (legend_position == "bottom") 6 else 0.6, "cm"),
                                                 barheight = unit(if (legend_position == "bottom") 0.6 else 6, "cm"))) +
    labs(x = "longitude", y = "latitude", title = legend_name) +
    theme_linedraw(base_size = pres_text_size) +
    theme(legend.position = legend_position, legend.title = element_blank())
}

wind_topdown_panels <- (wind_topdown_plot(chunk_wind, "cross_wind", "cross-wind (m/s)") + wind_topdown_plot(chunk_wind, "wind_support", "wind support (m/s)") +
                           wind_topdown_plot(chunk_wind, "wind_speed_10m", "wind speed (m/s)") + wind_topdown_plot(chunk_wind, "delta_t", "delta-T (°C)")) +
  plot_layout(ncol = 2) +
  plot_annotation(title = "D225_236: top-down views by wind conditions, Oct 3 14:30-15:30",
                   theme = theme(plot.title = element_text(size = pres_text_size)))

ggsave(file.path(fig_dir, "D225_236_topdown_chunk_wind.png"), wind_topdown_panels, width = 14, height = 12, dpi = 300)

#same wind conditions, zoomed in to 15:00-15:01: cross-wind, wind support, and raw heading side by side

chunk_zoom_wind <- chunk_wind %>%
  filter(between(hour(timestamp)*60 + minute(timestamp), 15*60+0, 15*60+1))

wind_topdown_zoom_panels <- (wind_topdown_plot(chunk_zoom_wind, "flight_direction", "heading (°)", palette = "plasma", legend_position = "bottom") +
                                wind_topdown_plot(chunk_zoom_wind, "cross_wind", "cross-wind (m/s)", legend_position = "bottom") +
                                wind_topdown_plot(chunk_zoom_wind, "wind_support", "wind support (m/s)", legend_position = "bottom")) +
  plot_layout(ncol = 3) +
  plot_annotation(title = "D225_236: top-down view by wind conditions, zoomed (15:00-15:01)",
                   theme = theme(plot.title = element_text(size = pres_text_size)))

ggsave(file.path(fig_dir, "D225_236_topdown_chunk_zoom_wind.png"), wind_topdown_zoom_panels, width = 20, height = 8, dpi = 300)
