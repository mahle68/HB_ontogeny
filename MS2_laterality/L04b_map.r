# script for making the map for Safi et al 2025
# Elham Nouani, PhD. 
# 11.12.2024, Konstanz, DE

library(tidyverse)
library(sf)
library(mapview)
library(ggridges)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

lat_zones <- seq(-30,65, by = 5)

#--------------------------------------------------------------------
# STEP 1: map tracking points ----------------------------------
#--------------------------------------------------------------------


# data prep
#open metadata to extract deployment date (the tracking data has some undeployed points in it)
deployment <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/EHB_metadata - Sheet1.csv") %>% 
  filter(nest_country == "Finland") %>% 
  mutate(deployment_dt_utc = as.POSIXct(deployment_dt_utc)) %>% 
  dplyr::select(ring_ID, deployment_dt_utc) %>% 
  rename(individual_local_identifier = ring_ID)

#life-cycle stages from L03a_tests_per_day.r AND append deployment
life_cycle <- readRDS("updated_life_cycle_nov24.rds") %>% 
  full_join(deployment, by = "individual_local_identifier")

#rewrite the life cycle file, to include deployment
saveRDS(life_cycle, file = "updated_life_cycle_nov24.rds")


data <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24.rds") 

cleaned_gps <- data %>% 
  #remove the points at (0,0) 
  filter(!(location_lat == 0 & location_long == 0)) %>% 
  mutate(unique_date = as.Date(timestamp)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, deployment_dt_utc, first_exploration, migration_start, migration_end), by = "individual_local_identifier") %>% 
  #remove data before deployment
  filter(timestamp >= deployment_dt_utc) %>% 
  #hourly subset to make the next step go faster!
  group_by(individual_local_identifier, yday(timestamp), hour(timestamp)) %>% #subset to hourly
  slice(1) %>% 
  ungroup() %>% 
  group_by(individual_local_identifier) %>% 
  rowwise() %>% 
  mutate(life_stage = case_when(
    between(unique_date, migration_start, migration_end) ~ "migration",
    between(unique_date, first_exploration, migration_start) ~ "post-fledging",
    unique_date < first_exploration ~ "pre-fledging",
    unique_date > migration_end ~ "wintering",
    TRUE ~ NA_character_
  )) %>% 
  ungroup() %>% 
  arrange(individual_local_identifier, timestamp) %>% 
  as.data.frame()

saveRDS(cleaned_gps, "cleaned_gps_for_laterality_map.rds")


#plot
cleaned_gps <- readRDS("cleaned_gps_for_laterality_map.rds")

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

#open the continent boundaries layer
world <- st_read("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
  st_crop(xmin = -17.5, xmax = 43, ymin = -35.6, ymax = 70) %>%
  st_union()

#create a rectangle to be the oceans
Polycoords <- data.frame(long = c(-17.5,43),
                         lat = c(-35.6,70))

pol <- st_polygon(
  list(
    cbind(
      Polycoords$lon[c(1,2,2,1,1)], 
      Polycoords$lat[c(1,1,2,2,1)])
  )
) %>% 
  st_sfc(crs = wgs)

lat_zones_for_map <- lat_zones[-20]

#make the flyway map
x11(height = 5.5, width = 3)
(flyway_map <- ggplot() +
    geom_sf(data = pol, fill = "#F0F8FF", col = "black") +
    geom_sf(data = world, fill = "white", col = "black", linewidth = 0.1) +
    # Use geom_segment instead of geom_hline to limit the horizontal lines
    geom_segment(aes(x = -17, xend = 42.5, y = lat_zones_for_map, yend = lat_zones_for_map),
                 linetype = "dashed", color = "gray75", linewidth = 0.5) +
    # Annotate the lat_zones
    annotate("text", x = rep(-16.9, length(lat_zones_for_map)), y = lat_zones_for_map + 1.5, #adjust the spacing based on the final plot size
             label = paste0(lat_zones_for_map, "°"), size = 2.3, hjust = 0, color = "gray75", fontface = "italic") +
    geom_path(data = subset(cleaned_gps, life_stage == "migration"), 
              aes(x = location_long, y = location_lat, 
                  group = individual_local_identifier), 
              linewidth = .4, lineend = "round", linetype = "solid") +
    geom_path(data = subset(cleaned_gps, life_stage == "post-fledging"), 
              aes(x = location_long, y = location_lat, 
                  group = individual_local_identifier), 
              linewidth = .4, lineend = "round", linetype = "dotted") +
    geom_path(data = subset(cleaned_gps, life_stage == "wintering"), 
              aes(x = location_long, y = location_lat, 
                  group = individual_local_identifier), 
              linewidth = .4, lineend = "round", linetype = "dotted") +
    ggtitle("  Hourly locations") +
    xlim(-17, 42.5) +
    ylim(-35, 65) +
    theme_void() +
    theme(text = element_text(size = 9))
)

ggsave(plot = flyway_map, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/flyway_map.png", 
       device = "png", width = 3, height = 5.5)


#--------------------------------------------------------------------
# STEP 2 whole trajectory: laterality ----------------------------------
#--------------------------------------------------------------------

#open laterality data (prepperd in L04_full_workflow.r). Same as the one used in the daily distribution plots
#add latitudinal zone category
filtered_w_LI <- readRDS("thinned_laterality_w_gps_wind_all_filters.rds") %>% 
  mutate(lat_zone = cut(location_lat_closest_gps_raw,
                        breaks = lat_zones,  # Use lat_zones directly
                        labels = lat_zones[-length(lat_zones)],  # Labels are the start of each zone
                        right = FALSE)) %>%  # Include the left bin edge, exclude the right
  #remove data with no associated latitude
  drop_na(lat_zone)

#violin pilot
x11(height = 5.5, width = 2.3)
(bank <- ggplot(filtered_w_LI, aes(x = mean_roll_mean, y = lat_zone)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.1, size = 0.2, color = "black") +
    geom_violin(trim = T, alpha = 0.7, color = "black", linewidth = 0.2) +
    geom_boxplot(width = 0.1, outliers = F, color = "black", linewidth = 0.2, fill = "white")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
    #annotate("text", x = -90, y = 20.1, label = "Bank angle (°)", size = 3.5, hjust = 0, color = "black") +
    ggtitle("Bank angle (°)") +
    xlim(-90, 90) +
    scale_y_discrete(expand = expansion(add = c(1.7, 1.5))) +  # Add space above and below. to match the lat zones of the map... place the violins in between the two lat lines
    theme_void() +
    theme(plot.margin = unit(c(.1, .5, .4, .5), "lines"),
          text = element_text(size = 9),
          axis.line.x = element_line(), 
          axis.ticks.x = element_line(), 
          axis.text.x = element_text(margin = margin(t = 1)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # Add this line
    )
)

ggsave(plot = bank, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/laterality_panel.png", 
       device = "png", width = 2.3, height = 5.5)



#--------------------------------------------------------------------
# STEP 3 whole trajectory: hourly wind speed ----------------------------------
#--------------------------------------------------------------------

#data prep
life_cycle <- readRDS("updated_life_cycle_nov24.rds")

#open wind speed annotated data from L04a_env_annotation.r
#filter out points before deployment
wind_speed <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24_wind.rds") %>% 
  filter(!(location_lat == 0 & location_long == 0)) %>% 
  mutate(unique_date = as.Date(timestamp)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, deployment_dt_utc, first_exploration, migration_start, migration_end), by = "individual_local_identifier") %>% 
  #remove data before deployment
  filter(timestamp >= deployment_dt_utc) %>% 
  group_by(individual_local_identifier, unique_hr) %>%  #unique hour has the date in it too 
  slice(1) %>%  #make sure there is one value for each hour for each location for each individual
  ungroup() %>% 
  mutate(lat_zone = cut(location_lat,
                        breaks = lat_zones,  # Use lat_zones directly
                        labels = lat_zones[-length(lat_zones)],  # Labels are the start of each zone
                        right = FALSE)) %>%  # Include the left bin edge, exclude the right
  #remove data with no associated latitude
  drop_na(lat_zone)

saveRDS(wind_speed, file = "wind_speed_for_the_map.rds")

#plot

wind_speed <- readRDS("wind_speed_for_the_map.rds")

#violin pilot
x11(height = 5.5, width = 2.3)
(wsp <- ggplot(wind_speed, aes(x = wind_speed, y = lat_zone)) +
    #geom_jitter(height = 0.2, width = 0, alpha = 0.1, size = 0.2, color = "black") + 
    geom_violin(trim = T, alpha = 0.7, color = "black", linewidth = 0.2) + 
    geom_boxplot(width = 0.1, outliers = F, color = "black", linewidth = 0.2, fill = "white")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
    #annotate("text", x = 0, y = 20.1, label = expression("Hourly wind speed (m s"^-1*")"), 
    #         size = 3.5, hjust = 0, color = "black") +
    scale_y_discrete(expand = expansion(add = c(1.7, 1.5))) +  # Add space above and below. to match the lat zones of the map
    #ggtitle(expression("Hourly wind speed (m s"^-1*")")) + #the superscript will make it difficult to line up the panels
    ggtitle("Hourly wind speed (m/s)") +
    theme_void() +
    theme(plot.margin = unit(c(.1, .5, .4, .5), "lines"),
          text = element_text(size = 9),
          axis.line.x = element_line(), 
          axis.ticks.x = element_line(), 
          axis.text.x = element_text(margin = margin(t = 1)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # Add this line
    )
)


ggsave(plot = wsp, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/wind_panel.png", 
       device = "png", width = 2.3, height = 5.5)


#--------------------------------------------------------------------
# STEP 4 migration: Daily distance ----------------------------------
#--------------------------------------------------------------------

#open migration data. prepared in L04_full_workflow......
migr_daily <- readRDS("data_migration_performance_models_2min_daily.rds") %>% 
  mutate(lat_zone = cut(location_lat,
                        breaks = lat_zones,  # Use lat_zones directly
                        labels = lat_zones[-length(lat_zones)],  # Labels are the start of each zone
                        right = FALSE)) %>%  # Include the left bin edge, exclude the right
  #remove data with no associated latitude
  drop_na(lat_zone)
  

#plot


#violin pilot
x11(height = 5.5, width = 2.3)
(distance <- ggplot(migr_daily, aes(x = daily_distance, y = lat_zone)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.1, size = 0.2, color = "black") +
    geom_violin(trim = T, alpha = 0.7, color = "black", linewidth = 0.2) +
    geom_boxplot(width = 0.1, outliers = F, color = "black", linewidth = 0.2, fill = "white")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
    #annotate("text", x = 0, y = 15, label = "Daily distance (km)", size = 3.5, hjust = 0, color = "black") +
    ggtitle("Daily distance (km)") +
    #xlim(-90, 90) +
    scale_y_discrete(expand = expansion(add = c(8.5, 1.5))) +  # Add space above and below. to match the lat zones of the map... place the violins in between the two lat lines
    theme_void() +
    theme(plot.margin = unit(c(.1, .5, .4, .5), "lines"),
          text = element_text(size = 9),
          axis.line.x = element_line(), 
          axis.ticks.x = element_line(), 
          axis.text.x = element_text(margin = margin(t = 1)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # Add this line
    )
)

ggsave(plot = distance, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/distance_panel.png", 
       device = "png", width = 2.3, height = 5.5)



#--------------------------------------------------------------------
# STEP 5 migration: vedba,  ----------------------------------
#--------------------------------------------------------------------

#data prep for all hourly variables

#open migration data. prepared in L04_full_workflow......
migr_hrly <- readRDS("data_migration_performance_models_2min_hrly.rds") %>% 
  mutate(lat_zone = cut(location_lat,
                        breaks = lat_zones,  # Use lat_zones directly
                        labels = lat_zones[-length(lat_zones)],  # Labels are the start of each zone
                        right = FALSE)) %>%  # Include the left bin edge, exclude the right
  #remove data with no associated latitude
  drop_na(lat_zone)


#plot

#vedba
x11(height = 5.5, width = 2.3)
(vedba <- ggplot(migr_hrly, aes(x = hrly_mean_vedba, y = lat_zone)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.1, size = 0.2, color = "black") +
    geom_violin(trim = T, alpha = 0.7, color = "black", linewidth = 0.2) +
    geom_boxplot(width = 0.1, outliers = F, color = "black", linewidth = 0.2, fill = "white")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
    #annotate("text", x = 0, y = 15, label = "Daily distance (km)", size = 3.5, hjust = 0, color = "black") +
    ggtitle("Hourly VeDBA (g)") +
    #xlim(-90, 90) +
    scale_y_discrete(expand = expansion(add = c(8.5, 1.5))) +  # Add space above and below. to match the lat zones of the map... place the violins in between the two lat lines
    theme_void() +
    theme(plot.margin = unit(c(.1, .5, .4, .5), "lines"),
          text = element_text(size = 9),
          axis.line.x = element_line(), 
          axis.ticks.x = element_line(), 
          axis.text.x = element_text(margin = margin(t = 1)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # Add this line
    )
)

ggsave(plot = vedba, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/vedba_panel.png", 
       device = "png", width = 2.3, height = 5.5)


#cumulative yaw
x11(height = 5.5, width = 2.3)
(yaw <- ggplot(migr_hrly, aes(x = hrly_mean_cum_yaw, y = lat_zone)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.1, size = 0.2, color = "black") +
    geom_violin(trim = T, alpha = 0.7, color = "black", linewidth = 0.2) +
    geom_boxplot(width = 0.1, outliers = F, color = "black", linewidth = 0.2, fill = "white")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
    #annotate("text", x = 0, y = 15, label = "Daily distance (km)", size = 3.5, hjust = 0, color = "black") +
    ggtitle("Hourly cumulative yaw (°)") +
    #xlim(-90, 90) +
    scale_y_discrete(expand = expansion(add = c(8.5, 1.5))) +  # Add space above and below. to match the lat zones of the map... place the violins in between the two lat lines
    theme_void() +
    theme(plot.margin = unit(c(.1, .5, .4, .5), "lines"),
          text = element_text(size = 9),
          axis.line.x = element_line(), 
          axis.ticks.x = element_line(), 
          axis.text.x = element_text(margin = margin(t = 1)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # Add this line
    )
)

ggsave(plot = yaw, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/yaw_panel.png", 
       device = "png", width = 2.3, height = 5.5)


#pitch
x11(height = 5.5, width = 2.3)
(pitch <- ggplot(migr_hrly, aes(x = hrly_mean_mean_pitch, y = lat_zone)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.1, size = 0.2, color = "black") +
    geom_violin(trim = T, alpha = 0.7, color = "black", linewidth = 0.2) +
    geom_boxplot(width = 0.1, outliers = F, color = "black", linewidth = 0.2, fill = "white")+
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
    #annotate("text", x = 0, y = 15, label = "Daily distance (km)", size = 3.5, hjust = 0, color = "black") +
    ggtitle("Hourly average pitch (°)") +
    #xlim(-90, 90) +
    scale_y_discrete(expand = expansion(add = c(8.5, 1.5))) +  # Add space above and below. to match the lat zones of the map... place the violins in between the two lat lines
    theme_void() +
    theme(plot.margin = unit(c(.1, .5, .4, .5), "lines"),
          text = element_text(size = 9),
          axis.line.x = element_line(), 
          axis.ticks.x = element_line(), 
          axis.text.x = element_text(margin = margin(t = 1)),
          panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5) # Add this line
    )
)

ggsave(plot = pitch, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/pitch_panel.png", 
       device = "png", width = 2.3, height = 5.5)


