# script for the plots in Safi et al 2025
# Elham Nouani, PhD. 
# 21.11.2024, Konstanz, DE

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(grid)
library(gridExtra)
library(mapview)


setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")
lat_zones <- seq(-30,65, by = 5)

###---------------------------------------------------------------###
### Fig 2: migration map and distribution of daily laterality     ###
###---------------------------------------------------------------###

##THE MAP -----------------------------------------------------------
cleaned_gps <- readRDS("cleaned_gps_for_laterality_map.rds")

wgs <- st_crs("+proj=longlat +datum=WGS84 +no_defs")

#open the continent boundaries layer
world <- st_read("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/GIS_files/continent_shapefile/World_Continents.shp") %>% 
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

stage_colors <- c(
  "post-fledging" = "#006400",  # Dark Green
  "migration" = "#8B0000",      # Dark Red
  "wintering" = "#008B8B"       # Dark Cyan
)

#make the flyway map
X11(height = 5.5, width = 3)
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
              color = stage_colors["migration"], linewidth = .4, lineend = "round", linetype = "solid") +
    geom_path(data = subset(cleaned_gps, life_stage == "post-fledging"), 
              aes(x = location_long, y = location_lat, 
                  group = individual_local_identifier), 
              color = stage_colors["post-fledging"], linewidth = .4, lineend = "round", linetype = "solid") +
    geom_path(data = subset(cleaned_gps, life_stage == "wintering"), 
              aes(x = location_long, y = location_lat, 
                  group = individual_local_identifier), 
              color = stage_colors["wintering"], linewidth = .4, lineend = "round", linetype = "solid") +
    ggtitle("  a") +
    xlim(-17, 42.5) +
    ylim(-35, 65) +
    theme_void() +
    theme(text = element_text(size = 9))
)


##Distribution of laterality -----------------------------------------------------------
#read in filtered data. this is not filtered for days since tagging
filtered_w_LI <- readRDS("thinned_laterality_w_gps_wind_all_filters.rds") %>% 
  mutate(life_stage = factor(life_stage, levels = c("post-fledging", "migration", "wintering"))) #reorder life stage

#reorder based on the value of handedness during post-fledging
# Determine the order of individual_local_identifier based on laterality_dir
ordered_identifiers <- filtered_w_LI %>%
  filter(life_stage == "post-fledging") %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  select(individual_local_identifier, laterality_dir_stage) %>% 
  arrange(laterality_dir_stage) %>% 
  ungroup() %>% 
  pull(individual_local_identifier)

#manually edit the orders so that the inds that dont change level are first
manual_orders <- c("D225_236", "D329_014","D163_696", "D225_232",  "D326_193", "D329_012",  "D329_015", 
                   "D225_231","D299_269", "D320_474", "D321_345", "D225_226", "D225_234", "D321_584", 
                   "D323_154", "D225_227", "D311_750", "D321_348", "D321_349", "D323_155", "D324_510","D225_228", 
                   "D299_270", "D299_271", "D320_475","D321_583" ,"D324_511", "D324_512", "D324_513", "D326_192", "D329_013")

# Reorder the factor levels of individual_local_identifier
filtered_w_LI$individual_local_identifier <- factor(filtered_w_LI$individual_local_identifier, 
                                                    levels = manual_orders)


#### ----------------------- plots for life stage and age
#### ridgelines for daily laterality for different life stages
#keep one row per day for each individual. to avoid overplotting
day_LI <- filtered_w_LI %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  slice(1)

#density distributions

# Assuming day_LI is your data frame
unique_life_stages <- unique(day_LI$life_stage)

# Define a mapping for the new titles
stage_labels <- c(
  "post-fledging" = "Post-fledging",
  "migration" = "Migration",
  "wintering" = "Wintering"
)

stage_colors <- c(
  "post-fledging" = "#006400",  # Dark Green
  "migration" = "#8B0000",      # Dark Red
  "wintering" = "#008B8B"       # Dark Cyan
)

# Create a list to store the plots
plots <- list()

# Loop over each life stage and create a plot
for (i in seq_along(unique_life_stages)) {
  stage <- unique_life_stages[i]
  
  # Filter the data for the current life stage
  stage_data <- subset(day_LI, life_stage == stage)
  
  # Create the plot for the current life stage
  p <- ggplot(stage_data, aes(x = laterality_bank_day, y = individual_local_identifier, 
                              color = laterality_dir_stage, fill = laterality_dir_stage)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, linewidth = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = c("right_handed" =  "#0d0887", "ambidextrous" = "#fb9f3a", "left_handed" = "#9c179e"),
                      labels = c("Right-handed \nLI = 0.25 to 1.0", 
                                 "Ambidextrous \nLI = -0.25 to 0.25",
                                 "Left-handed \nLI = -1.0 to -0.25")) +
    scale_color_manual(values = c("right_handed" =  "#0d0887", "ambidextrous" = "#fb9f3a", "left_handed" = "#9c179e"),
                       labels = c("Right-handed \nLI = 0.25 to 1.0", 
                                  "Ambidextrous \nLI = -0.25 to 0.25",
                                  "Left-handed \nLI = -1.0 to -0.25")) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    ggtitle(stage_labels[stage]) +
    labs(x = if (i == 2) "Laterality index" else "", # Only show x-axis label for the middle plot
         y = if (i == 1) "Individual ID" else NULL, # Only show y-axis label for the first plot
         fill = "Handedness",
         color = "Handedness") +
    theme_classic() +
    theme(text = element_text(size = 9),
          legend.text = element_text(size = 8.5),
          legend.title = element_text(size = 9),
          plot.title = element_text(hjust = 0.7,color = stage_colors[stage]), #center plot title and Set the title color
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)),
          axis.text.y = if (i == 1) element_text() else element_blank(), # Hide y-axis text for other plots
          legend.position = "none", # Remove legend from individual plots
          legend.key.spacing.y = unit(0.3, "cm"))
  
  # Add the plot to the list
  plots[[stage]] <- p
}

# Arrange the plots in one row
X11(height = 5.5, width = 4.2)
grid.arrange(grobs = plots, nrow = 1, widths = c(0.3, 0.2, 0.2),
             top = textGrob("  b", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 9)))


##also add flyway_map
X11(height = 4.7, width = 7)
p <- grid.arrange(
  flyway_map, 
  arrangeGrob(grobs = plots, nrow = 1, widths = c(0.3, 0.2, 0.2),
              top = textGrob("      b", vjust = 1, x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 10))),
  nrow = 1,
  widths = c(0.4, 0.6) # Adjust the widths to add space between flyway_map and the other plots
  
)

ggsave(plot = p, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/Fig_2.pdf", 
       device = "pdf", width = 7, height = 4.7)
