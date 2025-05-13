# script for plotting the laterality in the Finnish population of honey buzzards.
# Elham Nouani, PhD. 
# 21.11.2024, Konstanz, DE

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)


setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}


#######figure 2!!

##THE MAP -----------------------------------------------------------
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

lat_zones_for_map <-  seq(-30,65, by = 15)

stage_colors <- c(
  "post-fledging" = "#31688e",  # Dark Green
  "migration" = "#b5de2b",      # Dark Red
  "wintering" = "#1f9e89"       # Dark Cyan
)

#make the flyway map
X11(height = 4.09, width = 2.3)
(flyway_map <- ggplot() +
    geom_sf(data = pol, fill = "gray20", col = "gray20") +
    geom_sf(data = world, fill = "white", col = "gray20", linewidth = 0.1) +
    # Use geom_segment instead of geom_hline to limit the horizontal lines
    geom_segment(aes(x = -17, xend = 42.5, y = lat_zones_for_map, yend = lat_zones_for_map),
                 linetype = "dashed", color = "gray75", linewidth = 0.5) +
    # Annotate the lat_zones
    annotate("text", x = rep(-16.9, length(lat_zones_for_map)), y = lat_zones_for_map + 1.5, #adjust the spacing based on the final plot size
             label = paste0(lat_zones_for_map, "°"), size = 2, hjust = 0, color = "gray75", fontface = "italic") +
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
    theme(text = element_text(size = 8),
          plot.title = element_text(face = "bold", margin = margin(b = 5)), # Add margin below the title
          plot.margin = margin(0, 0, 0, 0, "pt"))
)

ggsave(plot = flyway_map, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/map_multi_panel.pdf", 
       device = "pdf", width = 2.3, height = 4.09, dpi = 600)

##DISTRIBUTION OF LATERALITY INDEX -----------------------------------------------------------

#read in filtered data. this is not filtered for days since tagging
filtered_w_LI <- readRDS("thinned_laterality_w_gps_wind_all_filters2.rds") %>% 
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
manual_orders <- c("D329_015", "D326_193", "D163_696", "D225_232", "D329_012",  
                   "D225_236", "D329_014", "D225_231","D299_269", "D320_474", "D321_345", "D225_226", "D225_234", "D321_584", 
                   "D323_154","D311_750", "D321_348", "D321_349", "D323_155", "D324_510","D225_227", "D225_228", 
                   "D299_270", "D299_271", "D320_475","D321_583" ,"D324_511", "D324_512", "D324_513", "D326_192", "D329_013")




# Reorder the factor levels of individual_local_identifier
filtered_w_LI$individual_local_identifier <- factor(filtered_w_LI$individual_local_identifier, 
                                                    levels = manual_orders)

#gray out the id label for individuals that don't have complete trajectories
text_col <- ifelse(filtered_w_LI$individual_local_identifier %in% c("D299_270", "D299_269", "D320_474", "D329_015", "D326_193", "D225_236") , "gray30", "black")

#### ----------------------- plots for life stage and age
#### ridgelines for daily laterality for different life stages
#keep one row per day for each individual. to avoid overplotting
day_LI <- filtered_w_LI %>% 
  group_by(individual_local_identifier, unique_date) %>% 
  slice(1)

#density distributions
X11(height = 4.09, width = 4.3)
(p <- ggplot(day_LI, aes(x = laterality_bank_day, y = individual_local_identifier, color = laterality_dir_stage, fill = laterality_dir_stage)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = c("right_handed" =  "#0d0887", "ambidextrous" = "#fb9f3a", "left_handed" = "#9c179e"),
                      labels = c("Right", 
                                 "No bias",
                                 "Left")) +
    scale_color_manual(values = c("right_handed" =  "#0d0887", "ambidextrous" = "#fb9f3a", "left_handed" = "#9c179e"),
                       labels = c("Right", 
                                  "No bias",
                                  "Left")) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    facet_wrap2(
      vars(life_stage), nrow = 1, scales = "free_x",
      strip = strip_themed(
        background_x = list(element_rect(color = "#31688e"),
                            element_rect(color = "#b5de2b"),
                            element_rect(color = "#1f9e89")),
        # text_x = list(
        #   "post-fledging" = element_text(color = "#31688e"),
        #   "migration" = element_text(color = "#b5de2b"),
        #   "wintering" = element_text(color = "#1f9e89")
      ),
      labeller = labeller(life_stage = c(
        "post-fledging" = "Post-fledging",
        "migration" = "Migration",
        "wintering" = "Wintering"
      ))) +
    labs(x = "Laterality index",
         y = "Individual ID",
         fill = "Lateralization",
         color = "Lateralization") +
    ggtitle("b") +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(face = "bold"), # make title bold
          panel.grid.major.y = element_line(color = "gray80"), # Add horizontal grid lines
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2)),
          axis.text.y = element_text(size = 6, face = "italic"), # Make y-axis text smaller and bold
          legend.key.width=unit(.3,"cm"),
          legend.key.height=unit(.15,"cm"),
          legend.position = "right")
)

ggsave(plot = p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/LI_distr_multi_panel.pdf", 
       device = "pdf", width = 4.3, height = 4.09, dpi = 600)
