
library(mapview)
library(sf)
library(terra)
library(lubridate)
library(tidyverse)
#library(rnaturalearth)
#library(rnaturalearthdata)

wgs <- crs("+proj=longlat +datum=WGS84 +no_defs")


# open a sample csv
sample <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/Pritish_collab_IMU/Sep29_23/flightMetrics_D320.csv") %>% 
  drop_na(location_long_closest_gps) %>% 
  st_as_sf(coords = c("location_long_closest_gps", "location_lat_closest_gps"))

st_crs(sample) <- "EPSG:4326"

mapview(sample %>% 
          mutate(abs_heading = abs(netHeadChange)), zcol = "abs_heading")


#try with ggplot

load("/home/enourani/Desktop/flyway_layer.RData")
#base <- ne_coastline(scale = 'medium', returnclass = 'sf') #flyway

ws <- ggplot(data = flyway) +
  geom_sf(col = "gray", fill = "gray") +
  #coord_sf(xlim = c(-10, 38), ylim = c(3, 64), expand = FALSE) +
  geom_path(data = sample, aes(x = location_long_closest_gps, y = location_lat_closest_gps, col = propFlap), 
            size = 2.5, lineend = "round") +
  scale_colour_viridis(option = "magma", na.value = "white", name = "m/s", alpha = 0.7) +
  theme_linedraw() +
  scale_x_continuous(breaks = c(0,30)) +
  scale_y_continuous(breaks = c(10,30,50)) +
  theme(axis.text = element_text(size = 10, colour = 1),
        legend.text = element_text(size = 10, colour = 1), 
        legend.title = element_text(size = 10, colour = 1),
        legend.position = "right",
        legend.background = element_rect(colour = NULL, fill = "white"))+
  labs(x = NULL, y = NULL, title = "Wind support") +
  facet_wrap(.~year, nrow = 1)
