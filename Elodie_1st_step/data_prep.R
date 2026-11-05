#prep data for 1st step project 2026
#Honey buzzard flight behavior along a gradient of human population density
#elham.nourani@unil.ch
#12.10.2026 Lausanne, CH

library(move2)
library(move2utils)
library(tidyverse)
library(lubridate)
library(sf)
library(mapview)

#---------------------------------------------------------------------------------
## Step 1: download GPS data from movebank and clean up                      #####
#---------------------------------------------------------------------------------

#individual IDs selected by Elodie
IDs <- c("D163_696", "D225_231", "D299_270", "D321_348", "D321_583", "D321_584", "D324_511")

#use gps data downloaded in Aug 2026

gps <- readRDS("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/data/raw_gps_aug2026.rds")

#clean up and sub-sample to 20 minute intervals

#clean wihtout using Kami's ftn
cleaned_gps <- gps |> 
  mt_filter_gps_quality() |> 
  mt_filter_movebank_visible() |> 
  #filter(!st_is_empty(geometry)) |> 
  mt_filter_unique() |> 
  mt_filter_unique(criterion = "first") #|> 
#arrange(mt_time()) |> 

saveRDS(cleaned_gps, "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/data/raw_gps_aug2026.rds")

#apply Kami's ftn... load from file if needed. 
#cleaned_gps <- cleaned_gps |> 
#  mt_clean_track(mass = 1.1, mode = "flying") #max mass in kg


#Elodie individuals and Europe only
Elodie_inds_EUR <- cleaned_gps |> 
  filter_track_data(.track_id = IDs) |> 
  st_crop(xmin = unname(st_bbox(cleaned_gps)[1]),
          ymin = 35.9,
          xmax = unname(st_bbox(cleaned_gps)[3]),
          ymax = unname(st_bbox(cleaned_gps)[4]))


#20 min filter
 gps_20min <- Elodie_inds_EUR |> 
   mt_thin_time(interval = "20 mins", tolerance = "5 min", remove = TRUE)  
 

#convert to sf
#firs move all track attributes to the event table
Elodie_inds_sf <- mt_as_event_attribute(gps_20min, names(mt_track_data(gps_20min))) |> 
  st_as_sf() |> 
  #keep only the relevant columns
  select(individual_local_identifier, tag_local_identifier, eobs_horizontal_accuracy_estimate, eobs_speed_accuracy_estimate, eobs_status, eobs_temperature, gps_dop, gps_satellite_count, ground_speed, heading, height_above_ellipsoid,
         timestamp, geometry, deploy_on_timestamp)

#drop the move2 component
class(Elodie_inds_sf) <- class(Elodie_inds_sf) |> 
  setdiff("move2")


Elodie_inds_sf$individual_local_identifier <- as.character(Elodie_inds_sf$individual_local_identifier)
Elodie_inds_sf$tag_local_identifier <- as.character(Elodie_inds_sf$tag_local_identifier)

#write one RDS file
#saveRDS(Elodie_inds_sf, "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Mentoring_Teaching/BEC_UNIL_MSc/1st_step/2026_HB_human_density/tracking_data_subset/honey_buzzard_subset.rds")
saveRDS(Elodie_inds_sf, "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Mentoring_Teaching/BEC_UNIL_MSc/1st_step/2026_HB_human_density/tracking_data_subset/honey_buzzard_subset_20min.rds")







# #for each individual, visualize to make sure there are no outliers, then write a csv to file
# 
# for(i in IDS){
#   
#   one_ind <- Elodie_inds_sf |> 
#     filter(individual_local_identifier == i)
#   
#   mapview(one_ind)
#   
#   saveRDS(one_ind, paste0("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Mentoring_Teaching/BEC_UNIL_MSc/1st_step/2026_HB_human_density/tracking_data_subset/",i,".rds" ))
# }
# 
# #20 min filter
# gps_20min <- cleaned_gps |> 
#   mt_thin_time(interval = "20 mins", tolerance = "5 min", remove = TRUE)  
# 
# saveRDS(gps_20min, "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/data/clean_gps_20min_aug2026.rds")
# 
# 
# #Elodie individuals and Europe only
# Elodie_inds <- gps_20min |> 
#   filter_track_data(.track_id = IDs) |> 
#   st_crop(xmin = unname(st_bbox(gps_20min)[1]),
#           ymin = 35.9,
#           xmax = unname(st_bbox(gps_20min)[3]),
#           ymax = unname(st_bbox(gps_20min)[4]))
#   
# 
# 
# 
#   
# 
# mapview(Elodie_inds_sf |>  filter(individual_local_identifier == "D163_696"))


##plot the individuals


library(ggplot2)
library(sf)
library(rnaturalearth)

# grab a basemap
world <- ne_countries(scale = "medium", returnclass = "sf")

# bbox of your data for zooming
bb <- st_bbox(Elodie_inds_sf)

ggplot() +
  geom_sf(data = world, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = Elodie_inds_sf, aes(color = individual_local_identifier),
          size = 0.4, alpha = 0.6) +
  coord_sf(xlim = c(bb["xmin"], bb["xmax"]),
           ylim = c(bb["ymin"], bb["ymax"]),
           expand = FALSE) +
  scale_color_brewer(palette = "Dark2", name = "Individual") +
  labs(x = NULL, y = NULL,
       title = "Elodie's individuals — Europe") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "aliceblue", color = NA),
        legend.position = "right")
