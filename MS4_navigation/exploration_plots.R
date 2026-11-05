# script for preparing plots of multi-year migration of honey buzzards
# Elham Nourani. 13.07.2026
#elham.nourani@unil.ch, Lausanne, CH


library(move2)
library(move2utils)
library(tidyverse)
library(lubridate)
library(sf)

#---------------------------------------------------------------------------------
## Step 1: download GPS data from movebank and clean up                      #####
#---------------------------------------------------------------------------------

gps <- movebank_download_study("European Honey Buzzard_Finland",
                                sensor_type_id = "gps",
                                remove_movebank_outliers = TRUE)

saveRDS(gps, file = "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/data/raw_gps_aug2026.rds")


#clean using Kami's ftn
cleaned_gps <- gps |> 
  mt_filter_gps_quality() |> 
  mt_filter_movebank_visible() |> 
  #filter(!st_is_empty(geometry)) |> 
  mt_filter_unique() |> 
  mt_filter_unique(criterion = "first") #|> 
  #arrange(mt_time()) |> 

cleaned_gps2 <- cleaned_gps |> 
  mt_clean_track(mass = 1.1, mode = "flying") #max mass in kg


saveRDS(cleaned_gps2, file = "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/data/cleaned_gps_aug2026.rds")


hrly_gps <- cleaned_gps |> 
  mt_thin_time(interval = "1 hour", tolerance = "5 min")  

hrly_gps$year <- year(mt_time(hrly_gps))

saveRDS(hrly_gps, file = "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/data/hrly_cleaned_gps_aug2026.rds")

#---------------------------------------------------------------------------------
## Step 2: summary info                                                       #####
#---------------------------------------------------------------------------------

africa <- rnaturalearth::ne_countries(continent = "Africa", returnclass = "sf", scale = 50) %>%
  st_union() %>%
  st_make_valid()

in_africa <- st_intersects(hrly_gps, africa, sparse = FALSE)[, 1]

year_status <- data.frame(individual_local_identifier = as.character(hrly_gps$individual_local_identifier),
                           year = hrly_gps$year,
                           in_africa = in_africa) %>%
  group_by(individual_local_identifier, year) %>%
  summarize(any_africa = any(in_africa), any_non_africa = any(!in_africa), .groups = "drop") %>%
  mutate(status = case_when(
    any_africa & any_non_africa ~ "migrated",
    any_africa & !any_non_africa ~ "stayed_in_africa",
    !any_africa & any_non_africa ~ "stayed_in_europe"
  ))

individual_summary <- year_status %>%
  group_by(individual_local_identifier) %>%
  summarize(n_years = n(),
            years_migrated = paste(year[status == "migrated"], collapse = ", "),
            years_stayed_in_africa = paste(year[status == "stayed_in_africa"], collapse = ", "),
            years_stayed_in_europe = paste(year[status == "stayed_in_europe"], collapse = ", "),
            .groups = "drop") %>%
  arrange(desc(n_years))

individual_summary

dir.create("figures", showWarnings = FALSE)
write.csv(individual_summary, "figures/individual_summary.csv", row.names = FALSE)

multi_year_ids <- year_status %>%
  filter(status == "migrated") %>%
  count(individual_local_identifier, name = "n_years_migrated") %>%
  filter(n_years_migrated > 1) %>%
  pull(individual_local_identifier)

#-----------------------------------------------------------------------------------
## Step 3: map of all trajectories, colored by individual                     ######
#-----------------------------------------------------------------------------------

all_tracks <- mt_track_lines(gps_clean)

(all_individuals_map <- ggplot() +
    geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)) +
    geom_sf(data = all_tracks, aes(color = individual_local_identifier), linewidth = 0.4) +
    coord_sf(xlim = c(-17, 53), ylim = c(-35.6, 70)) +
    labs(color = "Individual") +
    theme_linedraw() +
    theme(legend.position = "none"))

ggsave("figures/all_individuals_map.png", all_individuals_map, width = 9, height = 7, dpi = 300)

#---------------------------------------------------------------------------------
## Step 4: map of multi-year individuals, one panel per individual,          #####
## colored by year (consistent colors across panels)                        #####
#---------------------------------------------------------------------------------

multi_year_data <- hrly_gps %>%
  #filter(as.character(individual_local_identifier) %in% multi_year_ids)
  filter(as.character(individual_local_identifier) %in% c("D225_226","D225_228", "D225_231", "D320_475", "D321_584", "D323_154", "D324_510", "D324_512", "D329_014"))

multi_year_data <- mt_set_track_id(multi_year_data,
                                    paste(multi_year_data$individual_local_identifier, multi_year_data$year, sep = "::"))

multi_year_tracks <- mt_track_lines(multi_year_data) %>%
  separate(individual_local_identifier, into = c("individual_local_identifier", "year"), sep = "::") %>%
  mutate(year = factor(year, levels = sort(unique(year))))

multi_year_map <- ggplot() +
    geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)) +
    geom_sf(data = multi_year_tracks, aes(color = year), linewidth = 0.5) +
    coord_sf(xlim = c(-17, 53), ylim = c(-35.6, 70)) +
    facet_wrap(~ individual_local_identifier) +
    labs(color = "Year") +
    theme_linedraw()

ggsave("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/exploration_figs/maps/multi_year_map.png", multi_year_map, width = 10, height = 8, dpi = 300)


#---------------------------------------------------------------------------------
## Step 4b: map of multi-year for all individuals,                          #####
## colored by year (consistent colors across panels)                        #####
#---------------------------------------------------------------------------------


all_year_data <- mt_set_track_id(hrly_gps,
                                   paste(hrly_gps$individual_local_identifier, hrly_gps$year, sep = "::"))

all_ind_tracks <- mt_track_lines(all_year_data) %>%
  separate(individual_local_identifier, into = c("individual_local_identifier", "year"), sep = "::") %>%
  mutate(year = factor(year, levels = sort(unique(year))))

multi_year_map_all <- ggplot() +
    geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)) +
    geom_sf(data = all_ind_tracks, aes(color = year), linewidth = 0.5) +
    coord_sf(xlim = c(-17, 53), ylim = c(-35.6, 70)) +
    facet_wrap(~ individual_local_identifier) +
    labs(color = "Year") +
    theme_linedraw()


ggsave("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/exploration_figs/maps/multi_year_map_all.png", multi_year_map_all, width = 15, height = 8, dpi = 300)


#-----------------------------------------------------------------------------------
## Step 3b: map of all trajectories, colored by individual. One panel per year ######
#-----------------------------------------------------------------------------------

map_years <- ggplot() +
   geom_sf(data = rnaturalearth::ne_coastline(returnclass = "sf", scale = 50)) +
   geom_sf(data = all_ind_tracks, aes(color = individual_local_identifier), linewidth = 0.5) +
   coord_sf(xlim = c(-17, 53), ylim = c(-35.6, 70)) +
   facet_wrap(~ year) +
   labs(color = "Individual") +
   theme_linedraw()

ggsave("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MSxx_navigation/exploration_figs/maps/per_year_maps.png", map_years, width = 9, height = 7, dpi = 300)
