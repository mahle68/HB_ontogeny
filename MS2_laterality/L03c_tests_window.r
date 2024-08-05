#This code runs tests for whether honey buzzards have a preference for left or right bank angles over timw (using a moving window)
#at the scale of each 8-second burst
#Elham Nourani PhD.
#July 24. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)
library(viridis)
library(lme4)
library(mgcv)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#-------------------------------------------------------------------------------------
# STEP1: annotate data with life cycle stage
#-------------------------------------------------------------------------------------
#open data with laterality index calculated from L03b_tests_per_burst.r

#this file contains laterality index calculated per burst. NOT filtered for circling flight
laterality_1sec_LS <- readRDS("laterality_index_per_8sec_burst_LS.rds")

#identify circling cases with varying radii, first try to 45 degree 

#filter for circling flight only
laterality_circling <- laterality_1sec_LS %>% 
  filter(n_records >= 8) %>%  #remove bursts that are shorter than 8 seconds
  filter(cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45) %>%  #only keep thermalling flight. use a threshold of 45 degrees in 8 seconds.
  as.data.frame()

#assign a new column with every consecutive 24 instances of circling
laterality_24 <- laterality_circling %>% 
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = T) %>% 
  mutate(row_of_24 = ceiling(row_number() / 24)) %>% #set of 24 consecutive 8-sec bursts
  ungroup() %>% 
  mutate(row_of_24_id = paste0(individual_local_identifier, "_", row_of_24)) %>% 
  #recalculate laterality, but for each 24 consecutive bouts
  group_by(row_of_24_id) %>% #group by individual ID and day
  mutate(bank_direction = ifelse(cumulative_roll_8sec < 0, "left",
                                 ifelse(cumulative_roll_8sec > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_8sec < 0, "left",
                                    ifelse(cumulative_yaw_8sec > 0, "right", "straight"))) %>% 
  summarise(n_bursts = n(), #this should be 24 for all bursts
            median_timestamp = median(start_timestamp),
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            heading_left = sum(heading_direction == "left"),
            heading_right = sum(heading_direction == "right"),
            heading_straight = sum(heading_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
            laterality_heading = (heading_right - heading_left)/(heading_right + heading_left),
            row_of_24_timestamp = median(start_timestamp),
            row_of_24_tagging = median(day_since_tagging),
            individual_local_identifier = head(individual_local_identifier, 1),
            .groups = 'drop') %>% 
  as.data.frame()



#-------------------------------------------------------------------------------------
# STEP2: some plotting
#-------------------------------------------------------------------------------------

#### ridgelines for laterality in banking angle
ggplot(laterality_24, aes(x = laterality_bank, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()

ggplot(laterality_24, aes(x = laterality_bank, y = individual_local_identifier, fill = life_stage)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("post-fledging" = "yellow", "migration" = "#c4c4fc", "wintering" = "red")) +
  theme_minimal() +
  labs(fill = "Life Stage")


#### ridgelines for laterality in heading
ggplot(laterality_circling, aes(x = laterality_heading, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()


#### laterality over time
ggplot(data = laterality_24, aes(x = row_of_24_tagging, y = laterality_bank)) +
  geom_point() +
  geom_smooth( method = "loess") +
  facet_wrap(vars(individual_local_identifier))

#### gam
one_ind <- 
gam_model <- gam(laterality_bank ~ s(row_of_24_tagging, by = individual_local_identifier, bs = "fs", k = 10) + individual_local_identifier, 
                 data = laterality_24, method = "REML")
