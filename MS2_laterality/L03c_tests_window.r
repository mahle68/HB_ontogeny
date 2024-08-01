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

#-------------------------------------------------------------------------------------
# STEP1: annotate data with life cycle stage
#-------------------------------------------------------------------------------------
#open data with laterality index calculated from L03b_tests_per_burst.r

#this file contains laterality index calculated per burst. NOT filtered for circling flight
laterality_1sec_LS <- readRDS("laterality_index_per_8sec_burst_LS.rds")

#identify circling cases with varying radii, first try to 45 degree 


#filter for circling flight only
laterality_circling <- laterality_1sec_LS %>% 
  #mutate(burst_dur2 = difftime(end_timestamp, start_timestamp)) %>% ## OR use the n of records for this. that would be the number of rows
  #filter(burst_dur2 > 6.5 &  #remove bursts that are shorter than 6.5 seconds
  filter(n_records > 6.5 &  #remove bursts that are shorter than 6.5 seconds
           cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45) %>%  #only keep thermalling flight. use a threshold of 45 degrees in 8 seconds.
  as.data.frame()

#assign a new column with every consecutive 24 instances of circling
laterality_24 <- laterality_circling %>% 
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = T) %>% 
  mutate(rows_of_24 = paste0(individual_local_identifier, "_", imu_burst_id, "_", ceiling(row_number() / 24))) %>% 
  #recalculate laterality, but for each 24 consecutive bouts
  group_by(rows_of_24) %>% #group by individual ID and day
  mutate(bank_direction = ifelse(cumulative_roll_8sec < 0, "left",
                                 ifelse(cumulative_roll_8sec > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_8sec < 0, "left",
                                    ifelse(cumulative_yaw_8sec > 0, "right", "straight"))) %>% 
  summarise(n_bursts = n(),
            median_timestamp = median(start_timestamp),
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            heading_left = sum(heading_direction == "left"),
            heading_right = sum(heading_direction == "right"),
            heading_straight = sum(heading_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
            laterality_heading = (heading_right - heading_left)/(heading_right + heading_left),
            .groups = 'drop') %>% 
  as.data.frame()



#-------------------------------------------------------------------------------------
# STEP2: some plotting
#-------------------------------------------------------------------------------------

#### ridgelines for laterality in banking angle
ggplot(laterality_circling, aes(x = laterality_bank, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()

ggplot(laterality_circling, aes(x = laterality_bank, y = individual_local_identifier, fill = life_stage)) +
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


