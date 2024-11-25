#This code runs permutation tests to test whether honey buzzards have a preference for left or right bank angles during different life cycle stage
#Elham Nourani PhD.
#July 22. 2024. Konstanz, DE. 

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
#open data with laterality index calculated

#this file contains laterality index calculated per day. in the wintering ground and some points pre-migration, there was only one IMU recorded per day (n_bursts_in_day). filter as necessary
laterality_8sec <- readRDS("laterality_index_per_day.rds")


#use Ellen's life cycle stage estimations
life_cycle <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/HB_dates_Fin_all_M2.rds") %>% 
  mutate_at(c("migration_start", "migration_end"), as.Date) %>% 
  rename("individual_local_identifier" = "ind_id")

#identify the last day of data collection for individuals that have NA for end of migration (they died)
no_migr_end <- life_cycle %>% filter(is.na(migration_end)) %>% pull(ind_id)

missing_migration_ends <- laterality_8sec %>% 
  filter(individual_local_identifier %in% no_migr_end) %>% 
  mutate(unique_date = as.Date(unique_date)) %>% 
  group_by(individual_local_identifier) %>% 
  #assing end of data collection as migration end
  summarize(migration_end = max(unique_date)) %>% 
  ungroup() %>% 
  #add the other columns of life_cycle, so I can use rows-update in the next step to fill in the NAs
  left_join(life_cycle %>% select(-migration_end), by = "individual_local_identifier")

#update life_cycle to deal with missing dates
life_cycle <- life_cycle %>% 
  rows_update(missing_migration_ends, by = "individual_local_identifier") %>% 
  #one individual is also missing the start of migration, assign mid september
  mutate(migration_start = ifelse(is.na(migration_start), as.Date("2023-09-15"), migration_start) %>% as.Date())

saveRDS(life_cycle, file = "updated_life_cycle_nov24.rds")


laterality_8sec_LS <- laterality_8sec %>%
  mutate(unique_date = as.Date(unique_date)) %>% 
  full_join(life_cycle %>% select(individual_local_identifier, migration_start, migration_end), by = "individual_local_identifier") %>% 
  group_by(individual_local_identifier) %>% 
  rowwise() %>% 
  mutate(life_stage = case_when(
    between(unique_date, migration_start, migration_end) ~ "migration",
    unique_date < migration_start ~ "post-fledging",
    unique_date > migration_end ~ "wintering",
    TRUE ~ NA_character_
  )) %>% 
  ungroup() %>%
  filter(n_bursts_in_day >= 3) %>% #filter for number of IMU bursts per day.  
  as.data.frame()

saveRDS(laterality_8sec_LS, file = "laterality_daily_w_LS.rds")

#-------------------------------------------------------------------------------------
# STEP2: some plotting
#-------------------------------------------------------------------------------------

#### ridgelines for laterality in banking angle
ggplot(laterality_8sec_LS, aes(x = laterality_bank, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()

#### ridgelines for laterality in heading
ggplot(laterality_8sec_LS, aes(x = laterality_heading, y = individual_local_identifier)) +
  stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                      jittered_points = TRUE, 
                      point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  facet_wrap(vars(factor(life_stage, levels = c("post-fledging", "migration", "wintering")))) +
  theme_minimal()

#there is a pattern of higher laterality before migration, but the sample sizes are too small. So, calculate laterality for each 8-second burst and look at the distributions there. Continued in L03b_tests_per_burst.r


