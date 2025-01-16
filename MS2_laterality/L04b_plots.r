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


#open filtered data with LI calculated for each day and each life_stage
filtered_w_LI <- readRDS("thinned_laterality_w_gps_wind_all_filters.rds")



#use Ellen's life cycle stage estimations
life_cycle <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/HB_dates_Fin_all_M2.rds") %>% 
  mutate_at(c("migration_start", "migration_end"), as.Date) 

###----------------------------------###
### FOR DAYS   all flight types      ###
###----------------------------------###

#open data with daily laterality index and life cycle stages. from L03a_tests_per_day.r
#this is not filtered for circling flight only

daily_laterality <- readRDS("laterality_daily_w_LS.rds") %>%
  #assign handedness categories for coloring
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  #group by life-cycle stage and assign the mode of laterality as the laterality index for that stage
  group_by(individual_local_identifier, life_stage) %>% 
  mutate(mode_laterality_stage = getmode(laterality_dir) %>% factor(),
         mode_laterality_stage = fct_relevel(mode_laterality_stage, c("left_handed", "ambidextrous", "right_handed")),
         life_stage = factor(life_stage, levels = c("post-fledging", "migration", "wintering"))) %>% 
  ungroup()

#reorder based on the value of handedness during post-fledging
# Filter the data for the first level of life_stage
subset_first_stage <- daily_laterality %>%
  filter(life_stage == "post-fledging")

# Determine the order of individual_local_identifier based on laterality_dir
ordered_identifiers <- daily_laterality %>%
  filter(life_stage == "post-fledging") %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  select(individual_local_identifier, mode_laterality_stage) %>% 
  arrange(mode_laterality_stage) %>% 
  ungroup() %>% 
  pull(individual_local_identifier)


# Reorder the factor levels of individual_local_identifier
daily_laterality$individual_local_identifier <- factor(daily_laterality$individual_local_identifier, 
                                                       levels = ordered_identifiers)


#### ridgelines for laterality in banking angle
(p <- ggplot(daily_laterality, aes(x = laterality_bank, y = individual_local_identifier, color = mode_laterality_stage, fill = mode_laterality_stage)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("left_handed" = "#33638DFF" ,  "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    # labels = c("Left-handed \nLI = -1.0 to -0.25", "Ambidextrous \nLI = -0.25 to 0.25", "Right-handed \nLI = 0.25 to 1.0")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    ggtitle("Daily laterality index during different life cycle stages. Not filtered for flight type. \nright-side bias (LI = 1.0 to 0.25), left side bias (LI = −0.25 to -1), no bias (LI = −0.25 to 0.25)") +
    facet_wrap(vars(life_stage)) +
    labs(x = "Laterality index",
         y = "Individual ID",
         fill = "Handedness \n(mode of laterality category \nfor each life-cycle stage)",
         color = "Handedness \n(mode of laterality category \nfor each life-cycle stage)") +
    theme_minimal() +
    theme(legend.key.size = unit(.9, "cm"),
          plot.title = element_text(size = 11))
)

ggsave(plot = p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_lat_days_LS.jpg", 
       device = "jpg", width = 9, height = 9)

#conclusion was: we don't have that many days.

###----------------------------------###
### PER BURST. all flight types      ###
###----------------------------------###

#life-cycle stages from L03a_tests_per_day.r
life_cycle <- readRDS("updated_life_cycle_nov24.rds")

#add life stage to burst-level laterality data... 
laterality_per_burst <-  readRDS("laterality_index_per_8sec_burst_days_since.rds") %>%  #from L03b_tests_per_burst.r
  drop_na(individual_local_identifier) %>% 
  mutate(unique_date = as.Date(start_timestamp)) %>% 
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
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA),
    #add circling status
    circling_status = case_when(
      between(cumulative_yaw_8sec, -10, 10) ~ "straight",
      cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45 ~ "circling",
      .default = "shallow circling")) %>% 
  #group by life-cycle stage and assign the mode of laterality as the laterality index for that stage
  group_by(individual_local_identifier, life_stage) %>% 
  mutate(mode_laterality_stage = getmode(laterality_dir) %>% factor(),
         mode_laterality_stage = fct_relevel(mode_laterality_stage, c("left_handed", "ambidextrous", "right_handed")),
         life_stage = factor(life_stage, levels = c("post-fledging", "migration", "wintering"))) %>% 
  ungroup()


# Reorder the factor levels of individual_local_identifier based on the previous plot
laterality_per_burst$individual_local_identifier <- factor(laterality_per_burst$individual_local_identifier, 
                                                           levels = ordered_identifiers)


#### ridgelines for laterality in banking angle
(p2 <- ggplot(laterality_per_burst, aes(x = laterality_bank, y = individual_local_identifier, color = mode_laterality_stage, fill = mode_laterality_stage)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("left_handed" = "#33638DFF" ,  "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    # labels = c("Left-handed \nLI = -1.0 to -0.25", "Ambidextrous \nLI = -0.25 to 0.25", "Right-handed \nLI = 0.25 to 1.0")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    ggtitle("Laterality index for each 8-sec burst. Not filtered for flight type. \nright-side bias (LI = 1.0 to 0.25), left side bias (LI = −0.25 to -1), no bias (LI = −0.25 to 0.25)") +
    facet_wrap(vars(life_stage)) +
    labs(x = "Laterality index",
         y = "Individual ID",
         fill = "Handedness \n(mode of laterality category \nfor each life-cycle stage)",
         color = "Handedness \n(mode of laterality category \nfor each life-cycle stage)") +
    theme_minimal() +
    theme(legend.key.size = unit(.9, "cm"),
          plot.title = element_text(size = 11))
)


ggsave(plot = p2, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_lat_bursts_LS.jpg", 
       device = "jpg", width = 9, height = 9)


#### ridgelines for laterality in raw banking angle
(p3 <- ggplot(laterality_per_burst, aes(x = mean_roll_mean, y = individual_local_identifier, color = mode_laterality_stage, fill = mode_laterality_stage)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    xlim(-90,90) +
    scale_fill_manual(values = c("left_handed" = "#33638DFF" ,  "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    # labels = c("Left-handed \nLI = -1.0 to -0.25", "Ambidextrous \nLI = -0.25 to 0.25", "Right-handed \nLI = 0.25 to 1.0")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    ggtitle("Bank angle at each 8-sec burst. Not filtered for flight type. \nright-side bias (LI = 1.0 to 0.25), left side bias (LI = −0.25 to -1), no bias (LI = −0.25 to 0.25)") +
    facet_wrap(vars(life_stage)) +
    labs(x = "Bank angle",
         y = "Individual ID",
         fill = "Handedness \n(mode of laterality category \nfor each life-cycle stage)",
         color = "Handedness \n(mode of laterality category \nfor each life-cycle stage)") +
    theme_minimal() +
    theme(legend.key.size = unit(.9, "cm"),
          plot.title = element_text(size = 11))
)

ggsave(plot = p3, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_roll_bursts_LS.jpg", 
       device = "jpg", width = 9, height = 9)



###----------------------------------###
### PER BURST. only circling         ###
###----------------------------------###


circling_only <- laterality_per_burst %>% 
  filter(circling_status == "circling") %>% 
  #recalculate mode of handedness for this subset of the data
  group_by(individual_local_identifier, life_stage) %>% 
  mutate(mode_laterality_stage = getmode(laterality_dir) %>% factor(),
         mode_laterality_stage = fct_relevel(mode_laterality_stage, c("left_handed", "ambidextrous", "right_handed")),
         life_stage = factor(life_stage, levels = c("post-fledging", "migration", "wintering"))) %>% 
  ungroup() %>% 
  as.data.frame()

circling_data <

# Reorder the factor levels of individual_local_identifier based on the previous plots
circling_only$individual_local_identifier <- factor(circling_only$individual_local_identifier, 
                                                    levels = ordered_identifiers)



(p4 <- ggplot(circling_only, aes(x = laterality_bank, y = individual_local_identifier, color = mode_laterality_stage, fill = mode_laterality_stage)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("left_handed" = "#33638DFF" ,  "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    # labels = c("Left-handed \nLI = -1.0 to -0.25", "Ambidextrous \nLI = -0.25 to 0.25", "Right-handed \nLI = 0.25 to 1.0")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    ggtitle("Laterality index for each 8-sec burst. Circling flight. \nright-side bias (LI = 1.0 to 0.25), left side bias (LI = −0.25 to -1), no bias (LI = −0.25 to 0.25)") +
    facet_wrap(vars(life_stage)) +
    labs(x = "Laterality index",
         y = "Individual ID",
         fill = "Handedness \n(LI calculated for each life-cycle stage)",
         color = "Handedness \n(LI calculated for each life-cycle stage)") +
    theme_minimal() +
    theme(legend.key.size = unit(.9, "cm"),
          plot.title = element_text(size = 11))
)


ggsave(plot = p4, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_lat_bursts_LS_c.jpg", 
       device = "jpg", width = 9, height = 9)


#### ridgelines for laterality in raw banking angle
(p5 <- ggplot(circling_only, aes(x = mean_roll_mean, y = individual_local_identifier, color = mode_laterality_stage, fill = mode_laterality_stage)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    xlim(-90,90) +
    scale_fill_manual(values = c("left_handed" = "#33638DFF" ,  "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    # labels = c("Left-handed \nLI = -1.0 to -0.25", "Ambidextrous \nLI = -0.25 to 0.25", "Right-handed \nLI = 0.25 to 1.0")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    ggtitle("Bank angle at each 8-sec burst. Circling flight. \nright-side bias (LI = 1.0 to 0.25), left side bias (LI = −0.25 to -1), no bias (LI = −0.25 to 0.25)") +
    facet_wrap(vars(life_stage)) +
    labs(x = "Bank angle",
         y = "Individual ID",
         fill = "Handedness \n(mode of laterality category \nfor each life-cycle stage)",
         color = "Handedness \n(mode of laterality category \nfor each life-cycle stage)") +
    theme_minimal() +
    theme(legend.key.size = unit(.9, "cm"),
          plot.title = element_text(size = 11))
)

ggsave(plot = p5, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_roll_bursts_LS_c.jpg", 
       device = "jpg", width = 9, height = 9)



###-----------------------------------------------###
### PER BURST. only circling. LI for life stages  ###
###-----------------------------------------------###

#calculate overall laterality index for each individual at each life cycle stage
ind_ls_laterlaity_index <- laterality_per_burst %>%
  #filter for circling flight
  filter(circling_status == "circling") %>% 
  mutate(bank_direction = ifelse(mean_roll_mean < 0, "left",
                                 ifelse(mean_roll_mean > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_8sec < 0, "left",
                                    ifelse(cumulative_yaw_8sec > 0, "right", "straight"))) %>% 
  group_by(individual_local_identifier, life_stage) %>% 
  summarise(total_n = n(),
            bank_left = sum(bank_direction == "left"),
            bank_right = sum(bank_direction == "right"),
            bank_straight = sum(bank_direction == "straight"),
            heading_left = sum(heading_direction == "left"),
            heading_right = sum(heading_direction == "right"),
            heading_straight = sum(heading_direction == "straight"),
            laterality_bank = (bank_right - bank_left)/(bank_right + bank_left),
            laterality_heading = (heading_right - heading_left)/(heading_right + heading_left),
            .groups = 'drop') %>% 
  mutate(laterality_dir_LS = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  as.data.frame()

#add stage-specific laterality index as a new column to the data
circling_w_LI <- laterality_per_burst %>% 
  #filter for circling flight
  filter(circling_status == "circling") %>% 
  full_join(ind_ls_laterlaity_index %>% select(individual_local_identifier, life_stage, laterality_dir_LS), by = c("individual_local_identifier", "life_stage")) %>% 
  mutate(laterality_dir_LS = factor(laterality_dir_LS)) %>% 
  mutate(laterality_dir_LS = fct_relevel(laterality_dir_LS, c("left_handed", "ambidextrous", "right_handed"))) 


# Reorder the factor levels of individual_local_identifier based on the previous plots
circling_w_LI$individual_local_identifier <- factor(circling_w_LI$individual_local_identifier, 
                                                    levels = ordered_identifiers)



(p4 <- ggplot(circling_w_LI, aes(x = laterality_bank, y = individual_local_identifier, color = laterality_dir_LS, fill = laterality_dir_LS)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c("left_handed" = "#33638DFF" ,  "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    # labels = c("Left-handed \nLI = -1.0 to -0.25", "Ambidextrous \nLI = -0.25 to 0.25", "Right-handed \nLI = 0.25 to 1.0")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    ggtitle("Laterality index for each 8-sec burst. Circling flight. \nright-side bias (LI = 1.0 to 0.25), left side bias (LI = −0.25 to -1), no bias (LI = −0.25 to 0.25)") +
    facet_wrap(vars(life_stage)) +
    labs(x = "Laterality index",
         y = "Individual ID",
         fill = "Handedness \n(LI calculated for each \nlife-cycle stage)",
         color = "Handedness \n(LI calculated for each \nlife-cycle stage)") +
    theme_minimal() +
    theme(legend.key.size = unit(.9, "cm"),
          plot.title = element_text(size = 11))
)


ggsave(plot = p4, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_lat_bursts_LS_LI_c.jpg", 
       device = "jpg", width = 9, height = 9)


#### ridgelines for laterality in raw banking angle
(p5 <- ggplot(circling_w_LI, aes(x = mean_roll_mean, y = individual_local_identifier, color = laterality_dir_LS, fill = laterality_dir_LS)) +
    stat_density_ridges(quantile_lines = TRUE, rel_min_height = 0.01, alpha = 0.5,
                        jittered_points = TRUE, 
                        point_shape = "|", point_size = 1, point_alpha = 1, size = 0.2) + #Trailing tails can be cut off using the rel_min_height aesthetic.
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    xlim(-90,90) +
    scale_fill_manual(values = c("left_handed" = "#33638DFF" ,  "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    # labels = c("Left-handed \nLI = -1.0 to -0.25", "Ambidextrous \nLI = -0.25 to 0.25", "Right-handed \nLI = 0.25 to 1.0")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    ggtitle("Bank angle at each 8-sec burst. Circling flight. \nright-side bias (LI = 1.0 to 0.25), left side bias (LI = −0.25 to -1), no bias (LI = −0.25 to 0.25)") +
    facet_wrap(vars(life_stage)) +
    labs(x = "Bank angle",
         y = "Individual ID",
         fill = "Handedness \n(mode of laterality category \nfor each life-cycle stage)",
         color = "Handedness \n(mode of laterality category \nfor each life-cycle stage)") +
    theme_minimal() +
    theme(legend.key.size = unit(.9, "cm"),
          plot.title = element_text(size = 11))
)

ggsave(plot = p5, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_roll_bursts_LS_LI_c.jpg", 
       device = "jpg", width = 9, height = 9)


###-----------------------------------------------###
### PER BURST. only circling. NO LS. LI per ind   ###
###-----------------------------------------------###

#recreate figure 2, but with the same order as the plots in this script

#circling_w_LI_inds <- readRDS("circling_w_LI_population.rds")


# Reorder the factor levels of individual_local_identifier based on the previous plots
circling_w_LI_inds$individual_local_identifier <- factor(circling_w_LI_inds$individual_local_identifier, 
                                                    levels = ordered_identifiers)


X11(width = 7, height = 7)
(p_inds <- ggplot(data = circling_w_LI_inds, aes(x = mean_roll_mean, y = individual_local_identifier, 
                                            color = laterality_dir, fill = laterality_dir),
                  height = stat(density)) +
    geom_density_ridges(stat = "binline", bins = 150, scale = 0.98, alpha = 0.8, draw_baseline = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_x_continuous(breaks = seq(-90, 90, by = 10), 
                       labels = seq(-90, 90, by = 10),
                       limits = c(-92, 92)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.03))) +  # Add extra space above the highest level
    scale_fill_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#B8DE29FF", "right_handed" = "#20A387FF"),
                      labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    labs(x = "Bank angle (deg) averaged over each 8 second burst",
         y = "Individual ID",
         fill = "Handedness",
         color = "Handedness") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

ggsave(plot = p_inds, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/exploration_figs/ind_roll_fig2_reordered.jpg", 
       device = "jpg", width = 9, height = 9)


