# script for testing the hypotheses for Safi et al 2025.
# Elham Nouani, PhD. 
#17.09.2024, Konstanz, DE


#PNAS figure guidelines: width for 1 column = 3.42, 1.5 columns = 4.5 in, 2 columns = 7 in; max height = 9 in . Text should be at least 6-8 pts

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(mapview)
library(viridis)
library(lme4)
library(mgcv)
library(mgcViz)
library(INLA)
library(terra)
library(performance)
library(corrr)
library(gridExtra)

setwd("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de2/Work/Projects/HB_ontogeny_eobs/R_files/")

#------------------------------------------------------------
## Step 1: calculate euler angles, and laterality index #####
#------------------------------------------------------------

#done in script: L03b_tests_per_burst.r: data summarized for each 8-sec burst
#life-cycle stage needs to be updated.

laterality_1sec_days <- readRDS("laterality_index_per_8sec_burst_days_since.rds") %>% 
  filter(n_records >= 8) %>% #remove short bursts
  mutate(days_since_tagging = as.numeric(days_since_tagging),
         individual_local_identifier = as.factor(individual_local_identifier)) %>% 
  mutate(circling_status = case_when(
    between(cumulative_yaw_8sec, -10, 10) ~ "straight",
    cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45 ~ "circling",
    .default = "shallow circling"
  )) %>% 
  as.data.frame()

#check to see how the circling status category worked
ggplot(laterality_1sec_days, aes(x = cumulative_yaw_8sec, fill = circling_status)) +
  geom_histogram(bins = 250)

#-------------------------------------------
## Step 2: filter data for flight only #####
#-------------------------------------------
#use pitch.... but pitch can be negative when gliding!!

#-------------------------------------------------------
## Step 3.1: Is there laterality? Population-level #####
#-------------------------------------------------------

# #data: angle summaries for each second (from 01b_imu_processing.r)


X11(width = 7, height = 2.8) 
(p <- ggplot(data = laterality_1sec_days, aes(x = mean_roll_mean, fill = circling_status)) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
    scale_x_continuous(breaks =seq(-90, 90, by = 10), 
                       labels = seq(-90, 90, by = 10),
                       limits = c(-92, 92)) +
    scale_fill_manual(values = c("straight" = "#40e0d0ff", "circling" = "#8a2be2ff", "shallow circling" = "#ff7f50ff")) +
    labs(x = "Bank angle (deg) averaged over each 8 second burst",
         y = "Count",
         fill = "Circularity") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)
ggsave(plot = p, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distribution.pdf", 
       device = "pdf", width = 7, height = 2.8, dpi = 600)

#-------------------------------------------------------
## Step 3.2: Is there laterality? Individual-level #####
#-------------------------------------------------------

#calculate laterality index for each individual. overall.
#based on Bennison et al: 
#right-side bias (LI = 1.0 to 0.25), a left side bias (LI = −1 to −0.25), no bias (LI = −0.25 to 0.25)

#large circles only
circling_only <- laterality_1sec_days %>% 
  filter(circling_status == "circling") %>% #only keep thermaling flight (>45 deg)
  select(-c("n_records", "bank_left", "bank_right", "bank_straight", "heading_left", #remove laterality columns for the burst-specific indices, because i'll calculate this for each individual
            "heading_right", "heading_straight", "laterality_bank", "laterality_heading"))

#calculate overall laterality index for each individual
ind_laterlaity_index <- circling_only %>%
  mutate(bank_direction = ifelse(mean_roll_mean < 0, "left",
                                 ifelse(mean_roll_mean > 0, "right", "straight")),
         heading_direction = ifelse(cumulative_yaw_8sec < 0, "left",
                                    ifelse(cumulative_yaw_8sec > 0, "right", "straight"))) %>% 
  group_by(individual_local_identifier) %>% 
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
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  as.data.frame()

#add laterality as a new column to the data
circling_w_LI <- circling_only %>% 
  full_join(ind_laterlaity_index, by = "individual_local_identifier") %>% 
  mutate(laterality_dir = factor(laterality_dir)) %>% 
  mutate(laterality_dir = fct_relevel(laterality_dir, c("left_handed", "ambidextrous", "right_handed"))) 


saveRDS(circling_w_LI, file = "circling_w_LI_population.rds")

#re-order the individuals based on the laterality index
circling_w_LI$individual_local_identifier <- reorder(circling_w_LI$individual_local_identifier, desc(circling_w_LI$laterality_dir))

X11(width = 7, height = 7)
(p_inds <- ggplot(data = circling_w_LI, aes(x = mean_roll_mean, y = individual_local_identifier, 
                                            color = laterality_dir, fill = laterality_dir),
                  height = stat(density)) +
    geom_density_ridges(stat = "binline", bins = 150, scale = 0.98, alpha = 0.8, draw_baseline = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
    scale_x_continuous(breaks = seq(-90, 90, by = 10), 
                       labels = seq(-90, 90, by = 10),
                       limits = c(-92, 92)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.03))) +  # Add extra space above the highest level
    scale_fill_manual(values = c("left_handed" = "royalblue" , "ambidextrous" = "yellowgreen", "right_handed" = "lightcoral")) +
    scale_color_manual(values = c("left_handed" = "royalblue", "ambidextrous" = "yellowgreen", "right_handed" = "lightcoral")) +
    labs(x = "Bank angle (deg) averaged over each 8 second burst",
         y = "Individual",
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
#+
#guides(color = "none"))

ggsave(plot = p_inds, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distr_circling_inds.pdf", 
       device = "pdf", width = 7, height = 7, dpi = 600)

#---------------------------------------------------------------------
## Step 4: Is laterality more likely when the task is difficult? #####
#---------------------------------------------------------------------
#logistic regression: handedness ~ tightness of circles * age . maybe only for individuals with laterality

#### ----------------------- extract the IDs of individuals with handedness
handed_IDs <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de2/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  filter(laterality_dir != "ambidextrous") %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  distinct(individual_local_identifier)


#### ----------------------- use the 8-sec data and calculate handedness for each 8 second burst
laterality_circling_thin <- readRDS("laterality_index_per_8sec_burst_days_since.rds") %>% 
  filter(n_records >= 8) %>% # & #remove short bursts) 
  #individual_local_identifier %in% handed_IDs$individual_local_identifier) %>% #only keep individuals with handedness
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = TRUE) %>% 
  #thin the data, so that the consecutive bursts are at least 30 seconds apart
  mutate(time_lag_sec = if_else(row_number() == 1, 0, 
                                difftime(start_timestamp, lag(start_timestamp), units = "secs") %>% as.numeric())) %>% 
  filter(time_lag_sec >= 18) %>% 
  ungroup() %>% 
  mutate(days_since_tagging = as.numeric(days_since_tagging),
         individual_local_identifier = as.factor(individual_local_identifier),
         individual_local_identifier2 = individual_local_identifier,
         individual_local_identifier3 = individual_local_identifier) %>% #dublicate individual ID to use for inla random effects specification
  mutate(circling_status = case_when(
    between(cumulative_yaw_8sec, -10, 10) ~ "straight",
    cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45 ~ "circling",
    .default = "shallow circling"
  )) %>% 
  mutate(laterality_dir = case_when(
    between(laterality_bank, 0.25, 1.0) ~ "right_handed",
    between(laterality_bank, -1.0, -0.25) ~ "left_handed",
    between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
    .default = NA)) %>% 
  filter(circling_status == "circling") %>% # only keep circling flight
  mutate(laterality_bi = ifelse(laterality_dir == "ambidextrous", 0, 1), #create a binary variable for handedness vs not
         abs_cum_yaw = abs(cumulative_yaw_8sec)) %>% 
  #make sure to do all the filtering before scaling the variables!!!!
  mutate_at(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
              "abs_cum_yaw", "days_since_tagging"),
            list(z = ~scale(.))) %>% 
  as.data.frame()


#### ----------------------- environmental annotation

#6940 rows don't have a gps location associated with them.... so, redo GPS-matching and increase the time window to one hour, to match that of the env. data
#open gps data: segmented and annotated (original code in )1b_imu_processing.r
gps_ls <- list.files("gps_seg_apr24", full.names = T) %>% #these are sf files
  map(readRDS) %>% 
  bind_rows() %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier),
         row_id = row_number(), #assign a row id to be able to cbind the annotated data with the original data later on
         location_long = st_coordinates(.)[,1],
         location_lat = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  arrange(individual_local_identifier) %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame) %>% 
  head(-1) #remove the last element, the 32nd animal

#create a list of one element for each individual
or_ls <- laterality_circling_thin %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  arrange(individual_local_identifier) %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

# Define a function to find the closest GPS information and associate it with orientation data
find_closest_gps <- function(or_data, gps_data, time_tolerance = 60 * 60) {
  map_df(1:nrow(or_data), function(h) {
    or_row_time <- or_data[h, "start_timestamp"]
    gps_sub <- gps_data %>%
      filter(between(timestamp, or_row_time - time_tolerance, or_row_time + time_tolerance))
    
    if (nrow(gps_sub) >= 1) {
      time_diff <- abs(difftime(gps_sub$timestamp, or_row_time, units = "secs"))
      min_diff <- which.min(time_diff)
      or_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps", "ground_speed_closest_gps", 
                   "heading_closest_gps", "flight_type_sm2", "flight_type_sm3", "track_flight_seg_id")] <- 
        gps_sub[min_diff, c("timestamp", "location_long", "location_lat", "height_above_ellipsoid", "ground_speed", "heading", 
                            "flight_clust_sm2", "flight_clust_sm3", "track_flight_seg_id")]
    } else {
      or_data[h, c("timestamp_closest_gps", "location_long_closest_gps", "location_lat_closest_gps", "height_above_ellipsoid_closest_gps", "ground_speed_closest_gps", 
                   "heading_closest_gps", "row_id", "flight_type_sm2", "flight_type_sm3", "track_flight_seg_id")] <- NA
    }
    return(or_data[h, ])
  })
}

#make sure the order of individuals is the same in the two lists
# Create a list of data frames with or data and associated GPS information
(b <- Sys.time())
or_w_gps <- map2(or_ls, gps_ls, ~ find_closest_gps(or_data = .x, gps_data = .y))
Sys.time() - b # 3.8 hrs for orientation with one sec per row

#add a column comparing or and gps timestamps. then save one file per individual. also limit to migratory season!! 1.09 - 30.10
or_w_gps_df <- lapply(or_w_gps, function(x){
  x2 <- x %>% 
    mutate(imu_gps_timediff_sec = if_else(is.na(timestamp_closest_gps), NA, difftime(start_timestamp, timestamp_closest_gps, units = "secs") %>%  as.numeric()))
  
  #saveRDS(x2, 
  #        file = paste0("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/matched_gps_quat/",
  #                      x2$individual_local_identifier[1], "_quat_w_gps.rds"))
  x2
}) %>% 
  bind_rows()

sum(is.na(or_w_gps_df$start_timestamp_closest_gps))
sum(is.na(or_w_gps_df$timestamp_closest_gps))

#there are still many rows with no assinged gps
no_gps <- or_w_gps_df %>% 
  filter(is.na(timestamp_closest_gps))

#####
#how many hours of data need to be downloaded?
laterality_circling_thin %>% 
  group_by(year(start_timestamp), month(start_timestamp), day(start_timestamp), hour(start_timestamp)) %>% 
  slice(1)



#### ----------------------- look at multi-collinearity
laterality_circling_thin %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
                  "abs_cum_yaw", "days_since_tagging")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #sd of roll and pitch = 0.65; sd of yaw and cumulative yaw = 0.7

#### ----------------------- exploratory plot
ggplot(laterality_circling_thin, aes(x = factor(laterality_bi), y = abs(cumulative_yaw_8sec))) +
  geom_boxplot() +
  labs(x = "Laterality", y = "Absolute Cumulative Yaw (8 sec)") +
  theme_minimal()

ggplot(laterality_circling_thin, aes(x = days_since_tagging, y = abs(cumulative_yaw_8sec))) +
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(laterality_circling_thin, aes(x = days_since_tagging, y = abs(cumulative_yaw_8sec))) +
  geom_point() + 
  geom_smooth(method = "loess")

ggplot(circling_thinned_100, aes(x = days_since_tagging, y = abs(cumulative_yaw_8sec))) +
  geom_point() + 
  geom_smooth(method = "gam")


ggplot(circling_thinned_100, aes(x = days_since_tagging, y = mean_pitch_mean)) +
  geom_point() + 
  geom_smooth(method = "gam")

ggplot(laterality_circling_thin, aes(x = abs(cumulative_yaw_8sec), y = mean_pitch_mean)) +
  geom_point() + 
  geom_smooth(method = "gam")


#### ----------------------- model: binomial logistic regression with inla

#model with mean pitch for the strength of the thermal
# m_inla <- inla(laterality_bi ~ 1 + abs_cum_yaw_z + mean_pitch_mean_z + days_since_tagging_z + 
#                  f(individual_local_identifier, abs_cum_yaw_z, model = "iid") +  
#                  f(individual_local_identifier2, days_since_tagging_z, model = "iid") ,
#                  f(individual_local_identifier3, mean_pitch_mean_z, model = "iid"),
#                data = data, family = "binomial",
#                control.compute = list(cpo = TRUE),
#                control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#age as a smooth term
#first bin the age variable
data <- laterality_circling_thin %>% 
  mutate(age_group = inla.group(days_since_tagging_z, n = 200, method = "quantile"))

m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z + abs_cum_yaw_z + 
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                 f(age_group, model = "rw1"),
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#### ----------------------- model evaluation

#look at residuals
residuals <- data$laterality_bi - m_inla$summary.fitted.values$mean
plot(residuals)

#calculate variance explained: McFadden's pseudo R²
log_likelihood_full <- m_inla$mlik[1, 1]

null_model <- inla(laterality_bi ~ 1, family = "binomial",
                   control.compute = list(cpo = TRUE),
                   control.predictor = list(link = 1, compute = TRUE),
                   data = data)

log_likelihood_null <- null_model$mlik[1, 1]

pseudo_r_squared <- 1 - (log_likelihood_full / log_likelihood_null) #0.006

#model validation metrics
eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), # 0.52
                   Mlik = as.numeric(m_inla$mlik[,1][2])) # -15300.69

#### coefficients plot -----------------------------------------------------------------------------

# posterior means of coefficients
graph <- as.data.frame(summary(m_inla)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)

#remove weeks since dispersal
VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

#graph$Factor_n <- as.numeric(graph$Factor)

#plot the coefficients
X11(width = 3.42, height = 1.5)
(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray", linewidth = 0.5) +
    geom_point(color = "#8a2be2ff", size = 2)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Intercept", "Average pitch", "Absolute cumulative yaw"))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#8a2be2ff", linewidth = 0.5) +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

ggsave(plot = coefs, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_coeffs.pdf", 
       device = "pdf", width = 3.42, height = 1.5, dpi = 600)

#### ind_specific coefficients plot -----------------------------------------------------------------------------


#extract handedness for each individual. to use for coloring
handedness <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de2/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir)

#extract random effects,  ID is for the unique individuals
#mean_pitch_mean_z
random_effects_pitch <- m_inla$summary.random$individual_local_identifier

#abs_cum_yaw_z
random_effects_yaw <- m_inla$summary.random$individual_local_identifier2

#extract unique individual IDs from original data
ind_IDs <- unique(data$individual_local_identifier)

pitch <- random_effects_pitch %>% 
  mutate(coef = mean + graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Upper),
         variable = "mean_pitch_mean_z")

yaw <- random_effects_yaw %>% 
  mutate(coef = mean + graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Upper),
         variable = "abs_cum_yaw_z")

two_vars <- bind_rows(pitch, yaw) %>%
  full_join(handedness, by = c("ID" = "individual_local_identifier")) %>% 
  mutate(v_line = ifelse(variable == "mean_pitch_mean_z", graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate),
                         graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate)))

#re-order the individuals based on the laterality index
two_vars$ID <- reorder(two_vars$ID, desc(two_vars$laterality_dir))

X11(width = 7, height = 9)
(coefs_inds <- ggplot(two_vars, aes(x = coef, y = ID, color = laterality_dir)) +
    geom_vline(data = filter(two_vars, variable == "mean_pitch_mean_z"), 
               aes(xintercept = graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75",size = 0.5) + 
    geom_vline(data = filter(two_vars, variable == "abs_cum_yaw_z"), 
               aes(xintercept = graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", size = 0.5) +  
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = c("left_handed" = "royalblue", "ambidextrous" = "yellowgreen", "right_handed" = "lightcoral"),
                       name = "Handedness") +
    scale_y_discrete(labels = ind_IDs) +
    labs(x = "Estimate", y = "") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) + #increase distance between x-axis values and title
    facet_wrap(~ variable, scales = "free_x", labeller = as_labeller(c(
      "mean_pitch_mean_z" = "Average pitch",
      "abs_cum_yaw_z" = "Absolute cumulative yaw")
    )) # Separate panels for each variable
)

ggsave(plot = coefs_inds, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_coeffs_inds.pdf", 
       device = "pdf", width = 7, height = 9, dpi = 600)


#### plot the smooth term -----------------------------------------------------------------------------

# Extract the summary of the smooth term 
smooth_effects <- m_inla$summary.random$age_group %>% 
  #back transform the age values
  mutate(age = ID * attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "days_since_tagging_z"],'scaled:scale') +
           attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "days_since_tagging_z"],'scaled:center') )

# Plot the smooth term.... HIGHLIGHT MIGRATION period
X11(width = 3.42, height = 1.5)
(s <- ggplot(smooth_effects, aes(x = age, y = mean)) +
    geom_line(color = "#8a2be2ff", linewidth = 0.5) +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "#8a2be2ff", alpha = 0.12) +
    # xlim(1,300) +
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "gray", linewidth = 0.5) +
    labs(y = "Effect Size", x = "Days since tagging") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))
)


#combine the two plots for the linear and the non-linear terms
X11(width = 7, height = 1.7)
model_output_p <- grid.arrange(coefs, s, nrow = 1)

ggsave(plot = model_output_p, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_output.pdf", 
       device = "pdf", width = 7, height = 2, dpi = 600)



# Assuming 'data' is your dataset and 'm_inla' is your fitted model
# Extract predictions from the model
# predictions <- m_inla$summary.fitted.values
# 
# # Add predictions to your data
# data$predicted_prob <- predictions$mean

# # Plot the relationship between the probability of laterality and days since tagging
# ggplot(data, aes(x = days_since_tagging, y = predicted_prob)) +
#   geom_line(color = "blue") +
#   geom_point(aes(y = laterality_bi), alpha = 0.5, color = "red") +
#   labs(title = "Probability of Laterality vs Days Since Tagging",
#        x = "Days Since Tagging",
#        y = "Probability of Laterality") +
#   theme_minimal()


#-----------------------------------------------------------------------------------------------------------------------
## Step 5.1: Does laterality help with better performance when individuals are not experienced? Flight performance #####
#-----------------------------------------------------------------------------------------------------------------------
#flight performance ~ level of handedness * age




#-----------------------------------------------------------------------------------------------------------------------
## Step 5.2: Does laterality help with better performance when individuals are not experienced? migration performance #####
#-----------------------------------------------------------------------------------------------------------------------
#migration performance ~ level of handedness  (one row per individual)


