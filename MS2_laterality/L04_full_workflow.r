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

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

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

laterality_1sec_days <- laterality_1sec_days %>% 
  mutate(circling_status = as.factor(circling_status))

X11(width = 7, height = 2.8) 
(p <- ggplot(data = laterality_1sec_days, aes(x = mean_roll_mean, fill = circling_status)) +
    geom_histogram(binwidth = 1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_x_continuous(breaks =seq(-90, 90, by = 10), 
                       labels = seq(-90, 90, by = 10),
                       limits = c(-92, 92)) +
    scale_fill_manual(values = c("straight" = "#238A8DFF", "circling" = "#8a2be2ff", "shallow circling" = "#DCE319FF"),
                      labels = rev(c("Straight", "Shallow circling", "Circling"))) +
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

ggsave(plot = p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distribution.pdf", 
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
ggsave(plot = p_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distr_circling_inds.pdf", 
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

#### ----------------------- read in data on individuals' age at tagging (from Ellen's MSc thesis)
ages <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/aging_ellen-msc.csv")

#mean(ages$Estimated.age)
#1] 31.70968

#### ----------------------- use the 8-sec data and calculate handedness for each 8 second burst
laterality_circling_thin <- readRDS("laterality_index_per_8sec_burst_days_since.rds") %>% 
  filter(n_records >= 8) %>% # & #remove short bursts) 
  #full_join(ages, by = "individual_local_identifier") #do this when I have the correct ages from Ellen.
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
  as.data.frame()


#### ----------------------- environmental annotation

#6940 rows don't have a gps location associated with them.... so, redo GPS-matching and increase the time window to one hour, to match that of the env. data
#open raw gps data, matched with wind in L04a_env_annotation.r

gps_ls <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/all_gps_apr15_24_wind.rds") %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  arrange(individual_local_identifier) %>% 
  group_by(individual_local_identifier) %>% 
  group_split() %>% 
  map(as.data.frame)

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
      or_data[h, c("timestamp_closest_gps_raw", "location_long_closest_gps_raw", "location_lat_closest_gps_raw", "height_above_ellipsoid_closest_gps_raw", "ground_speed_closest_gps_raw", 
                   "heading_closest_gps",  "u_900", "v_900", "wind_speed")] <- 
        gps_sub[min_diff, c("timestamp", "location_long", "location_lat", "height_above_ellipsoid", "ground_speed", "heading", "u_900", "v_900", "wind_speed")]
    } else {
      or_data[h, c("timestamp_closest_gps_raw", "location_long_closest_gps_raw", "location_lat_closest_gps_raw", "height_above_ellipsoid_closest_gps_raw", "ground_speed_closest_gps_raw", 
                   "heading_closest_gps", "u_900", "v_900", "wind_speed")] <- NA
    }
    return(or_data[h, ])
  })
}

#make sure the order of individuals is the same in the two lists
# Create a list of data frames with or data and associated GPS information
(b <- Sys.time())
or_w_gps <- map2(or_ls, gps_ls, ~ find_closest_gps(or_data = .x, gps_data = .y))
Sys.time() - b # 2.5 mins

#add a column comparing or and gps timestamps. then save one file per individual.
or_w_gps_df <- lapply(or_w_gps, function(x){
  x2 <- x %>% 
    mutate(imu_gps_timediff_sec = if_else(is.na(timestamp_closest_gps_raw), NA, difftime(start_timestamp, timestamp_closest_gps_raw, units = "mins") %>%  as.numeric()))
  x2
}) %>% 
  bind_rows()

sum(is.na(or_w_gps_df$location_long_closest_gps_raw)) #406
sum(is.na(or_w_gps_df$start_location_long_closest_gps)) #6940

saveRDS(or_w_gps_df, file = "thinned_laterality_w_gps_wind.rds")

#there are still many rows with no assigned gps
no_gps <- or_w_gps_df %>% 
  filter(is.na(timestamp_closest_gps_raw))

#### ----------------------- filter for day since tagging and z-transform

quantile(or_w_gps_df$days_since_tagging, probs = 0.9) #271

circling_data <- or_w_gps_df %>% 
  #filter for days since tagging. only keep the 
  filter(days_since_tagging < quantile(days_since_tagging, probs = 0.9)) %>% 
  #make sure to do all the filtering before scaling the variables!!!!
  mutate_at(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
              "abs_cum_yaw", "wind_speed", "days_since_tagging"),
            list(z = ~scale(.))) %>% 
  as.data.frame()

saveRDS(circling_data, file = "thinned_laterality_for_modeling.rds")

#### ----------------------- look at multi-collinearity

circling_data %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
                  "abs_cum_yaw", "days_since_tagging", "wind_speed")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #sd of roll and pitch = 0.64; sd of yaw and cumulative yaw = 0.7

#### ----------------------- exploratory plot
ggplot(circling_data, aes(x = factor(laterality_bi), y = abs(cumulative_yaw_8sec))) +
  geom_boxplot() +
  labs(x = "Laterality", y = "Absolute Cumulative Yaw (8 sec)") +
  theme_minimal()

ggplot(circling_data, aes(x = days_since_tagging, y = abs(cumulative_yaw_8sec))) +
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(circling_data, aes(x = days_since_tagging, y = abs(cumulative_yaw_8sec))) +
  geom_point() + 
  geom_smooth(method = "loess")

ggplot(circling_data, aes(x = days_since_tagging, y = abs(cumulative_yaw_8sec))) +
  geom_point() + 
  geom_smooth(method = "gam")


ggplot(circling_data, aes(x = days_since_tagging, y = mean_pitch_mean)) +
  geom_point() + 
  geom_smooth(method = "gam")

ggplot(circling_data, aes(x = abs(cumulative_yaw_8sec), y = mean_pitch_mean)) +
  geom_point() + 
  geom_smooth(method = "gam")


#### ----------------------- model: binomial logistic regression with inla

circling_data <- readRDS("thinned_laterality_for_modeling.rds")


#create new data: make predictions for values that fall on a regular grid for optimal visualization 

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations. The max of the variables might be outliers, so use the 90% quantile instead
grd_pitch_yaw <- expand.grid(x = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = max(circling_data$mean_pitch_mean, na.rm = T),  length.out = 50),
                             y = seq(from = quantile(circling_data$abs_cum_yaw, .01, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(mean_pitch_mean = x,
         abs_cum_yaw = y)  %>% 
  mutate(wind_speed = attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "pitch_yaw")

grd_wind_yaw <- expand.grid(x = seq(from = quantile(circling_data$wind_speed, .01, na.rm = T), to = quantile(circling_data$wind_speed, .99, na.rm = T),  length.out = 50),
                            y = seq(from = quantile(circling_data$abs_cum_yaw, .01, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(wind_speed = x,
         abs_cum_yaw = y)  %>% 
  mutate(mean_pitch_mean = attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "wind_yaw")

grd_wind_pitch <- expand.grid(x = seq(from = quantile(circling_data$wind_speed, .01, na.rm = T), to = quantile(circling_data$wind_speed, .99, na.rm = T),  length.out = 50),
                              y = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = quantile(circling_data$mean_pitch_mean, .99, na.rm = T),  length.out = 50)) %>% 
  rename(wind_speed = x,
         mean_pitch_mean = y)  %>% 
  mutate(abs_cum_yaw = attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "wind_pitch")


#merge all together
grd_all <- bind_rows(grd_pitch_yaw, grd_wind_yaw, grd_wind_pitch) %>% 
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation for consistency
  mutate(mean_pitch_mean_z = (mean_pitch_mean - attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:scale'),
         abs_cum_yaw_z = (abs_cum_yaw - attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "abs_cum_yaw_z"],'scaled:scale'),
         wind_speed_z = (wind_speed - attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:scale'),
         days_since_tagging_z = (days_since_tagging - attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:scale'),
         individual_local_identifier = sample(circling_data$individual_local_identifier, nrow(.), replace = T) %>% as.factor(),
         laterality_bi = NA) 


#bin the age variable and append the new data
data <- circling_data %>% 
  dplyr::select(intersect(colnames(grd_all), colnames(.))) %>% 
  full_join(grd_all) %>% 
  mutate(individual_local_identifier2 = individual_local_identifier, #repeat individual ID column to be used in the model formula for random effects
         individual_local_identifier3 = individual_local_identifier,
         age_group = inla.group(days_since_tagging_z, n = 100, method = "quantile")) #age will be included as a smooth term

saveRDS(data, file = "thinned_laterality_for_modeling_w_new_data.rds")

m_inla1 <- inla(laterality_bi ~ 1 + mean_pitch_mean_z + abs_cum_yaw_z + wind_speed_z +
                  f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                  f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                  f(individual_local_identifier3, wind_speed_z, model = "iid") + 
                  f(age_group, model = "rw1"),
                data = data, family = "binomial",
                control.compute = list(cpo = TRUE),
                control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z * wind_speed_z +
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                 f(individual_local_identifier3, wind_speed_z, model = "iid") + 
                 f(age_group, model = "rw1"),
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#### model evaluation -----------------------

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
                   Mlik = as.numeric(m_inla$mlik[,1][2])) # -13888

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
X11(width = 7, height = 3)
(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#8a2be2ff", size = 2)  +
    labs(x = "Estimate", y = "") +
    #scale_y_discrete(labels = rev(c("Intercept", "Average pitch", "Absolute cumulative yaw"))) +
    scale_y_discrete(labels = rev(c("Intercept", "Average pitch", "Absolute cumulative yaw",
                                    "Wind speed", "Average pitch: Absolute cumulative yaw", "Average pitch: Wind speed",
                                    "Average cumulative yaw: Wind speed", "Average pitch : Average cumulative yaw: \nWind speed "))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#8a2be2ff", linewidth = 0.5) +
    ggtitle("a") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

ggsave(plot = coefs, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_coeffs.pdf", 
       device = "pdf", width = 7, height = 4, dpi = 600)

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

#wind_speed_z
random_effects_wind <- m_inla$summary.random$individual_local_identifier3

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

wind <- random_effects_wind %>% 
  mutate(coef = mean + graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "wind_speed_z") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "wind_speed_z") %>% pull(Upper),
         variable = "wind_speed_z")

three_vars <- bind_rows(pitch, yaw, wind) %>%
  full_join(handedness, by = c("ID" = "individual_local_identifier")) %>% 
  mutate(v_line = case_when(
    variable == "mean_pitch_mean_z" ~ graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate),
    variable == "abs_cum_yaw_z" ~ graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate),
    variable == "wind_speed_z" ~ graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate)
  ))


#re-order the individuals based on the laterality index
three_vars$ID <- reorder(three_vars$ID, desc(three_vars$laterality_dir))

X11(width = 9, height = 5)
(coefs_inds <- ggplot(three_vars, aes(x = coef, y = ID, color = laterality_dir)) +
    geom_vline(data = filter(three_vars, variable == "mean_pitch_mean_z"), 
               aes(xintercept = graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75",size = 0.5) + 
    geom_vline(data = filter(three_vars, variable == "abs_cum_yaw_z"), 
               aes(xintercept = graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", size = 0.5) +  
    geom_vline(data = filter(three_vars, variable == "wind_speed_z"), 
               aes(xintercept = graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", size = 0.5) +  
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed"),
                       name = "Handedness") +
    scale_y_discrete(labels = ind_IDs) +
    labs(x = "Estimate", y = "Individual ID") +
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
      "abs_cum_yaw_z" = "Absolute cumulative yaw",
      "wind_speed_z" = "Wind speed"
    )) # Separate panels for each variable
    ))

ggsave(plot = coefs_inds, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_coeffs_inds.pdf", 
       device = "pdf", width = 9, height = 5, dpi = 600)


#### plot the smooth term -----------------------------------------------------------------------------

# Extract the summary of the smooth term 
smooth_effects <- m_inla$summary.random$age_group %>% 
  #back transform the age values
  mutate(age = ID * attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:scale') +
           attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:center') )

# Plot the smooth term.... HIGHLIGHT MIGRATION period
X11(width = 3.42, height = 3)
(s <- ggplot(smooth_effects, aes(x = age, y = mean)) +
    geom_line(color = "#8a2be2ff", linewidth = 0.3) +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "#8a2be2ff", alpha = 0.12) +
    # xlim(1,300) +
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    labs(y = "Effect Size", x = "Days since tagging") +
    ggtitle("b") +
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
X11(width = 9, height = 3)
model_output_p <- grid.arrange(coefs, s, nrow = 1, widths = c(0.65, 0.35))

ggsave(plot = model_output_p, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_output.pdf", 
       device = "pdf", width = 9, height = 3, dpi = 600)


#### plot the interaction terms -----------------------------------------------------------------------------

terms <- unique(data$interaction)[-1] #remove the NA

for(i in terms){
  if(i == "pitch_yaw"){
    
    #extract information for rows that had NAs as response variables
    na_rows <- which(data$interaction == i)
    
    preds <- data.frame(yaw = data[na_rows,"abs_cum_yaw"],
                        pitch = data[na_rows,"mean_pitch_mean"],
                        preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
      #do some interpolation to smooth the plot
      terra::rast(type = "xyz") %>%
      focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
      as.data.frame(xy = T) %>%
      rename(preds = focal_mean)
    
    #plot
    X11(width = 7, height = 2)
    (pred_py <- preds %>% 
        ggplot() +
        geom_tile(aes(x = x, y = y, fill = preds)) +
        scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                             values = c(0.2, 0.5, 0.6, 0.7, 0.9),
                             limits = c(0, 1),
                             na.value = "white",
                             name = "Probability of Laterality") +
        guides(fill = guide_colourbar(title.vjust = .95)) + # the legend title needs to move up a bit
        labs(x = "Absolute cumulative yaw", y = "Average pitch") +
        ggtitle("a") +
        theme_classic() +
        theme(text = element_text(size = 11),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 9),
              panel.grid.minor = element_line(color = "white"),
              axis.title.x = element_text(margin = margin(t = 5))) + # increase distance between x-axis values and title
        scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
        scale_y_continuous(expand = c(0, 0))
    )
  }else if(i == "wind_yaw"){
    
    #extract information for rows that had NAs as response variables
    na_rows <- which(data$interaction == i)
    
    preds <- data.frame(yaw = data[na_rows,"abs_cum_yaw"],
                        pitch = data[na_rows,"mean_pitch_mean"],
                        preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
      #do some interpolation to smooth the plot
      terra::rast(type = "xyz") %>%
      focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
      as.data.frame(xy = T) %>%
      rename(preds = focal_mean)
    
    #plot
    X11(width = 7, height = 2)
    (pred_py <- preds %>% 
        ggplot() +
        geom_tile(aes(x = x, y = y, fill = preds)) +
        scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                             values = c(0.2, 0.5, 0.6, 0.7, 0.9),
                             limits = c(0, 1),
                             na.value = "white",
                             name = "Probability of Laterality") +
        guides(fill = guide_colourbar(title.vjust = .95)) + # the legend title needs to move up a bit
        labs(x = "Absolute cumulative yaw", y = "Average pitch") +
        ggtitle("a") +
        theme_classic() +
        theme(text = element_text(size = 11),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 9),
              panel.grid.minor = element_line(color = "white"),
              axis.title.x = element_text(margin = margin(t = 5))) + # increase distance between x-axis values and title
        scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
        scale_y_continuous(expand = c(0, 0))
    )
  }
}





#-----------------------------------------------------------------------------------------------------------------------
## Step 5.1: Does laterality help with better performance when individuals are not experienced? Flight performance #####
#-----------------------------------------------------------------------------------------------------------------------
#flight performance ~ level of handedness * age

#open data
laterality_circling_thin <- readRDS("thinned_laterality_w_gps.rds") #this file doesnt have the attributes... problematic for back-transforming

#I have mean of sd of yaw, pitch, and roll. consider calculating the maxes for this analysis.... will need to go back to the original
#code where I summarized the Quat for each 8-sec burst

#### ----------------------- look at multi-collinearity
laterality_circling_thin %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
                  "abs_cum_yaw", "days_since_tagging")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #sd of roll and pitch = 0.65; sd of yaw and cumulative yaw = 0.7

#### ----------------------- exploratory plot
ggplot(laterality_circling_thin, aes(x = factor(laterality_bi), y = mean_pitch_sd)) +
  geom_boxplot() +
  labs(x = "Laterality", y = "pitch wobble") +
  theme_minimal()

ggplot(laterality_circling_thin, aes(x = factor(laterality_bi), y = mean_roll_sd)) +
  geom_boxplot() +
  labs(x = "Laterality", y = "roll wobble") +
  theme_minimal()

ggplot(laterality_circling_thin, aes(x = factor(laterality_bi), y = mean_yaw_sd)) +
  geom_boxplot() +
  labs(x = "Laterality", y = "yaw wobble") +
  theme_minimal()

#model each flight performance separately, but because pitch and roll are correlated, maybe only do pitch and yaw? ...
#also because laterality is calculated using roll anyway.... maybe also control for thermal strength


#### ----------------------- model: regression with inla

#first bin the age variable
data <- laterality_circling_thin %>%  
  mutate(age_group = inla.group(days_since_tagging_z, n = 200, method = "quantile"),
         laterality_bi = as.factor(laterality_bi))

set.seed(777)
#then add new data to predict the interaction between age and laterality
#create a new dataset to use for making predictions for interaction of mean pitch and cumulative yaw
new_data <- expand.grid(x = unique(data$age_group), #range of age_group
                        y = c(0, 1)) %>% #laterality values
  rename(age_group = x,
         laterality_bi = y) %>%
  mutate(laterality_bi = as.factor(laterality_bi),
         #set the dependent variables to NA
         mean_yaw_sd_z = NA, 
         mean_pitch_sd_z = NA, 
         #randomly assign individual IDs
         individual_local_identifier = sample(laterality_circling_thin$individual_local_identifier, nrow(.), replace = T) %>% 
           as.factor())


#append the new datasets to original data
data <- data %>% 
  mutate(individual_local_identifier = as.factor(individual_local_identifier)) %>% 
  dplyr::select(names(new_data)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data)


####
m_inla <- inla(mean_yaw_sd_z ~ 1 + laterality_bi * age_group +
                 f(individual_local_identifier, age_group, model = "iid"),
               data = data,
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

m_inla <- inla(mean_pitch_sd_z ~ 1 + laterality_bi * age_group +
                 f(individual_local_identifier, age_group, model = "iid"),
               data = data,
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted


####  model evaluation -----------------------

#look at residuals
residuals <- data$mean_yaw_sd - m_inla$summary.fitted.values$mean
plot(residuals)

#calculate variance explained: McFadden's pseudo R²
log_likelihood_full <- m_inla$mlik[1, 1]

null_model <- inla(mean_yaw_sd ~ 1,
                   control.compute = list(cpo = TRUE),
                   control.predictor = list(link = 1, compute = TRUE),
                   data = data)

log_likelihood_null <- null_model$mlik[1, 1]

pseudo_r_squared <- 1 - (log_likelihood_full / log_likelihood_null) #0.002

#model validation metrics
eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), # 0.27 #improves when predicting the z-transformed variable
                   Mlik = as.numeric(m_inla$mlik[,1][2])) # -86459.69


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
X11(width = 3.42, height = 2.3)
(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#8a2be2ff", size = 2)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Intercept", "Laterality", "Days since tagging", "Laterality:Days since tagging"))) + 
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#8a2be2ff", linewidth = 0.5) +
    theme_classic() +
    ggtitle("a") +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5))) #increase distance between x-axis values and title
)

# ggsave(plot = coefs, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_coeffs.pdf", 
#        device = "pdf", width = 3.42, height = 1.5, dpi = 600)
# 
# ggsave(plot = coefs, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_pitch_model_coeffs.pdf", 
#        device = "pdf", width = 3.42, height = 1.5, dpi = 600)


## conclusions: horizontal wobble (sd_yaw) goes down with age and with laterality,
## vertical wobble (sd_pitch) goes up with age and laterality, BUT there is a negative interaction term
## TO DO: 1) plot the interaction term. 2) figure out if the response should be transformed or not. 3) add ind level
## random slope on age and plot ind-specific differences (but this isn't really the focus here... so maybe just add
## inds as a random effect on the slope, just to control for indivdiual variation. )


#### ind_specific coefficients plot -----------------------------------------------------------------------------

#extract handedness for each individual. to use for coloring
handedness <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir)

#extract random effects,  ID is for the unique individuals
#age
random_effects_age <- m_inla$summary.random$individual_local_identifier

#extract unique individual IDs from original data
ind_IDs <- unique(data$individual_local_identifier)

age <- random_effects_age %>% 
  mutate(coef = mean + graph %>% filter(Factor == "age_group") %>% pull(Estimate),
         lower = .[,4] +  graph %>% filter(Factor == "age_group") %>% pull(Lower),
         upper = .[,6] +  graph %>% filter(Factor == "age_group") %>% pull(Upper),
         variable = "age_group") %>% 
  full_join(handedness, by = c("ID" = "individual_local_identifier")) %>% 
  mutate(v_line =  graph %>% filter(Factor == "age_group") %>% pull(Estimate))

#re-order the individuals based on the laterality index
age$ID <- reorder(age$ID, desc(age$laterality_dir))

X11(width = 7, height = 9)
(coefs_inds <- ggplot(age, aes(x = coef, y = ID, color = laterality_dir)) +
    geom_vline(data = age, 
               aes(xintercept = graph %>% filter(Factor == "age_group") %>% pull(Estimate)), linetype = "dashed", color = "gray75",size = 0.5) + 
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    scale_color_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#20A387FF", "right_handed" = "#B8DE29FF"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed"),
                       name = Handedness) +
    #scale_color_manual(values = c("left_handed" = "#33638DFF", "ambidextrous" = "#20A387FF", "right_handed" = "#B8DE29FF"),
    #                   labels = c("Left-handed", "Ambidextrous", "Right-handed")) +
    #scale_color_manual(values = c("left_handed" = "royalblue", "ambidextrous" = "yellowgreen", "right_handed" = "lightcoral"),
    #                   name = "Handedness") +
    scale_y_discrete(labels = ind_IDs) +
    labs(x = "Estimate", y = "") +
    theme_classic() +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))  #increase distance between x-axis values and title
)

ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_coeffs_inds.pdf", 
       device = "pdf", width = 7, height = 9, dpi = 600)

ggsave(plot = coefs_inds, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_pitch_model_coeffs_inds.pdf", 
       device = "pdf", width = 7, height = 9, dpi = 600)


#### interaction plot -----------------------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(is.na(data$mean_yaw_sd_z))

preds <- data.frame(laterality_bi = data[na_rows,"laterality_bi"],
                    age_group = data[na_rows,"age_group"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"],
                    lower95 = m_inla$summary.fitted.values[na_rows,"0.025quant"],
                    upper95 = m_inla$summary.fitted.values[na_rows,"0.975quant"]) %>% 
  #back transform the age values
  mutate(age = age_group * sd(laterality_circling_thin$days_since_tagging) + mean(laterality_circling_thin$days_since_tagging),
         preds_yaw = preds * sd(laterality_circling_thin$mean_yaw_sd) + mean(laterality_circling_thin$mean_yaw_sd))
#preds_pitch = preds * sd(laterality_circling_thin$mean_pitch_sd) + mean(laterality_circling_thin$mean_pitch_sd),
#upper = upper95 * sd(laterality_circling_thin$mean_pitch_sd) + mean(laterality_circling_thin$mean_pitch_sd),
#lower = lower95 * sd(laterality_circling_thin$mean_pitch_sd) + mean(laterality_circling_thin$mean_pitch_sd))


#plot
(pred_p <- preds %>% 
    ggplot() +
    geom_point(aes(x = age, y = preds_yaw, color = laterality_bi), size = 0.2) +
    #geom_line(aes(x = age, y = preds_pitch, color = laterality_bi), linewidth = 0.5) +
    #geom_ribbon(aes(x = age, y = preds_pitch, fill = laterality_bi, ymin = lower, ymax = upper), alpha = 0.75) +
    geom_smooth(aes(x = age, y = preds_yaw, color = laterality_bi, fill = laterality_bi), method = "lm", se = T, alpha = 0.3, linewidth = 0.3) +
    scale_color_manual(values = c("0" = "gray40", "1" = "#8a2be2ff"),
                       name = "Laterality")+
    scale_fill_manual(values = c("0" = "gray40", "1" = "#8a2be2ff"),
                      guide = F)+
    labs(x = "Days since tagging", y = "Horizontal wobble") +
    #labs(x = "Age", y = "Vertical wobble") +
    ggtitle("b") +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(fill = c("0" = "gray40", "1" = "#8a2be2ff"),
                                                     alpha = 0.5))) +
    #xlim(1, 300) +
    theme(text = element_text(size = 11),
          legend.text = element_text(size = 10),
          legend.title = element_text(size = 9),
          #axis.text = element_text(color = "gray45"),
          #panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))
)


#combine the two plots for the coefficients and the interaction term
X11(width = 7, height = 2.3)
model_output_p <- grid.arrange(coefs, pred_p, nrow = 1, widths = c(0.4, 0.6))

#ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_output.pdf", 
#       device = "pdf", width = 7, height = 2.3, dpi = 600)

ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_pitch_model_output.pdf", 
       device = "pdf", width = 7, height = 2.3, dpi = 600)


#-----------------------------------------------------------------------------------------------------------------------
## Step 5.2: Does laterality help with better performance when individuals are not experienced? migration performance #####
#-----------------------------------------------------------------------------------------------------------------------
#migration performance ~ level of handedness  (one row per individual)

#extract handedness for each individual. to use for coloring
handedness <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir) %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier))

#data from Ellen
migr_dates <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/HB_dates_Fin_all_M2.rds")
migr_info <- readRDS("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/from_Ellen/mig_sum.rds") %>% 
  full_join(handedness) %>% 
  drop_na(individual_local_identifier) %>% 
  mutate(laterality_bi = ifelse(laterality_dir == "ambidextrous", 0, 1) %>% as.factor()) #create a binary variable for handedness vs not)


#exploratory plots
ggplot(migr_info, aes(x = laterality_dir, y = avg_speed_kmh)) +
  geom_boxplot() +
  geom_point() +
  labs(x = "Laterality", y = "avg_speed_kmh") +
  theme_classic()

ggplot(migr_info, aes(x = laterality_bi, y = avg_speed_kmh)) +
  geom_boxplot() +
  geom_point() +
  labs(x = "Laterality", y = "avg_speed_kmh") +
  theme_classic()


