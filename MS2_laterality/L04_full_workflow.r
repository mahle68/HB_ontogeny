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
#use pitch

#-------------------------------------------------------
## Step 3.1: Is there laterality? Population-level #####
#-------------------------------------------------------

# #data: angle summaries for each second (from 01b_imu_processing.r)
# sec <- readRDS("quat_summaries_1sec_Jul24.rds")
# 
# #plot the distribution of roll angles in polar histograms, one for each category of circling
# 
# ggplot(data = sec, aes(x = roll_mean)) +
#   geom_histogram(binwidth = 2, color = "blue", fill = "transparent") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
#   #coord_polar(start = -pi) +
#   #scale_x_continuous(breaks = seq(0, 330, 10)) +
#   theme_minimal()


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
    theme_linedraw() +
    theme(text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8.5),
          axis.text = element_text(color = "gray45"),
          panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))) #increase distance between x-axis values and title

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

X11(width = 7, height = 9)
(p_inds <- ggplot(data = circling_w_LI, aes(x = mean_roll_mean, y = individual_local_identifier, 
                                            color = laterality_dir, fill = laterality_dir),
                  height = stat(density)) +
    geom_density_ridges(stat = "binline", bins = 100, scale = 0.98, alpha = 0.8, draw_baseline = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
    scale_x_continuous(breaks = seq(-90, 90, by = 10), 
                       labels = seq(-90, 90, by = 10),
                       limits = c(-92, 92)) +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0.03))) +  # Add extra space above the highest level
    scale_fill_manual(values = c("left_handed" = "royalblue" , "ambidextrous" = "yellowgreen", "right_handed" = "lightcoral")) +
    scale_color_manual(values = c("left_handed" = "royalblue", "ambidextrous" = "yellowgreen", "right_handed" = "lightcoral")) +
    labs(x = "Bank angle (deg) averaged over each 8 second burst",
         y = "Individual",
         fill = "Handedness") +
    theme_linedraw() +
    theme(text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8.5),
          axis.text = element_text(color = "gray45"),
          panel.grid.major = element_line(color = "gray75"),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 5)))  +
    guides(color = "none"))

ggsave(plot = p_inds, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/roll_distr_circling_inds.pdf", 
       device = "pdf", width = 7, height = 9, dpi = 600)

#---------------------------------------------------------------------
## Step 4: Is laterality more likely when the task is difficult? #####
#---------------------------------------------------------------------
#logistic regression: handedness ~ tightness of circles * age . maybe only for individuals with laterality

#extract the IDs of individuals with handedness
handed_IDs <- readRDS("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de2/Work/Projects/HB_ontogeny_eobs/R_files/circling_w_LI_population.rds") %>% 
  filter(laterality_dir != "ambidextrous") %>% 
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  distinct(individual_local_identifier)


#measure of difficulty: whether they circle tighter or wider than the population average for that age
#i.e., model abs_cumulative_yaw ~ age. take the residuals as the proxy for difficulty.

#use the 8-sec data and calculate handedness for each 8 second burst
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
         # abs_cum_yaw_scale = scale(abs(cumulative_yaw_8sec)),
         # mean_pitch_scale = scale(mean_pitch_mean),
         # age_scale = scale(days_since_tagging)) %>%  #REPLACE WITH REAL AGE LTR
  as.data.frame()

#### look at multcollinearity
laterality_circling_thin %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_8sec", 
                  "abs_cum_yaw", "days_since_tagging")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #sd of roll and pitch = 0.65; sd of yaw and cumulative yaw = 0.7

# ###limit to the first 100 days
# laterality_circling_100 <- readRDS("laterality_index_per_8sec_burst_days_since.rds") %>% 
#   filter(n_records >= 8 & #remove short bursts
#            individual_local_identifier %in% handed_IDs$individual_local_identifier) %>% #only keep individuals with handedness
#   mutate(days_since_tagging = as.numeric(days_since_tagging),
#          individual_local_identifier = as.factor(individual_local_identifier),
#          individual_local_identifier2 = individual_local_identifier) %>% #dublicate individual ID to use for inla random effects specification
#   mutate(circling_status = case_when(
#     between(cumulative_yaw_8sec, -10, 10) ~ "straight",
#     cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45 ~ "circling",
#     .default = "shallow circling"
#   )) %>% 
#   mutate(laterality_dir = case_when(
#     between(laterality_bank, 0.25, 1.0) ~ "right_handed",
#     between(laterality_bank, -1.0, -0.25) ~ "left_handed",
#     between(laterality_bank, -0.25, 0.25) ~ "ambidextrous",
#     .default = NA)) %>% 
#   mutate(laterality_bi = ifelse(laterality_dir == "ambidextrous", 0, 1), #create a binary variable for handedness vs not
#          abs_cum_yaw = abs(cumulative_yaw_8sec),
#          abs_cum_yaw_scale = scale(abs(cumulative_yaw_8sec)),
#          mean_pitch_scale = scale(mean_pitch_mean),
#          age_scale = scale(days_since_tagging)) %>%  #REPLACE WITH REAL AGE LTR
#   filter(circling_status == "circling" & #only keep circling flight
#            days_since_tagging <= 100) %>% #3 months
#   as.data.frame()
# 
# #thin the data, so that the consecutive bursts are at least 18 seconds apart: 10 sec between the end of one and start of next
# circling_thinned_100 <- laterality_circling_100 %>% 
#   group_by(individual_local_identifier) %>% 
#   arrange(start_timestamp, .by_group = TRUE) %>% 
#   #thin the data, so that the consecutive bursts are at least 30 seconds apart
#   mutate(time_lag_sec = if_else(row_number() == 1, 0, 
#                                 difftime(start_timestamp, lag(start_timestamp), units = "secs") %>% as.numeric())) %>% 
#   filter(time_lag_sec >= 18) %>% 
#   as.data.frame()


#exploratory plot
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

#### GAM: predict abs_cumulative_yaw as a function of age

m_difficulty <- gam(abs_cum_yaw ~ s(days_since_tagging),
                    data = circling_thinned_100, method = "REML")
#gam.check(m_difficulty)

b <- getViz(m_difficulty)
o <- plot( sm(b, 1))
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

#try a model including pitch, to account for the strength of the thermal.
m_difficulty2 <- gam(abs_cum_yaw ~ mean_pitch_scale + s(age_scale, k = 10),
                     data = circling_thinned_100, method = "REML")
gam.check(m_difficulty2)

#plot the smooth term
b2 <- getViz(m_difficulty2)
o <- plot( sm(b2, 1))
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

#extract the residuals and append to the original dataset
circling_thinned_100 <- circling_thinned_100 %>% 
  mutate(difficulty = residuals(m_difficulty2))

#### try a linear model for difficulty
m_difficulty_ln <- lm(abs_cum_yaw ~ age_scale + mean_pitch_scale,
                      data = circling_thinned_100)


#check for temporal autocorrelation
#check_model(m_difficulty_ln)

par(mfrow =c(2,2)); plot(m_difficulty_ln)

pacf(residuals(m_difficulty_ln)) ####explore this later!!!!

plot(residuals(m_difficulty_ln), type = "b")
abline(h = 0,lty = 3)

### add residuals to the original data
circling_thinned <- circling_thinned %>% 
  mutate(difficulty_ln = residuals(m_difficulty_ln))

#### Model probability of laterality as a function of difficulty of the task
m_laterality <- glm(laterality_bi ~ scale(difficulty),
                    data = circling_thinned_100,
                    family = "binomial")

summary(m_laterality)
check_model(m_laterality)


#add ind random effect
m_laterality2 <- glmer(laterality_bi ~ scale(difficulty) + (1 | individual_local_identifier),
                       data = circling_thinned_100,
                       family = "binomial")

summary(m_laterality2)
check_model(m_laterality2)


#try a gam
m_laterality3 <- gam(laterality_bi ~ s(scale(difficulty)) + 
                       s(individual_local_identifier, bs = "re"),
                     data = circling_thinned_100, method = "REML", family = "binomial")

summary(m_laterality3)

#check for temporal autocorrelation
check.gamViz(getViz(m_laterality3))

b <- getViz(m_laterality3)
print(plot(b, allTerms = T), pages = 1) 


o <- plot( sm(b, 1))
o + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
  l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()


###################### effect plot
# Generate a sequence of difficulty values for prediction
difficulty_seq <- seq(min(laterality_circling_227$difficulty), max(laterality_circling_227$difficulty), length.out = 100)

# Create a data frame with the sequence of difficulty values
prediction_data <- data.frame(difficulty = difficulty_seq)

# Predict probabilities using the model
prediction_data$probability <- predict(m_laterality, newdata = prediction_data, type = "response")

# Plot the data
ggplot(prediction_data, aes(x = difficulty, y = probability)) +
  geom_line(color = "blue") +
  labs(title = "Probability of Laterality vs Difficulty Values",
       x = "Difficulty Values",
       y = "Probability of Laterality") +
  theme_minimal()


######################### old stuff with inla
#to make predictions with inla, the new data set can be appended to the original and the missing values will be predicted during the model building
set.seed(770)

#create a new dataset to use for making predictions for interaction of age and cumulative yaw
new_data_age_yaw <- expand.grid(x = seq(1,300, length.out = 50), #range of age
                        y = seq(1,360, length.out = 50)) %>% #range of cumulative yaw
  rename(abs_cum_yaw_z = y,
         days_since_tagging_z = x) %>%
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate (abs_cum_yaw_z = (abs_cum_yaw_z - attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "abs_cum_yaw_z"],'scaled:center'))/
            attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "abs_cum_yaw_z"],'scaled:scale'),
          days_since_tagging_z = (days_since_tagging_z - attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "days_since_tagging_z"],'scaled:center'))/
            attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "days_since_tagging_z"],'scaled:scale'),
          individual_local_identifier = rep(sample(laterality_circling_thin$individual_local_identifier, 10), nrow(.)/10), #randomly select 7 individual IDs
          individual_local_identifier2 = individual_local_identifier,
          individual_local_identifier3 = individual_local_identifier,
          laterality_bi = NA,
          #set other variables to their mean
          mean_pitch_mean_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "mean_pitch_mean_z"],'scaled:center'),
          mean_yaw_sd_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "mean_yaw_sd_z"],'scaled:center'),
          mean_pitch_sd_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "mean_pitch_sd_z"],'scaled:center'),
          cumulative_pitch_8sec_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "cumulative_pitch_8sec_z"],'scaled:center'),
          interaction = "age_yaw")

#create a new dataset to use for making predictions for interaction of mean pitch and cumulative yaw
new_data_pitch_yaw <- expand.grid(x = seq(-50, 50, length.out = 50), #range of pitch
                                y = seq(1, 360, length.out = 50)) %>% #range of cumulative yaw
  rename(abs_cum_yaw_z = y,
         mean_pitch_mean_z = x) %>%
  #calculate z-scores. get the mean(center) and sd(scale) from the previous z transformation. for consistency
  mutate (abs_cum_yaw_z = (abs_cum_yaw_z - attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "abs_cum_yaw_z"],'scaled:center'))/
            attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "abs_cum_yaw_z"],'scaled:scale'),
          mean_pitch_mean_z = (mean_pitch_mean_z - attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "mean_pitch_mean_z"],'scaled:center'))/
            attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "mean_pitch_mean_z"],'scaled:scale'),
          individual_local_identifier = rep(sample(laterality_circling_thin$individual_local_identifier, 10), nrow(.)/10), #randomly select 7 individual IDs
          individual_local_identifier2 = individual_local_identifier,
          individual_local_identifier3 = individual_local_identifier,
          laterality_bi = NA,
          #set other variables to their mean
          days_since_tagging_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "days_since_tagging_z"],'scaled:center'),
          mean_yaw_sd_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "mean_yaw_sd_z"],'scaled:center'),
          mean_pitch_sd_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "mean_pitch_sd_z"],'scaled:center'),
          cumulative_pitch_8sec_z = attr(laterality_circling_thin[,colnames(laterality_circling_thin) == "cumulative_pitch_8sec_z"],'scaled:center'),
          interaction = "pitch_yaw")


#append the new datasets to original data
data <- laterality_circling_thin %>% 
  mutate(interaction = "none") %>% 
  dplyr::select(names(new_data_pitch_yaw)) %>%  #only keep the columns that are necessary for the model
  bind_rows(new_data_pitch_yaw,new_data_age_yaw)





#### model: binomial logistic regression with inla -----------------------------------------------------------
#only age
# m_inla <- inla(laterality_bi ~ 1 + age_scale +  
#                  f(individual_local_identifier2, model = "iid"),  
#                data = data, family = "binomial",
#                control.compute = list(cpo = TRUE),
#                control.predictor = list(link = 1, compute = TRUE) #compute=t means that NA values will be predicted.
# )
# 
# #model validation metrics
# eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), #  0.5188965
#                    Mlik = as.numeric(m_inla$mlik[,1][2])) # -21795.92
# 
# #only age, but smooth
# data$age_scale_binned <- inla.group(data$age_scale, n = 50)
# 
# m_inla <- inla(laterality_bi ~ 1 + f(age_scale_binned, model = "rw1", constr = FALSE) +  
#                  f(individual_local_identifier, model = "iid"),  
#                data = data, family = "binomial",
#                control.compute = list(cpo = TRUE),
#                control.predictor = list(link = 1, compute = TRUE) #compute=t means that NA values will be predicted.
# )
# #*** warning *** iterative process seems to diverge, 'vb.correction' is aborted
# #*** Please (re-)consider your model, priors, confounding, etc.
# 
# (eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), # 0.5204054
#                     Mlik = as.numeric(m_inla$mlik[,1][2])) # -21794.79
# )
# 
# # age and cum yaw
# m_inla <- inla(laterality_bi ~ 1 + abs_cum_yaw_scale * age_scale +
#                  f(individual_local_identifier, abs_cum_yaw_scale, model = "iid") +  
#                  f(individual_local_identifier2, age_scale, model = "iid"),  
#                data = data, family = "binomial",
#                control.compute = list(cpo = TRUE),
#                control.predictor = list(link = 1, compute = TRUE) #compute=t means that NA values will be predicted.
# )
# 
# (eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), # 0.5184729
#                     Mlik = as.numeric(m_inla$mlik[,1][2])) # -21814.95
# )
# 
# ####
# # Bin the predictor variables
# data$abs_cum_yaw_scale_binned <- inla.group(data$abs_cum_yaw_scale, n = 100)
# data$age_scale_binned <- inla.group(data$age_scale, n = 100)
# 
# # Fit the model with binned variables
# m_inla <- inla(laterality_bi ~ 1 + abs_cum_yaw_scale_binned * age_scale_binned +
#                  f(individual_local_identifier, abs_cum_yaw_scale_binned, model = "rw2", constr = FALSE) + 
#                  f(individual_local_identifier2, age_scale_binned, model = "rw2", constr = FALSE),
#                data = data, family = "binomial",
#                control.compute = list(cpo = TRUE),
#                control.predictor = list(link = 1, compute = TRUE) #compute=t means that NA values will be predicted.
# )
# 
# data$interaction_term <- interaction(data$abs_cum_yaw_scale_binned, data$age_scale_binned)
# 
# # Create the interaction term
# data$interaction_term <- data$abs_cum_yaw_scale_binned * data$age_scale_binned
# 
# # Fit the model
# m_inla <- inla(
#   laterality_bi ~ 1 + 
#     f(abs_cum_yaw_scale_binned, model = "rw1") + 
#     f(age_scale_binned, model = "rw1") +
#     f(interaction_term, model = "rw1") + # Use rw1 for the interaction term
#     f(individual_local_identifier, model = "iid"),
#   data = data, 
#   family = "binomial",
#   control.compute = list(cpo = TRUE),
#   control.predictor = list(link = 1, compute = TRUE) # compute=TRUE means that NA values will be predicted.
# )


#model with mean pitch for the strength of the thermal
m_inla <- inla(laterality_bi ~ 1 + abs_cum_yaw_z * days_since_tagging_z + mean_pitch_mean_z +
                 f(individual_local_identifier, abs_cum_yaw_z, model = "iid") +  
                 f(individual_local_identifier2, days_since_tagging_z, model = "iid"),  
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted.

m_inla <- inla(laterality_bi ~ 1 + abs_cum_yaw_z + mean_pitch_mean_z + days_since_tagging_z + 
                 f(individual_local_identifier, abs_cum_yaw_z, model = "iid") +  
                 f(individual_local_identifier2, days_since_tagging_z, model = "iid") ,
                 f(individual_local_identifier3, mean_pitch_mean_z, model = "iid"),
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#age as a smooth term
#first bin the age variable
data <- data %>% 
  mutate(age_group = inla.group(days_since_tagging_z, n = 30, method = "quantile"))

m_inla <- inla(laterality_bi ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z + 
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                 f(age_group, model = "rw1"),
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#m_inla <- inla(laterality_bi ~ 1 + abs_cum_yaw_z + mean_pitch_mean_z +
#                 f(age_scale, model = "rw1") +
#                 f(individual_local_identifier, model = "iid"),  
#               data = data, family = "binomial",
#               control.compute = list(cpo = TRUE),
#               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted.



#### model evaluation -----------------------------------------------------------------------------

# calibration_data <- data.frame(
#   predicted = m_inla$summary.fitted.values$mean,
#   observed = data$laterality_bi
# )
# ggplot(calibration_data, aes(x = predicted, y = observed)) +
#   geom_point() +
#   geom_smooth(method = "loess") +
#   labs(x = "Predicted Probability", y = "Observed Proportion")
# 
# 
# #look at residuals
# residuals <- data$laterality_bi - m_inla$summary.fitted.values$mean
# plot(residuals)
# 
# #calculate variance explained: McFadden's pseudo R²
# log_likelihood_full <- m_inla$mlik[1, 1]
# 
# null_model <- inla(laterality_bi ~ 1, family = "binomial", 
#                    control.compute = list(cpo = TRUE),
#                    control.predictor = list(link = 1, compute = TRUE), 
#                    data = data)
# 
# log_likelihood_null <- null_model$mlik[1, 1]
# 
# pseudo_r_squared <- 1 - (log_likelihood_full / log_likelihood_null) #0.002
# 
# #model validation metrics
# eval <- data.frame(CPO = mean(m_inla_s$cpo$cpo, na.rm = T), #0.52
#                    Mlik = as.numeric(m_inla_s$mlik[,1][2])) # -21814.95

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

graph$Factor_n <- as.numeric(graph$Factor)

#plot the coefficients
#X11(width = 8, height = 6)

(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray", linewidth = 0.5) +
    geom_point(color = "cornflowerblue", size = 2)  +
    labs(x = "Estimate", y = "") +
    #scale_y_discrete(labels = rev(c("Intercept","Absolute cumulative yaw", "Age", "Average pitch", 
    #                                "Absolute cumulative yaw : age"))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "cornflowerblue", linewidth = 1) +
    theme_classic() +
    theme(text = element_text(size = 24)) 
)

#### interaction plot -----------------------------------------------------------------------------
#extract information for rows that had NAs as response variables

laterality_na <- which(data$interaction == "pitch_yaw")

preds <- data.frame(abs_cum_yaw_z = data[data$interaction == "pitch_yaw" ,"abs_cum_yaw_z"],
                    mean_pitch_mean_z = data[data$interaction == "pitch_yaw" ,"mean_pitch_mean_z"],
                    preds = m_inla$summary.fitted.values[laterality_na,"mean"]) %>% 
  mutate(prob_pres = exp(preds)/(1 + exp(preds))) #this should be between 0-1.

#create a raster and run a moving window to make the pattern easier to see.
pred_r <- preds %>% 
  dplyr::select(mean_pitch_mean_z, abs_cum_yaw_z, prob_pres) %>% #figure out if i need to do the transformation
  terra::rast(type = "xyz") %>%
  focal(w = 11, fun = median, na.policy = "all", na.rm = T) %>%
  as.data.frame(xy = T) %>% 
  rename(prob_pres = focal_median)

#plot
(  pred_p <- pred_r %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = prob_pres)) +
    scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                         na.value = "white", name = "Probability of laterality") +
    labs(x = "age", y = "circling intensity") +
    theme_classic()
)



###### old code for only a two-way interaction------------------------------------------------------------------------------------------

laterality_na <- which(is.na(data$laterality_bi))

preds <- data.frame(abs_cum_yaw_z = data[is.na(data$laterality_bi) ,"abs_cum_yaw_z"],
                    days_since_tagging_z = data[is.na(data$laterality_bi) ,"days_since_tagging_z"],
                    preds = m_inla$summary.fitted.values[laterality_na,"mean"]) %>% 
  mutate(prob_pres = exp(preds)/(1 + exp(preds))) #this should be between 0-1.

#create a raster and run a moving window to make the pattern easier to see.
pred_r <- preds %>% 
  dplyr::select(days_since_tagging_z, abs_cum_yaw_z, preds) %>% #figure out if i need to do the transformation
  terra::rast(type = "xyz") %>%
  focal(w = 11, fun = median, na.policy = "all", na.rm = T) %>%
  as.data.frame(xy = T) %>% 
  rename(prob_pres = focal_median)

#plot
(  pred_p <- pred_r %>% 
    ggplot() +
    geom_tile(aes(x = x, y = y, fill = prob_pres)) +
    scale_fill_gradient2(low = "#005AB5", mid = "seashell2", high = "#D41159",limits = c(0,1), midpoint = 0.5,
                         na.value = "white", name = "Probability of laterality") +
    labs(x = "age", y = "circling intensity") +
    theme_classic()
)

### odds ratio: likelihood of getting a presence
#https://cran.r-project.org/web/packages/sjPlot/vignettes/plot_model_estimates.html


################## even older stuff. with glmer

m1 <- glm(laterality_bi ~ abs(cumulative_yaw_8sec) * days_since_tagging,
          data = laterality_circling_227,
          family = "binomial")

# Calculate R-squared
deviance <- summary(m1)$deviance
null_deviance <- summary(m1)$null.deviance
rsquared <- 1 - (deviance / null_deviance)

# Print the R-squared value
print(rsquared) #0.002


# Add logistic fitted values back to dataframe as
# new column pred.g190
laterality_1sec_days$pred.g190 <- m1$fitted.values

head(laterality_1sec_days)

m2 <- glmer(laterality_bi ~ cumulative_pitch_8sec_z * abs_cum_yaw_z * days_since_tagging_z + 
              (1 | individual_local_identifier),
            data = laterality_circling_thin,
            family = "binomial")

summary(m2)

#check for temporal autocorrelation
library(performance)

check_model(m2)

#-----------------------------------------------------------------------------------------------------------------------
## Step 5.1: Does laterality help with better performance when individuals are not experienced? Flight performance #####
#-----------------------------------------------------------------------------------------------------------------------
#flight performance ~ level of handedness * age




#-----------------------------------------------------------------------------------------------------------------------
## Step 5.2: Does laterality help with better performance when individuals are not experienced? migration performance #####
#-----------------------------------------------------------------------------------------------------------------------
#migration performance ~ level of handedness  (one row per individual)


