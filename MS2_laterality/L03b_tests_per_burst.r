#This code runs tests for whether honey buzzards have a preference for left or right bank angles during different life cycle stage
#at the scale of each 8-second burst
#Elham Nourani PhD.
#July 24. 2024. Konstanz, DE. 

library(tidyverse)
library(lubridate)
library(sf)
library(ggridges)
library(viridis)
library(lme4)
library(mgcv)
library(mgcViz)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")


#-------------------------------------------------------------------------------------
# STEP1: annotate data with life cycle stage
#-------------------------------------------------------------------------------------
#open data with laterality index calculated

#this file contains laterality index calculated per burst. NOT filtered for circling flight
laterality_1sec <- readRDS("laterality_index_per_8sec_burst.rds")

#open meta-data to calculate day since tagging
meta_data <- read.csv("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/data/EHB_metadata - Sheet1.csv") %>% 
  mutate(deployment_dt_utc = as.POSIXct(deployment_dt_utc, tz = "UTC"))

laterality_1sec_days <- laterality_1sec %>%
  mutate(individual_local_identifier = as.character(individual_local_identifier)) %>% 
  full_join(meta_data %>% select(ring_ID, deployment_dt_utc), by = c("individual_local_identifier" = "ring_ID")) %>% 
  rowwise() %>% 
  mutate(days_since_tagging = floor(difftime(start_timestamp, deployment_dt_utc, unit = "days"))) %>% 
  ungroup() %>% 
  as.data.frame()

saveRDS(laterality_1sec_days, "laterality_index_per_8sec_burst_days_since.rds")
  
laterality_1sec_days <- readRDS("laterality_index_per_8sec_burst_days_since.rds") %>% 
  filter(n_records >= 8) %>% #remove short bursts
  #filter(cumulative_yaw_8sec >= 45 | cumulative_yaw_8sec <= -45) %>%  #remove straight movement. do i even need this??
  mutate(days_since_tagging = as.numeric(days_since_tagging),
         individual_local_identifier = as.factor(individual_local_identifier)) %>% 
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

#### HEATMAPS for 2D distributions
ggplot() + #these need to be standardized, because there are so much more data during migration
  stat_density_2d(data = laterality_1sec_days,
                  aes(x = days_since_tagging, y = cumulative_roll_8sec, (fill = after_stat(level)), geom = "polygon")) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

data_100 <- laterality_1sec_days %>% filter(days_since_tagging <= 100)

ggplot(laterality_1sec_days %>% filter(days_since_tagging <= 100), aes(x = days_since_tagging, y = cumulative_yaw_8sec) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw


#### scatterplots with smooth terms
ggplot(data = data_100, aes(x = days_since_tagging, y = abs(mean_roll_max))) +
  geom_point() +
  geom_smooth()

ggplot(data = data_100, aes(x = days_since_tagging, y = abs(laterality_bank))) +
  geom_point() +
  geom_smooth(method = "lm")

ggplot(data = data_100, aes(x = abs(cumulative_yaw_8sec), y = abs(laterality_bank))) +
  geom_point() +
  geom_smooth(method = "lm")

#-------------------------------------------------------------------------------------
# STEP3: GAM
#-------------------------------------------------------------------------------------

#the model will predict the strength of laterality as a function of age AND the tightness of the circle. We expect stronger laterality earlier in life and in tougher soaring conditions (narrower thermals)
#no need to have laterality bank between -1 and 1. I'm interested in the strength, so just take the absolute values
gam1 <- gam(laterality_bank ~ s(cumulative_yaw_8sec, days_since_tagging, by = individual_local_identifier), #cumulative yaw indicates difficulty. the tighter the circle, the more difficult the circling attenpt, the more likely to be lateral (hypothesis)
            data = laterality_1sec_days, method = "REML")
#try with absolute values of laterality and yaw
#test for termporal autocorrelation

saveRDS(gam1, file = "gam_full.rds")

#summary(gam1)

#plot the smooth interaction
#plot(gam1, page = 1, scheme = 2)

rm(gam1)
gc()
gc()

ctrl <- list(nthreads = 4) #this took a very long time to run. over 24 hrs. increase the n of cores next time

gam_abs <- gam(abs(laterality_bank) ~ s(abs(cumulative_yaw_8sec), days_since_tagging, by = individual_local_identifier), #cumulative yaw indicates difficulty. the tighter the circle, the more difficult the circling attenpt, the more likely to be lateral (hypothesis)
            data = laterality_1sec_days, method = "REML", control=ctrl)

summary(gam_abs) #R-sq.(adj) =  0.096   Deviance explained = 9.83%
#https://towardsdatascience.com/producing-insights-with-generalized-additive-models-gams-cf2b68b1b847

saveRDS(gam_abs, file = "gam_abs.rds")

#plot the smooth interaction
plot(gam_abs, page = 1, scheme = 2)


#use mgcViz for plotting: https://mfasiolo.github.io/mgcViz/articles/mgcviz.html
#create a mgcViz object
b <- getViz(gam_abs)
print(plot(b, allTerms = T), pages = 1) 
plot(sm(b, 1)) + l_fitRaster() + l_fitContour() + l_points()

### first 100 days

ctrl <- list(nthreads = 6)
(st_time <- Sys.time())
gam_100 <- gam(abs(laterality_bank) ~ s(days_since_tagging, abs(cumulative_yaw_8sec)), 
                      data = laterality_1sec_days %>% filter(days_since_tagging <= 100), method = "REML", family = betar(link = "probit"), control = ctrl)
Sys.time() - st_time #1.5 minutes

summary(gam_100) #R-sq.(adj) =  0.0316   Deviance explained = -1.08% ##worse than a null model!!!
plot(gam_100)
plot(gam_100, page = 1, scheme = 2, type = "response")



### try simpler GAMS 

ctrl <- list(nthreads = 6)
(st_time <- Sys.time())
gam_population <- gam(abs(laterality_bank) ~ te(days_since_tagging, abs(cumulative_yaw_8sec)), 
               data = laterality_1sec_days, method = "REML", family = betar(link = "probit"), control = ctrl)
Sys.time() - st_time #2.2 minutes
#model specs: The beta family (betar) with a logit link function ensures that the predicted values will always fall between 0 and 1
#"te" is the tensor product smooth, useful for interactions between two continuous variabels


summary(gam_population) #R-sq.(adj) =  0.0584   Deviance explained = 5.85%
plot(gam_population)
plot(gam_population, page = 1, scheme = 2, type = "response")

## gam with only days
gam_days <- gam(abs(laterality_bank) ~ s(days_since_tagging), 
                      data = laterality_1sec_days, method = "REML", family = betar(link = "probit"), control = ctrl)

plot(gam_days, page = 1, scheme = 2)

#gam with only tightness of circles
gam_yaw <- gam(abs(laterality_bank) ~ s(abs(cumulative_yaw_8sec)), 
                data = laterality_1sec_days, method = "REML", family = betar(link = "probit"), control = ctrl)

plot(gam_yaw, page = 1, scheme = 2)

#-------------------------------------------------------------------------------------
# STEP3B: binomial GAM
#-------------------------------------------------------------------------------------
#to ease model convergence, convert laterality index into a binomial response with 1 = laterality (LI = 1.0 to 0.25 OR -1 to -0.25 ), 2 = Ambidexterity (LI = -0.25 - 0.25). thresholds are taken from Bennison et al.

laterality_1sec_days <- readRDS("laterality_index_per_8sec_burst_days_since.rds") %>% 
  filter(n_records >= 8) %>% #remove short bursts
  mutate(days_since_tagging = as.numeric(days_since_tagging),
         individual_local_identifier = as.factor(individual_local_identifier),
         laterality_bank_dir = case_when(
           between(laterality_bank, 0.25, 1) ~ "right",
           between(laterality_bank, -1, -0.25) ~ "left",
           between(laterality_bank, -0.25, 0.25) ~ "ambidextrous"
         ),
         laterality_bank_bin = ifelse(laterality_bank_dir == "ambidextrous", 0, 1)) %>% #create a binomial variable to use as the response in the GAM
  as.data.frame()

#########run the binomial GAM ####
#https://noamross.github.io/gams-in-r-course/chapter4

#model without individual random effects
ctrl <- list(nthreads = 6)
(st_time <- Sys.time())
gam_population <- gam(laterality_bank_bin ~ te(days_since_tagging, abs(cumulative_yaw_8sec)), 
                      data = laterality_1sec_days, method = "REML", family = "binomial", control = ctrl)
Sys.time() - st_time #13 seconds

summary(gam_population) #using s for the smooth: R-sq.(adj) =  0.0251   Deviance explained = 2.05%
#using tensor smooth (useful when the two variables have different scales): R-sq.(adj) =  0.0174   Deviance explained = 1.42%

#convert the intercept from long-odds scale to a probability using the logistic function
plogis(gam_population$coefficients[1]) #0.6965779. interpretation: the model predicts a 67 percent baseline chance of a positive outcome (i.e. being lateral) if all predictors were at their average values. 

plot(gam_population, pages = 1, trans = plogis)
#contour plot with colors
plot(gam_population, scheme = 2, trans = plogis)

#3d plot
plot(gam_population, scheme = 1, trans = plogis)

#another plot
vis.gam(x = gam_population,                # GAM object
        view = c("days_since_tagging", "cumulative_yaw_8sec"),   # variables
        plot.type = "persp")    #3d plot

############## gam with tensor smooth, using ti: separate smooths for each variable in addition to the interaction. needs more data ####
ctrl <- list(nthreads = 6)
(st_time <- Sys.time())
gam_population2 <- gam(laterality_bank_bin ~ s(days_since_tagging) +
                         s(abs(cumulative_yaw_8sec)) +
                         ti(days_since_tagging, abs(cumulative_yaw_8sec)), 
                      data = laterality_1sec_days, method = "REML", family = "binomial", control = ctrl)
Sys.time() - st_time #42 seconds

summary(gam_population2) #-sq.(adj) =  0.0237   Deviance explained = 1.95%

plot(gam_population2, page = 1, trans = plogis, rug = T)
plot(gam_population2, scheme = 2, trans = plogis)

#3d plot
plot(gam_population2, scheme = 1, trans = plogis)


############## gam with tensor smooth: only first 100 days ####

data_100 <- laterality_1sec_days %>% 
  filter(days_since_tagging <= 100) %>% 
  as.data.frame()
  
ctrl <- list(nthreads = 10)
(st_time <- Sys.time())
gam_pop3 <- gam(laterality_bank_bin ~ s(days_since_tagging) +
                         s(abs(cumulative_yaw_8sec)) +
                         te(days_since_tagging, abs(cumulative_yaw_8sec)), 
                       data = data_100, method = "REML", family = "binomial", control = ctrl)
Sys.time() - st_time #24 seconds

summary(gam_pop3) #R-sq.(adj) =  0.0182   Deviance explained = 1.46%

plot(gam_pop3, page = 1, trans = plogis, rug = T, shift = coef(gam_pop3)[1], seWithMean = TRUE)
plot(gam_pop3, scheme = 2, trans = plogis, shift = coef(gam_pop3)[1], seWithMean = TRUE)

plot(gam_pop3, scheme = 2, trans = plogis, shift = coef(gam_pop3)[1], seWithMean = TRUE) #To get appropriately-scaled marginal effects plots, we can use the following specification: plot(binom_mod, pages = 1, trans = plogis, shift = coef(binom_mod)[1], seWithMean = TRUE)

b <- getViz(gam_pop3)
print(plot(b, allTerms = T, trans = plogis), pages = 1) 
plot(sm(b, 1)) + l_fitRaster() + l_fitContour() + l_points()

############## gam with tensor smooth: only first 100 days + individual effect ####

ctrl <- list(nthreads = 10)
(st_time <- Sys.time())
gam_ind <- gam(laterality_bank_bin ~ s(days_since_tagging) +
                  s(abs(cumulative_yaw_8sec)) +
                  ti(days_since_tagging, abs(cumulative_yaw_8sec)) +
                 s(individual_local_identifier, bs = "re"),
                data = data_100, method = "REML", family = "binomial", control = ctrl)
Sys.time() - st_time #56 seconds

summary(gam_ind) #R-sq.(adj) =  0.0216   Deviance explained = 1.74%

plot(gam_ind, page = 1, trans = plogis, rug = T)
plot(gam_ind, scheme = 2, trans = plogis)
plot(gam_ind, scheme = 1, trans = plogis)

vis.gam(x = gam_ind,                # GAM object
        view = c("days_since_tagging", "cumulative_yaw_8sec"),   # variables
        plot.type = "persp")    #3d plot

b <- getViz(gam_ind)
print(plot(b, allTerms = T, trans = plogis), pages = 1) 
plot(sm(b, 1)) + l_fitRaster() + l_fitContour() + l_points()

library(DHARMa) #https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html
simulationOutput <- simulateResiduals(fittedModel = gam_pop3, plot = F)
plot(simulationOutput)
plotResiduals(simulationOutput)
testDispersion(simulationOutput)


plot(gam_pop3, shade = TRUE, scale = 0, residuals = TRUE)

############## gam with circle radius as a categorical variable, time and ind ####

data_circle_r <- laterality_1sec_days %>% 
  mutate(turn_r = case_when(
    between(abs(cumulative_yaw_8sec), 0, 15) ~ "linear",
    between(abs(cumulative_yaw_8sec), 15, 90) ~ "shallow",
    abs(cumulative_yaw_8sec) >= 90 ~ "full",
  ) %>% as.factor())

#some exploration
ggplot(data_circle_r, aes(x = turn_r, y = abs(cumulative_roll_8sec))) +
  geom_boxplot()


ctrl <- list(nthreads = 10)
(st_time <- Sys.time())
gam_ind <- gam(laterality_bank_bin ~ s(days_since_tagging, by = turn_r) +
               s(individual_local_identifier, bs = "re"),
               data = data_circle_r, method = "REML", family = "binomial", control = ctrl)
Sys.time() - st_time #1.3 min

summary(gam_ind) #R-sq.(adj) =  0.0247   Deviance explained = 2.05%

############## gam with ONLY turning events ####
 
turns_only <- data_circle_r %>% 
  filter(turn_r != "linear") #PROVE the threshold makes sense using the gps data

ctrl <- list(nthreads = 10)
(st_time <- Sys.time())
gam_full <- gam(laterality_bank_bin ~ s(days_since_tagging) +
                 s(abs(cumulative_yaw_8sec)) +
                 ti(days_since_tagging, abs(cumulative_yaw_8sec)) +
                 s(individual_local_identifier, bs = "re"),
               data = turns_only, method = "REML", family = "binomial", control = ctrl)
Sys.time() - st_time #1 min

summary(gam_full)# R-sq.(adj) =  0.0407   Deviance explained = 3.54%

#use mgcViz
b <- getViz(gam_full)
print(plot(b, allTerms = T, trans = plogis, shift = coef(gam_full)[1], seWithMean = TRUE), pages = 1) 

############## gam with ONLY turning events AND ONLY 50% quantile of the data: i.e. until day 58 ####

turns_only <- data_circle_r %>% 
  filter(turn_r != "linear") %>%  #PROVE the threshold makes sense using the gps data
  filter(days_since_tagging <= quantile(days_since_tagging, probs = 0.5))

ctrl <- list(nthreads = 10)
(st_time <- Sys.time())
gam_58days <- gam(laterality_bank_bin ~ s(days_since_tagging) +
                  s(abs(cumulative_yaw_8sec)) +
                  ti(days_since_tagging, abs(cumulative_yaw_8sec)) +
                  s(individual_local_identifier, bs = "re"),
                data = turns_only, method = "REML", family = "binomial", control = ctrl)
Sys.time() - st_time #26 seconds

summary(gam_58days) #R-sq.(adj) =  0.0298   Deviance explained = 2.43%

#use mgcViz
b <- getViz(gam_58days)
print(plot(b, allTerms = T, trans = plogis, shift = coef(gam_58days)[1], seWithMean = TRUE), pages = 1) 
