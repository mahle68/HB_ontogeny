#update with BA: use the laterality index as a continuous variable to model the strength of lateralization, rather than the binary lateralized vs not.
#this is because with the calculation of lateralization in bank angle, i don't have many ambidextrous assignments within the 8-sec bursts

#---------------------------------------------------------------------
## Step 4: Is laterality more likely when the task is difficult? #####
#---------------------------------------------------------------------
#logistic regression: handedness ~ tightness of circles * age . maybe only for individuals with laterality
#read in filtered data. this is not filtered for days since tagging
filtered_w_LI <- readRDS("thinned_laterality_w_gps_wind_all_filters_BA.rds")


#### ----------------------- filter for day since tagging and z-transform
quantile(filtered_w_LI$days_since_tagging, probs = 0.9) #284 ... 

circling_data <- filtered_w_LI %>% 
  #filter for days since tagging. only keep the 
  filter(days_since_tagging < 284) %>% 
  #make sure to do all the filtering before scaling the variables!!!!
  mutate_at(c("mean_roll_sd", "bank_angle_deg_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_sum_8sec", 
              "abs_cum_yaw", "wind_speed", "days_since_tagging", "location_lat_closest_gps_raw"),
            list(z = ~scale(.))) %>%
  mutate(laterality_bi_8sec = ifelse(laterality_dir_8sec == "ambidextrous", 0, 1)) %>%  #create a binary variable for handedness vs not)
  as.data.frame()

saveRDS(circling_data, file = "thinned_laterality_for_modeling_BA.rds")

#### ----------------------- look at multi-collinearity

circling_data %>% 
  dplyr::select(c("mean_roll_sd", "bank_angle_deg_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_sum_8sec", 
                  "abs_cum_yaw", "wind_speed", "days_since_tagging", "location_lat_closest_gps_raw")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #sd of roll and pitch = 0.64; sd of roll and sd of bank = 0.56; sd of yaw and cumulative yaw = 0.68. age and latitude = -0.75

#### ----------------------- exploratory plot
ggplot(circling_data, aes(x = factor(laterality_bi_8sec), y = abs(cumulative_yaw_sum_8sec))) +
  geom_boxplot() +
  labs(x = "Laterality", y = "Absolute Cumulative Yaw (8 sec)") +
  theme_minimal()

ggplot(circling_data, aes(x = factor(laterality_dir_8sec), y = days_since_tagging)) +
  geom_boxplot() +
  labs(x = "Laterality", y = "") +
  theme_minimal()


#compare n and s hemisphere
ggplot(circling_data, aes(x = factor(laterality_dir), y = location_lat_closest_gps_raw)) +
  geom_boxplot() +
  labs(x = "Laterality", y = "latitude") +
  theme_minimal()

ggplot(circling_data, aes(x = days_since_tagging, y = abs(cumulative_yaw_sum_8sec))) +
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(circling_data, aes(x = days_since_tagging, y = location_lat_closest_gps_raw)) +
  geom_point() + 
  geom_smooth(method = "lm")

ggplot(circling_data, aes(x = days_since_tagging, y = laterality_bank)) +
  geom_point() + 
  geom_smooth(method = "lm")


#### ----------------------- model: binomial logistic regression with inla

circling_data <- readRDS("thinned_laterality_for_modeling_BA.rds")

#create new data: make predictions for values that fall on a regular grid for optimal visualization 

#to make sure the predictions cover the parameter space, create a dataset with all possible combinations. The max of the variables might be outliers, so use the 90% quantile instead
grd_pitch_yaw <- expand.grid(x = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = max(circling_data$mean_pitch_mean, na.rm = T),  length.out = 50),
                             y = seq(from = min(circling_data$abs_cum_yaw, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(mean_pitch_mean = x,
         abs_cum_yaw = y)  %>% 
  mutate(wind_speed = attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "pitch_yaw")

grd_wind_yaw <- expand.grid(x = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                            y = seq(from = min(circling_data$abs_cum_yaw, na.rm = T), to = quantile(circling_data$abs_cum_yaw, .99, na.rm = T),  length.out = 50)) %>% 
  rename(wind_speed = x,
         abs_cum_yaw = y)  %>% 
  mutate(mean_pitch_mean = attr(circling_data[,colnames(circling_data) == "mean_pitch_mean_z"],'scaled:center'), #set other variables to their mean
         days_since_tagging = attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'),
         interaction = "wind_yaw")

grd_wind_pitch <- expand.grid(x = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                              y = seq(from = quantile(circling_data$mean_pitch_mean, .01, na.rm = T), to = max(circling_data$mean_pitch_mean, na.rm = T),  length.out = 50)) %>% 
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
         laterality_bi_8sec = NA) 


#bin the age variable and append the new data
data <- circling_data %>% 
  dplyr::select(c(intersect(colnames(grd_all), colnames(.)), "life_stage", "location_lat_closest_gps_raw_z")) %>% 
  full_join(grd_all) %>% 
  mutate(individual_local_identifier2 = individual_local_identifier, #repeat individual ID column to be used in the model formula for random effects
         individual_local_identifier3 = individual_local_identifier,
         age_group = inla.group(days_since_tagging_z, n = 100, method = "quantile"), #age will be included as a smooth term
         lat_group = inla.group(location_lat_closest_gps_raw_z, n = 10, method = "quantile"))  


#re-order life stage, so that post-fledging is first
data$life_stage <- factor(data$life_stage, levels = c("post-fledging", "migration", "wintering"))


#saveRDS(data, file = "thinned_laterality_for_modeling_w_new_data_BA.rds")

data <- readRDS("thinned_laterality_for_modeling_w_new_data_BA.rds")

m_inla <- inla(laterality_bi_8sec ~ 1 + mean_pitch_mean_z * abs_cum_yaw_z * wind_speed_z + 
                 f(individual_local_identifier, mean_pitch_mean_z, model = "iid") +  
                 f(individual_local_identifier2, abs_cum_yaw_z, model = "iid") + 
                 f(individual_local_identifier3, wind_speed_z, model = "iid") + 
                 f(age_group, model = "rw1"),
               data = data, family = "binomial",
               control.compute = list(cpo = TRUE),
               control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#### model evaluation -----------------------

#look at residuals
residuals <- data$laterality_bi_8sec - m_inla$summary.fitted.values$mean
plot(residuals)

#calculate variance explained: McFadden's pseudo R²
log_likelihood_full <- m_inla$mlik[1, 1]

null_model <- inla(laterality_bi_8sec ~ 1, family = "binomial",
                   control.compute = list(cpo = TRUE),
                   control.predictor = list(link = 1, compute = TRUE),
                   data = data)

log_likelihood_null <- null_model$mlik[1, 1]

pseudo_r_squared <- 1 - (log_likelihood_full / log_likelihood_null) #0.006

#model validation metrics
eval <- data.frame(CPO = mean(m_inla$cpo$cpo, na.rm = T), # 0.52
                   Mlik = as.numeric(m_inla$mlik[,1][2])) # -13050.66

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
X11(width = 3.2, height = 2)
(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#0d0887", size = 1.5)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Intercept", "Average pitch", "Absolute cumulative yaw",
                                    "Wind speed", "Average pitch: Absolute cumulative yaw", "Average pitch: Wind speed",
                                    "Average cumulative yaw: Wind speed", "Average pitch : Average cumulative yaw: \nWind speed "))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#0d0887", linewidth = 0.5) +
    ggtitle("a") +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"), #top, right, bottom, left
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          plot.title = element_text(face = "bold"), # make title bold
          axis.title.x = element_text(margin = margin(t = 2))) #increase distance between x-axis values and title
)

ggsave(plot = coefs, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_coeffs_2min.pdf", 
       device = "pdf", width = 7, height = 4, dpi = 600)


#### plot the smooth term -----------------------------------------------------------------------------

# average start and end of migration... to add to the plot
population_migr_period <- readRDS("thinned_laterality_for_modeling4.rds") %>% 
  filter(life_stage == "migration") %>% 
  group_by(individual_local_identifier) %>% 
  arrange(start_timestamp, .by_group = T) %>% 
  summarize(migration_start = min(days_since_tagging),
            migration_end = max(days_since_tagging)) %>% 
  ungroup() %>% 
  summarize(avg_migration_start = mean(migration_start) %>% round(), 
            avg_migration_end = mean(migration_end) %>% round()) %>% 
  rename(xmin = avg_migration_start,
         xmax = avg_migration_end) %>% 
  mutate(ymin = -0.35,
         ymax = 0.48) %>% #format in case I want a shaded rectangle for the migration period.
  as.data.frame()



# Extract the summary of the smooth term 
smooth_effects <- m_inla$summary.random$age_group %>% 
  #back transform the age values
  mutate(age = ID * attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:scale') +
           attr(data[,colnames(data) == "days_since_tagging_z"],'scaled:center') )

# Plot the smooth term....
X11(width = 3.42, height = 3)
(s <- ggplot(smooth_effects, aes(x = age, y = mean)) +
    #geom_rect(data = population_migr_period, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
    #          fill="black", alpha=0.1, inherit.aes = FALSE) + #add migration period
    geom_line(color = "#0d0887", linewidth = 0.4) +
    geom_ribbon(aes(ymin = `0.025quant`, ymax = `0.975quant`), fill = "#0d0887", alpha = 0.12) +
    geom_hline(yintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_vline(xintercept = c(population_migr_period$xmin, population_migr_period$xmax),
               linetype = "dotted",  color = "black", linewidth = 0.5, show.legend = TRUE) +
    labs(y = "Effect Size", x = "Days since tagging") +
    ggtitle("b") +
    ylim(-0.35, 0.48) +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"), #top, right, bottom, left
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(face = "bold"), # make title bold
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2)))
)


#combine the two plots for the linear and the non-linear terms
#X11(width = 9, height = 3)
#model_output_p <- grid.arrange(coefs, s, nrow = 1, widths = c(0.6, 0.4))

#ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_output_2min.pdf", 
#       device = "pdf", width = 9, height = 3, dpi = 600)

#combine the two plots to make the multi-panel plot for all models
X11(width = 6.7, height = 2)
model_output_p <- grid.arrange(coefs, s, nrow = 1, widths = c(0.6, 0.4))

ggsave(plot = model_output_p, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_multi_panel.pdf", 
       device = "pdf", width = 6.7, height = 2, dpi = 600)

#### ind_specific coefficients plot -----------------------------------------------------------------------------

#extract handedness for each individual. to use for coloring. 
handedness <- circling_data %>% 
  group_by(individual_local_identifier) %>% 
  slice(1) %>% 
  ungroup() %>% 
  select(individual_local_identifier, laterality_dir_ind)


#extract random effects,  ID is for the unique individuals
#mean_pitch_mean_z
random_effects_pitch <- m_inla$summary.random$individual_local_identifier

#abs_cum_yaw_z
random_effects_yaw <- m_inla$summary.random$individual_local_identifier2

#wind_speed_z
random_effects_wind <- m_inla$summary.random$individual_local_identifier3

#extract unique individual IDs from original data
#ind_IDs <- unique(data$individual_local_identifier)

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
three_vars$laterality_dir_ind <- factor(three_vars$laterality_dir_ind, levels = c("left_handed", "ambidextrous", "right_handed"))
three_vars$ID <- reorder(three_vars$ID, desc(three_vars$laterality_dir_ind))


X11(width = 7, height = 8)
(coefs_inds <- ggplot(three_vars, aes(x = coef, y = ID, color = laterality_dir_ind)) +
    geom_vline(data = filter(three_vars, variable == "mean_pitch_mean_z"), 
               aes(xintercept = graph %>% filter(Factor == "mean_pitch_mean_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", linewidth = 0.5) + 
    geom_vline(data = filter(three_vars, variable == "abs_cum_yaw_z"), 
               aes(xintercept = graph %>% filter(Factor == "abs_cum_yaw_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", linewidth = 0.5) +  
    geom_vline(data = filter(three_vars, variable == "wind_speed_z"), 
               aes(xintercept = graph %>% filter(Factor == "wind_speed_z") %>% pull(Estimate)), linetype = "dashed", color = "gray75", linewidth = 0.5) +  
    geom_point(size = 2, position = position_dodge(width = .7))  +
    geom_linerange(aes(xmin = lower, xmax = upper), size = 0.8, position = position_dodge(width = .7)) +
    # scale_color_manual(values = c("left_handed" = "#33638DFF" , "ambidextrous" = "#B8DE29FF", "right_handed" =  "#20A387FF"),
    #                    labels = c("Left-handed", "Ambidextrous", "Right-handed"),
    #                    name = "Handedness") +
    scale_color_manual(values = c("left_handed" = "#9c179e" , "ambidextrous" = "#fb9f3a", "right_handed" =  "#0d0887"),
                       labels = c("Left-handed", "Ambidextrous", "Right-handed"),
                       name = "Handedness") +
    scale_y_discrete(labels = levels(three_vars$ID)) +
    labs(x = "Estimate", y = "Individual ID") +
    theme_classic() +
    theme(text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2))) + #increase distance between x-axis values and title
    facet_wrap(~ variable, scales = "free_x", labeller = as_labeller(c(
      "mean_pitch_mean_z" = "Average pitch",
      "abs_cum_yaw_z" = "Absolute cumulative yaw",
      "wind_speed_z" = "Wind speed"
    )) # Separate panels for each variable
    ))

ggsave(plot = coefs_inds, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_coeffs_inds_2min_plasma.pdf", 
       device = "pdf", width = 7, height = 8, dpi = 600)


#### plot the interaction terms -----------------------------------------------------------------------------

## yaw:pitch----------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "pitch_yaw")

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
    scale_fill_gradientn(#colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
      #colors = c("#0d0887", "#7e03a8" ,"#cc4778", "#f89540", "#f0f921"),
      colors = c("#0d0887", "#7e03a8", "white" ,"#d5546e", "#fdae32", "#f0f921"),
      values = c(0, 0.3, 0.5, 0.6, 0.7, 1),
      limits = c(0, 1),
      na.value = "white",
      name = "Probability\n of laterality") +
    guides(fill = guide_colourbar(title.vjust = 2.5)) + # the legend title needs to move up a bit
    labs(x = "Absolute cumulative yaw", y = "Average pitch") +
    ggtitle("c") +
    theme_classic() +
    theme(plot.margin = margin(0, 6, 0, 0, "pt"),
          text = element_text(size = 8),
          legend.direction="vertical",
          legend.position = "right",
          legend.key.width=unit(.2,"cm"),
          legend.key.height=unit(.4,"cm"),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          plot.title = element_text(face = "bold"), # make title bold
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2))) + # increase distance between x-axis values and title
    scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
    scale_y_continuous(expand = c(0, 0))
)

#ggsave(plot = pred_py, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/yaw_pitch.png", 
#       device = "png", width = 7, height = 2.1, dpi = 600)



## wind:yaw----------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "wind_yaw")

preds <- data.frame(yaw = data[na_rows,"abs_cum_yaw"],
                    wind = data[na_rows,"wind_speed"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 7, height = 2)
pred_wy <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#0d0887", "#7e03a8", "white" ,"#d5546e", "#fdae32", "#f0f921"),
                       values = c(0, 0.3, 0.5, 0.6, 0.7, 1),
                       limits = c(0, 1),
                       na.value = "white",
                       name = "Probability\n of laterality") +
  guides(fill = guide_colourbar(title.vjust = 2.5)) + # the legend title needs to move up a bit
  labs(x = "Absolute cumulative yaw", y = "Wind speed") +
  ggtitle("d") +
  theme_classic() +
  theme(plot.margin = margin(0, 6, 0, 0, "pt"),
        text = element_text(size = 8),
        legend.direction="vertical",
        legend.position = "right",
        legend.key.width=unit(.2,"cm"),
        legend.key.height=unit(.4,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(face = "bold"), # make title bold
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 2))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))

#ggsave(plot = pred_wy, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/wind_yaw.png", 
#       device = "png", width = 7, height = 2.1, dpi = 600)



##wind:pitch----------------------------------------------------------------

#extract information for rows that had NAs as response variables
na_rows <- which(data$interaction == "wind_pitch")

preds <- data.frame(pitch = data[na_rows,"mean_pitch_mean"],
                    wind = data[na_rows,"wind_speed"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"]) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 2.23, height = 1.67)
pred_wp <- preds %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#0d0887", "#7e03a8", "white" ,"#d5546e", "#fdae32", "#f0f921"),
                       values = c(0, 0.3, 0.5, 0.6, 0.7, 1),
                       limits = c(0, 1),
                       na.value = "white",
                       name = "Probability\n of laterality") +
  guides(fill = guide_colourbar(title.vjust = 2.5)) + # the legend title needs to move up a bit
  labs(x = "Average pitch", y = "Wind speed") +
  ggtitle("e") +
  theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
        text = element_text(size = 8),
        legend.direction="vertical",
        legend.position = "right",
        legend.key.width=unit(.2,"cm"),
        legend.key.height=unit(.4,"cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        plot.title = element_text(face = "bold"), # make title bold
        panel.grid.minor = element_line(color = "white"),
        axis.title.x = element_text(margin = margin(t = 2))) + # increase distance between x-axis values and title
  scale_x_continuous(expand = c(0, 0)) + #remove space between the raster and the axis
  scale_y_continuous(expand = c(0, 0))

#ggsave(plot = pred_wp, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/wind_pitch.png", 
#       device = "png", width = 7, height = 2.1, dpi = 600)



#combine the three plots for the coefficients and the interaction term-------------- VERTICAL
#X11(width = 4.5, height = 6.5)
#combined <- pred_py + pred_wy + pred_wp & theme(legend.position = "bottom")
#(p <- combined + plot_layout(guides = "collect", nrow = 3))

#ggsave(plot = p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/difficulty_model_ws_interactions_2min.pdf", 
#       device = "pdf", width = 4.5, height = 6.5, dpi = 600)

#combine the three plots. horizontal
X11(width = 6.7, height = 2)
combined <- pred_py + pred_wy + pred_wp & theme(legend.position = "right")
(p <- combined + plot_layout(guides = "collect", nrow = 1))

ggsave(plot = p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/interactions_multi_panel.pdf", 
       device = "pdf", width = 6.7, height = 2, dpi = 600)


#-----------------------------------------------------------------------------------------------------------------------
## Step 5.2: Does laterality help with better performance when individuals are not experienced? migration performance #####
#-----------------------------------------------------------------------------------------------------------------------

#### data prep -----------------------------------------------------------------------------

#read in thinned liberality data
circling_data <- readRDS("thinned_laterality_for_modeling4.rds")

##open migration data from L04c_migr_metrics.r
migr_hrly <- readRDS("hourly_migr_metrics_gps_vedba.rds")

#extract days of migration
migr_days <- migr_hrly %>% 
  distinct(ind_day)

#calculate hourly and daily summaries for wind and wobble using the laterality data, then append migration dataframe
migr_hrly_ww <- circling_data %>% 
  mutate(ind_day =  paste0(individual_local_identifier, "_", as.character(unique_date))) %>% 
  filter(ind_day %in% migr_days$ind_day) %>%  #subset for days in the migration data 
  mutate(dt_1hr = round_date(start_timestamp, "1 hour")) %>%  #assign unique hour
  group_by(individual_local_identifier, unique_date, dt_1hr) %>%  
  mutate(hrly_mean_wind_speed = mean(wind_speed, na.rm = T), #summarize wind for every hour. I have a unique wind value for each hour anyway
         hrly_mean_cum_yaw = mean(abs_cum_yaw, na.rm = T),
         hrly_max_cum_yaw = max(abs_cum_yaw, na.rm = T),
         hrly_max_mean_pitch = max(mean_pitch_mean, na.rm = T),
         hrly_mean_mean_pitch = mean(mean_pitch_mean, na.rm = T),) %>% 
  ungroup() %>% 
  group_by(individual_local_identifier, unique_date) %>% #summarize wind for every hour
  mutate(daily_max_wind = max(wind_speed, na.rm = T),
         daily_mean_wind = mean(wind_speed, na.rm = T),
         daily_max_cum_yaw = max(abs_cum_yaw, na.rm = T),
         daily_mean_cum_yaw = mean(abs_cum_yaw, na.rm = T),
         daily_max_mean_pitch = max(mean_pitch_mean, na.rm = T),
         daily_mean_mean_pitch = mean(mean_pitch_mean, na.rm = T),
         #I have already calculated daily LI, but in the previous round, I used mode of laterality in the models. so calculate the mode for making comparisons.
         mode_laterality = getmode(laterality_dir)) %>% 
  ungroup() %>% 
  select(4, 42, 59, 64, 69:74, 85:98) %>% 
  group_by(individual_local_identifier, unique_date, dt_1hr) %>% 
  slice(1) %>% #just keep one row per hour 
  #bind to the daily migration data. get rid of days in the migr_hrly data that don't have laterality data
  left_join(migr_hrly) %>% 
  as.data.frame()

saveRDS(migr_hrly_ww, file = "data_migration_performance_models_2min_hrly.rds") #this only has days that have laterality info

migr_daily_ww <- migr_hrly_ww %>% #ww: wind + wobble
  select(-contains("hrly")) %>%  #remove hourly data
  group_by(individual_local_identifier, unique_date) %>%  #keep one row per day
  slice(1) %>% 
  as.data.frame()

saveRDS(migr_daily_ww, file = "data_migration_performance_models_2min_daily.rds") #this only has days that have laterality info

#### exploratory plots-----------------------------------------------------------------------------

migr_long <- migr_daily_ww %>%
  pivot_longer(cols = c(daily_max_wind, daily_mean_wind, daily_max_cum_yaw, 
                        daily_mean_cum_yaw, daily_max_mean_pitch, 
                        daily_mean_mean_pitch, daily_distance, daily_avg_speed, daily_avg_altitude, daily_max_altitude,
                        daily_mean_vedba, daily_mean_vedba, daily_max_vedba, daily_min_vedba, daily_IQR_vedba),
               names_to = "measure", values_to = "value")

ggplot(migr_long, aes(x = mode_laterality, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 0.2) +
  labs(x = "Laterality", y = "Value") +
  theme_classic() +
  facet_wrap(~ measure, scales = "free_y")


ggplot(migr_long, aes(x = laterality_dir_ind, y = value)) +
  geom_violin(trim = FALSE, fill = "lightblue", alpha = 0.5) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.8, size = 0.2) +
  labs(x = "Laterality", y = "Value") +
  theme_classic() +
  facet_wrap(~ measure, scales = "free_y")


### Model migration performance as a function of laterality -----------------------------------------------------

migr_daily <- readRDS("data_migration_performance_models_2min_daily.rds") %>% 
  mutate(laterality_dir_day = as.character(laterality_dir_day), #convert to character, so that"ambidextrous" is the reference level
         laterality_bi_day = ifelse(laterality_dir_day == "ambidextrous", 0, 1),
         laterality_bi_day_mode =  ifelse(mode_laterality == "ambidextrous", 0, 1)) 

#replace very large values of max vertical speed with NAs
migr_daily[which(migr_daily$daily_max_vert_speed > 10), "daily_max_vert_speed"] <- NA

#z-transform the variables
data_m <- migr_daily %>%
  mutate(daily_abs_LI = abs(laterality_bank_day),
         laterality_bi_ind = ifelse(laterality_dir_ind == "ambidextrous", "ambidextrous", "handed"),
         laterality_bi_stage = ifelse(laterality_dir_stage == "ambidextrous", "ambidextrous", "handed")) %>% #because this is migration data only, there should be one assignment per individual
  #make sure to do all the filtering before scaling the variables!!!!
  mutate(across(contains("daily") , ~scale(.), .names = "{.col}_z")) %>%
  as.data.frame()

#check for autocorrelation
data_m %>% 
  dplyr::select(ends_with("_z")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #correlated: mean and max wind, max cum yaw and mean cum yaw, daily distance an average speed, mean vedba and IQR vedba; max and avg altitude;
# sd of vertical speed with avg altitude (0.57) and max altitude (0.64)

#model: separately for daily_mean_vedba, daily_max_altitude, daily_distance, daily_mean_cum_yaw, daily_mean_mean_pitch, daily_max_altitude
#vertical speed doesnt show a pattern. Omit because it is also unnecessary to have to explain the methodology.

response_vars <- c( "daily_distance", "daily_max_altitude",
                    "daily_mean_vedba", "daily_mean_cum_yaw", 
                    "daily_mean_mean_pitch")

response_names <- c(
  expression(atop(bold("f") * phantom("                             "), "Daily distance (km)")),
  expression(atop(bold("g") * phantom("                                  "), "Max flight altitude (m)")),
  expression(atop(bold("h") * phantom("                        "), "Avg VeDBA (g)")), 
  expression(atop(bold("i") * phantom("                                        "), "Avg cumulative yaw (deg)")), 
  expression(atop(bold("j") * phantom("                        "), "Avg pitch (deg)"))
)

plots_ls <- lapply(1:length(response_vars), function(response){
  
  #model
  #formula <- paste0(response_vars[response], " ~ 1 + daily_max_wind_z + laterality_bi_day") %>% formula()
  formula <- paste0(response_vars[response], " ~ 1 + daily_max_wind_z + daily_abs_LI_z") %>% formula()
  
  #don't include ind ID as a random intercept. some individuals don't have all levels of laterality
  m_inla <- inla(formula,
                 data = data_m,
                 control.compute = list(cpo = TRUE),
                 control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted
  
  # coefficients plot 
  
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
  
  #export the output as a latex table -----------
  # Convert to LaTeX table
  latex_table <- xtable(graph)
  
  # Specify the file path for each response variable
  # file_path <- paste0("/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/tables/latex_table_", response_vars[response], ".txt")
  
  # Open a connection to the file
  #sink(file_path)
  
  # Print the LaTeX code to the file
  #print(latex_table, type = "latex", include.rownames = FALSE)
  
  # Close the connection to the file
  #sink()
  
  #plot the coefficients -----------
  
  #remove intercept for better visualization. 
  graph <- graph[-1,]
  
  #X11(width = 4.5, height = 1.6)
  ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#0d0887", size = 1.5)  +
    labs(x = if (response == 5) "Estimate" else "", 
         y = "") +
    scale_y_discrete(labels = if (response %in% c(1,4)) rev(c("Intercept", "Max. wind speed", "Abs. Laterality Index")) else c("", "")) + 
    #scale_y_discrete(labels = if (response %in% c(1,4)) rev(c("Intercept", "Max. wind speed", "Laterality")) else c("", "")) + # Only add labels for first and fourth plot
    #scale_y_discrete(labels = rev(c("Intercept", "Max. wind speed", "Laterality"))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#0d0887", linewidth = 0.5) +
    ggtitle(response_names[response]) +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          text = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          axis.title.x = element_text(margin = margin(t = 2))) #increase distance between x-axis values and title
})

#put all plots together
X11(width = 6.7, height = 2.5)
combined <- reduce(plots_ls[1:5], `+`)
(p <- combined + plot_layout(ncol = 3))

ggsave(plot = p, filename = "/home/mahle68/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/migration_models_multi_panel_continuous.pdf", 
       device = "pdf", width = 6.7, height = 2.5, dpi = 600)