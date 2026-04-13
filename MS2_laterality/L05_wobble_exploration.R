#revisiting flight stability and its relationship with age and laterality
#April 8. 2026.Lausanne, CH
#Elham Nourani, elham.nourani@unil.ch 


library(tidyverse)
library(corrr)
library(INLA)
library(gridExtra)
library(patchwork)
library(terra)
library(ggridges)

#input data from the public folder

filtered_w_LI <- readRDS("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/public_data_for_Edmond/imu_wind_laterality_bursts.rds")


#repeat data filtering and z-transformation steps from 04_data_analysis.r in the public Github repo. Modify the missing data generation code to match what needs predicting: 

#### ----------------------- filter for day since tagging and z-transform
quantile(filtered_w_LI$days_since_tagging, probs = 0.9) #284 ... 

circling_data <- filtered_w_LI %>% 
  mutate(days_since_tagging = as.numeric(days_since_tagging)) %>% 
  #filter for days since tagging. only keep the 
  filter(days_since_tagging < 284) %>% 
  #make sure to do all the filtering before scaling the variables!!!!
  mutate_at(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_sum_8sec", 
              "abs_cum_yaw", "wind_speed", "days_since_tagging", "location_lat_closest_gps_raw"),
            list(z = ~scale(.))) %>% 
  as.data.frame()

#### ----------------------- look at multi-collinearity

circling_data %>% 
  dplyr::select(c("mean_roll_sd", "mean_yaw_mean", "mean_yaw_sd", "mean_pitch_mean" , "mean_pitch_sd", "cumulative_pitch_sum_8sec", 
                  "abs_cum_yaw", "days_since_tagging", "wind_speed", "location_lat_closest_gps_raw")) %>% 
  correlate() %>% 
  corrr::stretch() %>% 
  filter(abs(r) > 0.5) #correlated variables: sd of roll and pitch, sd of yaw and cumulative yaw

#### ----------------------- model: regression with inla

#first bin the age variable
data <- circling_data %>%  
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
         individual_local_identifier = sample(circling_data$individual_local_identifier, nrow(.), replace = T) %>% 
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


#### ----------------------- model evaluation

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
    scale_y_discrete(labels = rev(c("Intercept", "Laterality", "Age", "Laterality:Age"))) + 
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

#### interaction plot -----------------------------------------------------------------------------
#extract information for rows that had NAs as response variables
na_rows <- which(is.na(data$mean_yaw_sd_z))
na_rows <- which(is.na(data$mean_pitch_sd_z))

preds <- data.frame(laterality_bi = data[na_rows,"laterality_bi"],
                    age_group = data[na_rows,"age_group"],
                    preds = m_inla$summary.fitted.values[na_rows,"mean"],
                    lower95 = m_inla$summary.fitted.values[na_rows,"0.025quant"],
                    upper95 = m_inla$summary.fitted.values[na_rows,"0.975quant"]) %>% 
  #back transform the age values
  mutate(age = age_group * sd(circling_data$days_since_tagging) + mean(circling_data$days_since_tagging),
         #preds_yaw = preds * sd(circling_data$mean_yaw_sd) + mean(circling_data$mean_yaw_sd))
         preds_pitch = preds * sd(circling_data$mean_pitch_sd) + mean(circling_data$mean_pitch_sd))
#upper = upper95 * sd(circling_data$mean_pitch_sd) + mean(circling_data$mean_pitch_sd),
#lower = lower95 * sd(circling_data$mean_pitch_sd) + mean(circling_data$mean_pitch_sd))


#plot
(pred_p <- preds %>% 
    ggplot() +
    #geom_point(aes(x = age, y = preds_yaw, color = laterality_bi), size = 0.2) +
    #geom_smooth(aes(x = age, y = preds_yaw, color = laterality_bi, fill = laterality_bi), method = "lm", se = T, alpha = 0.3, linewidth = 0.3) +
    geom_point(aes(x = age, y = preds_pitch, color = laterality_bi), size = 0.2) +
    geom_smooth(aes(x = age, y = preds_pitch, color = laterality_bi, fill = laterality_bi), method = "lm", se = T, alpha = 0.3, linewidth = 0.3) +
    scale_color_manual(values = c("0" = "gray40", "1" = "#8a2be2ff"),
                       name = "Laterality")+
    scale_fill_manual(values = c("0" = "gray40", "1" = "#8a2be2ff"),
                      guide = F)+
    #labs(x = "Age", y = "Horizontal wobble") +
    labs(x = "Age", y = "Vertical wobble") +
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
model_output_p <- grid.arrange(coefs, pred_p, nrow = 1, widths = c(0.45, 0.55))

#ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_output.pdf", 
#       device = "pdf", width = 7, height = 2.3, dpi = 600)

ggsave(plot = model_output_p, filename = "/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/wobble_exploration/sd_pitch_model_output.png", 
       device = "png", width = 7, height = 2.3, dpi = 600)


#### wobble density plots -----------------------------------------------------------------------------


stage_colors <- c(
  "post-fledging" = "#31688e",  # Dark Green
  "migration" = "#b5de2b",      # Dark Red
  "wintering" = "#1f9e89"       # Dark Cyan
)


circling_data$life_stage <- factor(
  circling_data$life_stage,
  levels = names(stage_colors)
)



(h_wobble <- ggplot(circling_data, aes(x = mean_yaw_sd, y = life_stage, color = life_stage, fill = life_stage)) + 
  stat_density_ridges(jittered_points = F, rel_min_height = .01,
                      point_shape = "|", point_size = 3, point_alpha = 0.8, 
                      calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                      scale = 1.5, bandwidth = 0.2) +
  labs(y = "", x = "Horizontal wobble") +
  scale_color_manual(values = stage_colors, guide = FALSE) +
  scale_fill_manual(values = stage_colors, guide = FALSE) +
  coord_cartesian(ylim = c(1.5,4),
                  xlim = c(-0.2, 50)) +
  theme_classic() +
  theme(legend.text = element_text(size = 7),
        legend.position = "none",
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, margin = margin(t = 1)),  # Adjust the top margin for space
        axis.title.x = element_text(margin = margin(t = 20))) )



ggplot(
  circling_data,
  aes(
    x = mean_yaw_sd,
    color = life_stage,
    fill  = life_stage
  )
) +
  geom_density(
    alpha = 0.35,
    size = 1,
    position = "identity",
    adjust = 1
  ) +
  labs(
    x = "Horizontal wobble",
    y = "Density"
  ) +
  scale_color_manual(values = stage_colors, guide = FALSE) +
  scale_fill_manual(values = stage_colors, guide = FALSE) +
  coord_cartesian(
    xlim = c(0, 50)
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20, margin = margin(t = 15)),
    axis.title.x = element_text(margin = margin(t = 20))
  )


##v wobble
(v_wobble <- ggplot(circling_data, aes(x = mean_pitch_sd, y = life_stage, color = life_stage, fill = life_stage)) + 
    stat_density_ridges(jittered_points = F, rel_min_height = .01,
                        point_shape = "|", point_size = 3, point_alpha = 0.8, 
                        calc_ecdf = FALSE, panel_scaling = FALSE, alpha = 0.2, quantile_lines = TRUE,
                        scale = 1.5, bandwidth = 0.2) +
    labs(y = "", x = "Vertical wobble") +
    scale_color_manual(values = stage_colors, guide = FALSE) +
    scale_fill_manual(values = stage_colors, guide = FALSE) +
    coord_cartesian(ylim = c(1.5,4),
                    xlim = c(-0.2, 12)) +
    theme_classic() +
    theme(legend.text = element_text(size = 7),
          legend.position = "none",
          axis.text = element_text(size = 20),
          axis.title = element_text(size = 20, margin = margin(t = 1)),  # Adjust the top margin for space
          axis.title.x = element_text(margin = margin(t = 20))) )

