#revisiting flight stability and its relationship with age and laterality
#April 8. 2026.Lausanne, CH
#Elham Nourani, elham.nourani@unil.ch 


library(tidyverse)
library(corrr)
library(INLA)
library(gridExtra)
library(patchwork)
library(terra)


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

#### ----------------------- new data to look at interaction between laterality and age

new_data <- expand.grid(wind_speed = seq(from = min(circling_data$wind_speed,na.rm = T), to = max(circling_data$wind_speed, na.rm = T),  length.out = 50),
                        laterality_bi = c(0, 1), #laterality values
                        days_since_tagging = seq(from = min(circling_data$days_since_tagging, na.rm = T), to = max(circling_data$days_since_tagging, na.rm = T),  length.out = 50)) %>% 
  mutate(#set the dependent variables to NA. keep both the z-transormed and original columns. The original will be used for log-transformation of the response
    mean_yaw_sd_z = NA, 
    mean_pitch_sd_z = NA, 
    mean_yaw_sd = NA, 
    mean_pitch_sd = NA, 
    #calculate z scores
    wind_speed_z = (wind_speed - attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "wind_speed_z"],'scaled:scale'),
    days_since_tagging_z = (days_since_tagging - attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:center'))/attr(circling_data[,colnames(circling_data) == "days_since_tagging_z"],'scaled:scale'),
    #randomly assign individual IDs
    individual_local_identifier = sample(circling_data$individual_local_identifier, nrow(.), replace = T))

#### ----------------------- model: binomial logistic regression with inla

#bin the age variable and append the new data
data <- circling_data %>% 
  #only keep columns that exist in the new_data
  dplyr::select(intersect(colnames(new_data), colnames(.))) %>% 
  full_join(new_data) %>% 
  mutate(individual_local_identifier2 = individual_local_identifier, #repeat individual ID column to be used in the model formula for random effects
         individual_local_identifier3 = individual_local_identifier,
         age_group = inla.group(days_since_tagging_z, n = 100, method = "quantile"),
         laterality_bi = as.factor(laterality_bi),
         individual_local_identifier = as.factor(individual_local_identifier))

#### run the models one at a time! for pitch and yaw separately
m_inla_p <- inla(mean_pitch_sd_z ~ 1 + laterality_bi * age_group * wind_speed_z +
                   f(individual_local_identifier, age_group, model = "iid"),
                 data = data,
                 control.compute = list(cpo = TRUE),
                 control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

#log-transform instead of z-transform, to make sure response values are positive
m_inla_p <- inla(log(mean_pitch_sd) ~ 1 + laterality_bi * age_group * wind_speed_z +
                   f(individual_local_identifier, age_group, model = "iid"),
                 data = data,
                 control.compute = list(cpo = TRUE),
                 control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted

m_inla_y <- inla(mean_yaw_sd_z ~ 1 + laterality_bi * age_group * wind_speed_z +
                   f(individual_local_identifier, age_group, model = "iid"),
                 data = data,
                 control.compute = list(cpo = TRUE),
                 control.predictor = list(link = 1, compute = TRUE)) #compute=t means that NA values will be predicted


### ----------------------- coefficient plts
# posterior means of coefficients
graph <- as.data.frame(summary(m_inla_p)$fixed)
colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")

#graph$Model<-i
graph$Factor <- rownames(graph)


#### pitch_model: interaction plot -----------------------------------------------------------------------------

#extract information for rows that had NAs as response variables
na_rows <- which(is.na(data$mean_pitch_sd_z))

preds_p <- data.frame(days_since_tagging = data[na_rows,"days_since_tagging"],
                      wind = data[na_rows,"wind_speed"],
                      preds = m_inla_p$summary.fitted.values[na_rows,"mean"]) %>% 
  #back-transform the z-scores back to sd_pitch values
  #mutate(preds = preds * attr(circling_data[,colnames(circling_data) == "mean_pitch_sd_z"],'scaled:center') +
  #         attr(circling_data[,colnames(circling_data) == "mean_pitch_sd_z"],'scaled:scale')) %>% 
  #do some interpolation to smooth the plot
  terra::rast(type = "xyz") %>%
  focal(w = 3, fun = "mean", na.policy = "all", na.rm = T) %>% 
  as.data.frame(xy = T) %>%
  rename(preds = focal_mean)

#plot
#X11(width = 7, height = 2)
pred_py <- preds_p %>% 
  ggplot() +
  geom_tile(aes(x = x, y = y, fill = preds)) +
  scale_fill_gradientn(colors = c("#440154FF", "#39568CFF" ,"#20A387FF", "#83e22b", "#FDE725FF"),
                       values = c(0, 0.5, 0.6, 0.7, 1),
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


#combine the two plots for the coefficients and the interaction term
X11(width = 7, height = 2.3)
model_output_p <- grid.arrange(coefs, pred_py, nrow = 1, widths = c(0.4, 0.6))

#ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_yaw_model_output.pdf", 
#       device = "pdf", width = 7, height = 2.3, dpi = 600)

ggsave(plot = model_output_p, filename = "/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/figures/sd_pitch_model_output.pdf", 
       device = "pdf", width = 7, height = 2.3, dpi = 600)





#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
####### apr 8 stuff
data <- circling_data %>% 
  mutate(individual_local_identifier2 = individual_local_identifier, #repeat individual ID column to be used in the model formula for random effects
         individual_local_identifier3 = individual_local_identifier,
         age_group = inla.group(days_since_tagging_z, n = 100, method = "quantile"), 
         lat_group = inla.group(location_lat_closest_gps_raw, n = 10, method = "quantile"))  


m_yaw <- inla(mean_yaw_sd ~ laterality_bi * age_group + wind_speed +
                f(individual_local_identifier, age_group, model = "iid"),
              data = data,
              control.compute = list(cpo = T),
              control.predictor = list(link = 1, compute = T))

### ----------------------- coefficient plts
# posterior means of coefficients
graph <- as.data.frame(summary(m_yaw)$fixed)
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

#plot the coefficients
quartz(width = 6.2, height = 4)
(coefs <- ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype="dashed", 
               color = "gray75", linewidth = 0.5) +
    geom_point(color = "#0d0887", size = 1.5)  +
    labs(x = "Estimate", y = "") +
    scale_y_discrete(labels = rev(c("Intercept", "Laterality", "Age", "wind speed", "Laterality: Age"))) +
    geom_linerange(aes(xmin = Lower, xmax = Upper),color = "#0d0887", linewidth = 0.5) +
    #ggtitle("a") +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"), #top, right, bottom, left
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          legend.title = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          plot.title = element_text(face = "bold"), # make title bold
          axis.title.x = element_text(margin = margin(t = 2))) #increase distance between x-axis values and title
)


#interaction plot
na_rows <- which(is.na(data$mean_yaw_sd_z))

preds <- data.frame(laterality_bi = data[na_rows,"laterality_bi"],
                    age_group = data[na_rows,"age_group"],
                    wind_speed = data[na_rows, "wind_speed"],
                    preds = m_inla_y$summary.fitted.values[na_rows,"mean"],
                    lower95 = m_inla_y$summary.fitted.values[na_rows,"0.025quant"],
                    upper95 = m_inla_y$summary.fitted.values[na_rows,"0.975quant"]) %>% 
  #back transform the age values
  mutate(age = age_group * sd(circling_data$days_since_tagging) + mean(circling_data$days_since_tagging),
         preds_yaw = preds * sd(circling_data$mean_yaw_sd) + mean(circling_data$mean_yaw_sd))
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

model_output_p <- grid.arrange(coefs, pred_p, nrow = 1, widths = c(0.4, 0.6))
