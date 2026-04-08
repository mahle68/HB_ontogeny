#revisiting flight stability and its relationship with age and laterality
#April 8. 2026.Lausanne, CH
#Elham Nourani, elham.nourani@unil.ch 


library(tidyverse)
library(corrr)
library(INLA)


#input data from the public folder

filtered_w_LI <- readRDS("/Users/enourani/Library/CloudStorage/OneDrive-UniversitédeLausanne/WORK/Projects/HB_ontogeny_eobs/paper_prep/MS2_laterality/public_data_for_Edmond/imu_wind_laterality_bursts.rds")


#repeat data filtering and z-transformation steps from 04_data_analysis.r in the public Github repo. Modify the missing data generation code to match what needs predicting: 
#### ----------------------- filter for day since tagging and z-transform
quantile(filtered_w_LI$days_since_tagging, probs = 0.9) #284 ... 

circling_data <- filtered_w_LI %>% 
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

