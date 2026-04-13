# script for testing whether individual-level side bias (handed vs ambidextrous)
# is associated with migration performance in European honey buzzards.
# Individual-level bias classification is taken from laterality_dir_ind, which
# is the same classification used to color the distribution of laterality index figure (L04b_plots.r).
# Elham Nourani, PhD.
# Apr. 13. 2025. Konstanz, DE.

library(tidyverse)
library(INLA)
library(patchwork)
library(ggridges)

setwd("/home/enourani/ownCloud - enourani@ab.mpg.de@owncloud.gwdg.de/Work/Projects/HB_ontogeny_eobs/R_files/")

#---------------------------------------------------------------------------------
## Step 1: Load data and create bias group variable                          #####
#---------------------------------------------------------------------------------

#read in the daily migration performance dataset (from 04_data_analysis.r Step 2 / Edmond repository)
migr_daily <- readRDS("data_migration_performance_models_2min_daily2.rds")

#collapse left_handed and right_handed into a single "handed" group, retaining ambidextrous
#laterality_dir_ind is the individual-level classification (constant per individual)
#computed from the overall laterality index across all bursts, with thresholds at ±0.25 LI
data_m <- migr_daily %>%
  mutate(
    daily_abs_LI = abs(laterality_bank_day),
    bias_group = ifelse(laterality_dir_ind == "ambidextrous", "ambidextrous", "handed") %>%
      factor(levels = c("ambidextrous", "handed"))
  ) %>%
  mutate(across(c(daily_distance, daily_max_altitude, daily_mean_vedba,
                  daily_mean_cum_yaw, daily_mean_mean_pitch, daily_max_wind),
                ~scale(.), .names = "{.col}_z")) %>%
  as.data.frame()

#how many individuals per group?
data_m %>%
  group_by(individual_local_identifier) %>%
  slice(1) %>%
  ungroup() %>%
  count(bias_group)

#---------------------------------------------------------------------------------
## Step 2: Exploratory plots                                                  #####
#---------------------------------------------------------------------------------

#one row per individual-day for plotting
response_labels <- c(
  daily_distance     = "Daily distance (km)",
  daily_max_altitude = "Max flight altitude (m)",
  daily_mean_vedba   = "Avg VeDBA (g)",
  daily_mean_cum_yaw = "Avg total yaw (deg)",
  daily_mean_mean_pitch = "Avg pitch (deg)"
)

bias_colors <- c("ambidextrous" = "#fb9f3a", "handed" = "#0d0887")

plot_ls <- lapply(names(response_labels), function(rv) {
  ggplot(data_m %>% drop_na(all_of(rv)),
         aes(x = .data[[rv]], y = bias_group, fill = bias_group, color = bias_group)) +
    stat_density_ridges(quantile_lines = TRUE, quantiles = 2,
                        rel_min_height = 0.01, alpha = 0.4,
                        jittered_points = TRUE,
                        point_shape = "|", point_size = 1, point_alpha = 0.6) +
    scale_fill_manual(values = bias_colors) +
    scale_color_manual(values = bias_colors) +
    labs(x = response_labels[rv], y = "") +
    theme_classic() +
    theme(text = element_text(size = 8),
          legend.position = "none",
          axis.title.x = element_text(margin = margin(t = 2)))
})

X11(width = 10, height = 4)
reduce(plot_ls, `+`) + plot_layout(nrow = 1)


#---------------------------------------------------------------------------------
## Step 3: INLA models — migration performance ~ bias group + wind            #####
#---------------------------------------------------------------------------------

#bias_group is an individual-level (between-individual) factor.
#individual random effect captures the repeated daily observations within each individual.
#daily_max_wind_z is included as a covariate following 04_data_analysis.r Step 2.

response_vars <- c("daily_distance", "daily_max_altitude",
                   "daily_mean_vedba", "daily_mean_cum_yaw",
                   "daily_mean_mean_pitch")

response_names <- list(
  expression(atop(bold("a") * phantom("                       "), "Daily distance (km)")),
  expression(atop(bold("b") * phantom("                            "), "Max flight altitude (m)")),
  expression(atop(bold("c") * phantom("                 "), "Avg VeDBA (g)")),
  expression(atop(bold("d") * phantom("                        "), "Avg total yaw (deg)")),
  expression(atop(bold("e") * phantom("                  "), "Avg pitch (deg)"))
)

models_ls <- lapply(seq_along(response_vars), function(i) {

  formula <- as.formula(paste0(
    response_vars[i], "_z ~ 1 + daily_max_wind_z + bias_group + ",
    "f(individual_local_identifier, model = 'iid')"
  ))

  m <- inla(formula,
            data = data_m,
            control.compute = list(cpo = TRUE),
            control.predictor = list(link = 1, compute = TRUE))

  list(model = m, response = response_vars[i])
})

saveRDS(models_ls, file = "bias_group_migration_performance_models.rds")

#model fit summary
lapply(models_ls, function(x) {
  data.frame(
    response = x$response,
    CPO = mean(x$model$cpo$cpo, na.rm = TRUE),
    Mlik = as.numeric(x$model$mlik[, 1][2])
  )
}) %>%
  bind_rows()


#---------------------------------------------------------------------------------
## Step 4: Coefficient plots                                                  #####
#---------------------------------------------------------------------------------

coef_plots_ls <- lapply(seq_along(models_ls), function(i) {

  graph <- as.data.frame(summary(models_ls[[i]]$model)$fixed)
  colnames(graph)[which(colnames(graph) %in% c("0.025quant", "0.975quant"))] <- c("Lower", "Upper")
  colnames(graph)[which(colnames(graph) %in% c("mean"))] <- "Estimate"
  graph$Factor <- rownames(graph)

  graph <- graph %>%
    filter(Factor != "(Intercept)") %>%
    mutate(Factor = factor(Factor, levels = rev(unique(Factor))))

  ggplot(graph, aes(x = Estimate, y = Factor)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray75", linewidth = 0.5) +
    geom_point(color = "#0d0887", size = 1.5) +
    geom_linerange(aes(xmin = Lower, xmax = Upper), color = "#0d0887", linewidth = 0.5) +
    scale_y_discrete(labels = rev(c("Max wind speed", "Handed vs\nAmbidextrous"))) +
    labs(x = if (i == 5) "Estimate" else "",
         y = "") +
    ggtitle(response_names[[i]]) +
    theme_classic() +
    theme(plot.margin = margin(0, 0, 0, 0, "pt"),
          text = element_text(size = 8),
          panel.grid.minor = element_line(color = "white"),
          plot.title = element_text(face = "bold"),
          axis.title.x = element_text(margin = margin(t = 2)),
          axis.text.y = if (i == 1) element_text() else element_blank())
})

X11(width = 6.7, height = 2.5)
combined <- reduce(coef_plots_ls, `+`)
(p <- combined + plot_layout(nrow = 1))
