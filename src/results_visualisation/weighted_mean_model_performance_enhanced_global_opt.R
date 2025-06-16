source("code/ensemble_prediction.R")
library(dplyr)
library(tidyr)
library(ggplot2)


### preparing dataframes ####

#### Predicted ####

prediction_df_es_new <- read.csv("new_fitter/ES_fit/prediction_df_es.csv")
prediction_df_es_new <- prediction_df_es_new %>%
  mutate(Fit ="Location specific")

prediction_df_de_new <- read.csv ("new_fitter/DE_fit/prediction_df_de.csv")
prediction_df_de_new <- prediction_df_de_new %>%
  mutate(Fit ="Location specific")

prediction_df_apple_new <- read.csv("new_fitter/de_es_fit/prediction_df_apple.csv")
prediction_df_apple_new <- prediction_df_apple_new %>%
  mutate(Fit ="Apple fit")

prediction_df_new <- rbind(prediction_df_apple_new,
                       prediction_df_de_new,
                       prediction_df_es_new) %>%
  mutate(Cultivar = recode(Cultivar, "DelaRiega" = "De la Riega")) 

# add id 
prediction_df_id <- prediction_df_new %>%
  mutate(id = paste(Cultivar,location, Fit, sep = "_"))

# Ensure the Year is treated as numeric and clean any potential duplicates
predicted <- prediction_df_id %>%
  mutate(Year = as.numeric(Year)) %>%
  group_by(id, Year) %>%
  mutate(row_num = row_number()) %>%
  # Add distinct() to handle any potential duplicates
  distinct(id, Year, row_num, .keep_all = TRUE) %>%
  pivot_wider(
    id_cols = c(id, Year, Observed),
    names_from = row_num,
    values_from = Predicted,
    names_prefix = "R"
  ) %>%
  ungroup()

##### Confidence #####

# calculate iqr for calibration + validation because of small sample size
iqr_df <- prediction_df_id %>% 
  filter(round == 1) %>% 
  group_by(id) %>% 
  summarise(iqr = IQR(Observed))

confidence <- prediction_df_id %>% 
  # confidence should be based only on validation
  filter(data_type == 'Validation') %>% 
  group_by(id, round) %>% 
  summarise(rmse = chillR::RMSEP(predicted = Predicted, observed = Observed)) %>% 
  ungroup() %>% 
  # add the iqr information
  merge(iqr_df, by = 'id') %>% 
  mutate(rpiq = iqr / rmse) %>% 
  select(id, round, rpiq) %>% 
  # bring to wide format
  pivot_wider(names_from = round, values_from = rpiq, names_prefix = 'R')

weighted_mean <- weighted_mean_with_fail(predicted,confidence)

##### Plot Predicted vs. Observed using weighted mean.#####
# Restructure dataframe

weighted_mean_wider <- weighted_mean %>%
  separate(id, into = c("Cultivar", "Location", "Fit"), sep = "_") %>%
  select(Cultivar, Location, Fit, Year, Observed, weighted_pred, sd)

# factor location
weighted_mean_wider <- weighted_mean_wider %>%
  mutate(Cultivar = factor(Cultivar, levels = unique(Cultivar[order(Location, Cultivar)])),
         Location = factor(Location, levels = c("Germany", "Spain")))


performance_stats <- weighted_mean_wider %>%
  group_by(Cultivar,Fit) %>% 
  summarise(RMSE = chillR::RMSEP(predicted = weighted_pred, observed = Observed),
            RPIQ = chillR::RPIQ(weighted_pred,Observed),
            Bias = mean(weighted_pred - Observed)) %>% 
  reshape2::melt(id.vars = c("Cultivar","Fit")) %>% 
  reshape2::dcast(Cultivar~ Fit + variable,value.var = "value") %>%
  rename_with(~ gsub(" ", "_", .))


## summarize per location to compare overall performance between German and Spanish 
performance_stats <- weighted_mean_wider %>%
  group_by(Cultivar,Fit) %>% 
  summarise(RMSE = chillR::RMSEP(predicted = weighted_pred, observed = Observed),
            RPIQ = chillR::RPIQ(weighted_pred,Observed),
            Bias = mean(weighted_pred - Observed)) %>% 
  reshape2::melt(id.vars = c("Cultivar","Fit")) %>% 
  reshape2::dcast(Cultivar~ Fit + variable,value.var = "value") %>%
  rename_with(~ gsub(" ", "_", .))

performance_stats_location <- weighted_mean_wider %>%
  group_by(Location,Fit) %>% 
  summarise(RMSE = chillR::RMSEP(predicted = weighted_pred, observed = Observed),
            RPIQ = chillR::RPIQ(weighted_pred,Observed),
            Bias = mean(weighted_pred - Observed)) %>% 
  reshape2::melt(id.vars = c("Location","Fit")) %>% 
  reshape2::dcast(Location~ Fit + variable,value.var = "value") %>%
  rename_with(~ gsub(" ", "_", .))


#plot

p_weighted_mean <- weighted_mean_wider %>%
  ggplot(aes(x = Observed, y = weighted_pred)) +
  # Add reference line
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = "gray40") +
  # Add error bars with transparency
  geom_errorbar(aes(ymin = weighted_pred - sd, 
                    ymax = weighted_pred + sd, 
                    color = Fit, alpha = Fit), 
                #alpha = 0.7, 
                width = 0.5,
                linewidth = 1 ) +
  # Add points
  geom_point(aes(color = Fit, alpha = Fit), 
             size = 3) +
             #alpha = 0.3) +
  # Set colors
  scale_color_manual(name = "Calibration approaches:",
                     values = c("Apple fit" = "#E64B35", 
                                "Location specific" = "#4DBBD5"),
                     labels = c("Species-specific","Location-specific")) +
  
  # Set alpha: lower for blue, higher for red
  scale_alpha_manual(values = c("Apple fit" = 0.6, "Location specific" = 0.4)) +
  
  ylim(c(60,180))+ xlim(c(60,180))+
  
  geom_text(data = performance_stats, aes(x = 75, y = 150, label = paste0('RMSE: ', format(round(Apple_fit_RMSE, digits = 1), nsmall = 1),
                                                                          ' (', format(round(Location_specific_RMSE,1), nsmall = 1), ')')), hjust = 0,size = 3)+
  geom_text(data = performance_stats, aes(x = 75, y = 145, label = paste0('RPIQ: ', format(round(Apple_fit_RPIQ, digits = 1), nsmall = 1),
                                                                          ' (', format(round(Location_specific_RPIQ,1), nsmall = 1), ')')), hjust = 0,size = 3)+
  geom_text(data = performance_stats, aes(x = 75, y = 140, label = paste0('bias: ', format(round(Apple_fit_Bias, digits = 1), nsmall = 1),
                                                                          ' (', format(round(Location_specific_Bias,1), nsmall = 1), ')')), hjust = 0,size = 3)+
  
  #format axes
  scale_y_continuous(breaks = c(75, 96, 117, 138, 160), 
                     labels = c('15 Mar', '5 Apr', '26 Apr', '17 May', '8 Jun')) +
  scale_x_continuous(breaks = c(75, 96, 117, 138, 160), 
                     labels = c('15 Mar', '5 Apr', '26 Apr', '17 May', '8 Jun')) +
  # Labels
  labs(x = 'Observed Bloom Date',
       y = 'Predicted Bloom Date') +
  # Theming
  theme_bw(base_size = 14) +
  theme(
    #panel formatting
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "gray20"),
    
    # Axis formatting
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(),
    
    # legend formatting
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = -10, unit = "pt"),
    plot.margin = margin (1.5,1.5,1.5,2),
    
    # Facet formatting
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    
    # Add proper spacing
    panel.spacing = unit(0.5, "lines")
  ) +
  facet_wrap(~ Cultivar, ncol = 4, scales = "fixed") +
  guides(alpha = "none")


# save
ggsave('plots/model_performance/weighted_mean_prediction_dpi.jpeg',
       p_weighted_mean, height = 24, width = 19,
       dpi = 500, units = 'cm', device = 'jpeg')

##### calculate sd between splits #####

sd_splits <- weighted_mean_wider %>%
  group_by(Cultivar, Fit) %>%
  summarize(mean_sd = mean(sd),
            sd_cul = sd (sd))

