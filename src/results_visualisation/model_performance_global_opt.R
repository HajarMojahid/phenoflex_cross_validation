
# load libraries
library(chillR)
library(dplyr)
library(tidyr)
library(ggplot2)


# load prediction dataframes for all fits

val_es <- read.csv("old_fitter/05_summary_val_ES.csv") %>%
  mutate(fit = "location_fit",
         Cultivar =  Cultivar %>%
           sub("^data_", "", .) ) # remove "data_" prefix


val_de <- read.csv("old_fitter/06_summary_val_DE.csv") %>%
  mutate(
    fit = "location_fit",
    Cultivar = Cultivar %>%
      sub("^data_", "", .) %>%            # remove "data_" prefix
      gsub("([a-z])([A-Z])", "\\1 \\2", .) # add space to make like the apple dataframe 
)

val_apple <- read.csv("old_fitter/04_summary_val_both.csv") %>%
  mutate(fit = "apple_fit") %>%
  dplyr::select(-Species) # drop species column 

# bind together 

prediction_df <- rbind(val_apple,
                       val_de,
                       val_es) %>%
  mutate(Cultivar = recode(Cultivar, "DelaRiega" = "De la Riega"))

# order cultivars to have German cultivars together and Spanish ones together

spanish_cultivars <- c("Blanquina", "Clara", "Collaos", "Coloradona", "De la Riega",
                      "Perezosa", "Perico", "Raxao","Teorica", "Verdialona", "Xuanina")

german_cultivars <- c("Berlepsch", "Cox Orange", "Golden Delicious","James Grieve" , 
                       "Roter Boskoop")

ordered_cultivars <- c(german_cultivars, spanish_cultivars)

# Modify the data to use the ordered factor
prediction_df <- prediction_df %>%
  mutate(Cultivar = factor(Cultivar, levels = ordered_cultivars))


## calculate IQR for validation + calibration because of the small size sample

## load prediction dataframe with both calibration and validation datasets

# add fit column to distinguish between the different fits when calculating iqr

prediction_df_apple_old <- read.csv("old_fitter/prediction_df_apple.csv")%>%
  mutate(fit = "apple_fit") %>%
  dplyr::select(-Species) # drop species column 

prediction_df_de_old <- read.csv("old_fitter/prediction_df_de.csv")%>%
  mutate(fit = "location_fit") %>%
  mutate(Location = "Germany")

prediction_df_es_old <- read.csv("old_fitter/prediction_df_es.csv")%>%
  mutate(fit = "location_fit")  %>%
  mutate(Location = "Spain")

# put together 

prediction_df_old <- bind_rows(prediction_df_apple_old,
                               prediction_df_de_old,
                               prediction_df_es_old)  %>%
  mutate(Cultivar = recode(Cultivar, "DelaRiega" = "De la Riega")) 

## Add id to distinguish between the different combinations 
prediction_df_id_old <- prediction_df_old %>%
  mutate(id = paste(Cultivar,Location, fit, sep = "_"))


iqr_df <- prediction_df_old %>% 
  group_by(Cultivar,fit) %>% 
  summarise(iqr = IQR(pheno)) %>%
  ungroup()
  

# make performance dataframe 

performance_df<- prediction_df %>%
  group_by(Cultivar,fit) %>% 
  summarise(RMSE = chillR::RMSEP(Predicted,pheno),
            Bias = mean(Predicted-pheno)) %>%
  ungroup() %>% 
  #add the iqr information
  merge(iqr_df, by = c("Cultivar", "fit")) %>% 
  mutate(rpiq = iqr / RMSE) %>% 
  reshape2::melt(id.vars = c("Cultivar","fit")) %>% 
  reshape2::dcast(Cultivar~ fit + variable,value.var = "value") 

# summarize per location to compare overall performance between German and Spanish 

performance_by_location <- prediction_df %>%
  mutate(location = case_when(
    Cultivar %in% c("Berlepsch", "Cox Orange", "Golden Delicious", "James Grieve", "Roter Boskoop") ~ "Germany",
    TRUE ~ "Spain"
  )) %>%
  group_by(location, fit) %>%
  summarise(RMSE = RMSEP(Predicted,pheno), RPIQ = RPIQ(Predicted,pheno),
            Bias = mean(Predicted-pheno)) %>% 
  reshape2::melt(id.vars = c("location","fit")) %>% 
  reshape2::dcast(location~ fit + variable,value.var = "value")

#save
write.csv(performance_by_location,"results/model_performance/performance_per_location_old_fitter.csv", row.names = FALSE)

# plot 

p_performance_old_fitter <- prediction_df %>%
  
  ggplot(aes(x = pheno, y = Predicted)) +
  
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed', color = "gray40") +
  
  geom_point(aes(color = fit),
              size = 4,
              alpha = 0.7) +

  scale_color_manual(name = "Calibration approaches:",
                     values = c("apple_fit" = "#E64B35",
                                "location_fit" = "#4DBBD5"),
                     labels = c ("Species-specific", "Location-specific")) +
  
  ylim(c(60,180))+ xlim(c(60,180))+
  
  geom_text(data = performance_df, aes(x = 75, y = 145, label = paste0('RMSE: ', format(round(apple_fit_RMSE, digits = 1), nsmall = 1),
                                                                          ' (', format(round(location_fit_RMSE,1), nsmall = 1), ')')), hjust = 0, size = 3)+
  geom_text(data = performance_df, aes(x = 75, y = 140, label = paste0('RPIQ: ', format(round(apple_fit_rpiq, digits = 1), nsmall = 1),
                                                                          ' (', format(round(location_fit_rpiq,1), nsmall = 1), ')')), hjust = 0, size = 3)+
  geom_text(data = performance_df, aes(x = 75, y = 135, label = paste0('bias: ', format(round(apple_fit_Bias, digits = 1), nsmall = 1),
                                                                          ' (', format(round(location_fit_Bias,1), nsmall = 1), ')')), hjust = 0, size = 3)+
  
  # Format axes
  scale_y_continuous(breaks = c(75, 96, 117, 138, 160), 
                     labels = c('15 Mar', '5 Apr', '26 Apr', '17 May', '8 Jun'))+
  scale_x_continuous(breaks = c(75, 96, 117, 138, 160), 
                     labels = c('15 Mar', '5 Apr', '26 Apr', '17 May', '8 Jun'))+
  # Labels
  labs(x = 'Observed Bloom Date',
       y = 'Predicted Bloom Date') +
  # Theming
  theme_bw(base_size = 14) +
  theme(
    # Panel formatting
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "gray20"),
    
    # Axis formatting
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(),
    
    # Legend formatting
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key = element_rect(fill = "white", color = NA, size = 10),
    legend.text = element_text(size = 10),
    legend.margin = margin(t = -10, unit = "pt"),
    plot.margin = margin (1.5,1.5,1.5,2),
    
    # Facet formatting
    strip.background = element_blank(), #element_rect(fill = "gray95"),
    strip.text = element_text(size = 10, face = "bold"),
    
    # Add proper spacing
    panel.spacing = unit(0.5, "lines")
  ) +
  facet_wrap(~ Cultivar, ncol = 4, scales = "fixed")

ggsave('plots/model_performance/model_performance_old_fitter_last_run_rpiq_fixed_final.jpeg',
       p_performance_old_fitter, height = 24, width = 19,
       dpi = 300, units = 'cm', device = 'jpeg')


