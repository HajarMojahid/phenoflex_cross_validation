
# load libraries

library(ggplot2)
library(patchwork)
library(tidyverse)
library(readxl)
library(chillR)
source("code/get_temp_response_plot.R")


##### load temperature data ######
#we need hourly temperature data
#Germany
#read temperature data
temp_de <- read.csv('Edu/weather_data_final.csv')

#get hourly data from daily data
hourtemps_de <- stack_hourly_temps(temp_de, latitude = 50.62)

#only take dataframe with data
hourtemps_de <- hourtemps_de[[1]]

#Spain 
#read temperature data
temp_es <- read_excel('MaxMin1978and2020.xlsx')

#make day and month numeric
temp_es$Month <- as.numeric(temp_es$Month)
temp_es$Day <- as.numeric(temp_es$Day)

#check if the temperature data is complete
check_temperature_record(temp_es)

#get hourly data from daily data
hourtemps_es <- stack_hourly_temps(temp_es, latitude=43.48)

#only take dataframe with data
hourtemps_es <- hourtemps_es[[1]]


#temperature for the response plot
temperature <- seq(-10,60, by = 0.1)


##### Global optimization ######

# loading parameters
param_old_es <- read.csv("old_fitter/08_summary_pars_ES.csv") %>%
  mutate(fit ="Location-specific, Spain")

param_old_de <- read.csv("old_fitter/09_summary_pars_DE.csv") %>%
  mutate(fit ="Location-specific, Germany")

param_apple_old <- read.csv("old_fitter/07_summary_pars_both.csv") %>%
  mutate(fit ="Species-specific")

param_all_old <- bind_rows(param_apple_old,
                           param_old_de,
                           param_old_es) %>%
  mutate(Par = recode(Par, "s" = "slope"))

### keeping last run 

param_all_old_last_run <- param_all_old %>%
  mutate(fitter = "GO") %>%
  # picking one set of parameter for each location
  filter(Fit_res == 'Fit_res_05', Cultivar %in% c('data_Berlepsch', 'data_Clara')) %>% 
  filter(!(fit == "Species-specific" & Cultivar != "data_Berlepsch")) %>% # additional filter for the apple fit  
  dplyr::select(-Cultivar, -Fit_res)


##### Enhanced global opt #####

#load parameters 

param_new_es <- read.csv("new_fitter/ES_fit/param_sum_es.csv") %>%
  mutate(fit ="Location-specific, Spain")

param_new_de <- read.csv("new_fitter/DE_fit/param_sum_de.csv") %>%
  mutate(fit ="Location-specific, Germany")

param_apple_new <- read.csv("new_fitter/de_es_fit/param_sum_de_es.csv") %>%
  mutate(fit ="Species-specific")

param_all_new <- bind_rows(param_apple_new,
                           param_new_de,
                           param_new_es) %>%
  mutate(Par = recode(Par, "s" = "slope"))

## keep one set of parameter per row.
## and the best performing run in calibration 
## see script model_performance_all.R for best performing rounds. 

param_all_new <- param_all_new %>%
  mutate(fitter="EF") %>%
  filter(Cultivar %in% c("Berlepsch","Clara")) %>%
  filter(!(fit == "Species-specific" & Cultivar != "Berlepsch")) %>% # additional filter for the apple fit
  filter( !(fit == "Location-specific, Spain" & Fit_res != "Fit_res_01"), #round 1 for Spanish fit
          !(fit == "Location-specific, Germany" & Fit_res != "Fit_res_01"), # round 1 for German fit 
          !(fit == "Species-specific" & Fit_res != "Fit_res_03") ) %>% # round 3 for apple fit
  dplyr::select(- Cultivar,-Fit_res) 


###### putting fitting procedures together together ####

#bind 

param_all <- bind_rows (
  param_all_old_last_run,
  param_all_new
)

## making temperature response

temperature_responses<- param_all  %>%
  group_by(fit,fitter) %>%
  summarise(temperature = temperature,
            chill_response = gen_bell(par = Value, temp_values = temperature ),
            heat_response = GDH_response(temperature, Value)) %>%
  reshape2::melt(id.vars = c('fit','fitter', 'temperature'))



##### plot by facets chill response #####

chill_response_plot <- ggplot(filter(temperature_responses, variable == "chill_response"),
                              aes(temperature, value,  color = fit)) +
  
  geom_line(size = 0.75, aes(col = fit)) +
  
  scale_x_continuous(limits = c(-10, 27)) +
  
  xlab("Temperature (°C)") +
  
  ylab("Arbitrary units") +
  
  facet_grid(factor(fitter, levels = c("GO","EF"),labels = c("G.~optimization","E.G.~optimization")) ~ 
               factor(variable, labels = c("Chill response"))) + 
  
  scale_color_manual(name = "Calibration approaches:",
                     values = c("Location-specific, Germany" = "#4DBBD5",
                                "Location-specific, Spain" = "#BCD979",
                                "Species-specific" = "#E64B35"),
                     labels = c("German location-specific",
                                "Spanish location-specific",
                                "Species-specific")) +
  
  theme_bw() +
  theme(axis.title = element_text(),
        axis.text = element_text(size = 8),
        strip.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

##### plot by facets heat response #####
heat_response_plot <- ggplot(filter(temperature_responses, variable == "heat_response"),
                              aes(temperature, value,  color = fit)) +
  
  geom_line(size = 0.75, aes(col = fit)) +
  
  xlab("Temperature (°C)") +
  
  ylab("Arbitrary units") +
  
  facet_grid(factor(fitter, levels = c("GO","EF"), labels = c("G.~optimization","E.G.~optimization")) ~ 
               factor(variable, labels = c("Heat~response")), labeller = label_parsed) +
  
  scale_color_manual(name = "Calibration approaches:",
                     values = c("Location-specific, Germany" = "#4DBBD5",
                                "Location-specific, Spain" = "#BCD979",
                                "Species-specific" = "#E64B35"),
                     labels = c("German location-specific",
                                "Spanish location-specific",
                                "Species-specific")) +
  
  theme_bw() +
  theme(axis.title = element_text(),
        axis.text = element_text(size = 8),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.position = "bottom",
        legend.text = element_text(size = 8))


##### use patchwork to put the plots together #####

(chill_response_plot + heat_response_plot) +
  plot_layout(guides = "collect", axis_titles = "collect") & 
  theme(legend.title = element_text(size = 6),
        legend.position = "bottom",
        legend.key = element_rect(size = 6),
        legend.text = element_text(size = 6),
        legend.margin = margin(t = -10, unit = "pt"),
        plot.margin = margin (1,1,1,1))

## save 
ggsave("plots/temperature_response/temp_responses_final.margin_ready.png",
       width = 12, height = 10, units = "cm", dpi = 600)

